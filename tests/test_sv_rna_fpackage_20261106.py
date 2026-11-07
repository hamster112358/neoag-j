"""Package F tests: SV chain + RNA allele counting (B-E review B4-*, C3-*).

Each test names the review item it locks down:

B4-01  support FORMAT keys for GRIDSS/SvABA; absence of support keys is
       "not measured", never a measured zero.
B4-06  FILTER="." is "no decision made": excluded by default but ledgered,
       kept on request; genuinely-failed records stay excluded.
B4-07  tumor specificity distinguishes matched-normal-zero-support from
       unassessed normal (single-sample VCF cannot prove somatic).
B4-08  fusion joins honour the breakend orientation of the ALT column;
       orientation-impossible joins are rejected before peptide generation.
B4-09  a breakpoint snapped to the nearest CDS base is an APPROXIMATE
       projection, forcing low reconstruction confidence.
B4-10  SV clustering is the transitive closure of pairwise compatibility,
       identical for every input order; event VAF and confidence aggregate
       only over members compatible with the representative (round 2).
B4-14  predictor provenance follows the actual producer of the result file,
       not the global binding_stub flag.
B4-16  a gene missing from the expression table blanks event_expression
       instead of encoding absence as TPM 0.

(C3-01..C3-03 RNA pileup hardening is covered in
test_rna_allele_counts_script.py.)
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import pytest

from neoag.sv.phase1 import (
    _event_expression,
    _raw_event_from_sv,
    _support_rollup,
    _sv_event_row,
    build_sv_phase1_raw,
)
from neoag.sv.sv_filter import score_cluster_confidence
from neoag.sv.protein_reconstruct import reconstruct_cluster_protein
from neoag.sv.reference import FastaReference, GtfAnnotation
from neoag.sv.score_pipeline import run_sv_score
from neoag.sv.sv_callset import SVRecord, evidence_from_sample, parse_vcf_records
from neoag.sv.sv_merge import SVCluster, cluster_sv_records
from neoag.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"
MHCFLURRY_FIXTURE = ROOT / "data" / "fixtures" / "mhcflurry_predictions.csv"

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
)


def _mini_run(tmp_path, **overrides):
    kwargs = dict(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions_exact.tsv",
        expressed_products_tsv=FX / "expressed_products.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    kwargs.update(overrides)
    outdir = kwargs.pop("outdir", tmp_path / "adapter")
    return build_sv_phase1_raw(outdir=outdir, **kwargs)


def _bnd(alt: str, *, svtype: str = "BND", chrom1="chr1", pos1=60, chrom2="chr2", pos2=31) -> SVRecord:
    return SVRecord(
        record_id="bnd", caller="GRIDSS",
        chrom1=chrom1, pos1=pos1, chrom2=chrom2, pos2=pos2,
        svtype=svtype, alt=alt,
    )


# --------------------------------------------------------------------------- B4-01


def test_gridss_format_maps_to_measured_support():
    # SUR/SU split reads, PAI/PAIR spanning pairs, VD assembled reads,
    # PSC/PSR reference-side support bounding the local depth.
    sr, pe, alt, depth, observed = evidence_from_sample(
        {"GT": "1/1", "PAI": "3", "PAIR": "5", "PSC": "10", "PSR": "12", "SUR": "2", "SU": "4", "VD": "20"}
    )
    assert (sr, pe, alt, depth, observed) == (6, 8, 20, 36, True)


def test_svaba_support_wins_over_generic_pe():
    sr, pe, alt, depth, observed = evidence_from_sample(
        {"SPLIT_READS": "7", "PAIRED_READS": "9", "GENOME_DEPTH": "42"}
    )
    assert (sr, pe, alt, depth, observed) == (7, 9, 16, 42, True)
    # SvABA FORMATs also carry a PE field with a different meaning; the
    # explicit SPLIT_READS/PAIRED_READS counts must win.
    sr2, pe2, alt2, _, _ = evidence_from_sample(
        {"SPLIT_READS": "7", "PAIRED_READS": "9", "PE": "100"}
    )
    assert (sr2, pe2, alt2) == (7, 9, 16)


def test_absent_support_keys_are_not_measured_zero():
    sr, pe, alt, depth, observed = evidence_from_sample({"GT": "1/1", "DP": "30"})
    assert observed is False
    assert (sr, pe, alt) == (0, 0, 0)
    # A caller that reports support keys with zero values IS making a
    # negative measurement - that must stay observed.
    sr, pe, alt, depth, observed = evidence_from_sample({"SR": "5,0", "PR": "5,0", "DP": "20"})
    assert observed is True
    assert (sr, pe, alt, depth) == (0, 0, 0, 20)


def test_support_status_per_sample_from_vcf_layout(tmp_path):
    single = tmp_path / "single.vcf"
    single.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
        "chr1\t60\tbnd1\tN\tN]chr2:31]\t.\tPASS\tSVTYPE=BND\tPR:SR:DP\t12,6:8,4:60\n",
        encoding="utf-8",
    )
    recs = parse_vcf_records(single, "GRIDSS2", require_explicit_sample_names=False)
    assert recs[0].tumor_support_status == "MEASURED"
    # No normal column at all: "SAMPLE_ABSENT", never an implied normal-zero.
    assert recs[0].normal_support_status == "SAMPLE_ABSENT"

    no_fmt = tmp_path / "nofmt.vcf"
    no_fmt.write_text(
        VCF_HEADER + "chr1\t60\ta\tN\tN]chr2:31]\t.\tPASS\tSVTYPE=BND\t.\t.\t.\n",
        encoding="utf-8",
    )
    recs = parse_vcf_records(no_fmt, "C", tumor_sample_name="TUMOR", normal_sample_name="NORMAL")
    assert recs[0].tumor_support_status == "FORMAT_ABSENT"

    gt_only = tmp_path / "gt.vcf"
    gt_only.write_text(
        VCF_HEADER + "chr1\t60\tb\tN\tN]chr2:31]\t.\tPASS\tSVTYPE=BND\tGT\t./.\t./.\n",
        encoding="utf-8",
    )
    recs = parse_vcf_records(gt_only, "C", tumor_sample_name="TUMOR", normal_sample_name="NORMAL")
    assert recs[0].tumor_support_status == "NO_EVIDENCE_KEYS"


def test_cluster_rollup_prefers_measured():
    measured = SVRecord("m", "C", "chr1", 1, "chr1", 2, "DEL")
    measured.tumor_support_status = "NO_EVIDENCE_KEYS"
    other = SVRecord("o", "C", "chr1", 1, "chr1", 2, "DEL")
    other.tumor_support_status = "MEASURED"
    assert _support_rollup([measured, other], "tumor_support_status") == "MEASURED"
    assert _support_rollup([measured], "tumor_support_status") == "NO_EVIDENCE_KEYS"


# --------------------------------------------------------------------------- B4-06


def _mixed_filter_vcf(tmp_path) -> Path:
    path = tmp_path / "mix.vcf"
    path.write_text(
        VCF_HEADER
        + "chr1\t60\tbnd1\tN\tN]chr2:31]\t.\tPASS\tSVTYPE=BND\tPR:SR:DP\t10,0:5,0:50\t12,6:8,4:60\n"
        + "chr1\t61\tbnd2\tN\tN]chr2:32]\t.\t.\tSVTYPE=BND\tPR:SR:DP\t10,0:5,0:50\t12,6:8,4:60\n"
        + "chr1\t62\tbnd3\tN\tN]chr2:33]\t.\tLowQual\tSVTYPE=BND\tPR:SR:DP\t10,0:5,0:50\t12,6:8,4:60\n",
        encoding="utf-8",
    )
    return path


def test_unfiltered_records_excluded_but_ledgered(tmp_path):
    ledger: list[dict[str, str]] = []
    recs = parse_vcf_records(
        _mixed_filter_vcf(tmp_path), "GRIDSS2",
        tumor_sample_name="TUMOR", normal_sample_name="NORMAL", filter_ledger=ledger,
    )
    assert [r.record_id for r in recs] == ["bnd1"]
    assert [(row["record_id"], row["decision"]) for row in ledger] == [
        ("bnd2", "EXCLUDED_FILTER_NOT_ASSESSED"),
        ("bnd3", "EXCLUDED_FILTER_NOT_PASS"),
    ]


def test_unfiltered_include_policy_keeps_records(tmp_path):
    ledger: list[dict[str, str]] = []
    recs = parse_vcf_records(
        _mixed_filter_vcf(tmp_path), "GRIDSS2",
        tumor_sample_name="TUMOR", normal_sample_name="NORMAL",
        unfiltered_policy="include", filter_ledger=ledger,
    )
    # FILTER="." is kept for research review; a real FAIL still is not.
    assert [r.record_id for r in recs] == ["bnd1", "bnd2"]
    assert [(row["record_id"], row["decision"]) for row in ledger] == [
        ("bnd3", "EXCLUDED_FILTER_NOT_PASS"),
    ]


def test_phase1_writes_filter_exclusions_ledger(tmp_path):
    out = _mini_run(tmp_path)
    ledger_path = Path(out["sv_vcf_filter_exclusions"])
    assert ledger_path.is_file()
    # The all-PASS fixture yields a header-only ledger: the audit trail
    # exists whether or not anything was dropped.
    assert read_tsv(ledger_path) == []

    out2 = _mini_run(
        tmp_path,
        sv_vcfs=[_mixed_filter_vcf(tmp_path)],
        outdir=tmp_path / "adapter2",
        rna_junction_tsv=None,
        expressed_products_tsv=None,
    )
    rows = read_tsv(out2["sv_vcf_filter_exclusions"])
    assert [(r["record_id"], r["decision"]) for r in rows] == [
        ("bnd2", "EXCLUDED_FILTER_NOT_ASSESSED"),
        ("bnd3", "EXCLUDED_FILTER_NOT_PASS"),
    ]


# --------------------------------------------------------------------------- B4-07


def _svrow(**extra):
    row = {
        "event_id": "E1", "sample_id": "S", "callers": "GRIDSS2",
        "gene1": "GENE1", "gene2": "GENE2", "effect_class": "SV_Fusion",
        "sv_vaf_like": "0.2000", "tumor_alt_support": 10, "tumor_local_depth": 50,
        "normal_alt_support": 0, "normal_local_depth": 40,
        "chrom1": "chr1", "pos1": 60, "chrom2": "chr2", "pos2": 31,
    }
    row.update(extra)
    return row


def test_specificity_requires_matched_normal():
    # No normal evidence at all -> capped specificity, status says why.
    row = _raw_event_from_sv(_svrow(), None, {}, {}, "sv_wgs_phase1")
    assert (row["tumor_specificity"], row["tumor_specificity_status"]) == ("0.8500", "NORMAL_UNASSESSED")
    # A measured normal with zero support keeps full specificity.
    row = _raw_event_from_sv(_svrow(normal_support_status="MEASURED"), None, {}, {}, "sv_wgs_phase1")
    assert (row["tumor_specificity"], row["tumor_specificity_status"]) == ("1.0000", "MATCHED_NORMAL_ZERO_SUPPORT")
    # Support seen in the normal marks the event as possibly germline.
    row = _raw_event_from_sv(
        _svrow(normal_alt_support=3, normal_support_status="MEASURED"), None, {}, {}, "sv_wgs_phase1",
    )
    assert (row["tumor_specificity"], row["tumor_specificity_status"]) == ("0.7000", "MATCHED_NORMAL_SUPPORT_DETECTED")


def test_mini_pair_reports_measured_normal(tmp_path):
    out = _mini_run(tmp_path)
    events = read_tsv(out["sv_events_full"])
    assert events[0]["tumor_support_status"] == "MEASURED"
    assert events[0]["normal_support_status"] == "MEASURED"
    raw = read_tsv(out["raw_events"])
    assert raw[0]["tumor_specificity"] == "1.0000"
    assert raw[0]["tumor_specificity_status"] == "MATCHED_NORMAL_ZERO_SUPPORT"


# --------------------------------------------------------------------------- B4-08


def _annotation(tmp_path=None, *, mate_strand: str = "+", cds1_end: int | None = None):
    ref = FastaReference(FX / "mini_ref.fa")
    if mate_strand == "+" and cds1_end is None:
        return ref, GtfAnnotation(FX / "mini.gtf")
    lines = []
    for line in (FX / "mini.gtf").read_text().splitlines():
        if not line.strip():
            continue
        cells = line.split("\t")
        if mate_strand != "+" and cells[0] == "chr2":
            cells[6] = mate_strand
        if cds1_end is not None and cells[0] == "chr1" and cells[2] == "CDS":
            cells[4] = str(cds1_end)
        lines.append("\t".join(cells))
    path = tmp_path / "mini_variant.gtf"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return ref, GtfAnnotation(path)


def test_exact_geometry_accepts_compatible_join():
    # N]chr2:31] : local 5' side first, mate read rightward ("+"): both
    # fixture transcripts are "+", so the join is possible.
    ref, ann = _annotation()
    protein, meta = reconstruct_cluster_protein(SVCluster([_bnd("N]chr2:31]")]), "E1", "S", ann, ref)
    assert protein is not None
    assert meta["breakpoint_orientation_status"] == "ORIENTATION_EXACT"
    assert (meta["breakpoint1_cds_mapping"], meta["breakpoint2_cds_mapping"]) == ("EXACT", "EXACT")
    assert protein.gene == "GENE1::GENE2"
    assert protein.reconstruction_confidence == "medium"


def test_impossible_orientation_rejected_before_protein():
    # N[chr2:31[ demands the mate fragment read leftward ("-"), but GENE2's
    # transcript is "+": no spliced product can carry this junction.
    ref, ann = _annotation()
    protein, meta = reconstruct_cluster_protein(SVCluster([_bnd("N[chr2:31[")]), "E1", "S", ann, ref)
    assert protein is None
    assert meta["breakpoint_orientation_status"] == "JOIN_REJECTED_ORIENTATION"
    assert meta["reconstruction_reason"].startswith("impossible_join_orientation:alt=N[chr2:31[")


def test_mate_minus_strand_accepts_reverse_join(tmp_path):
    ref, ann = _annotation(tmp_path, mate_strand="-")
    protein, meta = reconstruct_cluster_protein(SVCluster([_bnd("N[chr2:31[")]), "E1", "S", ann, ref)
    assert protein is not None
    assert meta["breakpoint_orientation_status"] == "ORIENTATION_EXACT"
    assert protein.gene == "GENE1::GENE2"


def test_mate_first_alt_swaps_partner_roles():
    # "]chr1:60]N" at chr2:31 is the SAME junction as the mini record,
    # written mate-first: the 5' prefix must come from chr1 (GENE1), not
    # from the record's own chrom1 (GENE2).
    ref, ann = _annotation()
    protein, meta = reconstruct_cluster_protein(
        SVCluster([_bnd("]chr1:60]N", chrom1="chr2", pos1=31, chrom2="chr1", pos2=60)]),
        "E1", "S", ann, ref,
    )
    assert protein is not None
    assert protein.gene == "GENE1::GENE2"
    assert meta["breakpoint_orientation_status"] == "ORIENTATION_EXACT"


def test_symbolic_alt_orientation_unassessed():
    # <TRA> carries no coordinates, so no geometry can be checked; the join
    # proceeds with the old layout but can never claim medium confidence.
    ref, ann = _annotation()
    protein, meta = reconstruct_cluster_protein(SVCluster([_bnd("<TRA>", svtype="TRA")]), "E1", "S", ann, ref)
    assert protein is not None
    assert meta["breakpoint_orientation_status"] == "ORIENTATION_UNASSESSED"
    assert protein.reconstruction_confidence == "low"


# --------------------------------------------------------------------------- B4-09


def test_intronic_breakpoint_marks_approximate_projection(tmp_path):
    # GENE1's CDS ends at 100; DEL breakpoints at 120/130 are intronic and
    # get snapped to the last CDS base -> APPROXIMATE + low confidence.
    ref, ann = _annotation(tmp_path, cds1_end=100)
    record = SVRecord("del", "GRIDSS", "chr1", 120, "chr1", 130, "DEL")
    protein, meta = reconstruct_cluster_protein(SVCluster([record]), "E1", "S", ann, ref)
    assert meta["breakpoint1_cds_mapping"] == "APPROXIMATE_20BP"
    assert meta["breakpoint2_cds_mapping"] == "APPROXIMATE_30BP"
    assert meta["breakpoint_orientation_status"] == "NOT_APPLICABLE_LOCAL_SV"
    assert protein is not None
    assert protein.reconstruction_confidence == "low"
    assert "APPROXIMATE: hypothesis requires RNA evidence" in protein.reconstruction_reason


# --------------------------------------------------------------------------- B4-10


def _del(rid: str, p1: int, p2: int) -> SVRecord:
    return SVRecord(rid, "caller", "chr1", p1, "chr1", p2, "DEL")


def test_clustering_is_transitive_and_order_independent():
    # A~B and B~C within 200bp, A and C 300bp apart: every input order must
    # yield ONE cluster of all three (transitive closure), not an order-
    # dependent greedy partition.
    a, b, c = _del("A", 100, 150), _del("B", 250, 300), _del("C", 400, 450)
    for perm in itertools.permutations([a, b, c]):
        clusters = cluster_sv_records(list(perm), distance=200)
        assert len(clusters) == 1, perm
        assert sorted(clusters[0].record_ids) == ["A", "B", "C"]


def test_distant_records_stay_separate_clusters():
    a = _del("A", 100, 150)
    far = _del("D", 1000, 1050)
    clusters = cluster_sv_records([a, far], distance=200)
    assert len(clusters) == 2


# --------------------------------------------------------------------------- B4-14


@pytest.fixture()
def _skip_predictor_coverage(monkeypatch):
    # The SV fixtures' peptide space and the predictor stub fixtures were
    # never paired, so the core-predictor coverage cross-check would abort
    # before provenance is recorded.  These two tests pin WHERE the
    # mhcflurry provenance status comes from, not predictor coverage, so the
    # cross-check is silenced for them only.
    monkeypatch.setattr(
        "neoag.presentation._validate_core_predictor_coverage",
        lambda *args, **kwargs: None,
    )


def test_provided_mhcflurry_stays_real_under_binding_stub(tmp_path, _skip_predictor_coverage):
    adapter = _mini_run(tmp_path)
    out = run_sv_score(
        outdir=tmp_path / "scored",
        profile_name_or_path="sv_wgs_phase1",
        sample_id="SVMINI",
        raw_events=Path(adapter["raw_events"]),
        raw_peptides=Path(adapter["raw_peptides"]),
        mhcflurry=MHCFLURRY_FIXTURE,
        expression=FX / "expression.tsv",
        normal_expression=FX / "normal_expression.tsv",
        normal_hla_ligands=FX / "normal_hla_ligands.tsv",
        binding_stub=True,
        immunogenicity_stub=True,
    )
    assert out  # scoring completed
    prov = json.loads((tmp_path / "scored" / "provenance.sv_score.json").read_text(encoding="utf-8"))
    # A user-supplied evidence file stays real even while OTHER predictors
    # are stubbed by binding_stub=True (review B4-14).
    assert prov["tools"]["mhcflurry"]["status"] == "real"


def test_self_generated_stub_output_marked_stub(tmp_path, _skip_predictor_coverage):
    adapter = _mini_run(tmp_path)
    run_sv_score(
        outdir=tmp_path / "scored",
        profile_name_or_path="sv_wgs_phase1",
        sample_id="SVMINI",
        raw_events=Path(adapter["raw_events"]),
        raw_peptides=Path(adapter["raw_peptides"]),
        expression=FX / "expression.tsv",
        normal_expression=FX / "normal_expression.tsv",
        normal_hla_ligands=FX / "normal_hla_ligands.tsv",
        binding_stub=True,
        immunogenicity_stub=True,
    )
    prov = json.loads((tmp_path / "scored" / "provenance.sv_score.json").read_text(encoding="utf-8"))
    # The pipeline-generated stub CSV must be labelled stub, never dressed
    # up as a passthrough of a real prediction.
    assert prov["tools"]["mhcflurry"]["status"] == "stub"


# --------------------------------------------------------------------------- B4-16


def test_event_expression_vocabulary():
    both = {"GENE1": 12.0, "GENE2": 8.0}
    assert _event_expression("SV_Fusion", "GENE1", "GENE2", both, 0) == (8.0, "GENE_TPM_MEASURED")
    # A fusion is limited by its lower partner, but a partner with no row
    # contributes nothing - partial evidence, not zero.
    assert _event_expression("SV_Fusion", "GENE1", "GENE2", {"GENE1": 12.0}, 0) == (12.0, "GENE_ONLY_PARTIAL")
    assert _event_expression("SV_Fusion", "GENE1", "GENE2", {}, 0) == (None, "UNASSESSED_ID_NOT_MAPPED")
    assert _event_expression("SV_Noncoding", "", "", {}, 0) == (None, "UNASSESSED")


def test_missing_expression_rows_blank_the_value(tmp_path):
    empty_expr = tmp_path / "empty_expression.tsv"
    empty_expr.write_text("gene\tTPM\n", encoding="utf-8")
    out = _mini_run(tmp_path, expression_tsv=empty_expr)
    raw = read_tsv(out["raw_events"])
    # Absence is blank + status, never a presented TPM 0.0.
    assert raw[0]["event_expression"] == ""
    assert raw[0]["expression_evidence_status"] == "UNASSESSED_ID_NOT_MAPPED"


# ------------------------------------------------- B4-10 round-2 follow-up --

def test_b4_10_transitive_chain_counts_only_representative_compatible_evidence():
    """Operator test 1 (2026-11-07): a transitive SV cluster whose remote
    member has distinct breakpoints AND an incompatible orientation, with
    different read counts/depths.  Shows the representative, all member
    evidence, the calculated VAF, the confidence score and the event
    consequences; the decision is that only representative-compatible
    evidence contributes to the event."""
    # A: GRIDSS DEL chr1:100000-200000 (+/+), 40 alt / 50 depth, CI 10 bp.
    a = SVRecord(
        "A", "GRIDSS", "chr1", 100_000, "chr1", 200_000, "DEL",
        strand1="+", strand2="+", tumor_sr=12, tumor_pe=30,
        tumor_alt_support=40, tumor_local_depth=50,
        cipos="0,10", ciend="0,10", tumor_support_status="MEASURED",
    )
    # B: Manta at B:100150-200120, strand unknown (./.) - within 200 bp of
    # both A and C, so it chains them although A and C never agree.
    b = SVRecord(
        "B", "Manta", "chr1", 100_150, "chr1", 200_120, "DEL",
        strand1=".", strand2=".", tumor_sr=0, tumor_pe=1,
        tumor_alt_support=2, tumor_local_depth=30, tumor_support_status="MEASURED",
    )
    # C: SvABA DEL chr1:100300-200300 (-/-): 300 bp from A (beyond the 200 bp
    # window) with opposite orientation, but 150/180 bp from B - the classic
    # transitive chain member.  It carries the LARGEST raw counts so an
    # unrestricted aggregation would let it dominate the event.
    c = SVRecord(
        "C", "SvABA", "chr1", 100_300, "chr1", 200_300, "DEL",
        strand1="-", strand2="-", tumor_sr=18, tumor_pe=90,
        tumor_alt_support=200, tumor_local_depth=210,
        sniffles_support=25, cipos="0,3", ciend="0,3",
        tumor_support_status="MEASURED",
    )
    for permutation in itertools.permutations([a, b, c]):
        clusters = cluster_sv_records(list(permutation), distance=200)
        assert len(clusters) == 1
        cluster = clusters[0]
        assert cluster.representative.caller == "GRIDSS"  # priority + support
        assert sorted(cluster.record_ids) == ["A", "B", "C"]
        assert [r.caller for r in cluster.compatible_records] == sorted(
            [r.caller for r in cluster.compatible_records]
        )
    cluster = cluster_sv_records([b, c, a], distance=200)[0]
    assert [r.record_id for r in cluster.compatible_records] == ["A", "B"]
    assert [r.record_id for r in cluster.incompatible_records] == ["C"]

    # Member evidence as the operator asked to see it: A 40/50 (0.80),
    # B 2/30, C 200/210 (0.952).  Aggregated over ALL members (the old
    # behaviour) the event VAF would read 200/210 = 0.952 and the caller
    # count would claim 3 callers for an event only two callers described.
    old_style_vaf = max(a.tumor_alt_support, b.tumor_alt_support, c.tumor_alt_support) / max(
        a.tumor_local_depth, b.tumor_local_depth, c.tumor_local_depth
    )
    assert old_style_vaf == pytest.approx(0.952381, abs=1e-4)

    conf = score_cluster_confidence(cluster)
    # Compatible members: GRIDSS+Manta -> caller_score 2/3; GRIDSS split/
    # pair support (12/30) alone clamps support_score to 1.0; GRIDSS CI 10 bp
    # -> precision_score 1.0; no normal data, no RNA bonus.
    assert conf.score == pytest.approx(0.30 * (2 / 3) + 0.35 * 1.0 + 0.20 * 1.0, abs=1e-9)
    assert conf.score == pytest.approx(0.75)
    # Unrestricted aggregation would have scored C's 3-caller support and
    # long-read evidence into this event:
    assert conf.score < 0.30 * 1.0 + 0.35 * 1.0 + 0.20 * 1.0  # 0.85 old-style
    assert "incompatible_cluster_members_excluded" in conf.reason

    row = _sv_event_row(cluster, "SV1", "S1", {}, conf, None, 0)
    # Event coordinates/identity come from the single self-consistent
    # representative record.
    assert (row["chrom1"], row["pos1"], row["strand1"]) == ("chr1", 100_000, "+")
    assert (row["chrom2"], row["pos2"], row["strand2"]) == ("chr1", 200_000, "+")
    # VAF now describes the representative's event: 40/50, not C's 0.952.
    assert row["sv_vaf_like"] == "0.8000"
    assert row["tumor_alt_support"] == 40 and row["tumor_local_depth"] == 50
    assert row["tumor_sr"] == 12 and row["tumor_pe"] == 30  # not C's 18/90
    assert row["sniffles_support"] == 0  # C's 25 long reads excluded
    assert row["breakpoint_precision_bp"] == 10  # GRIDSS CI, not C's 3 bp
    assert row["callers"] == "GRIDSS,Manta" and row["caller_count"] == 2
    assert row["cluster_all_callers"] == "GRIDSS,Manta,SvABA"
    # Nothing hidden: all three members stay listed, C named as excluded.
    assert row["record_ids"] == "A,B,C"
    assert row["cluster_record_count"] == 3 and row["compatible_record_count"] == 2
    assert row["excluded_incompatible_record_ids"] == "C"
    assert row["excluded_incompatible_tumor_alt_max"] == 200
    assert row["event_evidence_basis"] == "REPRESENTATIVE_COMPATIBLE@v1"
    assert row["event_confidence_tier"] == conf.tier == "Tier1"

    # Event consequence layer: the phase-1 raw event inherits the corrected
    # VAF (clonality prior) and honest normal status (specificity cap).
    raw = _raw_event_from_sv(row, None, {}, {}, "default")
    assert float(raw["clonality"]) == 1.0  # 0.8/0.35 saturates
    assert raw["tumor_specificity"] == "0.8500"  # unassessed normal cap (B4-07)

    # A fully compatible cluster keeps the old basis label and discloses no
    # exclusions.
    plain = cluster_sv_records([a, b], distance=200)[0]
    plain_conf = score_cluster_confidence(plain)
    plain_row = _sv_event_row(plain, "SV0", "S1", {}, plain_conf, None, 0)
    assert plain_row["event_evidence_basis"] == "ALL_MEMBERS_COMPATIBLE"
    assert plain_row["excluded_incompatible_record_ids"] == ""
    assert "incompatible_cluster_members_excluded" not in plain_conf.reason
