"""Tests for work package G-d (2026-11-06): C4-01/03/04 + B4-02/03/05/11/12/13/15/17/18.

C4-01 find_cn no longer resolves overlapping segments by file order and
    records the chosen span plus an overlap-consistency status.
C4-04 _purity_bounds parses structured keys first and free-text ranges only
    when the whole string is exactly a numeric pair (no prose, no %).
B4-02 FastaReference fetches regions by byte-level random access (.fai or a
    streaming offset scan) instead of loading the genome into RAM.
B4-03 normal_hla_ligand_overlap reports "not_assessed" without a catalogue.
B4-05 the SV raw-event builder emits no duplicate keys and labels its
    heuristic event components via event_component_status.
B4-11 external CONFIRMED products blocked by junction QC land in a ledger.
B4-12 the SV scorer accepts the main branch's auxiliary predictors and
    records per-predictor evidence completeness.
B4-13 projected partner-junction matches of adjacent genes become a review
    flag instead of counted RNA support.
B4-15 capture interpretability is a joint two-sided tier (HIGH needs both).
B4-17 capture annotation is a declared typed SVCluster field, not setattr.
B4-18 gene-pair junction keys are exploratory only; read support transfers
    through exact event/coordinate keys alone.
"""

from __future__ import annotations

import ast
import dataclasses
import gzip
import json
from pathlib import Path

import pytest

from neoag.ccf_v2 import _purity_bounds, build_ccf_2, find_cn
from neoag.sv.evidence import JunctionReadEvidence, junction_positional_key, load_junction_reads
from neoag.sv.phase1 import build_sv_phase1_raw
from neoag.sv.score_pipeline import run_sv_score
from neoag.sv.peptide_builder import build_mhc1_peptides
from neoag.sv.protein_reconstruct import ProteinReconstruction
from neoag.sv.reference import FastaReference, GtfAnnotation, revcomp
from neoag.sv.sv_callset import SVRecord
from neoag.sv.sv_merge import SVCluster
from neoag.sv.wes_capture import CaptureRegions, annotate_cluster_capture, cluster_capture_annotation
from neoag.utils import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"
FXI = ROOT / "data" / "fixtures_sv_intron"

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _write_fasta(path: Path, seq: str, line_width: int = 70, name: str = "chr1") -> None:
    lines = [seq[i:i + line_width] for i in range(0, len(seq), line_width)]
    path.write_text(f">{name}\n" + "".join(line + "\n" for line in lines), encoding="utf-8")


def _write_fai(path: Path, seq_len: int, line_width: int = 70, name: str = "chr1") -> None:
    header = f">{name}\n"
    path.write_text(
        f"{name}\t{seq_len}\t{len(header.encode())}\t{line_width + 1}\t{line_width + 1}\n",
        encoding="utf-8",
    )


def _seq(n: int) -> str:
    return "".join("ACGT"[i % 4] for i in range(n))


# ---------------------------------------------------------------------------
# C4-01: find_cn segment selection
# ---------------------------------------------------------------------------


def _seg(start: int, end: int, cn: str, conf: str = "") -> dict[str, str]:
    row = {"chrom": "chr1", "start": str(start), "end": str(end), "total_cn": cn,
           "major_cn": cn, "minor_cn": "0"}
    if conf:
        row["cnv_confidence"] = conf
    return row


def test_find_cn_prefers_declared_quality_then_narrowness_not_file_order():
    wide_low = _seg(1, 100_000, "4", conf="low")
    narrow_high = _seg(140, 160, "2", conf="high")
    # File order puts the wrong (wide, low-quality) segment first.
    for rows in ([wide_low, narrow_high], [narrow_high, wide_low]):
        cn = find_cn("chr1", "150", rows)
        assert cn["total_cn"] == 2.0
        assert cn["cnv_segment_span"] == "140-160"
        assert cn["cnv_overlap_status"] == "DISCORDANT"
        # The genuine 2-vs-4 disagreement demotes the row even though the
        # high-quality segment won the selection.
        assert cn["cnv_confidence"] == "low"


def test_find_cn_equal_confidence_prefers_narrowest_segment():
    wide = _seg(1, 100_000, "2")
    mid = _seg(100, 500, "2")
    narrow = _seg(145, 155, "2")
    cn = find_cn("chr1", "150", [wide, mid, narrow])
    assert cn["cnv_segment_span"] == "145-155"
    assert cn["cnv_overlap_segments"] == "3"
    assert cn["cnv_overlap_status"] == "CONSISTENT"


def test_find_cn_reports_discordant_overlap_with_demoted_confidence():
    a = _seg(1, 200, "2", conf="high")
    b = _seg(100, 300, "6", conf="high")
    cn = find_cn("chr1", "150", [a, b])
    assert cn["cnv_overlap_status"] == "DISCORDANT"
    assert cn["cnv_overlap_segments"] == "2"
    assert cn["cnv_confidence"] == "low"  # genuine disagreement demotes the call
    # The winner is deterministic (same width -> coordinate order), not input order.
    assert cn == find_cn("chr1", "150", [b, a])


def test_find_cn_no_overlap_status_none_and_default_fields():
    cn = find_cn("chr1", "999", [_seg(1, 200, "2")])
    assert cn["matched"] is False
    assert cn["cnv_overlap_status"] == "NONE"
    assert cn["cnv_segment_span"] == ""
    assert cn["cnv_overlap_segments"] == "0"


def test_ccf_row_records_chosen_segment_and_overlap(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "mutation_source": "SNV",
        "chrom": "chr1", "pos": "150", "tumor_vaf": "0.25",
        "tumor_depth": "100", "tumor_alt_count": "25",
    }])
    write_tsv(purity, [{"purity": "0.5", "ploidy": "2"}])
    write_tsv(cnv, [
        {"chrom": "chr1", "start": "1", "end": "100000", "total_cn": "4",
         "major_cn": "4", "minor_cn": "0", "cnv_confidence": "low"},
        {"chrom": "chr1", "start": "140", "end": "160", "total_cn": "2",
         "major_cn": "1", "minor_cn": "1", "cnv_confidence": "high"},
    ])
    from neoag.config import load_profile
    rows = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv")
    assumptions = rows[0]["ccf_formula_assumptions"]
    assert "cnv_segment=140-160" in assumptions
    assert "cnv_overlap=DISCORDANT(n=2)" in assumptions
    # The narrow high-confidence segment won over the file-order-first one.
    assert rows[0]["total_cn"] == "2.0000"


# ---------------------------------------------------------------------------
# C4-04: purity bounds parsing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("info,expected", [
    ({"purity_low": "0.40", "purity_high": "0.60"}, (0.40, 0.60)),
    ({"purity_ci_low": "0.44", "purity_ci_high": "0.52"}, (0.44, 0.52)),
    ({"purity_range": "0.80-0.90"}, (0.80, 0.90)),
    ({"purity_range": "[0.45 - 0.55]"}, (0.45, 0.55)),
    ({"purity_range": "0.45 to 0.55"}, (0.45, 0.55)),
    ({"purity_range": "0.9-0.8"}, (0.8, 0.9)),  # reversed pair normalised
])
def test_purity_bounds_accepts_legit_forms(info, expected):
    assert _purity_bounds(info, 0.5) == expected


@pytest.mark.parametrize("info", [
    {"purity_range": "estimate 0.45 (Q30 base 0.98)"},   # prose
    {"purity_range": "90% tumor, purity ~0.4-0.6"},      # percent text
    {"purity_range": "0.5"},                              # single value
    {"purity_range": "n/a"},
    {},
])
def test_purity_bounds_falls_back_to_point_estimate(info):
    assert _purity_bounds(info, 0.5) == (0.5, 0.5)


def test_purity_bounds_structured_keys_beat_text():
    info = {"purity_low": "0.3", "purity_high": "0.4", "purity_range": "0.80-0.90"}
    assert _purity_bounds(info, 0.5) == (0.3, 0.4)


# ---------------------------------------------------------------------------
# B4-02: indexed FASTA access
# ---------------------------------------------------------------------------


def test_fasta_reference_plain_scan_and_fetch(tmp_path):
    seq = _seq(2000)
    fa = tmp_path / "ref.fa"
    _write_fasta(fa, seq)
    ref = FastaReference(fa)
    assert ref.contig_length("chr1") == 2000
    assert ref.fetch("chr1", 1, 30) == seq[:30]
    # Cross-line fetch (line width 70): spans two newline boundaries.
    assert ref.fetch("chr1", 65, 80) == seq[64:80]
    assert ref.fetch("chr1", 65, 80, "-") == revcomp(seq[64:80])
    # Clamping and alternate contig spelling.
    assert ref.fetch("chr1", 1995, 5000) == seq[1994:2000]
    assert ref.fetch("1", 1, 10) == seq[:10]
    with pytest.raises(KeyError):
        ref.fetch("chrZ", 1, 10)


def test_fasta_reference_fai_matches_scan(tmp_path):
    seq = _seq(1500)
    fa = tmp_path / "ref.fa"
    _write_fasta(fa, seq)
    scan = FastaReference(fa)
    _write_fai(fa.with_suffix(".fa.fai"), len(seq))
    indexed = FastaReference(fa)
    for start, end in [(1, 1), (70, 71), (700, 760), (1499, 1500)]:
        assert indexed.fetch("chr1", start, end) == scan.fetch("chr1", start, end)
        assert indexed.fetch("chr1", start, end) == seq[start - 1:end]


def test_fasta_reference_crlf_gz_and_missing_index_fallbacks(tmp_path):
    seq = _seq(300)
    crlf = tmp_path / "crlf.fa"
    lines = [seq[i:i + 70] for i in range(0, len(seq), 70)]
    crlf.write_text(">chr1\r\n" + "\r\n".join(lines) + "\r\n", encoding="ascii")
    # CRLF layout is only honest through the whole-file fallback; make sure the
    # scan path is not trusted with a .fai that would mis-stride it: without a
    # .fai the scan uses measured stride per line, which handles CRLF.
    ref = FastaReference(crlf)
    assert ref.fetch("chr1", 61, 85) == seq[60:85]

    gz = tmp_path / "ref.fa.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    gzref = FastaReference(gz)
    assert gzref.fetch("chr1", 1, 100) == seq[:100]


def test_fasta_reference_does_not_load_whole_genome(tmp_path):
    seq = _seq(50_000)
    fa = tmp_path / "big.fa"
    _write_fasta(fa, seq, line_width=100)
    ref = FastaReference(fa)
    assert ref._seqs is None  # scan mode: no materialised contigs
    ref.fetch("chr1", 1, 10)
    assert ref._seqs is None


# ---------------------------------------------------------------------------
# B4-03: normal HLA ligand overlap honesty
# ---------------------------------------------------------------------------


def _fusion_protein(seq: str, wt: str) -> ProteinReconstruction:
    return ProteinReconstruction(
        protein_sequence_id="PROT_T1", event_id="EVT_T1", sample_id="S1",
        gene="G1::G2", transcript_id="T1::T2", protein_type="SV_Fusion",
        protein_sequence=seq, wt_protein_sequence=wt,
        junction_aa_position=str(len(wt) - 4), novel_start_aa=str(len(wt) - 4),
        in_frame="yes",
        reconstruction_method="external_expressed_transcript:fixture",
        reconstruction_confidence="high",
    )


def _junction_peptides(normal_ligands=None):
    wt = "ABCDEFGHIJK"          # 11 aa
    novel = "LMNPQRS"
    prot = _fusion_protein(wt + novel, wt)
    return build_mhc1_peptides(
        [prot], ["HLA-A*02:01"], lengths=(8,),
        normal_ligands=normal_ligands,
        reference_peptides={"ZZZZZZZZ"},  # non-empty so ref-match is "assessed"
    )


def test_b4_03_absent_ligand_catalogue_is_not_assessed():
    peptides = _junction_peptides(normal_ligands=None)
    assert peptides and all(
        p.normal_hla_ligand_overlap == "not_assessed" for p in peptides
    )
    empty = _junction_peptides(normal_ligands=set())
    assert all(p.normal_hla_ligand_overlap == "not_assessed" for p in empty)


def test_b4_03_present_catalogue_reports_yes_no():
    peptides = _junction_peptides(normal_ligands={"ZZZZZZZZ"})
    assert peptides and all(p.normal_hla_ligand_overlap == "no" for p in peptides)
    some = peptides[0].peptide
    peptides2 = _junction_peptides(normal_ligands={some})
    hit = next(p for p in peptides2 if p.peptide == some)
    assert hit.normal_hla_ligand_overlap == "yes"


# ---------------------------------------------------------------------------
# B4-05: no duplicate dict keys + component provenance
# ---------------------------------------------------------------------------


def test_raw_event_builder_has_no_duplicate_dict_keys():
    tree = ast.parse(
        (ROOT / "src" / "neoag" / "sv" / "phase1.py").read_text(encoding="utf-8")
    )
    func = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "_raw_event_from_sv"
    )
    literals = [
        node for node in ast.walk(func)
        if isinstance(node, ast.Dict) and len(node.keys) > 20  # the base row dict
    ]
    assert literals, "base dict literal not found"
    for literal in literals:
        keys = [k.id for k in literal.keys if isinstance(k, ast.Name)]
        dupes = {k for k in keys if keys.count(k) > 1}
        assert not dupes, f"duplicate keys in _raw_event_from_sv base dict: {sorted(dupes)}"


def _run_mini_sv(tmp_path, **overrides):
    kwargs = dict(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions_exact.tsv",
        expressed_products_tsv=FX / "expressed_products.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    kwargs.update(overrides)
    return build_sv_phase1_raw(**kwargs)


def test_sv_event_labels_heuristic_components(tmp_path):
    out = _run_mini_sv(tmp_path)
    events = read_tsv(out["raw_events"])
    assert events
    status = events[0]["event_component_status"]
    assert "VAF_proxy" in status
    assert "fixed_prior" in status
    assert "not_assessed_pending_matched_population_data" in status
    # The clonality value stays the documented vaf/0.35 proxy.
    vaf = float(events[0]["tumor_vaf"])
    assert events[0]["clonality"] == f"{min(1.0, vaf / 0.35):.4f}"
    assert float(events[0]["persistence"]) == pytest.approx(0.70)


# ---------------------------------------------------------------------------
# B4-11: confirmed-product peptide exclusion ledger
# ---------------------------------------------------------------------------


def test_confirmed_product_without_junction_qc_is_ledged_not_silent(tmp_path):
    out = _run_mini_sv(tmp_path, rna_junction_tsv=None)
    rows = read_tsv(out["sv_peptide_exclusions"])
    assert len(rows) == 1
    row = rows[0]
    assert row["exclusion_reason"] == "rna_junction_qc_unavailable"
    assert row["junction_qc"] == "UNASSESSED"
    assert "Arriba" in row["product_provenance"]
    assert row["adjacency_key"]
    # The event stays evidence but produces no peptides.
    assert read_tsv(out["raw_peptides"]) == []
    prov = json.loads(Path(out["provenance"]).read_text(encoding="utf-8"))
    assert prov["confirmed_products_excluded_from_peptides"] == 1


def test_confirmed_product_with_qc_passing_junction_makes_peptides(tmp_path):
    out = _run_mini_sv(tmp_path)
    assert read_tsv(out["sv_peptide_exclusions"]) == []
    assert len(read_tsv(out["raw_peptides"])) > 0


# ---------------------------------------------------------------------------
# B4-12: shared predictor schema + missingness record
# ---------------------------------------------------------------------------


def _score_sv(tmp_path, *, stab=True, chop=True):
    adapter = tmp_path / "adapter"
    _run_mini_sv(adapter)
    peptides = read_tsv(adapter / "parsed" / "raw_peptides.tsv")
    pick = peptides[:3]
    # Binding evidence supplied for every candidate pair (provided files,
    # no stub runner) so the shared-schema join is exercised end to end.
    net_file = tmp_path / "net.xls"
    nlines = ["\t".join([
        "Pos", "Peptide", "HLA", "Score_EL", "%Rank_EL", "Score_BA", "%Rank_BA",
        "BindLevel", "AlleleSupport",
    ])]
    for p in peptides:
        nlines.append("\t".join(["1", p["peptide"], p["hla_allele"], "1.0", "1.0", "0.5", "2.0", "NB", "trained"]))
    net_file.write_text("\n".join(nlines) + "\n", encoding="utf-8")
    mhc_file = tmp_path / "mhc.csv"
    mlines = ["peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,mhcflurry_processing_score,mhcflurry_presentation_score"]
    mlines += [f"{p['peptide']},{p['hla_allele']},50.0,1.5,0.6,0.7" for p in peptides]
    mhc_file.write_text("\n".join(mlines) + "\n", encoding="utf-8")
    stab_file = tmp_path / "stab.tsv"
    rows = ["Peptide\tHLA\tscore\tpercentile_rank"]
    rows += [f"{p['peptide']}\t{p['hla_allele']}\t9.5\t1.5" for p in pick]
    stab_file.write_text("\n".join(rows) + "\n", encoding="utf-8")
    chop_file = tmp_path / "chop.tsv"
    crows = ["peptide\tnetchop_31d_max_score"]
    crows += [f"{p['peptide']}\t0.5" for p in pick]
    chop_file.write_text("\n".join(crows) + "\n", encoding="utf-8")
    out = run_sv_score(
        outdir=tmp_path / "scored",
        profile_name_or_path="sv_wgs_phase1",
        sample_id="SVMINI",
        raw_events=adapter / "parsed" / "raw_events.tsv",
        raw_peptides=adapter / "parsed" / "raw_peptides.tsv",
        netmhcpan=net_file,
        mhcflurry=mhc_file,
        netmhcstabpan=stab_file if stab else None,
        netchop=chop_file if chop else None,
        expression=FX / "expression.tsv",
        normal_expression=FX / "normal_expression.tsv",
        normal_hla_ligands=FX / "normal_hla_ligands.tsv",
        immunogenicity_stub=True,
    )
    return out, pick


def test_sv_score_records_predictor_missingness(tmp_path):
    out, _ = _score_sv(tmp_path, stab=True, chop=False)
    prov = json.loads(Path(out["provenance"]).read_text(encoding="utf-8"))
    status = prov["predictor_evidence_status"]
    assert status["netmhcstabpan"] == "PRESENT"
    assert status["netchop"] == "NOT_PROVIDED"
    assert status["netmhcpan"] == "PRESENT"
    assert status["mhcflurry"] == "PRESENT"


def test_sv_score_auxiliary_predictors_join_shared_schema(tmp_path):
    out, pick = _score_sv(tmp_path)
    pres = read_tsv(out["presentation_evidence"])
    wanted = {p["peptide"] for p in pick}
    matched = [r for r in pres if r["peptide"] in wanted]
    assert matched
    stab_rows = [r for r in matched if r.get("netmhcstabpan_rank") not in ("", None)]
    assert stab_rows, "provided NetMHCstabpan evidence did not join the SV peptides"
    assert all(float(r["netmhcstabpan_rank"]) == pytest.approx(1.5) for r in stab_rows)
    chop_rows = [r for r in matched if r.get("netchop_31d_max_score") not in ("", None)]
    assert chop_rows, "provided NetChop evidence did not join the SV peptides"
    # Provenance registry sees the passthrough predictor (not just the file).
    prov = json.loads(Path(out["provenance"]).read_text(encoding="utf-8"))
    tools = prov["tools"]
    names = {t.get("tool") for t in tools} if isinstance(tools, list) else set(tools)
    assert "netmhcstabpan" in names


# ---------------------------------------------------------------------------
# B4-13: adjacent-gene read-through projected match -> review flag
# ---------------------------------------------------------------------------

_GTF_ADJACENT = "\n".join([
    'chr1\tmini\tgene\t101\t750\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";',
    'chr1\tmini\texon\t101\t190\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1"; exon_number "1";',
    'chr1\tmini\texon\t301\t361\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1"; exon_number "2";',
    'chr1\tmini\texon\t601\t750\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1"; exon_number "3";',
    'chr1\tmini\tCDS\t101\t190\t.\t+\t0\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1";',
    'chr1\tmini\tCDS\t301\t361\t.\t+\t0\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1";',
    'chr1\tmini\tCDS\t601\t750\t.\t+\t2\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1";',
    'chr1\tmini\tgene\t851\t1500\t.\t+\t.\tgene_id "GENE2"; gene_name "GENE2";',
    'chr1\tmini\texon\t851\t950\t.\t+\t.\tgene_id "GENE2"; gene_name "GENE2"; transcript_id "TX2"; exon_number "1";',
    'chr1\tmini\texon\t1200\t1500\t.\t+\t.\tgene_id "GENE2"; gene_name "GENE2"; transcript_id "TX2"; exon_number "2";',
    'chr1\tmini\tCDS\t851\t950\t.\t+\t0\tgene_id "GENE2"; gene_name "GENE2"; transcript_id "TX2";',
    'chr1\tmini\tCDS\t1200\t1500\t.\t+\t0\tgene_id "GENE2"; gene_name "GENE2"; transcript_id "TX2";',
]) + "\n"


def test_readthrough_partner_detection(tmp_path):
    gtf = tmp_path / "adj.gtf"
    gtf.write_text(_GTF_ADJACENT, encoding="utf-8")
    ann = GtfAnnotation(gtf)
    assert ann.are_potential_readthrough_partners("chr1", "GENE1", "chr1", "GENE2")
    assert not ann.are_potential_readthrough_partners("chr1", "GENE1", "chr1", "GENE1")
    assert not ann.are_potential_readthrough_partners("chr1", "GENE1", "chr2", "GENE2")
    assert not ann.are_potential_readthrough_partners(
        "chr1", "GENE1", "chr1", "GENE2", max_gap_bp=10,
    )
    assert not ann.are_potential_readthrough_partners("chr1", "GENE1", "chr1", "MISSING")


def test_adjacent_gene_projected_match_is_review_flag(tmp_path):
    gtf = tmp_path / "adj.gtf"
    gtf.write_text(_GTF_ADJACENT, encoding="utf-8")
    ref = tmp_path / "mini_ref.fa"
    _write_fasta(ref, _seq(4000))
    vcf = tmp_path / "bnd_adj.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=SV>\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=End>\n"
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=CIPOS>\n"
        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=CIEND>\n"
        "##FORMAT=<ID=PR,Number=2,Type=Integer,Description=PR>\n"
        "##FORMAT=<ID=SR,Number=2,Type=Integer,Description=SR>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
        "chr1\t400\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1000;CIPOS=-5,5;CIEND=-5,5"
        "\tPR:SR:DP\t10,0:5,0:50\t14,8:10,6:60\n",
        encoding="utf-8",
    )
    junctions = tmp_path / "j.tsv"
    # Native read-through junction of GENE1->GENE2 (exon ends 361 -> 1200).
    write_tsv(junctions, [{
        "genome_build": "GRCh38", "chrom1": "chr1", "pos1": "361", "strand1": "+",
        "chrom2": "chr1", "pos2": "1200", "strand2": "+",
        "split_reads": "10", "spanning_reads": "6", "unique_start_count": "4",
        "min_anchor_bp": "30", "min_mapq": "60", "source_tool": "STAR",
    }], ["genome_build", "chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2",
         "split_reads", "spanning_reads", "unique_start_count", "min_anchor_bp",
         "min_mapq", "source_tool"])
    out = build_sv_phase1_raw(
        sample_id="SVMINI", sv_vcfs=[vcf], callers=["GRIDSS2"],
        tumor_sample_name="TUMOR", normal_sample_name="NORMAL",
        reference_fasta=ref, gencode_gtf=gtf, hla=FX / "hla.txt", outdir=tmp_path / "a",
        rna_junction_tsv=junctions,
        expression_tsv=FX / "expression.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    events = read_tsv(out["sv_events_full"])
    assert len(events) == 1
    assert events[0]["rna_evidence_match"] == "PROJECTED_PARTNER_JUNCTION_ADJACENT_GENES_REVIEW"
    assert events[0]["rna_support_status"] == "RNA_JUNCTION_REVIEW_ONLY"
    assert events[0]["rna_junction_reads"] == "0"
    # The review flag keeps the DNA hypothesis from being promoted to peptides.
    assert read_tsv(out["raw_peptides"]) == []


# ---------------------------------------------------------------------------
# B4-15: joint capture tier
# ---------------------------------------------------------------------------


def _regions(tmp_path, intervals: list[tuple[str, int, int]]) -> CaptureRegions:
    bed = tmp_path / "kit.bed"
    bed.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals), encoding="utf-8")
    return CaptureRegions.from_bed(bed)


def _cluster(bp1: int, bp2: int) -> SVCluster:
    rec = SVRecord("r1", "Manta", "chr1", bp1, "chr1", bp2, "BND")
    return SVCluster(records=[rec])


def test_capture_joint_tier(tmp_path):
    regions = _regions(tmp_path, [("chr1", 1000, 2000)])
    joint = regions.joint_interpretability
    assert joint("ON_TARGET", "ON_TARGET") == "HIGH"
    assert joint("ON_TARGET", "NEAR_TARGET_250BP") == "HIGH"
    assert joint("ON_TARGET", "NEAR_TARGET_1000BP") == "MEDIUM"   # one-sided partial
    assert joint("ON_TARGET", "OFF_TARGET") == "LOW"              # one-sided only
    assert joint("ON_TARGET", "NOT_ASSESSED") == "NOT_ASSESSED"
    assert joint("NEAR_TARGET_1000BP", "NEAR_TARGET_1000BP") == "MEDIUM"
    assert joint("OFF_TARGET", "OFF_TARGET") == "UNINTERPRETABLE"
    # End-to-end through classify_cluster: a single on-target breakend is no
    # longer reported as HIGH interpretability.
    cls = regions.classify_cluster(_cluster(1500, 9000))
    assert cls["breakend1_capture_status"] == "ON_TARGET"
    assert cls["breakend2_capture_status"] == "OFF_TARGET"
    assert cls["capture_interpretability"] == "LOW"
    both = regions.classify_cluster(_cluster(1500, 1600))
    assert both["capture_interpretability"] == "HIGH"
    ann = annotate_cluster_capture(_cluster(1500, 9000), regions)
    assert ann["capture_filter_status"] == "REVIEW"


# ---------------------------------------------------------------------------
# B4-17: capture annotation is a declared typed field
# ---------------------------------------------------------------------------


def test_capture_annotation_survives_object_copy(tmp_path):
    regions = _regions(tmp_path, [("chr1", 1000, 2000)])
    cluster = _cluster(1500, 1600)
    ann = annotate_cluster_capture(cluster, regions)
    assert "capture_interpretability" in ann
    # Declared field, not a setattr bolt-on.
    assert "capture_annotation" in {f.name for f in dataclasses.fields(SVCluster)}
    # Survives a constructor-based copy (the old setattr path lost this).
    clone = dataclasses.replace(cluster, records=list(cluster.records))
    assert cluster_capture_annotation(clone)["capture_interpretability"] == ann["capture_interpretability"]
    # Accessor always yields a dict, even for an unannotated cluster.
    assert cluster_capture_annotation(SVCluster(records=[])) == {}


# ---------------------------------------------------------------------------
# B4-18: exact junction keys vs exploratory gene-pair keys
# ---------------------------------------------------------------------------


def test_load_junction_reads_splits_exact_from_exploratory(tmp_path):
    tsv = tmp_path / "j.tsv"
    write_tsv(tsv, [
        {"chrom1": "chr1", "pos1": "1000", "chrom2": "chr2", "pos2": "2000",
         "split_reads": "8", "gene1": "ALK", "gene2": "EML4"},
        {"chrom1": "chr1", "pos1": "5000", "chrom2": "chr2", "pos2": "9000",
         "split_reads": "3", "gene1": "ALK", "gene2": "EML4"},
        {"event_id": "EVT_9", "junction_reads": "5"},
    ], ["chrom1", "pos1", "chrom2", "pos2", "split_reads", "gene1", "gene2",
        "event_id", "junction_reads"])
    data = load_junction_reads(tsv)
    assert isinstance(data, JunctionReadEvidence)
    assert data.exact[junction_positional_key("chr1", 1000, "chr2", 2000)] == 8
    assert data.exact[junction_positional_key("chr1", 5000, "chr2", 9000)] == 3
    assert data.exact["EVT_9"] == 5
    # Gene-pair aggregation exists but lives only in the exploratory bucket.
    assert data.gene_pair["ALK::EML4"] == 8
    assert "ALK::EML4" not in data.exact and "EML4::ALK" not in data.exact


def test_rna_junction_evidence_no_gene_pair_cross_attribution(tmp_path):
    from neoag.evidence_layer import build_rna_junction_evidence

    junctions = tmp_path / "junctions.tsv"
    write_tsv(junctions, [
        {"chrom1": "chr1", "pos1": "1000", "chrom2": "chr2", "pos2": "2000",
         "split_reads": "8", "gene1": "ALK", "gene2": "EML4"},
        {"chrom1": "chr1", "pos1": "5000", "chrom2": "chr2", "pos2": "9000",
         "split_reads": "3", "gene1": "ALK", "gene2": "EML4"},
    ], ["chrom1", "pos1", "chrom2", "pos2", "split_reads", "gene1", "gene2"])
    events = tmp_path / "events.tsv"
    write_tsv(events, [
        {"event_id": "EV_MATCH", "gene": "ALK::EML4", "chrom1": "chr1", "pos1": "1000",
         "chrom2": "chr2", "pos2": "2000", "mutation_source": "SV"},
        {"event_id": "EV_OTHER_BP", "gene": "ALK::EML4", "chrom1": "chr1", "pos1": "7777",
         "chrom2": "chr2", "pos2": "8888", "mutation_source": "SV"},
    ], ["event_id", "gene", "chrom1", "pos1", "chrom2", "pos2", "mutation_source"])
    peptides = tmp_path / "peptides.tsv"
    peptides.write_text("event_id\n", encoding="utf-8")
    rows = build_rna_junction_evidence(
        events, peptides, tmp_path / "out.tsv", junction_path=junctions,
    )
    by_id = {r["event_id"]: r for r in rows}
    # Exact-coordinate event keeps full support.
    assert by_id["EV_MATCH"]["junction_reads"] == "8"
    # Same gene pair, different breakpoints: the 8-read junction must NOT leak.
    assert by_id["EV_OTHER_BP"]["junction_reads"] == "0"
    # ...but the exploratory column still reports the pair-level support.
    assert by_id["EV_OTHER_BP"]["rna_junction_gene_pair_reads"] == "8"
