"""B-E review package E (2026-11-06): splice chain fixes B3b-01..13.

Each test mirrors one reviewed finding and demonstrates the behaviour change
directly at the adapter/prefilter level (plus two layer-level integrations).
"""

from pathlib import Path

import pytest

from neoag.splice.adapters.easyquant import parse_easyquant
from neoag.splice.adapters.irfinder import parse_irfinder
from neoag.splice.adapters.splice2neo import parse_splice2neo
from neoag.splice.adapters.spladder import parse_spladder_txt
from neoag.splice.adapters.immunopepper import parse_immunopepper_kmers, parse_immunopepper_meta
from neoag.splice.causal_status import CAUSAL_STATUS_PRIORITY
from neoag.splice.consensus import _favorable_presentation, build_consensus
from neoag.splice.coordinates import (
    JunctionNormalizationError,
    check_junction_format_vs_system,
    peptide_metadata,
    sniff_junction_table_format,
)
from neoag.splice.evidence_chains import _CAUSAL_PRIORITY as CHAINS_PRIORITY
from neoag.splice.pipeline import (
    _STATS_SELECTION_RULE,
    _STATUS_PRIORITY,
    SpliceLayer,
    _best,
    _numeric_pick,
    _stats_donor,
    build_splice_provenance_layer,
)
from neoag.splice_prefilter import prefilter_splice_peptides
from neoag.utils import read_tsv, write_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------- B3b-01 ----

def test_spladder_headerless_rejects_appended_integer_count_columns(tmp_path):
    # SplAdder appends integer sample-count columns after the exon run; the
    # bogus trailing "exon pair" (42-8) must be refused, not mis-parsed.
    txt = _write(
        tmp_path / "merge_graphs_exon_skip_C3.confirmed.txt",
        "chr3\t+\tEV42\tGENE3\t100\t150\t201\t250\t42\t8\n",
    )
    with pytest.raises(ValueError, match="review B3b-01"):
        parse_spladder_txt(txt, sample_id="S1")


def test_spladder_headerless_rejects_odd_coordinate_cells(tmp_path):
    txt = _write(
        tmp_path / "merge_graphs_exon_skip_C5.confirmed.txt",
        "chr3\t+\tEV42\tGENE3\t100\t150\t201\t250\t42\n",
    )
    with pytest.raises(ValueError, match="odd coordinate cell count"):
        parse_spladder_txt(txt, sample_id="S1")


def test_spladder_headerless_accepts_valid_event_row(tmp_path):
    txt = _write(
        tmp_path / "merge_graphs_exon_skip_C7.confirmed.txt",
        "chr3\t+\tEV42\tGENE3\t100\t150\t201\t250\tS1:1\tS1:8\n",
    )
    bundle = parse_spladder_txt(txt, sample_id="S1")
    assert bundle["events"][0]["junction_ids"] == "SJ|GRCh38|chr3|151|200|+"


# ---------------------------------------------------------------- B3b-02 ----

def test_best_numeric_conflict_picks_strongest_not_lexicographic():
    # Old behaviour: sorted()[0] -> "1.5" / "4.9" (string order).
    assert float(_best({"1.5", "2.0"})) == 2.0
    assert float(_best({"5.0", "4.9"})) == 5.0


def test_best_keeps_lexicographic_pick_for_text_values():
    assert _best({"BETA", "ALPHA"}) == "ALPHA"


def test_stats_donor_prefers_smallest_qvalue_then_strongest_effect():
    row_a = {"delta_psi": "-0.9", "qvalue": "0.4", "psi": "0.5", "source_record_ids": "R1"}
    row_b = {"delta_psi": "0.1", "qvalue": "0.00001", "psi": "0.6", "source_record_ids": "R2"}
    assert _stats_donor([row_a, row_b]) is row_b
    # Equal qvalues: the larger |delta_psi| row donates its whole tuple.
    row_d = {"delta_psi": "0.1", "qvalue": "0.4", "psi": "0.6", "source_record_ids": "R2"}
    assert _stats_donor([row_d, row_a]) is row_a
    # A row reporting no statistics at all never donates a tuple.
    row_e = {"qvalue": "", "delta_psi": "", "psi": "", "source_record_ids": "R3"}
    assert _stats_donor([row_e, row_a]) is row_a
    assert _stats_donor([row_e]) is None


def test_same_event_statistics_stay_a_coherent_tuple(tmp_path):
    # Operator test (round 2, review B3b-02): two splice2neo rows map to the
    # SAME event but carry different (delta_psi, qvalue) pairs.  The merged
    # row must preserve one row's complete pair - never the per-column
    # winners - and must name the exact source record it came from.
    txt = _write(
        tmp_path / "splice2neo.tsv",
        "junc_id\tevent_type\tgene\tpsi\tdelta_psi\tqvalue\n"
        "chr1:151-200:+\tSE\tGENE1\t0.5\t-0.9\t0.4\n"
        "chr1:151-200:+\tSE\tGENE1\t0.6\t0.1\t0.00001\n",
    )
    bundle = parse_splice2neo(txt, sample_id="S1")
    assert len(bundle["events"]) == 2
    assert bundle["events"][0]["splice_event_id"] == bundle["events"][1]["splice_event_id"]

    layer = SpliceLayer(sample_id="S1")
    layer.extend(bundle)
    layer.consolidate()
    events = layer.tables["events"]
    assert len(events) == 1
    merged = events[0]
    source = bundle["events"][1]
    # The statistically stronger pair belongs to row 2: q=0.00001 reported
    # FOR delta_psi=0.1.  The old per-column merge fabricated (-0.9, 0.00001)
    # - row 1's effect size paired with row 2's significance, a combination
    # no test ever produced.
    assert merged["delta_psi"] == source["delta_psi"]
    assert merged["qvalue"] == source["qvalue"]
    assert merged["psi"] == source["psi"]
    assert merged["stats_selection_rule"] == _STATS_SELECTION_RULE
    assert merged["stats_source_record_ids"] == source["source_record_ids"]
    avoided = [
        c for c in layer.tables["conflicts"] if c["conflict_type"] == "STATS_TUPLE_CROSS_PAIRING_AVOIDED"
    ]
    assert len(avoided) == 1
    assert avoided[0]["severity"] == "WARNING"
    assert avoided[0]["resolution_status"] == "COHERENT_TUPLE_SELECTED"
    assert "delta_psi" in avoided[0]["resolution_reason"]


def test_identical_stat_rows_merge_without_conflict(tmp_path):
    txt = _write(
        tmp_path / "splice2neo.tsv",
        "junc_id\tevent_type\tgene\tpsi\tdelta_psi\tqvalue\n"
        "chr1:151-200:+\tSE\tGENE1\t0.5\t-0.9\t0.4\n"
        "chr1:151-200:+\tSE\tGENE1\t0.5\t-0.9\t0.4\n",
    )
    layer = SpliceLayer(sample_id="S1")
    layer.extend(parse_splice2neo(txt, sample_id="S1"))
    layer.consolidate()
    merged = layer.tables["events"][0]
    assert merged["stats_selection_rule"] == "" or "stats_selection_rule" not in merged
    assert not any(
        c["conflict_type"] == "STATS_TUPLE_CROSS_PAIRING_AVOIDED" for c in layer.tables["conflicts"]
    )


def test_typed_aggregation_helpers():
    assert float(_numeric_pick({"0.0001", "0.002", "0.00002"}, mode="min")) == 0.00002
    assert float(_numeric_pick({"-0.9", "0.3"}, mode="maxabs")) == -0.9
    assert float(_numeric_pick({"7", "12"}, mode="max")) == 12.0


# ---------------------------------------------------------------- B3b-03 ----

_REGTOOLS = (
    "chrom\tstart\tend\tname\tscore\tstrand\n"
    "chr1\t150\t200\tJ1\t17\t+\n"
)
_SPLADDER_OTHER = "chr3\t+\tEV9\tGENE9\t100\t150\t201\t250\tS1:1\n"
_META_ON_REGTOOLS_JUNCTION = (
    "peptide\tid\treadFrame\treadFrameAnnotated\tgeneName\tgeneChr\tgeneStrand\t"
    "mutationMode\thasStopCodon\tisJunctionList\tvariantSegExpr\tmodifiedExonsCoord\n"
    "MARNDCEQGHILKFPSTWYV\tTX1\t0\tTrue\tGENE1\tchr1\t+\tsplice\tFalse\tTrue\t12.5\t"
    "chr1:100-150;chr1:201-250\n"
)


def test_immunopepper_join_is_caller_independent(tmp_path):
    # A RegTools-only junction (no SplAdder coverage) used to be outside the
    # ImmunoPepper target set, silently dropping its translation metadata.
    junctions = _write(tmp_path / "junctions.bed", _REGTOOLS)
    spladder = _write(tmp_path / "merge_graphs_exon_skip_A.confirmed.txt", _SPLADDER_OTHER)
    meta = _write(tmp_path / "meta.tsv", _META_ON_REGTOOLS_JUNCTION)
    build_splice_provenance_layer(
        sample_id="S1",
        outdir=tmp_path / "out",
        junctions=junctions,
        spladder_txt=[spladder],
        immunopepper_meta=[meta],
    )
    transcripts = read_tsv(tmp_path / "out" / "splice_transcript_hypotheses.tsv")
    chains = {row.get("junction_chain", "") for row in transcripts}
    assert any("SJ|GRCh38|chr1|151|200|+" in chain for chain in chains)


def test_kmer_orphans_get_one_explicit_summary(tmp_path):
    meta = _write(tmp_path / "meta.tsv", _META_ON_REGTOOLS_JUNCTION)
    bundle = parse_immunopepper_meta(meta, sample_id="S1")
    kmers = _write(
        tmp_path / "kmers.tsv",
        "kmer\tcoord\tisCrossJunction\tjunctionAnnotated\n"
        "QQQQQQQQ\t1:2:3:4\tTrue\tTrue\n",
    )
    kb = parse_immunopepper_kmers(
        kmers, sample_id="S1", meta_bundle={"orfs": bundle["orfs"], "transcripts": bundle["transcripts"]},
        record_unmapped_conflicts=False,
    )
    assert [c["conflict_type"] for c in kb["conflicts"]] == ["IMMUNOPEPPER_KMER_ORPHAN_SUMMARY"]
    assert kb["manifest"][0]["unmapped_conflict_rows"] == "suppressed_summary_only"


# ---------------------------------------------------------------- B3b-04 ----

_SJ_OUT_TAB = "chr1\t2000\t3000\t1\tGTAGC\t500\t120\t380\t400\n"
_STAR_BED = "chr1\t1999\t3000\tchr1:2000-3000/300-300\t255\t+\n"


def test_sniff_distinguishes_bed_from_sj_out_tab(tmp_path):
    assert sniff_junction_table_format(_write(tmp_path / "SJ.out.tab", _SJ_OUT_TAB)) == "star_sj"
    assert sniff_junction_table_format(_write(tmp_path / "junctions.bed", _STAR_BED)) == "bed"


def test_bed_input_with_star_identity_system_is_refused(tmp_path):
    bed = _write(tmp_path / "junctions.bed", _STAR_BED)
    with pytest.raises(JunctionNormalizationError, match="review B3b-04"):
        check_junction_format_vs_system(bed, source_tool="STAR", coordinate_system="auto")
    # The same BED is accepted under its real coordinate system.
    check_junction_format_vs_system(bed, source_tool="STAR", coordinate_system="bed0_half_open")


def test_sj_input_with_bed0_system_is_refused(tmp_path):
    sj = _write(tmp_path / "SJ.out.tab", _SJ_OUT_TAB)
    with pytest.raises(JunctionNormalizationError, match="review B3b-04"):
        check_junction_format_vs_system(sj, source_tool="STAR", coordinate_system="bed0_half_open")
    check_junction_format_vs_system(sj, source_tool="STAR", coordinate_system="star_sj")


# ---------------------------------------------------------------- B3b-05 ----

def _easyquant_inputs(tmp_path, junc: str, span: str):
    meta = _write(tmp_path / "meta.tsv", _META_ON_REGTOOLS_JUNCTION)
    bundle = parse_immunopepper_meta(meta, sample_id="S1")
    event_id = bundle["events"][0]["splice_event_id"]
    junction_id = bundle["events"][0]["junction_ids"]
    query_map = _write(
        tmp_path / "query_map.tsv",
        "query_id\tquery_name\tsplice_event_id\tjunction_id\tvariant_id\ttranscript_hypothesis_id\tposition_1based\n"
        f"Q1\tQ1\t{event_id}\t{junction_id}\tVAR|GRCh38|chr1|1000|A|G\tT1\t3\n",
    )
    quant = _write(
        tmp_path / "easyquant.tsv",
        "name\tjunc\tspan\tanch\ta\tb\tinterval\twithin_interval\tcoverage_perc\tcoverage_mean\tcoverage_median\n"
        f"Q1\t{junc}\t{span}\t4\t2\t2\t1-100\t50\t0.9\t20\t18\n",
    )
    return quant, query_map


def test_negative_targeted_requant_stays_prediction_only(tmp_path):
    quant, qmap = _easyquant_inputs(tmp_path, "0", "0")
    bundle = parse_easyquant(quant, sample_id="S1", query_map=qmap)
    link = bundle["causal_links"][0]
    assert link["causal_status"] == "DNA_PREDICTION_ONLY"
    assert link["targeted_requant_status"] == "TARGETED_REQUANT_NEGATIVE"


def test_spanning_only_requant_is_not_positive_cis_support(tmp_path):
    quant, qmap = _easyquant_inputs(tmp_path, "0", "5")
    bundle = parse_easyquant(quant, sample_id="S1", query_map=qmap)
    link = bundle["causal_links"][0]
    # Spanning pairs without junction-specific reads: honest spanning status,
    # but no positive DNA->RNA cis claim (review B3b-05).
    assert link["causal_status"] == "DNA_PREDICTION_ONLY"
    assert link["targeted_requant_status"] == "TARGETED_SPANNING_ONLY"


def test_positive_targeted_requant_still_supported(tmp_path):
    quant, qmap = _easyquant_inputs(tmp_path, "9", "5")
    bundle = parse_easyquant(quant, sample_id="S1", query_map=qmap)
    link = bundle["causal_links"][0]
    assert link["causal_status"] == "TARGETED_REQUANT_SUPPORTED"
    assert link["rna_junction_status"] == "EXACT_RNA_SUPPORTED"


# ---------------------------------------------------------------- B3b-06/07 ----

def _prefilter_rows(tmp_path, extra: dict):
    base = {
        "peptide_id": "P1", "event_id": "EV1", "sample_id": "S1",
        "event_type": "splice", "peptide_consequence": "splice_junction",
        "gene": "G1", "peptide": "ACDEFGHIK", "hla_allele": "HLA-A*02:01",
        "mhc_class": "I",
    }
    base.update({k: v for k, v in extra.items() if v is not None})
    peptides = tmp_path / "peptides.tsv"
    write_tsv(peptides, [base])
    return peptides


_PROFILE_REVIEW = {"gates": {}}
_PROFILE_REJECT = {"gates": {"junction_read_qc_incomplete_policy": "reject"}}


def _prefilter_status(tmp_path, extra, profile):
    peptides = _prefilter_rows(tmp_path, extra)
    return prefilter_splice_peptides(peptides, profile, tmp_path / "pf")


def test_incomplete_junction_qc_defaults_to_review_not_reject(tmp_path):
    result = _prefilter_status(tmp_path, {"junction_read_qc_status": "INCOMPLETE:MAPQ_MISSING"}, _PROFILE_REVIEW)
    row = next(r for r in result["decisions"]) if isinstance(result.get("decisions"), list) else None
    statuses = read_tsv(tmp_path / "pf" / "splice_prefilter_decisions.tsv")
    assert statuses[0]["prefilter_status"] == "REVIEW"
    assert "FORMAL_JUNCTION_READ_QC" in statuses[0]["unassessed_stages"]
    assert "FORMAL_JUNCTION_READ_QC" not in statuses[0]["failed_stages"]


def test_incomplete_junction_qc_rejects_only_under_explicit_policy(tmp_path):
    _prefilter_status(tmp_path, {"junction_read_qc_status": "INCOMPLETE:MAPQ_MISSING"}, _PROFILE_REJECT)
    statuses = read_tsv(tmp_path / "pf" / "splice_prefilter_decisions.tsv")
    assert statuses[0]["prefilter_status"] == "REJECT"
    assert "FORMAL_JUNCTION_READ_QC" in statuses[0]["failed_stages"]


def test_incomplete_junction_qc_never_reads_as_verified(tmp_path):
    # Review B3b-06 round-2 guard (HANDOFF acceptance): the REVIEW branch is
    # reasonable, but no downstream consumer may misread INCOMPLETE QC as
    # verified junction evidence.
    from neoag.splice_prefilter import _junction_read_qc_decision

    for value in ("INCOMPLETE", "INCOMPLETE:MAPQ_MISSING", "QC_INCOMPLETE_OVERHANG", "incomplete"):
        for policy in ("review", "reject", "anything_else"):
            verdict, _reason = _junction_read_qc_decision(
                [{"junction_read_qc_status": value}], incomplete_policy=policy
            )
            assert verdict is not True, (value, policy)

    _prefilter_status(tmp_path, {"junction_read_qc_status": "INCOMPLETE:MAPQ_MISSING"}, _PROFILE_REVIEW)
    statuses = read_tsv(tmp_path / "pf" / "splice_prefilter_decisions.tsv")
    row = statuses[0]
    assert row["prefilter_status"] == "REVIEW"
    # The pool is explicitly labelled exploratory, and the formal gate says
    # no - nothing here claims a high-confidence source chain.
    assert row["splice_candidate_pool"] == "EXPLORATION_EVIDENCE_INCOMPLETE"
    assert row["splice_formal_gate_pass"] == "no"
    # The stage detail keeps the raw status AND the insufficient marker, so
    # later steps see why the QC is unverified, never a bare "PASS".
    assert "INCOMPLETE:MAPQ_MISSING" in row["stage_details"]
    assert "insufficient_junction_read_qc_review" in row["stage_details"]


def test_o3_orf_grade_passes_credible_orf(tmp_path):
    _prefilter_status(
        tmp_path,
        {"junction_read_qc_status": "PASS", "orf_evidence_grade": "O3"},
        _PROFILE_REVIEW,
    )
    statuses = read_tsv(tmp_path / "pf" / "splice_prefilter_decisions.tsv")
    assert "CREDIBLE_ORF" not in statuses[0]["failed_stages"]
    assert "CREDIBLE_ORF" not in statuses[0]["unassessed_stages"]


def test_nmd_escape_status_passes_nmd_stage(tmp_path):
    # "NMD_ESCAPE" contains the FAIL token "NMD"; escape precedence must win.
    _prefilter_status(
        tmp_path,
        {"junction_read_qc_status": "PASS", "splice_nmd_status": "NMD_ESCAPE"},
        _PROFILE_REVIEW,
    )
    statuses = read_tsv(tmp_path / "pf" / "splice_prefilter_decisions.tsv")
    assert "NMD" not in statuses[0]["failed_stages"]
    assert "NMD" not in statuses[0]["unassessed_stages"]


def test_predicted_nmd_still_fails(tmp_path):
    _prefilter_status(
        tmp_path,
        {"junction_read_qc_status": "PASS", "splice_nmd_status": "PREDICTED_NMD"},
        _PROFILE_REVIEW,
    )
    statuses = read_tsv(tmp_path / "pf" / "splice_prefilter_decisions.tsv")
    assert "NMD" in statuses[0]["failed_stages"]


def test_top_pairs_truncation_follows_binder_rank_not_file_order(tmp_path):
    # 25 peptide/HLA rows for one event; only max_pairs=20 are kept.  The
    # kept set must be the best-ranked 20 regardless of input order.
    rows = []
    for i in range(25):
        rows.append({
            "peptide_id": f"P{i:02d}", "event_id": "EV1", "sample_id": "S1",
            "event_type": "splice", "peptide_consequence": "splice_junction",
            "gene": "G1", "peptide": "ACDEFGHIK", "hla_allele": "HLA-A*02:01",
            "mhc_class": "I", "el_rank": str(i),
        })
    peptides = tmp_path / "peptides.tsv"
    write_tsv(peptides, list(reversed(rows)))
    prefilter_splice_peptides(peptides, {"gates": {"max_splice_peptide_hla_per_event": 20}}, tmp_path / "pf")
    kept = {row["peptide_id"] for row in read_tsv(peptides) if row.get("event_id") == "EV1"}
    assert kept == {f"P{i:02d}" for i in range(20)}


# ---------------------------------------------------------------- B3b-09 ----

def test_missing_crosses_junction_is_unassessed_not_presumed(tmp_path):
    class _Rec:
        row: dict
    rec = _Rec()
    rec.row = {"peptide": "ACDEFGHIK", "hla_allele": "HLA-A*02:01"}
    meta = peptide_metadata(rec)
    assert meta["crosses_junction"] == "UNASSESSED"


# ---------------------------------------------------------------- B3b-10 ----

def _write_sj_tab(path: Path, chrom: str, start: int, end: int, reads: str) -> Path:
    return _write(
        path,
        f"{chrom}\t{start}\t{end}\t1\tGTAGC\t{reads}\t{reads}\t5\t400\n",
    )


def test_verified_reads_follow_actual_input_source(tmp_path):
    from neoag.splice.registry import build_support_index

    sj = _write_sj_tab(tmp_path / "SJ.out.tab", "chr1", 2000, 3000, "40")
    index = build_support_index(sj, sample_id="S1", source_tool="STAR-SJ")
    assert index.primary_tools == {"STAR-SJ"}
    rows = index.registry.entity_rows(sample_id="S1", primary_tools=set(index.primary_tools) or {"RegTools"})
    # Entity rows expose the tool-filtered max as "rna_junction_reads"
    # (registry.as_dict); unique=40, multi=5 -> record total 45.
    assert int(rows[0]["rna_junction_reads"]) == 45
    # The old hardcoded {"RegTools"} primary would have zeroed a STAR-only run.
    zeroed = index.registry.entity_rows(sample_id="S1", primary_tools={"RegTools"})
    assert int(zeroed[0]["rna_junction_reads"]) == 0


# ---------------------------------------------------------------- B3b-11 ----

def test_ir_site_coverage_is_not_split_reads(tmp_path):
    ir = _write(
        tmp_path / "ir.tsv",
        "Chr\tStart\tEnd\tName\tStrand\tIntronDepth\tSpliceLeft\tSpliceRight\tSpliceExact\tIRratio\tWarnings\tCNN\n"
        "chr1\t1000\t2000\tGENE1\t+\t35\t2\t1\t2\t0.42\tok\tpass\n",
    )
    bundle = parse_irfinder(ir, sample_id="S1")
    junction = bundle["junctions"][0]
    # Low-coverage intron retention must not satisfy the split-read gate.
    assert junction["unique_split_reads"] == ""
    assert junction["total_split_reads"] == ""
    assert junction["ir_site_coverage"] == "2"
    assert junction["ir_site_coverage_source"] == "IRFinder-S:SpliceExact"
    assert junction["split_read_support_status"] == "NOT_MEASURED_BY_IRFINDER"


# ---------------------------------------------------------------- B3b-12 ----

def _consensus_tables(percentile: str):
    return {
        "peptide_origins": [{
            "origin_peptide_id": "POR1", "peptide_id": "PEP1", "orf_id": "ORF1",
            "splice_event_id": "EV1", "sample_id": "S1", "gene": "G1",
            "peptide_sequence": "ACDEFGHIK", "crosses_junction": "true",
            "contains_novel_aa": "true",
        }],
        "peptide_origin_links": [],
        "presentation": [{
            "origin_peptide_id": "POR1",
            "best_presentation_percentile": percentile,
            "median_presentation_percentile": percentile,
        }],
    }


def test_weak_presentation_row_no_longer_earns_favorable_grade():
    rows = build_consensus(_consensus_tables("40.0"), sample_id="S1")
    row = rows[0]
    assert row["presentation_grade"] == "P1"
    assert "CAP_PRESENTATION_PREDICTED_NOT_FAVORABLE_R3" in row["cap_codes"]
    # This minimal fixture carries no ORF table, so HARD_ORF_INVALID pins the
    # tier at R4 before caps apply; the B3b-12 signal is the P1 grade + cap
    # (the old code graded percentile 40.0 as favorable P2 and earned R1/R2).
    assert row["final_evidence_tier"] == "R4"


def test_favorable_presentation_passes_the_gate():
    rows = build_consensus(_consensus_tables("1.2"), sample_id="S1")
    assert rows[0]["presentation_grade"] == "P2"
    assert "CAP_PRESENTATION_PREDICTED_NOT_FAVORABLE_R3" not in rows[0]["cap_codes"]


def test_favorable_helper_uses_best_value():
    assert _favorable_presentation([{"best_presentation_percentile": "1.5"}]) is True
    assert _favorable_presentation([{"best_presentation_percentile": "3.0", "median_presentation_percentile": "0.1"}]) is False
    assert _favorable_presentation([{"best_presentation_percentile": "not_a_number"}]) is False


# ---------------------------------------------------------------- B3b-13 ----

def test_one_shared_causal_ladder_everywhere():
    assert CHAINS_PRIORITY is CAUSAL_STATUS_PRIORITY
    assert _STATUS_PRIORITY["causal_status"] is CAUSAL_STATUS_PRIORITY
    # The old evidence_chains ladder tied these two; there must be no tie.
    assert (
        CAUSAL_STATUS_PRIORITY["PVACSPLICE_SUPPORTED"]
        > CAUSAL_STATUS_PRIORITY["TARGETED_SPANNING_ONLY"]
        > CAUSAL_STATUS_PRIORITY["DNA_RNA_CIS_SUPPORTED"]
    )


# ------------------------- B3b-02 round-3 follow-up (FURTHER_SUGGESTIONS WS1) --

def _missing_field_layer(tmp_path, swapped: bool):
    rows = [
        "chr1:151-200:+\tSE\tGENE1\t0.5\t0.4\t\n",
        "chr1:151-200:+\tSE\tGENE1\t0.6\t\t0.001\n",
    ]
    if swapped:
        rows.reverse()
    txt = _write(
        tmp_path / f"s2n_missing_{int(swapped)}.tsv",
        "junc_id\tevent_type\tgene\tpsi\tdelta_psi\tqvalue\n" + "".join(rows),
    )
    bundle = parse_splice2neo(txt, sample_id="S1")
    layer = SpliceLayer(sample_id="S1")
    layer.extend(bundle)
    layer.consolidate()
    return bundle, layer


def test_b3b_02_round3_missing_donor_field_is_never_backfilled(tmp_path):
    """Operator suggestion 2 (2026-11-08): row A reports delta_psi=0.4 with
    NO qvalue; row B reports qvalue=0.001 with NO delta_psi.  The merged
    primary tuple must not become (0.4, 0.001) - once a donor is selected
    its tuple stays intact, missing values included, and the complementary
    measurement survives only as a named disclosure with provenance."""
    for swapped in (False, True):  # donor choice invariant under permutation
        bundle, layer = _missing_field_layer(tmp_path, swapped)
        assert len(layer.tables["events"]) == 1
        merged = layer.tables["events"][0]
        # Row B owns the tuple (smallest reported qvalue wins).
        donor = [r for r in bundle["events"] if str(r.get("qvalue", "")).strip() == "0.001"][0]
        assert merged["qvalue"] == "0.001"
        # Row B never reported delta_psi, so the merged row carries NO effect
        # size - never row A's 0.4 dressed up as the donor's own measurement.
        assert merged["delta_psi"] == ""
        assert (merged["delta_psi"], merged["qvalue"]) != ("0.4", "0.001")
        # Missingness and donor identity are preserved in the output...
        assert merged["stats_missing_fields"] == "delta_psi"
        assert merged["stats_source_record_ids"] == donor["source_record_ids"]
        assert merged["stats_selection_rule"] == _STATS_SELECTION_RULE
        # ...and the complementary measurement is preserved separately.
        not_backfilled = [
            c for c in layer.tables["conflicts"]
            if c["conflict_type"] == "STATS_TUPLE_MISSING_FIELD_NOT_BACKFILLED"
        ]
        assert len(not_backfilled) == 1
        entry = not_backfilled[0]
        assert entry["field_name"] == "delta_psi"
        assert entry["observed_values"] == "0.4"
        assert entry["severity"] == "WARNING"
        assert entry["resolution_status"] == "DONOR_TUPLE_KEPT_INTACT"
        reporting = [r for r in bundle["events"] if str(r.get("delta_psi", "")).strip()][0]
        assert reporting["source_record_ids"] in entry["source_record_ids"]


def test_b3b_02_round3_complete_donor_tuple_has_no_missing_disclosure(tmp_path):
    # Both statistics present in the donor row: the tuple is intact by
    # construction, so no missing-field disclosure or conflict is emitted.
    txt = _write(
        tmp_path / "s2n_complete.tsv",
        "junc_id\tevent_type\tgene\tpsi\tdelta_psi\tqvalue\n"
        "chr1:151-200:+\tSE\tGENE1\t0.5\t-0.9\t0.4\n"
        "chr1:151-200:+\tSE\tGENE1\t0.6\t0.1\t0.00001\n",
    )
    layer = SpliceLayer(sample_id="S1")
    layer.extend(parse_splice2neo(txt, sample_id="S1"))
    layer.consolidate()
    merged = layer.tables["events"][0]
    assert float(merged["delta_psi"]) == 0.1 and float(merged["qvalue"]) == 0.00001
    assert merged.get("stats_missing_fields", "") == ""
    assert not any(
        c["conflict_type"] == "STATS_TUPLE_MISSING_FIELD_NOT_BACKFILLED"
        for c in layer.tables["conflicts"]
    )
