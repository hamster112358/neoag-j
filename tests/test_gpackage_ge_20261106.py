"""Package G-e tests: B-E review findings D1-01/02, D2-01..04, D3-02..05, D4-02/04/05.

Every test names the review item it pins. The theme of the package is
evidence-state inference honesty (BA vs EL endpoints, missing vs measured
zeros), a single total priority-cap ladder with profile-versioned cutoffs,
validated legacy baselines and sample-scoped event grouping, and canonical
HLA / empty-text identity vocabularies.
"""

from pathlib import Path
import json

from neoag import scoring
from neoag.candidate_identity import candidate_identity, normalize_hla
from neoag.config import load_profile
from neoag.evidence_consensus import (
    DEFAULT_RULES,
    _event_output,
    _legacy_baseline_order_status,
    _pareto_derived_grades,
    _rank_key,
    _row_text,
    build_evidence_consensus,
    load_consensus_rules,
    rank_evidence_consensus,
    score_evidence_consensus,
)
from neoag.evidence_states import derive_all_states, derive_presentation_consensus
from neoag.scoring import (
    apply_priority_cap_value,
    cap_rank,
    merge_escape_flags,
    priority,
    priority_thresholds,
    score,
)
from neoag.status_words import first_text, is_null_text
from neoag.utils import read_tsv, write_tsv
from neoag.vaccine_events import build_vaccine_event_tables


def complete_row(peptide_id="P1"):
    """Same shape as tests/test_evidence_consensus.py complete_row()."""
    return {
        "peptide_id": peptide_id,
        "event_id": "E" + peptide_id[1:],
        "event_type": "SNV",
        "efficacy_score": "0.5",
        "l3_hla_presentation_score": "0.9",
        "presentation_gate_status": "PASS",
        "netmhcpan_mt_rank_el": "0.2",
        "mhcflurry_presentation_score": "0.8",
        "l3_hla_binding_score": "0.8",
        "presentation_evidence_grade": "A",
        "l3_expression_score": "0.7",
        "expression_evidence_status": "GENE_AND_TRANSCRIPT_SUPPORTED",
        "rna_alt_reads": "8",
        "rna_vaf": "0.12",
        "l3_clonality_score": "0.6",
        "ccf_status": "clonal_like",
        "l3_tumor_specificity_score": "0.9",
        "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
        "mutant_specificity_gate_status": "PASS",
        "l3_normal_tissue_safety_score": "0.8",
        "safety_evidence_completeness": "COMPLETE",
        "safety_status": "PASS",
        "l3_apm_integrity_score": "0.7",
        "appm_evidence_completeness": "HIGH",
        "restricting_hla_lost": "false",
    }


# ---------------------------------------------------------------------------
# D1-01: BA rank must never be judged against EL thresholds
# ---------------------------------------------------------------------------

def test_ba_good_el_missing_is_unassessed_not_single_tool():
    # Review D1-01: previously binding_rank=0.3 fell into the EL slot and
    # produced PRESENTATION_SINGLE_TOOL grade 2 against EL thresholds.
    result = derive_presentation_consensus({"binding_rank": "0.3"}, DEFAULT_RULES)
    assert result["state"] == "PRESENTATION_UNASSESSED"
    assert result["assessed"] is False
    assert result["reason_code"] == "PRESENTATION_EL_UNASSESSED_BA_OBSERVED"
    assert "BA rank=0.3" in result["reason"]


def test_ba_bad_el_missing_is_unassessed_not_measured_worst():
    # Old code: binding_rank=80 became PRESENTATION_WEAK grade 0 — a fabricated
    # "measured worst" from an unassessed endpoint.
    result = derive_presentation_consensus({"binding_rank": "80.0"}, DEFAULT_RULES)
    assert result["state"] == "PRESENTATION_UNASSESSED"
    assert result["assessed"] is False


def test_ba_with_real_mhcflurry_keeps_single_tool_and_flags_el():
    result = derive_presentation_consensus(
        {"binding_rank": "0.3", "mhcflurry_presentation_score": "0.8"}, DEFAULT_RULES,
    )
    assert result["state"] == "PRESENTATION_SINGLE_TOOL"
    assert result["assessed"] is True
    assert "BA=0.3 not substituted" in result["reason"]


def test_el_present_drives_classification_regardless_of_ba():
    result = derive_presentation_consensus(
        {"netmhcpan_mt_rank_el": "3.0", "binding_rank": "0.1"}, DEFAULT_RULES,
    )
    # EL 3.0 > 2.0 supported threshold and no other core group -> weak.
    assert result["state"] == "PRESENTATION_WEAK"


def test_consensus_tiebreak_column_no_longer_receives_ba():
    # Review D3-02 (same endpoint separation): a binding_rank-only row must
    # leave the EL tiebreak slot as NA instead of 0.3.
    normalized = score_evidence_consensus(
        {"peptide_id": "P1", "event_id": "E1", "event_type": "SNV", "binding_rank": "0.3"},
    )
    assert normalized["netmhcpan_tiebreak_rank"] == "NA"
    assert normalized["presentation_consensus_state"] == "PRESENTATION_UNASSESSED"
    assert normalized["presentation_consensus_reason_code"] == "PRESENTATION_EL_UNASSESSED_BA_OBSERVED"


# ---------------------------------------------------------------------------
# D1-02: missing vs genuinely weak survives end-to-end into D3
# ---------------------------------------------------------------------------

def test_missing_presentation_is_unassessed_end_to_end():
    row = complete_row("P1")
    row["netmhcpan_mt_rank_el"] = ""
    row["mhcflurry_presentation_score"] = ""
    result = score_evidence_consensus(row)
    assert result["presentation_consensus_state"] == "PRESENTATION_UNASSESSED"
    assert "presentation_consensus" in result["evidence_missing_layers"]
    # Not-assessed must not materialise as the "worst measured" verdict:
    assert result["evidence_grade_uncapped"] != "R4"


def test_measured_worst_presentation_is_r4_end_to_end():
    row = complete_row("P2")
    row["netmhcpan_mt_rank_el"] = "99.0"  # a REAL measured poor prediction
    row["mhcflurry_presentation_score"] = "0.0"
    result = score_evidence_consensus(row)
    assert result["presentation_consensus_state"] == "PRESENTATION_WEAK"
    assert result["evidence_grade_uncapped"] == "R4"
    assert "presentation_consensus" not in result["evidence_missing_layers"]


# ---------------------------------------------------------------------------
# D2-01 / D2-02: one strictest-wins cap algebra
# ---------------------------------------------------------------------------

def test_unsafe_cap_helpers_are_gone():
    # Review D2-01/D2-02: the overwriting merge and the second ordering were
    # dead code; deleting them prevents re-wiring the wrong helper back in.
    assert not hasattr(scoring, "merge_peptide_escape")
    assert not hasattr(scoring, "apply_priority_cap")  # (apply_priority_cap_value stays)
    assert not hasattr(scoring, "_priority_rank")
    assert not hasattr(scoring, "_PRIORITY_ORDER")


def test_escape_merge_never_lifts_a_stricter_existing_cap():
    # Review D2-01: the live escape merge must keep the stricter earlier cap.
    merged = merge_escape_flags(
        {"peptide_id": "P1", "priority_cap": "C"},
        {"escape_status": "ESCAPE_CAUTION", "priority_cap": "A", "escape_multiplier": "0.9"},
    )
    assert merged["priority_cap"] == "C"
    merged_empty = merge_escape_flags(
        {"peptide_id": "P1", "priority_cap": "C"},
        {"escape_status": "ESCAPE_CAUTION", "priority_cap": "", "escape_multiplier": "0.9"},
    )
    assert merged_empty["priority_cap"] == "C"
    merged_stronger = merge_escape_flags(
        {"peptide_id": "P1", "priority_cap": "C"},
        {"escape_status": "ESCAPE_REJECT", "priority_cap": "D", "escape_multiplier": "0.0"},
    )
    assert merged_stronger["priority_cap"] == "D"


LADDER = ["D", "C_CAUTION", "C", "B_CAUTION", "B", "A"]


def test_cap_ladder_is_total_and_deterministic():
    # Review D2-02: exactly one ordering; a cap keeps any label at or below its
    # rank, and a plain label touched by a CAUTION cap becomes the caution
    # tier so the cap stays visible in final_priority.
    assert apply_priority_cap_value("C_CAUTION", "C") == "C_CAUTION"
    assert apply_priority_cap_value("C", "C_CAUTION") == "C_CAUTION"
    assert apply_priority_cap_value("B", "A") == "B"
    assert apply_priority_cap_value("B", "NONE") == "B"
    assert apply_priority_cap_value("B", "") == "B"
    assert apply_priority_cap_value("B", "SOMETHING_ELSE") == "B"  # off-ladder ignored
    for pri in LADDER:
        for cap in LADDER:
            expected = pri if cap_rank(pri) <= cap_rank(cap) else cap
            assert apply_priority_cap_value(pri, cap) == expected, (pri, cap)


def test_sequential_caps_are_order_independent():
    for base in LADDER:
        for first in LADDER:
            for second in LADDER:
                forward = apply_priority_cap_value(apply_priority_cap_value(base, first), second)
                backward = apply_priority_cap_value(apply_priority_cap_value(base, second), first)
                assert forward == backward, (base, first, second)


# ---------------------------------------------------------------------------
# D2-03: priority cutoffs are profile-locked, versioned, defaults unchanged
# ---------------------------------------------------------------------------

def test_priority_default_labels_are_pinned_at_boundaries():
    assert priority("PASS", 0.75) == "A"
    assert priority("PASS", 0.749) == "B"
    assert priority("PASS", 0.55) == "B"
    assert priority("PASS", 0.549) == "C"
    assert priority("PASS", 0.35) == "C"
    assert priority("PASS", 0.349) == "D"
    assert priority("CAUTION", 0.55) == "B_CAUTION"
    assert priority("CAUTION", 0.549) == "C_CAUTION"
    assert priority("SAFETY_PARTIAL", 0.35) == "C_CAUTION"
    assert priority("SAFETY_PARTIAL", 0.349) == "D"
    assert priority("FAIL", 0.9) == "D"


def test_shipped_profile_priority_thresholds_reproduce_history():
    profile = load_profile("default")
    section = profile["priority_thresholds"]
    assert section["score_a"] == 0.75
    assert section["score_b"] == 0.55
    assert section["score_c"] == 0.35
    assert str(section["version"])
    # With the shipped profile the labels are identical to the historical ones.
    assert priority("PASS", 0.75, profile) == "A"
    assert priority("PASS", 0.749, profile) == "B"


def test_priority_threshold_ablation_shifts_only_the_boundary():
    # Review D2-03 asks for a controlled ablation demonstration: raising no
    # constant, only a profile override, moves a score across a label.
    profile = {"priority_thresholds": {"score_a": 0.50, "version": "ablation"}}
    assert priority("PASS", 0.6) == "B"
    assert priority("PASS", 0.6, profile) == "A"
    assert priority_thresholds({"priority_thresholds": {"score_a": "oops"}})["score_a"] == 0.75


def test_peptide_output_records_priority_policy_version(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    presentation = tmp_path / "presentation.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV", "gene": "B2M",
        "event_confidence": "0.7", "event_expression": "10", "driver_relevance": "0.5",
        "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.7",
        "safety_status": "PASS",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "gene": "B2M", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "1.0", "el_rank": "1.0", "presentation_score": "0.8",
    }])
    write_tsv(presentation, [{"peptide_id": "P1"}])
    out_events = tmp_path / "ranked_events.tsv"
    out_peptides = tmp_path / "ranked_peptides.tsv"
    _, peps = score(
        events, peptides, presentation, None, None, None, None,
        load_profile("default"), out_events, out_peptides,
    )
    assert peps[0]["priority_policy_version"] == "priority-thresholds-v1-2026-11-06"


# ---------------------------------------------------------------------------
# D2-04: orphan peptide rows are ledgered, never silently dropped
# ---------------------------------------------------------------------------

def _write_scoring_inputs(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    presentation = tmp_path / "presentation.tsv"
    write_tsv(events, [
        {"event_id": "E1", "sample_id": "S1", "event_type": "SNV", "gene": "B2M",
         "event_confidence": "0.7", "event_expression": "10", "driver_relevance": "0.5",
         "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.7",
         "safety_status": "PASS"},
        {"event_id": "E2", "sample_id": "S1", "event_type": "Fusion", "gene": "FUS",
         "event_confidence": "0.7", "event_expression": "10", "driver_relevance": "0.5",
         "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.7",
         "safety_status": "PASS"},
    ])
    write_tsv(peptides, [
        {"peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
         "gene": "B2M", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I"},
        {"peptide_id": "P2", "event_id": "E_NOPE", "sample_id": "S1", "event_type": "SNV",
         "gene": "B2M", "peptide": "WWWWWWWWW", "hla_allele": "HLA-A*02:01", "mhc_class": "I"},
        {"peptide_id": "P3", "event_id": "E2", "sample_id": "S1", "event_type": "Fusion",
         "gene": "FUS", "peptide": "QQQQQQQQQ", "hla_allele": "HLA-B*07:02", "mhc_class": "I"},
    ])
    write_tsv(presentation, [{"peptide_id": "P1"}, {"peptide_id": "P2"}, {"peptide_id": "P3"}])
    return events, peptides, presentation


def test_score_ledgers_orphan_peptides(tmp_path):
    events, peptides, presentation = _write_scoring_inputs(tmp_path)
    profile = load_profile("default")
    profile = dict(profile)
    profile["enforce_enabled_sources"] = True
    profile["enabled_sources"] = {"types": ["SNV"]}
    out_events = tmp_path / "ranked_events.tsv"
    out_peptides = tmp_path / "ranked_peptides.tsv"
    _, peps = score(
        events, peptides, presentation, None, None, None, None,
        profile, out_events, out_peptides,
    )
    ledger = read_tsv(tmp_path / "ranked_peptides.unmatched_peptides.tsv")
    assert {row["peptide_id"] for row in ledger} == {"P2", "P3"}
    reasons = {row["peptide_id"]: row["reason"] for row in ledger}
    assert reasons["P2"] == "event_id_absent_from_raw_events"
    assert reasons["P3"] == "event_filtered_out_by_enabled_sources"
    assert [row["peptide_id"] for row in peps] == ["P1"]


def test_score_writes_header_only_ledger_when_nothing_is_lost(tmp_path):
    events, peptides, presentation = _write_scoring_inputs(tmp_path)
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "gene": "B2M", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    out_peptides = tmp_path / "ranked_peptides.tsv"
    score(
        events, peptides, presentation, None, None, None, None,
        load_profile("default"), tmp_path / "ranked_events.tsv", out_peptides,
    )
    ledger_path = tmp_path / "ranked_peptides.unmatched_peptides.tsv"
    assert ledger_path.is_file()
    assert read_tsv(ledger_path) == []  # header written, zero rows = visible count


# ---------------------------------------------------------------------------
# D3-02: placeholder numerics must not fake a favourable tie-break
# ---------------------------------------------------------------------------

def _tie_row(peptide_id, netmhcpan, mhcflurry):
    row = {
        "evidence_grade": "R3", "evidence_track": "MISSENSE", "pareto_front": "1",
        "fusion_orf_gate_grade": "0", "source_chain_confidence_grade": "0",
        "safety_completeness_grade": "0", "event_authenticity_grade": "0",
        "rna_support_grade": "0", "presentation_consensus_grade": "0",
        "mutant_specificity_grade": "0", "hla_appm_grade": "0",
        "ccf_confidence_grade": "0", "evidence_completeness_grade": "0",
        "netmhcpan_tiebreak_rank": netmhcpan, "mhcflurry_tiebreak_score": mhcflurry,
        "peptide_id": peptide_id,
    }
    return row


def test_tiebreak_orders_measured_values_then_missing_last():
    rows = [
        _tie_row("P_missing", "NA", "NA"),
        _tie_row("P_zero", "0.0", "NA"),
        _tie_row("P_low", "0.1", "NA"),
        _tie_row("P_worst", "1.9", "NA"),
    ]
    order = [row["peptide_id"] for row in sorted(rows, key=lambda row: _rank_key(row, DEFAULT_RULES))]
    assert order == ["P_zero", "P_low", "P_worst", "P_missing"]


def test_tiebreak_mhcflurry_missing_sorts_after_measured_zero():
    rows = [
        _tie_row("P_missing", "NA", "NA"),
        _tie_row("P_zero", "NA", "0.0"),
        _tie_row("P_good", "NA", "0.9"),
    ]
    order = [row["peptide_id"] for row in sorted(rows, key=lambda row: _rank_key(row, DEFAULT_RULES))]
    assert order == ["P_good", "P_zero", "P_missing"]


# ---------------------------------------------------------------------------
# D3-03: normal-junction sources evaluated separately, contradictions kept
# ---------------------------------------------------------------------------

def _fusion_normal_junction_grades(**fields):
    source = {"event_type": "Fusion", "peptide_consequence": "fusion", **fields}
    states = derive_all_states(source, DEFAULT_RULES)
    grades = _pareto_derived_grades(source, states, "FUSION", load_consensus_rules())
    return grades["normal_junction_safety_grade"]


def test_normal_junction_mixed_sources_are_conflict_not_safe():
    # Review D3-03: EXACT_MATCH in one field plus NOT_DETECTED in another used
    # to land in the safe branch (grade 3). It must not be a safe verdict.
    assert _fusion_normal_junction_grades(
        normal_junction_exact_match="EXACT_MATCH",
        normal_junction_assessment_status="NOT_DETECTED",
    ) == 1


def test_normal_junction_single_source_grades():
    assert _fusion_normal_junction_grades(normal_junction_exact_match="EXACT_MATCH") == 0
    assert _fusion_normal_junction_grades(normal_junction_status="NOT_DETECTED") == 3
    assert _fusion_normal_junction_grades(normal_junction_status="ABSENT") == 3
    assert _fusion_normal_junction_grades(normal_junction_status="UNASSESSED") == 1
    assert _fusion_normal_junction_grades() == 1
    assert _fusion_normal_junction_grades(normal_junction_status="REVIEW_PENDING") == 2


# ---------------------------------------------------------------------------
# D3-04: the legacy baseline must be validated before rank statistics
# ---------------------------------------------------------------------------

def test_legacy_order_validator_accepts_monotonic_inputs():
    assert _legacy_baseline_order_status([]) [0] == "NO_INPUT_ROWS"
    rows = [{"weighted_rank": "1"}, {"weighted_rank": "2"}, {"weighted_rank": "3"}]
    assert _legacy_baseline_order_status(rows)[0] == "VALIDATED"
    scores = [{"efficacy_score": "0.9"}, {"efficacy_score": "0.9"}, {"efficacy_score": "0.1"}]
    assert _legacy_baseline_order_status(scores)[0] == "VALIDATED"


def test_legacy_order_validator_rejects_unranked_and_ascending_inputs():
    status, detail = _legacy_baseline_order_status([{"efficacy_score": "0.4"}, {"efficacy_score": "0.9"}])
    assert status == "INVALID"
    assert "row 2" in detail
    status, detail = _legacy_baseline_order_status([{"peptide_id": "P1"}, {"peptide_id": "P2"}])
    assert status == "INVALID"
    assert "cannot be certified" in detail
    status, _ = _legacy_baseline_order_status([{"weighted_rank": "2"}, {"weighted_rank": "1"}])
    assert status == "INVALID"


def test_comparison_markdown_reports_validation_state(tmp_path):
    good = complete_row("P1")
    bad = complete_row("P2")
    good["efficacy_score"] = "0.9"
    bad["efficacy_score"] = "0.4"
    comprehensive = tmp_path / "comprehensive.tsv"
    write_tsv(comprehensive, [good, bad])
    rank_evidence_consensus(
        comprehensive,
        tmp_path / "ranked.tsv", tmp_path / "events.tsv",
        tmp_path / "states.tsv", tmp_path / "conflicts.tsv",
        load_consensus_rules(),
    )
    markdown = (tmp_path / "ranking_compare_weighted_vs_consensus.md").read_text()
    assert "VALIDATED" in markdown
    assert "NOT_COMPARABLE" not in markdown
    assert "Legacy baseline validation" in markdown

    comprehensive2 = tmp_path / "comprehensive_unsorted.tsv"
    write_tsv(comprehensive2, [dict(bad, efficacy_score="0.4"), dict(good, efficacy_score="0.9")])
    rank_evidence_consensus(
        comprehensive2,
        tmp_path / "ranked2.tsv", tmp_path / "events2.tsv",
        tmp_path / "states2.tsv", tmp_path / "conflicts2.tsv",
        load_consensus_rules(),
    )
    markdown2 = (tmp_path / "ranking_compare_weighted_vs_consensus.md").read_text()
    assert "INVALID" in markdown2
    assert "NOT_COMPARABLE (legacy row order is not a validated weighted ranking)" in markdown2


def test_missing_weighted_baseline_recorded_in_manifest(tmp_path):
    outdir = tmp_path / "consensus"
    comprehensive = tmp_path / "comprehensive.tsv"
    write_tsv(comprehensive, [complete_row("P1"), complete_row("P2")])
    result = build_evidence_consensus(comprehensive, outdir)
    assert result["weighted_baseline_status"] == "MISSING"
    manifest = json.loads((outdir / "evidence_consensus_run.json").read_text())
    assert manifest["weighted_baseline"]["status"] == "MISSING"
    assert manifest["legacy_baseline"]["order_status"] in {"VALIDATED", "INVALID"}


# ---------------------------------------------------------------------------
# D3-05: event grouping is sample-scoped as soon as sample_id exists
# ---------------------------------------------------------------------------

def _group_row(peptide_id, event_id, sample_id, rank):
    return {
        "peptide_id": peptide_id, "event_id": event_id, "sample_id": sample_id,
        "hla_allele": "HLA-A*02:01", "evidence_rank": str(rank),
        "evidence_track": "MISSENSE", "evidence_grade": "R3", "pareto_front": "1",
        "evidence_rank_key": f"key{rank}", "consensus_action": "X",
        "recommended_next_steps": "", "consensus_trace": "t",
        "hard_failure": "no", "manual_review_required": "no",
    }


def test_event_groups_are_isolated_by_sample():
    rows = [
        _group_row("P1", "E1", "S1", 1),
        _group_row("P2", "E1", "S2", 2),  # same event id, different sample
    ]
    groups = _event_output(rows, True)
    assert {row["event_group_id"] for row in groups} == {"EVENT:S1|E1", "EVENT:S2|E1"}
    assert all(row["member_event_count"] == "1" for row in groups)


def test_event_groups_keep_legacy_key_without_sample_id():
    rows = [_group_row("P1", "E1", "", 1), _group_row("P2", "E1", "", 2)]
    groups = _event_output(rows, True)
    assert [row["event_group_id"] for row in groups] == ["EVENT:E1"]


def test_unassessed_event_id_cannot_glue_peptides_together():
    # D4-05 vocabulary × D3-05 grouping: an "UNASSESSED" event id is not an
    # identity, so each peptide gets its own NO_EVENT group instead of sharing
    # one bogus event.
    rows = [
        _group_row("P1", "UNASSESSED", "S1", 1),
        _group_row("P2", "UNASSESSED", "S1", 2),
    ]
    groups = _event_output(rows, True)
    assert {row["event_group_id"] for row in groups} == {
        "EVENT:S1|NO_EVENT:P1", "EVENT:S1|NO_EVENT:P2",
    }


# ---------------------------------------------------------------------------
# D4-02: canonical HLA normalization in vaccine event grouping
# ---------------------------------------------------------------------------

def test_vaccine_grouping_collapses_hla_prefix_variants():
    events = [{"event_id": "E1", "event_group_id": "EVENT:E1", "event_type": "SNV",
               "gene": "GENE1", "event_evidence_rank": "1"}]
    peptides = [
        {"peptide_id": "P1", "event_id": "E1", "peptide": "SYFPEITHI",
         "hla_allele": "A*02:01", "evidence_rank": "1"},
        {"peptide_id": "P2", "event_id": "E1", "peptide": "SYFPEITHI",
         "hla_allele": "HLA-A*02:01", "evidence_rank": "2"},
    ]
    summaries, details = build_vaccine_event_tables(events, peptides)
    assert len(summaries) == 1
    summary = summaries[0]
    # Whitespace-only normalization used to report these as two unique
    # peptide-HLA units inside one event (D4-02 inflation).
    assert summary["unique_neoepitope_hla_count"] == "1"
    assert summary["raw_peptide_hla_row_count"] == "2"
    assert summary["hla_coverage_count"] == "1"
    assert summary["hla_alleles"] == "HLA-A*02:01"
    assert len(details) == 1
    detail = details[0]
    assert detail["hla_allele"] == "A*02:01"          # original representation kept
    assert detail["hla_allele_canonical"] == "HLA-A*02:01"
    assert detail["duplicate_evidence_row_count"] == "2"


def test_normalize_hla_canonical_forms():
    assert normalize_hla(" hla-a*02:01 ") == "HLA-A*02:01"
    assert normalize_hla("A*02:01") == "HLA-A*02:01"
    assert normalize_hla("A*02:01:01:02") == normalize_hla("a*02:01:01:02")
    assert normalize_hla("") == ""


# ---------------------------------------------------------------------------
# D4-04: empty-identity rows no longer collide on one hashed "|||||" id
# ---------------------------------------------------------------------------

def test_empty_coordinate_rows_get_distinct_unresolved_ids():
    row_a = {"peptide": "AAAAAAA", "hla_allele": "A*02:01", "source_note": "run1-well7"}
    row_b = {"peptide": "AAAAAAA", "hla_allele": "A*02:01", "source_note": "run1-well8"}
    id_a = candidate_identity(row_a)["event_identity_id"]
    id_b = candidate_identity(row_b)["event_identity_id"]
    assert id_a != id_b
    # Deterministic across calls for the same row content.
    assert candidate_identity(row_a)["event_identity_id"] == id_a
    # With coordinates present the coordinate fallback still collapses identically.
    coord = {"chrom": "chr1", "pos": "100", "ref": "A", "alt": "G", "event_type": "SNV", "gene": "G1"}
    assert (
        candidate_identity(coord)["event_identity_id"]
        == candidate_identity(dict(coord, peptide="ZZZZZZZ"))["event_identity_id"]
    )


# ---------------------------------------------------------------------------
# D4-05: one centralized null/UNASSESSED vocabulary
# ---------------------------------------------------------------------------

def test_shared_null_text_vocabulary():
    for value in ("", "  ", "NA", "n/a", "None", ".", "UNASSESSED"):
        assert is_null_text(value), value
    for value in ("E1", "0", "NOT_APPLICABLE", "HLA-A*02:01"):
        assert not is_null_text(value), value
    assert first_text({"a": "UNASSESSED", "b": "real"}, "a", "b") == "real"
    assert first_text({"a": "NA"}, "a", default="fallback") == "fallback"


def test_row_text_and_identity_share_the_vocabulary():
    # evidence_consensus._row_text used to ACCEPT "UNASSESSED" as a value while
    # candidate_identity._text rejected it (D4-05 divergence).
    assert _row_text({"event_id": "UNASSESSED"}, "event_id") == ""
    row = {
        "event_id": "UNASSESSED", "event_type": "SNV", "gene": "G1",
        "chrom": "chr1", "pos": "100", "ref": "A", "alt": "G",
        "peptide": "SYFPEITHI", "hla_allele": "A*02:01",
    }
    as_unassessed = candidate_identity(row)
    as_empty = candidate_identity(dict(row, event_id=""))
    assert as_unassessed == as_empty  # both fall through to coordinate identity
