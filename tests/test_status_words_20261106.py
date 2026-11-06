"""Regression tests for the shared status-polarity matcher (Package A).

Review findings covered: B5-01 (_truthy positive-first), B5-02 and D1-05
(dead classifiers removed), D1-03 ("ALT_SUPPORTED" in "NOT_ALT_SUPPORTED"),
D1-04 ("PASS" in "NOT_PASSED"), plus the same bug class found at the other
live positive-branch substring sites in source_chain/evidence_states.
"""

from neoag import evidence_states, source_chain
from neoag.evidence_states import (
    derive_event_authenticity,
    derive_mutant_specificity,
    derive_safety_state,
)
from neoag.source_chain import SUPPORTED, NEGATIVE, UNASSESSED, _somatic_status_requirement
from neoag.status_words import contains_positive, starts_negated


def test_contains_positive_rejects_negated_occurrences():
    assert contains_positive("ALT_SUPPORTED", "ALT_SUPPORTED")
    assert not contains_positive("NOT_ALT_SUPPORTED", "ALT_SUPPORTED")
    assert not contains_positive("NOT ALT SUPPORTED", "ALT_SUPPORTED")  # spaces normalise
    assert not contains_positive("UNSUPPORTED", "SUPPORTED")
    assert not contains_positive("NOT_PASSED", "PASS")
    assert not contains_positive("UNPHASED", "PHASED")
    assert not contains_positive("NOT_SOMATIC", "SOMATIC")
    assert not contains_positive("NOT_UNSAFE", "UNSAFE")  # double negative is not a reject
    assert contains_positive("UNSAFE", "UNSAFE")
    assert contains_positive("FILTER: PASS; REVIEW: none", "PASS")
    assert contains_positive("NOVEL_TAIL", "NOVEL")  # NO_ negation needs the underscore
    assert not contains_positive("", "PASS")
    assert not contains_positive("PASS", "")


def test_starts_negated_keeps_absence_semantics():
    assert starts_negated("NOT_SUPPORTED")
    assert starts_negated("NOT PASSED")
    assert starts_negated("NO_SUPPORT")
    assert not starts_negated("NOT_AVAILABLE")  # absence, not a negative verdict
    assert not starts_negated("NOT_APPLICABLE")
    assert not starts_negated("PASS")
    assert not starts_negated("")


def test_truthy_reads_negated_values_as_false():
    # B5-01: the old order matched SUPPORTED/PASS substrings first, so an
    # explicit negative flipped to True.
    assert source_chain._truthy({"crosses_junction": "NOT_SUPPORTED"}, "crosses_junction") is False
    assert source_chain._truthy({"crosses_junction": "NOT_PASSED"}, "crosses_junction") is False
    assert source_chain._truthy({"crosses_junction": "NOT CONFIRMED"}, "crosses_junction") is False
    assert source_chain._truthy({"crosses_junction": "PASS"}, "crosses_junction") is True
    assert source_chain._truthy({"crosses_junction": "CONFIRMED"}, "crosses_junction") is True
    assert source_chain._truthy({"crosses_junction": "NOT_AVAILABLE"}, "crosses_junction") is None
    assert source_chain._truthy({}, "crosses_junction") is None


def test_somatic_requirement_honours_negated_normal_status():
    # D1-03 class, live in source_chain: NOT_ALT_SUPPORTED must not mark the
    # candidate as matched-normal supported (fatal), while ALT_SUPPORTED must.
    row = {"matched_normal_status": "NOT_ALT_SUPPORTED", "normal_depth": "42", "normal_alt_count": "0"}
    assert _somatic_status_requirement(row, {}).status == SUPPORTED
    row_bad = {"matched_normal_status": "ALT_SUPPORTED", "normal_depth": "42"}
    assert _somatic_status_requirement(row_bad, {}).status == NEGATIVE
    # NOT_SOMATIC must not read as SOMATIC via the substring "SOMATIC".
    row_germ = {"matched_normal_status": "NOT_SOMATIC", "normal_depth": "42", "normal_alt_count": "1"}
    assert _somatic_status_requirement(row_germ, {}).status == UNASSESSED


def test_event_authenticity_negated_normal_status_not_hard_fail():
    # D1-03: "NOT_ALT_SUPPORTED" is the best tumor-specificity evidence and
    # must not trigger HARD_MATCHED_NORMAL_SUPPORT.
    row = {"matched_normal_status": "NOT_ALT_SUPPORTED", "normal_alt_count": "0", "normal_alt_vaf": "0.0"}
    res = derive_event_authenticity(row, {})
    assert not res["hard_fail"]
    row_bad = {"matched_normal_status": "ALT_SUPPORTED"}
    res_bad = derive_event_authenticity(row_bad, {})
    assert res_bad["hard_fail"] and res_bad["hard_code"] == "HARD_MATCHED_NORMAL_SUPPORT"


def test_mutant_specificity_negated_gate_is_not_mt_specific():
    # D1-04: gate "NOT_PASSED" contains "PASS" and previously produced the
    # top-grade MT_SPECIFIC state.
    row = {
        "peptide": "AXMX",
        "wildtype_peptide": "AXLX",
        "mutation_positions_in_peptide": "3",
        "mutant_specificity_gate_status": "NOT_PASSED",
    }
    res = derive_mutant_specificity(row, {})
    assert res["state"] != "MT_SPECIFIC"
    row_ok = dict(row, mutant_specificity_gate_status="PASS")
    assert derive_mutant_specificity(row_ok, {})["state"] == "MT_SPECIFIC"


def test_safety_state_requires_positive_pass():
    assert derive_safety_state({"safety_status": "PASS"}, {})["state"] == "SAFETY_PASS"
    assert derive_safety_state({"safety_status": "NOT_PASSED"}, {})["state"] != "SAFETY_PASS"
    assert derive_safety_state({"safety_status": "NOT_UNSAFE"}, {})["state"] != "SAFETY_REJECT"
    assert derive_safety_state({"safety_status": "UNSAFE"}, {})["hard_code"] == "HARD_SAFETY_REJECT"


def test_phasing_unphased_is_not_supported():
    # "UNPHASED" contains "PHASED"; it must not satisfy the read-backed
    # phasing requirement.
    row = {"phasing_required": "TRUE", "haplotype_status": "UNPHASED"}
    from neoag.source_chain import _phasing_requirement

    assert _phasing_requirement(row).status == UNASSESSED


def test_dead_classifiers_are_removed():
    # B5-02 / D1-05: both helpers had zero production call sites repo-wide;
    # the live derive_* functions replace them.
    assert not hasattr(evidence_states, "evidence_state")
    assert not hasattr(evidence_states, "STATES")
    assert not hasattr(evidence_states, "EvidenceState")
    assert not hasattr(evidence_states, "SUPPORTED_TOKENS")
    assert not hasattr(source_chain, "_status_from_text")
    assert not hasattr(source_chain, "_NEGATIVE_TOKENS")
