"""Typed scopes for hard-fail reason codes (review E-03).

A hard-fail gate must block exactly the descendants a defect invalidates.
An exact normal-proteome match for one peptide invalidates that peptide,
not its sibling peptides from the same event; a sample-identity failure
invalidates every candidate from the sample.  Each hard-fail code therefore
carries one of five scopes:

SAMPLE_GLOBAL
    Every candidate from the sample (identity/mix-up failures).
EVENT_GLOBAL
    Every candidate from the event (authenticity, artefact, germline
    matched-normal support, coordinate/build failures, normal junctions).
ORF_SPECIFIC
    Candidates translated from one defective reading frame (non-mutant ORF).
    Without an explicit ORF identifier on the row the scope is applied to
    the reported peptide row; it never escalates to the whole event.
PEPTIDE_SPECIFIC
    The reported peptide sequence, under every HLA restriction.
PEPTIDE_HLA_SPEC
    The reported peptide + HLA pair only (e.g. restricting allele confirmed
    lost — the same peptide restricted by a surviving allele is unaffected).

Unknown codes default to EVENT_GLOBAL: a new, unclassified reason keeps the
pre-E-03 conservative behaviour of blocking the whole event until someone
types it here deliberately.  Safety is never made looser by an omission.

Scopes are reported per event through the consensus event table
(``hard_failure_event_blocking`` / ``hard_failure_event_scope_codes`` /
``hard_failure_peptide_local_codes``) and proved through
review_integrity's ``hard_fail_not_promoted`` check.
"""

from __future__ import annotations

SAMPLE_GLOBAL = "SAMPLE_GLOBAL"
EVENT_GLOBAL = "EVENT_GLOBAL"
ORF_SPECIFIC = "ORF_SPECIFIC"
PEPTIDE_SPECIFIC = "PEPTIDE_SPECIFIC"
PEPTIDE_HLA_SPEC = "PEPTIDE_HLA_SPEC"

#: Scopes whose codes keep the event itself hard-excluded.
EVENT_BLOCKING_SCOPES = frozenset({SAMPLE_GLOBAL, EVENT_GLOBAL})

# Code -> scope.  Codes verified against their emission sites in
# evidence_states.py and source_chain.py (2026-11-07 review E-03).
HARD_FAIL_SCOPES: dict[str, str] = {
    # Sample identity / provenance: everything from this sample is affected.
    "SC_SAMPLE_IDENTITY_INVALID": SAMPLE_GLOBAL,
    # Event authenticity and coordinates: every peptide of the event.
    "HARD_EVENT_ARTIFACT": EVENT_GLOBAL,
    "SC_EVENT_ARTIFACT": EVENT_GLOBAL,
    "HARD_NORMAL_JUNCTION": EVENT_GLOBAL,
    "HARD_MATCHED_NORMAL_SUPPORT": EVENT_GLOBAL,
    "SC_COORDINATE_OR_ALLELE_INVALID": EVENT_GLOBAL,
    "SC_FORCE_C4_OVERRIDE": EVENT_GLOBAL,
    # Reading-frame defects: peptides from the defective ORF only.
    "HARD_NON_MUTANT_SEQUENCE": ORF_SPECIFIC,
    "SC_PEPTIDE_NOT_EVENT_DERIVED": ORF_SPECIFIC,
    # Peptide-sequence defects (this exact sequence only).
    "HARD_REFERENCE_PROTEOME_MATCH": PEPTIDE_SPECIFIC,
    "HARD_SAFETY_REJECT": PEPTIDE_SPECIFIC,
    "SC_FUSION_BOUNDARY_MAPPING_INCOMPLETE": PEPTIDE_SPECIFIC,
    # Peptide + HLA pair defects.
    "HARD_RESTRICTING_HLA_LOST": PEPTIDE_HLA_SPEC,
    "SC_PEPTIDE_HLA_TRACEABILITY_INVALID": PEPTIDE_HLA_SPEC,
    "SC_PEPTIDE_HLA_TRACEABILITY_CONFLICT": PEPTIDE_HLA_SPEC,
}


def hard_fail_scope(code: str) -> str:
    """Scope for one hard-fail code; unknown codes are treated as event-global."""
    return HARD_FAIL_SCOPES.get(str(code).strip().upper(), EVENT_GLOBAL)


def _split_codes(codes) -> list[str]:
    if isinstance(codes, str):
        items: list[str] = [code.strip().upper() for code in codes.replace(";", ",").split(",")]
    else:
        items = [str(code).strip().upper() for code in codes]
    return [code for code in items if code]


def event_blocking_codes(codes) -> list[str]:
    """Codes whose scope invalidates the whole event (order-preserving)."""
    return list(dict.fromkeys(
        code for code in _split_codes(codes) if hard_fail_scope(code) in EVENT_BLOCKING_SCOPES
    ))


def local_only_codes(codes) -> list[str]:
    """Codes that invalidate only the reported peptide/ORF/peptide-HLA row."""
    return list(dict.fromkeys(
        code for code in _split_codes(codes) if hard_fail_scope(code) not in EVENT_BLOCKING_SCOPES
    ))


def scopes_text(codes) -> str:
    """Human-auditable ``CODE:SCOPE`` pairs for the row's hard-fail codes."""
    return ";".join(f"{code}:{hard_fail_scope(code)}" for code in dict.fromkeys(_split_codes(codes)))
