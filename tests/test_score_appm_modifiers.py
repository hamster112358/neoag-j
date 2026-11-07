from pathlib import Path

from neoag.config import load_profile
from neoag.scoring import resolve_appm_peptide_modifiers_tsv, score
from neoag.utils import write_tsv


def test_resolve_appm_peptide_modifiers_from_summary_sibling(tmp_path):
    appm = tmp_path / "appm"
    appm.mkdir()
    summary = appm / "appm_summary.tsv"
    modifiers = appm / "appm_peptide_modifiers.tsv"
    summary.write_text("sample_id\nS1\n", encoding="utf-8")
    modifiers.write_text("peptide_id\tappm_multiplier\nP1\t0.5000\n", encoding="utf-8")
    resolved = resolve_appm_peptide_modifiers_tsv(appm_summary_tsv=summary)
    assert resolved == str(modifiers)


def test_score_applies_appm_peptide_modifiers(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    presentation = tmp_path / "presentation.tsv"
    appm = tmp_path / "appm"
    appm.mkdir()
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
        "immunogenicity_score": "0.6", "wildtype_binding_rank": "5.0",
        "self_similarity_score": "0.1", "normal_hla_ligand_overlap": "no",
    }])
    write_tsv(presentation, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_evidence_score": "0.8", "presentation_evidence_score": "0.7",
        "presentation_evidence_grade": "A",
    }])
    write_tsv(appm / "appm_summary.tsv", [{
        "mhc_i_integrity_score": "1.0", "mhc_ii_integrity_score": "1.0",
        "mhc_i_integrity_status": "MHC_I_INTACT", "mhc_ii_integrity_status": "MHC_II_INTACT",
    }])
    write_tsv(appm / "appm_peptide_modifiers.tsv", [{
        "peptide_id": "P1", "event_id": "E1", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "appm_multiplier": "0.2500", "appm_multiplier_reason": "mhc_i_defective",
        "appm_integrity_status": "MHC_I_DEFECTIVE", "appm_evidence_completeness": "LOW",
        "appm_review_required": "yes", "priority_cap": "C", "appm_action": "CAP",
    }])
    out_events = tmp_path / "ranked_events.tsv"
    out_peptides = tmp_path / "ranked_peptides.tsv"
    _, peps = score(
        events, peptides, presentation, appm / "appm_summary.tsv", None, None, None,
        load_profile("default"), out_events, out_peptides,
        appm_peptide_modifiers_tsv=appm / "appm_peptide_modifiers.tsv",
    )
    assert peps[0]["appm_multiplier"] == "0.2500"
    assert peps[0]["appm_multiplier_reason"] == "mhc_i_defective"
    assert peps[0]["appm_integrity_status"] == "MHC_I_DEFECTIVE"
    # Review D2-02: with one total cap ladder (D < C_CAUTION < C < ...) the
    # APPM "C" cap and the missing-safety-evidence "C_CAUTION" partial cap now
    # merge to the STRICTER label instead of silently keeping whichever was
    # applied first (the old code tied C with C_CAUTION).
    assert peps[0]["priority_cap"] == "C_CAUTION"


def test_d2_05_zero_composite_is_measured_and_only_source_none_falls_back():
    # Round-2 (HANDOFF section I, D2-05): availability, not magnitude,
    # decides the fallback.  A resolved composite of exactly 0.0 is a
    # measured result and stays authoritative; only the absence of component
    # scores (immunogenicity_source "none") may fall through to the legacy
    # immunogenicity_score field.
    from neoag.config import load_profile
    from neoag.scoring import normalize_immunogenicity_score

    profile = load_profile("default")
    presentation_zero = {"bigmhc_im_score": "0.0000", "hla_allele": "HLA-A*02:01"}
    peptide = {"immunogenicity_score": "0.9000"}
    assert normalize_immunogenicity_score(peptide, presentation_zero, profile) == 0.0
    # No component scores -> source "none" -> legacy value still honoured.
    assert normalize_immunogenicity_score(dict(peptide), {}, profile) == 0.9
    # Without a profile the composite path is not consulted at all.
    assert normalize_immunogenicity_score(dict(peptide), presentation_zero, None) == 0.9


def test_d2_05_sample_level_loh_penalty_names_its_allele_string_provenance():
    # Round-2 (HANDOFF section I, D2-05): when the LOH penalty comes from the
    # sample summary's comma-separated allele list rather than a per-candidate
    # APPM modifier row, the reason must say so and name the allele.
    from neoag.config import load_profile
    from neoag.scoring import appm_multiplier

    profile = load_profile("default")
    peptide = {"hla_allele": "HLA-A*02:01"}
    summary = {"mhc_i_integrity_score": "1.0", "hla_loh_alleles": "HLA-A*02:01,HLA-B*07:02"}
    value = appm_multiplier(peptide, summary, profile)
    assert value < 1.0
    reason = peptide["appm_multiplier_reason"]
    assert reason.startswith("sample_level_appm_summary+hla_loh_allele_string:")
    assert "HLA-A*02:01" in reason or "A*02:01" in reason
    untouched = {"hla_allele": "HLA-C*07:02"}
    appm_multiplier(untouched, summary, profile)
    assert "appm_multiplier_reason" not in untouched
