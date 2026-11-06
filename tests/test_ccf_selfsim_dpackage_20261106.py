"""B-E review package D (2026-11-06).

C4-02: unresolved clonality states must not be penalised below a *measured*
subclone, and clonality_indeterminate gets its own tunable multiplier.
C6-01: adapters must not fabricate self_similarity_score "0.0"/"0" for rows
where no self-similarity was computed; the value stays blank and the new
self_similarity_status column says UNASSESSED.
"""

import warnings

from neoag.adapters.easyfuse_adapter import easyfuse_row_to_peptide
from neoag.adapters.pvactools_parser import peptide_from_row
from neoag.adapters.variant_peptide_adapter import peptide_from_variant_row
from neoag.ccf_v2 import build_ccf_2, clonality_multiplier
from neoag.config import load_profile
from neoag.peptide_safety_gate import build_peptide_safety_gate
from neoag.schemas import PEPTIDE_FIELDS
from neoag.sv.peptide_builder import PeptideCandidate
from neoag.sv.phase1 import _raw_peptide_from_candidate
from neoag.utils import write_tsv


# ---------------------------------------------------------------- C4-02 ----

def test_unresolved_statuses_are_neutral_not_penalised():
    prof = load_profile("default")
    subclonal = clonality_multiplier("subclonal_compatible", prof)
    assert clonality_multiplier("unresolved", prof) == 1.0
    assert clonality_multiplier("vaf_only_unresolved", prof) == 1.0
    # The old fall-through charged 0.65 -- worse than a measured subclone.
    assert clonality_multiplier("unresolved", prof) >= subclonal


def test_clonality_indeterminate_has_its_own_multiplier():
    prof = load_profile("default")
    assert clonality_multiplier("clonality_indeterminate", prof) == 0.75
    custom = {"ccf_lite": {"indeterminate_multiplier": 0.5}}
    assert clonality_multiplier("clonality_indeterminate", custom) == 0.5
    # It is now separable from the subclonal bucket for tuning.
    assert clonality_multiplier("subclonal_compatible", custom) == 0.75


def test_unknown_status_still_falls_back_to_missing_ccf_multiplier():
    assert clonality_multiplier("totally_unknown", {"ccf_lite": {}}) == 0.65


def test_rna_only_multiplier_compat_chain():
    # An explicit rna_only key wins (sarcoma_rna_supported_v2_provisional
    # pins 0.85 and keeps its profile choice).
    assert clonality_multiplier(
        "RNA_ONLY_UNRESOLVED", {"ccf_lite": {"rna_only_unresolved_multiplier": 0.85}}
    ) == 0.85
    # Otherwise it follows the shared unresolved knob.
    assert clonality_multiplier(
        "RNA_ONLY_UNRESOLVED", {"ccf_lite": {"unresolved_multiplier": 0.9}}
    ) == 0.9
    assert clonality_multiplier("RNA_ONLY_UNRESOLVED", {"ccf_lite": {}}) == 1.0


def test_all_profiles_declare_the_new_multipliers():
    names = (
        "default",
        "sarcoma",
        "leukemia",
        "sv_wgs_phase1",
        "sv_wes_phase1_5",
        "sarcoma_rna_supported",
        "sarcoma_rna_supported_v2_provisional",
        "hcc1395_bigmhc_only",
    )
    for name in names:
        cfg = load_profile(name).get("ccf_lite", {})
        assert cfg.get("unresolved_multiplier") == 1.0, name
        assert cfg.get("indeterminate_multiplier") == 0.75, name


def _write_ccf_inputs(tmp_path):
    events = tmp_path / "events.tsv"
    purity = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(purity, [{"purity": "0.60", "ploidy": "2"}])
    write_tsv(cnv, [{
        "chrom": "chr1", "start": "1", "end": "1000",
        "total_cn": "2", "major_cn": "1", "minor_cn": "1",
    }])
    return events, purity, cnv


def test_missing_vaf_event_is_neutral_end_to_end(tmp_path):
    events, purity, cnv = _write_ccf_inputs(tmp_path)
    write_tsv(events, [{
        "event_id": "E_NOVAF", "sample_id": "S1",
        "mutation_source": "SNV", "chrom": "chr1", "pos": "100",
    }])
    row = build_ccf_2(events, purity, cnv, load_profile("default"), tmp_path / "ccf.tsv")[0]
    assert row["ccf_status"] == "unresolved"
    # Before C4-02 this row carried 0.6500 into the ranking.
    assert row["clonality_multiplier"] == "1.0000"


def test_missing_purity_event_is_neutral_end_to_end(tmp_path):
    events, _purity, cnv = _write_ccf_inputs(tmp_path)
    write_tsv(events, [{
        "event_id": "E_NOPUR", "sample_id": "S1", "mutation_source": "SNV",
        "chrom": "chr1", "pos": "100", "tumor_vaf": "0.30",
        "tumor_depth": "100", "tumor_alt_count": "30",
    }])
    row = build_ccf_2(events, None, cnv, load_profile("default"), tmp_path / "ccf.tsv")[0]
    assert row["ccf_status"] == "vaf_only_unresolved"
    assert row["clonality_multiplier"] == "1.0000"


# ---------------------------------------------------------------- C6-01 ----

def test_self_similarity_status_is_a_first_class_column():
    assert "self_similarity_status" in PEPTIDE_FIELDS
    assert (
        PEPTIDE_FIELDS.index("self_similarity_status")
        == PEPTIDE_FIELDS.index("self_similarity_score") + 1
    )


def test_variant_peptide_rows_leave_self_similarity_blank():
    row = peptide_from_variant_row(
        {"mutant_peptide": "AAAXAAAA", "wildtype_peptide": "AAAQAAAA"},
        "S1",
        {"event_id": "E1", "event_type": "SNV", "gene": "GENE"},
        "HLA-A*02:01",
    )
    assert row["self_similarity_score"] == ""
    assert row["self_similarity_status"] == "UNASSESSED"


def test_easyfuse_row_builder_no_longer_fabricates_zero():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        row = easyfuse_row_to_peptide(
            {"neo_peptide_sequence": "ABCDEFGHIJK", "frame": "neo_frame"},
            {"event_id": "EF1", "gene": "GENE1-GENE2"},
            sample_id="S1",
        )
    assert row is not None
    assert row["self_similarity_score"] == ""
    assert row["self_similarity_status"] == "UNASSESSED"


def test_phase1_sv_peptides_leave_self_similarity_blank():
    cand = PeptideCandidate(
        peptide_id="P1", event_id="SV1", sample_id="S1",
        protein_sequence_id="R1", peptide="ACDEFGHIK",
        wildtype_peptide="ACDEFGHIL", hla_allele="HLA-A*02:01",
        mhc_class="I", peptide_length=9, peptide_start_aa=1,
        peptide_end_aa=9, crosses_junction="yes", contains_novel_aa="yes",
        novel_aa_count=1, wildtype_match_status="NO_EXACT_MATCH",
    )
    row = _raw_peptide_from_candidate(cand, "SV_Fusion", "GENE1-GENE2")
    assert row["self_similarity_score"] == ""
    assert row["self_similarity_status"] == "UNASSESSED"


def test_pvactools_parser_maps_self_similarity_honestly():
    event = {"event_id": "E1", "event_type": "SNV", "gene": "GENE", "mutation_source": "SNV"}
    without = peptide_from_row(
        {"MT Epitope Seq": "AAGAIGGVLL", "HLA Allele": "HLA-A*02:01"},
        "S1", event, "pVACseq", source_file="ep.csv", source_row_number=2,
    )
    assert without["self_similarity_score"] == ""
    assert without["self_similarity_status"] == "UNASSESSED"
    with_sim = peptide_from_row(
        {"MT Epitope Seq": "AAGAIGGVLL", "HLA Allele": "HLA-A*02:01", "Self Similarity": "0.91"},
        "S1", event, "pVACseq", source_file="ep.csv", source_row_number=3,
    )
    assert with_sim["self_similarity_score"] == "0.91"
    assert with_sim["self_similarity_status"] == "COMPUTED"


def test_safety_gate_blank_self_similarity_is_not_a_zero_measurement(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "KRAS", "event_confidence": "0.9",
        "event_expression": "10", "driver_relevance": "0.8", "clonality": "0.8",
        "persistence": "0.8", "tumor_specificity": "0.8",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "SNV", "mutation_source": "SNV", "gene": "KRAS",
        "peptide": "VVVGAGGVGA", "wildtype_peptide": "VVVGAGGVGG",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2", "presentation_score": "0.9",
        "immunogenicity_score": "0.7", "wildtype_binding_rank": "0.2",
        "self_similarity_score": "",
        "self_similarity_status": "UNASSESSED",
    }])
    rows, _events = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        profile=load_profile("default"),
    )
    row = rows[0]
    # The audit trail shows "nothing was measured", not a fabricated zero.
    assert row["closest_self_similarity"] == ""
    assert "high_self_similarity_critical_tissue" not in row["safety_reason"]
