from neoag.candidate_identity import candidate_identity, normalize_hla
from neoag.reports_dual import _patient_representatives


def test_four_level_identity_keeps_events_but_collapses_same_peptide_hla():
    # The gene field carries the fusion pair: since the display-key redesign
    # (commits 7474afc/4a22add) the patient view collapses same-gene-pair
    # breakpoint hypotheses via FUSION_DISPLAY|GENE::PAIR.
    first = {
        "event_id": "FUSION_EWSR1_WT1_BP1",
        "protein_change": "ORF1",
        "peptide": "SSYGQQSEK",
        "hla_allele": "HLA-A*11:02",
        "event_type": "Fusion",
        "gene": "EWSR1::WT1",
    }
    second = {
        "event_id": "FUSION_EWSR1_WT1_BP2",
        "protein_change": "ORF2",
        "peptide": "ssygqqsek",
        "hla_allele": "A1102",
        "event_type": "Fusion",
        "gene": "EWSR1::WT1",
    }
    first_ids = candidate_identity(first)
    second_ids = candidate_identity(second)
    assert first_ids["event_identity_id"] != second_ids["event_identity_id"]
    assert first_ids["protein_change_identity_id"] != second_ids["protein_change_identity_id"]
    assert first_ids["peptide_sequence_id"] == second_ids["peptide_sequence_id"]
    assert first_ids["peptide_hla_id"] == second_ids["peptide_hla_id"]
    representatives = _patient_representatives([first, second], 10)
    assert len(representatives) == 1
    # The collapsed representative keeps both breakpoint hypotheses.
    assert representatives[0]["patient_display_hypothesis_count"] == "2"
    # Review D4-04: a caller-supplied event id is an exact primary identity.
    assert first_ids["event_identity_status"] == "EXACT_PRIMARY_ID"


def test_identity_status_names_the_derivation_layer():
    # Review D4-04: the status field must tell consumers, without prefix
    # sniffing, whether an event identity is exact, coordinate-derived, or a
    # row-digest fallback that is NOT a biological identity.
    exact = candidate_identity({"event_id": "E1"})
    assert exact["event_identity_status"] == "EXACT_PRIMARY_ID"
    coordinate = candidate_identity({"event_type": "snv", "chrom": "1", "pos": "100", "ref": "A", "alt": "T"})
    assert coordinate["event_identity_status"] == "COORDINATE_DERIVED_ID"
    unknown = candidate_identity({"peptide": "AAAA", "notes": "no id anywhere"})
    assert unknown["event_identity_status"] == "UNRESOLVED_ROW_FALLBACK"
    other_unknown = candidate_identity({"peptide": "MMMM", "notes": "different row"})
    assert unknown["event_identity_id"] != other_unknown["event_identity_id"]


def test_hla_normalization_is_stable():
    assert normalize_hla("A1102") == "HLA-A*11:02"
    assert normalize_hla("HLA-A*11:02") == "HLA-A*11:02"
