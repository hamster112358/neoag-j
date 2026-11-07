from neoag.config import load_profile
from neoag.peptide_safety_gate import build_peptide_safety_gate
from neoag.safety import load_normal_expression
from neoag.scoring import merge_safety_gate, score
from neoag.utils import write_tsv


def _inputs(tmp_path):
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
        "peptide": "VVVGAGGVG", "wildtype_peptide": "VVVGAGGVG",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2", "presentation_score": "0.9",
        "immunogenicity_score": "0.7", "wildtype_binding_rank": "0.2",
        "self_similarity_score": "0.1",
    }])
    return events, peptides


def test_layered_safety_fields_propagate_into_ranked_peptide():
    gate = {
        "safety_status": "SAFETY_PARTIAL",
        "normal_proteome_exact_match_status": "NOT_DETECTED",
        "normal_transcript_junction_match_status": "NOT_DETECTED",
        "normal_immunopeptidome_match_status": "NOT_DETECTED",
        "similar_peptide_cross_reactivity_status": "UNASSESSED",
        "source_gene_expression_context_status": "ASSESSED",
        "critical_organ_expression_context_status": "DETECTED_AUXILIARY",
        "hematopoietic_expression_context_status": "DETECTED_AUXILIARY",
        "tcr_contact_anchor_context_status": "NOT_APPLICABLE",
        "final_safety_conclusion": "PARTIAL_DIRECT_SAFETY_EVIDENCE",
    }
    row = merge_safety_gate({"peptide_id": "P1"}, gate)
    for field, value in gate.items():
        if field != "safety_status":
            assert row[field] == value


def test_missing_safety_resources_are_partial_not_pass(tmp_path):
    events, peptides = _inputs(tmp_path)
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        profile=load_profile("default"),
    )
    row = rows[0]
    assert row["safety_status"] == "SAFETY_PARTIAL"
    assert row["normal_tissue_max_tpm"] == ""
    assert row["reference_proteome_exact_match"] == "not_assessed"
    assert "normal_expression" in row["safety_missing_layers"]
    assert "reference_proteome" in row["safety_missing_layers"]
    assert event_rows[0]["event_safety_status"] == "SAFETY_PARTIAL"


def test_zero_expression_is_assessed_when_gene_is_present(tmp_path):
    events, peptides = _inputs(tmp_path)
    normal = tmp_path / "normal.tsv"
    write_tsv(normal, [{
        "gene": "KRAS", "normal_tissue_max_tpm": "0",
        "normal_hspc_tpm": "0", "critical_tissue_hit": "no",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_expression=normal,
    )
    assert rows[0]["normal_expression_status"] == "ASSESSED"
    assert rows[0]["normal_hspc_status"] == "ASSESSED"
    assert "normal_expression" not in rows[0]["safety_missing_layers"]


def test_event_caution_propagates_to_ranked_peptide(tmp_path):
    events, peptides = _inputs(tmp_path)
    presentation = tmp_path / "presentation.tsv"
    peptide_safety = tmp_path / "peptide_safety.tsv"
    event_safety = tmp_path / "event_safety.tsv"
    write_tsv(presentation, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "peptide": "VVVGAGGVG", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_evidence_score": "0.9", "presentation_evidence_score": "0.9",
        "presentation_evidence_grade": "A",
    }])
    write_tsv(peptide_safety, [{
        "peptide_id": "P1", "event_id": "E1", "safety_status": "PASS",
        "safety_reason": "no_major_signal", "safety_evidence_completeness": "1.0000",
    }])
    write_tsv(event_safety, [{
        "event_id": "E1", "event_safety_status": "CAUTION",
        "event_safety_reason": "critical_tissue_expression",
        "safety_evidence_completeness": "0.8000",
    }])
    _, ranked = score(
        events, peptides, presentation, None, None, None, None,
        load_profile("default"), tmp_path / "ranked_events.tsv", tmp_path / "ranked_peptides.tsv",
        peptide_safety_tsv=peptide_safety, event_safety_tsv=event_safety,
    )
    assert ranked[0]["safety_status"] == "CAUTION"
    assert "critical_tissue_expression" in ranked[0]["safety_reason"]


def test_peptide_local_proteome_match_does_not_fail_whole_event(tmp_path):
    events, peptides = _inputs(tmp_path)
    ref = tmp_path / "proteome.fa"
    ref.write_text(">normal\nXXVVVGAGGVGXX\n", encoding="utf-8")
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        profile=load_profile("default"),
        reference_proteome=ref,
    )
    assert rows[0]["safety_status"] == "FAIL"
    assert event_rows[0]["event_safety_status"] == "SAFETY_PARTIAL"


def test_indel_crosses_junction_does_not_require_normal_rna_junctions(tmp_path):
    events, peptides = _inputs(tmp_path)
    event_rows = [{
        "event_id": "E1", "sample_id": "S1", "event_type": "InDel",
        "mutation_source": "InDel", "gene": "GENE1",
    }]
    peptide_rows = [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "InDel", "mutation_source": "InDel", "gene": "GENE1",
        "peptide_consequence": "insertion", "crosses_junction": "yes",
        "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }]
    write_tsv(events, event_rows)
    write_tsv(peptides, peptide_rows)
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
    )
    assert rows[0]["normal_junction_assessment_status"] == "NOT_APPLICABLE"
    assert "normal_junction" not in rows[0]["safety_missing_layers"]


def test_ctat_normal_fusion_pair_matches_double_dash_catalog(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    junctions = tmp_path / "normal_fusions.tsv"
    write_tsv(events, [{
        "event_id": "EF1", "sample_id": "S1", "event_type": "Fusion",
        "mutation_source": "SV", "gene": "GENE1::GENE2",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "EF1", "sample_id": "S1",
        "event_type": "Fusion", "mutation_source": "SV", "gene": "GENE1::GENE2",
        "peptide_consequence": "fusion", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(junctions, [{
        "gene_pair": "GENE1--GENE2", "normal_sample_count": "4",
        "source": "GTEx_recurrent_StarF2019", "tissue": "Liver",
        "junction_class": "normal_recurrent_fusion",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        profile=load_profile("default"),
        normal_junctions=junctions,
    )
    assert rows[0]["normal_junction_assessment_status"] == "ASSESSED"
    assert rows[0]["normal_junction_seen"] == "yes"
    assert rows[0]["normal_junction_max_reads"] == "4"
    # Round-3 operator decision (FURTHER_SUGGESTIONS WS8): a recurrent normal
    # fusion of the same gene pair is CONTEXT, not proof that THIS candidate's
    # breakpoint exists in normal tissue - review-level by default, and the
    # exact-detection field stays honest about it.
    assert rows[0]["normal_junction_match_level"] == "GENE_PAIR_ONLY_CONTEXT"
    assert rows[0]["safety_status"] == "CAUTION"
    assert "gene_pair_only_normal_junction_context" in rows[0]["safety_reason"]
    assert rows[0]["normal_transcript_junction_match_status"] == "GENE_PAIR_CONTEXT_NOT_EXACT"

    # The conservative override remains available and is explicit.
    rows_reject, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety_reject.tsv",
        profile={"safety": {"normal_junction_gene_pair_only_action": "reject"}},
        normal_junctions=junctions,
    )
    assert rows_reject[0]["safety_status"] == "FAIL"
    assert "normal_junction_seen_gene_pair_only" in rows_reject[0]["safety_reason"]


def test_fusion_catalog_does_not_mark_splice_junction_assessed(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    junctions = tmp_path / "normal_fusions.tsv"
    write_tsv(events, [{
        "event_id": "S1", "sample_id": "S1", "event_type": "InDel",
        "mutation_source": "InDel", "gene": "GENE1",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "S1", "sample_id": "S1",
        "event_type": "InDel", "mutation_source": "InDel", "gene": "GENE1",
        "peptide_consequence": "splice_junction", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(junctions, [{
        "gene_pair": "OTHER1::OTHER2", "normal_sample_count": "4",
        "junction_class": "normal_recurrent_fusion",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_junctions=junctions,
    )
    assert rows[0]["normal_junction_assessment_status"] == "NOT_APPLICABLE"
    assert "normal_junction" not in rows[0]["safety_missing_layers"]


def test_observed_splice_junction_uses_gtex_coordinate_catalog(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    junctions = tmp_path / "normal_splice_junctions.tsv"
    write_tsv(events, [{
        "event_id": "SJ1", "sample_id": "S1", "event_type": "splice_junction",
        "mutation_source": "splice", "gene": "GENE1",
        "junction_id": "chr1:100-200:+",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "SJ1", "sample_id": "S1",
        "event_type": "splice_junction", "mutation_source": "splice", "gene": "GENE1",
        "junction_id": "chr1:100-200:+", "peptide_consequence": "splice_junction",
        "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(junctions, [{
        "junction_id": "chr1:100-200:+", "normal_sample_count": "1",
        "source": "GTEx_v11_exon_junctions", "junction_class": "normal_splice_junction",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_junctions=junctions,
    )
    assert rows[0]["normal_junction_assessment_status"] == "ASSESSED"
    assert rows[0]["normal_junction_seen"] == "yes"
    assert rows[0]["normal_junction_source"] == "GTEx_v11_exon_junctions"


def test_normal_expression_merges_gene_symbol_and_ensembl_alias_layers(tmp_path):
    normal = tmp_path / "normal.tsv"
    write_tsv(normal, [
        {
            "gene": "GENE1", "ensembl_gene_id": "ENSG000001",
            "normal_tissue_max_tpm": "12", "normal_expression_status": "ASSESSED",
            "normal_hspc_tpm": "0", "normal_hspc_status": "UNASSESSED",
        },
        {
            "gene": "ENSG000001", "ensembl_gene_id": "ENSG000001",
            "normal_tissue_max_tpm": "0", "normal_expression_status": "UNASSESSED",
            "normal_hspc_tpm": "23", "normal_hspc_status": "ASSESSED",
            "normal_hspc_unit": "HPA_nCPM",
        },
    ])
    loaded = load_normal_expression(normal)
    assert loaded["GENE1"]["normal_expression_status"] == "ASSESSED"
    assert loaded["GENE1"]["normal_hspc_status"] == "ASSESSED"
    assert loaded["GENE1"]["normal_hspc_tpm"] == 23
    assert loaded["ENSG000001"]["normal_tissue_max_tpm"] == 12


def test_splice_safety_recovers_exact_ensembl_id_from_source_provenance(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    normal = tmp_path / "normal.tsv"
    write_tsv(events, [{
        "event_id": "SJ|GRCh38|chr10|100|200|.",
        "sample_id": "S1",
        "event_type": "Splice",
        "mutation_source": "splice",
        "gene": "chr10:100",
        "source_records": "ENSG00000156395:E54.1-E62.1;JUNC00210279",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1",
        "event_id": "SJ|GRCh38|chr10|100|200|.",
        "sample_id": "S1",
        "event_type": "Splice",
        "mutation_source": "splice",
        "gene": "chr10:100",
        "peptide_consequence": "splice_junction",
        "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01",
    }])
    write_tsv(normal, [{
        "gene": "SORCS3",
        "ensembl_gene_id": "ENSG00000156395",
        "normal_tissue_max_tpm": "9.88",
        "normal_tissue_max_tissue": "Brain",
        "normal_expression_status": "ASSESSED",
        "normal_hspc_tpm": "1.0",
        "normal_hspc_status": "ASSESSED",
        "normal_hspc_unit": "HPA_nCPM",
    }])

    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        normal_expression=normal,
    )

    assert rows[0]["normal_expression_status"] == "ASSESSED"
    assert rows[0]["normal_hspc_status"] == "ASSESSED"
    assert rows[0]["normal_tissue_max_tpm"] == "9.8800"
    assert rows[0]["normal_hspc_tpm"] == "1.0000"
    assert "normal_expression" not in rows[0]["safety_missing_layers"]
    assert "normal_hspc" not in rows[0]["safety_missing_layers"]
    assert event_rows[0]["normal_expression_status"] == "ASSESSED"
    assert event_rows[0]["normal_hspc_status"] == "ASSESSED"
    assert rows[0]["source_gene_expression_context_status"] == "ASSESSED"
    assert rows[0]["normal_transcript_junction_match_status"] == "UNASSESSED"
    assert rows[0]["similar_peptide_cross_reactivity_status"] == "UNASSESSED"
    assert rows[0]["final_safety_conclusion"] == "PARTIAL_DIRECT_SAFETY_EVIDENCE"
    assert "normal_expression" not in rows[0]["safety_missing_layers"]
    assert "normal_hspc" not in rows[0]["safety_missing_layers"]


def test_missing_gene_in_available_normal_reference_has_specific_reason(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    normal = tmp_path / "normal.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "GENE_NOT_IN_REFERENCE",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "SNV", "mutation_source": "SNV", "gene": "GENE_NOT_IN_REFERENCE",
        "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01",
    }])
    write_tsv(normal, [{
        "gene": "OTHER_GENE", "normal_tissue_max_tpm": "1",
        "normal_expression_status": "ASSESSED", "normal_hspc_tpm": "0",
        "normal_hspc_status": "ASSESSED",
    }])

    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_expression=normal,
    )

    assert "normal_expression_gene_not_in_reference" in rows[0]["safety_reason"]
    assert "normal_hspc_gene_not_in_reference" in rows[0]["safety_reason"]


def _fusion_inputs(tmp_path):
    events = tmp_path / "f_events.tsv"
    peptides = tmp_path / "f_peptides.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "Fusion",
        "mutation_source": "SV", "gene": "GENE1::GENE2", "event_confidence": "0.9",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "Fusion", "mutation_source": "SV", "gene": "GENE1::GENE2",
        "peptide": "ABCDEFGHIJK", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2",
    }])
    return events, peptides


def test_c6_03_gene_pair_only_normal_junction_is_labelled_as_context(tmp_path):
    # Round-2 (HANDOFF section I, C6-03): a normal junction found under the
    # candidate's gene pair alone does not prove THIS candidate's breakpoint
    # occurs in normal tissue. The reject action may stay (conservative), but
    # the match level must say GENE_PAIR_ONLY_CONTEXT everywhere it is read.
    events, peptides = _fusion_inputs(tmp_path)
    junctions = tmp_path / "normal_junctions.tsv"
    write_tsv(junctions, [{"gene1": "GENE1", "gene2": "GENE2", "normal_reads": "9"}])
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        profile=load_profile("default"),
        normal_junctions=junctions,
    )
    # Round-3 (FURTHER_SUGGESTIONS WS8): review-level context is the DEFAULT.
    assert rows[0]["safety_status"] == "CAUTION"
    assert "gene_pair_only_normal_junction_context" in rows[0]["safety_reason"]
    assert rows[0]["normal_junction_match_level"] == "GENE_PAIR_ONLY_CONTEXT"
    assert event_rows[0]["normal_junction_match_level"] == "GENE_PAIR_ONLY_CONTEXT"
    assert event_rows[0]["event_safety_status"] != "FAIL"
    # Gene-pair context must never read as exact junction detection.
    assert rows[0]["normal_transcript_junction_match_status"] == "GENE_PAIR_CONTEXT_NOT_EXACT"

    # The explicit conservative override restores the old reject policy.
    profile = dict(load_profile("default"))
    safety = dict(profile.get("safety", {}))
    safety["normal_junction_gene_pair_only_action"] = "reject"
    profile["safety"] = safety
    rows2, event_rows2 = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety2.tsv",
        out_event_safety=tmp_path / "event_safety2.tsv",
        profile=profile,
        normal_junctions=junctions,
    )
    assert rows2[0]["safety_status"] == "FAIL"
    assert "normal_junction_seen_gene_pair_only" in rows2[0]["safety_reason"]
    assert event_rows2[0]["event_safety_status"] == "FAIL"


def _coord_junction_case(tmp_path, name, cand_strand, normal_strand):
    events = tmp_path / f"x_events_{name}.tsv"
    peptides = tmp_path / f"x_peptides_{name}.tsv"
    junctions = tmp_path / f"x_junctions_{name}.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "Fusion",
        "mutation_source": "SV", "gene": "GENE1::GENE2",
        "chrom": "chr1", "junction_start": "100", "junction_end": "200",
        **({"strand": cand_strand} if cand_strand else {}),
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "Fusion", "mutation_source": "SV", "gene": "GENE1::GENE2",
        "peptide": "ABCDEFGHIJK", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2",
    }])
    write_tsv(junctions, [{
        "chromosome": "chr1", "junction_start": "100", "junction_end": "200",
        "normal_reads": "9",
        **({"strand": normal_strand} if normal_strand else {}),
    }])
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / f"x_ps_{name}.tsv",
        out_event_safety=tmp_path / f"x_es_{name}.tsv",
        profile=load_profile("default"),
        normal_junctions=junctions,
    )
    return rows[0], event_rows[0]


def test_c6_03_exact_coordinate_match_keeps_its_own_level(tmp_path):
    # The honest split only matters if the exact match keeps its own label:
    # with matching coordinates AND a strand-aware match the level is
    # EXACT_NORMAL_JUNCTION_MATCH and the plain reject reason still applies.
    row, event_row = _coord_junction_case(tmp_path, "exact", "+", "+")
    assert row["safety_status"] == "FAIL"
    assert row["normal_junction_match_level"] == "EXACT_NORMAL_JUNCTION_MATCH"
    assert "normal_junction_seen" in row["safety_reason"]
    assert row["normal_transcript_junction_match_status"] == "DETECTED"
    assert event_row["normal_junction_match_level"] == "EXACT_NORMAL_JUNCTION_MATCH"

    # Review B3b-08: equal coordinates with OPPOSITE resolved strands is not
    # this junction at all - no match, not even a possible one.
    row_opp, _ = _coord_junction_case(tmp_path, "opposite", "+", "-")
    assert row_opp["normal_junction_match_level"] == ""
    assert row_opp["normal_junction_seen"] == "no"
    assert row_opp["normal_transcript_junction_match_status"] == "NOT_DETECTED"
    assert row_opp["safety_status"] != "FAIL"

    # One side with an unresolved strand: explicit ambiguous review state -
    # neither an exact positive (reject) nor a definitive negative.
    row_amb, event_amb = _coord_junction_case(tmp_path, "amb", "+", "")
    assert row_amb["normal_junction_match_level"] == "POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS"
    assert "possible_normal_junction_strand_ambiguous" in row_amb["safety_reason"]
    assert row_amb["safety_status"] == "CAUTION"
    assert row_amb["normal_transcript_junction_match_status"] == "POSSIBLE_STRAND_AMBIGUOUS"
    assert row_amb["final_safety_conclusion"] == "REVIEW_DIRECT_SAFETY_EVIDENCE"
    assert event_amb["event_safety_status"] == "CAUTION"
    # The candidate side being unstranded is symmetric.
    row_amb2, _ = _coord_junction_case(tmp_path, "amb2", "", "+")
    assert row_amb2["normal_junction_match_level"] == "POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS"
    # Both sides unstranded never earns "verified" either.
    row_both, _ = _coord_junction_case(tmp_path, "both", "", "")
    assert row_both["normal_junction_match_level"] == "POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS"
    assert row_both["safety_status"] == "CAUTION"


# --------------------- C6-04 / C6-05 round-3 follow-up (FURTHER_SUGGESTIONS WS5) --

def _normal_variant_event(tmp_path, name, **extra):
    events = tmp_path / f"events_{name}.tsv"
    peptides = tmp_path / f"peptides_{name}.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "KRAS", "event_confidence": "0.9",
        "event_expression": "10", "driver_relevance": "0.8", "clonality": "0.8",
        "persistence": "0.8", "tumor_specificity": "0.8", **extra,
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "SNV", "mutation_source": "SNV", "gene": "KRAS",
        "peptide": "VVVGAGGVG", "wildtype_peptide": "VVVGAGGVG",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2", "presentation_score": "0.9",
        "immunogenicity_score": "0.7", "wildtype_binding_rank": "0.2",
        "self_similarity_score": "0.1",
    }])
    return events, peptides


def _gate(tmp_path, name, profile=None, **extra):
    events, peptides = _normal_variant_event(tmp_path, name, **extra)
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / f"ps_{name}.tsv",
        out_event_safety=tmp_path / f"es_{name}.tsv",
        profile=profile,
    )
    return rows[0], event_rows[0]


def test_c6_04_missing_vaf_with_alt_reads_goes_to_review(tmp_path):
    # ALT-positive normal record with NO VAF: not a measured normal-negative
    # and not a confirmed normal-positive - a validation review.
    row, event_row = _gate(tmp_path, "review", normal_alt_count="5")
    assert row["matched_normal_status"] == "NORMAL_SUPPORT_REVIEW"
    assert "matched_normal_support_review" in row["safety_reason"]
    assert row["safety_status"] != "FAIL"
    assert row["review_required"] == "yes"
    assert row["final_safety_conclusion"] == "REVIEW_DIRECT_SAFETY_EVIDENCE"
    assert row["safety_decision_basis"] == "DIRECT_ASSESSED_EVIDENCE"
    assert event_row["event_safety_status"] == "CAUTION"
    assert "matched_normal_support_review" in event_row["event_safety_reason"]


def test_c6_04_missing_vaf_is_recalculable_from_same_allele_counts(tmp_path):
    # Counts and depth from the same allele allow an EXPLICIT recalculation;
    # the basis names the inputs so the recalculation is auditable.
    row, _ = _gate(tmp_path, "recalc", normal_alt_count="4", normal_depth="50")
    assert row["matched_normal_status"] == "normal_alt_support"
    assert row["normal_vaf"] == "0.0800"
    assert row["normal_support_evidence_basis"] == "VAF_RECALCULATED_ALT4_DEPTH50"
    assert row["safety_status"] == "FAIL"
    assert "matched_normal_support" in row["safety_reason"]


def test_c6_04_measured_zero_vaf_stays_distinct_from_missing(tmp_path):
    # Measured zero (0 ALT at real depth) is an ASSESSED normal-negative:
    # no review is raised for it.
    row, _ = _gate(tmp_path, "zero", normal_alt_count="0", normal_depth="30", normal_vaf="0.0")
    assert row["matched_normal_status"] == "no_clear_normal_support"
    assert "matched_normal_support_review" not in row["safety_reason"]
    assert row["normal_support_evidence_basis"] == "REPORTED"


def test_c6_04_contradictory_counts_and_vaf_warn_instead_of_deciding(tmp_path):
    # Reported VAF 0.5 with zero normal ALT reads: integrity warning, review,
    # and NO silent resolution into either a pass or a rejection.
    row, _ = _gate(
        tmp_path, "contradiction",
        normal_alt_count="0", normal_depth="100", normal_vaf="0.5000",
    )
    assert row["matched_normal_status"] == "NORMAL_VAF_CONTRADICTION"
    assert "normal_vaf_contradiction" in row["normal_support_integrity"]
    assert "normal_vaf_integrity_conflict" in row["safety_reason"]
    assert row["safety_status"] != "FAIL"
    assert row["review_required"] == "yes"


def test_c6_04_alt_positive_with_subthreshold_vaf_is_review_not_positive(tmp_path):
    # The old rule called ALT>=3 with VAF exactly 0.0 normal-POSITIVE.  A
    # read count and a VAF that disagree is ambiguity: review, never reject.
    row, _ = _gate(
        tmp_path, "subthresh",
        normal_alt_count="4", normal_depth="50", normal_vaf="0.0000",
    )
    assert row["matched_normal_status"] == "NORMAL_SUPPORT_REVIEW"
    assert row["safety_status"] != "FAIL"


def test_c6_04_inadequate_normal_coverage_is_not_a_measured_negative(tmp_path):
    row, _ = _gate(tmp_path, "thin", normal_alt_count="0", normal_depth="4", normal_vaf="0.0")
    assert row["matched_normal_status"] == "NORMAL_COVERAGE_INADEQUATE"
    assert row["matched_normal_status"] != "no_clear_normal_support"


def _anchor_gate(tmp_path, name, profile):
    events = tmp_path / f"events_{name}.tsv"
    peptides = tmp_path / f"peptides_{name}.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "KRAS",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "SNV", "mutation_source": "SNV", "gene": "KRAS",
        # Single difference at position 2 (an approximate MHC anchor slot)
        # while the wildtype peptide binds strongly.
        "peptide": "VAGVGAGGVG", "wildtype_peptide": "VVGVGAGGVG",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2",
        "wildtype_binding_rank": "0.0001",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / f"ps_{name}.tsv",
        profile=profile,
    )
    return rows[0]


def test_c6_05_anchor_heuristic_never_rejects_without_explicit_policy(tmp_path):
    row = _anchor_gate(tmp_path, "anchor_default", None)
    assert row["mutation_anchor_only"] == "yes"
    assert "anchor_only_mutation_with_wt_binding" in row["safety_reason"]
    # Caution only - the approximate cross-allele rule alone never rejects.
    assert row["safety_status"] == "CAUTION"
    assert row["safety_decision_basis"] == "HEURISTIC_ONLY"

    # Purely documentary mode keeps the evidence columns but raises no signal.
    row_record = _anchor_gate(
        tmp_path, "anchor_record", {"safety": {"anchor_only_action": "record_only"}},
    )
    assert row_record["safety_status"] != "CAUTION"
    assert row_record["mutation_anchor_only"] == "yes"
    assert row_record["safety_decision_basis"] == "HEURISTIC_ONLY"

    # Hard rejection is available ONLY as an explicit, separately justified
    # profile opt-in, and labels itself as such.
    row_reject = _anchor_gate(
        tmp_path, "anchor_reject", {"safety": {"anchor_only_action": "reject"}},
    )
    assert row_reject["safety_status"] == "FAIL"
    assert "anchor_only_reject_explicitly_enabled" in row_reject["safety_reason"]
