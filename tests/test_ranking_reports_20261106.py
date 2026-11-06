"""Two separate ranking reports (operator request, 2026-11-06).

The pipeline's consensus-based logic ranking and weight-based ranking must be
producible as two self-contained HTML reports from an existing case directory
(no pipeline change, no ranking change). These tests build a mini case dir and
verify both reports, their cross-references, derived columns, and the honest
behaviour when an artifact is missing.
"""
from __future__ import annotations

import json

import pytest

from neoag.ranking_reports import (
    CONSENSUS_REPORT_NAME, WEIGHT_REPORT_NAME, find_scoring_dir, write_ranking_reports)


def _write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    path.write_text("\n".join(
        ["\t".join(fields)] + ["\t".join(str(r.get(f, "")) for f in fields) for r in rows]
    ) + "\n", encoding="utf-8")


def _make_case(tmp_path):
    scoring = tmp_path / "final" / "scoring"
    scoring.mkdir(parents=True)
    ranked = [
        {"peptide_id": "P1", "peptide": "PEPTIDEAAAA", "gene": "G1", "hla_allele": "HLA-A*02:01",
         "mhc_class": "I", "efficacy_score": "0.4800", "immunology_composite_score": "0.6000",
         "l3_hla_binding_score": "0.7000", "l3_hla_presentation_score": "0.5000",
         "safety_multiplier": "1.0000", "presentation_gate_multiplier": "1.0000",
         "appm_multiplier": "0.8000", "ccf_multiplier": "1.0000", "escape_multiplier": "1.0000",
         "cross_platform_multiplier": "", "safety_status": "PASS", "priority": "A"},
        {"peptide_id": "P2", "peptide": "PEPTIDEBBBB", "gene": "G1", "hla_allele": "HLA-A*02:01",
         "mhc_class": "I", "efficacy_score": "0.2400", "immunology_composite_score": "0.6000",
         "l3_hla_binding_score": "0.5000", "l3_hla_presentation_score": "0.7000",
         "safety_multiplier": "1.0000", "presentation_gate_multiplier": "0.5000",
         "appm_multiplier": "0.8000", "ccf_multiplier": "1.0000", "escape_multiplier": "1.0000",
         "cross_platform_multiplier": "", "safety_status": "PASS", "priority": "B"},
        {"peptide_id": "P3", "peptide": "PEPTIDEBBBB", "gene": "G2", "hla_allele": "HLA-B*07:02",
         "mhc_class": "I", "efficacy_score": "0.0000", "immunology_composite_score": "0.6000",
         "l3_hla_binding_score": "0.9000", "l3_hla_presentation_score": "0.9000",
         "safety_multiplier": "0.0000", "presentation_gate_multiplier": "1.0000",
         "appm_multiplier": "0.8000", "ccf_multiplier": "1.0000", "escape_multiplier": "1.0000",
         "cross_platform_multiplier": "", "safety_status": "SAFETY_REJECT", "priority": "D"},
    ]
    _write_tsv(scoring / "ranked_peptides.tsv", ranked)
    consensus = [
        {"peptide_id": "P3", "evidence_rank": "1", "track_rank": "1",
         "peptide": "PEPTIDEBBBB", "gene": "G2", "hla_allele": "HLA-B*07:02",
         "evidence_grade": "R2", "evidence_grade_uncapped": "R1",
         "evidence_consensus_score": "0.90", "consensus_action": "ADVANCE",
         "source_chain_confidence_tier": "C2", "legacy_weighted_rank": "3",
         "event_authenticity_grade": "E2", "rna_support_grade": "N1",
         "presentation_consensus_grade": "P2", "mutant_specificity_grade": "S1",
         "safety_grade": "F2", "hla_appm_grade": "H1", "clonality_grade": "L1",
         "evidence_completeness_grade": "K2"},
        {"peptide_id": "P1", "evidence_rank": "2", "track_rank": "2",
         "peptide": "PEPTIDEAAAA", "gene": "G1", "hla_allele": "HLA-A*02:01",
         "evidence_grade": "R2", "evidence_grade_uncapped": "R2",
         "evidence_consensus_score": "0.80", "consensus_action": "ADVANCE",
         "source_chain_confidence_tier": "C2", "legacy_weighted_rank": "1",
         "event_authenticity_grade": "E2", "rna_support_grade": "N1",
         "presentation_consensus_grade": "P2", "mutant_specificity_grade": "S1",
         "safety_grade": "F1", "hla_appm_grade": "H1", "clonality_grade": "L1",
         "evidence_completeness_grade": "K2"},
        {"peptide_id": "P2", "evidence_rank": "3", "track_rank": "3",
         "peptide": "PEPTIDEBBBB", "gene": "G1", "hla_allele": "HLA-A*02:01",
         "evidence_grade": "R4", "evidence_grade_uncapped": "R4",
         "evidence_consensus_score": "0.30", "consensus_action": "MANUAL_REVIEW",
         "source_chain_confidence_tier": "C3", "legacy_weighted_rank": "2",
         "event_authenticity_grade": "E3", "rna_support_grade": "N2",
         "presentation_consensus_grade": "P3", "mutant_specificity_grade": "S2",
         "safety_grade": "F1", "hla_appm_grade": "H1", "clonality_grade": "L1",
         "evidence_completeness_grade": "K2"},
    ]
    _write_tsv(scoring / "ranked_peptides.evidence_consensus.tsv", consensus)
    _write_tsv(scoring / "weighted_vs_consensus_comparison.tsv", [
        {"peptide_id": "P1", "rank_shift_weighted_minus_consensus": "-1",
         "rank_shift_direction": "decreased"},
        {"peptide_id": "P2", "rank_shift_weighted_minus_consensus": "-1",
         "rank_shift_direction": "decreased"},
        {"peptide_id": "P3", "rank_shift_weighted_minus_consensus": "2",
         "rank_shift_direction": "improved"},
    ])
    _write_tsv(scoring / "ranked_events.evidence_consensus.tsv", [
        {"event_id": "EV1", "gene": "G1", "event_type": "SNV", "evidence_track": "CODING",
         "best_peptide_id": "P1", "best_peptide": "PEPTIDEAAAA",
         "best_hla_allele": "HLA-A*02:01", "best_peptide_evidence_rank": "2",
         "best_evidence_grade": "R2", "source_chain_confidence_tier": "C2",
         "source_chain_confidence_label": "ORTHOGONAL"}])
    (scoring / "evidence_consensus_run.json").write_text(json.dumps({
        "algorithm": "discrete_state_grade_track_pareto_source_chain_v3",
        "status": "PROVISIONAL_RESEARCH_ONLY",
        "rules_name": "test_rules", "rules_version": "0.1"}), encoding="utf-8")
    return scoring, ranked


def test_two_ranking_reports_written_with_cross_references(tmp_path):
    _scoring, _ranked = _make_case(tmp_path)
    paths = write_ranking_reports(tmp_path)
    assert paths["consensus"].name == CONSENSUS_REPORT_NAME
    assert paths["weight"].name == WEIGHT_REPORT_NAME
    consensus_html = paths["consensus"].read_text(encoding="utf-8")
    weight_html = paths["weight"].read_text(encoding="utf-8")

    # Each report presents exactly one ranking basis and points at the other.
    assert "Consensus-Based Logic Ranking" in consensus_html
    assert WEIGHT_REPORT_NAME in consensus_html
    assert "Weight-Based Ranking" in weight_html
    assert CONSENSUS_REPORT_NAME in weight_html
    # Consensus report: rules provenance and its provisional status are shown.
    assert "PROVISIONAL_RESEARCH_ONLY" in consensus_html
    assert "discrete_state_grade_track_pareto_source_chain_v3" in consensus_html

    # Consensus ordering is evidence_rank order (P3 first), NOT the weight
    # order. Match the leading rank cells of each table row (bare <td>P2</td>
    # would also occur as a synthetic layer-grade value).
    cells = ["<td>1</td><td>1</td><td>P3</td>",
             "<td>2</td><td>2</td><td>P1</td>",
             "<td>3</td><td>3</td><td>P2</td>"]
    positions = [consensus_html.index(c) for c in cells]
    assert positions == sorted(positions)
    # Layer grades and the rank-shift context ride along.
    assert "rank_shift_weighted_minus_consensus" in consensus_html
    assert "Consensus event ranking" in consensus_html


def test_weight_report_orders_by_efficacy_and_shows_chains_and_multipliers(tmp_path):
    _scoring, _ranked = _make_case(tmp_path)
    paths = write_ranking_reports(tmp_path)
    html = paths["weight"].read_text(encoding="utf-8")
    # Weight order by efficacy_score: P1 (0.48) then P2 (0.24) then P3 (0.00).
    positions = [html.index(f"<td>{pid}</td>") for pid in ("P1", "P2", "P3")]
    assert positions == sorted(positions)
    # Chain columns come from the profile's Layer-3 dimensions (derived, not
    # hardcoded): only the ones actually present in the table are listed.
    assert "l3_hla_binding_score" in html
    assert "l3_hla_presentation_score" in html
    # Multiplier columns are shown, plus the honest note that they are
    # decisions, not weights.
    assert "safety_multiplier" in html
    assert "presentation_gate_multiplier" in html
    assert "non-compensatory" in html
    # No weight_sensitivity.tsv yet -> the gap is stated, not skipped.
    assert "weight_sensitivity.tsv" in html


def test_weight_report_embeds_the_sensitivity_section(tmp_path):
    from neoag.sensitivity import write_weight_sensitivity
    scoring, ranked = _make_case(tmp_path)
    write_weight_sensitivity(ranked, {}, scoring / "weight_sensitivity.tsv", samples=40)
    paths = write_ranking_reports(tmp_path)
    html = paths["weight"].read_text(encoding="utf-8")
    assert "Weight Sensitivity (regional Dirichlet)" in html
    assert "chain_multiplier" in html
    assert "Regional draws" in html  # the verbatim formula is included
    assert "Evidence-chain contribution" in html


def test_missing_consensus_artifact_is_stated_not_faked(tmp_path):
    scoring = tmp_path / "final" / "scoring"
    scoring.mkdir(parents=True)
    _write_tsv(scoring / "ranked_peptides.tsv", [
        {"peptide_id": "P1", "peptide": "AAAAAAA", "efficacy_score": "0.1",
         "immunology_composite_score": "0.1"}])
    paths = write_ranking_reports(tmp_path)
    consensus_html = paths["consensus"].read_text(encoding="utf-8")
    assert "did not produce" in consensus_html


def test_find_scoring_dir_accepts_nested_and_bare_layouts(tmp_path):
    nested = tmp_path / "case" / "final" / "scoring"
    nested.mkdir(parents=True)
    (nested / "ranked_peptides.tsv").write_text("peptide_id\nP1\n", encoding="utf-8")
    assert find_scoring_dir(tmp_path / "case") == nested
    assert find_scoring_dir(nested.parent) == nested  # final/ passed directly
    assert find_scoring_dir(nested) == nested  # scoring dir passed directly
    with pytest.raises(FileNotFoundError):
        find_scoring_dir(tmp_path / "empty")


def test_cli_exposes_ranking_reports_command(tmp_path):
    import contextlib
    import io

    from neoag.cli import build_parser
    _scoring, _ranked = _make_case(tmp_path)
    buf = io.StringIO()
    with contextlib.suppress(SystemExit), contextlib.redirect_stdout(buf):
        args = build_parser().parse_args(["ranking-reports", "--case-dir", str(tmp_path)])
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        args.func(args)
    printed = out.getvalue()
    assert CONSENSUS_REPORT_NAME in printed and WEIGHT_REPORT_NAME in printed
    help_buf = io.StringIO()
    with contextlib.suppress(SystemExit), contextlib.redirect_stdout(help_buf):
        build_parser().parse_args(["ranking-reports", "--help"])
    assert "--out-dir" in help_buf.getvalue() and "--top" in help_buf.getvalue()
