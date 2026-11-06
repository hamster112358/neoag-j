"""Tests for the 2026-11-06 audit follow-up round (neoag-audit upload).

Each test maps to one uploaded document:
  purity fix.docx                 -> purity_confidence_route names the route,
                                     and it reaches the ccf_2 rows + reports
  CCF_multiplicity fix.docx       -> multiplicity_sensitivity is display-only
                                     pass-through into the ranked tables
  HLA fix.docx                    -> allele-level consensus as a parallel,
                                     review-only branch (pair table unchanged)
  sensitivity test on ranking.docx-> weight-sensitivity modes baseline/global/
                                     regional (SMAA acceptability, Dirichlet)
"""

from types import SimpleNamespace

from neoag.utils import write_tsv


# --- purity fix: the route that produced the confidence label is guaranteed ---

def test_purity_route_names_the_route_for_every_input_shape(tmp_path):
    from neoag.ccf_v2 import (
        ROUTE_CONSENSUS_MAP, ROUTE_EXPLICIT_COLUMN, ROUTE_NOT_ASSESSED,
        ROUTE_VALUE_THRESHOLDS, load_purity_info,
    )

    info = load_purity_info(None)
    assert info["purity_confidence_route"] == ROUTE_NOT_ASSESSED

    bare = tmp_path / "bare.tsv"
    write_tsv(bare, [{"purity": "0.40", "ploidy": "2"}])
    info = load_purity_info(bare)
    assert info["purity_confidence"] == "high"
    assert info["purity_confidence_route"] == ROUTE_VALUE_THRESHOLDS

    explicit = tmp_path / "explicit.tsv"
    write_tsv(explicit, [{"purity": "0.40", "ploidy": "2", "confidence": "low"}])
    info = load_purity_info(explicit)
    # The explicit column must win AND be named as the route.
    assert info["purity_confidence"] == "low"
    assert info["purity_confidence_route"] == ROUTE_EXPLICIT_COLUMN

    import json
    rec = tmp_path / "recommendation.json"
    rec.write_text(json.dumps({
        "status": "CONCORDANT", "recommended_purity": 0.55,
        "tool_values": {"FACETS": 0.55, "PURPLE": 0.57}, "n_tools": 2,
    }), encoding="utf-8")
    info = load_purity_info(rec)
    assert info["purity_confidence"] == "high"
    assert info["purity_confidence_route"] == ROUTE_CONSENSUS_MAP


def test_purity_route_reaches_every_ccf2_row(tmp_path):
    from neoag.ccf_v2 import ROUTE_EXPLICIT_COLUMN, build_ccf_2
    from neoag.config import load_profile

    events = tmp_path / "events.tsv"
    purity_f = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, [{"event_id": "E1", "sample_id": "S1", "mutation_source": "SNV",
                        "chrom": "chr9", "pos": "150", "tumor_vaf": "0.30",
                        "tumor_depth": "100", "tumor_alt_count": "30"}])
    write_tsv(purity_f, [{"purity": "0.80", "ploidy": "2", "source": "FACETS", "confidence": "high"}])
    write_tsv(cnv, [{"chrom": "chr9", "start": "1", "end": "1000", "total_cn": "2", "major_cn": "2", "minor_cn": "0"}])
    rows = build_ccf_2(events, purity_f, cnv, load_profile("default"), tmp_path / "ccf.tsv")
    assert rows[0]["purity_confidence_route"] == ROUTE_EXPLICIT_COLUMN


def test_reports_show_confidence_route_and_criteria():
    from neoag.reports_dual import _patient_purity_consensus

    bundle = SimpleNamespace(
        purity_consensus={
            "recommended_purity": "0.55", "recommended_ploidy": "2",
            "selected_tool": "多工具共识", "status": "CONCORDANT",
            "basis": "多工具结果已并列保留。",
            "confidence": "high",
            "confidence_route": "consensus_status_map",
            "confidence_criteria": "purity-confidence-v2: ...",
        },
        provenance={},
    )
    result, basis = _patient_purity_consensus(bundle)
    assert "置信度 high" in result
    assert "consensus_status_map" in result
    assert "purity-confidence-v2" in basis


# --- CCF multiplicity fix: display-only pass-through, never a rule ----------

def test_multiplicity_sensitivity_and_purity_fields_ride_into_scored_rows():
    from neoag.scoring import score_peptide

    prof = {"score_weights": {}, "scoring": {}, "gates": {"enabled": False}, "ccf_lite": {}}
    p = {"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I", "peptide_id": "P1"}
    e = {"event_id": "E1", "event_type": "SNV", "event_confidence": "0.9",
         "multiplicity_sensitivity": "material:intermediate..clonal",
         "purity_confidence": "medium",
         "purity_confidence_route": "explicit_confidence_column"}
    pres = {"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01"}
    out = score_peptide(p, e, prof, pres, {})
    assert out["multiplicity_sensitivity"] == "material:intermediate..clonal"
    assert out["purity_confidence_route"] == "explicit_confidence_column"
    # Display only: the same peptide scored without the display columns must
    # get the identical efficacy score (no ranking rule consumes them).
    out2 = score_peptide({"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01",
                          "mhc_class": "I", "peptide_id": "P1"},
                         {"event_id": "E1", "event_type": "SNV", "event_confidence": "0.9"},
                         prof, pres, {})
    assert out["efficacy_score"] == out2["efficacy_score"]


# --- HLA fix: allele level as a parallel branch, pair level unchanged -------

def _three_way_split_rows():
    return [
        {"tool": "OptiType", "locus": "A", "lowres1": "A*02:01", "lowres2": "A*30:01",
         "highres1": "", "highres2": "", "allele1": "A*02:01", "allele2": "A*30:01"},
        {"tool": "SpecHLA", "locus": "A", "lowres1": "A*02:01", "lowres2": "A*02:06",
         "highres1": "", "highres2": "", "allele1": "A*02:01", "allele2": "A*02:06"},
        {"tool": "HLA-LA", "locus": "A", "lowres1": "A*02:06", "lowres2": "A*30:01",
         "highres1": "", "highres2": "", "allele1": "A*02:06", "allele2": "A*30:01"},
    ]


def test_allele_level_consensus_resolves_the_three_way_pair_split():
    from neoag.agent_skills.hla_typing_compare import allele_consensus, consensus

    rows = _three_way_split_rows()
    pair = next(r for r in consensus(rows) if r["locus"] == "A")
    # Pair-level consensus stays exactly as before: three votes, no winner.
    assert pair["status"] == "DISCORDANT"
    assert pair["support"] == "1/3"

    by_allele = {r["allele_lowres"]: r for r in allele_consensus(rows) if r["locus"] == "A"}
    # Allele level sees what pair level cannot: every allele has two-tool
    # support; rank 1 (deterministic tie-break) is selected, the second pick
    # is flagged as an honest tie for the reviewer.
    assert by_allele["A*02:01"]["support"] == "2/3"
    assert by_allele["A*02:01"]["status"] == "ALLELE_CONSENSUS"
    assert by_allele["A*02:01"]["selected"] == "yes"
    assert by_allele["A*02:01"]["voting_tools"] == "OptiType;SpecHLA"
    assert by_allele["A*30:01"]["rank"] == "3"
    assert by_allele["A*30:01"]["selected"] == "no"
    assert by_allele["A*02:01"]["tie"] == "yes"


def test_allele_level_agrees_with_pair_level_when_tools_concur():
    from neoag.agent_skills.hla_typing_compare import allele_consensus

    rows = [
        {"tool": "OptiType", "locus": "A", "lowres1": "A*02:01", "lowres2": "A*30:01"},
        {"tool": "SpecHLA", "locus": "A", "lowres1": "A*02:01", "lowres2": "A*30:01"},
    ]
    alleles = [r for r in allele_consensus(rows) if r["locus"] == "A"]
    assert [(r["allele_lowres"], r["selected"], r["tie"]) for r in alleles] == [
        ("A*02:01", "yes", "no"), ("A*30:01", "yes", "no"),
    ]


# --- sensitivity test on ranking: baseline / global / regional modes --------

def _sens_rows(count=20):
    rows = []
    for i in range(count):
        rows.append({
            "peptide_id": f"P{i}", "peptide": "A" * 9, "hla_allele": "HLA-A*02:01",
            "mhc_class": "I",
            "netmhcpan_ba_rank": f"{0.1 + i * 0.3:.2f}",
            "netmhcpan_el_rank": f"{0.2 + i * 0.25:.2f}",
            "mhcflurry_presentation_score": f"{0.9 - i * 0.03:.3f}",
            "mhcflurry_processing_score": f"{0.7 - i * 0.02:.3f}",
            "netchop_31d_cterm_score": f"{0.8 - i * 0.025:.3f}",
            "netmhcstabpan_rank": f"{1.0 + i * 0.5:.2f}",
        })
    return rows


def test_sensitivity_modes_vector_counts_and_bounds():
    from neoag.presentation import weight_sensitivity_table

    rows = _sens_rows()
    baseline = weight_sensitivity_table(rows, {})
    # 1 + 2W with six active weights = 13 (the audit's own arithmetic).
    assert int(baseline[0]["n_vectors"]) == 13
    assert baseline[0]["mode"] == "baseline"
    assert baseline[0]["verdict"] in {"rank_stable", "rank_minor_shift", "rank_unstable"}

    glob = weight_sensitivity_table(rows, {}, mode="global")
    # Full +-perturb grid: 3^6 vectors (base included exactly once).
    assert int(glob[0]["n_vectors"]) == 3 ** 6

    regional = weight_sensitivity_table(rows, {}, mode="regional", samples=100)
    # Dirichlet draws plus the base vector.
    assert int(regional[0]["n_vectors"]) == 101
    for r in regional:
        for key in ("acceptability_top1", "acceptability_top10",
                    "acceptability_top50", "acceptability_top100"):
            assert 0.0 <= float(r[key]) <= 1.0
        assert 1.0 <= float(r["expected_rank"]) <= len(rows)
    # Acceptability is monotone in K.
    top = regional[0]
    assert float(top["acceptability_top1"]) <= float(top["acceptability_top10"])
    assert float(top["acceptability_top10"]) <= float(top["acceptability_top50"])
    assert float(top["acceptability_top50"]) <= float(top["acceptability_top100"])


def test_regional_mode_is_deterministic_and_concentration_controls_spread():
    from neoag.presentation import weight_sensitivity_table

    rows = _sens_rows()
    a = weight_sensitivity_table(rows, {}, mode="regional", samples=60, concentration=8.0)
    b = weight_sensitivity_table(rows, {}, mode="regional", samples=60, concentration=8.0)
    assert a == b  # fixed seed: same inputs, same table

    tight = weight_sensitivity_table(rows, {}, mode="regional", samples=60, concentration=200.0)
    loose = weight_sensitivity_table(rows, {}, mode="regional", samples=60, concentration=1.0)
    # A tighter Dirichlet region must not produce more rank dispersion than a
    # looser one, measured on the worst candidate.
    assert max(float(r["rank_std"]) for r in tight) <= max(float(r["rank_std"]) for r in loose)


def test_unknown_sensitivity_mode_is_rejected():
    from neoag.presentation import weight_sensitivity_table
    import pytest

    with pytest.raises(ValueError, match="unknown sensitivity mode"):
        weight_sensitivity_table(_sens_rows(3), {}, mode="everything")
