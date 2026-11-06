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


# --- 2026-11-06 operator round: sensitivity is ONE local mode; component
# count comes from the profile/spec, perturbation is wider, single-driver
# candidates are penalized, and the table feeds a two-part selection.


def test_local_vectors_track_active_components_not_a_hardcoded_count():
    from neoag.presentation import weight_sensitivity_table

    rows = _sens_rows(5)
    default = weight_sensitivity_table(rows, {})
    # Default profile activates six weights -> 1 + 2*6 vectors, and the table
    # SAYS the component count instead of the reader inferring it.
    assert default[0]["n_components"] == "6"
    assert default[0]["n_vectors"] == "13"

    thinned = weight_sensitivity_table(rows, {"presentation_weights": {
        "netmhcpan_ba": 0.25, "netmhcpan_el": 0.0,
        "mhcflurry_presentation": 0.30, "mhcflurry_processing": 0.0,
        "netchop_processing": 0.10, "netmhcstabpan_stability": 0.0,
    }})
    # Zeroing three weights must shrink the analysis to 1 + 2*3 = 7 vectors.
    assert thinned[0]["n_components"] == "3"
    assert thinned[0]["n_vectors"] == "7"


def test_default_perturbation_is_the_widened_range():
    from neoag.presentation import weight_sensitivity_table, SENSITIVITY_DEFAULT_PERTURB
    import inspect

    assert SENSITIVITY_DEFAULT_PERTURB == 0.5
    sig = inspect.signature(weight_sensitivity_table)
    assert sig.parameters["perturb"].default == 0.5
    # The removed mode machinery must be gone, not merely hidden.
    assert "mode" not in sig.parameters and "samples" not in sig.parameters
    row = weight_sensitivity_table(_sens_rows(3), {})[0]
    assert row["perturb"] == "0.50"


def test_single_driver_candidate_is_penalized_by_the_dominance_factor():
    from neoag.presentation import weight_sensitivity_table

    shared = {"peptide": "A" * 9, "hla_allele": "HLA-A*02:01", "mhc_class": "I",
              "mhcflurry_presentation_score": "", "mhcflurry_processing_score": "",
              "netchop_31d_cterm_score": "", "netmhcstabpan_rank": ""}
    balanced = {"peptide_id": "BAL", **shared,
                "netmhcpan_ba_rank": "0.0", "netmhcpan_el_rank": "0.0"}
    single = {"peptide_id": "ONE", **shared,
              "netmhcpan_ba_rank": "", "netmhcpan_el_rank": "",
              "netchop_31d_cterm_score": "1.0"}
    # Both score 1.0 under the base weights...
    out = {r["peptide_id"]: r for r in weight_sensitivity_table([balanced, single], {}, top_k=2)}
    assert out["BAL"]["base_score"] == "1.0000" and out["ONE"]["base_score"] == "1.0000"
    # ...but the single-driver candidate carries half the adjusted score.
    assert float(out["BAL"]["effective_components"]) > 1.9
    assert out["ONE"]["effective_components"] == "1.00"
    assert float(out["BAL"]["dominance_factor"]) > 0.95
    assert out["ONE"]["dominance_factor"] == "0.5000"
    assert int(out["BAL"]["adjusted_rank"]) < int(out["ONE"]["adjusted_rank"])
    assert out["ONE"]["adjusted_score"] == "0.5000"


def test_consistency_and_top_selection_mark_consistent_top_candidates():
    from neoag.presentation import weight_sensitivity_table

    base = {"peptide": "A" * 9, "hla_allele": "HLA-A*02:01", "mhc_class": "I",
            "mhcflurry_presentation_score": "", "mhcflurry_processing_score": "",
            "netchop_31d_cterm_score": "", "netmhcstabpan_rank": ""}
    r1 = {"peptide_id": "R1", **base, "netmhcpan_ba_rank": "0.0", "netmhcpan_el_rank": "9"}
    r2 = {"peptide_id": "R2", **base, "netmhcpan_ba_rank": "9", "netmhcpan_el_rank": "0.0"}
    # R2 leads on the default weights (el 0.35 > ba 0.25) but a +-50% ba-up or
    # el-down vector flips it out of the top 1: consistency must be < 1 and
    # top_selection must say no even for the leading candidate.
    out = {r["peptide_id"]: r for r in weight_sensitivity_table([r1, r2], {}, top_k=1)}
    assert out["R2"]["base_rank"] == "1"
    assert float(out["R2"]["consistency"]) < 1.0 and float(out["R2"]["consistency"]) > 0.8
    assert float(out["R1"]["consistency"]) < 0.2
    assert out["R1"]["top_selection"] == "no" and out["R2"]["top_selection"] == "no"
    # A candidate that can never fall out of the window (only row) is consistent.
    solo = weight_sensitivity_table([r2], {}, top_k=1)[0]
    assert solo["consistency"] == "1.0000" and solo["top_selection"] == "yes"


def test_technical_report_shows_both_tables_and_the_formula(tmp_path):
    from neoag.presentation import weight_sensitivity_table
    from neoag.reports_dual import ReportBundle, make_technical_report

    sens = weight_sensitivity_table(_sens_rows(5), {})
    bundle = ReportBundle(profile={}, events=[], peptides=[], weight_sensitivity=sens)
    make_technical_report(tmp_path / "tech.html", bundle)
    html = (tmp_path / "tech.html").read_text(encoding="utf-8")
    assert "Weight Sensitivity" in html
    assert "Table 1: ranking under the default weights" in html
    assert "Table 2: ranking after the sensitivity test" in html
    assert "adjusted_score = base_score" in html
    assert "dominance_factor = min(1, effective_components/2)" in html
    assert "robustness to the weight choice, not accuracy" in html


def test_invalid_sensitivity_arguments_are_rejected():
    from neoag.presentation import weight_sensitivity_table
    import pytest

    with pytest.raises(ValueError, match="perturb"):
        weight_sensitivity_table(_sens_rows(3), {}, perturb=1.0)
    with pytest.raises(ValueError, match="top_k"):
        weight_sensitivity_table(_sens_rows(3), {}, top_k=0)


# --- APPM integration note / scoring investigation 6.2: binding aggregation
# becomes a profile knob, and predictor disagreement is exported. The default
# must stay the historical max() so existing runs reproduce byte-identically.

def _binding_evidence(tmp_path, agg=None):
    from neoag.config import load_profile
    from neoag.presentation import build_presentation_evidence
    from neoag.utils import read_tsv

    hla = "HLA-A*02:01"
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [{
        "peptide_id": "P1", "event_id": "E1", "peptide": "AAGIGILTV",
        "hla_allele": hla, "mhc_class": "I",
    }])
    net = tmp_path / "net.tsv"
    write_tsv(net, [{"peptide": "AAGIGILTV", "hla_allele": hla,
                     "netmhcpan_ba_rank": "0.6", "netmhcpan_el_rank": "0.8"}])
    mhc = tmp_path / "mhc.tsv"
    write_tsv(mhc, [{"peptide": "AAGIGILTV", "hla_allele": hla,
                     "mhcflurry_affinity_percentile": "3.0",
                     "mhcflurry_processing_score": "0.7",
                     "mhcflurry_presentation_score": "0.8"}])
    if agg is None:
        profile = load_profile("default")
    else:
        path = tmp_path / "profile.toml"
        path.write_text(
            "[presentation_weights]\n"
            "netmhcpan_ba = 0.25\nnetmhcpan_el = 0.35\n"
            "mhcflurry_presentation = 0.30\nmhcflurry_processing = 0.10\n"
            "netchop_processing = 0.10\nnetmhcstabpan_stability = 0.05\n"
            f'binding_aggregation = "{agg}"\n',
            encoding="utf-8",
        )
        profile = load_profile(path)
    out = tmp_path / "pres.tsv"
    build_presentation_evidence(raw, net, mhc, profile, out)
    return read_tsv(out)[0]


def test_binding_evidence_default_stays_max(tmp_path):
    row = _binding_evidence(tmp_path)
    # BA rank 0.6 -> 1-0.6/10 = 0.94; MHCflurry percentile 3.0 -> 0.70.
    assert row["binding_evidence_score"] == "0.9400"
    assert row["binding_predictor_disagreement"] == "0.2400"
    # The shipped default profile explicitly pins the historical behaviour.
    from neoag.config import load_profile
    assert load_profile("default")["presentation_weights"]["binding_aggregation"] == "max"


def test_binding_aggregation_knob_mean_and_min(tmp_path):
    mean = _binding_evidence(tmp_path / "mean", agg="mean")
    mn = _binding_evidence(tmp_path / "min", agg="min")
    assert mean["binding_evidence_score"] == "0.8200"
    assert mn["binding_evidence_score"] == "0.7000"
    # Disagreement is reported identically for every aggregation.
    assert mean["binding_predictor_disagreement"] == mn["binding_predictor_disagreement"] == "0.2400"


def test_binding_disagreement_blank_with_one_predictor(tmp_path):
    from neoag.config import load_profile
    from neoag.presentation import build_presentation_evidence
    from neoag.utils import read_tsv

    hla = "HLA-A*02:01"
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [{"peptide_id": "P1", "peptide": "AAGIGILTV",
                     "hla_allele": hla, "mhc_class": "I"}])
    net = tmp_path / "net.tsv"
    write_tsv(net, [{"peptide": "AAGIGILTV", "hla_allele": hla,
                     "netmhcpan_ba_rank": "0.6", "netmhcpan_el_rank": "0.8"}])
    out = tmp_path / "pres.tsv"
    build_presentation_evidence(raw, net, None, load_profile("default"), out)
    row = read_tsv(out)[0]
    assert row["binding_evidence_score"] == "0.9400"
    assert row["binding_predictor_disagreement"] == ""


def test_unknown_binding_aggregation_is_rejected(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="binding_aggregation"):
        _binding_evidence(tmp_path, agg="optimistic")


def test_disagreement_column_in_presentation_schema():
    from neoag.schemas import PRESENTATION_FIELDS
    assert "binding_predictor_disagreement" in PRESENTATION_FIELDS
