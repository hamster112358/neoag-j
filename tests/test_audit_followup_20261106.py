"""Tests for the 2026-11-06 audit follow-up round (neoag-audit upload).

Each test maps to one uploaded document:
  purity fix.docx                 -> purity_confidence_route names the route,
                                     and it reaches the ccf_2 rows + reports
  CCF_multiplicity fix.docx       -> multiplicity_sensitivity is display-only
                                     pass-through into the ranked tables
  HLA fix.docx                    -> allele-level consensus as a parallel,
                                     review-only branch (pair table unchanged)
  sensitivity test on ranking.docx-> regional Dirichlet weight sensitivity on
                                     the FINAL ranking (operator correction:
                                     one mode, seeded at the default weights)
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


# --- sensitivity test on ranking: ONE regional mode (Dirichlet seeded at the
# default vector, operator correction 2026-11-06: "start from a seed of
# default, and use Dirichlet distribution to select the value for each
# weights"; wider region + more draws = the expanded range; run on the FINAL
# ranking's evidence chains; reporting window covers at least the top 50).

def _final_rows(count=20):
    """Final-ranking-shaped rows: l3_* chain columns + composite + efficacy."""
    from neoag.model_layers import default_l3_weights, l3_dimension_field

    weights = default_l3_weights({})
    total = sum(weights.values())
    names = list(weights)
    rows = []
    for i in range(count):
        values = {}
        for idx, name in enumerate(names):
            # Deterministic per-row/per-chain variation so weight draws can
            # reorder candidates, over a range that stays inside [0.05, 0.95].
            values[name] = 0.90 - 0.013 * i + 0.05 * (((i * 23 + idx * 7) % 17) / 16.0 - 0.5)
            values[name] = max(0.05, min(0.95, values[name]))
        composite = sum(weights[n] * values[n] for n in weights) / total
        row = {
            "peptide_id": f"P{i}", "peptide": "A" * 9, "hla_allele": "HLA-A*02:01",
            "mhc_class": "I", "gene": f"GENE{i}",
            "immunology_composite_score": f"{composite:.4f}",
            "efficacy_score": f"{composite * 0.8:.4f}",  # 0.8 = modifier chain M
        }
        for name in names:
            row[l3_dimension_field(name)] = f"{values[name]:.4f}"
        rows.append(row)
    return rows


def test_chains_and_draw_counts_derive_from_the_profile():
    from neoag.model_layers import default_l3_weights
    from neoag.sensitivity import (
        SENSITIVITY_DEFAULT_SAMPLES, weight_sensitivity_table,
    )

    window, chains = weight_sensitivity_table(_final_rows(5), {})
    # Chain count/names come from default_l3_weights at run time, and the
    # table SAYS them instead of the reader inferring a number.
    assert sorted(c["chain"] for c in chains) == sorted(default_l3_weights({}))
    assert len(chains) == len(default_l3_weights({}))
    assert window[0]["n_chains"] == str(len(default_l3_weights({})))
    assert window[0]["n_draws"] == str(SENSITIVITY_DEFAULT_SAMPLES + 1)

    thinned_weights = {k: (0.0 if k == "expression" else v)
                       for k, v in default_l3_weights({}).items()}
    window, chains = weight_sensitivity_table(
        _final_rows(5), {"l3_weights": thinned_weights})
    # Zeroing one weight must shrink the ACTIVE chain count by one (the
    # diagnostic keeps the zeroed row with weight 0.0000).
    assert window[0]["n_chains"] == str(len(thinned_weights) - 1)
    assert next(c for c in chains if c["chain"] == "expression")["default_weight"] == "0.0000"


def test_regional_dirichlet_is_the_only_mode_with_a_wide_default_region():
    import contextlib
    import inspect
    import io

    import neoag.presentation as presentation
    from neoag.cli import build_parser
    from neoag.sensitivity import (
        SENSITIVITY_DEFAULT_CONCENTRATION, SENSITIVITY_DEFAULT_SAMPLES,
        SENSITIVITY_DEFAULT_TOP_K, weight_sensitivity_table,
    )

    # One mode: the parameters are region width + draw count + window + seed.
    # There is no perturb parameter and no mode switch anywhere.
    sig = inspect.signature(weight_sensitivity_table)
    assert {"samples", "concentration", "top_k", "seed"} <= set(sig.parameters)
    assert "perturb" not in sig.parameters and "mode" not in sig.parameters
    assert (SENSITIVITY_DEFAULT_SAMPLES, SENSITIVITY_DEFAULT_CONCENTRATION,
            SENSITIVITY_DEFAULT_TOP_K) == (600, 3.0, 50)
    # The superseded local mode is gone from the code, not merely hidden.
    assert not hasattr(presentation, "weight_sensitivity_table")
    assert not hasattr(presentation, "SENSITIVITY_FORMULA")
    buf = io.StringIO()
    with contextlib.suppress(SystemExit), contextlib.redirect_stdout(buf):
        build_parser().parse_args(["build-presentation-evidence", "--help"])
    assert "--sensitivity-perturb" not in buf.getvalue()
    # The dedicated command owns the analysis now.
    buf2 = io.StringIO()
    with contextlib.suppress(SystemExit), contextlib.redirect_stdout(buf2):
        build_parser().parse_args(["weight-sensitivity", "--help"])
    assert "--samples" in buf2.getvalue() and "--concentration" in buf2.getvalue()


def test_dirichlet_region_is_seeded_unbiased_and_deterministic():
    from neoag.sensitivity import dirichlet_weight_vectors, sensitivity_chains

    base = {name: wt for name, _, wt in sensitivity_chains({})}
    total = sum(base.values())
    vectors = dirichlet_weight_vectors(base, 1500, 3.0, 7)
    assert vectors[0] == base  # the default vector is the anchor/seed
    for vec in vectors:
        assert abs(sum(vec.values()) - total) < 1e-9
        assert all(vec[name] >= 0.0 for name in vec)
    # Unbiased: the mean of every drawn weight sits at its default.
    for name, wt in base.items():
        mean_w = sum(vec[name] for vec in vectors) / len(vectors)
        assert abs(mean_w - wt) < 0.04, (name, mean_w, wt)
    # Determinism: the same seed reproduces every draw exactly.
    assert dirichlet_weight_vectors(base, 1500, 3.0, 7) == vectors
    assert dirichlet_weight_vectors(base, 1500, 3.0, 8) != vectors
    # A chain with default weight 0 is never switched on by a draw.
    zero = dirichlet_weight_vectors({"a": 0.5, "b": 0.0, "c": 0.5}, 10, 3.0, 1)
    assert all(vec["b"] == 0.0 for vec in zero)


def test_anchor_vector_reproduces_the_stored_efficacy_exactly():
    from neoag.sensitivity import weight_sensitivity_table

    profile = {"l3_weights": {"alpha": 0.5, "beta": 0.5}}
    row = {
        "peptide_id": "ANCH", "peptide": "A" * 9, "hla_allele": "HLA-A*02:01",
        "l3_alpha_score": "0.6000", "l3_beta_score": "0.8000",
        "immunology_composite_score": "0.7000",
        "efficacy_score": "0.5600",  # = 0.7000 * 0.8 multiplier chain
    }
    window, _chains = weight_sensitivity_table([row], profile, samples=40)
    # score(default vector) == stored efficacy_score; a lone candidate can
    # never fall out of any window, so it is trivially consistent.
    assert window[0]["base_score"] == "0.5600"
    assert window[0]["base_rank"] == "1"
    assert window[0]["verdict"] == "rank_stable"
    assert window[0]["consistency"] == "1.0000"
    assert window[0]["top_selection"] == "yes"


def test_single_chain_driver_is_penalized_by_the_dominance_factor():
    from neoag.sensitivity import weight_sensitivity_table

    profile = {"l3_weights": {"alpha": 0.45, "beta": 0.45, "gamma": 0.10}}
    common = {"peptide": "A" * 9, "hla_allele": "HLA-A*02:01", "mhc_class": "I",
              "l3_gamma_score": "0.0000"}
    bal = {"peptide_id": "BAL", **common,
           "l3_alpha_score": "1.0000", "l3_beta_score": "1.0000",
           "immunology_composite_score": "0.9000", "efficacy_score": "0.9000"}
    one = {"peptide_id": "ONE", **common,
           "l3_alpha_score": "1.0000", "l3_beta_score": "0.0000",
           "immunology_composite_score": "0.4500", "efficacy_score": "0.4500"}
    rows, _chains = weight_sensitivity_table([bal, one], profile, samples=40)
    out = {r["peptide_id"]: r for r in rows}
    assert out["ONE"]["dominance"] == "1.0000"
    assert out["ONE"]["effective_chains"] == "1.00"
    assert out["ONE"]["dominance_factor"] == "0.5000"
    assert float(out["BAL"]["dominance_factor"]) == 1.0
    # expected_score is the mean over region draws (unbiased, small sampling
    # error), so compare with a tolerance: BAL ~ 0.90, ONE ~ 0.45 * 0.5.
    assert abs(float(out["BAL"]["adjusted_score"]) - 0.90) < 0.10
    assert abs(float(out["ONE"]["adjusted_score"]) - 0.225) < 0.05
    assert int(out["BAL"]["adjusted_rank"]) < int(out["ONE"]["adjusted_rank"])


def test_window_reports_at_least_the_top_50_but_selection_stays_downstream():
    from neoag.sensitivity import SENSITIVITY_DEFAULT_TOP_K, weight_sensitivity_table

    window, _chains = weight_sensitivity_table(_final_rows(60), {})
    # The whole pool (60) was re-ranked; the REPORT window keeps every
    # candidate inside the default or median top-50 - the 13th and 21st are
    # covered, and so is the 50th; nothing about the window selects a final
    # candidate list - that stays a review decision.
    assert SENSITIVITY_DEFAULT_TOP_K == 50
    assert window[0]["n_scorable"] == "60"
    assert len(window) >= 50
    by_rank = {int(r["base_rank"]): r for r in window}
    assert 13 in by_rank and 21 in by_rank
    assert any(int(r["base_rank"]) >= 50 for r in window) or len(window) >= 50


def test_borderline_candidates_have_their_flips_marked():
    from neoag.sensitivity import weight_sensitivity_table

    profile = {"l3_weights": {"alpha": 0.5, "beta": 0.5}}
    common = {"peptide": "A" * 9, "hla_allele": "HLA-A*02:01", "mhc_class": "I"}
    r1 = {"peptide_id": "R1", **common, "l3_alpha_score": "1.0000", "l3_beta_score": "0.0000",
          "immunology_composite_score": "1.0000", "efficacy_score": "1.0000"}
    r2 = {"peptide_id": "R2", **common, "l3_alpha_score": "0.0000", "l3_beta_score": "1.0000",
          "immunology_composite_score": "1.0000", "efficacy_score": "1.0000"}
    # Perfect tie under the default; a Dirichlet draw favours one side or the
    # other, so neither candidate may sit in the top-1 under every draw.
    rows, _chains = weight_sensitivity_table([r1, r2], profile, samples=400, top_k=1)
    out = {r["peptide_id"]: r for r in rows}
    for rid in ("R1", "R2"):
        assert 0.15 < float(out[rid]["consistency"]) < 0.85
        assert out[rid]["top_selection"] == "no"
        assert int(out[rid]["rank_spread"]) >= 1
    solo = weight_sensitivity_table([r1], profile, samples=40, top_k=1)[0][0]
    assert solo["consistency"] == "1.0000" and solo["top_selection"] == "yes"


def test_technical_report_shows_both_tables_the_chains_and_the_formula(tmp_path):
    from neoag.reports_dual import ReportBundle, make_technical_report
    from neoag.sensitivity import weight_sensitivity_table

    sens, chains = weight_sensitivity_table(_final_rows(8), {})
    bundle = ReportBundle(profile={}, events=[], peptides=[],
                          weight_sensitivity=sens, weight_sensitivity_chains=chains)
    make_technical_report(tmp_path / "tech.html", bundle)
    html = (tmp_path / "tech.html").read_text(encoding="utf-8")
    assert "Weight Sensitivity (regional Dirichlet)" in html
    assert "Table 1: default-weight ranking" in html
    assert "Table 2: stability-adjusted ranking" in html
    assert "adjusted_score = expected_score" in html
    assert "dominance_factor = min(1, effective_chains/2)" in html
    assert "robustness to the weight choice, not accuracy" in html
    # The per-chain diagnostic lists the real chain names.
    assert "hla_presentation" in html


def test_invalid_sensitivity_arguments_are_rejected():
    import pytest

    from neoag.sensitivity import weight_sensitivity_table

    with pytest.raises(ValueError, match="samples"):
        weight_sensitivity_table(_final_rows(3), {}, samples=0)
    with pytest.raises(ValueError, match="concentration"):
        weight_sensitivity_table(_final_rows(3), {}, concentration=0.0)
    with pytest.raises(ValueError, match="top_k"):
        weight_sensitivity_table(_final_rows(3), {}, top_k=0)


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
