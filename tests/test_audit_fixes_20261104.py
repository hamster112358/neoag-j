"""Regression tests for the 2026-11-04 audit findings (neoag-audit repo).

APPM finding 2: appm_v2.load_hla_loh() ignored consensus_status and counted every
allele of a status-bearing table as lost (missing status == lost), which capped the
MHC-I integrity score at 0.45 for a report whose consensus had one lost allele.

Mutant-RNA finding 1: evidence_states.derive_rna_support() returned RNA_CONFIRMED
grade 3 from the RNA_ALT_SUPPORTED status string alone, while evidence_layer marks
any alt>0 as RNA_ALT_SUPPORTED - one read reached the highest level, bypassing the
snv_min_alt_reads / snv_min_vaf / junction thresholds.

CCF C2/C4: _probability_from_interval fabricated 0/1 certainty without a usable
interval, and the probability_* column names overstated an interval fraction.

APPM F6: immune_escape returned a hardcoded "LOW" for cytotoxic-killing resistance.

HLA A2/A3: netmhcstabpan was read and written through the presentation table but
excluded from the weighted score; the weight is now profile-tunable so saved score
tables can be reranked without rerunning predictors.
"""
import csv

from neoag.utils import write_tsv

from neoag.appm_v2 import load_hla_loh
from neoag.evidence_consensus import load_consensus_rules
from neoag.evidence_states import derive_rna_support


def _write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_consensus_status_table_marks_only_lost_alleles(tmp_path):
    tsv = _write_tsv(
        tmp_path / "hla_loh_consensus.tsv",
        ["sample_id", "hla_allele", "locus", "consensus_status"],
        [
            ["S1", "HLA-A*02:06", "A", "Retained"],
            ["S1", "HLA-A*30:01", "A", "Discordant"],
            ["S1", "HLA-B*13:02", "B", "Lost"],
            ["S1", "HLA-B*48:01", "B", "Retained"],
            ["S1", "HLA-C*06:02", "C", "Discordant"],
            ["S1", "HLA-C*08:01", "C", "Retained"],
        ],
    )
    assert load_hla_loh(tsv) == {"HLA-B*13:02"}


def test_consensus_vocabulary_blank_and_single_tool_not_lost(tmp_path):
    # F3 policy (audit): a single-tool loss claim is not a consensus loss.
    tsv = _write_tsv(
        tmp_path / "consensus2.tsv",
        ["hla_allele", "consensus_status"],
        [
            ["HLA-A*02:01", "CONSENSUS_LOST"],
            ["HLA-A*11:01", "SINGLE_TOOL_LOST"],
            ["HLA-B*07:02", "CONSENSUS_RETAINED"],
            ["HLA-C*07:01", ""],
        ],
    )
    assert load_hla_loh(tsv) == {"HLA-A*02:01"}


def test_legacy_allele_only_table_still_lists_lost(tmp_path):
    tsv = _write_tsv(
        tmp_path / "lost_only.tsv",
        ["sample_id", "hla_allele"],
        [["S1", "HLA-A*24:02"], ["S1", "HLA-B*08:01"]],
    )
    assert load_hla_loh(tsv) == {"HLA-A*24:02", "HLA-B*08:01"}


def test_single_alt_read_is_not_rna_confirmed():
    rules = load_consensus_rules()
    row = {
        "event_type": "SNV",
        "rna_support_status": "RNA_ALT_SUPPORTED",
        "rna_alt_reads": "1",
        "rna_depth": "14",
        "rna_vaf": "0.07",
    }
    state = derive_rna_support(row, rules)
    assert state["state"] == "RNA_LOW_SUPPORT"
    assert state["grade"] == 1


def test_status_string_cannot_beat_thresholds_but_numbers_confirm():
    rules = load_consensus_rules()
    strong = derive_rna_support(
        {"event_type": "SNV", "rna_support_status": "RNA_ALT_SUPPORTED",
         "rna_alt_reads": "8", "rna_depth": "20", "rna_vaf": "0.40"},
        rules,
    )
    assert strong["state"] == "RNA_CONFIRMED" and strong["grade"] == 3
    medium = derive_rna_support(
        {"event_type": "SNV", "rna_support_status": "RNA_ALT_SUPPORTED",
         "rna_alt_reads": "4", "rna_depth": "20", "rna_vaf": "0.20"},
        rules,
    )
    assert medium["state"] == "RNA_CONFIRMED" and medium["grade"] == 2


def test_junction_status_with_one_read_is_low_support():
    rules = load_consensus_rules()
    weak = derive_rna_support(
        {"event_type": "FUSION", "rna_support_status": "RNA_JUNCTION_SUPPORTED",
         "rna_junction_reads": "1"},
        rules,
    )
    assert weak["state"] == "RNA_LOW_SUPPORT"
    ok = derive_rna_support(
        {"event_type": "FUSION", "rna_support_status": "RNA_JUNCTION_SUPPORTED",
         "rna_junction_reads": "5"},
        rules,
    )
    assert ok["state"] == "RNA_CONFIRMED"


def test_curated_status_without_numbers_still_confirms():
    rules = load_consensus_rules()
    state = derive_rna_support(
        {"event_type": "SNV", "rna_support_status": "RNA_ALT_SUPPORTED"}, rules
    )
    assert state["state"] == "RNA_CONFIRMED" and state["grade"] == 3


def test_interval_without_usable_ci_is_unassessed_not_certain():
    from neoag.ccf_v2 import _interval_fraction_over_threshold
    assert _interval_fraction_over_threshold(0.9, None, None, 0.8, "gt") == ""
    assert _interval_fraction_over_threshold(0.1, None, None, 0.3, "lt") == ""
    assert _interval_fraction_over_threshold(0.9, 0.95, 0.8, 0.8, "gt") == ""
    assert _interval_fraction_over_threshold(0.8, 0.7, 0.9, 0.8, "gt") == "0.5000"
    assert _interval_fraction_over_threshold(0.9, 0.85, 0.95, 0.8, "gt") == "1.0000"


def test_output_columns_use_honest_interval_fraction_names():
    from neoag.ccf_v2 import CCF_FIELDS
    assert "interval_fraction_ccf_above_0_8" in CCF_FIELDS
    assert "probability_ccf_gt_0_8" not in CCF_FIELDS


def test_cytotoxic_resistance_is_not_assessed_until_wired(tmp_path):
    from neoag.immune_escape import build_immune_escape_evidence
    peps = _write_tsv(
        tmp_path / "peptides.tsv",
        ["peptide_id", "event_id", "peptide", "hla_allele", "mhc_class"],
        [["P1", "E1", "AAAAAAAAA", "HLA-A*02:01", "I"]],
    )
    out = tmp_path / "escape"
    build_immune_escape_evidence(sample_id="S1", raw_peptides=peps, outdir=out)
    rows = []
    for f in sorted(out.glob("*.tsv")):
        with f.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if "cytotoxic_killing_resistance_status" in row:
                    rows.append(row)
    assert rows and all(
        row["cytotoxic_killing_resistance_status"] == "NOT_ASSESSED" for row in rows
    )


def _presentation_scores(tmp_path, stability_weight=None):
    """Saved score tables only -> presentation_evidence_score per peptide id.

    Two peptides are identical in NetMHCpan / MHCflurry evidence and differ only
    in NetMHCstabpan rank (P1 = 0.5 stable, P2 = 8 unstable). Without a weight
    override the default netmhcstabpan_stability weight applies; passing one
    writes a throwaway profile (the audit's sensitivity rerun, HLA A3).
    """
    from neoag.presentation import build_presentation_evidence
    from neoag.config import load_profile
    from neoag.utils import read_tsv

    hla = "HLA-A*02:01"
    pep_rows = [["P1", "E1", "AAGIGILTV", hla, "I"], ["P2", "E2", "SLYYRRSLNV", hla, "I"]]
    raw = _write_tsv(tmp_path / "raw.tsv",
                     ["peptide_id", "event_id", "peptide", "hla_allele", "mhc_class"], pep_rows)
    net = _write_tsv(tmp_path / "net.tsv",
                     ["peptide", "hla_allele", "netmhcpan_ba_rank", "netmhcpan_el_rank"],
                     [[r[2], hla, "0.6", "0.8"] for r in pep_rows])
    mhc = _write_tsv(
        tmp_path / "mhc.tsv",
        ["peptide", "hla_allele", "mhcflurry_affinity_percentile",
         "mhcflurry_processing_score", "mhcflurry_presentation_score"],
        [[r[2], hla, "0.5", "0.7", "0.8"] for r in pep_rows])
    stab = _write_tsv(tmp_path / "stab.tsv", ["peptide", "hla_allele", "netmhcstabpan_rank"],
                      [["AAGIGILTV", hla, "0.5"], ["SLYYRRSLNV", hla, "8"]])
    if stability_weight is None:
        profile = load_profile("default")
    else:
        profile_path = tmp_path / "profile.toml"
        profile_path.write_text(
            "[presentation_weights]\n"
            "netmhcpan_ba = 0.25\nnetmhcpan_el = 0.35\n"
            "mhcflurry_presentation = 0.30\nmhcflurry_processing = 0.10\n"
            "netchop_processing = 0.10\n"
            "netmhcstabpan_stability = " + str(stability_weight) + "\n",
            encoding="utf-8",
        )
        profile = load_profile(profile_path)
    out = tmp_path / "pres.tsv"
    build_presentation_evidence(raw, net, mhc, profile, out, stab)
    return {r["peptide_id"]: float(r["presentation_evidence_score"]) for r in read_tsv(out)}


def test_stabpan_rank_now_moves_presentation_score(tmp_path):
    scores = _presentation_scores(tmp_path)
    # Before the fix the stability table was read and passed through but never
    # scored: identical evidence gave identical presentation scores.
    assert scores["P1"] > scores["P2"]


def test_stability_weight_is_profile_tunable(tmp_path):
    base = _presentation_scores(tmp_path / "zero", stability_weight=0.0)
    assert abs(base["P1"] - base["P2"]) < 1e-9
    default = _presentation_scores(tmp_path / "low", stability_weight=0.05)
    heavy = _presentation_scores(tmp_path / "high", stability_weight=0.5)
    assert 0 < default["P1"] - default["P2"] < heavy["P1"] - heavy["P2"]


def test_build_presentation_cli_exposes_netchop():
    from neoag.cli import build_parser
    parser = build_parser()
    args = parser.parse_args([
        "build-presentation-evidence", "--raw-peptides", "r.tsv",
        "--out", "o.tsv", "--netchop", "chop.tsv",
    ])
    assert args.netchop == "chop.tsv"


# ---------------------------------------------------------------------------
# Second audit batch (2026-11-04): HLA A1, QC A3, B3, C1/C3/C5, F5
# ---------------------------------------------------------------------------

def test_a1_consensus_preserves_ranked_votes_and_tie():
    from neoag.agent_skills.hla_typing_compare import consensus

    def row(tool, a1, a2):
        return {"tool": tool, "locus": "A", "lowres1": a1, "lowres2": a2,
                "highres1": "", "highres2": ""}

    discordant = consensus([
        row("OptiType", "HLA-A*02:01", "HLA-A*30:01"),
        row("SpecHLA", "HLA-A*02:06", "HLA-A*30:01"),
        row("HLA-LA", "HLA-A*02:01", "HLA-A*68:01"),
    ])
    row_a = next(r for r in discordant if r["status"] == "DISCORDANT")
    assert row_a["tie"] == "yes"
    votes = row_a["ranked_votes"].split(";")
    assert len(votes) == 3 and all(v.endswith("=1") for v in votes)
    # Deterministic: alphabetical among ties, so the "best" pick is reproducible.
    assert row_a["consensus_2field"] == min(v.split("=")[0] for v in votes)

    agree = consensus([
        row("OptiType", "HLA-A*02:01", "HLA-A*30:01"),
        row("SpecHLA", "HLA-A*02:01", "HLA-A*30:01"),
        row("HLA-LA", "HLA-A*02:01", "HLA-A*68:01"),
    ])
    row_b = next(r for r in agree if r["status"] == "CONSENSUS")
    assert row_b["tie"] == "no"
    assert row_b["ranked_votes"] == "HLA-A*02:01 / HLA-A*30:01=2;HLA-A*02:01 / HLA-A*68:01=1"


def test_a3_weight_sensitivity_exposes_weight_driven_ranks():
    from neoag.presentation import weight_sensitivity_table

    rows = [
        {"peptide_id": "P1", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
         "netmhcpan_ba_rank": "1", "netmhcpan_el_rank": "99", "mhcflurry_presentation_score": "",
         "mhcflurry_processing_score": "", "netchop_31d_cterm_score": "", "netmhcstabpan_rank": ""},
        {"peptide_id": "P2", "peptide": "BBBBBBBBB", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
         "netmhcpan_ba_rank": "99", "netmhcpan_el_rank": "1", "mhcflurry_presentation_score": "",
         "mhcflurry_processing_score": "", "netchop_31d_cterm_score": "", "netmhcstabpan_rank": ""},
    ]
    # With the default weights the el weight (0.35) beats ba (0.25); at a +-50%
    # sweep the ba-up vector flips them - the table must reveal that.
    out = weight_sensitivity_table(rows, {}, perturb=0.5)
    assert out and all(r["n_vectors"] == "13" for r in out)
    assert {r["peptide_id"] for r in out} == {"P1", "P2"}
    p1 = next(r for r in out if r["peptide_id"] == "P1")
    assert int(p1["rank_spread"]) >= 1
    assert p1["verdict"] in {"rank_minor_shift", "rank_unstable"}
    assert float(p1["score_range"]) > 0


def test_b3_unresolved_immunogenicity_downgrades_instead_of_silent_redistribution():
    from neoag.scoring import compute_peptide_efficacy

    weights = {"immunogenicity": 0.15, "presentation_evidence": 0.25, "binding_evidence": 0.20, "event_score": 0.35}
    prof = {"score_weights": weights, "scoring": {"immunogenicity_unresolved_multiplier": 0.90},
            "gates": {"enabled": False}, "ccf_lite": {}}
    prof100 = {"score_weights": weights, "scoring": {"immunogenicity_unresolved_multiplier": 1.0},
               "gates": {"enabled": False}, "ccf_lite": {}}
    pep = {"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
           "binding_rank": "1", "el_rank": "1", "wildtype_peptide": "TTTTTTTTT"}
    ev = {"event_id": "E1", "event_type": "SNV", "event_confidence": "0.9"}
    pres = {"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01",
            "presentation_evidence_score": "0.9", "binding_evidence_score": "0.9",
            "evidence_completeness": "1.0"}
    no_immuno = compute_peptide_efficacy(pep, ev, pres, prof)
    no_immuno_full = compute_peptide_efficacy(pep, ev, pres, prof100)
    assert no_immuno["immunogenicity_resolved"] == "no"
    assert float(no_immuno["efficacy_score"]) == round(0.90 * float(no_immuno_full["efficacy_score"]), 4)
    resolved = compute_peptide_efficacy({**pep, "prime_score": "0.8000"}, ev, pres, prof)
    resolved_full = compute_peptide_efficacy({**pep, "prime_score": "0.8000"}, ev, pres, prof100)
    assert resolved["immunogenicity_resolved"] == "yes"
    assert resolved["efficacy_score"] == resolved_full["efficacy_score"]


def _ccf_rows(tmp_path, cnv_rows, events_rows, purity="0.70"):
    from neoag.ccf_v2 import build_ccf_2
    from neoag.config import load_profile

    events = tmp_path / "events.tsv"
    purity_f = tmp_path / "purity.tsv"
    cnv = tmp_path / "cnv.tsv"
    write_tsv(events, events_rows)
    write_tsv(purity_f, [{"purity": purity, "ploidy": "2", "source": "FACETS", "confidence": "high"}])
    write_tsv(cnv, cnv_rows)
    return build_ccf_2(events, purity_f, cnv, load_profile("default"), tmp_path / "ccf.tsv")


def _event(eid, vaf, depth=100):
    return {"event_id": eid, "sample_id": "S1", "mutation_source": "SNV",
            "chrom": "chr9", "pos": "150", "tumor_vaf": vaf,
            "tumor_depth": str(depth), "tumor_alt_count": str(int(float(vaf) * depth))}


def test_c1_diploid_fallback_never_pretends_to_be_well_supported(tmp_path):
    # No CNV segment overlaps chr9 -> diploid fallback even at clonal-looking VAF.
    rows = _ccf_rows(tmp_path,
                     [{"chrom": "chr1", "start": "1", "end": "1000", "total_cn": "3", "major_cn": "2", "minor_cn": "1"}],
                     [_event("E1", "0.35")])
    row = rows[0]
    assert "DIPLOID_FALLBACK" in row["local_cnv_status"]
    assert row["ccf_confidence"] == "low"
    assert row["clonality_status"] != "clonal_compatible"
    assert float(row["clonality_multiplier"]) < 1.0


def test_c3_multiplicity_sensitivity_flags_material_candidates(tmp_path):
    rows = _ccf_rows(
        tmp_path,
        [{"chrom": "chr9", "start": "1", "end": "1000", "total_cn": "2", "major_cn": "2", "minor_cn": "0"}],
        [_event("MAT", "0.25", 200), _event("STB", "0.05", 200)],
        purity="0.50",
    )
    by_id = {r["event_id"]: r for r in rows}
    assert by_id["MAT"]["multiplicity_sensitivity"].startswith("material:")
    assert by_id["STB"]["multiplicity_sensitivity"] == "stable"


def test_c5_purity_confidence_criteria_are_explicit(tmp_path):
    from neoag.ccf_v2 import PURITY_CONFIDENCE_CRITERIA, load_purity_info

    purity = tmp_path / "purity.tsv"
    write_tsv(purity, [{"purity": "0.80", "ploidy": "2"}])
    info = load_purity_info(purity)
    assert info["purity_confidence"] == "high"
    assert info["purity_confidence_criteria"] == PURITY_CONFIDENCE_CRITERIA
    assert "0.20<=purity<0.35=medium" in PURITY_CONFIDENCE_CRITERIA


def test_f5_low_expression_suggests_but_does_not_prove_biallelic_loss():
    from neoag.appm_v2 import _biallelic_status, _functional_status, _gene_integrity_score, _gene_integrity_status

    status, reason = _biallelic_status("TP53", {"TP53": ["stop_gained"]}, {}, {"TP53": 0.01}, True, 1.0)
    assert status == "SUSPECTED_BIALLELIC_LOSS"
    assert "no_second_hit_evidence" in reason
    assert _functional_status(status) == "caution"
    integrity = _gene_integrity_status(biallelic_status=status, damaging=True, copy_status="",
                                       loh_status="", expression_status="low")
    assert integrity == "suspected_biallelic_loss"
    assert _gene_integrity_score(integrity) == 0.25
    # Hard copy-level evidence still gives the hard call.
    hard, _ = _biallelic_status("TP53", {"TP53": ["stop_gained"]},
                                {"TP53": {"copy_number_status": "homozygous deletion", "loh_status": ""}},
                                {"TP53": 8.0}, True, 1.0)
    assert hard == "BIALLELIC_LOSS"
