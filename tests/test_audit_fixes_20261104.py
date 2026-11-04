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


def test_consensus_vocabulary_and_blank_status_not_lost(tmp_path):
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
    assert load_hla_loh(tsv) == {"HLA-A*02:01", "HLA-A*11:01"}


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
