"""Regression tests for the G-b package (B-E review C1-02/03/04/05/06/08/09/11).

Each test names the audit item it locks:

* C1-02 / C1-09  - absence is persisted as "" everywhere, never as the
  old 99.0 / 0.0 / "99" / "0" sentinels; measured zeros survive.
* C1-03          - the local IEDB-style heuristic never enters the
  tool-named evidence column; it is off by default and, when explicitly
  enabled, lands in separate *_local_estimate columns.
* C1-04          - annotation missingness is "" (no "NA" placeholders).
* C1-05          - one non-standard residue blanks its own IEDB row with a
  reason instead of crashing the batch.
* C1-06          - a stored tool output is only reused when a sidecar
  fingerprint of peptides + version + mode still matches.
* C1-08          - no predictor result means the restricting HLA stays
  unknown; the first requested allele is never fabricated as the binder.
* C1-11          - trained vs pseudoallele vs unassessed allele support is
  recorded at predict time and carried to the presentation columns.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from neoag.adapters.iedb_immunogenicity import (
    nonstandard_residues,
    predict_immunogenicity,
    predict_immunogenicity_checked,
    predict_pairs,
)
from neoag.adapters.mhcflurry import parse_mhcflurry
from neoag.adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_standard_xls
from neoag.adapters.peptide_netmhcpan import annotate_variant_peptide_row
from neoag.config import load_profile
from neoag.presentation import build_presentation_evidence
from neoag.utils import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# C1-05: peptide alphabet validation in the IEDB-style heuristic
# ---------------------------------------------------------------------------

def test_nonstandard_residues_helper():
    assert nonstandard_residues("AAGIGILTV") == ""
    # Sorted, unique, uppercase: X, U, B, Z and the stop-marker "*" all listed.
    assert nonstandard_residues("AXU*BZ") == "*BUXZ"


def test_predict_immunogenicity_raises_clean_error_for_unreadable_residue():
    # Previously this was a bare KeyError("X") from IMMUNOSCALE[aa].
    with pytest.raises(ValueError, match="residue"):
        predict_immunogenicity("AXGIGILTV", "HLA-A*02:01")


def test_predict_checked_reports_reason_without_assessment():
    score, status = predict_immunogenicity_checked("AXGIGILTV", "HLA-A*02:01")
    assert score is None
    assert status == "UNASSESSED_NONSTANDARD_RESIDUE:X"
    score, status = predict_immunogenicity_checked("", "HLA-A*02:01")
    assert score is None
    assert status == "UNASSESSED_EMPTY_PEPTIDE"
    score, status = predict_immunogenicity_checked("AAGIGILTV", "HLA-A*02:01")
    assert score is not None
    assert status == "ASSESSED"


def test_predict_pairs_blank_one_bad_row_and_continue_others():
    # C1-05: "one X in an SV/fusion ORF peptide kills the whole annotation
    # run" - after the fix it blanks exactly its own row and carries why.
    rows = predict_pairs([
        ("AAGIGILTV", "HLA-A*02:01"),
        ("AXGIGILTV", "HLA-A*02:01"),
        ("SIINFEKL", "H-2-Kb"),
    ])
    assert len(rows) == 3
    assert rows[0]["iedb_assessment_status"] == "ASSESSED"
    assert rows[0]["iedb_immunogenicity_score"] != ""
    assert rows[1]["iedb_immunogenicity_score"] == ""
    assert rows[1]["iedb_assessment_status"] == "UNASSESSED_NONSTANDARD_RESIDUE:X"
    assert rows[2]["iedb_assessment_status"] == "ASSESSED"
    assert rows[2]["iedb_immunogenicity_score"] != ""


# ---------------------------------------------------------------------------
# C1-09: MHCflurry absence is typed missingness, not "99"/"0" predictions
# ---------------------------------------------------------------------------

def test_parse_mhcflurry_missing_metrics_blank_measured_zero_preserved(tmp_path):
    csv_path = tmp_path / "mhcflurry.csv"
    csv_path.write_text(
        "peptide,allele,mhcflurry_affinity_percentile,mhcflurry_processing_score\n"
        "MEASUREDZERO,HLA-A*02:01,0.0,0\n"
        "GAPS,HLA-A*02:01,,\n"
        "GARBAGE,HLA-A*02:01,abc,xyz\n",
        encoding="utf-8",
    )
    rows = {r["peptide"]: r for r in parse_mhcflurry(csv_path, "S1")}
    # C1-02 disposition: "protecting legitimate zero measurements".
    assert rows["MEASUREDZERO"]["mhcflurry_affinity_percentile"] == "0.0"
    assert rows["MEASUREDZERO"]["mhcflurry_processing_score"] == "0.0"
    # Absence (and unparseable cells) are blank, not 99.0 / 0.0.
    assert rows["GAPS"]["mhcflurry_affinity_percentile"] == ""
    assert rows["GAPS"]["mhcflurry_processing_score"] == ""
    assert rows["GAPS"]["mhcflurry_affinity"] == ""
    assert rows["GARBAGE"]["mhcflurry_affinity_percentile"] == ""
    # C1-11: MHCflurry reports no per-allele training coverage, so every row
    # records the coverage explicitly as unassessed.
    assert rows["MEASUREDZERO"]["mhcflurry_allele_support_status"] == "MODEL_COVERAGE_UNASSESSED"


# ---------------------------------------------------------------------------
# C1-02 + C1-11: NetMHCpan parse/xls round-trip
# ---------------------------------------------------------------------------

def test_parse_netmhcpan_blank_cells_stay_blank_not_sentinels(tmp_path):
    xls = tmp_path / "netmhcpan.xls"
    xls.write_text(
        "Pos\tPeptide\tHLA\tScore_EL\t%Rank_EL\tScore_BA\t%Rank_BA\tBindLevel\n"
        "0\tAAGIGILTV\tHLA-A*02:01\t0.1\t\t50.0\t\t\n",
        encoding="utf-8",
    )
    row = parse_netmhcpan(xls, "S1")[0]
    # Empty rank cells used to be written as 99.0/0.0; now they stay "".
    assert row["netmhcpan_el_rank"] == ""
    assert row["netmhcpan_ba_rank"] == ""
    assert row["netmhcpan_ba_score"] == "50.0"
    assert row["netmhcpan_el_score"] == "0.1"
    # A third-party table cannot say anything about allele support.
    assert row["netmhcpan_allele_support_status"] == ""


def test_standard_xls_allele_support_roundtrip(tmp_path):
    # C1-11: the trained/pseudoallele statement recorded at run time must
    # survive the xls -> evidence-TSV hop.
    out = tmp_path / "netmhcpan.xls"
    write_netmhcpan_standard_xls(out, [{
        "Peptide": "AAGIGILTV", "HLA": "HLA-A*02:01", "Score_EL": "0.1",
        "%Rank_EL": "0.5", "Score_BA": "50", "%Rank_BA": "1.0",
        "AlleleSupport": "PSEUDOALLELE_EXTRAPOLATED",
    }])
    row = parse_netmhcpan(out, "S1")[0]
    assert row["netmhcpan_allele_support_status"] == "PSEUDOALLELE_EXTRAPOLATED"


# ---------------------------------------------------------------------------
# C1-11: runner mode detection
# ---------------------------------------------------------------------------

def _peplist_line(allele: str, peptide: str) -> str:
    return f"1 {allele} {peptide} PEPLIST 0.1 0.5 0.2 1.0 50"


def test_local_by_allele_records_pseudoallele_notice(tmp_path, monkeypatch):
    from neoag.tools import runner as runner_mod

    def fake_run_cmd(cmd, cwd, env=None):
        # The -a argument carries netmhcpan_allele_string()'s compact form:
        # HLA-A*02:01 -> HLA-A0201, HLA-B*07:02 -> HLA-B0702.
        allele = cmd[cmd.index("-a") + 1]
        peptide = "AAGIGILTV" if allele == "HLA-A0201" else "SIINFEKL"
        stdout = _peplist_line(allele, peptide) + "\n"
        if allele == "HLA-B0702":
            stdout = (
                "WARNING: Allele HLA-B*07:02 is not in the set of trained "
                "alleles; results are computed with pseudoallele HLA-B*07:03P\n"
                + stdout
            )
        return SimpleNamespace(returncode=0, stdout=stdout, stderr="")

    monkeypatch.setattr(runner_mod, "_run_cmd", fake_run_cmd)
    out_xls = tmp_path / "tools" / "netmhcpan.xls"
    runner_mod._run_netmhcpan_local_by_allele(
        [("AAGIGILTV", "HLA-A*02:01"), ("SIINFEKL", "HLA-B*07:02")],
        out_xls,
        netmhcpan_exe="netmhcpan_real",
    )
    rows = {r["peptide"]: r for r in parse_netmhcpan(out_xls, "S1")}
    assert rows["AAGIGILTV"]["netmhcpan_allele_support_status"] == "TRAINED_ALLELE"
    assert rows["SIINFEKL"]["netmhcpan_allele_support_status"] == "PSEUDOALLELE_EXTRAPOLATED"


def test_local_pmhc_mode_records_unassessed_allele_support(tmp_path, monkeypatch):
    from neoag.tools import runner as runner_mod

    def fake_run_cmd(cmd, cwd, env=None):
        return SimpleNamespace(
            returncode=0,
            stdout=_peplist_line("HLA-A*02:01", "AAGIGILTV") + "\n",
            stderr="",
        )

    monkeypatch.setattr(runner_mod, "_run_cmd", fake_run_cmd)
    out_xls = tmp_path / "tools" / "netmhcpan.xls"
    runner_mod._run_netmhcpan_local(
        [("AAGIGILTV", "HLA-A*02:01")], out_xls, netmhcpan_exe="netmhcpan_real",
    )
    row = parse_netmhcpan(out_xls, "S1")[0]
    # The pan-mode output gives no per-allele trained/pseudoallele statement;
    # the runner must say UNASSESSED rather than claim support it cannot see.
    assert row["netmhcpan_allele_support_status"] == "ALLELE_SUPPORT_UNASSESSED"


# ---------------------------------------------------------------------------
# C1-11 + C1-02: presentation copies the statuses and persists blanks
# ---------------------------------------------------------------------------

def _presentation_fixture(tmp_path, *, net_rows, mhc_rows):
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "peptide": "AAGIGILTV", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    net = tmp_path / "net.tsv"
    net.write_text(net_rows, encoding="utf-8")
    mhc = tmp_path / "mhc.tsv"
    mhc.write_text(mhc_rows, encoding="utf-8")
    out = tmp_path / "pres.tsv"
    build_presentation_evidence(raw, net, mhc, load_profile("default"), out)
    return read_tsv(out)[0]


def test_presentation_records_status_from_new_style_evidence(tmp_path):
    row = _presentation_fixture(
        tmp_path,
        net_rows=(
            "peptide\thla_allele\tnetmhcpan_ba_rank\tnetmhcpan_el_rank\t"
            "netmhcpan_allele_support_status\n"
            "AAGIGILTV\tHLA-A*02:01\t0.6\t0.8\tPSEUDOALLELE_EXTRAPOLATED\n"
        ),
        mhc_rows=(
            "peptide\thla_allele\tmhcflurry_affinity_percentile\t"
            "mhcflurry_allele_support_status\n"
            "AAGIGILTV\tHLA-A*02:01\t3.0\tMODEL_COVERAGE_UNASSESSED\n"
        ),
    )
    assert row["netmhcpan_allele_support_status"] == "PSEUDOALLELE_EXTRAPOLATED"
    assert row["mhcflurry_allele_support_status"] == "MODEL_COVERAGE_UNASSESSED"
    # The rollup takes the most pessimistic statement, and extrapolation is
    # judged on the NetMHCpan chain alone.
    assert row["predictor_allele_support_status"] == "PSEUDOALLELE_EXTRAPOLATED"
    assert row["allele_extrapolation_status"] == "EXTRAPOLATED"


def test_presentation_records_trained_allele_status(tmp_path):
    row = _presentation_fixture(
        tmp_path,
        net_rows=(
            "peptide\thla_allele\tnetmhcpan_ba_rank\tnetmhcpan_el_rank\t"
            "netmhcpan_allele_support_status\n"
            "AAGIGILTV\tHLA-A*02:01\t0.6\t0.8\tTRAINED_ALLELE\n"
        ),
        mhc_rows="peptide\thla_allele\n",
    )
    assert row["predictor_allele_support_status"] == "TRAINED_ALLELE"
    assert row["allele_extrapolation_status"] == "NOT_EXTRAPOLATED"


def test_presentation_legacy_evidence_says_unassessed_not_trained(tmp_path):
    # Evidence tables written before the C1-11 producers have no status
    # column: "show unknown when unsupported", never a guessed TRAINED.
    row = _presentation_fixture(
        tmp_path,
        net_rows="peptide\thla_allele\tnetmhcpan_ba_rank\tnetmhcpan_el_rank\n"
                 "AAGIGILTV\tHLA-A*02:01\t0.6\t0.8\n",
        mhc_rows="peptide\thla_allele\n",
    )
    assert row["netmhcpan_allele_support_status"] == "ALLELE_SUPPORT_UNASSESSED"
    assert row["mhcflurry_allele_support_status"] == ""
    assert row["predictor_allele_support_status"] == "MODEL_COVERAGE_UNASSESSED"
    assert row["allele_extrapolation_status"] == "UNASSESSED"


def test_presentation_persists_blanks_instead_of_sentinels(tmp_path):
    # C1-02/C1-09 end-to-end: an evidence hit with an EMPTY rank cell used
    # to land in the table as "99.0"; scoring already ignored it, so the
    # column now agrees with what was actually scored (blank).
    row = _presentation_fixture(
        tmp_path,
        net_rows="peptide\thla_allele\tnetmhcpan_ba_rank\tnetmhcpan_el_rank\n"
                 "AAGIGILTV\tHLA-A*02:01\t\t\n",
        mhc_rows=(
            "peptide\thla_allele\tmhcflurry_affinity_percentile\t"
            "mhcflurry_processing_score\tmhcflurry_presentation_score\n"
            "AAGIGILTV\tHLA-A*02:01\t3.0\t0.7\t0.8\n"
        ),
    )
    assert row["netmhcpan_ba_rank"] == ""
    assert row["netmhcpan_el_rank"] == ""
    # The remaining predictor still scores: absence excluded, not zeroed in.
    assert row["binding_evidence_score"] == "0.7000"  # 1 - 3.0/10
    assert row["mhcflurry_affinity_percentile"] == "3.0"


def test_presentation_missing_predictor_columns_blank(tmp_path):
    row = _presentation_fixture(
        tmp_path,
        net_rows="peptide\thla_allele\n",
        mhc_rows="peptide\thla_allele\n",
    )
    assert row["netmhcpan_ba_rank"] == ""
    assert row["mhcflurry_affinity_percentile"] == ""
    assert row["mhcflurry_processing_score"] == ""
    assert row["netmhcstabpan_rank"] == ""
    assert row["predictor_allele_support_status"] == ""


# ---------------------------------------------------------------------------
# C1-08: no fabricated restricting allele
# ---------------------------------------------------------------------------

def test_no_prediction_leaves_restricting_hla_unknown(tmp_path):
    row = annotate_variant_peptide_row(
        {"mutant_peptide": "AAGIGILTV", "wildtype_peptide": "SLYYRRSLNV"},
        ["HLA-A*02:01", "HLA-B*07:02"],
        index={},
    )
    # The old code stamped hla_alleles[0] here, inventing a binder.
    assert row["hla_allele"] == ""
    # The requested pool is still preserved, separately, unambiguously.
    assert row["sample_hla_alleles"] == "HLA-A*02:01,HLA-B*07:02"
    # Every dependent lookup stayed blank (C1-04 vocabulary).
    for col in ("netmhcpan_mt_ic50", "netmhcpan_wt_ic50", "prime_wt_score",
                "netmhcstabpan_wt_rank", "bigmhc_im_wt_score",
                "iedb_immunogenicity_score"):
        assert row[col] == "", col


# ---------------------------------------------------------------------------
# C1-06: reuse requires an input fingerprint match
# ---------------------------------------------------------------------------

def _fake_mhcflurry_runner(counter):
    def runner(ctx, out_path: Path) -> Path:
        counter["calls"] += 1
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(f"peptide,allele\nR{counter['calls']},HLA-A*02:01\n", encoding="utf-8")
        return out_path
    return runner


def test_run_tool_reuses_only_while_input_fingerprint_matches(tmp_path, monkeypatch):
    from neoag.tools import runner as runner_mod
    from neoag.tools.registry import RunContext

    peptides = tmp_path / "raw_peptides.tsv"
    peptides.write_text("peptide\thla_allele\nAAGIGILTV\tHLA-A*02:01\n", encoding="utf-8")
    ctx = RunContext(sample_id="S1", outdir=tmp_path / "run", raw_peptides=peptides)
    out_path = tmp_path / "tools" / "mhcflurry.csv"

    counter = {"calls": 0}
    monkeypatch.setitem(runner_mod.RUNNERS, "mhcflurry", _fake_mhcflurry_runner(counter))
    monkeypatch.setattr(
        runner_mod, "check_tool",
        lambda name, exe: SimpleNamespace(available=True, message="fake", name=name),
    )

    # 1. First run executes and records the sidecar fingerprint.
    runner_mod.run_tool("mhcflurry", ctx, out_path)
    assert counter["calls"] == 1
    sidecar = Path(str(out_path) + ".inputs.json")
    assert sidecar.is_file()

    # 2. Identical inputs -> reused, runner not invoked.
    runner_mod.run_tool("mhcflurry", ctx, out_path)
    assert counter["calls"] == 1

    # 3. ONE added peptide (C1-06's requested test) -> fingerprint mismatch,
    #    the stale table is recomputed instead of reused.
    peptides.write_text(
        "peptide\thla_allele\nAAGIGILTV\tHLA-A*02:01\nSIINFEKL\tHLA-B*07:02\n",
        encoding="utf-8",
    )
    runner_mod.run_tool("mhcflurry", ctx, out_path)
    assert counter["calls"] == 2

    # 4. A legacy output with no sidecar at all is never blindly reused.
    sidecar.unlink()
    runner_mod.run_tool("mhcflurry", ctx, out_path)
    assert counter["calls"] == 3

    # 5. The manual escape hatch keeps working.
    monkeypatch.setenv("NEOAG_FORCE_TOOL_RERUN", "1")
    runner_mod.run_tool("mhcflurry", ctx, out_path)
    assert counter["calls"] == 4
