"""G-package work package G-f tests (B-E review findings E-01/02/04/05/06/07/08).

Every test here demonstrates the behaviour *after* the fix and, where the
review named a concrete failure mode, the old behaviour is stated in the
assertions or comments so the effect of the change is visible.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from neoag.agent_skills.patient_report import main as patient_report_main
from neoag.benchmark_system import (
    _join_ranked_to_ligandome,
    _labelled_ligand_rows,
    run_appm_multiplier_delta,
)
from neoag.evidence_consensus import _derived_grade_caps
from neoag.open_neo.review import _place_fusion_section, _priority
from neoag.open_neo.review_integrity import audit_review_inputs
from neoag.report_from_final import _purity_qc_status
from neoag.reports_dual import (
    ReportBundle,
    _patient_private_path_pattern,
    _patient_release_audit,
    load_report_bundle,
    make_patient_report,
)
from neoag.scoring import recommended as recommended_text
from neoag.status_words import phasing_needs_resolution
from neoag.utils import read_tsv
from neoag.validation_design import classify_validation_mode


# ---------------------------------------------------------------- E-01 -----

def test_patient_report_skill_escapes_markup_from_data_fields(tmp_path):
    """E-01: TSV-derived text used to be embedded raw into patient HTML."""
    rec = tmp_path / "recommendation.tsv"
    rec.write_text(
        "event_id\tpeptide_id\tpeptide\thla_allele\tgene\tevent_type\tfinal_priority\n"
        "E1\tP1\tAAAAAAAAA\tHLA-A*02:01\t<script>alert(1)</script>\tSNV\tB\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "out"
    rc = patient_report_main([
        "--recommendation", str(rec),
        "--outdir", str(outdir),
        "--title", "报告<img src=x onerror=alert(1)>",
    ])
    assert rc == 0
    html_text = (outdir / "patient_report.html").read_text(encoding="utf-8")
    # No executable markup may survive; the injected payload must appear
    # only as escaped text.
    assert "<script>" not in html_text
    assert "<img" not in html_text
    assert "onerror=alert(1)>" not in html_text.replace("&gt;", "")
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html_text
    assert "&lt;img src=x onerror=alert(1)&gt;" in html_text
    # The markdown sibling keeps the raw string (it is not an HTML sink).
    assert "<script>alert(1)</script>" in (outdir / "patient_report.md").read_text(encoding="utf-8")


# ---------------------------------------------------------------- E-02 -----

def test_private_path_pattern_covers_deployment_roots():
    """E-02: /data1 and /data2 (and /opt, /srv) used to escape the scan."""
    pattern = _patient_private_path_pattern({})
    for leak in (
        "/mnt/data/x.tsv", "/root/x", "/home/u/x",
        "/data1/neoag/final/out.tsv", "/data2/caojingqi/case/out.tsv",
        "/opt/tools/bin/x", "/srv/refs/x.fa",
        "nohup: ignoring input",
        "Traceback (most recent call last):",
    ):
        assert pattern.search(leak), leak
    assert not pattern.search("候选位点中位深度 120x；逻辑文件名 ranked_events.evidence_consensus.tsv")


def test_private_path_pattern_includes_manifest_roots():
    """E-02: roots seen in the run's own provenance are scanned too."""
    pattern = _patient_private_path_pattern(
        {"input_files": {"somatic_vcf": "/weirdroot42/case/t.vcf.gz"}}
    )
    assert pattern.search("/weirdroot42/leaked/path")


def test_release_audit_check10_fails_on_deployment_path(tmp_path):
    bundle = ReportBundle(
        profile={"_profile_name": "t"}, events=[], peptides=[], validation_rows=[],
        provenance={"sample_id": "S1"}, sample_id="S1",
    )
    rows, _status = _patient_release_audit(
        bundle, [], [], {}, {}, "see /data2/caojingqi/production/final/rank_peptides.tsv",
    )
    check10 = next(row for row in rows if row["编号"] == "10")
    assert check10["状态"] == "失败"
    rows_ok, _ = _patient_release_audit(bundle, [], [], {}, {}, "逻辑文件名 ranked_events.tsv")
    assert next(row for row in rows_ok if row["编号"] == "10")["状态"] == "通过"


# ---------------------------------------------------------------- E-04 -----

def _integrity_artifacts(tmp_path, grades, cap_reasons):
    scoring = tmp_path / "result" / "scoring"
    scoring.mkdir(parents=True)
    peptide_header = (
        "peptide_id\tevent_id\tevidence_grade\thard_failure\t"
        "netmhcpan_mt_rank_el\tmhcflurry_presentation_score\tevidence_grade_cap_reasons\n"
    )
    peptide_rows = "".join(
        f"P{index}\tE{index}\t{grades[index]}\tno\t0.5\t0.9\t{cap_reasons}\n"
        for index in range(len(grades))
    )
    peptide_path = scoring / "ranked_peptides.evidence_consensus.tsv"
    peptide_path.write_text(peptide_header + peptide_rows, encoding="utf-8")
    (scoring / "all_tool_results.tsv").write_text(peptide_header + peptide_rows, encoding="utf-8")
    (scoring / "ranked_events.evidence_consensus.tsv").write_text("event_id\n", encoding="utf-8")
    (scoring / "ranked_peptides.weighted_baseline.tsv").write_text("peptide_id\n", encoding="utf-8")
    (scoring / "validation_plan.tsv").write_text("peptide_id\n", encoding="utf-8")
    (scoring / "run_manifest.json").write_text(json.dumps({"run_id": "R1"}), encoding="utf-8")
    return {
        "run_manifest": str(scoring / "run_manifest.json"),
        "consensus_events": str(scoring / "ranked_events.evidence_consensus.tsv"),
        "consensus_peptides": str(peptide_path),
        "weighted_baseline": str(scoring / "ranked_peptides.weighted_baseline.tsv"),
        "all_tool_results": str(scoring / "all_tool_results.tsv"),
        "validation_plan": str(scoring / "validation_plan.tsv"),
    }


def test_integrity_warns_on_no_top_grade_with_missing_asset_caps(tmp_path):
    """E-04: the all-R3 / no-R1-R2 collapse used to pass silently."""
    artifacts = _integrity_artifacts(
        tmp_path, ["R3"] * 20, "CAP_SAFETY_PARTIAL,CAP_NORMAL_JUNCTION_UNASSESSED",
    )
    result = audit_review_inputs(artifacts, tmp_path / "audit")
    checks = {row["check"]: row for row in result["checks"]}
    concentration = checks["no_top_grade_concentration"]
    assert concentration["status"] == "WARN"
    assert "no R1/R2 candidates" in concentration["detail"]
    assert "capped=20/20" in concentration["detail"]
    assert "capped_by_absent_resource=20 (100%)" in concentration["detail"]
    # A warning must never be a hard failure: it triggers review, not rejection.
    assert "no_top_grade_concentration" not in {row["check"] for row in result["blocking_issues"]}
    # The all-R4 alert still exists unchanged next to the new check.
    assert checks["evidence_grade_distribution"]["status"] == "PASS"


def test_integrity_top_grade_check_passes_when_one_r1_present(tmp_path):
    artifacts = _integrity_artifacts(
        tmp_path, ["R1"] + ["R3"] * 19, "CAP_NORMAL_JUNCTION_UNASSESSED",
    )
    result = audit_review_inputs(artifacts, tmp_path / "audit")
    checks = {row["check"]: row for row in result["checks"]}
    assert checks["no_top_grade_concentration"]["status"] == "PASS"


# ---------------------------------------------------------------- E-05 -----

@pytest.mark.parametrize("value", [
    "PHASING_REQUIRED",
    "UNPHASED",
    "PHASE_UNRESOLVED",
    "PHASING_UNRESOLVED",
    "LOW_CONFIDENCE",
    "PHASED_CIS LOW_CONFIDENCE",
])
def test_phasing_helper_fires_on_shared_tokens(value):
    assert phasing_needs_resolution(value), value


@pytest.mark.parametrize("value", [
    "",
    "PHASED_CIS",
    "PHASED_NOT_CIS",
    "PHASED_NOT_CIS low",
    "PHASING_NOT_REQUIRED",
    "NOT_LOW_CONFIDENCE",
    "REPAIR_REQUIRED",
])
def test_phasing_helper_stays_silent_on_resolved_or_negated(value):
    # "PHASED_NOT_CIS low" documents why bare LOW is not a token: every
    # resolved trans pair carries confidence "low".
    assert not phasing_needs_resolution(value), value


@pytest.mark.parametrize("phasing_value", [
    "UNPHASED",
    "PHASING_REQUIRED",
    "PHASE_UNRESOLVED",
    "LOW_CONFIDENCE",
])
def test_d3_caps_and_review_gate_share_one_phasing_vocabulary(phasing_value):
    """E-05: UNPHASED used to escape review, UNRESOLVED escaped the D3 cap."""
    states = {
        "clonality": {"state": "CLONAL", "assessed": True},
        "safety": {"state": "SAFETY_PASS"},
        "hla_appm": {"state": "HLA_LOH_PASSED"},
        "mutant_specificity": {"state": "MT_SPECIFIC"},
        "presentation_consensus": {"state": "PRESENTATION_CONCORDANT"},
        "event_authenticity": {"state": "EVENT_STRONG"},
        "rna_support": {"state": "RNA_CONFIRMED"},
    }
    caps = _derived_grade_caps({"haplotype_status": phasing_value}, states, "MISSENSE")
    assert ("R3", "CAP_PHASING_REQUIRED") in caps, phasing_value
    assert _priority("R2", {"phase_status": phasing_value})[:2] == (
        "COMPLETE_EVIDENCE", "PHASING_FIRST",
    ), phasing_value
    # The scoring layer and validation planner use the same vocabulary:
    assert "read-backed phasing" in str(recommended_text(
        {"safety_status": "PASS", "haplotype_status": phasing_value}, 1.0, 1.0,
    ))
    assert classify_validation_mode({
        "evidence_grade": "R2", "priority": "B", "safety_status": "PASS",
        "haplotype_status": phasing_value,
    }) == "phasing_required", phasing_value


def test_resolved_phasing_not_gated_by_any_stage():
    states = {
        "clonality": {"state": "CLONAL", "assessed": True},
        "safety": {"state": "SAFETY_PASS"},
        "hla_appm": {"state": "HLA_LOH_PASSED"},
        "mutant_specificity": {"state": "MT_SPECIFIC"},
        "presentation_consensus": {"state": "PRESENTATION_CONCORDANT"},
        "event_authenticity": {"state": "EVENT_STRONG"},
        "rna_support": {"state": "RNA_CONFIRMED"},
    }
    caps = _derived_grade_caps(
        {"haplotype_status": "PHASED_CIS", "phase_confidence": "high"}, states, "MISSENSE",
    )
    assert ("R3", "CAP_PHASING_REQUIRED") not in caps
    assert _priority("R2", {"phase_status": "PHASED_CIS"})[1] != "PHASING_FIRST"


# ---------------------------------------------------------------- E-06 -----

def test_purity_qc_status_distinguishes_missing_qc_from_assessed(tmp_path):
    """E-06: absence of a QC file used to be reported as ASSESSED."""
    missing = tmp_path / "nothing"
    assert _purity_qc_status("PURPLE", missing / "out") == "RESULT_PATH_MISSING"

    purple_dir = tmp_path / "purple"
    purple_dir.mkdir()
    assert _purity_qc_status("PURPLE", purple_dir) == "QC_NOT_PROVIDED"

    marker = purple_dir / "sample.purple.qc"
    marker.write_text("something\telse\n", encoding="utf-8")
    assert _purity_qc_status("PURPLE", purple_dir) == "QC_MARKER_WITHOUT_STATUS"

    marker.write_text("sample\tvalue\nQCStatus\tPASS\n", encoding="utf-8")
    assert _purity_qc_status("PURPLE", purple_dir) == "PASS"

    facets_dir = tmp_path / "facets"
    facets_dir.mkdir()
    assert _purity_qc_status("FACETS", facets_dir) == "QC_NOT_PROVIDED"
    (facets_dir / "facets_omni2p5_summary.tsv").write_text(
        "metric\tvalue\nstatus\tPASS\n", encoding="utf-8",
    )
    assert _purity_qc_status("FACETS", facets_dir) == "PASS"

    other_dir = tmp_path / "sequenza"
    other_dir.mkdir()
    (other_dir / "seq.txt").write_text("purity 0.5\n", encoding="utf-8")
    assert _purity_qc_status("SEQUENZA", other_dir) == "QC_NOT_PROVIDED"


# ---------------------------------------------------------------- E-07 -----

def _audit_bundle(provenance):
    return ReportBundle(
        profile={"_profile_name": "requested_profile"},
        events=[], peptides=[], validation_rows=[],
        provenance=provenance, sample_id="S1",
    )


def _release_audit(path):
    return json.loads(path.with_suffix(".release_audit.json").read_text(encoding="utf-8"))


def test_release_audit_fails_on_silent_profile_fallback(tmp_path):
    """E-07: a fallback profile must block the formal release, loudly."""
    bundle = _audit_bundle({
        "sample_id": "S1",
        "profile_status": {
            "requested": "sarcoma_rna_supported_v2_provisional",
            "status": "FALLBACK_DEFAULTS",
            "error": "FileNotFoundError: profile not found",
        },
    })
    out = tmp_path / "patient.html"
    make_patient_report(out, bundle)
    payload = _release_audit(out)
    check12 = next(row for row in payload["checks"] if row["编号"] == "12")
    assert check12["状态"] == "失败"
    assert "FALLBACK_DEFAULTS" in check12["说明"]
    assert payload["status"] != "PASS"
    # E-08c: the release banner is anchored right before the first section
    # instead of a fixed list index, so it stays above the fold.
    html_text = out.read_text(encoding="utf-8")
    assert html_text.index("发布状态：") < html_text.index("class='section'")


def test_release_audit_passes_check12_when_profile_loaded(tmp_path):
    bundle = _audit_bundle({
        "sample_id": "S1",
        "profile_status": {"requested": "p1", "status": "LOADED"},
    })
    out = tmp_path / "patient.html"
    make_patient_report(out, bundle)
    payload = _release_audit(out)
    check12 = next(row for row in payload["checks"] if row["编号"] == "12")
    assert check12["状态"] == "通过"


def test_release_audit_check12_passes_when_status_not_recorded(tmp_path):
    # Reports rebuilt without going through write_reports_from_final have no
    # profile_status; that keeps the historical behaviour (not recorded).
    out = tmp_path / "patient.html"
    make_patient_report(out, _audit_bundle({"sample_id": "S1"}))
    payload = _release_audit(out)
    check12 = next(row for row in payload["checks"] if row["编号"] == "12")
    assert check12["状态"] == "通过"
    assert "NOT_RECORDED" in check12["说明"]


def test_rebuild_labels_profile_fallback_in_provenance(tmp_path, monkeypatch):
    """E-07 end-to-end: the rebuild path labels the stub profile."""
    import neoag.report_from_final as rff

    production = tmp_path / "prod"
    final = production / "final"
    (final / "scoring").mkdir(parents=True)
    (final / "scoring" / "ranked_events.evidence_consensus.tsv").write_text(
        "event_id\nE1\n", encoding="utf-8",
    )
    (final / "scoring" / "ranked_peptides.evidence_consensus.tsv").write_text(
        "peptide_id\tevent_id\nP1\tE1\n", encoding="utf-8",
    )

    calls = {"n": 0}

    def boom(name):
        calls["n"] += 1
        raise RuntimeError("profile missing on this host")

    monkeypatch.setattr(rff, "load_profile", boom)
    result = rff.write_reports_from_final(final)
    assert calls["n"] == 1
    provenance = json.loads(Path(result["provenance"]).read_text(encoding="utf-8"))
    status = provenance["profile_status"]
    assert status["status"] == "FALLBACK_DEFAULTS"
    assert status["requested"] == "sarcoma_rna_supported_v2_provisional"
    assert "profile missing on this host" in status["error"]
    payload = _release_audit(Path(result["evidence_report_patient"]))
    check12 = next(row for row in payload["checks"] if row["编号"] == "12")
    assert check12["状态"] == "失败"

    monkeypatch.setattr(rff, "load_profile", lambda name: {"gates": {}})
    result = rff.write_reports_from_final(final)
    provenance = json.loads(Path(result["provenance"]).read_text(encoding="utf-8"))
    assert provenance["profile_status"]["status"] == "LOADED"


# ---------------------------------------------------------------- E-08 -----

@pytest.mark.parametrize("markup,expected_placement", [
    ("<html><body><p>x</p></body></html>", "INSERTED_BEFORE_BODY_END"),
    ("<html><body onload='x'><p>x</p></BODY>", "INSERTED_BEFORE_BODY_END"),
    ("<html><body><p>x</p>\n</body>\n</html>", "INSERTED_BEFORE_BODY_END"),
    ("<html><div>x</div></html>", "INSERTED_BEFORE_HTML_END"),
    ("<div>fragment without closing tags</div>", "APPENDED_NO_CLOSING_TAG"),
])
def test_fusion_section_placement_survives_markup_drift(markup, expected_placement):
    """E-08: the old replace("</body></html>") silently dropped the section."""
    section = "<div class='section'><h2>融合基因检测结果</h2></div>"
    updated, placement = _place_fusion_section(markup, section)
    assert placement == expected_placement
    assert section in updated
    if expected_placement == "INSERTED_BEFORE_BODY_END":
        assert updated.index(section) < updated.lower().rindex("</body>")
    elif expected_placement == "INSERTED_BEFORE_HTML_END":
        assert updated.index(section) < updated.lower().rindex("</html>")
    else:
        assert updated.startswith(markup) and updated.endswith(section)


def _purity_bundle(tmp_path, tools):
    return load_report_bundle(
        profile={"_profile_name": "t"},
        events=[], peptides=[],
        outdir=tmp_path,
        provenance={"sample_id": "S1", "purity_cnv_tools": tools},
    )


def test_multi_tool_purity_fallback_uses_median_not_first(tmp_path):
    """E-08: the multi-tool fallback used to hand back purities[0]."""
    bundle = _purity_bundle(tmp_path, [
        {"tool": "PURPLE", "purity": "0.20", "ploidy": "2.8", "status": "FOUND"},
        {"tool": "FACETS", "purity": "0.80", "ploidy": "3.1", "status": "FOUND"},
        {"tool": "SEQUENZA", "purity": "0.55", "ploidy": "", "status": "FOUND"},
    ])
    consensus = bundle.provenance["purity_cnv_consensus"]
    assert consensus["recommended_purity"] == "0.55"
    assert consensus["status"] == "MULTI_TOOL_ARBITRATION_REQUIRED"
    assert "中位数" in consensus["basis"] and "不静默选择单一工具" in consensus["basis"]

    concordant = _purity_bundle(tmp_path, [
        {"tool": "PURPLE", "purity": "0.50", "ploidy": "", "status": "FOUND"},
        {"tool": "FACETS", "purity": "0.52", "ploidy": "", "status": "FOUND"},
    ]).provenance["purity_cnv_consensus"]
    assert concordant["recommended_purity"] == "0.52"
    assert concordant["status"] == "MULTI_TOOL_REVIEW"


def test_pairing_status_records_its_evidence_level(tmp_path):
    """E-08: pairing statements must state what evidence they rest on."""
    empty = load_report_bundle(
        profile={"_profile_name": "t"}, events=[], peptides=[],
        outdir=tmp_path, provenance={"sample_id": "S1"},
    )
    assert empty.provenance["pairing_evidence_level"] == "NO_PAIRING_EVIDENCE"

    with_depth = load_report_bundle(
        profile={"_profile_name": "t"}, events=[], peptides=[],
        outdir=tmp_path,
        provenance={"sample_id": "S1", "normal_dna_depth": "候选位点中位深度 120x（n=50）"},
    )
    assert with_depth.provenance["pairing_evidence_level"] == "NORMAL_DNA_DEPTH_EVIDENCE"
    assert "BAM Matcher" in with_depth.provenance["pairing_status"]


def test_join_ranked_keeps_measured_zero_and_fills_only_missing(tmp_path):
    """E-08: a genuine 0.0 ligandome score used to be replaced by the ranked score."""
    lig = tmp_path / "ligandome.tsv"
    lig.write_text(
        "peptide_id\tpeptide\thla_allele\tobserved_label\tefficacy_score\n"
        "P1\tAAAA\tHLA-A*02:01\t1\t0.0\n"
        "P2\tBBBB\tHLA-B*07:02\t0\t\n",
        encoding="utf-8",
    )
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text(
        "peptide_id\tpeptide\thla_allele\tefficacy_score\tappm_multiplier\n"
        "P1\tAAAA\tHLA-A*02:01\t0.9\t2.0\n"
        "P2\tBBBB\tHLA-B*07:02\t0.4\t\n",
        encoding="utf-8",
    )
    labelled, status = _labelled_ligand_rows(lig)
    assert status == "completed"
    joined = _join_ranked_to_ligandome(ranked, labelled)
    by_id = {row["peptide_id"]: row for row in joined}
    assert by_id["P1"]["_score"] == "0.0"   # measured zero survives
    assert by_id["P2"]["_score"] == "0.4"   # missing score falls back to ranked


def test_multiplier_delta_keeps_zero_multiplier_visible(tmp_path):
    """E-08: `or 1.0` silently rewrote a real 0.0 multiplier as neutral."""
    lig = tmp_path / "ligandome.tsv"
    lig.write_text(
        "peptide_id\tpeptide\thla_allele\tobserved_label\tefficacy_score\n"
        "P1\tAAAA\tHLA-A*02:01\t1\t0.8\n"
        "P2\tBBBB\tHLA-B*07:02\t0\t0.2\n",
        encoding="utf-8",
    )
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text(
        "peptide_id\tpeptide\thla_allele\tefficacy_score\tappm_multiplier\n"
        "P1\tAAAA\tHLA-A*02:01\t0.8\t0.0\n"
        "P2\tBBBB\tHLA-B*07:02\t0.2\t2.0\n",
        encoding="utf-8",
    )
    out = run_appm_multiplier_delta(outdir=tmp_path, ligandome_ms=lig, ranked_peptides=ranked)
    rows = read_tsv(out)
    assert rows and all(row["benchmark_status"] == "completed" for row in rows)
    assert all("non-positive APPM multiplier" in row["interpretation"] for row in rows)
    assert all("1 row(s)" in row["interpretation"] for row in rows)
