"""Round-3 review C6-02 (FURTHER_SUGGESTIONS WS2): pre-run safety-resource validation.

Acceptance properties exercised here:
- all required resources valid -> release mode PROCEED, with paths, hashes and
  declared builds recorded per resource;
- a required file missing -> exploratory continues with an explicit
  PROCEED_WITH_UNASSESSED; release mode BLOCKs unless a documented waiver
  (resource + reason + authorized_by) applies;
- a file declaring the wrong genome build -> BUILD_MISMATCH, which is a
  correctness defect and is never waiverable;
- a missing normal-junction reference for a run without junction candidates
  -> NOT_APPLICABLE, never a false blocking condition;
- a valid EMPTY table is distinguished from MISSING and from UNREADABLE and
  from PARSE_FAILURE;
- locus coverage over the candidate gene set reported separately (COVERAGE_LOW);
- the gate records the preflight beside its outputs without changing any gate
  row, and the release-integrity audit traces the recorded decision.
"""

import hashlib
import json
import os
from pathlib import Path

import pytest

from neoag.cli import main as cli_main
from neoag.open_neo.review_integrity import audit_review_inputs
from neoag.peptide_safety_gate import build_peptide_safety_gate
from neoag.safety_preflight import (
    SAFETY_PREFLIGHT_SCHEMA_VERSION,
    SAFETY_RELEASE_CRITERIA,
    run_safety_preflight,
    write_safety_preflight,
)


def _rows_by_resource(report):
    return {row["resource"]: row for row in report["rows"]}


def _proteome(tmp_path, name="proteome.fasta"):
    p = tmp_path / name
    p.write_text(">sp|P00001|A_HUMAN\nAAAAAAAAA\n>sp|P00002|B_HUMAN\nBBBBBBBBB\n", encoding="utf-8")
    return p


def _ligands(tmp_path, name="ligands.tsv"):
    p = tmp_path / name
    p.write_text("peptide\thla_allele\nAAAAAAAAA\tHLA-A*02:01\n", encoding="utf-8")
    return p


def _expression(tmp_path, genes, name="expression.tsv"):
    p = tmp_path / name
    p.write_text("gene\tmax_tpm\n" + "".join(f"{g}\t1.0\n" for g in genes), encoding="utf-8")
    return p


def _junctions(tmp_path, name="junctions.tsv"):
    p = tmp_path / name
    p.write_text("junction_id\tcoverage\nchr1:100-200\t9\n", encoding="utf-8")
    return p


def _all_valid(tmp_path):
    return {
        "reference_proteome": _proteome(tmp_path),
        "normal_hla_ligands": _ligands(tmp_path),
        "normal_expression": _expression(tmp_path, ["EML4", "ALK"]),
        "normal_junctions": _junctions(tmp_path),
    }


def test_c6_02_all_valid_resources_pass_release_and_record_provenance(tmp_path):
    paths = _all_valid(tmp_path)
    report = run_safety_preflight(
        paths=paths,
        event_types=["Fusion Fusion"],
        candidate_genes=["EML4", "ALK"],
        mode="release",
    )
    assert report["schema_version"] == SAFETY_PREFLIGHT_SCHEMA_VERSION
    assert report["release_criteria"] == SAFETY_RELEASE_CRITERIA
    assert report["release_decision"] == "PROCEED"
    assert report["blocking_resources"] == ""
    assert report["unresolved_resources"] == ""
    rows = _rows_by_resource(report)
    assert all(row["status"] == "OK" for row in rows.values())
    for name, path in paths.items():
        row = rows[name]
        assert row["path"] == str(path)
        assert row["sha256"] == hashlib.sha256(Path(path).read_bytes()).hexdigest()
    assert rows["normal_junctions"]["scope"] == "junction_events"


def test_c6_02_missing_required_resource_exploratory_tolerates_release_blocks(tmp_path):
    paths = _all_valid(tmp_path)
    missing = tmp_path / "absent.tsv"  # never created
    paths["normal_expression"] = missing

    exploratory = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"], mode="exploratory",
    )
    assert exploratory["release_decision"] == "PROCEED_WITH_UNASSESSED"
    assert "normal_expression" in exploratory["unresolved_resources"]
    row = _rows_by_resource(exploratory)["normal_expression"]
    assert row["status"] == "MISSING"
    assert "UNASSESSED" in row["evidence_consequence"]  # consequence is stated up front

    release = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"], mode="release",
    )
    assert release["release_decision"] == "BLOCK_RELEASE"
    assert release["blocking_resources"] == "normal_expression"


def test_c6_02_documented_waiver_frees_release_but_incomplete_waiver_does_not(tmp_path):
    paths = _all_valid(tmp_path)
    paths["normal_expression"] = None  # not provided at all

    waiver = [{"resource": "normal_expression", "reason": "cohort lacks normal RNA", "authorized_by": "PI"}]
    release = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"],
        mode="release", waivers=waiver,
    )
    assert release["release_decision"] == "PROCEED"
    assert release["blocking_resources"] == ""
    row = _rows_by_resource(release)["normal_expression"]
    assert row["status"] == "WAIVED"
    assert "cohort lacks normal RNA" in row["detail"] and "authorized_by=PI" in row["detail"]

    incomplete = [{"resource": "normal_expression", "reason": "no normal RNA"}]  # no authorized_by
    release2 = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"],
        mode="release", waivers=incomplete,
    )
    assert release2["release_decision"] == "BLOCK_RELEASE"
    assert _rows_by_resource(release2)["normal_expression"]["status"] == "MISSING"


def test_c6_02_wrong_genome_build_is_never_waiverable(tmp_path):
    paths = _all_valid(tmp_path)
    proteome = paths["reference_proteome"]
    Path(proteome).with_name(proteome.name + ".build").write_text("GRCh37\n", encoding="utf-8")

    report = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"],
        genome_build="GRCh38", mode="release",
        waivers=[{"resource": "reference_proteome", "reason": "only build on disk", "authorized_by": "PI"}],
    )
    row = _rows_by_resource(report)["reference_proteome"]
    assert row["status"] == "BUILD_MISMATCH"
    assert row["declared_build"] == "GRCh37"
    assert "waiver REJECTED" in row["detail"]
    assert report["release_decision"] == "BLOCK_RELEASE"
    assert "reference_proteome" in report["blocking_resources"]


def test_c6_02_matching_declared_build_passes_and_undeclared_is_disclosed(tmp_path):
    paths = _all_valid(tmp_path)
    proteome = paths["reference_proteome"]
    Path(proteome).with_name(proteome.name + ".build").write_text("GRCh38\n", encoding="utf-8")
    report = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"],
        genome_build="GRCh38", mode="release",
    )
    rows = _rows_by_resource(report)
    assert rows["reference_proteome"]["status"] == "OK"
    # expression does not declare a build: disclosed, but never blocking
    assert rows["normal_expression"]["status"] == "OK"
    assert "genome build undeclared" in rows["normal_expression"]["detail"]
    assert report["release_decision"] == "PROCEED"


def test_c6_02_inapplicable_missing_junction_reference_never_blocks(tmp_path):
    paths = _all_valid(tmp_path)
    paths["normal_junctions"] = None

    snv_run = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=["EML4"], mode="release",
    )
    row = _rows_by_resource(snv_run)["normal_junctions"]
    assert row["status"] == "NOT_APPLICABLE"
    assert "normal_junctions" not in snv_run["unresolved_resources"]
    assert snv_run["release_decision"] == "PROCEED"

    fusion_run = run_safety_preflight(
        paths=paths, event_types=["Fusion Fusion"], candidate_genes=["EML4", "ALK"], mode="release",
    )
    assert _rows_by_resource(fusion_run)["normal_junctions"]["status"] == "MISSING"
    assert fusion_run["release_decision"] == "BLOCK_RELEASE"
    assert fusion_run["blocking_resources"] == "normal_junctions"


@pytest.mark.skipif(os.geteuid() == 0, reason="chmod 000 is not enforced for root")
def test_c6_02_empty_valid_missing_unreadable_parse_failure_are_distinguished(tmp_path):
    paths = _all_valid(tmp_path)

    paths["normal_expression"] = _expression(tmp_path, [], name="empty_expr.tsv")
    report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
    assert _rows_by_resource(report)["normal_expression"]["status"] == "EMPTY_VALID"
    assert report["release_decision"] == "BLOCK_RELEASE"  # valid but evidentially empty

    paths["normal_expression"] = tmp_path / "absent_expr.tsv"
    report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
    assert _rows_by_resource(report)["normal_expression"]["status"] == "MISSING"

    unreadable = tmp_path / "unreadable.tsv"
    unreadable.write_text("gene\tmax_tpm\nEML4\t1\n", encoding="utf-8")
    unreadable.chmod(0o000)
    try:
        paths["normal_expression"] = unreadable
        report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
        assert _rows_by_resource(report)["normal_expression"]["status"] == "UNREADABLE"
        # a correctness defect: the waiver is rejected
        report = run_safety_preflight(
            paths=paths, event_types=["SNV SNV"], mode="release",
            waivers=[{"resource": "normal_expression", "reason": "r", "authorized_by": "a"}],
        )
        assert _rows_by_resource(report)["normal_expression"]["status"] == "UNREADABLE"
        assert "waiver REJECTED" in _rows_by_resource(report)["normal_expression"]["detail"]
    finally:
        unreadable.chmod(0o644)

    wrong_schema = tmp_path / "wrong_schema.tsv"
    wrong_schema.write_text("sample\tvalue\nS1\t5\n", encoding="utf-8")
    paths["normal_expression"] = wrong_schema
    report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
    expr = _rows_by_resource(report)["normal_expression"]
    assert expr["status"] == "PARSE_FAILURE"
    assert "required column" in expr["detail"]

    bare_fasta = tmp_path / "bare.fasta"
    bare_fasta.write_text("AAAAAAAAA\nBBBBBBBBB\n", encoding="utf-8")  # sequences, no headers
    paths["reference_proteome"] = bare_fasta
    report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
    assert _rows_by_resource(report)["reference_proteome"]["status"] == "PARSE_FAILURE"

    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("", encoding="utf-8")
    paths["reference_proteome"] = empty_fasta
    report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
    assert _rows_by_resource(report)["reference_proteome"]["status"] == "EMPTY_VALID"


def test_c6_02_low_locus_coverage_is_separate_and_waiverable(tmp_path):
    paths = _all_valid(tmp_path)
    paths["normal_expression"] = _expression(tmp_path, ["EML4"])  # only 1 of 3 wanted genes

    wanted = ["EML4", "ALK", "TP53"]
    report = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=wanted, mode="release",
    )
    expr = _rows_by_resource(report)["normal_expression"]
    assert expr["status"] == "COVERAGE_LOW"
    assert expr["gene_coverage"] == "0.3333"
    assert report["release_decision"] == "BLOCK_RELEASE"

    report = run_safety_preflight(
        paths=paths, event_types=["SNV SNV"], candidate_genes=wanted, mode="release",
        waivers=[{"resource": "normal_expression", "reason": "panel cohort", "authorized_by": "PI"}],
    )
    assert _rows_by_resource(report)["normal_expression"]["status"] == "WAIVED"
    assert report["release_decision"] == "PROCEED"

    # without a candidate gene set the coverage axis does not fire
    report = run_safety_preflight(paths=paths, event_types=["SNV SNV"], mode="release")
    expr = _rows_by_resource(report)["normal_expression"]
    assert expr["status"] == "OK"
    assert expr["gene_coverage"] == ""


def test_c6_02_unknown_mode_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="unknown safety preflight mode"):
        run_safety_preflight(paths=_all_valid(tmp_path), mode="turbo")


def test_c6_02_gate_records_preflight_sidecar_without_changing_gate_rows(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    from neoag.utils import write_tsv

    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "KRAS", "event_confidence": "0.9",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "SNV", "mutation_source": "SNV", "gene": "KRAS",
        "peptide": "VVVGAGGVG", "wildtype_peptide": "VVVGAGGVG",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I", "binding_rank": "0.2",
    }])
    out_gate = tmp_path / "peptide_safety.tsv"
    rows, _event_rows = build_peptide_safety_gate(
        raw_events=events, raw_peptides=peptides, out_peptide_safety=out_gate,
    )
    sidecar = tmp_path / "safety_preflight.json"
    assert sidecar.is_file() and (tmp_path / "safety_preflight.tsv").is_file()
    report = json.loads(sidecar.read_text(encoding="utf-8"))
    assert report["mode"] == "exploratory"
    assert report["release_decision"] == "PROCEED_WITH_UNASSESSED"
    expr_rows = _rows_by_resource(report)
    assert expr_rows["reference_proteome"]["status"] == "MISSING"
    # SNV-only run: the junction reference must not appear as a blocking gap
    assert expr_rows["normal_junctions"]["status"] == "NOT_APPLICABLE"
    # recording only: the gate rows are untouched by the preflight
    assert rows[0]["safety_status"] == "SAFETY_PARTIAL"

    # release mode through the profile is still record-only for the gate itself
    rows2, _ = build_peptide_safety_gate(
        raw_events=events, raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety_release.tsv",
        profile={"safety": {"preflight_mode": "release"}},
    )
    report2 = json.loads((tmp_path / "safety_preflight.json").read_text(encoding="utf-8"))
    assert report2["mode"] == "release"
    assert report2["release_decision"] == "BLOCK_RELEASE"
    assert rows2[0]["safety_status"] == "SAFETY_PARTIAL"


def _audit_fixture(root: Path) -> dict:
    scoring = root / "scoring"
    scoring.mkdir(parents=True)
    (scoring / "ranked_peptides.weighted_baseline.tsv").write_text(
        "peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\nP2\t0.7\tC\n", encoding="utf-8"
    )
    (scoring / "ranked_peptides.evidence_consensus.tsv").write_text(
        "evidence_rank\tpeptide_id\tevent_id\tevent_type\tgene\tpeptide\thla_allele\tevidence_grade\tpareto_front\trna_support_state\trna_support_status\tsafety_state\tpresentation_consensus_state\tmutant_specificity_state\tclonality_state\tccf_confidence\thla_appm_state\tevidence_completeness_state\thard_failure\n"
        "1\tP1\tE1\tSNV\tGENE1\tSYFPEITHI\tHLA-A*02:01\tR1\t1\tRNA_CONFIRMED\tRNA_ALT_DETECTED\tSAFETY_PASS\tPRESENTATION_CONSISTENT_STRONG\tMT_SPECIFIC\tCLONAL_LIKE\thigh\tHLA_APPM_RETAINED\tCOMPLETE\tno\n"
        "2\tP2\tE2\tfusion\tGENE2::GENE3\tABCDEFGHI\tHLA-B*07:02\tR3\t1\tRNA_UNASSESSED\tRNA_ALT_NOT_DETECTED\tSAFETY_PARTIAL\tPRESENTATION_SINGLE_TOOL\tNOVEL_JUNCTION\tunresolved\tlow\tHLA_LOH_UNASSESSED\tPARTIAL\tno\n",
        encoding="utf-8",
    )
    (scoring / "all_tool_results.tsv").write_text(
        (scoring / "ranked_peptides.evidence_consensus.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    (scoring / "validation_plan.tsv").write_text(
        "event_id\trecommended_validation\nE1\tMT/WT pair\nE2\ttargeted RNA\n", encoding="utf-8")
    (scoring / "run_manifest.json").write_text(json.dumps(
        {"schema_version": "1.0", "algorithm": "discrete_state_grade_track_pareto_v2",
         "status": "PROVISIONAL_RESEARCH_ONLY", "outputs": {}}) + "\n", encoding="utf-8")
    (scoring / "ranked_events.evidence_consensus.tsv").write_text(
        "event_evidence_rank\tevent_group_id\tevent_id\tgene\tevent_type\tbest_evidence_grade\t"
        "best_pareto_front\tmanual_review_required\trecommended_next_steps\tevent_consensus_trace\t"
        "representative_1_peptide_id\trepresentative_1_peptide\trepresentative_1_hla_allele\n"
        "1\tEVENT:E1\tE1\tGENE1\tSNV\tR1\t1\tno\tMT/WT peptide assay\tcomplete\tP1\tSYFPEITHI\tHLA-A*02:01\n"
        "2\tEVENT:E2\tE2\tGENE2::GENE3\tFusion\tR3\t1\tyes\ttargeted RNA and second caller\tRNA-only\tP2\tABCDEFGHI\tHLA-B*07:02\n",
        encoding="utf-8",
    )
    return {
        "run_manifest": str(scoring / "run_manifest.json"),
        "consensus_events": str(scoring / "ranked_events.evidence_consensus.tsv"),
        "consensus_peptides": str(scoring / "ranked_peptides.evidence_consensus.tsv"),
        "weighted_baseline": str(scoring / "ranked_peptides.weighted_baseline.tsv"),
        "all_tool_results": str(scoring / "all_tool_results.tsv"),
        "validation_plan": str(scoring / "validation_plan.tsv"),
    }


def test_c6_02_release_audit_traces_recorded_preflight_decision(tmp_path):
    artifacts = _audit_fixture(tmp_path / "result")
    scoring = Path(artifacts["run_manifest"]).parent
    sidecar = scoring / "safety_preflight.json"

    def _audit(name):
        result = audit_review_inputs(dict(artifacts), tmp_path / name)
        checks = {row["check"]: row for row in result["checks"]}
        return result, checks

    # absence is disclosed at INFO level: legacy runs must not flip status
    result, checks = _audit("audit-absent")
    assert checks["safety_resource_preflight"]["status"] == "NOT_RECORDED"
    assert checks["safety_resource_preflight"]["severity"] == "INFO"
    baseline_status = result["status"]
    assert baseline_status != "BLOCKED"

    sidecar.write_text(json.dumps({
        "schema_version": SAFETY_PREFLIGHT_SCHEMA_VERSION, "mode": "exploratory",
        "release_decision": "PROCEED_WITH_UNASSESSED",
        "unresolved_resources": "normal_expression", "blocking_resources": "",
    }), encoding="utf-8")
    _result, checks = _audit("audit-unassessed")
    row = checks["safety_resource_preflight"]
    assert row["status"] == "WARN" and "normal_expression" in row["detail"]

    sidecar.write_text(json.dumps({
        "schema_version": SAFETY_PREFLIGHT_SCHEMA_VERSION, "mode": "release",
        "release_decision": "BLOCK_RELEASE",
        "unresolved_resources": "normal_junctions", "blocking_resources": "normal_junctions",
    }), encoding="utf-8")
    result, checks = _audit("audit-blocked")
    assert checks["safety_resource_preflight"]["status"] == "FAIL"
    assert result["status"] == "BLOCKED"


def test_c6_02_cli_reports_rows_and_exit_code(tmp_path):
    paths = _all_valid(tmp_path)
    events = tmp_path / "raw_events.tsv"
    from neoag.utils import write_tsv

    write_tsv(events, [{
        "event_id": "E1", "event_type": "Fusion", "mutation_source": "SV",
        "gene": "EML4::ALK",
    }])
    out_json = tmp_path / "preflight" / "safety_preflight.json"
    base = [
        "safety-preflight", "--mode", "release", "--raw-events", str(events),
        "--reference-proteome", str(paths["reference_proteome"]),
        "--normal-hla-ligands", str(paths["normal_hla_ligands"]),
        "--normal-expression", str(paths["normal_expression"]),
        "--normal-junctions", str(paths["normal_junctions"]),
        "--out-json", str(out_json),
    ]
    cli_main(base)  # valid run exits normally
    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["release_decision"] == "PROCEED"
    # genes were parsed from "EML4::ALK" and coverage computed against them
    expr = {row["resource"]: row for row in report["rows"]}["normal_expression"]
    assert expr["status"] == "OK" and expr["gene_coverage"] == "1.0000"
    assert (tmp_path / "preflight" / "safety_preflight.json.tsv").is_file()  # default tsv name

    blocked = list(base)
    flag = blocked.index("--normal-junctions")
    del blocked[flag:flag + 2]
    with pytest.raises(SystemExit) as excinfo:
        cli_main(blocked)
    assert excinfo.value.code == 2

    waived = blocked + [
        "--waiver", "normal_junctions:no normal RNA library:PI",
        "--out-tsv", str(tmp_path / "preflight" / "waived.tsv"),
    ]
    # --out-json was already consumed above; the waived run reuses the same json path
    cli_main(waived)
    report = json.loads(out_json.read_text(encoding="utf-8"))
    assert report["release_decision"] == "PROCEED"
