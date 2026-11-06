"""Package G-c tests (B-E review C2-01/03/04/06, C5-01/02/05), 2026-11-06.

C2-01  prime component is rank-only; score-only rows stay unassessed.
C2-03  evidence merges use canonical (peptide, HLA) keys on both sides.
C2-04  existing *_evidence.tsv files are reused only with a matching input
       fingerprint sidecar; stale evidence is re-run or registered missing.
C2-06  neutral-sentinel matching is numeric; mode="mean" is arithmetic.
C5-01  only loss-flavoured LOH states cap peptides; UNASSESSED /
       SINGLE_TOOL_RETAINED record a reason without penalising.
C5-02  provided-but-unmapped proteomics/ligandome inputs are labelled
       input_only_coverage_unresolved instead of counted as gene evidence.
C5-05  paths advertises explainability outputs only after they exist.
"""

from pathlib import Path

import pytest

from neoag.appm_v2 import build_appm_2
from neoag.evidence_provenance import (
    ProvenanceRegistry,
    STATUS_MISSING,
    STATUS_REAL,
)
from neoag.immune_escape import build_immune_escape_evidence
from neoag.immunogenicity_composite import (
    _evidence_fingerprint,
    _evidence_fingerprint_matches,
    _is_neutral_or_blank,
    _write_evidence_fingerprint,
    apply_immunogenicity_evidence,
    combine_component_scores,
    component_scores,
    run_immunogenicity_predictors,
)
from neoag.tools.registry import RunContext
from neoag.utils import read_tsv, write_tsv


def _profile(sources=None, **immuno):
    cfg = {"enabled": True, "sources": sources or ["prime", "bigmhc_im"]}
    cfg.update(immuno)
    return {"immunogenicity": cfg}


# --------------------------------------------------------------------------
# C2-03: normalized merge keys
# --------------------------------------------------------------------------

def test_merge_joins_across_hla_notation_and_case(tmp_path):
    prime = tmp_path / "prime_evidence.tsv"
    prime.write_text(
        "sample_id\tpeptide\thla_allele\tprime_score\tprime_rank\tevidence_status\n"
        "S\tMTPEPTID\tHLA-A*02:01\t0.90\t1.0\treal\n"
        "S\tWTPEPTID\tHLA-A*02:01\t0.20\t50.0\treal\n",
        encoding="utf-8",
    )
    bigmhc = tmp_path / "bigmhc_im_evidence.tsv"
    bigmhc.write_text(
        "sample_id\tpeptide\thla_allele\tbigmhc_im_score\tevidence_status\n"
        "S\tMTPEPTID\tB*58:01\t0.80\treal\n",
        encoding="utf-8",
    )
    rows = [{
        "peptide": "mtpeptid",
        "wildtype_peptide": "wtpeptid",
        "hla_allele": "A*02:01",
    }]
    apply_immunogenicity_evidence(
        rows, {"prime": prime, "bigmhc_im": bigmhc, "deepimmuno": None, "iedb": None},
        _profile(),
    )
    assert rows[0]["prime_score"] == "0.90"
    assert rows[0]["prime_wt_score"] == "0.20"
    # Same row asked with a second allele (B*58:01 written vs HLA-B*58:01
    # typed) still merges; before C2-03 both sides joined raw strings.
    rows_b = [{"peptide": "MTPEPTID", "hla_allele": "HLA-B*58:01"}]
    apply_immunogenicity_evidence(
        rows_b, {"prime": None, "bigmhc_im": bigmhc, "deepimmuno": None, "iedb": None},
        _profile(),
    )
    assert rows_b[0]["bigmhc_im_score"] == "0.80"


def test_iedb_merge_joins_across_notation(tmp_path):
    iedb = tmp_path / "iedb_immunogenicity.tsv"
    iedb.write_text(
        "peptide\thla_allele\tiedb_immunogenicity_score\n"
        "AAAAAAAAA\tHLA-A*02:01\t0.30000\n",
        encoding="utf-8",
    )
    rows = [{"peptide": "aaaaaaaaa", "hla_allele": "A*02:01"}]
    apply_immunogenicity_evidence(
        rows, {"prime": None, "bigmhc_im": None, "deepimmuno": None, "iedb": iedb},
        _profile(sources=["iedb"]),
    )
    assert rows[0]["iedb_immunogenicity_score"] == "0.30000"


# --------------------------------------------------------------------------
# C2-06: numeric neutral sentinel + strict mean mode
# --------------------------------------------------------------------------

@pytest.mark.parametrize("raw,neutral", [
    ("", True), ("0.5", True), ("0.5000", True), ("0.50000", True),
    ("0.50001", False), ("0.4", False), ("abc", True),
])
def test_neutral_sentinel_is_numeric_not_string_form(raw, neutral):
    assert _is_neutral_or_blank(raw) is neutral


def test_formatted_neutral_iedb_does_not_enter_composite():
    # The heuristic writes "%.5f": its neutral output used to slip through the
    # {"", "0.5", "0.5000"} string sentinel and count as a real component.
    profile = _profile(sources=["iedb"])
    assert component_scores({"iedb_immunogenicity_score": "0.50000"}, profile) == {}
    parts = component_scores({"iedb_immunogenicity_score": "0.25000"}, profile)
    assert parts == {"iedb": pytest.approx(0.5)}


def test_mean_mode_is_arithmetic_even_with_hla_multipliers():
    profile = {"immunogenicity": {
        "enabled": True,
        "composite": "mean",
        "hla_source_multipliers": {"HLA-A*02:01": {"prime": 2.0, "bigmhc_im": 1.0}},
    }}
    parts = {"prime": 0.2, "bigmhc_im": 0.8}
    # Old behaviour with these multipliers: weighted (2*0.2+0.8)/3 = 0.4.
    assert combine_component_scores(parts, profile, hla_allele="HLA-A*02:01") == pytest.approx(0.5)
    weighted_profile = {"immunogenicity": {
        "enabled": True,
        "composite": "weighted",
        "hla_source_multipliers": {"HLA-A*02:01": {"prime": 2.0, "bigmhc_im": 1.0}},
    }}
    assert combine_component_scores(parts, weighted_profile, hla_allele="HLA-A*02:01") == pytest.approx(0.4)


def test_zero_multiplier_still_excludes_source_in_mean_mode():
    profile = {"immunogenicity": {
        "enabled": True,
        "composite": "mean",
        "hla_source_multipliers": {"HLA-A*02:01": {"bigmhc_im": 0.0}},
    }}
    parts = {"prime": 0.2, "bigmhc_im": 0.8}
    assert combine_component_scores(parts, profile, hla_allele="HLA-A*02:01") == pytest.approx(0.2)


# --------------------------------------------------------------------------
# C2-04: fingerprint-guarded evidence reuse
# --------------------------------------------------------------------------

def _prime_tool_table(path: Path, rows=("MTPEPTID\tHLA-A*02:01\t0.90\t1.0",)):
    body = "".join(f"{r}\n" for r in rows)
    path.write_text(
        "peptide\thla_allele\tprime_score\tprime_rank\n" + body, encoding="utf-8",
    )


def _peptides(path: Path, peptide="MTPEPTID"):
    path.write_text(
        "peptide_id\tevent_id\tpeptide\thla_allele\n"
        f"P1\tE1\t{peptide}\tHLA-A*02:01\n",
        encoding="utf-8",
    )


def test_fingerprint_helpers_require_sidecar(tmp_path):
    peptides = tmp_path / "raw_peptides.tsv"
    _peptides(peptides)
    evidence = tmp_path / "prime_evidence.tsv"
    evidence.write_text("peptide\nX\n", encoding="utf-8")
    fp = _evidence_fingerprint(peptides, "prime")
    assert fp["tool"] == "prime" and len(fp["peptide_input_sha256"]) == 64
    assert not _evidence_fingerprint_matches(evidence, fp)  # absent sidecar = no reuse
    _write_evidence_fingerprint(evidence, fp)
    assert _evidence_fingerprint_matches(evidence, fp)
    _peptides(peptides, peptide="OTHERPEP")
    assert not _evidence_fingerprint_matches(evidence, _evidence_fingerprint(peptides, "prime"))


def test_evidence_reused_only_while_peptide_input_matches(tmp_path):
    peptides = tmp_path / "raw_peptides.tsv"
    _peptides(peptides)
    outdir = tmp_path / "run"
    ctx = RunContext(sample_id="S", outdir=outdir, raw_peptides=peptides)
    tools_dir = outdir / "tools"
    tools_dir.mkdir(parents=True, exist_ok=True)
    tool_tsv = tools_dir / "prime.tsv"
    _prime_tool_table(tool_tsv)
    # What run_prime/run_tool leaves behind for a completed run (C1-06 sidecar).
    _write_evidence_fingerprint(tool_tsv, _evidence_fingerprint(peptides, "prime"))

    reg = ProvenanceRegistry()
    paths = run_immunogenicity_predictors(
        peptides, outdir, _profile(sources=["prime"]), ctx, skip=True, provenance_registry=reg,
    )
    evidence = outdir / "presentation" / "prime_evidence.tsv"
    assert paths["prime"] == evidence
    sidecar = Path(str(evidence) + ".inputs.json")
    assert sidecar.is_file(), "fresh evidence write must record its input fingerprint"
    assert reg._records["prime"].status == STATUS_REAL

    # Same peptides: reuse is accepted.
    reg2 = ProvenanceRegistry()
    paths2 = run_immunogenicity_predictors(
        peptides, outdir, _profile(sources=["prime"]), ctx,
        skip=True, reuse_sources={"prime"}, provenance_registry=reg2,
    )
    assert paths2["prime"] == evidence
    assert reg2._records["prime"].mode == "passthrough"

    # Peptide input changed: the sidecar no longer matches, the stale file is
    # refused and prime is registered missing (nothing to recompute from).
    _peptides(peptides, peptide="OTHERPEP")
    reg3 = ProvenanceRegistry()
    paths3 = run_immunogenicity_predictors(
        peptides, outdir, _profile(sources=["prime"]), ctx,
        skip=True, reuse_sources={"prime"}, provenance_registry=reg3,
    )
    assert paths3["prime"] is None
    assert reg3._records["prime"].status == STATUS_MISSING


def test_reuse_decline_recomputes_when_possible(tmp_path, monkeypatch):
    peptides = tmp_path / "raw_peptides.tsv"
    _peptides(peptides)
    outdir = tmp_path / "run"
    ctx = RunContext(sample_id="S", outdir=outdir, raw_peptides=peptides)
    tools_dir = outdir / "tools"
    tools_dir.mkdir(parents=True, exist_ok=True)
    _prime_tool_table(tools_dir / "prime.tsv")
    _write_evidence_fingerprint(tools_dir / "prime.tsv", _evidence_fingerprint(peptides, "prime"))
    run_immunogenicity_predictors(
        peptides, outdir, _profile(sources=["prime"]), ctx, skip=True,
    )

    _peptides(peptides, peptide="NEWPEPTID")
    tool_tsv = tools_dir / "prime.tsv"

    def fake_run_prime(run_ctx, out_path):
        _prime_tool_table(Path(out_path), rows=("NEWPEPTID\tHLA-A*02:01\t0.50\t2.0",))
        # Mirror what run_prime/run_tool records after a completed run.
        _write_evidence_fingerprint(Path(out_path), _evidence_fingerprint(peptides, "prime"))

    monkeypatch.setattr("neoag.tools.runner.run_prime", fake_run_prime)
    monkeypatch.setattr(
        "neoag.immunogenicity_composite._immunogenicity_tool_available",
        lambda name, run_ctx: True,
    )
    reg = ProvenanceRegistry()
    paths = run_immunogenicity_predictors(
        peptides, outdir, _profile(sources=["prime"]), ctx,
        reuse_sources={"prime"}, provenance_registry=reg,
    )
    evidence = outdir / "presentation" / "prime_evidence.tsv"
    assert paths["prime"] == evidence
    assert reg._records["prime"].status == STATUS_REAL
    merged_rows = [r for r in read_tsv(evidence) if r.get("peptide") == "NEWPEPTID"]
    assert merged_rows and merged_rows[0]["prime_score"] == "0.50"
    assert _evidence_fingerprint_matches(evidence, _evidence_fingerprint(peptides, "prime"))


# --------------------------------------------------------------------------
# C5-01: absence states no longer cap peptides
# --------------------------------------------------------------------------

def _escape_flags(tmp_path, peptide_rows, loh_rows):
    peptides = tmp_path / "peptides.tsv"
    write_tsv(peptides, peptide_rows)
    paths = {"outdir": tmp_path / "out"}
    kwargs = {}
    if loh_rows is not None:
        loh = tmp_path / "loh.tsv"
        write_tsv(loh, loh_rows)
        kwargs["hla_loh_tsv"] = loh
    result = build_immune_escape_evidence(
        sample_id="S", raw_peptides=peptides, outdir=paths["outdir"], **kwargs,
    )
    return {r["peptide_id"]: r for r in read_tsv(result["peptide_escape_flags"])}


def _row(pid, allele="HLA-A*02:01"):
    return {
        "peptide_id": pid, "event_id": pid, "peptide": "AAAAAAAAA",
        "hla_allele": allele, "mhc_class": "I",
    }


def test_no_loh_input_does_not_cap_peptides(tmp_path):
    flags = _escape_flags(tmp_path, [_row("P1")], None)
    assert flags["P1"]["escape_status"] == "ESCAPE_PASS"
    assert flags["P1"]["priority_cap"] == ""
    assert flags["P1"]["escape_reason"] == "restricting_hla_unassessed"


def test_allele_missing_from_loh_table_does_not_cap(tmp_path):
    loh = [{"hla_allele": "HLA-B*07:02", "consensus_status": "CONSENSUS_RETAINED"}]
    flags = _escape_flags(tmp_path, [_row("P1")], loh)
    assert flags["P1"]["escape_status"] == "ESCAPE_PASS"
    assert flags["P1"]["priority_cap"] == ""
    assert "restricting_hla_unassessed" in flags["P1"]["escape_reason"]


@pytest.mark.parametrize("status", ["SINGLE_TOOL_LOST", "DISCORDANT"])
def test_loss_flavored_loh_still_caps(tmp_path, status):
    loh = [{"hla_allele": "HLA-A*02:01", "consensus_status": status}]
    flags = _escape_flags(tmp_path, [_row("P1")], loh)
    assert flags["P1"]["escape_status"] == "HLA_LOH_INCOMPLETE"
    assert flags["P1"]["priority_cap"] == "C_CAUTION"
    assert flags["P1"]["escape_reason"] == f"restricting_hla_{status.lower()}"


# --------------------------------------------------------------------------
# C5-02: provided assay inputs become input_only, never fake gene evidence
# --------------------------------------------------------------------------

def _appm(tmp_path, *, proteomics=None, ligandome=None):
    kwargs = {}
    if proteomics is not None:
        p = tmp_path / "proteomics.tsv"
        write_tsv(p, proteomics)
        kwargs["proteomics_tsv"] = p
    if ligandome is not None:
        l = tmp_path / "ligandome.tsv"
        write_tsv(l, ligandome)
        kwargs["hla_ligandome_tsv"] = l
    out = tmp_path / "appm"
    paths = build_appm_2(sample_id="S", outdir=out, **kwargs)
    completeness = read_tsv(paths["appm_evidence_completeness"])[0]
    summary = read_tsv(paths["appm_summary"])[0]
    genes = {r["gene"]: r for r in read_tsv(paths["appm_gene_status"])}
    return paths, completeness, summary, genes


def test_unmapped_proteomics_input_is_labelled_input_only(tmp_path):
    paths, completeness, summary, genes = _appm(
        tmp_path, proteomics=[{"gene": "B2M", "protein": "present"}],
    )
    assert genes["B2M"]["protein_status"] == "input_only_coverage_unresolved"
    assert genes["B2M"]["ligandome_support"] == "not_assessed"
    assert completeness["protein_assessed"] == "no"
    assert completeness["functional_validation_status"] == "input_only_coverage_unresolved"
    assert completeness["validation_evidence_source"] == "protein_provided_without_gene_mapping"
    assert "protein_gene_mapping" in completeness["missing_evidence"].split(";")
    assert summary["claim_strength"] == "no_gene_level_functional_validation"
    assert summary["functional_validation_status"] == "input_only_coverage_unresolved"


def test_unmapped_ligandome_input_is_labelled_input_only(tmp_path):
    _, completeness, _, genes = _appm(
        tmp_path, ligandome=[{"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01"}],
    )
    assert genes["B2M"]["ligandome_support"] == "input_only_coverage_unresolved"
    assert "ligandome" in completeness["missing_evidence"]
    assert "ligandome_gene_mapping" in completeness["missing_evidence"].split(";")


def test_assay_inputs_without_files_keep_not_assessed(tmp_path):
    _, completeness, summary, genes = _appm(tmp_path)
    assert genes["B2M"]["protein_status"] == "not_assessed"
    assert genes["B2M"]["ligandome_support"] == "not_assessed"
    assert completeness["functional_validation_status"] == "computational_proxy"
    assert summary["claim_strength"] == "computational_triage_only"
    assert "gene_mapping" not in completeness["missing_evidence"]


# --------------------------------------------------------------------------
# C5-05: paths never advertise files the explainability step did not write
# --------------------------------------------------------------------------

def test_submodule_scores_absent_when_explainability_fails(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("injected explainability failure")

    monkeypatch.setattr("neoag.appm_explainability.enhance_appm_outputs_v042", boom)
    paths = build_appm_2(sample_id="S", outdir=tmp_path / "appm")
    assert "appm_submodule_scores" not in paths
    assert not Path(paths["appm_gene_status"]).name.startswith("appm_submodule")
    assert (tmp_path / "appm" / "appm_explainability_warning.tsv").is_file()
    for key, value in paths.items():
        assert Path(value).is_file(), f"advertised path {key} does not exist"


def test_submodule_scores_advertised_only_after_successful_write(tmp_path):
    paths = build_appm_2(sample_id="S", outdir=tmp_path / "appm")
    assert "appm_submodule_scores" in paths
    assert Path(paths["appm_submodule_scores"]).is_file()
