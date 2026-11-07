"""Pre-run safety-resource validation gate (review C6-02, round 3).

One explicit check of every safety reference BEFORE the expensive safety
assessment and ranking, so a missing or wrong resource surfaces as a named
preflight row instead of an unexplained concentration of capped/unassessed
candidates at the end of a run.

Each resource is validated on six axes, kept as DISTINCT statuses:

- MISSING            the file was never provided / does not exist
- UNREADABLE         it exists but cannot be opened
- PARSE_FAILURE      it opens but is not consumable (missing header/required
                     columns, no FASTA records, undecodable content)
- EMPTY_VALID        it is well-formed with a header but has zero data rows:
                     structurally fine, evidentially worthless
- BUILD_MISMATCH     it declares a genome build different from the run
- COVERAGE_LOW       it parses, but its locus coverage over the candidate
                     gene set is below the declared minimum

OK is the only passing status.  Two operating modes, chosen explicitly:

- exploratory (default): missing/deficient resources are TOLERATED, but every
  affected evidence layer must remain UNASSESSED downstream and the priority
  caps stay in force (the report states exactly which consequences apply).
- release: contract-required resources must be OK or covered by a documented
  waiver (resource, reason, authorized_by).  Waivers may excuse MISSING,
  EMPTY_VALID and COVERAGE_LOW only - an unreadable, unparseable or wrong-
  build resource is a correctness defect, never a waiverable gap.

The criteria are versioned (SAFETY_RELEASE_CRITERIA) so release behavior
cannot silently drift between runs, and the whole report is machine-readable
(JSON + TSV) for the run manifest and the release-integrity audit.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

SAFETY_PREFLIGHT_SCHEMA_VERSION = "safety_preflight@v1"
SAFETY_RELEASE_CRITERIA = "safety_release_criteria@v1"

# Statuses a documented waiver may excuse in release mode.  Anything else
# (UNREADABLE, PARSE_FAILURE, BUILD_MISMATCH) is a correctness defect.
WAIVERABLE_STATUSES = {"MISSING", "EMPTY_VALID", "COVERAGE_LOW"}
OK_STATUSES = {"OK", "NOT_APPLICABLE", "WAIVED"}

# What each resource must look like to parse, and what it is required for.
# ``scope`` states which candidate classes consume the resource: a resource
# for a candidate class that is absent from the run is NOT_APPLICABLE and can
# never produce a false blocking condition.
RESOURCE_SPECS: dict[str, dict[str, Any]] = {
    "reference_proteome": {
        "kind": "fasta",
        "scope": "peptides",
        "required": True,
        "evidence_consequence": (
            "reference_proteome_exact_match and novelty stay UNASSESSED; "
            "candidates keep the missing-resource priority cap"
        ),
    },
    "normal_hla_ligands": {
        "kind": "tsv",
        "scope": "peptides",
        "required": False,
        "required_for_release": True,
        "required_columns_any": [["peptide"], ["sequence"], ["peptide_seq"]],
        "evidence_consequence": "normal immunopeptidome match layer stays UNASSESSED",
    },
    "normal_expression": {
        "kind": "tsv",
        "scope": "all_events",
        "required": True,
        "required_columns_any": [["gene"], ["symbol"], ["gene_id"]],
        "locus_column_any": ["gene", "symbol", "gene_id"],
        "evidence_consequence": (
            "normal-tissue and HSPC context stay UNASSESSED; safety_evidence_"
            "incomplete partial status and its priority cap stay in force"
        ),
    },
    "normal_junctions": {
        "kind": "tsv",
        "scope": "junction_events",
        "required": False,
        "required_for_release": True,
        "evidence_consequence": "normal-junction safety layer stays UNASSESSED for junction candidates",
    },
}

_JUNCTION_EVENT_TOKENS = ("fusion", "splice", "junction")


def _declared_build(path: Path) -> str:
    sidecar = path.with_name(path.name + ".build")
    if sidecar.is_file():
        try:
            return sidecar.read_text(encoding="utf-8").strip()
        except OSError:
            return ""
    try:
        opener = __import__("gzip").open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                low = line.lower()
                if "genome_build" in low or "build" in low:
                    for token in ("GRCh38", "GRCh37", "hg38", "hg19", "T2T"):
                        if token.lower() in low:
                            return token
    except OSError:
        return ""
    return ""


def _has_junction_scope(event_types: Iterable[str]) -> bool:
    return any(tok in str(t).lower() for t in event_types for tok in _JUNCTION_EVENT_TOKENS)


def _scan_tsv(path: Path, required_columns_any: list[list[str]] | None, locus_column: str | None,
              wanted_genes: set[str]) -> tuple[str, int, str, float]:
    """Return (status, rows, parse_detail, gene_coverage_fraction)."""
    opener = __import__("gzip").open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            header = fh.readline()
            if not header:
                # Not even a header line: not a "valid empty" table.
                return "PARSE_FAILURE", 0, "file has no header line", 0.0
            header_cols = header.rstrip("\r\n").split("\t")
            if required_columns_any is not None:
                if not any(all(col in header_cols for col in option) for option in required_columns_any):
                    return "PARSE_FAILURE", 0, f"header lacks any required column group {required_columns_any}", 0.0
            locus_idx = header.index(locus_column) if locus_column and locus_column in header_cols else None
            rows = 0
            seen_genes: set[str] = set()
            for line in fh:
                if not line.strip():
                    continue
                rows += 1
                if locus_idx is not None and wanted_genes:
                    parts = line.rstrip("\r\n").split("\t")
                    if locus_idx < len(parts) and parts[locus_idx].strip():
                        seen_genes.add(parts[locus_idx].strip())
            coverage = 0.0
            if locus_idx is not None and wanted_genes:
                coverage = len(seen_genes & wanted_genes) / len(wanted_genes)
            if rows == 0:
                return "EMPTY_VALID", 0, "", coverage
            return "OK", rows, "", coverage
    except (OSError, UnicodeError):
        return "UNREADABLE", 0, "", 0.0


def _scan_fasta(path: Path) -> tuple[str, int, str]:
    opener = __import__("gzip").open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            saw_record_header = False
            saw_content = False
            records = 0
            for line in fh:
                if line.startswith(">"):
                    saw_record_header = True
                    records += 1
                elif line.strip():
                    saw_content = True
            if not saw_record_header:
                if saw_content:
                    return "PARSE_FAILURE", 0, "sequence data without FASTA headers"
                return "EMPTY_VALID", 0, ""
            return "OK", records, ""
    except (OSError, UnicodeError):
        return "UNREADABLE", 0, ""


def inspect_resource(
    name: str,
    path: str | Path | None,
    *,
    applicable: bool = True,
    candidate_genes: set[str] | None = None,
    genome_build: str = "",
    min_coverage: float = 0.5,
) -> dict[str, str]:
    spec = RESOURCE_SPECS[name]
    row: dict[str, str] = {
        "resource": name,
        "path": str(path or ""),
        "scope": str(spec["scope"]),
        "required": "yes" if spec.get("required") else "no",
        "required_for_release": "yes" if spec.get("required_for_release", spec.get("required")) else "no",
        "status": "",
        "rows": "",
        "sha256": "",
        "declared_build": "",
        "gene_coverage": "",
        "detail": "",
        "evidence_consequence": str(spec["evidence_consequence"]),
    }
    if not applicable:
        row["status"] = "NOT_APPLICABLE"
        row["detail"] = f"no {spec['scope']} candidates in this run; absence of this resource never blocks"
        return row
    if not path or not Path(path).is_file():
        row["status"] = "MISSING"
        return row
    p = Path(path)
    try:
        row["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
    except OSError:
        row["status"] = "UNREADABLE"
        return row
    row["declared_build"] = _declared_build(p)
    if genome_build and row["declared_build"] and row["declared_build"].lower() != str(genome_build).lower():
        row["status"] = "BUILD_MISMATCH"
        row["detail"] = f"resource declares {row['declared_build']}, run uses {genome_build}"
        return row
    if spec["kind"] == "fasta":
        status, rows, detail = _scan_fasta(p)
        row.update(status=status, rows=str(rows), detail=detail)
    else:
        locus_column = None
        for col in spec.get("locus_column_any", []):
            locus_column = col
            break
        status, rows, detail, coverage = _scan_tsv(
            p, spec.get("required_columns_any"), locus_column, candidate_genes or set(),
        )
        row.update(status=status, rows=str(rows), detail=detail)
        if locus_column and candidate_genes:
            row["gene_coverage"] = f"{coverage:.4f}"
            if status == "OK" and coverage < min_coverage:
                row["status"] = "COVERAGE_LOW"
                row["detail"] = f"covers {row['gene_coverage']} of {len(candidate_genes)} candidate genes (minimum {min_coverage:g})"
    if row["status"] == "OK" and not row["declared_build"] and genome_build:
        row["detail"] = (row["detail"] + "; " if row["detail"] else "") + "genome build undeclared by the resource (disclosed, not blocking)"
    return row


def run_safety_preflight(
    *,
    paths: Mapping[str, str | Path | None],
    event_types: Iterable[str] = (),
    candidate_genes: Iterable[str] = (),
    genome_build: str = "",
    mode: str = "exploratory",
    waivers: Iterable[Mapping[str, str]] = (),
    min_coverage: float = 0.5,
) -> dict[str, Any]:
    mode = str(mode or "exploratory").lower()
    if mode not in {"exploratory", "release"}:
        raise ValueError(f"unknown safety preflight mode {mode!r}; expected exploratory|release")
    event_types = [str(t) for t in event_types]
    genes = {str(g) for g in candidate_genes if str(g).strip()}
    have_peptides = bool(event_types)  # resource applicability is event-type aware
    junction_run = _has_junction_scope(event_types) if event_types else True
    applicability = {
        "reference_proteome": have_peptides,
        "normal_hla_ligands": have_peptides,
        "normal_expression": True,
        "normal_junctions": junction_run,
    }
    waiver_map = {
        str(w.get("resource", "")): w for w in waivers if str(w.get("resource", ""))
    }
    rows: list[dict[str, str]] = []
    blocking: list[str] = []
    unresolved: list[str] = []
    for name in RESOURCE_SPECS:
        row = inspect_resource(
            name, paths.get(name),
            applicable=applicability.get(name, True),
            candidate_genes=genes, genome_build=genome_build, min_coverage=min_coverage,
        )
        waiver = waiver_map.get(name)
        if row["status"] != "OK" and row["status"] != "NOT_APPLICABLE" and waiver:
            if row["status"] in WAIVERABLE_STATUSES and str(waiver.get("reason", "")).strip() and str(waiver.get("authorized_by", "")).strip():
                row["status"] = "WAIVED"
                row["detail"] = (row["detail"] + "; " if row["detail"] else "") + (
                    f"waived: {str(waiver.get('reason')).strip()} (authorized_by={str(waiver.get('authorized_by')).strip()})"
                )
            elif row["status"] not in WAIVERABLE_STATUSES:
                row["detail"] = (row["detail"] + "; " if row["detail"] else "") + (
                    f"waiver REJECTED for {row['status']} (correctness defects are never waiverable)"
                )
        if row["status"] not in OK_STATUSES:
            unresolved.append(name)
            required_in_release = row["required_for_release"] == "yes"
            if mode == "release" and required_in_release:
                blocking.append(name)
        rows.append(row)
    if mode == "release":
        decision = "BLOCK_RELEASE" if blocking else "PROCEED"
    else:
        decision = "PROCEED_WITH_UNASSESSED" if unresolved else "PROCEED"
    return {
        "schema_version": SAFETY_PREFLIGHT_SCHEMA_VERSION,
        "release_criteria": SAFETY_RELEASE_CRITERIA,
        "mode": mode,
        "genome_build": genome_build,
        "event_types_scopes": {
            "peptides": have_peptides, "junction_events": junction_run,
        },
        "release_decision": decision,
        "blocking_resources": ";".join(blocking),
        "unresolved_resources": ";".join(unresolved),
        "rows": rows,
    }


def write_safety_preflight(report: Mapping[str, Any], out_json: str | Path, out_tsv: str | Path | None = None) -> tuple[Path, ...]:
    from .utils import write_tsv

    out = Path(out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    written: tuple[Path, ...] = (out,)
    if out_tsv:
        tsv_path = Path(out_tsv)
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        cols = [
            "resource", "path", "scope", "required", "required_for_release", "status",
            "rows", "sha256", "declared_build", "gene_coverage", "detail", "evidence_consequence",
        ]
        write_tsv(tsv_path, list(report.get("rows", [])), cols)
        written = (out, tsv_path)
    return written


def summarize_preflight_for_manifest(report: Mapping[str, Any]) -> dict[str, str]:
    """Compact form for a run manifest / provenance record."""
    return {
        "safety_preflight_schema": str(report.get("schema_version", "")),
        "safety_preflight_mode": str(report.get("mode", "")),
        "safety_preflight_decision": str(report.get("release_decision", "")),
        "safety_preflight_unresolved": str(report.get("unresolved_resources", "")),
        "safety_preflight_blocking": str(report.get("blocking_resources", "")),
    }
