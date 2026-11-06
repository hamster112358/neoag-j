"""Composite immunogenicity from PRIME, BigMHC_IM, and optional IEDB fallback."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .utils import clamp, norm_rank, to_float
from .evidence_provenance import (
    ProvenanceRegistry,
    STATUS_INVALID_FOR_PRODUCTION,
    STATUS_REAL,
    STATUS_STUB,
)

IMMUNOGENICITY_FIELDS = [
    "prime_score",
    "prime_rank",
    "bigmhc_im_score",
    "deepimmuno_score",
    "iedb_immunogenicity_score",
    "immunogenicity_composite_score",
    "immunogenicity_source",
]


def immunogenicity_config(profile: Mapping[str, Any]) -> dict[str, Any]:
    cfg = profile.get("immunogenicity", {})
    weights = cfg.get("weights") or {}
    if not isinstance(weights, dict):
        weights = {}
    hla_mult = cfg.get("hla_source_multipliers") or {}
    if not isinstance(hla_mult, dict):
        hla_mult = {}
    sources = cfg.get("sources") or ["prime", "bigmhc_im"]
    return {
        "enabled": bool(cfg.get("enabled", True)),
        "sources": list(sources),
        "composite": str(cfg.get("composite", "mean")).lower(),
        "weights": {str(k): float(v) for k, v in weights.items()},
        "hla_source_multipliers": hla_mult,
        "use_iedb_fallback": bool(cfg.get("use_iedb_fallback", False)),
    }


def _canonical_hla(hla_allele: str) -> str:
    from .adapters.peptide_input import normalize_hla_allele

    return normalize_hla_allele(hla_allele) if hla_allele.strip() else ""


def _merge_key(peptide: Any, hla: Any) -> tuple[str, str]:
    """Canonical (peptide, HLA) join key for evidence merges (review C2-03).

    The adapter side already indexes evidence rows by upper-cased peptide and
    normalized HLA (adapters/prime._pair_key and friends); the merge side must
    apply the identical normalization or notation drift such as A*02:01 versus
    HLA-A*02:01 silently produces empty scores.  Original strings stay in the
    rows themselves for traceability.
    """
    return str(peptide or "").strip().upper(), _canonical_hla(str(hla or ""))


def effective_source_weights(
    hla_allele: str,
    parts: Mapping[str, float],
    profile: Mapping[str, Any],
) -> dict[str, float]:
    """Base source weights with optional per-HLA multipliers (0 = exclude source)."""
    cfg = immunogenicity_config(profile)
    base = cfg["weights"]
    multipliers_table = cfg["hla_source_multipliers"]
    canon = _canonical_hla(hla_allele)
    hla_mult: Mapping[str, Any] = {}
    if canon and isinstance(multipliers_table.get(canon), Mapping):
        hla_mult = multipliers_table[canon]

    out: dict[str, float] = {}
    for name in parts:
        w = float(base.get(name, 1.0))
        mult = float(hla_mult.get(name, 1.0)) if hla_mult else 1.0
        out[name] = w * mult
    return out


def normalize_iedb_score(raw: Any) -> float:
    return clamp((to_float(raw, 0.0) + 0.5) / 1.5)


def normalize_prime_component(row: Mapping[str, Any]) -> float | None:
    # Review C2-01: only %Rank is on a scale comparable with the other
    # components.  Raw PRIME scores are log-odds; clamping them saturated every
    # positive value at 1.0, and no uncalibrated transform may replace that
    # (review disposition: require empirical calibration or leave the
    # component unassessed).  A row without a rank therefore yields no prime
    # component instead of a fabricated one.
    rank_raw = str(row.get("prime_rank", "")).strip()
    if rank_raw:
        return norm_rank(rank_raw)
    return None


def _is_neutral_or_blank(raw: str) -> bool:
    """Blank or numeric-neutral (0.5) score string (review C2-06).

    The old sentinel set {"", "0.5", "0.5000"} was string-compared, so the
    heuristic's "%.5f" output ("0.50000") slipped through as a genuine call.
    Interpretation must not depend on decimal formatting.  Unparseable values
    count as unassessed too: converting them through to_float(..., 0.0) would
    fabricate a component score from garbage.
    """
    if not raw:
        return True
    try:
        return abs(float(raw) - 0.5) < 1e-9
    except ValueError:
        return True


def normalize_bigmhc_im_component(row: Mapping[str, Any]) -> float | None:
    raw = str(row.get("bigmhc_im_score", "")).strip()
    if not raw:
        return None
    return clamp(to_float(raw, 0.0))


def normalize_deepimmuno_component(row: Mapping[str, Any]) -> float | None:
    raw = str(row.get("deepimmuno_score", "")).strip()
    if not raw:
        return None
    return clamp(to_float(raw, 0.0))


def component_scores(row: Mapping[str, Any], profile: Mapping[str, Any]) -> dict[str, float]:
    if _is_stub_evidence(row):
        return {}
    cfg = immunogenicity_config(profile)
    if not cfg["enabled"]:
        return {}

    parts: dict[str, float] = {}
    for src in cfg["sources"]:
        if src == "prime":
            val = normalize_prime_component(row)
            if val is not None:
                parts["prime"] = val
        elif src == "bigmhc_im":
            val = normalize_bigmhc_im_component(row)
            if val is not None:
                parts["bigmhc_im"] = val
        elif src == "deepimmuno":
            val = normalize_deepimmuno_component(row)
            if val is not None:
                parts["deepimmuno"] = val
        elif src == "iedb":
            raw = str(row.get("iedb_immunogenicity_score", "")).strip()
            if not _is_neutral_or_blank(raw):
                parts["iedb"] = normalize_iedb_score(raw)

    if not parts and cfg["use_iedb_fallback"]:
        raw = str(row.get("iedb_immunogenicity_score", "")).strip()
        if not _is_neutral_or_blank(raw):
            parts["iedb"] = normalize_iedb_score(raw)
    return parts


def _is_stub_evidence(row: Mapping[str, Any]) -> bool:
    # Review C2-02: the row-level status wins first — a real tool file can
    # still contain placeholder rows the tool never returned (C1-01), which
    # the file-level column cannot see.  The file-level fallback also treats
    # invalid_for_production as unusable, not just plain stub.
    row_status = str(row.get("evidence_row_status", "")).strip().lower()
    if row_status:
        return row_status != STATUS_REAL
    return str(row.get("evidence_status", "")).strip().lower() in {STATUS_STUB, STATUS_INVALID_FOR_PRODUCTION}


def _all_stub_source_rows(rows: list[dict[str, str]]) -> bool:
    return bool(rows) and all(
        str(row.get("source_file", "")).strip().lower() == "stub" for row in rows
    )


def _non_stub_evidence_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if not _is_stub_evidence(row)]


def _immunogenicity_tool_available(name: str, ctx: Any) -> bool:
    if name == "bigmhc_im":
        import os
        from pathlib import Path

        bigmhc_dir = Path(ctx.executables.get("bigmhc_dir") or os.environ.get("BIGMHC_DIR", ""))
        return bool(str(bigmhc_dir)) and (bigmhc_dir / "src" / "predict.py").is_file()
    from .tools.runner import check_tool

    try:
        return check_tool(name, ctx.exe(name)).available
    except Exception:
        return False


def combine_component_scores(
    parts: Mapping[str, float],
    profile: Mapping[str, Any],
    *,
    hla_allele: str = "",
) -> float:
    if not parts:
        return 0.0
    cfg = immunogenicity_config(profile)
    mode = cfg["composite"]
    weights = effective_source_weights(hla_allele, parts, profile)
    active = {name: val for name, val in parts.items() if weights.get(name, 0.0) > 0.0}
    if not active:
        return 0.0

    if mode == "max":
        return clamp(max(active.values()))

    if mode == "weighted":
        num = 0.0
        den = 0.0
        for name, val in active.items():
            w = weights[name]
            num += w * val
            den += w
        return clamp(num / den) if den > 0 else 0.0

    # Review C2-06: "mean" is a strict arithmetic mean over the active
    # components.  It previously switched to weighting whenever HLA multipliers
    # existed, violating the mode-name contract; weighted combination now only
    # happens under its own named mode above.  Multipliers still gate sources
    # in or out (a zero multiplier excludes the source via `active`).
    return clamp(sum(active.values()) / len(active))


def resolve_immunogenicity_score(
    peptide: Mapping[str, Any],
    presentation: Mapping[str, Any] | None,
    profile: Mapping[str, Any],
) -> tuple[float, dict[str, str]]:
    """Return 0–1 immunogenicity score and metadata fields for output rows."""
    row = {**peptide, **(presentation or {})}
    parts = component_scores(row, profile)
    hla = str(row.get("hla_allele", "")).strip()
    composite = combine_component_scores(parts, profile, hla_allele=hla)
    meta = {
        "immunogenicity_composite_score": f"{composite:.4f}" if parts else "",
        "immunogenicity_source": "+".join(sorted(parts)) if parts else "none",
    }
    return composite, meta


def has_resolved_immunogenicity(
    peptide: Mapping[str, Any],
    presentation: Mapping[str, Any] | None,
    profile: Mapping[str, Any],
) -> bool:
    row = {**peptide, **(presentation or {})}
    if component_scores(row, profile):
        return True
    raw = str(peptide.get("immunogenicity_score", "")).strip()
    if not _is_neutral_or_blank(raw):
        val = to_float(raw, 0.5)
        return abs(val - 0.5) >= 1e-6 and val > 0.0
    return False


def enrich_row_immunogenicity(row: dict[str, str], profile: Mapping[str, Any]) -> None:
    composite, meta = resolve_immunogenicity_score(row, row, profile)
    row.update(meta)


def merge_prime(
    presentation_rows: list[dict[str, str]],
    prime_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = _merge_key(row.get("peptide", ""), row.get("hla_allele", ""))
        prime = prime_by_pair.get(key, {})
        row["prime_score"] = prime.get("prime_score", "")
        row["prime_rank"] = prime.get("prime_rank", "")
        wt_key = _merge_key(row.get("wildtype_peptide", ""), row.get("hla_allele", ""))
        wt_prime = prime_by_pair.get(wt_key, {}) if wt_key[0] else {}
        row["prime_wt_score"] = wt_prime.get("prime_score", "")
        row["prime_wt_rank"] = wt_prime.get("prime_rank", "")


def merge_bigmhc_im(
    presentation_rows: list[dict[str, str]],
    bigmhc_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = _merge_key(row.get("peptide", ""), row.get("hla_allele", ""))
        hit = bigmhc_by_pair.get(key, {})
        row["bigmhc_im_score"] = hit.get("bigmhc_im_score", "")
        wt_key = _merge_key(row.get("wildtype_peptide", ""), row.get("hla_allele", ""))
        wt_hit = bigmhc_by_pair.get(wt_key, {}) if wt_key[0] else {}
        row["bigmhc_im_wt_score"] = wt_hit.get("bigmhc_im_score", "")


def merge_deepimmuno(
    presentation_rows: list[dict[str, str]],
    deepimmuno_by_pair: dict[tuple[str, str], dict[str, str]],
) -> None:
    for row in presentation_rows:
        key = _merge_key(row.get("peptide", ""), row.get("hla_allele", ""))
        hit = deepimmuno_by_pair.get(key, {})
        row["deepimmuno_score"] = hit.get("deepimmuno_score", "")


def apply_immunogenicity_composite(
    presentation_rows: list[dict[str, str]],
    profile: Mapping[str, Any],
) -> None:
    for row in presentation_rows:
        enrich_row_immunogenicity(row, profile)


def _evidence_fingerprint(raw_peptides: Path, name: str) -> dict[str, str]:
    """What must not change for a stored *_evidence.tsv to stay reusable (C2-04).

    Mirrors the C1-06 runner-sidecar contract: the peptide table bytes cover
    the peptide set and HLA notation, the producer version covers the software
    side.  A file-exists check alone silently applied stale predictions to
    changed peptides or an upgraded predictor.
    """
    from .evidence_provenance import detect_tool_version

    digest = hashlib.sha256()
    pep = Path(raw_peptides)
    if pep.is_file():
        digest.update(pep.read_bytes())
    return {
        "tool": name,
        "peptide_input_sha256": digest.hexdigest(),
        "tool_version": str(detect_tool_version(name) or ""),
    }


def _evidence_fingerprint_path(evidence_path: Path) -> Path:
    return evidence_path.with_name(evidence_path.name + ".inputs.json")


def _evidence_fingerprint_matches(evidence_path: Path, fingerprint: dict[str, str]) -> bool:
    """Stored sidecar must agree on every field; absent sidecar = no reuse."""
    try:
        stored = json.loads(_evidence_fingerprint_path(evidence_path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if not isinstance(stored, dict):
        return False
    return all(str(stored.get(key, "")) == value for key, value in fingerprint.items())


def _write_evidence_fingerprint(evidence_path: Path, fingerprint: dict[str, str]) -> None:
    _evidence_fingerprint_path(evidence_path).write_text(
        json.dumps(fingerprint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def run_immunogenicity_predictors(
    raw_peptides: Path,
    outdir: Path,
    profile: Mapping[str, Any],
    ctx: Any,
    *,
    skip: bool = False,
    skip_sources: set[str] | None = None,
    reuse_sources: set[str] | None = None,
    provenance_registry: ProvenanceRegistry | None = None,
) -> dict[str, Path | None]:
    """Run configured immunogenicity tools and write presentation evidence TSVs."""
    from .adapters.bigmhc_im import parse_bigmhc_im, write_bigmhc_im_evidence
    from .adapters.deepimmuno import parse_deepimmuno, write_deepimmuno_evidence
    from .adapters.iedb_immunogenicity import predict_pairs as predict_iedb_pairs, write_immunogenicity_evidence
    from .adapters.prime import parse_prime, write_prime_evidence
    from .tools.runner import run_bigmhc_im, run_deepimmuno, run_prime

    cfg = immunogenicity_config(profile)
    skip_sources = skip_sources or set()
    reuse_sources = reuse_sources or set()
    registry = provenance_registry or ProvenanceRegistry()
    tools_dir = outdir / "tools"
    pres_dir = outdir / "presentation"
    pres_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path | None] = {
        "prime": None,
        "bigmhc_im": None,
        "deepimmuno": None,
        "iedb": None,
    }

    prime_tool = tools_dir / "prime.tsv"
    prime_reuse = pres_dir / "prime_evidence.tsv"
    prime_fp = _evidence_fingerprint(raw_peptides, "prime")
    if "prime" in cfg["sources"] and "prime" not in skip_sources:
        if "prime" in reuse_sources and prime_reuse.is_file() and prime_reuse.stat().st_size > 0:
            if _evidence_fingerprint_matches(prime_reuse, prime_fp):
                paths["prime"] = prime_reuse
                registry.register_passthrough("prime", prime_reuse)
            else:
                print(
                    f"[neoag] prime: reuse declined for {prime_reuse} (no matching input "
                    "fingerprint); recomputing instead of trusting stale evidence (review C2-04)",
                    flush=True,
                )
        prime_can_use_tool_output = skip or ctx.stub or _immunogenicity_tool_available("prime", ctx)
        if paths["prime"] is not None:
            pass
        elif not skip and not prime_can_use_tool_output:
            registry.register_missing("prime")
        elif not skip:
            run_prime(ctx, prime_tool)
        if (
            paths["prime"] is None and prime_can_use_tool_output and prime_tool.is_file()
            and _evidence_fingerprint_matches(prime_tool, prime_fp)
        ):
            prime_path = pres_dir / "prime_evidence.tsv"
            prime_rows = parse_prime(prime_tool, ctx.sample_id)
            # A non-strict fallback may have stubbed the whole tool output;
            # register that as stub, not real (review C1-01).
            prime_is_stub = ctx.stub or _all_stub_source_rows(prime_rows)
            write_prime_evidence(
                prime_path,
                prime_rows,
                provenance=registry.register_stub("prime") if prime_is_stub else registry.register_real("prime", prime_tool),
            )
            _write_evidence_fingerprint(prime_path, prime_fp)
            paths["prime"] = prime_path
        elif prime_reuse.is_file() and _evidence_fingerprint_matches(prime_reuse, prime_fp):
            paths["prime"] = prime_reuse
            registry.register_passthrough("prime", paths["prime"])
        elif prime_reuse.is_file():
            print(
                f"[neoag] prime: existing {prime_reuse} does not match this run's peptide input; "
                "not reusing stale evidence (review C2-04)",
                flush=True,
            )
            registry.register_missing("prime")
        elif prime_can_use_tool_output and prime_tool.is_file():
            print(
                f"[neoag] prime: tool output {prime_tool} lacks a matching input fingerprint; "
                "refusing to re-parse stale tool output (review C2-04)",
                flush=True,
            )
            registry.register_missing("prime")
        elif ctx.stub:
            registry.register_stub("prime")
        else:
            registry.register_missing("prime")
    else:
        registry.register_not_used("prime")

    bigmhc_tool = tools_dir / "bigmhc_im.tsv"
    bigmhc_reuse = pres_dir / "bigmhc_im_evidence.tsv"
    bigmhc_fp = _evidence_fingerprint(raw_peptides, "bigmhc_im")
    if "bigmhc_im" in cfg["sources"] and "bigmhc_im" not in skip_sources:
        if "bigmhc_im" in reuse_sources and bigmhc_reuse.is_file() and bigmhc_reuse.stat().st_size > 0:
            if _evidence_fingerprint_matches(bigmhc_reuse, bigmhc_fp):
                paths["bigmhc_im"] = bigmhc_reuse
                registry.register_passthrough("bigmhc_im", bigmhc_reuse)
            else:
                print(
                    f"[neoag] bigmhc_im: reuse declined for {bigmhc_reuse} (no matching input "
                    "fingerprint); recomputing instead of trusting stale evidence (review C2-04)",
                    flush=True,
                )
        bigmhc_can_use_tool_output = skip or ctx.stub or _immunogenicity_tool_available("bigmhc_im", ctx)
        if paths["bigmhc_im"] is not None:
            pass
        elif not skip and not bigmhc_can_use_tool_output:
            registry.register_missing("bigmhc_im")
        elif not skip:
            run_bigmhc_im(ctx, bigmhc_tool)
        if (
            paths["bigmhc_im"] is None and bigmhc_can_use_tool_output and bigmhc_tool.is_file()
            and _evidence_fingerprint_matches(bigmhc_tool, bigmhc_fp)
        ):
            bigmhc_path = pres_dir / "bigmhc_im_evidence.tsv"
            bigmhc_rows = parse_bigmhc_im(bigmhc_tool, ctx.sample_id)
            bigmhc_is_stub = ctx.stub or _all_stub_source_rows(bigmhc_rows)
            write_bigmhc_im_evidence(
                bigmhc_path,
                bigmhc_rows,
                provenance=registry.register_stub("bigmhc_im") if bigmhc_is_stub else registry.register_real("bigmhc_im", bigmhc_tool),
            )
            _write_evidence_fingerprint(bigmhc_path, bigmhc_fp)
            paths["bigmhc_im"] = bigmhc_path
        elif bigmhc_reuse.is_file() and _evidence_fingerprint_matches(bigmhc_reuse, bigmhc_fp):
            paths["bigmhc_im"] = bigmhc_reuse
            registry.register_passthrough("bigmhc_im", paths["bigmhc_im"])
        elif bigmhc_reuse.is_file():
            print(
                f"[neoag] bigmhc_im: existing {bigmhc_reuse} does not match this run's peptide input; "
                "not reusing stale evidence (review C2-04)",
                flush=True,
            )
            registry.register_missing("bigmhc_im")
        elif bigmhc_can_use_tool_output and bigmhc_tool.is_file():
            print(
                f"[neoag] bigmhc_im: tool output {bigmhc_tool} lacks a matching input fingerprint; "
                "refusing to re-parse stale tool output (review C2-04)",
                flush=True,
            )
            registry.register_missing("bigmhc_im")
        elif ctx.stub:
            registry.register_stub("bigmhc_im")
        else:
            registry.register_missing("bigmhc_im")
    else:
        registry.register_not_used("bigmhc_im")

    deepimmuno_tool = tools_dir / "deepimmuno.tsv"
    deepimmuno_reuse = pres_dir / "deepimmuno_evidence.tsv"
    deepimmuno_fp = _evidence_fingerprint(raw_peptides, "deepimmuno")
    if "deepimmuno" in cfg["sources"] and "deepimmuno" not in skip_sources:
        if "deepimmuno" in reuse_sources and deepimmuno_reuse.is_file() and deepimmuno_reuse.stat().st_size > 0:
            if _evidence_fingerprint_matches(deepimmuno_reuse, deepimmuno_fp):
                paths["deepimmuno"] = deepimmuno_reuse
                registry.register_passthrough("deepimmuno", deepimmuno_reuse)
            else:
                print(
                    f"[neoag] deepimmuno: reuse declined for {deepimmuno_reuse} (no matching input "
                    "fingerprint); recomputing instead of trusting stale evidence (review C2-04)",
                    flush=True,
                )
        deepimmuno_can_use_tool_output = skip or ctx.stub or _immunogenicity_tool_available("deepimmuno", ctx)
        if paths["deepimmuno"] is not None:
            pass
        elif not skip and not deepimmuno_can_use_tool_output:
            registry.register_missing("deepimmuno")
        elif not skip:
            run_deepimmuno(ctx, deepimmuno_tool)
        if (
            paths["deepimmuno"] is None and deepimmuno_can_use_tool_output and deepimmuno_tool.is_file()
            and _evidence_fingerprint_matches(deepimmuno_tool, deepimmuno_fp)
        ):
            deepimmuno_path = pres_dir / "deepimmuno_evidence.tsv"
            deepimmuno_rows = parse_deepimmuno(deepimmuno_tool, ctx.sample_id)
            deepimmuno_is_stub = ctx.stub or _all_stub_source_rows(deepimmuno_rows)
            write_deepimmuno_evidence(
                deepimmuno_path,
                deepimmuno_rows,
                provenance=registry.register_stub("deepimmuno") if deepimmuno_is_stub else registry.register_real("deepimmuno", deepimmuno_tool),
            )
            _write_evidence_fingerprint(deepimmuno_path, deepimmuno_fp)
            paths["deepimmuno"] = deepimmuno_path
        elif deepimmuno_reuse.is_file() and _evidence_fingerprint_matches(deepimmuno_reuse, deepimmuno_fp):
            paths["deepimmuno"] = deepimmuno_reuse
            registry.register_passthrough("deepimmuno", paths["deepimmuno"])
        elif deepimmuno_reuse.is_file():
            print(
                f"[neoag] deepimmuno: existing {deepimmuno_reuse} does not match this run's peptide input; "
                "not reusing stale evidence (review C2-04)",
                flush=True,
            )
            registry.register_missing("deepimmuno")
        elif deepimmuno_can_use_tool_output and deepimmuno_tool.is_file():
            print(
                f"[neoag] deepimmuno: tool output {deepimmuno_tool} lacks a matching input fingerprint; "
                "refusing to re-parse stale tool output (review C2-04)",
                flush=True,
            )
            registry.register_missing("deepimmuno")
        elif ctx.stub:
            registry.register_stub("deepimmuno")
        else:
            registry.register_missing("deepimmuno")
    else:
        registry.register_not_used("deepimmuno")

    iedb_path = pres_dir / "iedb_immunogenicity.tsv"
    if "iedb" in cfg["sources"] or cfg["use_iedb_fallback"]:
        iedb_fp = _evidence_fingerprint(raw_peptides, "iedb")
        iedb_rewritten = False
        if not skip and "iedb" in cfg["sources"]:
            from .tools.prep import unique_peptide_hla_pairs
            pairs = unique_peptide_hla_pairs(raw_peptides)
            write_immunogenicity_evidence(
                iedb_path,
                predict_iedb_pairs(pairs),
                provenance=registry.register_derived("iedb", iedb_path),
            )
            _write_evidence_fingerprint(iedb_path, iedb_fp)
            iedb_rewritten = True
        if iedb_path.is_file() and (iedb_rewritten or _evidence_fingerprint_matches(iedb_path, iedb_fp)):
            paths["iedb"] = iedb_path
            if not registry.has("iedb"):
                registry.register_derived("iedb", iedb_path)
        elif iedb_path.is_file():
            print(
                f"[neoag] iedb: existing {iedb_path} does not match this run's peptide input; "
                "not reusing stale derived evidence (review C2-04)",
                flush=True,
            )
        elif "iedb" in cfg["sources"]:
            registry.register_not_used("iedb")
    else:
        registry.register_not_used("iedb")
    return paths


def apply_immunogenicity_evidence(
    pres_rows: list[dict[str, str]],
    immuno_paths: dict[str, Path | None],
    profile: Mapping[str, Any],
) -> None:
    from .adapters.bigmhc_im import bigmhc_by_pair
    from .adapters.deepimmuno import deepimmuno_by_pair
    from .adapters.prime import prime_by_pair
    from .utils import read_tsv

    if immuno_paths.get("prime"):
        merge_prime(pres_rows, prime_by_pair(_non_stub_evidence_rows(read_tsv(immuno_paths["prime"]))))
    if immuno_paths.get("bigmhc_im"):
        merge_bigmhc_im(pres_rows, bigmhc_by_pair(_non_stub_evidence_rows(read_tsv(immuno_paths["bigmhc_im"]))))
    if immuno_paths.get("deepimmuno"):
        merge_deepimmuno(pres_rows, deepimmuno_by_pair(_non_stub_evidence_rows(read_tsv(immuno_paths["deepimmuno"]))))
    if immuno_paths.get("iedb"):
        iedb_by_pair = {
            _merge_key(r.get("peptide", ""), r.get("hla_allele", "")): r.get("iedb_immunogenicity_score", "")
            for r in read_tsv(immuno_paths["iedb"])
        }
        for row in pres_rows:
            key = _merge_key(row.get("peptide", ""), row.get("hla_allele", ""))
            row["iedb_immunogenicity_score"] = iedb_by_pair.get(key, "")
    apply_immunogenicity_composite(pres_rows, profile)
