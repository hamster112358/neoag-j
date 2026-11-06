"""Two separate ranking reports for one case (operator request, 2026-11-06).

The pipeline produces two independent rankings of the same candidates:

- the consensus-based LOGIC ranking (ranked_peptides.evidence_consensus.tsv):
  discrete evidence states and grades combined by the consensus rules
  (evidence_rank / evidence_grade), no numeric weights involved;
- the WEIGHT-based ranking (ranked_peptides.tsv): efficacy_score, a weighted
  mean over the Layer-3 evidence chains times the decision multipliers.

Both currently share one technical report. The operator asked for two separate
reports so each ranking can be read, reviewed and circulated on its own
terms. This module writes exactly that pair from an existing case directory;
it changes no ranking, no pipeline step and no existing report.

Design rules kept from the weight-sensitivity work: every quantity shown is
derived from the files/code at report time (chain columns come from
sensitivity_chains(), column lists are filtered against the actual headers),
and anything a case did not produce is stated as missing instead of faked.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .reports_dual import REPORT_CSS, _table, esc, weight_sensitivity_html
from .sensitivity import _dense_ranks, sensitivity_chains
from .utils import read_tsv, to_float

CONSENSUS_REPORT_NAME = "ranking_report_consensus.html"
WEIGHT_REPORT_NAME = "ranking_report_weight.html"

# Looked up under the case root, in order; a bare scoring dir also works.
_SCORING_SUBDIRS = ("final/scoring", "scoring", "")

# Columns curated for each report; entries a given file lacks are dropped at
# report time, so these are preference lists, not contracts.
_CONSENSUS_MAIN = [
    "evidence_rank", "track_rank", "peptide_id", "peptide", "gene", "hla_allele",
    "evidence_grade", "evidence_grade_uncapped", "evidence_consensus_score",
    "consensus_action", "source_chain_confidence_tier", "legacy_weighted_rank",
]
_CONSENSUS_LAYERS = [
    "event_authenticity_grade", "rna_support_grade", "presentation_consensus_grade",
    "mutant_specificity_grade", "safety_grade", "hla_appm_grade", "clonality_grade",
    "evidence_completeness_grade",
]
_CONSENSUS_EVENTS = [
    "event_id", "gene", "event_type", "evidence_track", "best_peptide_id",
    "best_peptide", "best_hla_allele", "best_peptide_evidence_rank",
    "best_evidence_grade", "source_chain_confidence_tier",
    "source_chain_confidence_label",
]
_WEIGHT_MAIN = [
    "rank", "peptide_id", "peptide", "gene", "hla_allele", "mhc_class",
    "efficacy_score", "immunology_composite_score",
]
_WEIGHT_MULTIPLIERS = [
    "safety_multiplier", "presentation_gate_multiplier", "appm_multiplier",
    "ccf_multiplier", "escape_multiplier", "cross_platform_multiplier",
    "mutant_specificity_multiplier",
]
_WEIGHT_TAIL = ["safety_status", "priority"]


def _page(title: str, body: list[str]) -> str:
    return ("<!DOCTYPE html><html><head><meta charset='utf-8'>"
            f"<title>{esc(title)}</title>{REPORT_CSS}</head><body>"
            + "\n".join(body) + "</body></html>")


def _present(rows: list[Mapping[str, Any]], columns: list[str]) -> list[str]:
    if not rows:
        return []
    fields = set(rows[0])
    return [c for c in columns if c in fields]


def _counts(rows: list[Mapping[str, Any]], column: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for row in rows:
        key = str(row.get(column) or "(blank)")
        out[key] = out.get(key, 0) + 1
    return out


def _counts_table(counts: dict[str, int], label: str) -> str:
    return _table([{"value": k, "rows": str(v)} for k, v in sorted(counts.items())],
                  [label, "rows"])


def find_scoring_dir(case_root: str | Path) -> Path:
    """Locate the scoring artifact dir under a case root (final/, scoring/, or bare)."""
    root = Path(case_root)
    for rel in _SCORING_SUBDIRS:
        cand = (root / rel) if rel else root
        if (cand / "ranked_peptides.tsv").exists() or \
                (cand / "ranked_peptides.evidence_consensus.tsv").exists():
            return cand
    raise FileNotFoundError(
        f"no ranked_peptides.tsv / ranked_peptides.evidence_consensus.tsv found under {root} "
        f"(looked in: {', '.join(d or '.' for d in _SCORING_SUBDIRS)})")


def _read_optional(path: Path) -> list[dict[str, str]]:
    return read_tsv(path) if path.exists() else []


def _consensus_meta(scoring_dir: Path) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    path = scoring_dir / "evidence_consensus_run.json"
    if path.exists():
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}
    return payload if isinstance(payload, dict) else {}


def build_consensus_report(scoring_dir: Path, top: int, case_label: str = "") -> str:
    rows = _read_optional(scoring_dir / "ranked_peptides.evidence_consensus.tsv")
    events = _read_optional(scoring_dir / "ranked_events.evidence_consensus.tsv")
    run = _consensus_meta(scoring_dir)
    label = f" - {esc(case_label)}" if case_label else ""
    body = [f"<h1>Consensus-Based Logic Ranking{label}</h1>"]
    body.append(
        "<p>This report presents ONE of the pipeline's two independent rankings: the "
        "<b>consensus-based logic ranking</b>. Candidates are ordered by <code>evidence_rank</code>, "
        "produced by discrete evidence states and grades combined through the consensus rules - "
        "no numeric weights are involved. The weight-based ranking (efficacy_score) is presented "
        f"separately in <code>{WEIGHT_REPORT_NAME}</code>, in the same directory.</p>")
    rules = {
        "algorithm": run.get("algorithm", "(run manifest not found)"),
        "rules_name": run.get("rules_name", ""),
        "rules_version": run.get("rules_version", ""),
        "status": run.get("status", ""),
    }
    body.append(_table([rules], ["algorithm", "rules_name", "rules_version", "status"]))
    if str(rules["status"]).upper() == "PROVISIONAL_RESEARCH_ONLY":
        body.append("<div class='warn'>The consensus rules are "
                    "PROVISIONAL_RESEARCH_ONLY: this ranking is an independent research view; "
                    "it does not modify the weighted baseline.</div>")
    if not rows:
        body.append("<div class='warn'>This case did not produce "
                    "<code>ranked_peptides.evidence_consensus.tsv</code> (the evidence-consensus "
                    "step did not run or its output was not shipped), so no consensus ranking "
                    "exists to report.</div></body></html>")
        return _page("Consensus-Based Logic Ranking", body)

    ordered = sorted(rows, key=lambda r: (to_float(r.get("evidence_rank"), 1e9),
                                          str(r.get("peptide_id") or "")))
    grades = _counts(ordered, "evidence_grade")
    actions = _counts(ordered, "consensus_action") if "consensus_action" in ordered[0] else {}
    tiers = _counts(ordered, "source_chain_confidence_tier") if "source_chain_confidence_tier" in ordered[0] else {}
    summary = [{"metric": "candidate peptide-HLA rows", "value": str(len(ordered))}]
    for label, counts in (("evidence_grade", grades), ("consensus_action", actions),
                          ("source_chain tier", tiers)):
        if counts:
            summary.append({"metric": label,
                            "value": "; ".join(f"{k}={v}" for k, v in sorted(counts.items()))})
    body.append("<div class='section'><h2>Summary</h2>")
    body.append(_table(summary, ["metric", "value"]))
    body.append("</div>")

    # Rank-shift context from the weighted-vs-consensus comparison, when present.
    comparison = _read_optional(scoring_dir / "weighted_vs_consensus_comparison.tsv")
    shift_by_id = {str(r.get("peptide_id") or ""): r for r in comparison}
    enriched = []
    for i, row in enumerate(ordered[:top], start=1):
        item = dict(row)
        shift = shift_by_id.get(str(row.get("peptide_id") or ""), {})
        for key in ("rank_shift_weighted_minus_consensus", "rank_shift_direction"):
            if key in (comparison[0] if comparison else {}):
                item[key] = shift.get(key, "")
        enriched.append(item)
    headers = _present(enriched, _CONSENSUS_MAIN) + _present(enriched, _CONSENSUS_LAYERS)
    extra = [c for c in ("rank_shift_weighted_minus_consensus", "rank_shift_direction")
             if comparison and c in comparison[0]]
    headers += extra
    body.append(f"<div class='section'><h2>Consensus ranking - top {min(top, len(ordered))} "
                f"of {len(ordered)} candidates</h2>")
    body.append("<p>Layer grade columns are the per-question verdicts the consensus rules "
                "combined; <code>legacy_weighted_rank</code> and the rank-shift columns show "
                "where the same candidate sits in the weight-based ranking, for orientation "
                "only.</p>")
    body.append(_table(enriched, headers))
    body.append("</div>")

    if events:
        event_headers = _present(events, _CONSENSUS_EVENTS)
        body.append(f"<div class='section'><h2>Consensus event ranking "
                    f"({len(events)} event rows)</h2>")
        body.append(_table(events, event_headers))
        body.append("</div>")
    return _page("Consensus-Based Logic Ranking", body)


def build_weight_report(scoring_dir: Path, top: int, case_label: str = "") -> str:
    rows = _read_optional(scoring_dir / "ranked_peptides.tsv")
    label = f" - {esc(case_label)}" if case_label else ""
    body = [f"<h1>Weight-Based Ranking{label}</h1>"]
    body.append(
        "<p>This report presents ONE of the pipeline's two independent rankings: the "
        "<b>weight-based ranking</b>. Candidates are ordered by <code>efficacy_score</code> = "
        "weighted mean over the Layer-3 evidence chains (weights from the deployed profile "
        "defaults) times the decision multipliers. The consensus-based logic ranking is "
        f"presented separately in <code>{CONSENSUS_REPORT_NAME}</code>, in the same directory.</p>")
    if not rows:
        body.append("<div class='warn'>This case did not produce "
                    "<code>ranked_peptides.tsv</code>, so no weight-based ranking exists to "
                    "report.</div>")
        return _page("Weight-Based Ranking", body)

    ordered = sorted(rows, key=lambda r: (-(to_float(r.get("efficacy_score"), -1.0)),
                                          str(r.get("peptide_id") or "")))
    ranks = _dense_ranks([to_float(r.get("efficacy_score"), 0.0) for r in ordered])
    for row, rank in zip(ordered, ranks):
        row["rank"] = str(rank)

    chain_cols = [(name, col) for name, col, _w in sensitivity_chains(None)]
    active_chains = [(name, col) for name, col in chain_cols if col in ordered[0]]
    summary = [{"metric": "candidate rows ranked", "value": str(len(ordered))}]
    eff_counts = _counts(ordered, "efficacy_score")
    summary.append({"metric": "distinct efficacy_score values", "value": str(len(eff_counts))})
    if "safety_status" in ordered[0]:
        summary.append({"metric": "safety_status",
                        "value": "; ".join(f"{k}={v}" for k, v in sorted(_counts(ordered, "safety_status").items()))})
    body.append("<div class='section'><h2>Summary</h2>")
    body.append(_table(summary, ["metric", "value"]))
    body.append("</div>")

    headers = (_present(ordered, _WEIGHT_MAIN)
               + [col for _name, col in active_chains]
               + _present(ordered, _WEIGHT_MULTIPLIERS)
               + _present(ordered, _WEIGHT_TAIL))
    body.append(f"<div class='section'><h2>Weight ranking - top {min(top, len(ordered))} "
                f"of {len(ordered)} candidates</h2>")
    body.append(
        "<p>The l3_* columns are the Layer-3 evidence-chain inputs, listed because the score "
        "is a weighted mean over exactly those columns (chain set derived from the code and "
        "profile at report time, not hardcoded). The multiplier columns are the decision "
        "factors applied after the weighted mean (safety, presentation gate, APPM, CCF, "
        "escape, cross-platform, mutant specificity): they are non-compensatory verdicts, not "
        "weights - a safety rejection must not be tradable against a high binding score. When "
        "every candidate shares the same verdicts they are run-constant and only rescale the "
        "score; they reorder candidates as soon as any verdict differs. The weight-sensitivity "
        "section below prints each candidate's multiplier product (chain_multiplier) so "
        "constancy is answerable from the table.</p>")
    body.append(_table(ordered[:top], headers))
    body.append("</div>")

    sens = _read_optional(scoring_dir / "weight_sensitivity.tsv")
    chains = _read_optional(scoring_dir / "weight_sensitivity_chains.tsv")
    if sens:
        body.append(weight_sensitivity_html(sens, chains))
    else:
        body.append("<div class='info'>No <code>weight_sensitivity.tsv</code> in this case: "
                    "run the <code>weight-sensitivity</code> command against this case's "
                    "ranked_peptides.tsv to add the regional-Dirichlet stability section "
                    "to this report.</div>")
    return _page("Weight-Based Ranking", body)


def write_ranking_reports(case_dir: str | Path, out_dir: str | Path | None = None,
                          top: int = 50) -> dict[str, Path]:
    """Write the two ranking reports for one case; returns {consensus, weight} paths."""
    scoring_dir = find_scoring_dir(case_dir)
    out = Path(out_dir) if out_dir else Path(case_dir) / "reports"
    out.mkdir(parents=True, exist_ok=True)
    consensus_path = out / CONSENSUS_REPORT_NAME
    weight_path = out / WEIGHT_REPORT_NAME
    case_label = Path(str(case_dir).rstrip("/")).name
    consensus_path.write_text(
        build_consensus_report(scoring_dir, top, case_label), encoding="utf-8")
    weight_path.write_text(
        build_weight_report(scoring_dir, top, case_label), encoding="utf-8")
    return {"consensus": consensus_path, "weight": weight_path}
