from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Any

from ..utils import read_tsv, first, to_float
from .identity import normalize_contig


def load_expression(path: str | Path | None) -> dict[str, float]:
    if not path or not Path(path).exists():
        return {}
    out: dict[str, float] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "gene_name", "symbol", "SYMBOL"], "")
        if gene:
            out[gene] = to_float(first(r, ["TPM", "tpm", "expression", "expr", "gene_tpm"], "0"), 0.0)
    return out


def load_normal_expression(path: str | Path | None) -> dict[str, dict[str, Any]]:
    if not path or not Path(path).exists():
        return {}
    out: dict[str, dict[str, Any]] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "gene_name", "symbol", "SYMBOL"], "")
        if not gene:
            continue
        out[gene] = {
            "normal_tissue_max_tpm": to_float(first(r, ["normal_tissue_max_tpm", "max_tpm", "TPM"], "0"), 0.0),
            "normal_hspc_tpm": to_float(first(r, ["normal_hspc_tpm", "hspc_tpm"], "0"), 0.0),
            "critical_tissue_hit": first(r, ["critical_tissue_hit", "critical"], "no"),
        }
    return out


def normal_expr_for_genes(genes: list[str], normal_expr: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    nt = 0.0
    nh = 0.0
    crit = False
    for g in genes:
        e = normal_expr.get(g, {})
        nt = max(nt, to_float(e.get("normal_tissue_max_tpm"), 0.0))
        nh = max(nh, to_float(e.get("normal_hspc_tpm"), 0.0))
        crit = crit or str(e.get("critical_tissue_hit", "")).strip().lower() in {"1", "true", "yes", "y"}
    return {"normal_tissue_max_tpm": nt, "normal_hspc_tpm": nh, "critical_tissue_hit": "yes" if crit else "no"}


def load_normal_ligands(path: str | Path | None) -> set[str]:
    if not path or not Path(path).exists():
        return set()
    return {
        first(r, ["peptide", "Peptide", "sequence"], "").strip().upper()
        for r in read_tsv(path)
        if first(r, ["peptide", "Peptide", "sequence"], "").strip()
    }


def junction_positional_key(chrom1: str, pos1: Any, chrom2: str, pos2: Any) -> str:
    """Strand-insensitive sorted junction coordinate key (mirrors the exact
    positional index used by the SV pipeline)."""
    left = (normalize_contig(chrom1), int(float(pos1)))
    right = (normalize_contig(chrom2), int(float(pos2)))
    if right < left:
        left, right = right, left
    return f"{left[0]}:{left[1]}|{right[0]}:{right[1]}"


@dataclass
class JunctionReadEvidence:
    """Review B4-18: a gene-pair label alone must never transfer read support.

    Two different breakpoints can involve the same gene pair; joining on
    "G1::G2" silently attributed one junction's reads to the other event.
    Evidence transfer now happens only through ``exact`` keys (event-level
    identifiers or junction coordinates).  ``gene_pair`` is retained solely
    for exploratory grouping/reporting and is never consumed as support.
    """

    exact: dict[str, int] = field(default_factory=dict)
    gene_pair: dict[str, int] = field(default_factory=dict)


def load_junction_reads(path: str | Path | None) -> JunctionReadEvidence:
    """Load optional RNA junction support, split by join strictness.

    Exact keys: event_id, sv_event_id, or the junction's own coordinates
    (chrom1/pos1 + chrom2/pos2).  Gene-pair keys (gene_pair/fusion column or
    gene1/gene2) go to the exploratory bucket only.
    Value columns: junction_reads, split_reads, spanning_reads, reads.
    """
    out = JunctionReadEvidence()
    if not path or not Path(path).exists():
        return out
    for r in read_tsv(path):
        val = int(to_float(first(r, ["junction_reads", "split_reads", "spanning_reads", "reads"], "0"), 0.0))
        exact_keys = [first(r, ["event_id"], ""), first(r, ["sv_event_id"], "")]
        c1 = first(r, ["chrom1", "chr1"], "")
        p1 = first(r, ["pos1", "breakpoint1"], "")
        c2 = first(r, ["chrom2", "chr2"], "")
        p2 = first(r, ["pos2", "breakpoint2"], "")
        if c1 and c2 and p1 and p2:
            try:
                exact_keys.append(junction_positional_key(c1, p1, c2, p2))
            except (TypeError, ValueError):
                pass
        for k in exact_keys:
            if k:
                out.exact[k] = max(out.exact.get(k, 0), val)
        pair_keys = [first(r, ["gene_pair", "fusion"], "")]
        g1 = first(r, ["gene1", "left_gene"], "")
        g2 = first(r, ["gene2", "right_gene"], "")
        if g1 and g2:
            pair_keys.extend([f"{g1}::{g2}", f"{g2}::{g1}"])
        for k in pair_keys:
            if k:
                out.gene_pair[k] = max(out.gene_pair.get(k, 0), val)
    return out
