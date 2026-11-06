from __future__ import annotations
from pathlib import Path
from .utils import read_tsv, write_tsv, to_float, norm_rank, clamp, safe_id
from .schemas import PRESENTATION_FIELDS
from .evidence_provenance import ProvenanceRecord, ProvenanceRegistry, provenance_derived, attach_provenance
from .adapters.peptide_input import normalize_hla_allele

def by_key(rows):
    result = {}
    for row in rows:
        peptide = str(row.get("peptide") or "").strip().upper()
        hla = normalize_hla_allele(str(row.get("hla_allele") or ""))
        # A caller-local peptide_hla_key may refer to an older peptide ID.
        # Exact sequence + normalized HLA is the canonical biological key.
        key = safe_id(f"{peptide}_{hla}") if peptide and hla else str(row.get("peptide_hla_key") or "")
        if key:
            result[key] = row
    return result


def _validate_core_predictor_coverage(peptides, evidence, tool):
    """Reject a supplied core predictor table that cannot map to this run."""
    if not evidence or not peptides:
        return
    candidate_keys = {
        safe_id(
            f"{str(row.get('peptide') or '').strip().upper()}_"
            f"{normalize_hla_allele(str(row.get('hla_allele') or ''))}"
        )
        for row in peptides
        if str(row.get("peptide") or "").strip()
        and normalize_hla_allele(str(row.get("hla_allele") or ""))
    }
    overlap = len(candidate_keys.intersection(evidence))
    minimum = 1 if len(candidate_keys) < 100 else max(5, int(len(candidate_keys) * 0.01))
    if overlap < minimum:
        raise ValueError(
            f"{tool} evidence maps to only {overlap}/{len(candidate_keys)} candidate peptide-HLA keys; "
            "check that the file is the normalized evidence table for this run rather than an "
            "incompatible raw output, stale copied result, or table with mismatched HLA notation"
        )

def grade(binding, presentation, complete):
    if complete <= 0:
        return "MISSING"
    if binding >= 0.85 and presentation >= 0.85:
        return "A"
    if binding >= 0.75 and presentation >= 0.65:
        return "B"
    if binding >= 0.65:
        return "C_BINDING_ONLY"
    return "D_WEAK"

def build_presentation_evidence(
    raw_peptides,
    netmhcpan,
    mhcflurry,
    profile,
    out=None,
    netmhcstabpan=None,
    netchop=None,
    provenance_registry: ProvenanceRegistry | None = None,
):
    peptides = read_tsv(raw_peptides)
    net = by_key(read_tsv(netmhcpan)) if netmhcpan else {}
    mhc = by_key(read_tsv(mhcflurry)) if mhcflurry else {}
    stab = by_key(read_tsv(netmhcstabpan)) if netmhcstabpan else {}
    chop_rows = read_tsv(netchop) if netchop else []
    chop_by_id = {row.get("peptide_id", ""): row for row in chop_rows if row.get("peptide_id")}
    chop_by_peptide = {row.get("peptide", ""): row for row in chop_rows if row.get("peptide")}
    _validate_core_predictor_coverage(peptides, net, "NetMHCpan")
    _validate_core_predictor_coverage(peptides, mhc, "MHCflurry")
    w = profile.get("presentation_weights", {})
    # Defaults come from PRESENTATION_WEIGHT_DEFAULTS - the same table the
    # weight-sensitivity module uses, so scoring and sensitivity can never
    # disagree about which components exist or what their defaults are.
    def _w(key: str) -> float:
        return float(w.get(key, PRESENTATION_WEIGHT_DEFAULTS[key]))
    w_ba = _w("netmhcpan_ba")
    w_el = _w("netmhcpan_el")
    w_mhcf = _w("mhcflurry_presentation")
    w_proc = _w("mhcflurry_processing")
    w_chop = _w("netchop_processing")
    w_stab = _w("netmhcstabpan_stability")
    # Audit follow-up 2026-11-06 (binding aggregation): the binding evidence
    # score historically took max(normalized NetMHCpan BA, normalized
    # MHCflurry affinity), which trusts whichever predictor is most optimistic
    # (audit "APPM integration" note; presentation-scoring investigation 6.2).
    # The aggregation is now a profile knob; the default stays "max" so every
    # existing run reproduces byte-identically, and mean/min are available for
    # A/B without code changes. |BA - MHCflurry| is exported either way as a
    # disagreement quality feature.
    binding_agg = str(w.get("binding_aggregation", "max")).strip().lower()
    if binding_agg not in {"max", "mean", "min"}:
        raise ValueError(
            f"presentation_weights.binding_aggregation must be max|mean|min, got {binding_agg!r}"
        )
    rows = []
    for p in peptides:
        hla = normalize_hla_allele(str(p.get("hla_allele") or ""))
        key = safe_id(f"{str(p.get('peptide','')).strip().upper()}_{hla}")
        wt_peptide = p.get("wildtype_peptide", "")
        wt_key = safe_id(f"{str(wt_peptide).strip().upper()}_{hla}") if wt_peptide else ""
        n = net.get(key, {})
        m = mhc.get(key, {})
        wt_n = net.get(wt_key, {}) if wt_key else {}
        wt_m = mhc.get(wt_key, {}) if wt_key else {}
        s = stab.get(key, {})
        wt_s = stab.get(wt_key, {}) if wt_key else {}
        c = chop_by_id.get(p.get("peptide_id", ""), {}) or chop_by_peptide.get(p.get("peptide", ""), {})
        ba = n.get("netmhcpan_ba_rank", "")
        el = n.get("netmhcpan_el_rank", "")
        pct = m.get("mhcflurry_affinity_percentile", "")
        proc = m.get("mhcflurry_processing_score", "")
        pres = m.get("mhcflurry_presentation_score", "")
        def evidence_value(key, default=""):
            return n.get(key) or m.get(key) or p.get(key) or default
        ba_s = norm_rank(ba) if ba != "" else None
        el_s = norm_rank(el) if el != "" else None
        pct_s = norm_rank(pct) if pct != "" else None
        proc_s = clamp(to_float(proc, -1)) if proc != "" else None
        pres_s = clamp(to_float(pres, -1)) if pres != "" else None
        chop_raw = c.get("netchop_31d_cterm_score") or c.get("netchop_31d_max_score", "")
        chop_s = clamp(to_float(chop_raw, -1)) if chop_raw != "" else None
        binding_parts = [x for x in [ba_s, pct_s] if x is not None]
        if binding_parts:
            if binding_agg == "mean":
                binding = sum(binding_parts) / len(binding_parts)
            elif binding_agg == "min":
                binding = min(binding_parts)
            else:
                binding = max(binding_parts)
        else:
            binding = norm_rank(p.get("binding_rank", 99))
        # Quality feature requested by the scoring investigation (6.2): how far
        # apart the two binding predictors disagree, in normalized-score units.
        binding_disagreement = (
            f"{abs(ba_s - pct_s):.4f}" if ba_s is not None and pct_s is not None else ""
        )
        num = den = 0.0
        is_mhci = str(p.get("mhc_class") or "I").upper() not in {"II", "MHC-II", "CLASSII"}
        weighted_evidence = [(ba_s,w_ba),(el_s,w_el),(pres_s,w_mhcf),(proc_s,w_proc)]
        if is_mhci:
            weighted_evidence.append((chop_s, w_chop))
            # NetMHCstabpan rank (percentile, lower = more stable) was read, matched
            # and written through this table but never reached the score. It enters as
            # a small signal: stability is conditional on binding, so it corroborates
            # rather than drives; weight is profile-tunable (audit items A2/A3).
            stab_rank = s.get("netmhcstabpan_rank", "")
            stab_s = norm_rank(stab_rank) if stab_rank != "" else None
            weighted_evidence.append((stab_s, w_stab))
        for val, wt in weighted_evidence:
            if val is not None:
                num += val * wt; den += wt
        if den:
            presentation = num / den
            complete = min(1.0, den / sum(wt for _, wt in weighted_evidence))
        else:
            presentation = clamp(to_float(p.get("presentation_score"), 0.0))
            complete = 0.25 if p.get("presentation_score") else 0.0
        rows.append({
            "peptide_id": p.get("peptide_id",""),
            "event_id": p.get("event_id",""),
            "sample_id": p.get("sample_id",""),
            "peptide": p.get("peptide",""),
            # Retained in memory so downstream immunogenicity merges can look up
            # the matched WT peptide. The canonical presentation schema stores
            # the resulting WT scores rather than this helper field.
            "wildtype_peptide": wt_peptide,
            "hla_allele": p.get("hla_allele",""),
            "mhc_class": p.get("mhc_class",""),
            "netmhcpan_ba_rank": str(to_float(ba, 99.0)),
            "netmhcpan_el_rank": str(to_float(el, 99.0)),
            "netmhcpan_mt_rank_ba": evidence_value("netmhcpan_mt_rank_ba", ba),
            "netmhcpan_mt_rank_el": evidence_value("netmhcpan_mt_rank_el", el),
            "netmhcpan_wt_rank_ba": wt_n.get("netmhcpan_ba_rank") or evidence_value("netmhcpan_wt_rank_ba"),
            "netmhcpan_wt_rank_el": wt_n.get("netmhcpan_el_rank") or evidence_value("netmhcpan_wt_rank_el"),
            "netmhcstabpan_score": str(to_float(s.get("netmhcstabpan_score"), 0.0)) if s else "",
            "netmhcstabpan_rank": str(to_float(s.get("netmhcstabpan_rank"), 99.0)) if s else "",
            "netmhcstabpan_wt_score": str(to_float(wt_s.get("netmhcstabpan_score"), 0.0)) if wt_s else "",
            "netmhcstabpan_wt_rank": str(to_float(wt_s.get("netmhcstabpan_rank"), 99.0)) if wt_s else "",
            "netchop_31d_max_score": c.get("netchop_31d_max_score", ""),
            "netchop_31d_mean_score": c.get("netchop_31d_mean_score", ""),
            "netchop_31d_cterm_score": c.get("netchop_31d_cterm_score", ""),
            "netchop_31d_cleavage_sites": c.get("netchop_31d_cleavage_sites", ""),
            "netchop_processing_status": c.get("netchop_processing_status", "ASSESSED" if c else "UNASSESSED"),
            "mhcflurry_affinity_percentile": str(to_float(pct, 99.0)),
            "mhcflurry_processing_score": str(to_float(proc, 0.0)),
            "mhcflurry_presentation_score": str(to_float(pres, 0.0)),
            "mhcflurry_wt_affinity_percentile": wt_m.get("mhcflurry_affinity_percentile") or evidence_value("mhcflurry_wt_affinity_percentile"),
            "mhcflurry_wt_processing_score": wt_m.get("mhcflurry_processing_score") or evidence_value("mhcflurry_wt_processing_score"),
            "mhcflurry_wt_presentation_score": wt_m.get("mhcflurry_presentation_score") or evidence_value("mhcflurry_wt_presentation_score"),
            "prime_wt_score": evidence_value("prime_wt_score"),
            "prime_wt_rank": evidence_value("prime_wt_rank"),
            "bigmhc_im_wt_score": evidence_value("bigmhc_im_wt_score"),
            "binding_evidence_score": f"{binding:.4f}",
            "binding_predictor_disagreement": binding_disagreement,
            "presentation_evidence_score": f"{presentation:.4f}",
            "evidence_completeness": f"{complete:.4f}",
            "presentation_evidence_grade": grade(binding, presentation, complete),
        })
    registry = provenance_registry or ProvenanceRegistry()
    if netmhcpan and not registry.has("netmhcpan"):
        registry.register_passthrough("netmhcpan", netmhcpan)
    if mhcflurry and not registry.has("mhcflurry"):
        registry.register_passthrough("mhcflurry", mhcflurry)
    if netmhcstabpan and not registry.has("netmhcstabpan"):
        registry.register_passthrough("netmhcstabpan", netmhcstabpan)
    if netchop and not registry.has("netchop"):
        registry.register_passthrough("netchop", netchop)
    if out:
        summary = registry.tool_summary_fields()
        composite = provenance_derived(
            "presentation_composite",
            out,
            upstream="netmhcpan+mhcflurry+netmhcstabpan+netchop+immunogenicity",
        )
        for row in rows:
            row.update(summary)
            row.update(composite.as_fields())
        write_tsv(out, rows, PRESENTATION_FIELDS)
    return rows


WEIGHT_SENSITIVITY_FIELDS = [
    "peptide_id", "peptide", "hla_allele", "mhc_class",
    "base_score", "base_rank", "best_rank", "worst_rank", "rank_spread",
    "score_min", "score_max", "score_range", "verdict",
    "n_vectors", "n_components", "perturb", "top_k",
    "expected_rank", "rank_std",
    "dominance", "effective_components", "dominance_factor",
    "consistency", "adjusted_score", "adjusted_rank", "top_selection",
]

# Single source of truth for the class-I presentation aggregate. The scoring
# builder and the sensitivity re-rank below both enumerate THIS spec and these
# defaults, so "which components are used and how many" is one fact in code -
# never a number restated (or mis-stated) in a report. The default profile
# activates six components; a profile that zeroes any of them changes the
# vector count, the component count and the tables automatically.
# (evidence column, profile weight key, value style, class-I-only)
PRESENTATION_COMPONENTS = (
    ("netmhcpan_ba_rank", "netmhcpan_ba", "rank", False),
    ("netmhcpan_el_rank", "netmhcpan_el", "rank", False),
    ("mhcflurry_presentation_score", "mhcflurry_presentation", "score", False),
    ("mhcflurry_processing_score", "mhcflurry_processing", "score", False),
    ("netchop_31d_cterm_score", "netchop_processing", "score", True),
    ("netmhcstabpan_rank", "netmhcstabpan_stability", "rank", True),
)

PRESENTATION_WEIGHT_DEFAULTS = {
    "netmhcpan_ba": 0.25,
    "netmhcpan_el": 0.35,
    "mhcflurry_presentation": 0.30,
    "mhcflurry_processing": 0.10,
    "netchop_processing": 0.10,
    "netmhcstabpan_stability": 0.05,
}

# Operator round of 2026-11-06: the three modes (baseline/global/regional)
# were options, not a bundle - keep ONE analysis, the local +-perturbation
# test, with a wider default perturbation. global/regional are removed.
SENSITIVITY_DEFAULT_PERTURB = 0.5
SENSITIVITY_DEFAULT_TOP_K = 10

SENSITIVITY_FORMULA = (
    "Components and weights: v_i = norm_rank(rank) = max(0, 1 - min(r,10)/10) for rank-style "
    "inputs or v_i = clamp(score) for score-style inputs; score = sum(w_i*v_i) / sum(w_i) over the "
    "components PRESENT in the row (class-I-only components skipped for class II; weight 0 = off). "
    "Perturbation vectors: the base vector plus, for each active component i (base weight > 0), the "
    "base vector with w_i scaled by (1+P) and by (1-P) -> 1 + 2n vectors, n = number of active "
    "components, P = perturb (default 0.5). "
    "Dominance: contribution share p_i = (w_i*v_i) / sum(w_j*v_j) under the BASE weights; "
    "dominance D = sum(p_i^2) (equals 1/k when k components share the score equally, approaches 1 "
    "when a single component drives); effective_components = 1/D; "
    "dominance_factor = min(1, effective_components/2) - full credit needs the score carried by at "
    "least two roughly equal drivers, a single-driver candidate keeps half credit. "
    "Selection: consistency = fraction of the 1+2n vectors under which the candidate stays in the "
    "top K (K = top_k, default 10); adjusted_score = base_score * consistency * dominance_factor; "
    "top_selection = yes only when the candidate is in the default top K AND consistency = 1. "
    "Stability under this test is robustness to the weight choice, NOT accuracy."
)


def base_sensitivity_weights(profile) -> dict[str, float]:
    w = (profile or {}).get("presentation_weights", {}) or {}
    return {key: float(w.get(key, PRESENTATION_WEIGHT_DEFAULTS[key]))
            for key in PRESENTATION_WEIGHT_DEFAULTS}


def _sensitivity_row_parts(row, vec):
    """(value, weight) for every component this row actually contributes."""
    is_mhci = str(row.get("mhc_class") or "I").upper() not in {"II", "MHC-II", "CLASSII"}
    parts = []
    for col, key, style, class_i_only in PRESENTATION_COMPONENTS:
        if class_i_only and not is_mhci:
            continue
        wt = float(vec.get(key, 0.0))
        if wt <= 0:
            continue
        raw = str(row.get(col, "")).strip()
        if raw == "":
            continue
        val = norm_rank(raw) if style == "rank" else clamp(to_float(raw, -1.0))
        if val is None:
            continue
        parts.append((float(val), wt))
    return parts


def weight_sensitivity_table(evidence_rows, profile,
                             perturb: float = SENSITIVITY_DEFAULT_PERTURB,
                             top_k: int = SENSITIVITY_DEFAULT_TOP_K):
    """Audit A3, rebuilt 2026-11-06 per operator: local weight sensitivity only.

    No predictors are rerun: every component value is already in the
    presentation evidence table, so the weighted sum is recombined under the
    base weight vector and, for each ACTIVE component (weight > 0 - the
    component set and count come from PRESENTATION_COMPONENTS and the profile,
    never a hardcoded number), under that component scaled by (1+-perturb).
    The default perturbation is 0.5 (weights stretched to half and one-and-a-
    half times), the wider range the operator asked for.

    Per candidate: the raw picture (best/worst rank, spread, verdict), an
    expected-rank/standard-deviation summary, and the selection pair the
    operator specified - consistency (share of vectors keeping the candidate
    in the top K) and a dominance penalty (a candidate whose score is driven
    by one component is punished: dominance_factor = min(1, (1/sum p_i^2)/2)).
    adjusted_score = base_score * consistency * dominance_factor orders the
    post-sensitivity table; top_selection=yes marks the consistent
    top-selection candidates (default top K AND never falling out of it).
    SENSITIVITY_FORMULA states every equation verbatim; the technical report
    renders both tables and that formula. This is robustness, not accuracy,
    and it never touches the production score.
    """
    perturb = float(perturb)
    if not (0.0 < perturb < 1.0):
        raise ValueError(f"perturb must be in (0, 1), got {perturb!r}")
    top_k = int(top_k)
    if top_k < 1:
        raise ValueError(f"top_k must be >= 1, got {top_k!r}")

    base = base_sensitivity_weights(profile)
    active = [key for key, val in base.items() if val > 0]
    vectors: dict[str, dict[str, float]] = {"base": dict(base)}
    for key in active:
        for factor, tag in ((1.0 + perturb, "up"), (1.0 - perturb, "down")):
            vectors[f"{key}_{tag}"] = {**base, key: base[key] * factor}

    def score_row(row, vec):
        parts = _sensitivity_row_parts(row, vec)
        den = sum(wt for _, wt in parts)
        if not den:
            return None
        return sum(val * wt for val, wt in parts) / den

    scores = {name: [score_row(r, vec) for r in evidence_rows] for name, vec in vectors.items()}

    def dense_ranks(vals):
        ordered = sorted({v for v in vals if v is not None}, reverse=True)
        pos = {v: i + 1 for i, v in enumerate(ordered)}
        return [pos[v] if v is not None else len(ordered) + 1 for v in vals]

    rank_map = {name: dense_ranks(vals) for name, vals in scores.items()}
    n = len(evidence_rows)
    minor_band = max(3, int(0.05 * n))

    out = []
    adjusted_values: list[float | None] = []
    for i, row in enumerate(evidence_rows):
        base_score = scores["base"][i]
        if base_score is None:
            continue
        rs = [rank_map[name][i] for name in vectors]
        spread = max(rs) - min(rs)
        svals = [vals[i] for vals in scores.values() if vals[i] is not None]
        verdict = "rank_stable" if spread == 0 else (
            "rank_minor_shift" if spread <= minor_band else "rank_unstable")
        expected = sum(rs) / len(rs)
        variance = sum((r - expected) ** 2 for r in rs) / len(rs)
        consistency = sum(1 for r in rs if r <= top_k) / len(rs)

        # Dominance under the BASE weights: how concentrated the score's
        # contributions are. Simpson concentration D = sum p_i^2; the penalty
        # punishes one or a few weights driving the score.
        contributions = [val * wt for val, wt in _sensitivity_row_parts(row, base)]
        total = sum(contributions)
        if total > 0:
            dominance = sum((c / total) ** 2 for c in contributions)
            effective = 1.0 / dominance if dominance > 0 else 0.0
        else:
            dominance = 1.0
            effective = 0.0
        dominance_factor = min(1.0, effective / 2.0)

        adjusted = base_score * consistency * dominance_factor
        adjusted_values.append(adjusted)
        out.append({
            "peptide_id": row.get("peptide_id", ""),
            "peptide": row.get("peptide", ""),
            "hla_allele": row.get("hla_allele", ""),
            "mhc_class": row.get("mhc_class", ""),
            "base_score": f"{base_score:.4f}",
            "base_rank": str(rank_map["base"][i]),
            "best_rank": str(min(rs)),
            "worst_rank": str(max(rs)),
            "rank_spread": str(spread),
            "score_min": f"{min(svals):.4f}",
            "score_max": f"{max(svals):.4f}",
            "score_range": f"{max(svals) - min(svals):.4f}",
            "verdict": verdict,
            "n_vectors": str(len(vectors)),
            "n_components": str(len(active)),
            "perturb": f"{perturb:.2f}",
            "top_k": str(top_k),
            "expected_rank": f"{expected:.2f}",
            "rank_std": f"{variance ** 0.5:.2f}",
            "dominance": f"{dominance:.4f}",
            "effective_components": f"{effective:.2f}",
            "dominance_factor": f"{dominance_factor:.4f}",
            "consistency": f"{consistency:.4f}",
            "adjusted_score": f"{adjusted:.4f}",
            "adjusted_rank": "",  # filled after the loop, once every value exists
            "top_selection": "yes" if (rank_map["base"][i] <= top_k and consistency >= 1.0 - 1e-12) else "no",
            "_adjusted_raw": adjusted,
        })

    adjusted_ranks = dense_ranks(adjusted_values)
    for item, adj_rank in zip(out, adjusted_ranks):
        item["adjusted_rank"] = str(adj_rank)
        item.pop("_adjusted_raw", None)
    return out


def write_weight_sensitivity(evidence_rows, profile, path,
                             perturb: float = SENSITIVITY_DEFAULT_PERTURB,
                             top_k: int = SENSITIVITY_DEFAULT_TOP_K):
    rows = weight_sensitivity_table(evidence_rows, profile,
                                    perturb=perturb, top_k=top_k)
    write_tsv(path, rows, WEIGHT_SENSITIVITY_FIELDS)
    return path
