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
    w_ba = float(w.get("netmhcpan_ba", 0.25))
    w_el = float(w.get("netmhcpan_el", 0.35))
    w_mhcf = float(w.get("mhcflurry_presentation", 0.30))
    w_proc = float(w.get("mhcflurry_processing", 0.10))
    w_chop = float(w.get("netchop_processing", 0.10))
    w_stab = float(w.get("netmhcstabpan_stability", 0.05))
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
    "peptide_id", "peptide", "hla_allele", "mhc_class", "base_score", "base_rank",
    "best_rank", "worst_rank", "rank_spread", "score_min", "score_max", "score_range",
    "verdict", "n_vectors",
    # Audit follow-up 2026-11-06 (ranking sensitivity modes): SMAA-style
    # acceptability columns; mode names baseline|global|regional.
    "mode", "expected_rank", "rank_std",
    "acceptability_top1", "acceptability_top10", "acceptability_top50", "acceptability_top100",
]

SENSITIVITY_MODES = ("baseline", "global", "regional")
SENSITIVITY_MAX_VECTORS = 4000


def _weight_vectors(mode: str, base: dict[str, float], perturb: float,
                    samples: int, concentration: float, seed: int) -> dict[str, dict[str, float]]:
    """Weight vectors to re-rank under, per sensitivity mode.

    baseline (audit A3): the base vector plus one up/down perturbation per
    active weight - 1 + 2W configurations, cheap enough for every run.
    global: the full +-perturb grid over all active weights (3^W), capped
    deterministically at SENSITIVITY_MAX_VECTORS.
    regional (SMAA-style): Dirichlet draws centred on the base weights - the
    audit's "Dirichlet distribution can be one of options in Regional mode" -
    concentrated around the scientifically motivated vector (higher
    concentration = tighter region), sampled with a fixed seed so a rerun of
    the same inputs produces the same table.
    """
    import random as _random

    keys = [k for k, v in base.items() if v > 0]
    if mode not in SENSITIVITY_MODES:
        raise ValueError(f"unknown sensitivity mode {mode!r}; expected one of {SENSITIVITY_MODES}")
    vectors: dict[str, dict[str, float]] = {"base": dict(base)}
    if mode == "baseline":
        for key, val in base.items():
            if val <= 0:
                continue
            for factor, tag in ((1.0 + perturb, "up"), (1.0 - perturb, "down")):
                vec = dict(base)
                vec[key] = val * factor
                vectors[f"{key}_{tag}"] = vec
        return vectors
    if mode == "global":
        combos = []
        for mask in range(3 ** len(keys)):
            rest, combo = mask, []
            for _ in keys:
                combo.append(rest % 3)
                rest //= 3
            combos.append(combo)
        if len(combos) > SENSITIVITY_MAX_VECTORS:
            stride = -(-len(combos) // SENSITIVITY_MAX_VECTORS)
            combos = combos[::stride]
        for index, combo in enumerate(combos):
            if all(value == 1 for value in combo):
                continue  # this is the base vector
            vec = dict(base)
            for key, state in zip(keys, combo):
                vec[key] = base[key] * (1.0 - perturb if state == 0 else 1.0 + perturb)
            vectors[f"grid{index}"] = vec
        return vectors
    rng = _random.Random(seed)
    alpha = [max(base[key] * concentration, 1e-6) for key in keys]
    total = sum(alpha)
    for index in range(samples):
        draws = [rng.gammavariate(a, 1.0) for a in alpha]
        drawn = sum(draws) or 1.0
        # Keep the Dirichlet draw at the base weight total so the perturbed
        # vector stays comparable (the weighted mean renormalizes anyway).
        scale = total / drawn
        vectors[f"dirichlet{index}"] = {
            **{key: base[key] for key in base if key not in keys},
            **{key: value * scale for key, value in zip(keys, draws)},
        }
    return vectors


def weight_sensitivity_table(evidence_rows, profile, perturb: float = 0.25,
                             mode: str = "baseline", samples: int = 300,
                             concentration: float = 10.0, seed: int = 20261106):
    """Audit item A3 (+2026-11-06 follow-up): weight-choice sensitivity.

    No predictors are rerun: every component score is already in the
    presentation evidence table, so the weighted sum is simply recombined under
    a set of weight vectors and the candidates re-ranked.  The result is a
    per-patient QC table - a candidate whose rank moves only because of the
    weight choice is not a robust pick.  Default weights come from the profile;
    the mechanism is exactly the sensitivity rerun the audit suggested, and
    NetMHCstabpan/NetChop participate with their real (tunable) weights.

    mode selects the vector set (see _weight_vectors):
      baseline  1 + 2W configurations (the original audit-A3 test),
      global    the full +-perturb grid (far more combinations),
      regional  SMAA-style Dirichlet draws around the base weights.
    Every mode additionally reports the stochastic acceptability of each
    candidate - the fraction of weight vectors under which it lands in the
    top-K - and its expected rank / rank dispersion across vectors (SMAA).
    This is robustness, not accuracy: it never changes the production score,
    it tells the reviewer which ranks survive the weight choice.
    """

    w = profile.get("presentation_weights", {}) or {}
    base = {
        "netmhcpan_ba": float(w.get("netmhcpan_ba", 0.25)),
        "netmhcpan_el": float(w.get("netmhcpan_el", 0.35)),
        "mhcflurry_presentation": float(w.get("mhcflurry_presentation", 0.30)),
        "mhcflurry_processing": float(w.get("mhcflurry_processing", 0.10)),
        "netchop_processing": float(w.get("netchop_processing", 0.10)),
        "netmhcstabpan_stability": float(w.get("netmhcstabpan_stability", 0.05)),
    }
    vectors = _weight_vectors(mode, base, perturb, samples, concentration, seed)

    def score_row(row, vec):
        is_mhci = str(row.get("mhc_class") or "I").upper() not in {"II", "MHC-II", "CLASSII"}
        parts = []
        for col, key, style in (
            ("netmhcpan_ba_rank", "netmhcpan_ba", "rank"),
            ("netmhcpan_el_rank", "netmhcpan_el", "rank"),
            ("mhcflurry_presentation_score", "mhcflurry_presentation", "score"),
            ("mhcflurry_processing_score", "mhcflurry_processing", "score"),
            ("netchop_31d_cterm_score", "netchop_processing", "score"),
            ("netmhcstabpan_rank", "netmhcstabpan_stability", "rank"),
        ):
            if key in ("netchop_processing", "netmhcstabpan_stability") and not is_mhci:
                continue
            raw = str(row.get(col, "")).strip()
            if raw == "":
                continue
            val = norm_rank(raw) if style == "rank" else clamp(to_float(raw, -1.0))
            if val is None:
                continue
            parts.append((float(val), vec[key]))
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
    for i, row in enumerate(evidence_rows):
        base_score = scores["base"][i]
        if base_score is None:
            continue
        rs = [rank_map[name][i] for name in vectors]
        spread = max(rs) - min(rs)
        svals = [vals[i] for vals in scores.values() if vals[i] is not None]
        verdict = "rank_stable" if spread == 0 else ("rank_minor_shift" if spread <= minor_band else "rank_unstable")
        expected = sum(rs) / len(rs)
        variance = sum((r - expected) ** 2 for r in rs) / len(rs)
        acceptability = {
            f"acceptability_top{k}": f"{sum(1 for r in rs if r <= k) / len(rs):.4f}"
            for k in (1, 10, 50, 100)
        }
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
            "mode": mode,
            "expected_rank": f"{expected:.2f}",
            "rank_std": f"{variance ** 0.5:.2f}",
            **acceptability,
        })
    return out


def write_weight_sensitivity(evidence_rows, profile, path, perturb: float = 0.25,
                             mode: str = "baseline", samples: int = 300,
                             concentration: float = 10.0, seed: int = 20261106):
    rows = weight_sensitivity_table(evidence_rows, profile, perturb=perturb,
                                    mode=mode, samples=samples,
                                    concentration=concentration, seed=seed)
    write_tsv(path, rows, WEIGHT_SENSITIVITY_FIELDS)
    return path
