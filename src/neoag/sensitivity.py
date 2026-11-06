"""Final-ranking weight sensitivity (regional Dirichlet, operator design 2026-11-06).

ONE analysis mode: the regional test. The default weight vector is the seed;
every draw pulls the weights from a Dirichlet distribution centred on that
default (Gamma draws rescaled to the default total, fixed RNG seed). The
region width is the concentration parameter (lower = wider); the number of
tests is the sample count. "Wider possible perturbation" therefore means a
wider sampling region with more draws - not a bigger perturbation in any
single run.

The analysis runs on the FINAL peptide ranking, not on the presentation
sub-score: the evidence chains are the Layer-3 dimensions of
default_l3_weights(profile), derived at run time from the code and profile -
the chain count and names are never hardcoded (n_chains column). The stored
per-candidate modifier chain M = efficacy_score / immunology_composite_score
(safety, APPM, presentation gate, CCF, escape, cross-platform, specificity
multipliers - constant across draws) is recombined with every perturbed
composite, so at the default vector the recomputed score equals the stored
efficacy_score exactly.

Stability is not accuracy: the ranking is also penalised when a single chain
or a few chains drive the score (Simpson concentration dominance factor), and
candidates whose rank moves wildly across the region are penalised by their
rank inter-quartile range. The production score is never modified; this is a
QC view written by the separate ``weight-sensitivity`` command.
"""
from __future__ import annotations

import random
from pathlib import Path
from statistics import median, quantiles

from .utils import read_tsv, write_tsv, to_float, clamp
from .model_layers import default_l3_weights, l3_dimension_field

SENSITIVITY_DEFAULT_SAMPLES = 600
SENSITIVITY_DEFAULT_CONCENTRATION = 3.0
SENSITIVITY_DEFAULT_TOP_K = 50
SENSITIVITY_DEFAULT_SEED = 20261106

WEIGHT_SENSITIVITY_FIELDS = [
    "peptide_id", "peptide", "hla_allele", "mhc_class", "gene",
    "base_score", "base_rank", "expected_rank", "median_rank", "rank_iqr",
    "best_rank", "worst_rank", "rank_spread",
    "score_min", "score_max", "score_range", "expected_score",
    "consistency", "dominance", "effective_chains", "dominance_factor",
    "rank_stability", "verdict",
    "adjusted_score", "adjusted_rank", "top_selection",
    "n_draws", "n_chains", "n_scorable", "concentration", "top_k", "seed",
]

CHAIN_DIAG_FIELDS = ["chain", "default_weight", "mean_share", "n_candidates"]

# Operator correction of 2026-11-06: keep only the regional mode (Dirichlet
# draws seeded at the profile default). The earlier local +-perturbation mode
# is removed. The formula below is rendered verbatim in the technical report.
SENSITIVITY_FORMULA = (
    "What is re-ranked: the FINAL peptide ranking. Its evidence chains are the Layer-3 "
    "dimensions of default_l3_weights(profile); the chain names and the count are derived "
    "from the code and profile at run time and printed in the n_chains column (deployed "
    "default profile: ten chains), never hardcoded here. The six predictor components of "
    "the presentation table are not separate final-ranking chains: they are averaged into "
    "the single hla_presentation chain first, and that inner weighting stays at its "
    "default - this test perturbs only the top-level chain weights. "
    "Scores: s_k = stored l3_<chain>_score (missing or blank = 0, exactly as production "
    "scoring treats it); under a weight vector w: composite(w) = sum(w_k*s_k) / sum(w_k); "
    "M = efficacy_score / immunology_composite_score is the stored per-candidate chain of "
    "multipliers (safety, APPM, presentation gate, CCF, escape, cross-platform, mutant "
    "specificity - constant across draws); score(w) = composite(w) * M. When the stored "
    "composite is 0 every score is 0 regardless of M, so M is taken as 1 there. M is "
    "deliberately NOT sampled: "
    "its factors are discrete decisions (safety rejection, presentation-gate tier, CCF "
    "state, escape flag), not uncertain weights - perturbing a decision belongs in a "
    "separate what-if scenario (re-rank with the chosen decision flipped), not in the "
    "weight region. At the default "
    "vector score(w) equals the stored efficacy_score exactly. "
    "Regional draws: the default vector is the seed; each draw takes w_k = W * g_k / "
    "sum(g) with g_k ~ Gamma(base_k * concentration, 1) - a Dirichlet distribution centred "
    "on the defaults (W = sum of base weights; the draw mean of every w_k is its default, "
    "so the region is unbiased). concentration (default 3.0) sets the region width - lower "
    "draws are spread wider over the weight simplex; samples (default 600) is the number "
    "of draws, plus the default vector itself as an anchor, all reproducible from the "
    "fixed seed. Chains whose default weight is 0 stay 0 in every draw. "
    "Per candidate over all draws (dense ranks, ties share a rank): median_rank, "
    "expected_rank, rank_std omitted for brevity, rank_iqr = inclusive quartiles Q3-Q1, "
    "best/worst/spread, and consistency = fraction of draws with rank <= K "
    "(the SMAA rank-acceptability index at K = top_k, default 50). "
    "Dominance penalty (operator rule: one or a few weights driving the score is "
    "punished): p_k = base_k*s_k / sum(base_j*s_j) under the DEFAULT vector; dominance "
    "D = sum(p_k^2) (= 1/k when k chains share the score equally, -> 1 for a single "
    "driver); effective_chains = 1/D; dominance_factor = min(1, effective_chains/2) - "
    "full credit needs the score carried by at least two roughly equal chains. "
    "Stability ranking: rank_stability = 1 - rank_iqr / n_candidates; "
    "adjusted_score = expected_score * dominance_factor * rank_stability; adjusted_rank "
    "orders the second table. top_selection = yes only when base_rank <= K AND "
    "consistency = 1. Verdict: rank_stable (spread 0), rank_minor_shift (spread <= 5% of "
    "n_candidates, floor 3), rank_unstable otherwise. "
    "Stability under this test is robustness to the weight choice, NOT accuracy."
)


def sensitivity_chains(profile) -> list[tuple[str, str, float]]:
    """(chain_name, column, default_weight) for every final-ranking chain.

    Derived from default_l3_weights(profile) - the same function production
    scoring uses - so "which chains and how many" is one fact in code.
    """
    weights = default_l3_weights(profile or {})
    return [(name, l3_dimension_field(name), float(wt))
            for name, wt in weights.items()]


def dirichlet_weight_vectors(base: dict[str, float], samples: int,
                             concentration: float, seed: int) -> list[dict[str, float]]:
    """Dirichlet region around the default vector, deterministic from seed.

    Draw k: g_i = Gamma(base_i * concentration, 1) for the active (base > 0)
    chains, w_i = total * g_i / sum(g); zero-weight chains stay exactly 0.
    The first vector is the untouched default (the anchor).
    """
    total = sum(base.values())
    active = [key for key, wt in base.items() if wt > 0]
    vectors: list[dict[str, float]] = [dict(base)]
    rng = random.Random(seed)
    alphas = [(key, base[key] * concentration) for key in active]
    for _ in range(int(samples)):
        gammas = [(key, rng.gammavariate(alpha, 1.0)) for key, alpha in alphas]
        gsum = sum(g for _, g in gammas)
        vec = dict(base)
        if gsum > 0:
            for key, g in gammas:
                vec[key] = total * g / gsum
        vectors.append(vec)
    return vectors


def _dense_ranks(values: list[float]) -> list[int]:
    ordered = sorted(set(values), reverse=True)
    pos = {v: i + 1 for i, v in enumerate(ordered)}
    return [pos[v] for v in values]


def weight_sensitivity_table(peptide_rows, profile,
                             samples: int = SENSITIVITY_DEFAULT_SAMPLES,
                             concentration: float = SENSITIVITY_DEFAULT_CONCENTRATION,
                             top_k: int = SENSITIVITY_DEFAULT_TOP_K,
                             seed: int = SENSITIVITY_DEFAULT_SEED):
    """Regional Dirichlet stability of the final ranking; returns (rows, chain_diag).

    rows covers every candidate in the default top-K or the median top-K
    (K = top_k, default 50, so the 13th, 21st, ... all get a systematic
    verdict - the window is a REPORTING window, not a selection cap);
    chain_diag gives each chain's default weight and its mean contribution
    share over that window. A candidate with no stored efficacy_score or
    immunology_composite_score cannot be re-scored and is excluded from the
    rank pool (counted in n_scorable).
    """
    samples = int(samples)
    concentration = float(concentration)
    top_k = int(top_k)
    seed = int(seed)
    if samples < 1:
        raise ValueError(f"samples must be >= 1, got {samples!r}")
    if concentration <= 0:
        raise ValueError(f"concentration must be > 0, got {concentration!r}")
    if top_k < 1:
        raise ValueError(f"top_k must be >= 1, got {top_k!r}")

    chains = sensitivity_chains(profile)
    base = {name: wt for name, _, wt in chains}
    total_w = sum(base.values())
    if not chains or total_w <= 0:
        return [], []

    scorable = []
    for row in peptide_rows:
        eff = to_float(str(row.get("efficacy_score") or "").strip(), None)
        comp = to_float(str(row.get("immunology_composite_score") or "").strip(), None)
        if eff is None or comp is None:
            continue
        values = {name: clamp(to_float(row.get(col), 0.0)) for name, col, _ in chains}
        # composite == 0 implies every chain value is 0 (all weights >= 0 and
        # the composite is their weighted mean), so score(w) = 0 for every w
        # and the multiplier chain value is irrelevant; 1.0 keeps it finite.
        multiplier = eff / comp if comp > 0 else 1.0
        scorable.append((row, values, multiplier))

    n = len(scorable)
    if n == 0:
        return [], []

    vectors = dirichlet_weight_vectors(base, samples, concentration, seed)
    names = list(base)
    score_cols = []
    for vec in vectors:
        score_cols.append([
            sum(vec[name] * values[name] for name in names) / total_w * mult
            for _, values, mult in scorable
        ])
    rank_cols = [_dense_ranks(col) for col in score_cols]

    minor_band = max(3, int(0.05 * n))
    out = []
    for i, (row, values, _mult) in enumerate(scorable):
        rs = [col[i] for col in rank_cols]
        ss = [col[i] for col in score_cols]
        expected_rank = sum(rs) / len(rs)
        sorted_rs = sorted(rs)
        q1, _q2, q3 = quantiles(sorted_rs, n=4, method="inclusive")
        rank_iqr = q3 - q1
        consistency = sum(1 for r in rs if r <= top_k) / len(rs)
        rank_stability = 1.0 - rank_iqr / n

        contributions = [base[name] * values[name] for name in names]
        total_c = sum(contributions)
        if total_c > 0:
            dominance = sum((c / total_c) ** 2 for c in contributions)
            effective_chains = 1.0 / dominance if dominance > 0 else 0.0
        else:
            dominance = 1.0
            effective_chains = 0.0
        dominance_factor = min(1.0, effective_chains / 2.0)

        expected_score = sum(ss) / len(ss)
        adjusted = expected_score * dominance_factor * rank_stability
        spread = max(rs) - min(rs)
        verdict = "rank_stable" if spread == 0 else (
            "rank_minor_shift" if spread <= minor_band else "rank_unstable")
        out.append({
            "peptide_id": row.get("peptide_id", ""),
            "peptide": row.get("peptide", ""),
            "hla_allele": row.get("hla_allele", ""),
            "mhc_class": row.get("mhc_class", ""),
            "gene": row.get("gene", ""),
            "base_score": f"{ss[0]:.4f}",
            "base_rank": str(rs[0]),
            "expected_rank": f"{expected_rank:.2f}",
            "median_rank": f"{median(rs):.1f}",
            "rank_iqr": f"{rank_iqr:.1f}",
            "best_rank": str(min(rs)),
            "worst_rank": str(max(rs)),
            "rank_spread": str(spread),
            "score_min": f"{min(ss):.4f}",
            "score_max": f"{max(ss):.4f}",
            "score_range": f"{max(ss) - min(ss):.4f}",
            "expected_score": f"{expected_score:.4f}",
            "consistency": f"{consistency:.4f}",
            "dominance": f"{dominance:.4f}",
            "effective_chains": f"{effective_chains:.2f}",
            "dominance_factor": f"{dominance_factor:.4f}",
            "rank_stability": f"{rank_stability:.4f}",
            "verdict": verdict,
            "adjusted_score": f"{adjusted:.4f}",
            "adjusted_rank": "",  # filled after the loop, once every value exists
            "top_selection": "",  # filled after base ranks are final
            "n_draws": str(len(vectors)),
            "n_chains": str(sum(1 for wt in base.values() if wt > 0)),
            "n_scorable": str(n),
            "concentration": f"{concentration:.2f}",
            "top_k": str(top_k),
            "seed": str(seed),
            "_adjusted_raw": adjusted,
            "_shares": [c / total_c for c in contributions] if total_c > 0 else None,
        })

    adjusted_ranks = _dense_ranks([item["_adjusted_raw"] for item in out])
    for item, adj_rank in zip(out, adjusted_ranks):
        item["adjusted_rank"] = str(adj_rank)

    # Reporting window: every candidate inside the default top-K or inside the
    # median top-K (ties in a median of .5 use <= so the count stays >= K when
    # available). The rank pool above always covered ALL candidates.
    window = [item for item in out
              if int(item["base_rank"]) <= top_k or float(item["median_rank"]) <= top_k]
    window.sort(key=lambda item: int(item["base_rank"]))
    for item in window:
        item["top_selection"] = (
            "yes" if (int(item["base_rank"]) <= top_k
                      and float(item["consistency"]) >= 1.0 - 1e-12) else "no")

    shares = [item.pop("_shares") for item in window if item.get("_shares")]
    for item in out:
        item.pop("_adjusted_raw", None)
        item.pop("_shares", None)
    chain_diag = []
    n_window = len(shares)
    for idx, (name, _col, wt) in enumerate(chains):
        mean_share = sum(sh[idx] for sh in shares) / n_window if n_window else 0.0
        chain_diag.append({
            "chain": name,
            "default_weight": f"{wt:.4f}",
            "mean_share": f"{mean_share:.4f}",
            "n_candidates": str(n_window),
        })
    return window, chain_diag


def write_weight_sensitivity(peptide_rows, profile, path,
                             samples: int = SENSITIVITY_DEFAULT_SAMPLES,
                             concentration: float = SENSITIVITY_DEFAULT_CONCENTRATION,
                             top_k: int = SENSITIVITY_DEFAULT_TOP_K,
                             seed: int = SENSITIVITY_DEFAULT_SEED):
    """Write weight_sensitivity.tsv plus the per-chain diagnostic alongside."""
    rows, chain_diag = weight_sensitivity_table(
        peptide_rows, profile, samples=samples, concentration=concentration,
        top_k=top_k, seed=seed)
    path = Path(path)
    write_tsv(path, rows, WEIGHT_SENSITIVITY_FIELDS)
    chains_path = path.with_name(f"{path.stem}_chains{path.suffix}")
    write_tsv(chains_path, chain_diag, CHAIN_DIAG_FIELDS)
    return path, chains_path
