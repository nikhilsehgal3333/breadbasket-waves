"""Likelihood multiplication factor (LMF) of concurrent low (or high) yields in two regions.

Kornhuber et al. (2023), Methods, Eqs. 1-3: p_w>1 = N(both regions below trend | more than one
wave event) / N(more than one wave event); p_0 = the same for years with no wave event;
LMF = p_w>1 / p_0. Concurrent high yields are treated in the same way.

Added by this study, and labelled as its own: per cell, a Katz log 95 % interval and an exact
one-sided hypergeometric p-value; across cells, a permutation test that shuffles whole years of
wave counts against whole years of yields, with four primary tests fixed before any result.
"""
from __future__ import annotations

import math
from itertools import combinations

import numpy as np
import pandas as pd

from . import params
from .waves import year_class

Z95 = 1.959963984540054


def lmf(hits_multi: int, n_multi: int, hits_zero: int, n_zero: int) -> float:
    if n_multi == 0 or n_zero == 0:
        return math.nan
    p1, p0 = hits_multi / n_multi, hits_zero / n_zero
    if p0 == 0:
        return math.inf if p1 > 0 else math.nan
    return p1 / p0


def katz_ci(hits_multi: int, n_multi: int, hits_zero: int, n_zero: int, z: float = Z95) -> tuple[float, float]:
    """Katz log interval for a ratio of two proportions; adds 0.5 to the cells when one is zero."""
    a, n1, c, n2 = float(hits_multi), float(n_multi), float(hits_zero), float(n_zero)
    if n1 == 0 or n2 == 0:
        return (math.nan, math.nan)
    if a == 0 or c == 0:
        a, c, n1, n2 = a + 0.5, c + 0.5, n1 + 1.0, n2 + 1.0
    ratio = (a / n1) / (c / n2)
    se = math.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2)
    return (ratio * math.exp(-z * se), ratio * math.exp(z * se))


def _hypergeom_pmf(x: int, total: int, successes: int, draws: int) -> float:
    return math.comb(successes, x) * math.comb(total - successes, draws - x) / math.comb(total, draws)


def p_greater(hits_multi: int, n_multi: int, hits_zero: int, n_zero: int) -> float:
    """Exact P(X >= hits_multi) for X ~ Hypergeometric(n_multi + n_zero, hits, n_multi)."""
    total, successes, draws = n_multi + n_zero, hits_multi + hits_zero, n_multi
    top = min(successes, draws)
    return float(sum(_hypergeom_pmf(x, total, successes, draws) for x in range(hits_multi, top + 1)))


def p_less(hits_multi: int, n_multi: int, hits_zero: int, n_zero: int) -> float:
    """Exact P(X <= hits_multi), for the test that concurrent high yields are rarer."""
    total, successes, draws = n_multi + n_zero, hits_multi + hits_zero, n_multi
    low = max(0, draws - (total - successes))
    return float(sum(_hypergeom_pmf(x, total, successes, draws) for x in range(low, hits_multi + 1)))


def holm(pvalues: list[float]) -> list[float]:
    order = np.argsort(pvalues)
    m = len(pvalues)
    adjusted = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvalues[i]))
        adjusted[i] = running
    return adjusted.tolist()


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out.tolist()


def classes_from_counts(counts: pd.Series, minimum: int = 2) -> pd.Series:
    return counts.map(lambda n: year_class(int(n), minimum))


def pair_counts(anom: pd.DataFrame, classes: pd.Series, r1: str, r2: str, direction: str) -> tuple[int, int, int, int]:
    years = anom.index.intersection(classes.index)
    a, b, cls = anom.loc[years, r1], anom.loc[years, r2], classes.loc[years]
    ok = a.notna() & b.notna()
    both = ((a < 0) & (b < 0)) if direction == "low" else ((a > 0) & (b > 0))
    multi, zero = ok & (cls == "multi"), ok & (cls == "zero")
    return int((both & multi).sum()), int(multi.sum()), int((both & zero).sum()), int(zero.sum())


def lmf_table(anom: pd.DataFrame, classes: pd.Series, regions: tuple[str, ...] = params.REGIONS) -> pd.DataFrame:
    """One row per region pair and direction, with counts, LMF, interval and exact p-value."""
    rows = []
    for r1, r2 in combinations(regions, 2):
        for direction in ("low", "high"):
            a, n1, c, n2 = pair_counts(anom, classes, r1, r2, direction)
            lo, hi = katz_ci(a, n1, c, n2)
            p = p_greater(a, n1, c, n2) if direction == "low" else p_less(a, n1, c, n2)
            rows.append((f"{r1} x {r2}", r1, r2, direction, a, n1, c, n2, lmf(a, n1, c, n2), lo, hi, p))
    return pd.DataFrame(rows, columns=["pair", "region_1", "region_2", "direction", "hits_multi",
                                       "n_multi", "hits_zero", "n_zero", "lmf", "ci_low", "ci_high",
                                       "p_exact"])


def _pair_masks(anom: pd.DataFrame, years: np.ndarray, regions: tuple[str, ...]):
    pairs = list(combinations(regions, 2))
    low = np.zeros((len(pairs), len(years)), dtype=bool)
    high = np.zeros_like(low)
    valid = np.zeros_like(low)
    for i, (r1, r2) in enumerate(pairs):
        a, b = anom.loc[years, r1].to_numpy(), anom.loc[years, r2].to_numpy()
        valid[i] = np.isfinite(a) & np.isfinite(b)
        low[i] = valid[i] & (a < 0) & (b < 0)
        high[i] = valid[i] & (a > 0) & (b > 0)
    return low, high, valid


def _lmf_vector(both: np.ndarray, valid: np.ndarray, multi: np.ndarray, zero: np.ndarray) -> np.ndarray:
    n1 = (valid & multi).sum(axis=1).astype(float)
    n0 = (valid & zero).sum(axis=1).astype(float)
    h1 = (both & multi).sum(axis=1).astype(float)
    h0 = (both & zero).sum(axis=1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p1 = np.where(n1 > 0, h1 / n1, np.nan)
        p0 = np.where(n0 > 0, h0 / n0, np.nan)
        out = np.where(p0 > 0, p1 / p0, np.where(p1 > 0, np.inf, np.nan))
    return out


def primary_statistics(anom: pd.DataFrame, wave_counts: pd.DataFrame, years: np.ndarray,
                       regions: tuple[str, ...] = params.REGIONS, minimum: int = 2) -> dict[str, int]:
    """Per wave: the number of pairs with LMF(low) > 1 and with LMF(high) < 1."""
    low, high, valid = _pair_masks(anom, years, regions)
    out = {}
    for k in wave_counts.columns:
        n = wave_counts.loc[years, k].to_numpy()
        multi, zero = n >= minimum, n == 0
        out[f"w{k}_low_gt1"] = int(np.sum(_lmf_vector(low, valid, multi, zero) > 1))
        out[f"w{k}_high_lt1"] = int(np.sum(_lmf_vector(high, valid, multi, zero) < 1))
    return out


def permutation_test(anom: pd.DataFrame, wave_counts: pd.DataFrame, years: np.ndarray,
                     n_perm: int = params.N_PERMUTATIONS, rng: np.random.Generator | None = None,
                     regions: tuple[str, ...] = params.REGIONS, minimum: int = 2) -> pd.DataFrame:
    """Shuffle whole years of wave counts (all wavenumbers together) against whole years of
    yields (all regions together). p = (1 + #{T* >= T}) / (1 + B). Holm-adjusted over the four."""
    rng = rng or params.rng_for("permutation")
    observed = primary_statistics(anom, wave_counts, years, regions, minimum)
    low, high, valid = _pair_masks(anom, years, regions)
    counts = wave_counts.loc[years].to_numpy()            # (year, wavenumber)
    exceed = {key: 0 for key in observed}
    for _ in range(n_perm):
        shuffled = counts[rng.permutation(len(years))]
        for j, k in enumerate(wave_counts.columns):
            multi, zero = shuffled[:, j] >= minimum, shuffled[:, j] == 0
            if int(np.sum(_lmf_vector(low, valid, multi, zero) > 1)) >= observed[f"w{k}_low_gt1"]:
                exceed[f"w{k}_low_gt1"] += 1
            if int(np.sum(_lmf_vector(high, valid, multi, zero) < 1)) >= observed[f"w{k}_high_lt1"]:
                exceed[f"w{k}_high_lt1"] += 1
    keys = list(observed)
    p = [(1 + exceed[key]) / (1 + n_perm) for key in keys]
    return pd.DataFrame({"test": keys, "observed_pairs": [observed[k] for k in keys],
                         "of_pairs": len(list(combinations(regions, 2))), "p_perm": p,
                         "p_holm": holm(p)})
