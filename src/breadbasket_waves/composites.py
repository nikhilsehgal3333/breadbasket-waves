"""Composite yield anomalies for years with more than one wave event and years with none.

Kornhuber et al. (2023), Methods, "Composite analysis": "the values are divided by the average
yield values for the three years preceding and following the event, generating standardized
composites (in terms of % of regional production)", with "500 replicates of the composites using
~90% of the full sample". The study reports that spread as a subsample spread, not a confidence
interval, and adds a stratified bootstrap interval with replacement.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import params


def neighbour_ratio(series: pd.Series, half: int = params.NEIGHBOUR_YEARS,
                    breaks: tuple[int, ...] = ()) -> pd.Series:
    """Percent difference from the mean of the `half` years before and after (the year itself
    is left out). NaN when the window leaves the series or crosses a break year."""
    series = series.sort_index()
    years = series.index.to_numpy()
    out = pd.Series(np.nan, index=series.index, dtype=float)
    for year in years:
        window = [y for y in range(year - half, year + half + 1) if y != year]
        if any(y not in series.index for y in window):
            continue
        if any(year - half < b <= year + half for b in breaks):
            continue
        base = series.loc[window].mean()
        out[year] = 100.0 * (series[year] / base - 1.0)
    return out


def composite(std: pd.Series, classes: pd.Series) -> dict[str, float]:
    cls = classes.reindex(std.index)
    ok = std.notna()
    multi = std[ok & (cls == "multi")]
    zero = std[ok & (cls == "zero")]
    return {"multi_mean": float(multi.mean()), "zero_mean": float(zero.mean()),
            "n_multi": int(multi.size), "n_zero": int(zero.size)}


def subsample_spread(std: pd.Series, classes: pd.Series, rng: np.random.Generator,
                     frac: float = params.SUBSAMPLE_FRACTION, n: int = params.N_SUBSAMPLES) -> dict[str, float]:
    """The paper's resampling: n draws of about 90 % of each sample, without replacement."""
    cls = classes.reindex(std.index)
    out = {}
    for name in ("multi", "zero"):
        values = std[std.notna() & (cls == name)].to_numpy()
        size = max(1, int(round(frac * values.size)))
        draws = [rng.choice(values, size=size, replace=False).mean() for _ in range(n)]
        out[f"{name}_median"] = float(np.median(draws))
        out[f"{name}_q05"], out[f"{name}_q95"] = (float(q) for q in np.quantile(draws, [0.05, 0.95]))
    return out


def bootstrap_interval(std: pd.Series, classes: pd.Series, rng: np.random.Generator,
                       n: int = params.N_BOOTSTRAP) -> dict[str, float]:
    """Stratified bootstrap with replacement: 95 % interval of each class mean and of the gap."""
    cls = classes.reindex(std.index)
    multi = std[std.notna() & (cls == "multi")].to_numpy()
    zero = std[std.notna() & (cls == "zero")].to_numpy()
    m = np.array([rng.choice(multi, multi.size).mean() for _ in range(n)])
    z = np.array([rng.choice(zero, zero.size).mean() for _ in range(n)])
    out = {}
    for name, draws in (("multi", m), ("zero", z), ("gap", m - z)):
        out[f"{name}_lo"], out[f"{name}_hi"] = (float(q) for q in np.quantile(draws, [0.025, 0.975]))
    return out
