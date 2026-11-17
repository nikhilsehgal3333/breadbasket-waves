"""Singular spectrum analysis (SSA) trends for detrending yields.

Kornhuber et al. (2023), Methods: "We then simply detrend modeled and observed crop yields using
singular spectrum analysis (SSA) ... This isolates interannual yield variability from longer-term
trends, and thus quantifies yield variations relative to ~10-year moving windows."

The paper gives neither the window nor the number of components. The study plan fixes a window
of L = 10 years and the first reconstructed component as the trend. A series can be split at
break years (EEU at 1992), and each segment is detrended on its own.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import params


def _hankel_average(m: np.ndarray) -> np.ndarray:
    rows, cols = m.shape
    out = np.zeros(rows + cols - 1)
    count = np.zeros(rows + cols - 1)
    for i in range(rows):
        out[i:i + cols] += m[i]
        count[i:i + cols] += 1
    return out / count


def components(x: np.ndarray, window: int) -> np.ndarray:
    """All reconstructed components, shape (n_components, len(x)); they sum back to x."""
    x = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError("SSA needs a complete series")
    if not 2 <= window <= len(x) - 1:
        raise ValueError(f"window {window} does not fit a series of length {len(x)}")
    trajectory = np.lib.stride_tricks.sliding_window_view(x, window).T     # (window, K)
    u, s, vt = np.linalg.svd(trajectory, full_matrices=False)
    return np.array([_hankel_average(s[i] * np.outer(u[:, i], vt[i])) for i in range(len(s))])


def trend(x: np.ndarray, window: int = params.SSA_WINDOW, k: int = params.SSA_COMPONENTS) -> np.ndarray:
    return components(x, window)[:k].sum(axis=0)


def detrend(series: pd.Series, window: int = params.SSA_WINDOW, k: int = params.SSA_COMPONENTS,
            breaks: tuple[int, ...] = ()) -> pd.Series:
    """Year-indexed series -> year-indexed trend. Each break year starts a new segment."""
    series = series.sort_index()
    cuts = [series.index.min()] + [b for b in breaks if series.index.min() < b <= series.index.max()]
    cuts.append(series.index.max() + 1)
    parts = []
    for start, stop in zip(cuts[:-1], cuts[1:]):
        seg = series[(series.index >= start) & (series.index < stop)]
        parts.append(pd.Series(trend(seg.values, window, k), index=seg.index))
    return pd.concat(parts)


def moving_average_trend(series: pd.Series, span: int = 11, breaks: tuple[int, ...] = ()) -> pd.Series:
    """Sensitivity: centred moving average, shrinking at the ends, segment by segment."""
    series = series.sort_index()
    cuts = [series.index.min()] + [b for b in breaks if series.index.min() < b <= series.index.max()]
    cuts.append(series.index.max() + 1)
    parts = []
    for start, stop in zip(cuts[:-1], cuts[1:]):
        seg = series[(series.index >= start) & (series.index < stop)]
        parts.append(seg.rolling(span, center=True, min_periods=1).mean())
    return pd.concat(parts)
