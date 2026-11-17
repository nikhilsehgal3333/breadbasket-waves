"""Wave events from the 250 hPa meridional wind.

Follows Kornhuber et al. (2023), Methods, "Wave analysis and wave event definition": weekly means
of the mid-latitude (37.5-57.5 N) mean of v250 in June-August, decomposed by a fast Fourier
transform along longitude. A week is an event for wavenumber k when its amplitude exceeds the mean
plus 1.5 standard deviations of all JJA weeks of k in the base period (1960-2014).

Choices the paper does not state, fixed in the study plan: cos-latitude weights, amplitude
2|X_k|/N, weeks from 1 June (days 1-91; 31 August is dropped), thresholds pooled over all JJA
weeks of each wavenumber, standard deviation with ddof = 1. Consecutive event weeks are counted
separately, never merged.
"""
from __future__ import annotations

import datetime as dt

import numpy as np
import pandas as pd

from . import params
from .fetch import NCEP_LATS


def lat_weighted_mean(v: np.ndarray, lats: np.ndarray = NCEP_LATS) -> np.ndarray:
    """(day, lat, lon) -> (day, lon), weighted by cos(latitude)."""
    w = np.cos(np.deg2rad(np.asarray(lats, dtype=float)))
    return np.tensordot(np.asarray(v, dtype=float), w / w.sum(), axes=([1], [0]))


def weekly_means(x: np.ndarray, n_weeks: int = params.N_WEEKS) -> np.ndarray:
    """(day, lon) -> (week, lon): non-overlapping 7-day means of the first n_weeks * 7 days."""
    days = n_weeks * 7
    if x.shape[0] < days:
        raise ValueError(f"need {days} days, got {x.shape[0]}")
    return x[:days].reshape(n_weeks, 7, x.shape[1]).mean(axis=1)


def wave_amplitudes(x: np.ndarray, ks: tuple[int, ...] = params.WAVENUMBERS) -> np.ndarray:
    """(..., lon) -> (..., len(ks)): amplitude 2|X_k| / N of the zonal Fourier series."""
    n = x.shape[-1]
    spectrum = np.fft.rfft(x, axis=-1)
    return np.stack([2.0 * np.abs(spectrum[..., k]) / n for k in ks], axis=-1)


def year_amplitudes(v: np.ndarray) -> np.ndarray:
    """(day, lat, lon) wind for one JJA -> (week, k) amplitudes."""
    return wave_amplitudes(weekly_means(lat_weighted_mean(v)))


def thresholds(amps: dict[int, np.ndarray], base: tuple[int, int] = params.WAVE_BASE,
               n_sd: float = params.THRESHOLD_SD) -> np.ndarray:
    """amps maps year -> (week, k). Return mean + n_sd * SD per k over all base-period weeks."""
    pooled = np.concatenate([amps[y] for y in range(base[0], base[1] + 1)], axis=0)
    return pooled.mean(axis=0) + n_sd * pooled.std(axis=0, ddof=1)


def week_dates(year: int, week: int) -> tuple[dt.date, dt.date]:
    """First and last day of week 1..13 of JJA."""
    start = dt.date(year, 6, 1) + dt.timedelta(days=7 * (week - 1))
    return start, start + dt.timedelta(days=6)


def week_month(year: int, week: int) -> int:
    """Calendar month of the week's middle day (used to pick the CRU month)."""
    start, _ = week_dates(year, week)
    return (start + dt.timedelta(days=3)).month


def event_table(amps: dict[int, np.ndarray], thr: np.ndarray,
                ks: tuple[int, ...] = params.WAVENUMBERS) -> pd.DataFrame:
    """Long table with one row per year, week and wavenumber."""
    rows = []
    for year in sorted(amps):
        a = amps[year]
        for w in range(a.shape[0]):
            for j, k in enumerate(ks):
                rows.append((year, w + 1, week_month(year, w + 1), k, float(a[w, j]),
                             bool(a[w, j] > thr[j])))
    return pd.DataFrame(rows, columns=["year", "week", "month", "wavenumber", "amplitude", "event"])


def year_counts(events: pd.DataFrame) -> pd.DataFrame:
    """Rows: year. Columns: wavenumber. Values: number of event weeks."""
    table = events.pivot_table(index="year", columns="wavenumber", values="event", aggfunc="sum")
    return table.astype(int).sort_index()


def year_class(n_event_weeks: int, minimum: int = 2) -> str:
    """'multi' for at least `minimum` event weeks, 'zero' for none, 'one' otherwise (excluded)."""
    if n_event_weeks >= minimum:
        return "multi"
    return "zero" if n_event_weeks == 0 else "one"
