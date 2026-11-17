"""Surface maps from CRU TS 4.10: maximum temperature and precipitation in the months that
contain wave event weeks.

Simplified relative to Kornhuber et al. (2023), which used weekly ERA5 anomalies centred on each
event: CRU TS is monthly, so each event week is assigned to the month of its middle day, and the
composite averages those month-fields. Anomalies are taken from a per-cell, per-calendar-month
linear trend fitted on 1961-2014 (the paper also removes a linear trend per grid point). Cells
with no station within reach in that month (stn = 0) are masked.

Only the figure is published. The CRU-derived grids stay in local/ (Open Database License).
"""
from __future__ import annotations

import gzip
import shutil

import numpy as np
import pandas as pd
import xarray as xr

from . import fetch, params
from .paths import LOCAL

LAT_RANGE = (20.0, 75.0)
MONTHS = (6, 7, 8)


def _unzipped(var: str, a: int, b: int):
    """Unpack one CRU file to a scratch path; the caller deletes it after reading (about 370 MB
    each unpacked, so keeping all 16 would cost about 6 GB)."""
    src = fetch.cru_path(var, a, b)
    out = LOCAL / "cru" / "unz" / src.name[:-3]
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    with gzip.open(src) as f, open(tmp, "wb") as g:
        shutil.copyfileobj(f, g)
    tmp.replace(out)
    return out


def load(var: str, first: int = params.WAVE_BASE[0], last: int = params.LAST_YEAR) -> tuple[xr.DataArray, xr.DataArray]:
    """JJA fields for 20-75 N: (values, station counts), time = month starts."""
    vals, stns = [], []
    for a, b in fetch.CRU_DECADES:
        if b < first or a > last:
            continue
        path = _unzipped(var, a, b)
        try:
            with xr.open_dataset(path) as ds:
                sel = ds.sel(lat=slice(*LAT_RANGE))
                sel = sel.isel(time=np.isin(sel.time.dt.month, MONTHS) &
                               (sel.time.dt.year >= first) & (sel.time.dt.year <= last))
                vals.append(sel[var].load())
                stns.append(sel["stn"].load())
        finally:
            path.unlink(missing_ok=True)
    v, s = xr.concat(vals, "time"), xr.concat(stns, "time")
    return v, s


def detrend_linear(field: xr.DataArray, base: tuple[int, int] = params.YIELD_BASE) -> xr.DataArray:
    """Anomaly from a per-cell, per-calendar-month linear trend fitted on the base years."""
    out = []
    for m in MONTHS:
        f = field.isel(time=field.time.dt.month == m)
        years = f.time.dt.year.values.astype(float)
        in_base = (years >= base[0]) & (years <= base[1])
        x = years[in_base] - years[in_base].mean()
        y = f.values[in_base]
        slope = np.nansum(x[:, None, None] * (y - np.nanmean(y, axis=0)), axis=0) / np.sum(x ** 2)
        intercept = np.nanmean(y, axis=0)
        fit = intercept + slope * (years - years[in_base].mean())[:, None, None]
        out.append(f - fit)
    return xr.concat(out, "time").sortby("time")


def composite_maps(anom: xr.DataArray, stn: xr.DataArray, months: list[tuple[int, int]]) -> xr.DataArray:
    """Mean anomaly over the given (year, month) fields, masking cells with stn = 0."""
    keys = set(months)
    pick = np.array([(int(t.dt.year), int(t.dt.month)) in keys for t in anom.time])
    a = anom.isel(time=pick).where(stn.isel(time=pick) > 0)
    return a.mean("time", skipna=True)


def event_months(events: pd.DataFrame, k: int, last: int = params.YIELD_BASE[1]) -> list[tuple[int, int]]:
    e = events[(events.wavenumber == k) & events.event & (events.year <= last)]
    return sorted(set(zip(e.year.astype(int), e.month.astype(int))))
