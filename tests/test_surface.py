import numpy as np
import pandas as pd
import xarray as xr

from breadbasket_waves import surface


def _field(values_by_time, times):
    return xr.DataArray(np.asarray(values_by_time, dtype=float), dims=("time", "lat", "lon"),
                        coords={"time": times, "lat": [30.0, 40.0], "lon": [0.0, 1.0]})


def test_per_cell_detrend_removes_trend():
    times = pd.to_datetime([f"{y}-{m:02d}-16" for y in range(1961, 2015) for m in (6, 7, 8)])
    years = times.year.to_numpy()
    vals = 20 + 0.03 * (years - 1961)[:, None, None] + np.zeros((len(times), 2, 2))
    vals[times.month == 7] += 5.0                       # a different level per calendar month
    anom = surface.detrend_linear(_field(vals, times), base=(1961, 2014))
    assert np.allclose(anom.values, 0.0, atol=1e-9)


def test_composite_masks_cells_without_stations():
    times = pd.to_datetime(["2000-06-16", "2000-07-16", "2001-06-16"])
    anom = _field(np.stack([np.full((2, 2), 1.0), np.full((2, 2), 3.0), np.full((2, 2), 9.0)]), times)
    stn = _field(np.ones((3, 2, 2)), times)
    stn[1, 0, 0] = 0                                    # no station in July 2000 at one cell
    comp = surface.composite_maps(anom, stn, [(2000, 6), (2000, 7)])
    assert comp.values[0, 0] == 1.0                     # July 2000 masked there
    assert comp.values[1, 1] == 2.0


def test_event_months_are_unique_pairs():
    ev = pd.DataFrame({"year": [2003, 2003, 2003, 2020], "week": [1, 2, 6, 1], "month": [6, 6, 7, 6],
                       "wavenumber": [7, 7, 7, 7], "amplitude": 1.0, "event": [True, True, True, True]})
    assert surface.event_months(ev, 7, last=2014) == [(2003, 6), (2003, 7)]
