import numpy as np

from breadbasket_waves import waves
from breadbasket_waves.fetch import NCEP_LATS, NCEP_LONS


def _zonal_wave(k, amplitude, phase=0.3):
    return amplitude * np.cos(2 * np.pi * k * NCEP_LONS / 360.0 + phase)


def test_fft_recovers_known_wavenumber_amplitude():
    x = _zonal_wave(5, 7.5)[None, :]
    amps = waves.wave_amplitudes(x, ks=(5, 7))
    assert np.isclose(amps[0, 0], 7.5)
    assert amps[0, 1] < 1e-9


def test_two_waves_are_separated():
    x = (_zonal_wave(5, 3.0) + _zonal_wave(7, 4.0, phase=1.1))[None, :]
    amps = waves.wave_amplitudes(x, ks=(5, 7))
    assert np.allclose(amps[0], [3.0, 4.0])


def test_lat_weighted_mean_of_constant_is_constant():
    v = np.full((92, len(NCEP_LATS), len(NCEP_LONS)), 4.2)
    assert np.allclose(waves.lat_weighted_mean(v), 4.2)


def test_lat_weights_favour_lower_latitudes():
    v = np.zeros((1, len(NCEP_LATS), 1))
    v[0, -1, 0] = 1.0          # 37.5 N
    low = waves.lat_weighted_mean(v)[0, 0]
    v[0, -1, 0], v[0, 0, 0] = 0.0, 1.0   # 57.5 N
    high = waves.lat_weighted_mean(v)[0, 0]
    assert low > high


def test_weekly_means_give_13_weeks():
    x = np.arange(92 * 2, dtype=float).reshape(92, 2)
    weekly = waves.weekly_means(x)
    assert weekly.shape == (13, 2)
    assert np.allclose(weekly[0], x[:7].mean(axis=0))
    assert np.allclose(weekly[12], x[84:91].mean(axis=0))   # 31 August (day 92) is dropped


def test_threshold_is_mean_plus_1p5_sd():
    rng = np.random.default_rng(0)
    amps = {y: rng.gamma(4.0, 1.0, size=(13, 2)) for y in range(1960, 2015)}
    pooled = np.concatenate([amps[y] for y in range(1960, 2015)])
    thr = waves.thresholds(amps, base=(1960, 2014))
    assert np.allclose(thr, pooled.mean(axis=0) + 1.5 * pooled.std(axis=0, ddof=1))


def test_scaling_the_wind_leaves_events_unchanged():
    rng = np.random.default_rng(1)
    amps = {y: rng.gamma(4.0, 1.0, size=(13, 2)) for y in range(1960, 2015)}
    events = waves.event_table(amps, waves.thresholds(amps, base=(1960, 2014)))
    scaled = {y: 3.0 * a for y, a in amps.items()}
    events_scaled = waves.event_table(scaled, waves.thresholds(scaled, base=(1960, 2014)))
    assert (events["event"].values == events_scaled["event"].values).all()


def test_adjacent_event_weeks_count_as_two():
    amps = {2000: np.zeros((13, 2)), 2001: np.zeros((13, 2))}
    amps[2000][3, 1] = amps[2000][4, 1] = 10.0          # two adjacent wave-7 weeks
    events = waves.event_table(amps, thr=np.array([5.0, 5.0]))
    counts = waves.year_counts(events)
    assert counts.loc[2000, 7] == 2
    assert counts.loc[2001, 7] == 0


def test_single_event_years_are_excluded():
    assert waves.year_class(0) == "zero"
    assert waves.year_class(1) == "one"
    assert waves.year_class(2) == "multi"
    assert waves.year_class(1, minimum=1) == "multi"


def test_week_months():
    assert waves.week_month(2018, 1) == 6
    assert waves.week_month(2018, 5) == 7      # 29 June to 5 July, middle day 2 July
    assert waves.week_month(2018, 13) == 8
