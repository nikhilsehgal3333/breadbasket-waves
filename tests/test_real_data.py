"""Mechanics checks on the downloaded data. The whole file is skipped when local/ is empty.

These check that the data and the pipeline behave as the method implies. None of them asserts
agreement with the paper: that is reported in RESULTS.md, never tested."""
import numpy as np
import pytest

from breadbasket_waves import fetch, params, waves, yields

HAVE_NCEP = all(fetch.ncep_path(y).exists() for y in range(params.WAVE_BASE[0], params.LAST_YEAR + 1))
HAVE_FAO = fetch.fao_path().exists()
needs_ncep = pytest.mark.skipif(not HAVE_NCEP, reason="NCEP slices not downloaded")
needs_fao = pytest.mark.skipif(not HAVE_FAO, reason="FAOSTAT file not downloaded")


@needs_ncep
@pytest.mark.parametrize("year", [1964, 2000, 2015])          # two leap years and one other
def test_wind_slices_have_92_days_and_the_band(year):
    assert np.load(fetch.ncep_path(year)).shape == (92, len(fetch.NCEP_LATS), len(fetch.NCEP_LONS))


@needs_ncep
def test_event_share_is_what_the_threshold_implies():
    # mean + 1.5 SD leaves roughly 5-11 % of weeks above it for skewed, positive amplitudes.
    amps = {y: waves.year_amplitudes(np.load(fetch.ncep_path(y)))
            for y in range(params.WAVE_BASE[0], params.WAVE_BASE[1] + 1)}
    events = waves.event_table(amps, waves.thresholds(amps))
    share = events.groupby("wavenumber").event.mean()
    assert ((share > 0.05) & (share < 0.11)).all(), share.to_dict()


@needs_fao
def test_region_names_match_this_fao_release():
    fao = yields.load_fao(fetch.fao_path())
    names = fao.drop_duplicates("code").set_index("code").area
    regions = yields.load_regions()
    wrong = {int(r.code): (r["name"], names.get(r.code)) for _, r in regions.iterrows()
             if names.get(r.code) != r["name"]}
    assert not wrong, wrong


@needs_fao
def test_fao_spans_match_the_region_map():
    fao = yields.load_fao(fetch.fao_path())
    wheat = fao[(fao.crop == "wheat") & (fao.element == "production") & fao.value.notna()]
    first = wheat.groupby("code").year.min()
    last = wheat.groupby("code").year.max()
    for _, r in yields.load_regions().iterrows():
        assert first[r.code] == r.first_year, (r["name"], first[r.code])
        assert last[r.code] == r.last_year, (r["name"], last[r.code])


@needs_fao
def test_every_listed_country_reports_wheat_every_year():
    fao = yields.load_fao(fetch.fao_path())
    yields.coverage(fao, yields.load_regions(), range(params.YIELD_BASE[0], params.LAST_YEAR + 1))
