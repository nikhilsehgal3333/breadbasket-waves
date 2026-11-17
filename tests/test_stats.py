import math

import numpy as np
import pandas as pd
import pytest

from breadbasket_waves import composites, lmf, params, ssa


# --- SSA -----------------------------------------------------------------------------------

def test_ssa_components_sum_back_to_the_series():
    x = np.random.default_rng(0).normal(size=40).cumsum()
    assert np.allclose(ssa.components(x, 10).sum(axis=0), x)


def test_ssa_recovers_linear_trend():
    years = np.arange(1961, 2015)
    x = 2.0 + 0.05 * (years - 1961)
    assert np.allclose(ssa.trend(x, 10, 2), x, atol=1e-8)
    assert params.SSA_COMPONENTS == 2


def test_one_component_leaves_a_trend_residual():
    # Why the plan uses k = 2: with k = 1 a plain linear trend is not removed, and the residual
    # (up to about 2.5 % of the level here) is as large as the effects being measured.
    x = 2.0 + 0.05 * np.arange(54)
    residual = np.abs(x - ssa.trend(x, 10, 1)).max() / x.mean()
    assert residual > 0.01


def test_ssa_constant_stays_constant_and_scales():
    x = np.full(30, 3.0)
    assert np.allclose(ssa.trend(x), 3.0)
    y = np.random.default_rng(1).normal(5, 1, 30)
    assert np.allclose(ssa.trend(7 * y), 7 * ssa.trend(y))


def test_detrend_splits_at_breaks():
    s = pd.Series(np.r_[np.full(31, 1.0), np.full(33, 5.0)], index=range(1961, 2025))
    t = ssa.detrend(s, breaks=(1992,))
    assert np.allclose(t.loc[:1991], 1.0) and np.allclose(t.loc[1992:], 5.0)


# --- LMF -----------------------------------------------------------------------------------

def test_worked_example_gives_1_6():
    # Kornhuber et al. (2023), Methods: 6 of 15 multi-event years, 5 of 20 zero-event years.
    assert lmf.lmf(6, 15, 5, 20) == pytest.approx(1.6)


def test_worked_example_interval_and_p():
    lo, hi = lmf.katz_ci(6, 15, 5, 20)
    assert (round(lo, 2), round(hi, 2)) == (0.60, 4.26)
    assert lmf.p_greater(6, 15, 5, 20) == pytest.approx(0.281, abs=0.002)
    # With the control fixed at 5/20, a pair needs 9/15 to pass 0.05.
    assert lmf.p_greater(9, 15, 5, 20) < 0.05 < lmf.p_greater(8, 15, 5, 20)


def test_zero_denominator_gives_nan_and_zero_control_gives_inf():
    assert math.isnan(lmf.lmf(1, 0, 2, 10))
    assert math.isinf(lmf.lmf(2, 10, 0, 10))
    assert math.isnan(lmf.lmf(0, 10, 0, 10))
    lo, hi = lmf.katz_ci(2, 10, 0, 10)
    assert math.isfinite(lo) and math.isfinite(hi)


def test_hypergeometric_tails_sum_to_one_plus_point():
    a, n1, c, n2 = 4, 12, 6, 20
    point = lmf._hypergeom_pmf(a, n1 + n2, a + c, n1)
    assert lmf.p_greater(a, n1, c, n2) + lmf.p_less(a, n1, c, n2) == pytest.approx(1 + point)


def test_holm_and_bh():
    assert lmf.holm([0.01, 0.04, 0.03, 0.2]) == pytest.approx([0.04, 0.09, 0.09, 0.2])
    q = lmf.benjamini_hochberg([0.01, 0.02, 0.03, 0.5])
    assert q == pytest.approx([0.04, 0.04, 0.04, 0.5])


def _toy(seed=0, n_years=54):
    rng = np.random.default_rng(seed)
    years = np.arange(1961, 1961 + n_years)
    anom = pd.DataFrame(rng.normal(size=(n_years, 5)), index=years, columns=params.REGIONS)
    counts = pd.DataFrame({5: rng.poisson(1.1, n_years), 7: rng.poisson(1.0, n_years)}, index=years)
    return anom, counts, years


def test_lmf_table_has_ten_pairs_per_direction():
    anom, counts, years = _toy()
    table = lmf.lmf_table(anom, lmf.classes_from_counts(counts[7]))
    assert len(table) == 20 and set(table.direction) == {"low", "high"}


def test_permutation_keeps_class_sizes_and_p_is_positive():
    anom, counts, years = _toy(1)
    result = lmf.permutation_test(anom, counts, years, n_perm=200, rng=np.random.default_rng(5))
    assert (result.p_perm > 0).all() and (result.p_perm <= 1).all()
    shuffled = counts.to_numpy()[np.random.default_rng(5).permutation(len(years))]
    assert sorted(shuffled[:, 1]) == sorted(counts[7])          # same class sizes


def test_primary_statistics_match_the_table():
    anom, counts, years = _toy(2)
    stats = lmf.primary_statistics(anom, counts, years)
    table = lmf.lmf_table(anom, lmf.classes_from_counts(counts[7]))
    assert stats["w7_low_gt1"] == int((table[table.direction == "low"].lmf > 1).sum())
    assert stats["w7_high_lt1"] == int((table[table.direction == "high"].lmf < 1).sum())


# --- Composites ----------------------------------------------------------------------------

def test_ratio_to_neighbours_on_toy_series():
    s = pd.Series(100.0, index=range(1990, 2001))
    s[1995] = 90.0
    r = composites.neighbour_ratio(s)
    assert r[1995] == pytest.approx(-10.0)
    assert np.isnan(r[1990]) and np.isnan(r[2000])          # windows leave the series
    assert np.isnan(r[1992])                                 # needs 1989, before the series
    assert r[1993] == pytest.approx(100 * (100 / ((100 * 5 + 90) / 6) - 1))


def test_ratio_windows_do_not_cross_breaks():
    s = pd.Series(1.0, index=range(1980, 2005))
    r = composites.neighbour_ratio(s, breaks=(1992,))
    assert r.loc[1989:1994].isna().all() and r.loc[1983:1988].notna().all()


def test_composite_means_by_class():
    std = pd.Series([-5.0, -3.0, 1.0, 2.0, np.nan], index=[1, 2, 3, 4, 5])
    cls = pd.Series(["multi", "multi", "zero", "zero", "multi"], index=[1, 2, 3, 4, 5])
    c = composites.composite(std, cls)
    assert c == {"multi_mean": -4.0, "zero_mean": 1.5, "n_multi": 2, "n_zero": 2}
