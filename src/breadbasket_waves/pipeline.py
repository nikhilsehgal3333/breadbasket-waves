"""The analysis steps, in plan order. Each step reads verified inputs from local/ and writes
byte-stable outputs to results/, and records its key numbers in the manifest."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import fetch, figures, params, waves, yields
from .io_utils import Manifest, write_csv
from .paths import DATA as DATA_DIR, RESULTS


# --- S3: wave events ---------------------------------------------------------------------------

def wave_step(manifest: Manifest) -> dict:
    years = range(params.WAVE_BASE[0], params.LAST_YEAR + 1)
    amps = {y: waves.year_amplitudes(np.load(fetch.ncep_path(y))) for y in years}
    thr = waves.thresholds(amps)
    events = waves.event_table(amps, thr)
    counts = waves.year_counts(events)
    write_csv(events, RESULTS / "wave_weeks.csv", ["year", "week", "wavenumber"])
    out_counts = counts.reset_index().rename(columns={5: "wave5_event_weeks", 7: "wave7_event_weeks"})
    for k in params.WAVENUMBERS:
        out_counts[f"wave{k}_class"] = counts[k].map(waves.year_class).to_numpy()
    write_csv(out_counts, RESULTS / "wave_year_counts.csv", ["year"])
    write_csv(pd.DataFrame({"wavenumber": list(params.WAVENUMBERS), "threshold_m_s": thr}),
              RESULTS / "wave_thresholds.csv", ["wavenumber"])

    base = counts.loc[params.WAVE_BASE[0]:params.WAVE_BASE[1]]
    n_weeks = len(base) * params.N_WEEKS
    src = "results/wave_year_counts.csv"
    manifest.add("base_weeks", n_weeks, src, "JJA weeks 1960-2014")
    for k in params.WAVENUMBERS:
        manifest.add(f"events_w{k}_base", int(base[k].sum()), src, f"wave-{k} event weeks, 1960-2014")
        manifest.add(f"event_share_w{k}_base", round(100 * base[k].sum() / n_weeks, 1), src, "% of weeks")
        cls = base[k].loc[params.YIELD_BASE[0]:].map(waves.year_class)
        for c in ("multi", "zero", "one"):
            manifest.add(f"years_{c}_w{k}_1961_2014", int((cls == c).sum()), src, f"wave-{k} year class")
        ext = counts.loc[params.WAVE_BASE[1] + 1:, k]
        manifest.add(f"events_w{k}_2015_2024", int(ext.sum()), src, f"wave-{k} event weeks, 2015-2024")
        manifest.add(f"threshold_w{k}", round(float(thr[list(params.WAVENUMBERS).index(k)]), 3),
                     "results/wave_thresholds.csv", "m/s")

    # Events per decade (the reanalysis rests on fewer observations before the satellite era, 1979).
    decade = (counts.index // 10) * 10
    per_decade = counts.groupby(decade).sum()
    per_decade.index.name = "decade"
    per_decade = per_decade.reset_index().rename(columns={5: "wave5_event_weeks", 7: "wave7_event_weeks"})
    per_decade["years"] = counts.groupby(decade).size().to_numpy()
    write_csv(per_decade, RESULTS / "wave_events_per_decade.csv", ["decade"])

    # Timing: share of event weeks in June (weeks 1-4), the rest of the season.
    for k in params.WAVENUMBERS:
        e = events[(events.wavenumber == k) & events.event & (events.year <= params.WAVE_BASE[1])]
        manifest.add(f"share_june_w{k}", round(100 * (e.month == 6).mean(), 1), "results/wave_weeks.csv",
                     f"% of wave-{k} event weeks (1960-2014) whose middle day is in June")

    # Known case, reported, not tested: summers named by Kornhuber et al. (2019).
    named = pd.DataFrame({"year": list(params.WAVE7_NAMED_SUMMERS),
                          "wave7_event_weeks": [int(counts.loc[y, 7]) for y in params.WAVE7_NAMED_SUMMERS]})
    named["wave7_class"] = named.wave7_event_weeks.map(waves.year_class)
    write_csv(named, RESULTS / "wave7_named_summers.csv", ["year"])
    for _, r in named.iterrows():
        manifest.add(f"named_w7_{int(r.year)}", int(r.wave7_event_weeks), "results/wave7_named_summers.csv")

    figures.fig_wave_events(counts, events, params.WAVE_BASE[1], params.WAVE7_NAMED_SUMMERS)
    return {"counts": counts, "events": events, "thresholds": thr}


# --- S4a: FAO coverage and regional yields ------------------------------------------------------

def yield_step(manifest: Manifest) -> dict:
    regions = yields.load_regions()
    yields.check_successions(regions)
    fao = yields.load_fao(fetch.fao_path())
    years = range(params.YIELD_BASE[0], params.LAST_YEAR + 1)
    names = fao.drop_duplicates("code").set_index("code").area
    mismatched = [(int(r.code), r["name"], names.get(r.code)) for _, r in regions.iterrows()
                  if names.get(r.code) != r["name"]]
    cov = yields.coverage(fao, regions, years)
    write_csv(cov, RESULTS / "fao_coverage.csv", ["region", "year", "code"])
    table = yields.region_table(fao, regions, years)
    write_csv(table, RESULTS / "regional_yields.csv", ["region", "year"])
    variants = {
        "maize_only": yields.region_table(fao, regions, years, crops=("maize",)),
        "wheat_only": yields.region_table(fao, regions, years, crops=("wheat",)),
        "china_351": yields.region_table(fao, regions, years, overrides={
            "EAS": [(351, 1961, 2024), (141, 1961, 2024)]}),
    }
    for name, t in variants.items():
        write_csv(t, RESULTS / f"regional_yields_{name}.csv", ["region", "year"])
    manifest.add("fao_latest_year", int(fao.year.max()), "results/regional_yields.csv")
    manifest.add("fao_name_mismatches", len(mismatched), "data/regions.csv",
                 "; ".join(f"{c}: {a!r} vs {b!r}" for c, a, b in mismatched))
    est = table[table.region != "WORLD"].estimated_share
    manifest.add("estimated_share_max", round(100 * float(est.max()), 1), "results/regional_yields.csv",
                 "% of regional production flagged estimated or imputed, worst region-year")
    return {"fao": fao, "regions": regions, "table": table, "variants": variants}


def spot_checks(fao: pd.DataFrame, regions: pd.DataFrame, table: pd.DataFrame, n: int = 10) -> pd.DataFrame:
    """CP1: recompute n seeded region-years with a plain loop over the raw rows."""
    rng = params.rng_for("cp1-spot-checks")
    pool = table[table.region != "WORLD"][["region", "year"]].to_numpy()
    picks = pool[rng.choice(len(pool), size=n, replace=False)]
    rows = []
    for region, year in picks:
        year = int(year)
        prod = area = 0.0
        codes = yields.active_codes(regions, region, year)
        for code in codes:
            for crop in ("maize", "wheat"):
                p = fao[(fao.code == code) & (fao.crop == crop) & (fao.element == "production") & (fao.year == year)].value
                a = fao[(fao.code == code) & (fao.crop == crop) & (fao.element == "area") & (fao.year == year)].value
                if len(p) == 1 and len(a) == 1 and np.isfinite(p.iloc[0]) and np.isfinite(a.iloc[0]) and a.iloc[0] > 0:
                    prod += float(p.iloc[0])
                    area += float(a.iloc[0])
        pipe = table[(table.region == region) & (table.year == year)].iloc[0]
        rows.append({"region": region, "year": year, "codes": " ".join(map(str, codes)),
                     "loop_production_t": prod, "loop_area_ha": area,
                     "loop_yield_t_ha": prod / area, "pipeline_yield_t_ha": float(pipe.yield_t_ha),
                     "match": bool(np.isclose(prod / area, pipe.yield_t_ha, rtol=1e-12))})
    out = pd.DataFrame(rows)
    write_csv(out, RESULTS / "cp1_fao_spot_checks.csv", ["region", "year"])
    return out


# --- S4b: anomalies and the year audit ----------------------------------------------------------
# (Written before checkpoint 1; run only after it.)

from . import composites, lmf, ssa  # noqa: E402

ALL_REGIONS = params.REGIONS + ("ALL", "WORLD")


def region_series(table: pd.DataFrame) -> pd.DataFrame:
    """Year x region yields, with ALL = the five regions pooled (production over area)."""
    wide_p = table.pivot(index="year", columns="region", values="production_t")
    wide_a = table.pivot(index="year", columns="region", values="area_ha")
    out = (wide_p / wide_a)[list(params.REGIONS) + ["WORLD"]]
    out["ALL"] = wide_p[list(params.REGIONS)].sum(axis=1) / wide_a[list(params.REGIONS)].sum(axis=1)
    return out[list(ALL_REGIONS)]


def breaks_for(region: str) -> tuple[int, ...]:
    return (params.EEU_BREAK,) if region in ("EEU", "ALL") else ()


def anomalies(series: pd.DataFrame, first: int, last: int, k: int = params.SSA_COMPONENTS,
              method: str = "ssa") -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    """Relative anomaly y / trend - 1 per region, with the trend fitted on first..last only."""
    anom, trends = {}, {}
    for region in series.columns:
        s = series.loc[first:last, region]
        if method == "ssa":
            t = ssa.detrend(s, params.SSA_WINDOW, k, breaks_for(region))
        else:
            t = ssa.moving_average_trend(s, 11, breaks_for(region))
        trends[region] = t
        anom[region] = s / t - 1.0
    return pd.DataFrame(anom), trends


def year_audit(anom: pd.DataFrame, counts: pd.DataFrame, path_tag: str, minimum: int = 2) -> pd.DataFrame:
    years = anom.index
    out = pd.DataFrame({"year": years})
    for k in params.WAVENUMBERS:
        out[f"wave{k}_event_weeks"] = counts.loc[years, k].to_numpy()
        out[f"wave{k}_class"] = counts.loc[years, k].map(lambda n: waves.year_class(int(n), minimum)).to_numpy()
    for region in anom.columns:
        out[f"{region}_anomaly_pct"] = 100 * anom[region].to_numpy()
    write_csv(out, RESULTS / f"year_audit_{path_tag}.csv", ["year"])
    return out


# --- S5: LMF tables and the four primary tests ---------------------------------------------------

def lmf_tables(anom: pd.DataFrame, counts: pd.DataFrame, years, minimum: int = 2) -> pd.DataFrame:
    parts = []
    for k in params.WAVENUMBERS:
        classes = lmf.classes_from_counts(counts.loc[years, k], minimum)
        t = lmf.lmf_table(anom.loc[years, list(params.REGIONS)], classes)
        t.insert(0, "wavenumber", k)
        parts.append(t)
    out = pd.concat(parts, ignore_index=True)
    out["q_bh"] = lmf.benjamini_hochberg(out.p_exact.tolist())
    return out


def lmf_step(manifest: Manifest, anom: pd.DataFrame, counts: pd.DataFrame, years, tag: str,
             n_perm: int = params.N_PERMUTATIONS, minimum: int = 2, record: bool = True) -> dict:
    years = np.asarray(list(years))
    table = lmf_tables(anom, counts, years, minimum)
    write_csv(table, RESULTS / f"lmf_{tag}.csv", ["wavenumber", "direction", "pair"])
    primary = lmf.permutation_test(anom.loc[years, list(params.REGIONS)], counts.loc[years], years,
                                   n_perm=n_perm, rng=params.rng_for(f"permutation-{tag}"),
                                   minimum=minimum)
    write_csv(primary, RESULTS / f"primary_tests_{tag}.csv", ["test"])
    if record:
        src = f"results/lmf_{tag}.csv"
        cell = table[(table.wavenumber == 7) & (table.direction == "low") & (table.pair == "NA x WEU")].iloc[0]
        for field in ("hits_multi", "n_multi", "hits_zero", "n_zero"):
            manifest.add(f"{tag}_w7_NA_WEU_{field}", int(cell[field]), src)
        for field in ("lmf", "ci_low", "ci_high", "p_exact"):
            manifest.add(f"{tag}_w7_NA_WEU_{field}", round(float(cell[field]), 3), src)
        psrc = f"results/primary_tests_{tag}.csv"
        for _, r in primary.iterrows():
            manifest.add(f"{tag}_{r.test}_pairs", int(r.observed_pairs), psrc, "of 10 pairs")
            manifest.add(f"{tag}_{r.test}_p_holm", round(float(r.p_holm), 3), psrc)
        manifest.add(f"{tag}_cells_q_below_0_1", int((table.q_bh < 0.1).sum()), src, "of 40 cells")
    return {"table": table, "primary": primary}


# --- S6: composites -----------------------------------------------------------------------------

def composite_step(manifest: Manifest, series: pd.DataFrame, counts: pd.DataFrame, first: int,
                   last: int, tag: str, record: bool = True) -> pd.DataFrame:
    rows = []
    for region in series.columns:
        s = series.loc[first:last, region]
        std = composites.neighbour_ratio(s, params.NEIGHBOUR_YEARS, breaks_for(region))
        for k in params.WAVENUMBERS:
            classes = lmf.classes_from_counts(counts.loc[std.index, k])
            c = composites.composite(std, classes)
            spread = composites.subsample_spread(std, classes, params.rng_for(f"subsample-{tag}-{region}-{k}"))
            boot = composites.bootstrap_interval(std, classes, params.rng_for(f"bootstrap-{tag}-{region}-{k}"))
            rows.append({"region": region, "wavenumber": k, **c, **spread, **boot,
                         "gap": c["multi_mean"] - c["zero_mean"]})
    out = pd.DataFrame(rows)
    write_csv(out, RESULTS / f"composites_{tag}.csv", ["wavenumber", "region"])
    if record:
        for _, r in out.iterrows():
            manifest.add(f"{tag}_comp_w{r.wavenumber}_{r.region}", round(float(r.multi_mean), 2),
                         f"results/composites_{tag}.csv", "% vs neighbouring years, multi-event years")
    return out


# --- S5-S7 and sensitivities, run after checkpoint 1 -------------------------------------------

def fao_extract(fao: pd.DataFrame, regions: pd.DataFrame) -> None:
    """Filtered FAOSTAT rows used by the study (CC BY 4.0, FAO). Committed only on Nikhil's go."""
    codes = set(regions.code)
    ext = fao[fao.code.isin(codes)][["code", "area", "crop", "element", "year", "value", "flag"]]
    write_csv(ext, DATA_DIR / "fao_extract.csv", ["code", "crop", "element", "year"], float_format="%.15g")


def _cell(table: pd.DataFrame, k: int, direction: str, pair: str) -> pd.Series:
    return table[(table.wavenumber == k) & (table.direction == direction) & (table.pair == pair)].iloc[0]


def after_cp1(manifest: Manifest, w: dict, y: dict, n_perm: int = params.N_PERMUTATIONS) -> dict:
    counts = w["counts"]
    series = region_series(y["table"])
    out = {}

    # Reproduction period 1961-2014 and full period 1961-2024.
    for tag, last in (("1961_2014", params.YIELD_BASE[1]), ("1961_2024", params.LAST_YEAR)):
        anom, trends = anomalies(series, params.YIELD_BASE[0], last)
        year_audit(anom, counts, tag)
        res = lmf_step(manifest, anom, counts, anom.index, tag, n_perm)
        comp = composite_step(manifest, series, counts, params.YIELD_BASE[0], last, tag)
        out[tag] = {"anom": anom, "trends": trends, **res, "composites": comp}
    figures.fig_regional_yields(y["table"], out["1961_2014"]["trends"], params.REGIONS)
    figures.fig_lmf(out["1961_2014"]["table"], "fig3_lmf_1961_2014.png", "1961-2014")
    figures.fig_lmf(out["1961_2024"]["table"], "fig3b_lmf_1961_2024.png", "1961-2024")
    figures.fig_composites(out["1961_2014"]["composites"], "fig4_composites_1961_2014.png")

    # Extension: the ten later years, and earlier years whose sign flips when the trend is refitted.
    full, base = out["1961_2024"]["anom"], out["1961_2014"]["anom"]
    ext = year_audit(full.loc[params.YIELD_BASE[1] + 1:], counts, "2015_2024")
    flips = []
    for region in params.REGIONS:
        b, f = base[region], full.loc[base.index, region]
        for yr in b.index[(np.sign(b) != np.sign(f))]:
            flips.append({"region": region, "year": int(yr), "anomaly_1961_2014_pct": 100 * b[yr],
                          "anomaly_1961_2024_pct": 100 * f[yr]})
    write_csv(pd.DataFrame(flips, columns=["region", "year", "anomaly_1961_2014_pct", "anomaly_1961_2024_pct"]),
              RESULTS / "sign_flips_after_refit.csv", ["region", "year"])
    manifest.add("sign_flips_after_refit", len(flips), "results/sign_flips_after_refit.csv",
                 "region-years before 2015 whose anomaly sign changes when the trend includes 2015-2024")
    for k in params.WAVENUMBERS:
        cls = ext[f"wave{k}_class"]
        for c in ("multi", "zero", "one"):
            manifest.add(f"ext_years_{c}_w{k}", int((cls == c).sum()), "results/year_audit_2015_2024.csv")

    # Sensitivities on 1961-2014: the four primary statistics and the worked-example cell.
    variants = {
        "primary": dict(),
        "ssa_k1": dict(k=1),
        "moving_average": dict(method="ma"),
        "maize_only": dict(table="maize_only"),
        "wheat_only": dict(table="wheat_only"),
        "china_351": dict(table="china_351"),
        "eeu_from_1992": dict(eeu_from=1992),
        "one_or_more_events": dict(minimum=1),
    }
    rows = []
    for name, opt in variants.items():
        tab = y["variants"][opt["table"]] if "table" in opt else y["table"]
        ser = region_series(tab)
        anom, _ = anomalies(ser, params.YIELD_BASE[0], params.YIELD_BASE[1], k=opt.get("k", params.SSA_COMPONENTS),
                            method=opt.get("method", "ssa"))
        if "eeu_from" in opt:
            anom.loc[:opt["eeu_from"] - 1, "EEU"] = np.nan
        minimum = opt.get("minimum", 2)
        years = np.asarray(anom.index)
        table = lmf_tables(anom, counts, years, minimum)
        prim = lmf.permutation_test(anom.loc[years, list(params.REGIONS)], counts.loc[years], years,
                                    n_perm=n_perm, rng=params.rng_for(f"permutation-sens-{name}"), minimum=minimum)
        cell = _cell(table, 7, "low", "NA x WEU")
        row = {"variant": name, "w7_NA_WEU_lmf": cell.lmf, "w7_NA_WEU_counts":
               f"{int(cell.hits_multi)}/{int(cell.n_multi)} vs {int(cell.hits_zero)}/{int(cell.n_zero)}"}
        for _, r in prim.iterrows():
            row[f"{r.test}_pairs"] = int(r.observed_pairs)
            row[f"{r.test}_p_holm"] = r.p_holm
        rows.append(row)
    sens = pd.DataFrame(rows)
    write_csv(sens, RESULTS / "sensitivities_1961_2014.csv", ["variant"])

    # Extension without Ukraine 2022-2024 (war years).
    tab = yields.region_table(y["fao"], y["regions"], range(params.YIELD_BASE[0], params.LAST_YEAR + 1),
                              overrides={"EEU": [(c, a, (2021 if c == 230 else b)) for c, a, b in
                                                 y["regions"][y["regions"].region == "EEU"][["code", "first_year", "last_year"]].itertuples(index=False)]})
    anom, _ = anomalies(region_series(tab), params.YIELD_BASE[0], params.LAST_YEAR)
    t = lmf_tables(anom, counts, anom.index)
    write_csv(t, RESULTS / "lmf_1961_2024_without_ukraine_2022_2024.csv", ["wavenumber", "direction", "pair"])
    out["sensitivities"] = sens
    return out


# --- S8: surface maps (CRU TS) ------------------------------------------------------------------

def surface_step(manifest: Manifest, events: pd.DataFrame) -> None:
    from . import surface
    maps, n_months = {}, {}
    for var in fetch.CRU_VARS:
        values, stn = surface.load(var)
        anom = surface.detrend_linear(values)
        for k in params.WAVENUMBERS:
            months = surface.event_months(events, k)
            n_months[k] = len(months)
            maps[(k, var)] = surface.composite_maps(anom, stn, months)
    for k in params.WAVENUMBERS:
        manifest.add(f"surface_event_months_w{k}", n_months[k], "results/wave_weeks.csv",
                     "distinct JJA months containing event weeks, 1960-2014")
    figures.fig_surface(maps, n_months, fetch.ne_path())
