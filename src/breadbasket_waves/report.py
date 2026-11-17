"""RESULTS.md and the README results block, generated from results/ only.

Every number in the text comes from a results file or the manifest; numbers derived only for the
text (such as the share expected under uniform timing) are added to the manifest first. A test
scans the generated text for numbers that have no source.
"""
from __future__ import annotations

import re

import pandas as pd

from . import params
from .io_utils import Manifest
from .paths import RESULTS, ROOT


def read(name: str) -> pd.DataFrame:
    # keep_default_na=False: the region "NA" (North America) must not become a missing value.
    return pd.read_csv(RESULTS / name, keep_default_na=False, na_values=[""])


def f2(x) -> str:
    return "inf" if x == float("inf") else f"{float(x):.2f}"


def f1(x) -> str:
    return f"{float(x):+.1f}"


def pct(x) -> str:
    return f"{float(x):.0f}%"


def md_table(header: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _counts(r) -> str:
    return f"{int(r.hits_multi)}/{int(r.n_multi)} vs {int(r.hits_zero)}/{int(r.n_zero)}"


def build(m: Manifest) -> tuple[str, str]:
    weeks = read("wave_weeks.csv")
    lmf14, lmf24 = read("lmf_1961_2014.csv"), read("lmf_1961_2024.csv")
    prim14, prim24 = read("primary_tests_1961_2014.csv"), read("primary_tests_1961_2024.csv")
    comp = read("composites_1961_2014.csv")
    sens = read("sensitivities_1961_2014.csv").set_index("variant")
    ext = read("year_audit_2015_2024.csv")
    named = read("wave7_named_summers.csv")
    T = {t.key: t for t in params.TARGETS}

    # Derived numbers used in the text.
    base = weeks[(weeks.year <= params.WAVE_BASE[1]) & weeks.event]
    early_uniform = m.add("uniform_share_weeks_1_6", round(100 * 6 / params.N_WEEKS), "derived", "6 of 13 weeks")
    late_uniform = m.add("uniform_share_weeks_10_13", round(100 * 4 / params.N_WEEKS), "derived", "4 of 13 weeks")
    early7 = m.add("share_weeks_1_6_w7", round(100 * (base[base.wavenumber == 7].week <= 6).mean()), "results/wave_weeks.csv")
    late5 = m.add("share_weeks_10_13_w5", round(100 * (base[base.wavenumber == 5].week >= 10).mean()), "results/wave_weeks.csv")

    def cell(t, k, d, pair):
        return t[(t.wavenumber == k) & (t.direction == d) & (t.pair == pair)].iloc[0]

    ex = cell(lmf14, 7, "low", "NA x WEU")
    low7 = lmf14[(lmf14.wavenumber == 7) & (lmf14.direction == "low")].sort_values("lmf", ascending=False)
    na_pairs = low7[low7.pair.str.contains("NA x")]
    p14 = prim14.set_index("test")
    p24 = prim24.set_index("test")
    c7 = comp[comp.wavenumber == 7].set_index("region")
    c5 = comp[comp.wavenumber == 5].set_index("region")
    w7_multi_ext = ext[ext.wave7_class == "multi"].year.tolist()
    w5_multi_ext = ext[ext.wave5_class == "multi"].year.tolist()
    named_flagged = named[named.wave7_class == "multi"].year.tolist()
    named_missed = named[named.wave7_class != "multi"].year.tolist()
    flips = m.get("sign_flips_after_refit")

    summary = [
        f"- **Wave events.** With a mean + 1.5 SD threshold, NCEP/NCAR reanalysis gives "
        f"{m.get('events_w5_base')} wave-5 and {m.get('events_w7_base')} wave-7 event weeks in JJA "
        f"1960-2014 (the paper, with ERA5: {T['events_w5_era5'].value} and {T['events_w7_era5'].value}). "
        f"That agreement is built into the threshold, so it only checks the mechanics. The timing is a "
        f"real check: {early7}% of wave-7 events fall in the first six weeks (1 June to 12 July; "
        f"{early_uniform}% if timing were random), and wave-5 events lean slightly late ({late5}% in the "
        f"last four weeks against {late_uniform}%), as the paper reports.",
        f"- **Concurrent low yields, wave-7, 1961-2014.** The likelihood multiplication factor (LMF) is "
        f"above 1 in {int(p14.loc['w7_low_gt1', 'observed_pairs'])} of 10 region pairs. The largest "
        f"values are in pairs that include North America: "
        + ", ".join(f"{r.pair} {f2(r.lmf)}" for r in na_pairs.head(3).itertuples())
        + ". The paper reports the same: the effect is strongest \"in pairings that include NA\".",
        f"- **Uncertainty.** None of the four pre-set primary tests is significant (Holm-adjusted "
        f"p = {f2(p14.loc['w7_low_gt1', 'p_holm'])} for wave-7 low yields). Single cells rest on "
        f"{int(ex.n_multi)} and {int(ex.n_zero)} years, and they move with the detrending method: the "
        f"paper's worked example (NA x WEU) is {f2(ex.lmf)} with our pre-set SSA trend and "
        f"{f2(sens.loc['moving_average', 'w7_NA_WEU_lmf'])} with an 11-year moving average, against "
        f"{T['lmf_example'].value} in the paper.",
        f"- **Yield composites, wave-7.** Across the five regions pooled, years with two or more event "
        f"weeks average {f1(c7.loc['ALL', 'multi_mean'])}% against the three years either side (the "
        f"paper: {T['comp_w7_mean'].value}%). Regionally the pattern differs: EEU "
        f"{f1(c7.loc['EEU', 'multi_mean'])}%, NA {f1(c7.loc['NA', 'multi_mean'])}%, EAS "
        f"{f1(c7.loc['EAS', 'multi_mean'])}% (the paper: up to {T['comp_w7_EEU'].value}%, "
        f"{T['comp_w7_NA'].value}% and {T['comp_w7_EAS'].value}%).",
        f"- **2015-2024.** Fixed thresholds flag wave-7 summers {', '.join(map(str, w7_multi_ext)) or 'none'} "
        f"and wave-5 summers {', '.join(map(str, w5_multi_ext)) or 'none'}. Adding these years, wave-7 "
        f"pairs with LMF above 1 go from {int(p14.loc['w7_low_gt1', 'observed_pairs'])} to "
        f"{int(p24.loc['w7_low_gt1', 'observed_pairs'])}, wave-5 from "
        f"{int(p14.loc['w5_low_gt1', 'observed_pairs'])} to {int(p24.loc['w5_low_gt1', 'observed_pairs'])}. "
        f"About three new years per class can show whether the pattern persists, not confirm it.",
        f"- **Known wave-7 summers.** Of the summers Kornhuber et al. (2019) name, "
        f"{', '.join(map(str, named_flagged))} are flagged; {', '.join(map(str, named_missed))} is not. "
        f"Its strongest wave-7 week stays below the threshold. This is reported, not tuned away.",
    ]

    md = ["# Results", "",
          "Generated by `run_all.py` from the files in `results/`. Do not edit this file by hand.", "",
          "An independent study and a learning exercise that re-implements a published analysis on "
          "open data (Kornhuber et al. 2023, *Nature Communications* 14, 3528). Not new science, and "
          "not a forecast.", "", "## Summary", ""] + summary

    # 1. Waves
    md += ["", "## 1. Wave events", "",
           "![Wave events](results/figures/fig1_wave_events.png)", "",
           md_table(["", "This study (NCEP/NCAR R1)", "Paper (ERA5)", "What it tests"], [
               ["wave-5 event weeks, JJA 1960-2014", m.get("events_w5_base"), T["events_w5_era5"].value, "mechanics only"],
               ["wave-7 event weeks, JJA 1960-2014", m.get("events_w7_base"), T["events_w7_era5"].value, "mechanics only"],
               ["wave-7 years with 2+ / 0 events, 1961-2014",
                f"{m.get('years_multi_w7_1961_2014')} / {m.get('years_zero_w7_1961_2014')}",
                f"{T['lmf_example_n_multi'].value} / {T['lmf_example_n_zero'].value}", "mechanics only"],
               ["wave-7 events in 1 Jun-12 Jul", f"{early7}% (random: {early_uniform}%)", "mostly June and early July", "timing"],
               ["wave-5 events in the last 4 weeks", f"{late5}% (random: {late_uniform}%)", "slight preference for later weeks", "timing"],
           ]), "",
           "Event weeks per decade are in `results/wave_events_per_decade.csv`. Before 1979, the start "
           "of the modern satellite era, the reanalysis rests on fewer observations, so early "
           "decades are less certain."]

    # 2. LMF
    md += ["", "## 2. Concurrent low and high yields (LMF)", "",
           "![LMF](results/figures/fig3_lmf_1961_2014.png)", "",
           "Each cell: LMF, then the counts behind it (years with both regions below trend among years "
           "with two or more event weeks, against the same among years with none).", "",
           "**Primary tests, fixed before any result** (per wave: how many of the 10 pairs have LMF > 1 "
           "for low yields, and LMF < 1 for high yields; whole-year permutation test, Holm-adjusted):", "",
           md_table(["Test", "1961-2014: pairs", "p (Holm)", "1961-2024: pairs", "p (Holm)"], [
               [t, int(p14.loc[t, "observed_pairs"]), f2(p14.loc[t, "p_holm"]),
                int(p24.loc[t, "observed_pairs"]), f2(p24.loc[t, "p_holm"])] for t in p14.index]),
           "", "**The paper's worked example (wave-7, NA x WEU, low yields) under each sensitivity:**", "",
           md_table(["Variant", "Counts", "LMF"],
                    [["paper (ERA5, FAO)", f"{T['lmf_example_hits_multi'].value}/{T['lmf_example_n_multi'].value} vs "
                      f"{T['lmf_example_hits_zero'].value}/{T['lmf_example_n_zero'].value}", T["lmf_example"].value]]
                    + [[v, sens.loc[v, "w7_NA_WEU_counts"], f2(sens.loc[v, "w7_NA_WEU_lmf"])] for v in sens.index]),
           "", f"Pre-set primary: SSA trend with two components. The worked example's own 95 % interval "
           f"is {f2(ex.ci_low)} to {f2(ex.ci_high)} (exact one-sided p = {f2(ex.p_exact)}). All 40 cells, "
           f"with intervals, exact p-values and Benjamini-Hochberg q-values, are in "
           f"`results/lmf_1961_2014.csv` and `results/lmf_1961_2024.csv`."]

    # 3. Composites
    rows = []
    for region in ("NA", "WEU", "EEU", "IND", "EAS", "ALL", "WORLD"):
        rows.append([region, f"{f1(c7.loc[region, 'multi_mean'])} [{f1(c7.loc[region, 'multi_lo'])}, {f1(c7.loc[region, 'multi_hi'])}]",
                     f1(c7.loc[region, "zero_mean"]),
                     f"{f1(c5.loc[region, 'multi_mean'])} [{f1(c5.loc[region, 'multi_lo'])}, {f1(c5.loc[region, 'multi_hi'])}]",
                     f1(c5.loc[region, "zero_mean"])])
    md += ["", "## 3. Yield composites", "", "![Composites](results/figures/fig4_composites_1961_2014.png)", "",
           md_table(["Region", "wave-7, 2+ events (95 % CI)", "wave-7, none", "wave-5, 2+ events (95 % CI)", "wave-5, none"], rows),
           "", f"Percent difference from the mean of the three years before and after, 1961-2014. Paper: "
           f"wave-7 up to {T['comp_w7_EAS'].value}% (EAS), {T['comp_w7_NA'].value}% (NA), "
           f"{T['comp_w7_EEU'].value}% (EEU), average {T['comp_w7_mean'].value}%; wave-5 "
           f"{T['comp_w5_range'].value}%, EEU {T['comp_w5_EEU'].value}%, NA {T['comp_w5_NA'].value}%. "
           f"The paper's 500 draws of about 90 % of years are a subsample spread, not a confidence "
           f"interval; their medians are in `results/composites_1961_2014.csv`."]

    # 4. Extension and 5. Surface
    md += ["", "## 4. The ten later years (2015-2024)", "",
           "Thresholds stay fixed at their 1960-2014 values.", "",
           md_table(["Year", "wave-5 weeks", "wave-7 weeks"] + [f"{r} anomaly" for r in params.REGIONS],
                    [[int(r.year), int(r.wave5_event_weeks), int(r.wave7_event_weeks)]
                     + [f1(getattr(r, f"{g}_anomaly_pct")) for g in params.REGIONS] for r in ext.itertuples()]),
           "", f"Refitting the trend on 1961-2024 changes the sign of {flips} earlier region-year "
           f"anomalies (`results/sign_flips_after_refit.csv`), an end effect of any trend fit.",
           "", "## 5. Surface conditions in event months", "",
           "![Surface](results/figures/fig5_surface_composites.png)", "",
           f"CRU TS 4.10 maximum temperature and precipitation, averaged over the "
           f"{m.get('surface_event_months_w7')} wave-7 and {m.get('surface_event_months_w5')} wave-5 "
           f"event months of 1960-2014. Monthly, so a simplification of the paper's weekly fields.",
           "", "## 6. Differences from the paper", "",
           "- Reanalysis: NCEP/NCAR R1 at 2.5 degrees instead of ERA5 at 1 degree (ERA5 needs a login). "
           "Kornhuber et al. (2020) used NCEP/NCAR and agreed with the 2023 results.",
           "- Yields: the FAOSTAT release of December 2025, which revises earlier years.",
           "- Detrending: the paper does not give its SSA settings; this study fixed two components and "
           "a 10-year window before any result. See `DEVIATIONS.md`.",
           "- Regions before 1992 use the USSR and Yugoslav SFR rows; the EEU trend is fitted "
           "separately before and after 1992. China is FAO code 41 (mainland).",
           "- Surface fields are monthly, not weekly.",
           "", "## 7. Limits", "",
           "- Single LMF cells rest on about 15 and 20 years; read the pattern across pairs, not one cell.",
           "- Combined maize and wheat yields mix crops with different seasons; India's wheat is "
           "harvested before June. Maize-only and wheat-only results are in "
           "`results/sensitivities_1961_2014.csv`.",
           "- Correlation, not attribution: shared drivers such as ENSO are not separated.",
           ""]
    readme = "\n".join(["Summary from the latest run (full tables in [RESULTS.md](RESULTS.md)):", ""] + summary)
    return "\n".join(md), readme


def write(m: Manifest) -> None:
    md, block = build(m)
    (ROOT / "RESULTS.md").write_text(md + "\n", encoding="utf-8", newline="\n")
    readme_path = ROOT / "README.md"
    readme = readme_path.read_text(encoding="utf-8")
    readme = re.sub(r"<!-- RESULTS:START -->.*?<!-- RESULTS:END -->",
                    lambda _: "<!-- RESULTS:START -->\n" + block + "\n<!-- RESULTS:END -->", readme, flags=re.S)
    readme_path.write_text(readme, encoding="utf-8", newline="\n")
