# breadbasket-waves

An independent study by Nikhil Sehgal. It is a learning exercise that re-implements a published
analysis on open data: the observed link between amplified summer wave patterns in the jet stream
and low maize and wheat yields in several breadbasket regions in the same year, reported by
Kornhuber et al. (2023). It then adds the ten years after the paper's period ends (2015-2024).

It is not new science, and it is not a forecast.

## The question

Using only open data, does the observed rise in concurrent low maize-and-wheat yields in summers
with amplified wave-5 and wave-7 patterns reproduce, and does it hold when the years after 2014
are added?

## Results

<!-- RESULTS:START -->
Summary from the latest run (full tables in [RESULTS.md](RESULTS.md)):

- **Wave events.** With a mean + 1.5 SD threshold, NCEP/NCAR reanalysis gives 57 wave-5 and 56 wave-7 event weeks in JJA 1960-2014 (the paper, with ERA5: 62 and 54). That agreement is built into the threshold, so it only checks the mechanics. The timing is a real check: 64% of wave-7 events fall in the first six weeks (1 June to 12 July; 46% if timing were random), and wave-5 events lean slightly late (35% in the last four weeks against 31%), as the paper reports.
- **Concurrent low yields, wave-7, 1961-2014.** The likelihood multiplication factor (LMF) is above 1 in 7 of 10 region pairs. The largest values are in pairs that include North America: NA x WEU 5.60, NA x EEU 5.60, NA x IND 3.50. The paper reports the same: the effect is strongest "in pairings that include NA".
- **Uncertainty.** None of the four pre-set primary tests is significant (Holm-adjusted p = 1.00 for wave-7 low yields). Single cells rest on 15 and 21 years, and they move with the detrending method: the paper's worked example (NA x WEU) is 5.60 with our pre-set SSA trend and 1.68 with an 11-year moving average, against 1.6 in the paper.
- **Yield composites, wave-7.** Across the five regions pooled, years with two or more event weeks average -2.6% against the three years either side (the paper: -2 to -3%). Regionally the pattern differs: EEU -6.2%, NA -3.3%, EAS +1.0% (the paper: up to -3%, -6% and -7%).
- **2015-2024.** Fixed thresholds flag wave-7 summers 2015, 2022 and wave-5 summers 2021, 2022. Adding these years, wave-7 pairs with LMF above 1 go from 7 to 6, wave-5 from 7 to 9. About three new years per class can show whether the pattern persists, not confirm it.
- **Known wave-7 summers.** Of the summers Kornhuber et al. (2019) name, 2003, 2006, 2015 are flagged; 2018 is not. Its strongest wave-7 week stays below the threshold. This is reported, not tuned away.
<!-- RESULTS:END -->

## How to run

```
uv sync
uv run pytest
uv run python run_all.py
```

`run_all.py` downloads the data into `local/` (about 600 MB, checked against the hashes in
`data/sources.csv`), then rebuilds every table, figure, `RESULTS.md` and `site/index.html`.

## Data and credits

- NCEP/NCAR Reanalysis 1 daily meridional wind: data provided by the NOAA PSL, Boulder, Colorado,
  USA, from their website at https://psl.noaa.gov (Kalnay et al. 1996, *Bulletin of the American
  Meteorological Society* 77, 437-470).
- FAO. 2025. FAOSTAT: Crops and livestock products. Accessed on 17 November 2026.
  https://www.fao.org/faostat/en/#data/QCL. Licence: CC-BY-4.0. `data/fao_extract.csv` is a
  filtered extract (the study's countries, maize and wheat, production and harvested area,
  1961-2024), kept because FAO replaces its bulk file in place. FAO does not endorse this study.
- CRU TS 4.10 (Harris et al. 2020, *Scientific Data* 7, 109), Open Database License. Only figures
  derived from it are published here.
- Made with Natural Earth (public domain).

Full source records, with download dates and SHA-256 hashes, are in `data/sources.csv`. Method
choices and every change to them are in `DEVIATIONS.md`.
