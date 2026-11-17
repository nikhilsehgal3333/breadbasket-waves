"""Regional maize and wheat yields from the FAOSTAT bulk file.

Regions follow Kornhuber et al. (2023), Methods: national production is summed and divided by the
summed harvested area in each year. Countries are matched on FAO area codes through
data/regions.csv, never on names, because names change between releases. Predecessor states
(USSR, Yugoslav SFR, Belgium-Luxembourg, Serbia and Montenegro) stand in for their successors
before the successors' first year.
"""
from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from .paths import DATA

ITEMS = {"15": "wheat", "56": "maize"}             # never 446 (green maize)
ELEMENTS = {"5510": ("production", "t"), "5312": ("area", "ha")}
MEMBER = "Production_Crops_Livestock_E_All_Data_(Normalized).csv"
ESTIMATED_FLAGS = {"E", "I"}                       # estimated or imputed by FAO
REGIONS_CSV = DATA / "regions.csv"
# (predecessor, successor) pairs: their year spans in regions.csv must never overlap.
SUCCESSIONS = ((15, 255), (228, 185), (228, 230), (228, 146), (248, 154), (248, 186), (248, 272),
               (186, 272))


class CoverageError(ValueError):
    pass


def load_regions(path: Path = REGIONS_CSV) -> pd.DataFrame:
    # keep_default_na=False: otherwise pandas reads the region "NA" (North America) as missing.
    return pd.read_csv(path, dtype={"region": str, "code": int, "first_year": int, "last_year": int},
                       keep_default_na=False)


def load_fao(zip_path: Path) -> pd.DataFrame:
    """Rows for wheat and maize production and harvested area, with units checked."""
    rows = []
    with zipfile.ZipFile(zip_path) as z, z.open(MEMBER) as f:
        reader = csv.reader(io.TextIOWrapper(f, encoding="utf-8", errors="strict"))
        header = next(reader)
        i = {h: n for n, h in enumerate(header)}
        for r in reader:
            item, element = r[i["Item Code"]], r[i["Element Code"]]
            if item not in ITEMS or element not in ELEMENTS:
                continue
            name, unit = ELEMENTS[element]
            if r[i["Unit"]] != unit:
                raise ValueError(f"element {element}: unit {r[i['Unit']]!r}, expected {unit!r}")
            value = r[i["Value"]]
            rows.append((int(r[i["Area Code"]]), r[i["Area Code (M49)"]].lstrip("'"), r[i["Area"]],
                         ITEMS[item], name, int(r[i["Year"]]),
                         float(value) if value != "" else np.nan, r[i["Flag"]]))
    df = pd.DataFrame(rows, columns=["code", "m49", "area", "crop", "element", "year", "value", "flag"])
    return df.sort_values(["code", "crop", "element", "year"]).reset_index(drop=True)


def active_codes(regions: pd.DataFrame, region: str, year: int) -> list[int]:
    sel = regions[(regions.region == region) & (regions.first_year <= year) & (regions.last_year >= year)]
    return sorted(sel.code.tolist())


def check_successions(regions: pd.DataFrame) -> None:
    """No year may count a predecessor together with one of its successors."""
    span = {r.code: (r.first_year, r.last_year, r["name"]) for _, r in regions.iterrows()}
    for pred, succ in SUCCESSIONS:
        if pred in span and succ in span:
            (a0, a1, an), (b0, b1, bn) = span[pred], span[succ]
            if not (a1 < b0 or b1 < a0):
                raise CoverageError(f"{an} ({a0}-{a1}) and {bn} ({b0}-{b1}) overlap")


def coverage(fao: pd.DataFrame, regions: pd.DataFrame, years: range) -> pd.DataFrame:
    """One row per region, year and listed code: does FAO report wheat production and area?

    Raises CoverageError when a listed country has no wheat production or area in a year it is
    listed for (wheat is grown by every listed country; maize is not)."""
    have = set(zip(fao.code, fao.crop, fao.element, fao.year))
    out, missing = [], []
    for region in sorted(set(regions.region) - {"WORLD"}):
        for year in years:
            for code in active_codes(regions, region, year):
                ok = all((code, "wheat", e, year) in have for e in ("production", "area"))
                maize = all((code, "maize", e, year) in have for e in ("production", "area"))
                out.append((region, year, code, ok, maize))
                if not ok:
                    missing.append(f"{region} {year} code {code}")
    if missing:
        raise CoverageError("missing wheat rows: " + ", ".join(missing[:20]))
    return pd.DataFrame(out, columns=["region", "year", "code", "wheat", "maize"])


def region_table(fao: pd.DataFrame, regions: pd.DataFrame, years: range,
                 crops: tuple[str, ...] = ("maize", "wheat"),
                 overrides: dict[str, list[tuple[int, int, int]]] | None = None) -> pd.DataFrame:
    """Production, area and yield per region and year for the chosen crops.

    `overrides` replaces a region's codes, e.g. {"EAS": [(351, 1961, 2024), (141, 1961, 2024)]}
    for the China 351 sensitivity."""
    sub = fao[fao.crop.isin(crops)]
    per_crop = sub.pivot_table(index=["code", "year", "crop"], columns="element", values="value",
                               aggfunc="first")
    # A crop counts in a country-year only when both production and area are reported.
    per_crop = per_crop[np.isfinite(per_crop.get("production")) & np.isfinite(per_crop.get("area"))
                        & (per_crop.get("area") > 0)]
    flags = sub[sub.element == "production"].set_index(["code", "year", "crop"]).flag
    per_crop = per_crop.join(flags.rename("flag"), how="left")
    per_crop["estimated"] = np.where(per_crop.flag.isin(ESTIMATED_FLAGS), per_crop.production, 0.0)
    totals = per_crop.groupby(level=["code", "year"])[["production", "area", "estimated"]].sum()
    rows = []
    for region in list(dict.fromkeys(regions.region)):
        for year in years:
            if overrides and region in overrides:
                codes = [c for c, a, b in overrides[region] if a <= year <= b]
            else:
                codes = active_codes(regions, region, year)
            prod = area = est = 0.0
            for code in codes:
                if (code, year) in totals.index:
                    t = totals.loc[(code, year)]
                    prod += float(t.production)
                    area += float(t.area)
                    est += float(t.estimated)
            rows.append((region, year, prod, area, prod / area if area > 0 else np.nan,
                         est / prod if prod > 0 else np.nan))
    return pd.DataFrame(rows, columns=["region", "year", "production_t", "area_ha",
                                       "yield_t_ha", "estimated_share"])
