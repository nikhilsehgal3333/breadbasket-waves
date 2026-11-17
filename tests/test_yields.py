import csv
import io
import zipfile

import numpy as np
import pandas as pd
import pytest

from breadbasket_waves import yields

HEADER = ["Area Code", "Area Code (M49)", "Area", "Item Code", "Item Code (CPC)", "Item",
          "Element Code", "Element", "Year Code", "Year", "Unit", "Value", "Flag", "Note"]


def _row(code, area, item, element, year, value, flag="A", unit=None):
    item_name = {"15": "Wheat", "56": "Maize (corn)", "446": "Green corn (maize)"}[item]
    el_name, el_unit = {"5510": ("Production", "t"), "5312": ("Area harvested", "ha")}[element]
    return [str(code), f"'{code:03d}", area, item, "'0111", item_name, element, el_name,
            str(year), str(year), unit or el_unit, str(value), flag, ""]


def _zip(tmp_path, rows):
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(HEADER)
    w.writerows(rows)
    path = tmp_path / "fao.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr(yields.MEMBER, buf.getvalue())
    return path


def _regions(rows):
    return pd.DataFrame(rows, columns=["region", "code", "name", "first_year", "last_year", "role"])


def test_codes_not_names_and_green_maize_ignored(tmp_path):
    rows = [
        _row(41, "China, mainland", "15", "5510", 2000, 100), _row(41, "China, mainland", "15", "5312", 2000, 50),
        _row(351, "China", "15", "5510", 2000, 999), _row(351, "China", "15", "5312", 2000, 1),
        _row(41, "China, mainland", "446", "5510", 2000, 777),
        _row(5000, "World", "15", "5510", 2000, 5000), _row(5000, "World", "15", "5312", 2000, 2500),
    ]
    fao = yields.load_fao(_zip(tmp_path, rows))
    assert set(fao.crop) == {"wheat"}                    # item 446 dropped
    assert fao.m49.str.startswith("'").sum() == 0        # apostrophes stripped
    regions = _regions([("EAS", 41, "renamed in this test", 1961, 2024, "country"),
                        ("WORLD", 5000, "World", 1961, 2024, "aggregate")])
    table = yields.region_table(fao, regions, range(2000, 2001))
    eas = table[table.region == "EAS"].iloc[0]
    assert eas.production_t == 100 and eas.area_ha == 50   # 351 is not added to 41


def test_combined_yield_is_ratio_of_sums(tmp_path):
    rows = [_row(231, "USA", "15", "5510", 2000, 60), _row(231, "USA", "15", "5312", 2000, 20),
            _row(231, "USA", "56", "5510", 2000, 300), _row(231, "USA", "56", "5312", 2000, 30),
            _row(33, "Canada", "15", "5510", 2000, 40), _row(33, "Canada", "15", "5312", 2000, 20)]
    fao = yields.load_fao(_zip(tmp_path, rows))
    regions = _regions([("NA", 231, "USA", 1961, 2024, "country"), ("NA", 33, "Canada", 1961, 2024, "country")])
    t = yields.region_table(fao, regions, range(2000, 2001)).iloc[0]
    assert t.yield_t_ha == pytest.approx((60 + 300 + 40) / (20 + 30 + 20))
    wheat_only = yields.region_table(fao, regions, range(2000, 2001), crops=("wheat",)).iloc[0]
    assert wheat_only.yield_t_ha == pytest.approx((60 + 40) / (20 + 20))


def test_crop_without_area_is_left_out(tmp_path):
    rows = [_row(231, "USA", "15", "5510", 2000, 60), _row(231, "USA", "15", "5312", 2000, 20),
            _row(231, "USA", "56", "5510", 2000, 300)]          # maize area missing
    fao = yields.load_fao(_zip(tmp_path, rows))
    regions = _regions([("NA", 231, "USA", 1961, 2024, "country")])
    assert yields.region_table(fao, regions, range(2000, 2001)).iloc[0].yield_t_ha == pytest.approx(3.0)


def test_wrong_unit_raises(tmp_path):
    rows = [_row(231, "USA", "15", "5510", 2000, 60, unit="1000 t")]
    with pytest.raises(ValueError):
        yields.load_fao(_zip(tmp_path, rows))


def test_predecessor_then_successor(tmp_path):
    rows = []
    for year in (1990, 1991):
        rows += [_row(228, "USSR", "15", "5510", year, 90), _row(228, "USSR", "15", "5312", year, 60)]
    for year in (1992, 1993):
        rows += [_row(185, "Russian Federation", "15", "5510", year, 50), _row(185, "Russian Federation", "15", "5312", year, 25)]
    fao = yields.load_fao(_zip(tmp_path, rows))
    regions = _regions([("EEU", 228, "USSR", 1961, 1991, "predecessor"),
                        ("EEU", 185, "Russian Federation", 1992, 2024, "country")])
    t = yields.region_table(fao, regions, range(1990, 1994)).set_index("year")
    assert list(t.yield_t_ha.round(2)) == [1.5, 1.5, 2.0, 2.0]


def test_region_map_has_no_double_counting():
    regions = yields.load_regions()
    yields.check_successions(regions)
    bad = pd.concat([regions, _regions([("EEU", 185, "Russian Federation", 1990, 2024, "country")])])
    with pytest.raises(yields.CoverageError):
        yields.check_successions(bad.drop_duplicates(subset=["code"], keep="last"))


def test_region_map_matches_the_plan():
    regions = yields.load_regions()
    assert set(regions.region) == {"NA", "WEU", "EEU", "IND", "EAS", "WORLD"}
    assert 351 not in set(regions.code)                  # China = 41 (mainland); 351 is a sensitivity
    belgium = regions.set_index("code")
    assert belgium.loc[15, "last_year"] == 1999 and belgium.loc[255, "first_year"] == 2000
    assert belgium.loc[79, "first_year"] == 1961        # Germany continuous from 1961


def test_missing_country_raises(tmp_path):
    rows = [_row(231, "USA", "15", "5510", 2000, 60), _row(231, "USA", "15", "5312", 2000, 20)]
    fao = yields.load_fao(_zip(tmp_path, rows))
    regions = _regions([("NA", 231, "USA", 1961, 2024, "country"), ("NA", 33, "Canada", 1961, 2024, "country")])
    with pytest.raises(yields.CoverageError):
        yields.coverage(fao, regions, range(2000, 2001))
