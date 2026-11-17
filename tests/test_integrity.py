"""Integrity checks on the published text: every number has a source, and the text is clean."""
import json
import re
from pathlib import Path

import pandas as pd
import pytest

from breadbasket_waves import params
from breadbasket_waves.paths import RESULTS, ROOT

HAVE_RESULTS = (RESULTS / "manifest.json").exists() and (ROOT / "RESULTS.md").exists()
needs_results = pytest.mark.skipif(not HAVE_RESULTS, reason="run run_all.py first")

TEXT_FILES = ["README.md", "RESULTS.md", "DEVIATIONS.md", "LICENSE"]


def _numbers(text: str) -> set[str]:
    text = re.sub(r"https?://\S+|doi\S*|`[^`]*`|\([^)]*\.png\)", " ", text)      # links, code, image paths
    text = re.sub(r"\*Nature Communications\* 14, 3528|BAMS 77, 437-470|Kalnay et al\. 1996|"
                  r"Harris et al\. 2020|CRU TS 4\.10", " ", text)
    return {n.lstrip("+-").rstrip("%") for n in re.findall(r"[+-]?\d+(?:\.\d+)?%?", text)}


def _known() -> set[str]:
    known: set[str] = set()

    def add(v):
        try:
            x = float(v)
        except (TypeError, ValueError):
            known.update(n.lstrip("+-") for n in re.findall(r"[+-]?\d+(?:\.\d+)?", str(v)))
            return
        for s in (f"{x:g}", f"{x:.0f}", f"{x:.1f}", f"{x:.2f}", f"{abs(x):.1f}", f"{abs(x):.2f}", f"{abs(x):g}"):
            known.add(s.lstrip("+-"))

    for item in json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8")).values():
        add(item["value"])
    for t in params.TARGETS:
        add(t.value)
    for csv in RESULTS.glob("*.csv"):
        for v in pd.read_csv(csv, keep_default_na=False).to_numpy().ravel():
            add(v)
            for part in re.findall(r"\d+", str(v)):
                known.add(part)
    known.update(str(y) for y in range(1900, 2031))           # years
    known.update(str(i) for i in range(0, 14))                # weeks, figure and section numbers
    known.update({"1.5", "37.5", "57.5", "250", "2.5", "40", "500", "90", "95", "10", "20", "75"})
    return known


@needs_results
@pytest.mark.parametrize("name", ["README.md", "RESULTS.md"])
def test_numbers_in_text_come_from_results(name):
    text = (ROOT / name).read_text(encoding="utf-8")
    if name == "README.md":                                   # only the generated block carries results
        text = text.split("<!-- RESULTS:START -->")[1].split("<!-- RESULTS:END -->")[0]
    unknown = _numbers(text) - _known()
    assert not unknown, f"numbers in {name} with no source: {sorted(unknown)}"


def test_text_scan_clean():
    files = [ROOT / f for f in TEXT_FILES if (ROOT / f).exists()]
    files += list((ROOT / "src").rglob("*.py")) + list((ROOT / "tests").rglob("*.py"))
    for path in files:
        text = path.read_text(encoding="utf-8")
        assert chr(0x2014) not in text, f"em dash in {path.name}"
        assert "Nik" + "khil" not in text, f"old spelling in {path.name}"
        assert "Co-" + "Authored-By" not in text, f"co-author line in {path.name}"
