"""Byte-stable outputs and the results manifest.

Every CSV is written with sorted rows, "\\n" line endings and a fixed float format, so two runs on
any machine give identical files. Every number that appears in README.md, RESULTS.md or the site
is first recorded in results/manifest.json with the file it comes from.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .paths import RESULTS


def write_csv(df: pd.DataFrame, path: Path, sort_by: list[str] | None = None,
              float_format: str = "%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if sort_by:
        df = df.sort_values(sort_by, kind="mergesort")
    df.to_csv(path, index=False, lineterminator="\n", float_format=float_format, encoding="utf-8")
    return path


def _plain(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


class Manifest:
    """Key numbers with their source file, in one sorted JSON file."""

    def __init__(self, path: Path = RESULTS / "manifest.json", load: bool = False):
        """A build starts from an empty manifest; load=True reads the saved one (for reports)."""
        self.path = path
        self.items: dict[str, dict[str, Any]] = {}
        if load and path.exists():
            self.items = json.loads(path.read_text(encoding="utf-8"))

    def add(self, key: str, value: Any, source: str, note: str = "") -> Any:
        self.items[key] = {"value": _plain(value), "source": source, "note": note}
        return value

    def get(self, key: str) -> Any:
        return self.items[key]["value"]

    def save(self) -> Path:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(self.items, indent=1, sort_keys=True, ensure_ascii=False) + "\n"
        self.path.write_text(text, encoding="utf-8", newline="\n")
        return self.path
