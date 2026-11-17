"""Rebuild everything: fetch and verify the data, then every table, figure and report.

    uv run python run_all.py              fetch what is missing, verify every hash, rebuild
    uv run python run_all.py --offline    no network; stops if a file is missing
    uv run python run_all.py --repin      accept changed sources (log the change in DEVIATIONS.md)
    uv run python run_all.py --until cp1  stop after the outputs for checkpoint 1
"""
from __future__ import annotations

import argparse
import sys
import warnings

from breadbasket_waves import fetch, pipeline, report, sitebuild
from breadbasket_waves.io_utils import Manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--repin", action="store_true")
    parser.add_argument("--no-surface", action="store_true", help="skip the CRU TS downloads and maps")
    parser.add_argument("--until", choices=["cp1", "all"], default="all")
    args = parser.parse_args(argv)
    warnings.filterwarnings("ignore", message=".*dap2.*")

    print("S2  fetch and verify sources")
    fetch.run(offline=args.offline, repin=args.repin, surface=not args.no_surface)
    manifest = Manifest()
    print("S3  wave events")
    w = pipeline.wave_step(manifest)
    print("S4a regional yields and coverage")
    y = pipeline.yield_step(manifest)
    print("CP1 spot checks")
    checks = pipeline.spot_checks(y["fao"], y["regions"], y["table"])
    manifest.add("cp1_spot_checks_matching", int(checks.match.sum()), "results/cp1_fao_spot_checks.csv",
                 f"of {len(checks)} seeded region-years recomputed with a plain loop")
    manifest.save()
    if args.until == "cp1":
        print("stopped after checkpoint 1 outputs")
        return 0
    pipeline.fao_extract(y["fao"], y["regions"])
    print("S5-S7 LMF, permutation tests, composites, extension, sensitivities")
    pipeline.after_cp1(manifest, w, y)
    if not args.no_surface:
        print("S8  surface maps (CRU TS)")
        pipeline.surface_step(manifest, w["events"])
    print("S9  RESULTS.md, README and site")
    report.write(manifest)
    manifest.save()
    sitebuild.build()
    return 0


if __name__ == "__main__":
    sys.exit(main())
