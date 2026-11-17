"""Download every source into local/ and verify it against the lock file data/sources.csv.

Trust on first use: the first download of a source writes its SHA-256 into the lock file. Every
later run recomputes the hash and stops on a mismatch. Re-pinning a changed source needs
`run_all.py --repin` and an entry in DEVIATIONS.md.

Only the reanalysis is subset on the server: one OPeNDAP request per year returns the June-August
250 hPa meridional wind between 37.5 and 57.5 degrees north, saved as a NumPy array (.npy files
are byte-stable; netCDF/HDF5 files are not).
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import io
import time
from pathlib import Path

import numpy as np

from . import net, params
from .paths import DATA, LOCAL

LOCK = DATA / "sources.csv"
FIELDS = ("key", "url", "retrieved", "bytes", "sha256", "licence", "citation")

NCEP_URL = "dap2://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/Dailies/pressure/vwnd.{year}.nc"
FAO_URL = ("https://bulks-faostat.fao.org/production/"
           "Production_Crops_Livestock_E_All_Data_(Normalized).zip")
CRU_URL = ("https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.10/cruts.2604091129.v4.10/"
           "{var}/cru_ts4.10.{a}.{b}.{var}.dat.nc.gz")
CRU_DECADES = ((1951, 1960), (1961, 1970), (1971, 1980), (1981, 1990), (1991, 2000),
               (2001, 2010), (2011, 2020), (2021, 2025))
CRU_VARS = ("tmx", "pre")
NE_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/v5.1.2/geojson/"
          "ne_110m_admin_0_countries.geojson")

# The NCEP/NCAR R1 grid inside the latitude band, north to south, and all longitudes.
NCEP_LATS = np.arange(57.5, 37.4, -2.5)
NCEP_LONS = np.arange(0.0, 360.0, 2.5)

LICENCE = {
    "ncep": "No usage restrictions; acknowledge the NOAA PSL",
    "fao": "CC BY 4.0 with the FAO Statistical Database Terms of Use",
    "cru": "Open Database License 1.0 with the Database Contents License",
    "ne": "Public domain",
}
CITATION = {
    "ncep": ("Kalnay et al. (1996), Bulletin of the American Meteorological Society 77, 437-470. "
             "NCEP Reanalysis data provided by the NOAA PSL, Boulder, Colorado, USA"),
    "fao": "FAO. FAOSTAT: Crops and livestock products (QCL). https://www.fao.org/faostat/en/#data/QCL",
    "cru": ("Harris, I., Osborn, T. J., Jones, P. & Lister, D. (2020). Version 4 of the CRU TS "
            "monthly high-resolution gridded multivariate climate dataset. Scientific Data 7, 109"),
    "ne": "Natural Earth, 1:110m Admin 0 Countries, release v5.1.2. https://www.naturalearthdata.com",
}


class LockMismatch(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_lock() -> dict[str, dict[str, str]]:
    if not LOCK.exists():
        return {}
    with open(LOCK, encoding="utf-8", newline="") as f:
        return {row["key"]: row for row in csv.DictReader(f)}


def write_lock(rows: dict[str, dict[str, str]]) -> None:
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with open(LOCK, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        for key in sorted(rows):
            w.writerow({k: rows[key].get(k, "") for k in FIELDS})


def verify(key: str, path: Path, url: str, kind: str, repin: bool = False) -> str:
    """Check path against the lock; add the entry on first use. Return the hash."""
    digest = sha256(path)
    lock = read_lock()
    entry = lock.get(key)
    if entry is None or (repin and entry["sha256"] != digest):
        lock[key] = {
            "key": key, "url": url, "retrieved": dt.datetime.now(dt.UTC).date().isoformat(),
            "bytes": str(path.stat().st_size), "sha256": digest,
            "licence": LICENCE[kind], "citation": CITATION[kind],
        }
        write_lock(lock)
    elif entry["sha256"] != digest:
        raise LockMismatch(
            f"{key}: the file in local/ does not match the lock (lock {entry['sha256'][:12]}, "
            f"file {digest[:12]}). The source may have changed. Re-pin only with --repin and a "
            f"DEVIATIONS.md entry.")
    return digest


def _save(path: Path, body: bytes) -> None:
    """Write atomically, so an interrupted download never leaves a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_bytes(body)
    tmp.replace(path)


def _ncep_session():
    import requests
    s = requests.Session()
    s.headers["User-Agent"] = net.USER_AGENT
    return s


def ncep_year_array(year: int, session=None) -> np.ndarray:
    """Return the JJA 250 hPa meridional wind for one year as float32 (day, lat, lon)."""
    import xarray as xr
    ds = xr.open_dataset(NCEP_URL.format(year=year), engine="pydap", session=session or _ncep_session())
    v = ds["vwnd"].sel(level=params.LEVEL_HPA, lat=slice(params.LAT_BAND[1], params.LAT_BAND[0]),
                       time=slice(f"{year}-06-01", f"{year}-08-31")).load()
    days = [str(t)[:10] for t in v.time.values]
    if len(days) != 92 or days[0] != f"{year}-06-01" or days[-1] != f"{year}-08-31":
        raise ValueError(f"{year}: unexpected JJA days {days[0]}..{days[-1]} ({len(days)})")
    if not (np.allclose(v.lat.values, NCEP_LATS) and np.allclose(v.lon.values, NCEP_LONS)):
        raise ValueError(f"{year}: unexpected grid")
    return np.asarray(v.values, dtype=np.float32)


def ncep_path(year: int) -> Path:
    return LOCAL / "ncep" / f"vwnd250_jja_{year}.npy"


def ensure_ncep(year: int, offline: bool = False, repin: bool = False, session=None) -> Path:
    path = ncep_path(year)
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        # pydap makes 3-4 requests per year (dds, das, data), so the per-request gap of net.py
        # does not apply; after HTTP 429 wait 30 s and retry, at most 3 times (hard rule).
        for attempt in range(4):
            try:
                array = ncep_year_array(year, session)
                break
            except Exception as error:
                if "429" not in str(error) or attempt == 3:
                    raise
                time.sleep(30)
        buf = io.BytesIO()
        np.save(buf, array, allow_pickle=False)
        _save(path, buf.getvalue())
        time.sleep(3 * net.MIN_GAP)
    url = NCEP_URL.format(year=year) + " [vwnd, level 250, lat 57.5..37.5, JJA]"
    verify(f"ncep_vwnd250_jja_{year}", path, url, "ncep", repin)
    return path


def ensure_url(key: str, url: str, path: Path, kind: str, offline: bool = False,
               repin: bool = False) -> Path:
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        _save(path, net.get(url))
    verify(key, path, url, kind, repin)
    return path


def fao_path() -> Path:
    return LOCAL / "fao" / "Production_Crops_Livestock_E_All_Data_(Normalized).zip"


def cru_path(var: str, a: int, b: int) -> Path:
    return LOCAL / "cru" / f"cru_ts4.10.{a}.{b}.{var}.dat.nc.gz"


def ne_path() -> Path:
    return LOCAL / "ne" / "ne_110m_admin_0_countries.geojson"


def run(offline: bool = False, repin: bool = False, surface: bool = True, log=print) -> None:
    """Fetch and verify every source. Failures are logged and re-raised at the end."""
    failures: list[str] = []
    session = None if offline else _ncep_session()
    for year in range(params.WAVE_BASE[0], params.LAST_YEAR + 1):
        try:
            ensure_ncep(year, offline, repin, session)
        except LockMismatch:
            raise
        except Exception as error:  # logged, then reported together
            failures.append(f"ncep {year}: {error}")
            log(f"  failed: ncep {year}: {error}")
    jobs = [("fao_qcl_bulk", FAO_URL, fao_path(), "fao"), ("ne_110m_admin0", NE_URL, ne_path(), "ne")]
    if surface:
        for var in CRU_VARS:
            for a, b in CRU_DECADES:
                jobs.append((f"cru_{var}_{a}_{b}", CRU_URL.format(var=var, a=a, b=b),
                             cru_path(var, a, b), "cru"))
    for key, url, path, kind in jobs:
        try:
            ensure_url(key, url, path, kind, offline, repin)
        except LockMismatch:
            raise
        except Exception as error:
            failures.append(f"{key}: {error}")
            log(f"  failed: {key}: {error}")
    if failures:
        raise RuntimeError("fetch incomplete:\n" + "\n".join(failures))
