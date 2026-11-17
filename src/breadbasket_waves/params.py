"""Method settings (fixed in the study plan before any result) and the published numbers this
study compares against, each with its source and the source's own words.

The targets feed the comparison table in RESULTS.md. No test asserts agreement with them:
agreement or disagreement is reported, never forced.
"""
from __future__ import annotations

import zlib
from dataclasses import dataclass

import numpy as np

# --- Sources -----------------------------------------------------------------------------------

KORNHUBER_2023 = (
    "Kornhuber, K., Lesk, C., Schleussner, C.-F., Jägermeyr, J., Pfleiderer, P. & Horton, R. M. "
    "(2023). Risks of synchronized low yields are underestimated in climate and crop model "
    "projections. Nature Communications 14, 3528. https://doi.org/10.1038/s41467-023-38906-7"
)
KORNHUBER_2020 = (
    "Kornhuber, K., Coumou, D., Vogel, E., Lesk, C., Donges, J. F., Lehmann, J. & Horton, R. M. "
    "(2020). Amplified Rossby waves enhance risk of concurrent heatwaves in major breadbasket "
    "regions. Nature Climate Change 10, 48-53. https://doi.org/10.1038/s41558-019-0637-z"
)
KORNHUBER_2019 = (
    "Kornhuber, K., Osprey, S., Coumou, D., Petri, S., Petoukhov, V., Rahmstorf, S. & Gray, L. "
    "(2019). Extreme weather events in early summer 2018 connected by a recurrent hemispheric "
    "wave-7 pattern. Environmental Research Letters 14, 054002. "
    "https://doi.org/10.1088/1748-9326/ab13bf"
)

# --- Method settings (from the plan; change only through DEVIATIONS.md) -------------------------

LEVEL_HPA = 250.0
LAT_BAND = (37.5, 57.5)          # degrees north
WAVENUMBERS = (5, 7)
N_WEEKS = 13                     # 13 non-overlapping 7-day means: 1 June to 30 August
THRESHOLD_SD = 1.5               # event: weekly amplitude above mean + 1.5 SD (ddof = 1)
WAVE_BASE = (1960, 2014)         # the paper's wave climatology period
YIELD_BASE = (1961, 2014)        # FAOSTAT starts in 1961
LAST_YEAR = 2024                 # latest FAOSTAT year (checked 17 Nov 2026)
REGIONS = ("NA", "WEU", "EEU", "IND", "EAS")
EEU_BREAK = 1992                 # USSR and Yugoslav SFR rows end in 1991
SSA_WINDOW = 10
SSA_COMPONENTS = 2               # 2, not 1: one component leaves up to 2.5 % of a linear trend (DEVIATIONS.md)
NEIGHBOUR_YEARS = 3              # composites: ratio to the mean of 3 years before and 3 after
N_PERMUTATIONS = 10_000
N_SUBSAMPLES = 500
SUBSAMPLE_FRACTION = 0.9
N_BOOTSTRAP = 2_000
SEED = 20261117


def rng_for(name: str) -> np.random.Generator:
    """A random stream that depends only on the base seed and the name (never on hash())."""
    return np.random.default_rng([SEED, zlib.crc32(name.encode("utf-8"))])


# --- Published targets -------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    key: str
    value: float | int | str
    unit: str
    what: str
    source: str
    where: str
    quote: str


_M = "Methods"
_R = "Results"
_Q_EVENTS = "For ERA5 that is 54 wave-7 events and 62 wave-5 for JJA 1960–2014"
_Q_LMF = ("During years with more than one wave-7 event the probability of concurrent low yields "
          "is increased by a factor of 1.6 or 60% compared to years with no wave-7 event")
_Q_LMF_N = ("Throughout the dataset Nw>1=15 years are identified in which more than one wave event "
            "occurs")
_Q_COMP7 = ("years with more than one wave event are associated with regional crop yield anomalies "
            "of up to −7% in EAS for, −6% to in NA and −3% in EEU for wave-7 with an average "
            "response across the selected regions of −2% to −3%")
_Q_COMP5 = ("For wave-5 we observe values of −1% to −3% for respectively with highest regional "
            "anomalies observed in EEU (−3.7%) and NA (−2%)")

TARGETS: tuple[Target, ...] = (
    Target("events_w7_era5", 54, "event weeks", "High-amplitude wave-7 weeks, JJA, ERA5, 1960-2014", KORNHUBER_2023, _M, _Q_EVENTS),
    Target("events_w5_era5", 62, "event weeks", "High-amplitude wave-5 weeks, JJA, ERA5, 1960-2014", KORNHUBER_2023, _M, _Q_EVENTS),
    Target("lmf_example_n_multi", 15, "years", "Years with more than one wave-7 event (ERA5, FAO)", KORNHUBER_2023, _M, _Q_LMF_N),
    Target("lmf_example_hits_multi", 6, "years", "Of those, years with below-trend yields in both NA and WEU", KORNHUBER_2023, _M, _Q_LMF_N),
    Target("lmf_example_n_zero", 20, "years", "Years with no wave-7 event", KORNHUBER_2023, _M, _Q_LMF_N),
    Target("lmf_example_hits_zero", 5, "years", "Of those, years with below-trend yields in both NA and WEU", KORNHUBER_2023, _M, _Q_LMF_N),
    Target("lmf_example", 1.6, "ratio", "LMF for concurrent low yields, wave-7, NA x WEU", KORNHUBER_2023, _M, _Q_LMF),
    Target("comp_w7_EAS", -7, "% (up to)", "Composite yield anomaly, years with more than one wave-7 event, EAS", KORNHUBER_2023, _R, _Q_COMP7),
    Target("comp_w7_NA", -6, "% (up to)", "Composite yield anomaly, wave-7, NA", KORNHUBER_2023, _R, _Q_COMP7),
    Target("comp_w7_EEU", -3, "% (up to)", "Composite yield anomaly, wave-7, EEU", KORNHUBER_2023, _R, _Q_COMP7),
    Target("comp_w7_mean", "-2 to -3", "%", "Composite yield anomaly, wave-7, average across regions", KORNHUBER_2023, _R, _Q_COMP7),
    Target("comp_w5_range", "-1 to -3", "%", "Composite yield anomalies, wave-5, range across regions", KORNHUBER_2023, _R, _Q_COMP5),
    Target("comp_w5_EEU", -3.7, "%", "Composite yield anomaly, wave-5, EEU", KORNHUBER_2023, _R, _Q_COMP5),
    Target("comp_w5_NA", -2, "%", "Composite yield anomaly, wave-5, NA", KORNHUBER_2023, _R, _Q_COMP5),
)

# Summers that Kornhuber et al. (2019) name as showing the amplified wave-7 pattern (abstract):
# "This pattern was also observed during the European heatwaves of 2003, 2006 and 2015 among
# others", plus 2018, the subject of the paper. Reported as a known case, never tested.
WAVE7_NAMED_SUMMERS = (2003, 2006, 2015, 2018)


def target(key: str) -> Target:
    for t in TARGETS:
        if t.key == key:
            return t
    raise KeyError(key)
