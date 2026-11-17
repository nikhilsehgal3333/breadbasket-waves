"""Figures, written as PNG with no build date or software tag, so they are byte-stable."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .paths import FIGURES  # noqa: E402

INK = "#1f2933"
MUTED = "#7b8794"
WAVE_COLOURS = {5: "#2f6f9f", 7: "#c2571a"}
plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "svg.hashsalt": "breadbasket-waves",
})


def _save(fig, name: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)
    return path


def fig_wave_events(counts: pd.DataFrame, events: pd.DataFrame, base_end: int, named: tuple[int, ...]) -> Path:
    """Event weeks per JJA for wave-5 and wave-7, and the week-of-summer timing of events."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 5.2), gridspec_kw={"width_ratios": [3.2, 1]})
    for row, k in enumerate((5, 7)):
        ax = axes[row, 0]
        years = counts.index.to_numpy()
        ax.bar(years, counts[k].to_numpy(), color=WAVE_COLOURS[k], width=0.8)
        ax.axvline(base_end + 0.5, color=MUTED, lw=0.8, ls="--")
        ax.set_ylabel(f"wave-{k} event weeks")
        ax.set_xlim(years.min() - 1, years.max() + 1)
        if k == 7:
            ax.scatter(list(named), [-0.35] * len(named), marker="^", s=22, color=INK, clip_on=False, zorder=4)
            ax.set_ylim(bottom=0)
        ax.text(base_end + 1, ax.get_ylim()[1] * 0.92, "extension", color=MUTED, fontsize=8)
        timing = events[(events.wavenumber == k) & events.event].groupby("week").size()
        timing = timing.reindex(range(1, 14), fill_value=0)
        axes[row, 1].bar(timing.index, timing.values, color=WAVE_COLOURS[k], width=0.8)
        axes[row, 1].set_xticks([1, 5, 9, 13], ["1 Jun", "29 Jun", "27 Jul", "24 Aug"])
        axes[row, 1].set_ylabel("events, all years")
    axes[0, 0].set_title("High-amplitude wave weeks per summer (NCEP/NCAR R1, 250 hPa, 37.5-57.5 N)")
    axes[0, 1].set_title("Week of summer")
    axes[1, 0].set_xlabel("year (triangles: wave-7 summers named by Kornhuber et al. 2019)")
    fig.tight_layout()
    return _save(fig, "fig1_wave_events.png")


def fig_regional_yields(table: pd.DataFrame, trends: dict[str, pd.Series], regions: tuple[str, ...]) -> Path:
    fig, axes = plt.subplots(len(regions), 1, figsize=(8, 1.7 * len(regions)), sharex=True)
    for ax, region in zip(axes, regions):
        s = table[table.region == region].set_index("year").yield_t_ha
        ax.plot(s.index, s.values, color=INK, lw=1)
        t = trends[region]
        ax.plot(t.index, t.values, color="#c2571a", lw=1)
        below = s[s < t.reindex(s.index)]
        ax.scatter(below.index, below.values, s=6, color="#c2571a", zorder=3)
        ax.set_ylabel(f"{region}\nt/ha")
    axes[0].set_title("Combined maize and wheat yield (FAOSTAT) with SSA trend; dots: below trend")
    axes[-1].set_xlabel("year")
    fig.tight_layout()
    return _save(fig, "fig2_regional_yields.png")


def fig_lmf(table: pd.DataFrame, name: str, period: str) -> Path:
    """Heatmaps as in Kornhuber et al. (2023), Fig. 4: concurrent low yields above the diagonal,
    concurrent high yields below. Cells show the LMF and the counts behind it."""
    regions = ("NA", "WEU", "EEU", "IND", "EAS")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    cmap = plt.get_cmap("RdBu_r")
    for ax, k in zip(axes, (7, 5)):
        grid = np.full((5, 5), np.nan)
        labels = [["" for _ in regions] for _ in regions]
        t = table[table.wavenumber == k]
        for _, r in t.iterrows():
            i, j = regions.index(r.region_1), regions.index(r.region_2)
            a, b = (i, j) if r.direction == "low" else (j, i)
            value = r.lmf
            grid[a, b] = np.log2(value) if np.isfinite(value) and value > 0 else np.nan
            shown = "inf" if np.isinf(value) else f"{value:.2f}"
            labels[a][b] = f"{shown}\n{int(r.hits_multi)}/{int(r.n_multi)} v {int(r.hits_zero)}/{int(r.n_zero)}"
        ax.imshow(grid, cmap=cmap, vmin=-1.5, vmax=1.5)
        for a in range(5):
            for b in range(5):
                if labels[a][b]:
                    strong = np.isfinite(grid[a, b]) and abs(grid[a, b]) > 0.9
                    ax.text(b, a, labels[a][b], ha="center", va="center", fontsize=7,
                            color="white" if strong else INK)
        ax.set_xticks(range(5), regions)
        ax.set_yticks(range(5), regions)
        ax.set_title(f"wave-{k}, {period}\nabove diagonal: low yields; below: high yields", fontsize=9)
    fig.suptitle("Likelihood multiplication factor (LMF); colour log2 scale, red = more likely in "
                 "multi-event years", fontsize=9)
    fig.tight_layout()
    return _save(fig, name)


def fig_composites(comp: pd.DataFrame, name: str) -> Path:
    regions = ("NA", "WEU", "EEU", "IND", "EAS", "ALL", "WORLD")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, k in zip(axes, (7, 5)):
        c = comp[comp.wavenumber == k].set_index("region").loc[list(regions)]
        x = np.arange(len(regions))
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.bar(x - 0.18, c.multi_mean, 0.36, color=WAVE_COLOURS[k], label="2 or more event weeks")
        ax.bar(x + 0.18, c.zero_mean, 0.36, color="#c9ced6", label="no event week")
        ax.errorbar(x - 0.18, c.multi_mean, yerr=[c.multi_mean - c.multi_lo, c.multi_hi - c.multi_mean],
                    fmt="none", ecolor=INK, lw=0.8, capsize=2)
        ax.scatter(x - 0.18, c.multi_median, marker="_", s=80, color="white", zorder=3)
        ax.set_xticks(x, regions)
        ax.set_title(f"wave-{k}: yield vs the 3 years either side (%)")
    axes[0].set_ylabel("% difference")
    axes[0].legend(fontsize=7, frameon=False)
    fig.text(0.01, -0.02, "Bars: mean. Whiskers: 95 % bootstrap interval. White tick: median of the "
             "paper-style 90 % subsamples (a spread, not an interval).", fontsize=7, color=MUTED)
    fig.tight_layout()
    return _save(fig, name)


def _outlines(ax, geojson_path: Path) -> None:
    import json
    data = json.loads(Path(geojson_path).read_text(encoding="utf-8"))
    for feature in data["features"]:
        geom = feature["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            ring = np.asarray(poly[0])
            ax.plot(ring[:, 0], ring[:, 1], color=INK, lw=0.35)


def fig_surface(maps: dict, n_months: dict, geojson_path: Path) -> Path:
    """maps[(k, var)] -> 2-D DataArray (lat, lon)."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 4.6))
    specs = {"tmx": ("RdBu_r", 1.5, "max temperature anomaly (deg C)"),
             "pre": ("BrBG", 15.0, "precipitation anomaly (mm/month)")}
    for row, k in enumerate((7, 5)):
        for col, var in enumerate(("tmx", "pre")):
            ax = axes[row, col]
            cmap, lim, label = specs[var]
            m = maps[(k, var)]
            im = ax.pcolormesh(m.lon, m.lat, m.values, cmap=cmap, vmin=-lim, vmax=lim, shading="auto")
            _outlines(ax, geojson_path)
            ax.set_xlim(-170, 180)
            ax.set_ylim(20, 75)
            ax.axhspan(37.5, 57.5, color="none", ec=MUTED, ls="--", lw=0.6)
            ax.set_title(f"wave-{k} event months (n = {n_months[k]}): {label}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([30, 45, 60])
            fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    fig.text(0.01, 0.0, "CRU TS 4.10, JJA months containing event weeks, 1960-2014; anomalies from a "
             "per-cell linear trend (1961-2014); cells with no station masked. Dashed: the 37.5-57.5 N "
             "band of the wave analysis. Outlines: Natural Earth.", fontsize=7, color=MUTED)
    fig.tight_layout()
    return _save(fig, "fig5_surface_composites.png")
