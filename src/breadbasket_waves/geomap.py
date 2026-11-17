"""An interactive breadbasket map for the site, as plain inline SVG (no JavaScript).

The five regions of the paper are shaded on a Northern Hemisphere map drawn from Natural Earth
outlines. Arcs join each pair of regions: red when concurrent low yields are more likely in
summers with two or more wave event weeks (LMF > 1), blue when less likely, thicker for larger
effects. Hovering an arc shows the counts behind it. A CSS-only toggle switches wave-7 and wave-5.
All numbers come from the LMF results file.
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import pandas as pd

from . import params

REGION_COUNTRIES = {
    "NA": {"USA", "CAN"},
    "WEU": {"FRA", "CHE", "ESP", "PRT", "GBR", "BEL", "NLD", "DEU"},
    "EEU": {"GRC", "BGR", "MDA", "MKD", "UKR", "ROU", "SRB", "ALB", "RUS"},
    "IND": {"IND"},
    "EAS": {"CHN", "MNG"},
}
REGION_NAMES = {"NA": "North America", "WEU": "Western Europe", "EEU": "Eastern Europe",
                "IND": "India", "EAS": "Eastern Asia"}
# Label and arc anchor points (lon, lat), placed in the main growing areas.
ANCHORS = {"NA": (-97, 42), "WEU": (2, 47), "EEU": (36, 50), "IND": (78, 24), "EAS": (113, 36)}
LON0, LON1, LAT0, LAT1 = -170.0, 150.0, 8.0, 74.0
SCALE = 3.0


def _xy(lon: float, lat: float) -> tuple[float, float]:
    return (round((lon - LON0) * SCALE, 1), round((LAT1 - lat) * SCALE, 1))


def _path(rings) -> str:
    parts = []
    for ring in rings:
        pts = [_xy(lon, lat) for lon, lat in ring]
        parts.append("M" + " L".join(f"{x},{y}" for x, y in pts) + " Z")
    return " ".join(parts)


def _countries(geojson: Path) -> str:
    data = json.loads(geojson.read_text(encoding="utf-8"))
    region_of = {c: r for r, cs in REGION_COUNTRIES.items() for c in cs}
    out = []
    for feat in sorted(data["features"], key=lambda f: f["properties"].get("ADM0_A3", "")):
        code = feat["properties"].get("ADM0_A3", "")
        geom = feat["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        rings = [p[0] for p in polys if any(LAT0 - 10 <= lat <= LAT1 + 10 and LON0 - 20 <= lon <= LON1 + 30
                                            for lon, lat in p[0])]
        if not rings:
            continue
        region = region_of.get(code)
        cls = f"land r-{region}" if region else "land"
        name = html.escape(feat["properties"].get("NAME", code))
        out.append(f'<path class="{cls}" d="{_path(rings)}"><title>{name}</title></path>')
    return "\n".join(out)


def _arc(a: str, b: str, row: pd.Series) -> str:
    (x1, y1), (x2, y2) = _xy(*ANCHORS[a]), _xy(*ANCHORS[b])
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    lift = 0.22 * math.hypot(x2 - x1, y2 - y1)
    cx, cy = round(mx, 1), round(my - lift, 1)
    value = float(row.lmf)
    if not math.isfinite(value) or value <= 0:
        width, cls = 1.0, "arc arc-na"
    else:
        width = round(1.0 + 3.2 * min(abs(math.log2(value)), 2.5), 1)
        cls = "arc arc-up" if value > 1 else "arc arc-down"
    tip = (f"{REGION_NAMES[a]} x {REGION_NAMES[b]}: both below trend in {int(row.hits_multi)} of "
           f"{int(row.n_multi)} wave summers vs {int(row.hits_zero)} of {int(row.n_zero)} quiet "
           f"summers. LMF {value:.2f} (95% CI {row.ci_low:.2f}-{row.ci_high:.2f}, exact p {row.p_exact:.2f})")
    return (f'<path class="{cls}" stroke-width="{width}" d="M{x1},{y1} Q{cx},{cy} {x2},{y2}">'
            f"<title>{html.escape(tip)}</title></path>")


def _layer(table: pd.DataFrame, k: int) -> str:
    t = table[(table.wavenumber == k) & (table.direction == "low")]
    arcs = []
    for r in t.itertuples():
        arcs.append(_arc(r.region_1, r.region_2, pd.Series(r._asdict())))
    return f'<g class="layer layer-w{k}">' + "\n".join(arcs) + "</g>"


def svg(table: pd.DataFrame, geojson: Path) -> str:
    w, h = round((LON1 - LON0) * SCALE), round((LAT1 - LAT0) * SCALE)
    labels = []
    for region, (lon, lat) in ANCHORS.items():
        x, y = _xy(lon, lat)
        labels.append(f'<circle class="node" cx="{x}" cy="{y}" r="5"/>'
                      f'<text class="lbl" x="{x}" y="{round(y + 17, 1)}">{REGION_NAMES[region]}</text>')
    band = (_xy(LON0, params.LAT_BAND[1]), _xy(LON1, params.LAT_BAND[0]))
    return (f'<svg class="geomap" viewBox="0 0 {w} {h}" role="img" aria-label="Breadbasket map">'
            f'<rect class="sea" width="{w}" height="{h}"/>'
            f'<rect class="band" x="0" y="{band[0][1]}" width="{w}" height="{round(band[1][1] - band[0][1], 1)}">'
            f"<title>37.5-57.5 N: the latitude band of the wave analysis</title></rect>"
            f"{_countries(geojson)}{_layer(table, 7)}{_layer(table, 5)}{''.join(labels)}</svg>")


CSS = """.mapcard { background:var(--card); border:1px solid var(--line); border-radius:16px; padding:14px; margin:1.6em 0; }
.mapcard input { position:absolute; opacity:0; pointer-events:none; }
.toggle { display:flex; gap:8px; margin:0 0 10px; flex-wrap:wrap; }
.toggle label { cursor:pointer; padding:5px 14px; border-radius:99px; border:1px solid var(--line); font-size:0.9rem; color:var(--muted); }
#w7:checked ~ .toggle label[for=w7], #w5:checked ~ .toggle label[for=w5] { background:var(--accent); color:var(--bg); border-color:transparent; }
.layer-w5 { display:none; } #w5:checked ~ svg .layer-w5 { display:inline; } #w5:checked ~ svg .layer-w7 { display:none; }
.geomap { width:100%; height:auto; display:block; border-radius:10px; }
.geomap .sea { fill:var(--soft); } .geomap .band { fill:var(--accent); opacity:0.07; }
.geomap .land { fill:var(--card); stroke:var(--line); stroke-width:0.6; }
.geomap .r-NA, .geomap .r-WEU, .geomap .r-EEU, .geomap .r-IND, .geomap .r-EAS { fill:var(--accent); fill-opacity:0.28; stroke:var(--accent); stroke-opacity:0.5; }
.geomap .land:hover { fill-opacity:0.5; }
.geomap .arc { fill:none; stroke-linecap:round; opacity:0.85; cursor:help; }
.geomap .arc:hover { opacity:1; filter:drop-shadow(0 0 3px rgba(0,0,0,0.35)); }
.geomap .arc-up { stroke:#c2410c; } .geomap .arc-down { stroke:#2563eb; } .geomap .arc-na { stroke:var(--muted); stroke-dasharray:3 3; }
.geomap .node { fill:var(--ink); } .geomap .lbl { font-size:12px; fill:var(--ink); text-anchor:middle; font-weight:600;
  paint-order:stroke; stroke:var(--card); stroke-width:3px; }
.legend { color:var(--muted); font-size:0.86rem; margin:8px 2px 0; }
.legend .up { color:#c2410c; font-weight:600; } .legend .down { color:#2563eb; font-weight:600; }"""


def card(table: pd.DataFrame, geojson: Path, period: str) -> str:
    return (f'<div class="mapcard"><input type="radio" name="wave" id="w7" checked>'
            f'<input type="radio" name="wave" id="w5">'
            f'<div class="toggle"><label for="w7">wave-7</label><label for="w5">wave-5</label></div>'
            f"{svg(table, geojson)}"
            f'<p class="legend">Arcs join breadbasket regions, {period}. <span class="up">Red</span>: '
            f'both regions below their yield trend more often in summers with two or more wave event '
            f'weeks (LMF above 1). <span class="down">Blue</span>: less often. Thicker means a larger '
            f'factor. Hover an arc for the counts behind it; single arcs rest on about 15 and 20 years. '
            f'Shaded band: 37.5-57.5 N, where the waves are measured.</p></div>')
