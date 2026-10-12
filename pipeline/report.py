"""Build the survey report page (results/report/index.html) from the results files.

  python report.py

Numbers come from results/survey/*/summary.json, results/candidates/*, results/report/transit_fits.json and
results/report/extras.json; the narrative for the headline candidates lives in CURATED below. Figures are
drawn as inline SVG styled by the page's colour tokens, so they follow the viewer's light/dark theme.
"""
import glob
import html
import json
import re
import math
from datetime import datetime, timezone

import env  # noqa: F401
import numpy as np
import pandas as pd

import prior
from env import CATALOGS, RESULTS

OUT = RESULTS / "report"
CAND = RESULTS / "candidates"
BTJD0 = 2457000.0
R_EARTH_PER_RSUN = 109.08
RE = "R<sub>⊕</sub>"
RSUN = "R<sub>☉</sub>"
MSUN = "M<sub>☉</sub>"
RHOSUN = "ρ<sub>☉</sub>"
REPO = "https://github.com/specbug/ad-astra/blob/claude/gifted-albattani-szay6x"


def esc(x):
    return html.escape(str(x), quote=True)


def fmt(x, nd=2, dash="–"):
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return dash
    except TypeError:
        return dash
    if isinstance(x, (int, np.integer)):
        return f"{int(x):,}"
    return f"{x:,.{nd}f}"


def btjd_to_utc(t):
    jd = t + BTJD0
    return datetime.fromtimestamp((jd - 2440587.5) * 86400.0, tz=timezone.utc)


# ----------------------------------------------------------------------------------------------- data

def survey_stats():
    rows = []
    for f in sorted(glob.glob(str(RESULTS / "survey" / "s*" / "summary.json"))):
        s = json.load(open(f))
        name = f.split("/")[-2]
        rows.append({"sector": s["sector"], "product": "spoc" if "spoc" in name else "qlp",
                     "n_real": s["rates"]["False"]["n_lc"], "n_flip": s["rates"]["True"]["n_lc"],
                     "n_events": s["n_events"], "n_clean": s["n_clean"], "n_trackb": s["n_trackb"],
                     "n_trackc": s["n_trackc"], "n_trackd": s["n_trackd"]})
    return pd.DataFrame(rows)


def stage2_stats(sv):
    import stage2
    p = CAND / "trackb_stage2.parquet"
    if not p.exists():
        return None
    d = stage2.classify_shape(pd.read_parquet(p))
    vetted = set(zip(d.sector, d["product"]))
    lc = sv[[(s, pr) in vetted for s, pr in zip(sv.sector, sv["product"])]]
    pl = d[d.planet_like & (d.Teff >= 4000) & (d.fit_t14_h > 2)]
    out = {"n_vetted_real": int((~d.invert.astype(bool)).sum()), "n_vetted_flip": int(d.invert.astype(bool).sum()),
           "n_pl_real": int((d.planet_like & ~d.invert.astype(bool)).sum()),
           "n_pl_flip": int((d.planet_like & d.invert.astype(bool)).sum()), "by_product": {}, "by_sector": []}
    for prod in ("qlp", "spoc"):
        n = lc[lc["product"] == prod][["n_real", "n_flip"]].sum()
        x = pl[pl["product"] == prod]
        r, f = int((~x.invert.astype(bool)).sum()), int(x.invert.astype(bool).sum())
        out["by_product"][prod] = {"n_real_lc": int(n.n_real), "n_flip_lc": int(n.n_flip), "real": r, "flip": f}
    q = pl[pl["product"] == "qlp"]
    for s in sorted(set(lc[lc["product"] == "qlp"].sector)):
        n = lc[(lc["product"] == "qlp") & (lc.sector == s)].iloc[0]
        x = q[q.sector == s]
        out["by_sector"].append({"sector": int(s), "real": int((~x.invert.astype(bool)).sum()),
                                 "flip": int(x.invert.astype(bool).sum()), "n_real_lc": int(n.n_real),
                                 "n_flip_lc": int(n.n_flip)})
    return out


def completeness():
    p = CAND / "completeness_toi.csv"
    if not p.exists():
        return None
    c = pd.read_csv(p)
    rows = []
    for lo, hi in ((1500, 3000), (3000, 6000), (6000, 12000), (12000, 1e7)):
        x = c[(c.depth_ppm >= lo) & (c.depth_ppm < hi)]
        if len(x):
            f = x[x.found]
            rows.append({"bin": f"{lo/1e3:g}–{hi/1e3:g} ppt" if hi < 1e7 else f"> {lo/1e3:g} ppt", "n": len(x),
                         "found": float(x.found.mean()),
                         "pl": float(f.planet_like.mean()) if len(f) and f.planet_like.notna().any() else float("nan")})
    return {"n": len(c), "found": float(c.found.mean()), "rows": rows}


# ------------------------------------------------------------------------------------------- figures

def svg_lightcurves(panels, labels, ylim=None, height=210, panel_w=250, id_="lc", xr=None, xstep=None):
    """Small multiples of binned transit photometry (ppm vs hours from mid-transit) with the model."""
    n = len(panels)
    padl, padr, padt, padb = 46, 10, 22, 34
    W = n * panel_w
    H = height
    allx = np.concatenate([p["x"] for p in panels])
    if xr is None:
        xr = float(np.ceil(max(abs(allx.min()), abs(allx.max()))))
    if ylim is None:
        ally = np.concatenate([np.asarray(p["y"])[np.abs(np.asarray(p["x"])) <= xr] for p in panels])
        lo, hi = np.percentile(ally, [1, 99])
        ylim = (min(lo, min(min(p["my"]) for p in panels)) - 150, hi + 150)
    y0, y1 = ylim
    step = 250 if (y1 - y0) < 2200 else 500
    yt = np.arange(np.ceil(y0 / step) * step, y1 + 1, step)
    xstep = xstep or (6 if xr > 10 else 4 if xr > 6 else 2)
    xt = np.arange(-int(xr // xstep) * xstep, int(xr) + 1, xstep)
    parts = [f'<svg class="lc" id="{id_}" viewBox="0 0 {W} {H}" role="img" '
             f'aria-label="Binned TESS photometry around each transit, with the fitted transit model">']
    for i, (p, lab) in enumerate(zip(panels, labels)):
        x0 = i * panel_w + padl
        pw = panel_w - padl - padr
        ph = H - padt - padb

        def X(v):
            return x0 + (v + xr) / (2 * xr) * pw

        def Y(v):
            return padt + (y1 - v) / (y1 - y0) * ph
        parts.append(f'<g class="panel" data-i="{i}">')
        parts.append(f'<text class="ptitle" x="{x0}" y="{padt - 8}">{esc(lab)}</text>')
        for v in yt:
            parts.append(f'<line class="grid" x1="{x0}" x2="{x0 + pw}" y1="{Y(v):.1f}" y2="{Y(v):.1f}"/>')
            if i == 0:
                parts.append(f'<text class="tick" x="{x0 - 6}" y="{Y(v) + 4:.1f}" text-anchor="end">{int(v)}</text>')
        for v in xt:
            parts.append(f'<text class="tick" x="{X(v):.1f}" y="{H - padb + 16}" text-anchor="middle">{int(v):+d}</text>'
                         .replace("+0<", "0<"))
        parts.append(f'<line class="axis" x1="{x0}" x2="{x0 + pw}" y1="{padt + ph}" y2="{padt + ph}"/>')
        for xv, yv, ev in zip(p["x"], p["y"], p["e"]):
            if not (-xr <= xv <= xr):
                continue
            yv_c = min(max(yv, y0), y1)
            parts.append(f'<line class="err" x1="{X(xv):.1f}" x2="{X(xv):.1f}" y1="{Y(min(yv + ev, y1)):.1f}" '
                         f'y2="{Y(max(yv - ev, y0)):.1f}"/>')
            parts.append(f'<circle class="pt" cx="{X(xv):.1f}" cy="{Y(yv_c):.1f}" r="2.6" '
                         f'data-x="{xv:.2f}" data-y="{yv:.0f}" data-e="{ev:.0f}"/>')
        pts = " ".join(f"{X(a):.1f},{Y(min(max(b, y0), y1)):.1f}" for a, b in zip(p["mx"], p["my"]) if -xr <= a <= xr)
        parts.append(f'<polyline class="model" points="{pts}"/>')
        parts.append("</g>")
    parts.append(f'<text class="axlabel" x="{W / 2}" y="{H - 4}" text-anchor="middle">hours from mid-transit</text>')
    parts.append(f'<text class="axlabel" transform="translate(11,{padt + (H - padt - padb) / 2}) rotate(-90)" '
                 f'text-anchor="middle">flux − 1 (ppm)</text>')
    parts.append("</svg>")
    return "".join(parts)


# ------------------------------------------------------------------------------------------ helpers

def chip(text, kind):
    return f'<span class="chip chip-{kind}">{esc(text)}</span>'


def check_rows(items):
    icons = {"pass": "✓", "warn": "!", "fail": "✕", "info": "i"}
    out = ['<ul class="checks">']
    for kind, title, body in items:
        out.append(f'<li class="check check-{kind}"><span class="ci" aria-hidden="true">{icons[kind]}</span>'
                   f'<span class="sr">{ {"pass": "Passed", "warn": "Caution", "fail": "Failed", "info": "Note"}[kind]}: </span>'
                   f'<div><strong>{title}</strong> {body}</div></li>')
    out.append("</ul>")
    return "".join(out)


def table(df, cols, cls="data"):
    """cols: list of (header, column or callable(row) -> html, align)."""
    h = [f'<div class="tablewrap"><table class="{cls}"><thead><tr>']
    for head, _, al in cols:
        h.append(f'<th class="{al}">{head}</th>')
    h.append("</tr></thead><tbody>")
    for _, r in df.iterrows():
        h.append("<tr>")
        for _, c, al in cols:
            v = c(r) if callable(c) else esc(r[c])
            h.append(f'<td class="{al}">{v}</td>')
        h.append("</tr>")
    h.append("</tbody></table></div>")
    return "".join(h)


def pm(v, nd=2, unit=""):
    """value with asymmetric errors [med, +, -]"""
    if v is None:
        return "–"
    m, up, lo = v
    if abs(up - lo) / max(abs(up), abs(lo), 1e-12) < 0.25:
        return f"{m:.{nd}f} ± {max(up, lo):.{nd}f}{unit}"
    return f"{m:.{nd}f}<sup>+{up:.{nd}f}</sup><sub>−{lo:.{nd}f}</sub>{unit}"


# ---------------------------------------------------------------------------------------------- style

CSS = r"""
/* Layout: one reading column (~70ch) with figure plates that widen to 1120px; facts sit in a
   two-column spec sheet beside each candidate's evidence on wide screens. */
:root {
  --bg: #f3f5f8; --surface: #ffffff; --ink: #121924; --ink-2: #3a4555; --muted: #657184;
  --line: #d8dee6; --line-2: #e8ecf1;
  --accent: #8f5300; --mark: #c47a12; --data: #4b5768; --blue: #245a8e;
  --good: #23704a; --warn: #9a5c00; --bad: #ad3434;
  --good-bg: #e6f2eb; --warn-bg: #f8eedd; --bad-bg: #f6e3e3; --blue-bg: #e3ecf6; --neutral-bg: #eceff3;
  --font-display: "Archivo", "Helvetica Neue", Arial, sans-serif;
  --font-body: "Source Serif 4", "Iowan Old Style", Georgia, serif;
  --font-mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
  --step--1: 0.84rem; --step-0: 1.06rem; --step-1: 1.25rem; --step-2: 1.6rem; --step-3: 2.2rem; --step-4: 3.1rem;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #0b1017; --surface: #111823; --ink: #e7ebf1; --ink-2: #c2cad5; --muted: #8b97a8;
  --line: #243142; --line-2: #1a2431;
  --accent: #f0b35e; --mark: #f0b35e; --data: #a3afc0; --blue: #86b6e6;
  --good: #73c79a; --warn: #f0b35e; --bad: #f08a8a;
  --good-bg: #13291f; --warn-bg: #2b2112; --bad-bg: #2d1717; --blue-bg: #132235; --neutral-bg: #182130;
  color-scheme: dark; } }
:root[data-theme="dark"] {
  --bg: #0b1017; --surface: #111823; --ink: #e7ebf1; --ink-2: #c2cad5; --muted: #8b97a8;
  --line: #243142; --line-2: #1a2431;
  --accent: #f0b35e; --mark: #f0b35e; --data: #a3afc0; --blue: #86b6e6;
  --good: #73c79a; --warn: #f0b35e; --bad: #f08a8a;
  --good-bg: #13291f; --warn-bg: #2b2112; --bad-bg: #2d1717; --blue-bg: #132235; --neutral-bg: #182130;
  color-scheme: dark; }

* { box-sizing: border-box; }
body { background: var(--bg); color: var(--ink); font-family: var(--font-body); font-size: var(--step-0);
  line-height: 1.62; padding-inline: 16px; padding-block: 0 64px; }
.wrap { max-width: 1120px; margin: 0 auto; }
.col { max-width: 44rem; }
a { color: var(--blue); text-underline-offset: 2px; }
a:focus-visible, button:focus-visible { outline: 2px solid var(--mark); outline-offset: 2px; }
h1, h2, h3, .kicker, .label, th, .chip, .stat .v, nav { font-family: var(--font-display); }
h1, h2, h3 { text-wrap: balance; line-height: 1.15; margin: 0; }
h1 { font-size: var(--step-4); font-weight: 800; font-stretch: 125%; letter-spacing: -0.01em; }
h2 { font-size: var(--step-3); font-weight: 750; font-stretch: 112%; }
h3 { font-size: var(--step-1); font-weight: 700; }
p { margin: 0; }
.mono, code, .num { font-family: var(--font-mono); font-size: 0.92em; }
.num, td.r, .stat .v { font-variant-numeric: tabular-nums; }
.kicker, .label { text-transform: uppercase; letter-spacing: 0.09em; font-size: var(--step--1); color: var(--muted);
  font-weight: 600; }
.muted { color: var(--muted); }

header.mast { padding-block: 56px 28px; display: grid; gap: 14px; border-bottom: 1px solid var(--line); }
header.mast .dek { font-size: var(--step-1); color: var(--ink-2); max-width: 40rem; }
.strip { display: flex; flex-wrap: wrap; gap: 8px 28px; margin-top: 8px; }
.strip div { display: grid; }
.strip .v { font-family: var(--font-mono); font-size: var(--step-1); color: var(--ink); }
nav.toc { display: flex; flex-wrap: wrap; gap: 6px 18px; padding-block: 14px; font-size: var(--step--1);
  border-bottom: 1px solid var(--line); }
nav.toc a { color: var(--ink-2); text-decoration: none; }
nav.toc a:hover { color: var(--accent); }

section { padding-block: 44px 8px; display: grid; gap: 18px; }
section + section { border-top: 1px solid var(--line-2); }
.lede { font-size: var(--step-1); color: var(--ink-2); }

.findings { display: grid; gap: 14px; margin: 0; padding: 0; list-style: none; }
.findings li { display: grid; grid-template-columns: 9.5rem 1fr; gap: 6px 18px; align-items: start;
  padding-block: 12px; border-top: 1px solid var(--line-2); }
.findings li > div { min-width: 0; }
@media (max-width: 640px) { .findings li { grid-template-columns: 1fr; } }

.chip { display: inline-flex; align-items: center; gap: 6px; font-size: 0.74rem; font-weight: 650;
  letter-spacing: 0.04em; text-transform: uppercase; padding: 3px 9px; border-radius: 999px; white-space: nowrap; }
.chip-new { background: var(--warn-bg); color: var(--warn); }
.chip-tce { background: var(--blue-bg); color: var(--blue); }
.chip-known { background: var(--neutral-bg); color: var(--ink-2); }
.chip-eb { background: var(--bad-bg); color: var(--bad); }
.chip-ok { background: var(--good-bg); color: var(--good); }

.cand { display: grid; gap: 22px; }
.cand-head { display: grid; gap: 8px; }
.cand-head .chips { display: flex; flex-wrap: wrap; gap: 8px; }
.cand-grid { display: grid; grid-template-columns: minmax(0, 1fr) 17rem; gap: 28px; align-items: start; }
.cand-grid > div { min-width: 0; display: grid; gap: 16px; }
@media (max-width: 900px) { .cand-grid { grid-template-columns: 1fr; } }
dl.spec { display: grid; grid-template-columns: auto 1fr; gap: 7px 14px; margin: 0; padding: 16px 18px;
  background: var(--surface); border: 1px solid var(--line); border-radius: 10px; font-size: 0.93rem; }
dl.spec dt { color: var(--muted); font-family: var(--font-display); font-size: 0.8rem; text-transform: uppercase;
  letter-spacing: 0.06em; padding-top: 2px; }
dl.spec dd { margin: 0; font-family: var(--font-mono); font-size: 0.86rem; font-variant-numeric: tabular-nums; }
dl.spec dd sup, dl.spec dd sub, td sup, td sub { font-size: 0.7em; }
sub, sup { font-size: 0.62em; line-height: 0; }
.plate { background: var(--surface); border: 1px solid var(--line); border-radius: 10px; padding: 14px 14px 10px;
  display: grid; gap: 8px; position: relative; }
.plate figcaption { font-size: var(--step--1); color: var(--muted); line-height: 1.5; }
.scroll-x { overflow-x: auto; }
svg.lc { width: 100%; min-width: 620px; height: auto; display: block; }
svg.lc .grid { stroke: var(--line-2); stroke-width: 1; }
svg.lc .axis { stroke: var(--line); stroke-width: 1; }
svg.lc .tick { fill: var(--muted); font: 10.5px var(--font-mono); }
svg.lc .ptitle { fill: var(--ink-2); font: 600 11.5px var(--font-display); letter-spacing: 0.02em; }
svg.lc .axlabel { fill: var(--muted); font: 11px var(--font-display); }
svg.lc .err { stroke: var(--data); stroke-opacity: 0.35; stroke-width: 1; }
svg.lc .pt { fill: var(--data); stroke: var(--surface); stroke-width: 1; }
svg.lc .pt.hot { fill: var(--accent); }
svg.lc .model { fill: none; stroke: var(--mark); stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }
.tip { position: absolute; pointer-events: none; background: var(--ink); color: var(--bg); font: 12px var(--font-mono);
  padding: 5px 8px; border-radius: 6px; white-space: nowrap; transform: translate(-50%, -130%); }
.tip b { font-weight: 600; }

ul.checks { list-style: none; margin: 0; padding: 0; display: grid; gap: 10px; }
.check { display: grid; grid-template-columns: 1.6rem 1fr; gap: 10px; align-items: start; font-size: 0.98rem; }
.check .ci { width: 1.4rem; height: 1.4rem; border-radius: 50%; display: grid; place-items: center;
  font: 700 0.8rem var(--font-display); margin-top: 2px; }
.check-pass .ci { background: var(--good-bg); color: var(--good); }
.check-warn .ci { background: var(--warn-bg); color: var(--warn); }
.check-fail .ci { background: var(--bad-bg); color: var(--bad); }
.check-info .ci { background: var(--blue-bg); color: var(--blue); }
.sr { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); }

.tablewrap { overflow-x: auto; border: 1px solid var(--line); border-radius: 10px; background: var(--surface); }
table { border-collapse: collapse; width: 100%; font-size: 0.86rem; }
th { text-align: left; font-size: 0.74rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted);
  font-weight: 600; padding: 10px 12px; border-bottom: 1px solid var(--line); white-space: nowrap; background: var(--surface); }
td { padding: 8px 12px; border-bottom: 1px solid var(--line-2); vertical-align: top; }
td.r, th.r { text-align: right; }
td.m { font-family: var(--font-mono); font-size: 0.82rem; }
tr:last-child td { border-bottom: 0; }
td .chip { font-size: 0.66rem; }

ul.nexts { list-style: none; margin: 0; padding: 0; display: grid; gap: 8px; }
ul.nexts li { display: grid; gap: 1px; font-size: 0.9rem; padding-bottom: 8px; border-bottom: 1px solid var(--line-2); }
ul.nexts li:last-child { border-bottom: 0; padding-bottom: 0; }
.plates { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 16px; }
.plates .plate img { width: 100%; height: auto; border-radius: 6px; background: #fff; }
.callout { border-left: 3px solid var(--mark); padding: 4px 0 4px 16px; color: var(--ink-2); }
.two { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 18px; }
.two > div { min-width: 0; display: grid; gap: 10px; align-content: start; }
footer { padding-block: 40px 0; border-top: 1px solid var(--line); margin-top: 40px; color: var(--muted);
  font-size: var(--step--1); display: grid; gap: 8px; }
@media (max-width: 520px) { h1 { font-size: 2.3rem; } h2 { font-size: 1.7rem; } }
@media (prefers-reduced-motion: no-preference) { nav.toc a { transition: color .15s; } }
"""

JS = r"""
(() => {
  for (const svg of document.querySelectorAll('svg.lc')) {
    const plate = svg.closest('.plate');
    const tip = document.createElement('div'); tip.className = 'tip'; tip.hidden = true; plate.appendChild(tip);
    const pts = [...svg.querySelectorAll('circle.pt')];
    let last = null;
    const show = (ev) => {
      const r = svg.getBoundingClientRect(); const vb = svg.viewBox.baseVal;
      const sx = vb.width / r.width, sy = vb.height / r.height;
      const mx = (ev.clientX - r.left) * sx, my = (ev.clientY - r.top) * sy;
      let best = null, bd = 1e9;
      for (const p of pts) { const dx = p.cx.baseVal.value - mx, dy = p.cy.baseVal.value - my;
        const d = dx * dx + dy * dy * 0.25; if (d < bd) { bd = d; best = p; } }
      if (!best || bd > 900) { tip.hidden = true; if (last) last.classList.remove('hot'); last = null; return; }
      if (last && last !== best) last.classList.remove('hot');
      best.classList.add('hot'); last = best;
      const pr = plate.getBoundingClientRect();
      tip.textContent = '';
      const b = document.createElement('b'); b.textContent = `${best.dataset.y} ± ${best.dataset.e} ppm`;
      tip.append(b, document.createTextNode(`  at ${(+best.dataset.x).toFixed(1)} h`));
      tip.style.left = `${best.cx.baseVal.value / sx + r.left - pr.left}px`;
      tip.style.top = `${best.cy.baseVal.value / sy + r.top - pr.top}px`;
      tip.hidden = false;
    };
    svg.addEventListener('pointermove', show);
    svg.addEventListener('pointerleave', () => { tip.hidden = true; if (last) last.classList.remove('hot'); });
  }
})();
"""


# ----------------------------------------------------------------------------------------- sections

def section_hip71639(fit, ex):
    f = fit["HIP71639"]
    e = ex["HIP71639"]
    rs, rs_err = 1.596, 0.063
    rp = f["rp_rs"][0] * rs * R_EARTH_PER_RSUN
    rp_err = rp * math.hypot(max(f["rp_rs"][1], f["rp_rs"][2]) / f["rp_rs"][0], rs_err / rs)
    depth = -min(min(p["my"]) for p in f["plot"])
    labels = [f"Sector {s} · {btjd_to_utc(p['tc']).strftime('%b %Y')}" for s, p in zip((11, 38, 65, 102), f["plot"])]
    svg = svg_lightcurves(f["plot"], labels, id_="lc-hip71639", xr=13, xstep=6)
    nxt = "".join(f"<li><span class='num'>{btjd_to_utc(r['btjd']).strftime('%-d %b %Y, %H:%M')} UTC</span>"
                  f"<span class='muted'>{esc(r['note'])}</span></li>" for r in e["next"])
    checks = check_rows([
        ("pass", "Four transits, one period.",
         f"Transits in sectors 11, 38, 65 and 102 fit a single linear ephemeris. Every shorter alias "
         f"(P/2 to P/8) predicts extra transits inside the observed data that are not there, ruled out at "
         f"more than {e['alias_min_sigma']:.0f}σ."),
        ("pass", "The shape fits a planet on this star.",
         f"Flat-bottomed, {pm(f['t14_h'], 2)} h long, with short ingress and egress. The stellar density the transit "
         f"shape implies at this period, {pm(f['rho_solar'], 2)} {RHOSUN}, is {e['rho_verdict']} with the catalogue "
         f"value for this F star ({e['rho_tic']:.2f} {RHOSUN})."),
        ("pass", "Depths agree.",
         "The four transits have the same depth within the noise; SPOC's odd-versus-even test found the same."),
        ("pass", "The signal is on this star.",
         "The dip is present in the summed flux of all pixels (680, 454, 665 and 623 ppm in the four sectors). "
         "The only Gaia source within 2′ bright enough to mimic it is TIC 128703032 (T = 10.6, 80″ away), which "
         "would need eclipses about 7% deep. Its brightest pixel changes by less than 0.25% at the transit times."),
        ("pass", "No sign of a stellar companion.",
         "Gaia DR3: RUWE 0.94, no non-single-star solution, radial velocity stable to 0.13 km/s over 15 epochs."),
        ("warn", "SPOC's centroid test flagged an offset.",
         "The 2023 DV report gives a 32″ offset of the difference-image centroid, but rates 0 of its 3 difference "
         "images as good quality: the star saturates the detector and bleeds along columns, which breaks "
         "centroiding. The pixel test above does not rely on centroids."),
        ("info", "Not yet ruled out:",
         "an eclipsing companion closer to the star than Gaia can separate from a G = 5.9 star. Eisner et al. "
         "(2021) report Gemini speckle imaging with no nearby companion and spectra consistent with a planet; "
         "those radial velocities, phased on the 23.88-day period, would also bound the planet's mass."),
    ])
    return f"""
<section id="hip71639" class="cand">
  <div class="cand-head">
    <p class="kicker">Independent confirmation · Track B (long-period transits) · SPOC 2-minute data</p>
    <h2>HIP 71639: a warm Neptune every 23.9 days</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("Community TOI since 2020", "known")}{chip("Period posted May 2026", "known")}{chip("SPOC TCE with a 31× wrong period", "tce")}</div>
    <p class="lede col">An F star at the edge of naked-eye visibility (V ≈ 6.1), 34 parsecs away in Lupus, dims by
    {depth:.0f} parts per million for 7.2 hours every {f['period_d'][0]:.4f} days. A planet about {rp:.1f} times
    the size of Earth would do this.</p>
  </div>
  <div class="cand-grid">
    <div>
      <figure class="plate" style="margin:0">
        <div class="scroll-x">{svg}</div>
        <figcaption>TESS 2-minute photometry (PDCSAP) in 20-minute bins, each window normalised by a quadratic
        fitted outside the transit. The orange line is one transit model fitted to all four windows together.
        Hover a point for its value.</figcaption>
      </figure>
      <div class="col" style="display:grid;gap:12px">
        <h3>Who saw it before</h3>
        <p>Volunteers of the citizen-science project Planet Hunters TESS spotted the sector-11 transit. Eisner et
        al. (2021, <a href="https://arxiv.org/abs/2011.13944">MNRAS 501, 4669</a>) list the star as a single-transit
        planet candidate, and the project lists it among the community TOIs it sent to ExoFOP. The paper reports
        43 MINERVA-Australis and 2 LCO/NRES spectra consistent with a planet and Gemini speckle imaging showing no
        nearby companion; the ground-based data were left for a later paper. CHIRON spectra followed in 2023.</p>
        <p>On 15 May 2026 another ExoFOP user posted the period, P = 23.8843518 ± 0.0000098 d, from the same four
        transits in sectors 11, 38, 65 and 102. This survey found it independently in October and agrees within
        0.1σ ({f['period_d'][0]:.7f} ± {max(f['period_d'][1:]):.7f} d). The planet is still not a TOI.</p>
        <p>NASA's SPOC pipeline also detected the signal. Its multi-sector searches of sectors 1–39 and then 1–65
        found the transits in sectors 11, 38 and 65 and reported a threshold-crossing event (multiple-event
        statistic 26.6). One transit fell in each sector, two years apart, so the pipeline took the gap between
        them as the period: <span class="num">740.42 d</span>. At that period the transit is far too short for the
        star (its fit implied a stellar density of 8.6 ± 6.5 {RHOSUN}), and the centroid test was flagged. In the
        TOI lists available to us it never became a TOI.</p>
        <p>The fourth transit, in sector 102 (April 2026), arrives {abs(e['s102_offset_min']):.0f} minutes from
        where 106 cycles of a {f['period_d'][0]:.3f}-day period put it, and 740.42 days is exactly 31 of those
        cycles. What this page adds is the vetting below and the transit fit.</p>
      </div>
      <div class="col" style="display:grid;gap:12px">
        <h3>Vetting</h3>
        {checks}
      </div>
    </div>
    <div>
      <dl class="spec">
        <dt>Star</dt><dd>HIP 71639 · HD 128582 · TIC 128703021</dd>
        <dt>Position</dt><dd>14h39m10.8s −46°35′04″</dd>
        <dt>Brightness</dt><dd>T 5.59 · G 5.94</dd>
        <dt>Host</dt><dd>F7 · 6280 K · 1.60 {RSUN} · 1.23 {MSUN}</dd>
        <dt>Distance</dt><dd>34.0 pc</dd>
        <dt>Period</dt><dd>{f['period_d'][0]:.5f} ± {max(f['period_d'][1:]):.5f} d</dd>
        <dt>Mid-transit</dt><dd>BTJD {f['tref_btjd'][0]:.4f} ± {max(f['tref_btjd'][1:]):.4f}</dd>
        <dt>Depth</dt><dd>{depth:.0f} ppm</dd>
        <dt>Duration</dt><dd>{pm(f['t14_h'], 2)} h</dd>
        <dt>Impact b</dt><dd>{pm(f['b'], 2)}</dd>
        <dt>Rp / R★</dt><dd>{pm(f['rp_rs'], 4)}</dd>
        <dt>Radius</dt><dd>{rp:.1f} ± {rp_err:.1f} {RE}</dd>
        <dt>Orbit</dt><dd>{e['a_au']:.3f} AU · {e['s_earth']:.0f}× Earth's insolation</dd>
        <dt>T<sub>eq</sub></dt><dd>≈ {e['teq']:.0f} K</dd>
        <dt>RV signal</dt><dd>≈ {e['k_ms']:.1f} m/s for {e['mass_assumed']:.0f} M<sub>⊕</sub></dd>
        <dt>Transits</dt><dd>S11 · S38 · S65 · S102</dd>
      </dl>
      <div class="plate">
        <p class="label">Next transits seen from Chile</p>
        <ul class="nexts">{nxt}</ul>
        <p class="muted" style="font-size:var(--step--1)">Fraction of the 7.2-hour transit at night from La Silla
        (Sun below −12°, star above 30°). Timing uncertainty is under 2 minutes.</p>
      </div>
    </div>
  </div>
</section>"""


def section_hip56966(fit, ex):
    if "HIP56966" not in fit or "HIP56966" not in ex:
        return ""
    f = fit["HIP56966"]
    e = ex["HIP56966"]
    rs, rs_err = 1.069, 0.05
    rp = f["rp_rs"][0] * rs * R_EARTH_PER_RSUN
    rp_err = rp * math.hypot(max(f["rp_rs"][1], f["rp_rs"][2]) / f["rp_rs"][0], rs_err / rs)
    depth = -min(min(p["my"]) for p in f["plot"])
    labels = [f"Sector {s} · {btjd_to_utc(p['tc']).strftime('%b %Y')}" for s, p in zip((10, 99, 100), f["plot"])]
    svg = svg_lightcurves(f["plot"], labels, id_="lc-hip56966", panel_w=330, xr=10, xstep=4)
    nxt = "".join(f"<li><span class='num'>{btjd_to_utc(r['btjd']).strftime('%-d %b %Y, %H:%M')} UTC</span>"
                  f"<span class='muted'>{esc(r['note'])}</span></li>" for r in e["next"])
    checks = check_rows([
        ("pass", "Four transits, one period.",
         f"Dips in sectors 10, 64, 99 and 100 (710, 651, 694 and 799 ppm, each ±60) fall on one ephemeris with "
         f"4 minutes rms: 74 and 51 cycles separate the earlier ones from the sector-99 transit. Every shorter alias "
         f"predicts transits inside the data that are missing, ruled out at more than {e['alias_min_sigma']:.0f}σ. "
         f"The only other predicted transit since 2019 fell in the mid-sector gap of sector 37."),
        ("pass", "No secondary eclipse.",
         "Half an orbit after each transit the light curves are flat to about 150 ppm, as for a planet and unlike most "
         "eclipsing binaries."),
        ("pass", "The signal is on this star.",
         "A source test on the 2-minute pixels of three transits: the difference image is fitted as the star's "
         "image placed at each candidate position. The target fits best. Every Gaia source within 150″ bright "
         "enough to cause the dip with a 50% eclipse fits worse by Δχ² ≥ 25 summed over the three; the nearest "
         "plausible one (G = 13.9, 37″ away, which would need 10%-deep eclipses) by Δχ² = 79."),
        ("pass", "No sign of a stellar companion.",
         "Gaia DR3: RUWE 0.81, no non-single-star solution, radial velocity stable to 0.18 km/s over 12 epochs."),
        ("warn", "Slightly short transit.",
         f"A circular orbit at this period gives a {pm(f['t14_h'], 1)}-hour transit if the planet crosses at impact "
         f"parameter ≈ {e['b_circ']:.2f}; the shape prefers {pm(f['b'], 2)}. Consistent within 1σ, or a modest "
         f"eccentricity."),
        ("info", "Not yet ruled out:", "an eclipsing companion within about 2″, which imaging and radial velocities "
         "would exclude."),
        ("info", "On record:",
         "Salinas et al. (2025) posted the sector-10 transit to ExoFOP as a single-transit community candidate "
         "(778 ppm, 2.2 h, no period). Gemini speckle imaging in March 2020 found no companion brighter than "
         "Δm = 6.9 at 0.5″. The other three transits and the period appear nowhere else."),
    ])
    return f"""
<section id="hip56966" class="cand">
  <div class="cand-head">
    <p class="kicker">Candidate 1 · Track B, found during vetting · SPOC 2-minute data</p>
    <h2>HIP 56966: a sub-Neptune every 19.8 days</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("Period new", "new")}{chip("Community TOI since 2025, no period", "known")}</div>
    <p class="lede col">A Sun-like G0 star (V = 8.8), 73 parsecs away in Centaurus, dims by {depth:.0f} ppm for
    2.8 hours every {f['period_d'][0]:.4f} days. A planet about {rp:.1f} times the size of Earth would do this.
    SPOC never saw two of its transits in one search, so it raised no alert.</p>
  </div>
  <div class="cand-grid">
    <div>
      <figure class="plate" style="margin:0">
        <div class="scroll-x">{svg}</div>
        <figcaption>TESS 2-minute photometry around the transits in sectors 10, 99 and 100 (20-minute bins), with
        one model fitted to all three. The sector-64 transit, seen only in full-frame images, agrees in depth.</figcaption>
      </figure>
      <div class="col" style="display:grid;gap:12px">
        <h3>Vetting</h3>
        {checks}
      </div>
    </div>
    <div>
      <dl class="spec">
        <dt>Star</dt><dd>HIP 56966 · TIC 266376971</dd>
        <dt>Position</dt><dd>11h40m40.7s −58°06′16″</dd>
        <dt>Brightness</dt><dd>V 8.84 · T 8.28</dd>
        <dt>Host</dt><dd>G0 · 5970 K · 1.07 {RSUN} · 1.09 {MSUN}</dd>
        <dt>Distance</dt><dd>72.6 pc</dd>
        <dt>Period</dt><dd>{f['period_d'][0]:.5f} ± {max(f['period_d'][1:]):.5f} d</dd>
        <dt>Mid-transit</dt><dd>BTJD {f['tref_btjd'][0]:.4f} ± {max(f['tref_btjd'][1:]):.4f}</dd>
        <dt>Depth</dt><dd>{depth:.0f} ppm</dd>
        <dt>Duration</dt><dd>{pm(f['t14_h'], 2)} h</dd>
        <dt>Rp / R★</dt><dd>{pm(f['rp_rs'], 4)}</dd>
        <dt>Radius</dt><dd>{rp:.1f} ± {rp_err:.1f} {RE}</dd>
        <dt>Orbit</dt><dd>{e['a_au']:.3f} AU · {e['s_earth']:.0f}× Earth's insolation</dd>
        <dt>T<sub>eq</sub></dt><dd>≈ {e['teq']:.0f} K</dd>
        <dt>RV signal</dt><dd>≈ {e['k_ms']:.1f} m/s for {e['mass_assumed']:.0f} M<sub>⊕</sub></dd>
        <dt>Transits</dt><dd>S10 · S64 · S99 · S100</dd>
      </dl>
      <div class="plate">
        <p class="label">Next transits seen from Chile</p>
        <ul class="nexts">{nxt}</ul>
      </div>
    </div>
  </div>
</section>"""


def _newcands():
    p = CAND / "new_candidates.json"
    return json.load(open(p)) if p.exists() else {}


def section_tic402710280(fit, ex):
    key = "TIC402710280"
    nc = _newcands().get("402710280")
    if key not in fit or key not in ex or not nc:
        return ""
    f, e = fit[key], ex[key]
    rs, rs_err = 0.817, 0.05
    rp = f["rp_rs"][0] * rs * R_EARTH_PER_RSUN
    rp_err = rp * math.hypot(max(f["rp_rs"][1], f["rp_rs"][2]) / f["rp_rs"][0], rs_err / rs)
    depth = -min(min(p["my"]) for p in f["plot"][:4])
    labels = [f"S{s} · {btjd_to_utc(p['tc']).strftime('%b %Y')}" for s, p in zip(f["sectors"], f["plot"])]
    svg = svg_lightcurves(f["plot"], labels, id_="lc-tic402710280", panel_w=250, xr=12, xstep=6)
    tr = nc["transits"]
    c = nc["centroid"]
    near = c["closest_neighbours"]
    dens = nc["density_test"]["dchi2_rho_fixed"]
    nxt = "".join(f"<li><span class='num'>{btjd_to_utc(r['btjd']).strftime('%-d %b %Y, %H:%M')} UTC</span>"
                  f"<span class='muted'>{esc(r['note'])}</span></li>" for r in e["next"])
    checks = check_rows([
        ("pass", "Five transits, one period, none missing.",
         "Dips in sectors 39, 66, 93, 100 and 102 (" + ", ".join(f"{t['depth_ppm']:.0f}" for t in tr) +
         " ppm; signal-to-noise " + ", ".join(f"{t['snr']:.0f}" for t in tr) + ") fall on one ephemeris. "
         "Every shorter alias predicts transits inside the data that are not there. Every transit this period "
         "predicts inside the data is there: sectors 12 and 101 contain none, and the one due in sector 103 falls "
         "in a data gap."),
        ("pass", "The signal is on this star.",
         f"Difference-image centroids of the five transits, measured from the star's own image on each frame: "
         f"{c['offset_arcsec']:.1f}″ ± {c['centroid']['err']:.1f}″ from the target. Every Gaia star able to fake "
         f"the dip lies {near[1]['sigma']:.0f}σ or more from that position, except one: a G = 16.2 background star "
         f"2″ away (parallax 0.69 mas against the target's 15.0)."),
        ("pass", "Not that background star.",
         "To make a 0.2% dip through 6.4 magnitudes of dilution it would need a 73%-deep eclipse, which takes a "
         "companion nearly as large as the star and a long, V-shaped ingress. Fitting the transits as a diluted "
         "eclipse of that kind is worse than a planet by Δχ² = 79 (radius ratio ≥ 0.7), and still by 65 at 0.5 and "
         "27 at 0.3."),
        ("pass", "A transit this star can make.",
         f"Without the sector-102 transit the fitted stellar density is {pm(f['rho_solar'], 2)} solar; the "
         f"catalogue K dwarf's 1.61 costs Δχ² = {dens['without_S102']:.1f}. The sector-102 transit, 1.4 days from "
         f"the start of an orbit where our light curve has a ramp, alone pulls the density down (Δχ² = "
         f"{dens['all_five']:.0f} with it); it is shown but left out of the shape fit."),
        ("pass", "No secondary eclipse, even and odd alike.",
         f"Half an orbit later the light curve is flat (37 to 120 ± 48 ppm); even and odd transits have depths "
         f"{nc['odd_even_ppm']['even'][0]} ± {nc['odd_even_ppm']['even'][1]} and {nc['odd_even_ppm']['odd'][0]} ± "
         f"{nc['odd_even_ppm']['odd'][1]} ppm (2.3σ)."),
        ("pass", "No sign of a stellar companion.", "Gaia DR3: RUWE 1.00, no non-single-star solution, a single "
         "image in every scan."),
        ("info", "Two transits only in new data.", "The sector 100 and 102 transits come from our own photometry of "
         "the full-frame images; no public light curve of those sectors exists yet (the public QLP light curves end "
         "at sector 98)."),
        ("info", "On record:", "nothing. Not a TOI or community TOI, no planet parameters, imaging or spectra on "
         "ExoFOP (checked live on 12 October 2026); not in SIMBAD's object types, the NASA archive or the literature "
         "searches. SPOC never observed it at 2-minute cadence, so it has no SPOC TCE."),
    ])
    return f"""
<section id="tic402710280" class="cand">
  <div class="cand-head">
    <p class="kicker">New candidate · period hunt, sectors 99–108 full-frame images</p>
    <h2>TIC 402710280: a sub-Neptune every 36.3 days</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("Not a community TOI", "new")}{chip("Period new", "new")}{chip("Nothing on record", "new")}</div>
    <p class="lede col">A K dwarf (T = 9.2, V ≈ 10) 67 parsecs away in Triangulum Australe dims by {depth:.0f} ppm
    for {f['t14_h'][0]:.1f} hours every {f['period_d'][0]:.4f} days: a planet about {rp:.1f} times the size of
    Earth. The survey first saw two of its transits, 1,454 days apart; our own light curves of the newest
    full-frame images (sectors 100 and 102) added two more, and the period they fix finds a fifth in sector 66.</p>
  </div>
  <div class="cand-grid">
    <div>
      <figure class="plate" style="margin:0">
        <div class="scroll-x">{svg}</div>
        <figcaption>Every transit (20-minute bins) with one model fitted to the first four: TESS-SPOC (sectors 39
        and 66), QLP (93) and our own full-frame photometry (100 and 102). The sector-102 window sits on a ramp at
        the start of an orbit.</figcaption>
      </figure>
      <div class="col" style="display:grid;gap:12px">
        <h3>Vetting</h3>
        {checks}
      </div>
    </div>
    <div>
      <dl class="spec">
        <dt>Star</dt><dd>TIC 402710280 · TYC 9041-03677-1</dd>
        <dt>Position</dt><dd>16h34m05.7s −63°32′59″</dd>
        <dt>Brightness</dt><dd>G 9.78 · T 9.24</dd>
        <dt>Host</dt><dd>K dwarf · 5180 K · 0.82 {RSUN} · 0.88 {MSUN}</dd>
        <dt>Distance</dt><dd>66.8 pc</dd>
        <dt>Period</dt><dd>{f['period_d'][0]:.5f} ± {max(f['period_d'][1:]):.5f} d</dd>
        <dt>Mid-transit</dt><dd>BTJD {f['tref_btjd'][0]:.4f} ± {max(f['tref_btjd'][1:]):.4f}</dd>
        <dt>Depth</dt><dd>{depth:.0f} ppm</dd>
        <dt>Duration</dt><dd>{pm(f['t14_h'], 2)} h</dd>
        <dt>Rp / R★</dt><dd>{pm(f['rp_rs'], 4)}</dd>
        <dt>Radius</dt><dd>{rp:.1f} ± {rp_err:.1f} {RE}</dd>
        <dt>Orbit</dt><dd>{e['a_au']:.3f} AU · {e['s_earth']:.0f}× Earth's insolation</dd>
        <dt>T<sub>eq</sub></dt><dd>≈ {e['teq']:.0f} K</dd>
        <dt>RV signal</dt><dd>≈ {e['k_ms']:.1f} m/s for {e['mass_assumed']:.0f} M<sub>⊕</sub></dd>
        <dt>Transits</dt><dd>S39 · S66 · S93 · S100 · S102</dd>
        <dt>Next in TESS</dt><dd>sectors 113–114</dd>
      </dl>
      <div class="plate">
        <p class="label">Next transits seen from Chile</p>
        <ul class="nexts">{nxt}</ul>
      </div>
    </div>
  </div>
</section>"""


def section_toi2493(fit, ex):
    key = "TOI2493"
    nc = _newcands().get("123664207")
    if key not in fit or key not in ex or not nc:
        return ""
    f, e = fit[key], ex[key]
    rs, rs_err = 0.652, 0.059
    rp = f["rp_rs"][0] * rs * R_EARTH_PER_RSUN
    rp_err = rp * math.hypot(max(f["rp_rs"][1], f["rp_rs"][2]) / f["rp_rs"][0], rs_err / rs)
    depth = -min(min(p["my"]) for p in f["plot"])
    labels = [f"S{s} · {btjd_to_utc(p['tc']).strftime('%b %Y')}" for s, p in zip(f["sectors"], f["plot"])]
    svg = svg_lightcurves(f["plot"], labels, id_="lc-toi2493", panel_w=300, xr=12, xstep=6)
    tr = nc["transits"]
    c = nc["centroid"]
    near = c["closest_neighbours"]
    nxt = "".join(f"<li><span class='num'>{btjd_to_utc(r['btjd']).strftime('%-d %b %Y, %H:%M')} UTC</span>"
                  f"<span class='muted'>{esc(r['note'])}</span></li>" for r in e["next"])
    checks = check_rows([
        ("pass", "Three transits, one period.",
         "The TOI's transit (sector 33, 2021) and two more found by the period hunt in sectors 98 and 107 (" +
         ", ".join(f"{t['depth_ppm']:.0f}" for t in tr) + " ppm; signal-to-noise " + ", ".join(f"{t['snr']:.1f}" for t in tr) +
         ") share one period; every shorter alias predicts transits inside sectors 6, 33, 87, 98 or 107 that are not there."),
        ("pass", "The signal is on this star.",
         f"Difference-image centroids of the three transits (2-minute pixels): {c['offset_arcsec']:.1f}″ ± "
         f"{c['centroid']['err']:.1f}″ from the target; every Gaia star able to fake the dip is {near[0]['sigma']:.0f}σ "
         f"or more away."),
        ("pass", "A transit this star can make.",
         f"The fitted stellar density, {pm(f['rho_solar'], 2)} solar, matches the K dwarf's 2.4; fixing it changes "
         f"χ² by less than 0.1."),
        ("warn", "Modest signal.", "Each transit alone is a 7–9σ event, the sector-98 one the weakest. TESS should "
         "see the next transit in sector 110, and a 1-metre telescope can catch the 1.5-ppt, 6-hour dip from the ground."),
        ("info", "On record:", "TOI-2493.01 since February 2021, as a single transit (\"single at TBJD~2226\"), "
         "period 0. SPOC's only detection is a 15.1-day TCE in sector 33 that the transits in sectors 98 and 107 "
         "rule out. No period on ExoFOP or in the NASA archive, and no paper found."),
    ])
    return f"""
<section id="toi2493" class="cand">
  <div class="cand-head">
    <p class="kicker">New period for a TOI · period hunt, SPOC 2-minute data</p>
    <h2>TOI-2493.01: a temperate sub-Neptune every 79 days</h2>
    <div class="chips">{chip("TOI since 2021, no period", "known")}{chip("Period new", "new")}</div>
    <p class="lede col">A K5 dwarf (T = 11.4) 93 parsecs away dims by {depth:.0f} ppm for {f['t14_h'][0]:.1f} hours
    every {f['period_d'][0]:.3f} days. The planet would be about {rp:.1f} times the size of Earth and receive
    {e['s_earth']:.1f} times Earth's insolation, on the inner edge of the star's habitable zone.</p>
  </div>
  <div class="cand-grid">
    <div>
      <figure class="plate" style="margin:0">
        <div class="scroll-x">{svg}</div>
        <figcaption>SPOC 2-minute photometry of the three transits (20-minute bins) with one model fitted to
        all three.</figcaption>
      </figure>
      <div class="col" style="display:grid;gap:12px">
        <h3>Vetting</h3>
        {checks}
      </div>
    </div>
    <div>
      <dl class="spec">
        <dt>Star</dt><dd>TOI-2493 · TIC 123664207</dd>
        <dt>Position</dt><dd>05h58m59.1s −21°34′30″</dd>
        <dt>Brightness</dt><dd>T 11.43</dd>
        <dt>Host</dt><dd>K5 dwarf · 4240 K · 0.65 {RSUN} · 0.66 {MSUN}</dd>
        <dt>Distance</dt><dd>93.1 pc</dd>
        <dt>Period</dt><dd>{f['period_d'][0]:.4f} ± {max(f['period_d'][1:]):.4f} d</dd>
        <dt>Mid-transit</dt><dd>BTJD {f['tref_btjd'][0]:.4f} ± {max(f['tref_btjd'][1:]):.4f}</dd>
        <dt>Depth</dt><dd>{depth:.0f} ppm</dd>
        <dt>Duration</dt><dd>{pm(f['t14_h'], 2)} h</dd>
        <dt>Radius</dt><dd>{rp:.1f} ± {rp_err:.1f} {RE}</dd>
        <dt>Orbit</dt><dd>{e['a_au']:.3f} AU · {e['s_earth']:.1f}× Earth's insolation</dd>
        <dt>T<sub>eq</sub></dt><dd>≈ {e['teq']:.0f} K</dd>
        <dt>Transits</dt><dd>S33 · S98 · S107</dd>
        <dt>Next in TESS</dt><dd>sector 110 (15 Nov 2026)</dd>
      </dl>
      <div class="plate">
        <p class="label">Next transits seen from Chile</p>
        <ul class="nexts">{nxt}</ul>
      </div>
    </div>
  </div>
</section>"""


def section_hd11506(fit, ex):
    f = fit["HD11506"]
    e = ex["HD11506"]
    rs, rs_err = 1.350, 0.06
    rp = f["rp_rs"][0] * rs * R_EARTH_PER_RSUN
    rp_err = rp * math.hypot(max(f["rp_rs"][1], f["rp_rs"][2]) / f["rp_rs"][0], rs_err / rs)
    labels = ["Sector 30 · Oct 2020", "Sector 107 · Aug 2026"]
    svg = svg_lightcurves(f["plot"], labels, id_="lc-hd11506", panel_w=330, xr=10, xstep=4)
    per = f["per_transit"]
    depth = -min(min(p["my"]) for p in f["plot"])
    d_i = [depth * (x["rp_rs"] / f["rp_rs"][0]) ** 2 for x in per]
    al = e["aliases"]
    allowed = [a for a in al if a["status"] == "allowed"]
    alias_rows = "".join(
        f"<tr><td class='r m'>{a['n']}</td><td class='r m'>{a['P']:.2f}</td><td class='r m'>{a['e_min_b50']:.2f}</td>"
        f"<td class='r m'>{a['b_circ']:.2f}</td></tr>" for a in allowed[::-1][:12])
    checks = check_rows([
        ("pass", "Two transits, same shape.",
         f"Durations agree; depths {d_i[0]:.0f} and {d_i[1]:.0f} ppm, "
         f"which differ by about {e['depth_diff_sigma']:.1f}σ."),
        ("pass", "Seen in two independent data sets.",
         "Sector 30 in both the 2-minute and the full-frame (QLP) light curves; sector 107 in 2-minute data. "
         "Both are in the raw (SAP) flux as well as the corrected flux."),
        ("pass", "On target.",
         "Difference-image centroids shift by 0.27 and 0.39 pixels, consistent with the star. No other Gaia source "
         "within 88″ is bright enough to cause a 700 ppm dip, even if fully eclipsed."),
        ("pass", "Not one of the known planets.",
         f"The three radial-velocity planets have periods of 223.9, 1618 and ~26,500 days. Neither of the first "
         f"two fits the {e['dt']:.2f}-day gap between the transits."),
        ("warn", "The transit is short for this star.",
         f"At the shortest allowed period ({allowed[-1]['P']:.1f} d) a circular orbit gives a 3.2-hour transit "
         f"only if the planet crosses near the stellar limb (impact parameter ≈ {allowed[-1]['b_circ']:.2f}), "
         f"while the shape prefers a more central crossing (b = {pm(f['b'], 2)}). An eccentric orbit "
         f"(e ≳ {allowed[-1]['e_min_b50']:.1f}, transit near periastron) also fits. Longer periods need more "
         "eccentricity."),
        ("info", "On record:", "the TESS single-transit working group posted the sector-30 transit to ExoFOP in "
         "August 2024 as a community TOI, but with a depth of 7,210 ± 3,331 ppm and a duration of 26.9 ± 7.4 hours, "
         "ten and eight times what the 2-minute data show. Nothing on record mentions the second transit or the "
         "period constraints."),
        ("info", "Not yet ruled out:", "an eclipsing companion closer than ~1″, or two different planets each "
         "transiting once."),
    ])
    return f"""
<section id="hd11506" class="cand">
  <div class="cand-head">
    <p class="kicker">Candidate 2 · Track B (duo transit) · SPOC 2-minute data</p>
    <h2>HD 11506: a fourth planet in a known planetary system?</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("Second transit new", "new")}{chip("Community TOI since 2024", "known")}{chip("Host has 3 RV planets", "known")}</div>
    <p class="lede col">HD 11506 is a bright Sun-like star in Cetus with three giant planets found by the
    radial-velocity method, none known to transit. TESS caught two matching 3.2-hour dips, almost six years apart,
    of the size a planet {rp:.1f} times the size of Earth would make.</p>
  </div>
  <div class="cand-grid">
    <div>
      <figure class="plate" style="margin:0">
        <div class="scroll-x">{svg}</div>
        <figcaption>TESS 2-minute photometry around the two transits (20-minute bins), with one transit shape
        fitted to both. The period is unknown, so only the shape and the two mid-times are fitted.</figcaption>
      </figure>
      <div class="col" style="display:grid;gap:12px">
        <h3>What the period could be</h3>
        <p>The transits are {e['dt']:.2f} days apart, so the period is that gap divided by a whole number. Each
        choice predicts other transits; where TESS was watching at those times (sectors 3, 30 and 107) and saw
        nothing, that choice is ruled out. {len(allowed)} of {len(al)} choices survive, the shortest being
        {allowed[-1]['P']:.1f} days. The table lists the shortest ones with the impact parameter a circular
        orbit would need to give a 3.2-hour transit, and the least eccentricity that works with the fitted shape.</p>
        <div class="tablewrap"><table><thead><tr><th class="r">n</th><th class="r">Period (d)</th>
        <th class="r">min. e</th><th class="r">b if circular</th></tr></thead><tbody>{alias_rows}</tbody></table></div>
      </div>
      <div class="col" style="display:grid;gap:12px">
        <h3>Vetting</h3>
        {checks}
      </div>
    </div>
    <div>
      <dl class="spec">
        <dt>Star</dt><dd>HD 11506 · HIP 8770 · TIC 266680498</dd>
        <dt>Position</dt><dd>01h52m50.5s −19°30′25″</dd>
        <dt>Brightness</dt><dd>T 6.98 · G 7.39</dd>
        <dt>Host</dt><dd>G0 · 6030 K · 1.35 {RSUN} · 1.12 {MSUN}</dd>
        <dt>Distance</dt><dd>51.3 pc</dd>
        <dt>Transit 1</dt><dd>BTJD {f['tc_0'][0]:.4f} ± {max(f['tc_0'][1:]):.4f}</dd>
        <dt>Transit 2</dt><dd>BTJD {f['tc_1'][0]:.4f} ± {max(f['tc_1'][1:]):.4f}</dd>
        <dt>Depth</dt><dd>{depth:.0f} ppm</dd>
        <dt>Duration</dt><dd>{pm(f['t14_h'], 2)} h</dd>
        <dt>Rp / R★</dt><dd>{pm(f['rp_rs'], 4)}</dd>
        <dt>Radius</dt><dd>{rp:.1f} ± {rp_err:.1f} {RE}</dd>
        <dt>Period</dt><dd>{e['dt']:.2f} d / n, ≥ {allowed[-1]['P']:.1f} d</dd>
        <dt>Gaia DR3</dt><dd>RUWE 0.92 · RV σ 0.13 km/s</dd>
      </dl>
    </div>
  </div>
</section>"""


def trackb_shortlist(tb):
    """Single/duo candidates worth listing: not TOIs, planet-sized (<= 16 R_E), long period, not a much-flagged SPOC signal."""
    x = tb[tb.verdict.isin(["single", "repeating"])].copy()
    if "is_toi" in x:
        x = x[~x.is_toi.astype(bool)]
    longp = (x.verdict == "single") | x.period_guess.isna() | (x.period_guess >= 15)
    return x[(x.rp_rearth <= 16) & (x.n_spoc_tce.fillna(0) < 5) & longp & ~x.tic.isin([128703021, 266680498, 266376971])]


def verdict_chip(v):
    m = {"single": ("single transit", "new"), "repeating": ("repeats", "new"), "duo": ("duo", "new"),
         "EB (deeper eclipses)": ("eclipsing binary", "eb"), "TOI on neighbour": ("TOI nearby", "known"),
         "anomaly candidate": ("anomaly", "new"), "EB (deeper dips)": ("eclipsing binary", "eb"),
         "periodic (EB/planet)": ("periodic", "known")}
    t, k = m.get(v, (v, "known"))
    return chip(t, k)


# Track B stars that web searches found already confirmed (names as published); the NASA confirmed-planet
# cross-match (known.flag_confirmed) adds the rest.
CONFIRMED_LIT = {65910228: "NGTS-38 b", 279401253: "TIC 279401253 b", 4672985: "TIC 4672985 b",
                 137966491: "HD 278555 b", 318707987: "HD 3167 c", 221567884: "HD 148797 b, c",
                 453147896: "NGTS-39 b", 65612701: "TIC 65612701 b", 70803960: "K2-290 c"}
# Confirmed planets without TOIs that Track B found blind (checked against the published periods), and one it lost.
RECOVERED = [("NGTS-38 b", "180.5 d, one transit"), ("TIC 279401253 b", "76.8 d, one transit"),
             ("TIC 4672985 b", "69.0 d, one transit"), ("NGTS-39 b", "58.2 d"), ("K2-290 c", "48.4 d, one transit"),
             ("HD 148797 b and c", "42.1 and 68.2 d"), ("HD 3167 c", "29.85 d, one transit"),
             ("HD 278555 b", "28.84 d, period recovered"), ("TIC 65612701 b", "10.22 d, period recovered")]
MISSED = [("TIC 245076932 b", "21.6 d: all nine transits found, then rejected as an eclipsing binary")]


def period_hunt_summary():
    """Counts and the table of new periods from the period hunt and its review (period_hunt*.json, period_review.json)."""
    rv = _json("period_review.json")
    if not rv:
        return None
    rv = rv if isinstance(rv, list) else list(rv.values())
    hunt = []
    for name in ("period_hunt.json", "period_hunt_ffi.json"):
        x = _json(name)
        hunt += x if isinstance(x, list) else []
    best = {}
    for r in hunt:
        k = (r["tic"], r["id"])
        if k not in best or r.get("ffi_sectors"):
            best[k] = r
    targets = [r for r in best.values() if not str(r["id"]).startswith("single")]
    n_t = len(targets)
    n_data = sum(1 for r in targets if r.get("status") not in ("no data", "error"))
    n_rec = sum(1 for r in targets if r.get("recovered"))
    n_ffi = sum(1 for r in targets if r.get("ffi_sectors"))
    rv = [o for o in rv if not str(o.get("id", "")).startswith("single") and "error" not in o]
    cred = [o for o in rv if o.get("n_credible", 0) >= 1]
    def known_P(o):
        ps = []
        for txt in (o.get("nasa") or "", o.get("exofop_periods") or ""):
            ps += [float(v) for v in re.findall(r"P=([0-9.]+)", txt)]
        return ps
    def agrees(o, ps):
        al = o.get("aliases") or []
        return any(abs(a / p - 1) < 2e-3 for a in al for p in ps)
    on_rec = [o for o in cred if known_P(o)]
    n_agree = sum(1 for o in on_rec if agrees(o, known_P(o)))
    new = [o for o in cred if not known_P(o) and o.get("n_aliases") and o["n_aliases"] <= 3]
    return {"n_targets": n_t, "n_data": n_data, "n_rec": n_rec, "n_ffi": n_ffi, "n_review": len(rv), "n_cred": len(cred),
            "n_onrec": len(on_rec), "n_agree": n_agree, "new": new,
            "n_unique": sum(1 for o in cred if o.get("n_aliases") == 1)}


def section_periodhunt(summary_only=False):
    S = period_hunt_summary()
    if S is None:
        return ""
    if summary_only:
        return (f'<li>{chip(f"{len(S["new"])} periods", "new")}<div><strong>Periods for TOIs and community TOIs that '
                f'have none.</strong> All {S["n_targets"]} were searched in every TESS sector, including our own light curves '
                f'of sectors 99–108; {S["n_cred"]} gained transits that pass the checks. Where a period was already on '
                f'record, the hunt agrees in {S["n_agree"]} of {S["n_onrec"]} cases; {len(S["new"])} periods, narrowed to '
                f'three or fewer choices, are on no list (TOI-2493.01 is the best of them).</div></li>')
    pv = _json("period_vet.json")
    rows = []
    for o in S["new"]:
        t = int(o["tic"])
        ev = [e for e in o.get("events", []) if e.get("ok")]
        v = pv.get(f"{t}|{o['id']}", {})
        rows.append(dict(id=o["id"], tic=t, disp=o.get("disp", ""), depth=o.get("depth_ref_ppm"), dur=o.get("dur_h"),
                         n=1 + len(ev), sectors=", ".join(str(e["sector"]) for e in ev), al=o.get("aliases") or [],
                         pmin=o.get("p_min_circular"), pix=cv_text(v) if v else "not tested",
                         verdict=v.get("verdict", "")))
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("n", ascending=False)
    cols = [("TOI / CTOI", lambda r: f"<span class='m'>{esc(r.id)}</span><div class='muted' style='font-size:0.78rem'>TIC {r.tic}</div>", "l"),
            ("Depth", lambda r: f"{r.depth:.0f} ppm" if r.depth == r.depth else "–", "r"),
            ("Transits", lambda r: f"{r.n} <span class='muted'>(new: S{esc(r.sectors)})</span>", "l"),
            ("Periods allowed (d)", lambda r: esc(", ".join(f"{a:.3f}" for a in r.al[:3])), "l"),
            ("Shortest circular", lambda r: f"{r.pmin:.0f} d" if r.pmin == r.pmin and r.pmin else "–", "r"),
            ("Pixels", lambda r: f"<span class='muted' style='font-size:0.8rem'>{esc(r.pix)}</span>", "l")]
    return f"""
<section id="periodhunt">
  <p class="kicker">Period hunt · TOIs and community TOIs without a period</p>
  <h2>Periods for single-transit candidates, from every sector</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>ExoFOP lists {S['n_targets']} TOIs and community TOIs (not false positives) with no period, most of them single
  transits. For each, every sector of TESS data was loaded (SPOC 2-minute, TESS-SPOC and QLP light curves, plus our
  own light curves of sectors 99–108 for {S['n_ffi']} of them), the catalogued transit was re-found ({S['n_rec']} of
  {S['n_data']} with data), other planets on the star were masked, and the rest of the data was searched for dips of
  the same duration and similar depth. Each new dip was then re-measured on its own sector's light curve against a
  local baseline and kept only if it is significant there, fully covered, away from sector edges and gaps, free of
  background jumps, present in the uncorrected (SAP) flux, and of the right depth. The periods that put all kept
  transits on one ephemeris without predicting transits the data do not show are listed.</p>
  <p>{S['n_cred']} targets gained at least one such transit, and {S['n_unique']} are left with a single possible
  period. As a test: where a period is already on record (a confirmed planet in NASA's archive, or one posted to
  ExoFOP), it is among the hunt's allowed periods in {S['n_agree']} of {S['n_onrec']} cases. The table lists the
  periods that are on no list, narrowed to three or fewer, with the pixel test of their transits. TOI-2493.01,
  above, is the strongest; the others are automated results that deserve the same one-by-one check before anyone
  relies on them.</p>
  </div>
  {table(df, cols) if len(df) else ""}
  <p class="muted" style="font-size:var(--step--1)">Shortest circular: the shortest period at which a central
  transit across the catalogued star would last as long as the one observed; periods far below it need an
  eccentric orbit or a different star.</p>
</section>"""


def section_trackb(tb, singles, lit=None):
    if tb is None or not len(tb):
        return ""
    litm = {int(x["tic"]): x["status"] for x in (lit or [])}
    show = trackb_shortlist(tb).sort_values("snr_max", ascending=False)
    n_eb = int((tb.verdict == "EB (deeper eclipses)").sum())
    n_nb = int((tb.verdict == "TOI on neighbour").sum())
    n_toi = int(tb.is_toi.sum()) if "is_toi" in tb else 0
    rest = tb[tb.verdict.isin(["single", "repeating"]) & ~tb.tic.isin(show.tic) & ~tb.is_toi.astype(bool)]
    kp_all = (show.known_planet.fillna("").astype(str) != "") if "known_planet" in show else pd.Series(False, index=show.index)
    lit_known = show.tic.map(lambda t: litm.get(int(t)) == "known in this role")
    conf_mask = kp_all | show.tic.isin(CONFIRMED_LIT)
    conf, cand = show[conf_mask], show[lit_known & ~conf_mask]

    def conf_name(r):
        if int(r.tic) in CONFIRMED_LIT:  # the planet the survey actually saw
            return CONFIRMED_LIT[int(r.tic)]
        kp = r.get("known_planet")
        return kp.split("; ")[0] if isinstance(kp, str) and kp else f"TIC {int(r.tic)}"
    ps_txt = (f", and all {len(show)} were matched by position against NASA's confirmed-planet table of 5 October 2026"
              if "known_planet" in show else "")
    n_on_record = int(sum(prior.prior(int(t))["kind"] is not None for t in show.tic))
    pub_txt = ""
    if len(conf):
        pub_txt += (f" {len(conf)} are confirmed planets or brown dwarfs that never became TOIs "
                    f"({esc(', '.join(sorted(conf_name(r) for _, r in conf.iterrows())))}), so the TOI filter could not "
                    "catch them; that the blind search found them is a useful check of its reach.")
    if len(cand):
        pub_txt += (f" {len(cand)} {'was' if len(cand) == 1 else 'were'} already published as "
                    f"{'a candidate' if len(cand) == 1 else 'candidates'} "
                    f"({esc(', '.join(f'TIC {int(t)}' for t in cand.tic))}).")
    if pub_txt:
        pub_txt += " All are marked \"published\" where they appear in the table."
    def published(r):
        kp = r.get("known_planet")
        return (kp if isinstance(kp, str) and kp else None), litm.get(int(r.tic))

    def tcell(r):
        pr = prior.prior(int(r.tic))
        label = {"blend": ("neighbour's planet", "eb"), "confirmed": ("published", "known"), "ctoi": ("CTOI", "known"),
                 "multiple": ("published", "known"), "binary": ("binary", "eb")}.get(pr["kind"])
        mark = " " + chip(*label) if label else ""
        if pr["notes"]:
            mark += f"<div class='muted' style='font-size:0.78rem;margin-top:2px'>{esc('; '.join(pr['notes']))}</div>"
        return f"<span class='m'>TIC {int(r.tic)}</span>{mark}"

    def p_cell(r):
        s = singles.get(str(int(r.tic))) if singles else None
        if s and s.get("p_min_allowed"):
            return f"≥ {s['p_min_allowed']:.0f}"
        p = r.get("period_guess")
        return fmt(p, 1) if p == p and r.verdict == "repeating" else "–"

    cols = [("Star", tcell, "l"), ("T mag", lambda r: fmt(r.Tmag, 1), "r"), ("T<sub>eff</sub>", lambda r: fmt(r.Teff, 0), "r"),
            (f"R★ ({RSUN})", lambda r: fmt(r.rad, 2), "r"), ("Sectors with dips", lambda r: esc(r.sectors), "l"),
            ("Depth (ppt)", lambda r: fmt(r.depth * 1e3, 1), "r"), ("T14 (h)", lambda r: fmt(r.t14_h, 1), "r"),
            (f"Rp ({RE})", lambda r: fmt(r.rp_rearth, 1), "r"), ("Period (d)", p_cell, "r"),
            ("SPOC TCEs", lambda r: fmt(int(r.n_spoc_tce)) if r.n_spoc_tce == r.n_spoc_tce and r.n_spoc_tce >= 0 else "–", "r"),
            ("Status", lambda r: verdict_chip(r.verdict), "l")]
    return f"""
<section id="trackb">
  <p class="kicker">Track B · single and duo transits</p>
  <h2>Other transit candidates</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>Every planet-sized dip on a dwarf star that passed the shape tests and is not a known TOI was followed up
  star by star: the dip finder was run on every TESS sector of the star, SPOC's own detections were read from
  the data-validation reports, and TOIs within 150″ were checked for a matching period. Of {len(tb)} such
  stars, {n_eb} turned out to be eclipsing binaries (deeper eclipses elsewhere in their data), {n_nb} carry
  a TOI on a neighbouring star with the same period, and {n_toi} became TOIs during 2026. Of the rest, {len(rest)}
  are left out here because they are larger than 16 Earth radii, repeat on periods shorter than 15 days, or were
  flagged in five or more SPOC searches without becoming TOIs (most likely rejected as binaries).
  {len(show)} remain; the 60 with the strongest dips are listed below. Most are giant-sized and could be brown dwarfs or small stars; all need the pixel and period work
  done for the two headline candidates. The 60 listed were also searched on the web by TIC and Gaia number{ps_txt};
  each is labelled with what ExoFOP, SIMBAD, VizieR and NASA's archive already hold. {n_on_record} of the
  {len(show)} are on record, most as single-transit community TOIs posted by the TESS Single Transit Planet
  Candidate working group or the Salinas et al. (2025) transformer search.{pub_txt}</p>
  </div>
  {table(show.head(60), cols)}
  <p class="muted" style="font-size:var(--step--1)">Rp assumes the catalogue stellar radius and corrects for
  contaminating flux. "Period" is the shortest period the other TESS data allow for a single transit, or the
  ephemeris through repeated events.</p>
</section>"""


def section_plates(items):
    if not items:
        return ""
    cards = []
    for path, title, text in items:
        cards.append(f'<figure class="plate" style="margin:0"><a href="{esc(path)}"><img src="{esc(path)}" '
                     f'alt="{esc(title)} light-curve dossier" decoding="async"></a><figcaption><strong>{esc(title)}</strong> '
                     f'{text}</figcaption></figure>')
    return f'<div class="plates">{"".join(cards)}</div>'


TRACKC_NOTES = {
    168789840: "Known sextuple (Powell et al. 2021); recovered blind.",
    309658221: "Published pulsating eclipsing binary (Lee et al. 2019); our second period is twice its δ Scuti pulsation period, so it is not an eclipse.",
    321471064: "Listed as a quadruple in the Multiple Star Catalog (WDS 11421-6216).",
    293954660: "Published as a doubly eclipsing pair (2.814 and 4.904 d; our second period is twice the latter) and found there to be a blend of two separate binaries.",
    389836747: "Published as a doubly eclipsing candidate (EA+EA, 2.567 and 2.730 d) in an A&A table of TESS multiply eclipsing systems (2022).",
}


def _json(name):
    p = CAND / name
    return json.load(open(p)) if p.exists() else {}


def eclipse_source(pix, ztf, j):
    """Where eclipse set j comes from: the centroid verdict, with ZTF's identification when it has one."""
    v = ((pix.get("centroid") or {}).get(str(j)) or {}).get("verdict")
    if not v:
        return "untested", None
    z = (ztf or {}).get(f"set{j}_eclipser") or []
    if v == "on target":
        txt, who = "target", "target"
    elif v.startswith("on G"):
        txt, who = v[3:].replace('"', "″"), v[3:]
    elif v.startswith("target or"):
        txt, who = v.replace('"', "″"), None
    elif v == "no usable centroids":
        return "untested", None
    else:
        txt, who = "offset, no catalogued star", None
    if z:
        zz = z[0].replace('"', "″")
        txt += f" · ZTF sees it on {zz}" if who == "target" else " · ZTF ✓" if who and z[0] == who else f" · ZTF: {zz}"
    return txt, who


def trackc_calibration():
    """How the centroid test does on Kostov et al. (2022) quadruples, whose authors traced both eclipse sets to one star."""
    k = _json("eb_pixels_kostov.json")
    n = {"both": 0, "neighbour": 0, "other": 0}
    for o in k.values():
        v = o.get("centroid_verdict")
        if v == "both sets on target":
            n["both"] += 1
        elif v == "a set on a neighbour":
            n["neighbour"] += 1
        elif o.get("centroid"):
            n["other"] += 1
    return n


def section_trackc(tc, lit=None):
    if tc is None or not len(tc):
        return ""
    pix = _json("eb_pixels.json")
    ztf = _json("ztf_check.json")
    como = pd.read_csv(CAND / "trackc_neighbours_comoving.csv") if (CAND / "trackc_neighbours_comoving.csv").exists() else None
    x = tc[(tc.P2_recurs >= 2)].copy()
    x = x.sort_values(["P2_recurs", "snr2_max"], ascending=False).head(30)
    rows = []
    for _, r in x.iterrows():
        t = int(r.tic)
        pr = prior.prior(t)
        p_ = pix.get(str(t), {})
        z = ztf.get(str(t))
        s1, w1 = eclipse_source(p_, z, 1)
        s2, w2 = eclipse_source(p_, z, 2)
        published = pr["kind"] == "multiple" or any(str(l.get("tic")) == str(t) and l.get("status") == "known in this role"
                                                     for l in (lit or []))
        if w1 == "target" and w2 == "target":
            verdict = "both on target" + (" (published)" if published else "")
        elif (w1 and w1 != "target") or (w2 and w2 != "target"):
            verdict = "blend"
        else:
            verdict = "unresolved"
        cm = ""
        if como is not None:
            c = como[(como.tic == t) & como.comoving]
            if len(c):
                cm = (f"The star the other set comes from moves with this one (proper motions within "
                      f"{c.dpm_masyr.iloc[0]:.1f} mas/yr; parallaxes differ by {abs(c.plx_diff_sigma.iloc[0]):.0f}σ): "
                      f"possibly a wide pair of binaries, {c.sep_au.iloc[0]:,.0f} AU apart.")
        rows.append(dict(tic=t, Tmag=r.Tmag, P1=r.P1, d1=r.depth1, s1=s1, P2=r.P2, d2=r.depth2, s2=s2,
                         rec=int(r.P2_recurs), verdict=verdict, kind=pr["kind"], note=prior.summary(t),
                         extra=" ".join(v for v in (TRACKC_NOTES.get(t, ""), cm) if v)))
    df = pd.DataFrame(rows)
    order = {"both on target (published)": 0, "both on target": 1, "unresolved": 2, "blend": 3}
    df = df.sort_values(["verdict", "Tmag"], key=lambda q: q.map(order) if q.name == "verdict" else q)
    kinds = {"both on target (published)": "known", "both on target": "new", "unresolved": "known", "blend": "eb"}
    cols = [("Star", lambda r: f"<span class='m'>TIC {r.tic}</span>", "l"), ("T", lambda r: fmt(r.Tmag, 1), "r"),
            ("P<sub>1</sub> (d)", lambda r: fmt(r.P1, 4), "r"), ("Depth", lambda r: f"{r.d1 * 100:.1f}%", "r"),
            ("Set 1 is on", lambda r: esc(r.s1), "l"),
            ("P<sub>2</sub> (d)", lambda r: fmt(r.P2, 4), "r"), ("Depth", lambda r: f"{r.d2 * 100:.2f}%", "r"),
            ("Set 2 is on", lambda r: esc(r.s2), "l"),
            ("Verdict", lambda r: chip(r.verdict, kinds.get(r.verdict, "known")) +
             (f"<div class='muted' style='font-size:0.78rem;margin-top:4px;max-width:22rem'>{esc(r.extra)}</div>" if r.extra else ""), "l"),
            ("On record", lambda r: f"<span class='muted' style='font-size:0.8rem'>{esc(r.note)}</span>", "l")]
    vc = df.verdict.value_counts().to_dict()
    n_pub_on = vc.get("both on target (published)", 0)
    n_new_on = vc.get("both on target", 0)
    n_blend = vc.get("blend", 0)
    n_unres = vc.get("unresolved", 0)
    n_ztf = sum(1 for _, r in df.iterrows() for sv in (r.s1, r.s2) if "ZTF ✓" in sv)
    cal = trackc_calibration()
    n_cal = sum(cal.values())
    cal_txt = (f"On {n_cal} of Kostov et al.'s (2022) quadruple candidates, whose authors traced both eclipse "
               f"sets to one star with their own pixel analysis, the test puts both sets on the target for "
               f"{cal['both']}, a set on a neighbour for {cal['neighbour']}, and is unresolved for {cal['other']}. "
               if n_cal else "")
    n_all = len(tc)
    return f"""
<section id="trackc">
  <p class="kicker">Track C · eclipsing triples and quadruples</p>
  <h2>Stars with two independent sets of eclipses</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>For every light curve with several deep dips, the main eclipsing period was found, a template of those
  eclipses was subtracted, and the residuals were searched for a second eclipse train. {n_all} stars showed a
  strong second signal that is not a harmonic of the first and is eclipse-shaped (box depth at least three times
  the swing of a sine wave at that period, which rejects pulsations and spots). The table lists those whose second
  period recurs in two or more sectors.</p>
  <p>A TESS pixel is 21″ wide, so two unrelated eclipsing binaries an arcminute apart blend into one light curve.
  Each eclipse set was therefore traced to its star: for every eclipse, the difference image (just before and
  after, minus during) shows only the light that disappeared, and a fitted point-spread function says where it
  sits, measured from the target's own image on that frame (Gaia positions moved to the date). Up to three sectors
  per star, every eclipse, combined with a 2″ systematic floor per sector. {cal_txt}Where ZTF's 1″ pixels reach
  (north of −30°), its light curves of the target and of every neighbour able to make the eclipses settle it
  directly: in {n_ztf} cases ZTF sees the eclipses on exactly the neighbour the centroids point to.</p>
  <p><strong>The outcome is clean and negative.</strong> All {n_pub_on} systems that come out with both eclipse
  sets on the target are already published as quadruples, multiply eclipsing systems or a sextuple; every
  unpublished one ({n_blend}) is a blend of two binaries on different stars, with the second star named in the
  table{f", and {n_unres} could not be traced" if n_unres else ""}. {f"{n_new_on} unpublished system{'s' if n_new_on != 1 else ''} keep both sets on the target. " if n_new_on else ""}Two of the blends are pairs of stars that move together
  on the sky, possibly wide binaries of binaries; their parallaxes disagree at 4–5σ, so this is not settled.</p>
  </div>
  {table(df, cols)}
</section>"""


def section_sources():
    """What the data beyond TESS's public light curves contributed (12 October 2026 extension)."""
    import glob as _g
    n_ffi = len(_g.glob(str(CATALOGS.parent / "cache" / "ffi_lc" / "*.npz")))
    kelt = _json("kelt_check.json")
    n_kelt = sum(1 for v in kelt.values() if v.get("status") == "searched")
    ztf = _json("ztf_check.json")
    n_ztf_sys = sum(1 for v in ztf.values() if v.get("status") == "checked")
    n_ztf_id = sum(len(v.get(f"set{j}_eclipser") or []) for v in ztf.values() for j in (1, 2))
    n_ztf_out = sum(len(v.get(f"set{j}_ruled_out") or []) for v in ztf.values() for j in (1, 2))
    g = pd.read_csv(CAND / "gaia_checks.csv") if (CAND / "gaia_checks.csv").exists() else pd.DataFrame()
    n_g = int((g.get("gaia_match") == "ok").sum()) if len(g) else 0
    n_orb = int(g.nss_orbits.notna().sum()) if len(g) else 0
    rows = [
        ("TESS full-frame images, sectors 99–108", "NASA · our own photometry",
         f"{n_ffi} light curves of {n_ffi} star-sectors that have no public full-frame light curve yet (QLP's end at "
         "sector 98, TESS-SPOC's near 80): two of TIC 402710280's five transits, and the extra sectors behind the "
         "period hunt's second pass."),
        ("Gaia DR3 (sources, orbits, eclipsing binaries)", "ESA · Gaia archive",
         f"{n_g} named stars matched: RUWE, resolved close companions, radial-velocity scatter and {n_orb} orbital "
         "solutions, quoted under “On record” throughout. In Track C, four spectroscopic orbits have the period of "
         "the eclipse set the pixels put on the target; for the blends, Gaia shows which neighbours move with the "
         "target. Positions with proper motions anchor every pixel test."),
        ("ZTF DR23 light curves", "NSF / Caltech · IRSA's copy on AWS",
         f"{n_ztf_sys} Track C systems north of −30°: the eclipses seen directly on {n_ztf_id} neighbours, and "
         f"{n_ztf_out} star–eclipse-set pairings ruled out (flat where a neighbour would need a deep eclipse). Every "
         "identification agrees with the TESS centroids."),
        ("KELT survey (2006–2018)", "NASA Exoplanet Archive",
         f"An anchored period search for the 28 bright single transits: KELT covers {n_kelt} of them, and injected "
         "copies of their TESS transits are not recovered, so its silence constrains nothing."),
        ("SuperWASP DR1", "NASA Exoplanet Archive", "Checked for the headline stars; too sparse to help."),
        ("Hipparcos–Gaia accelerations (Brandt 2021)", "ESA missions, via VizieR",
         "No unexplained acceleration for any named Hipparcos star (HIP 56966: χ² = 1.1)."),
        ("ExoFOP, SIMBAD, VizieR, NASA archive, ADS", "", "The novelty checks; ExoFOP re-read live for the new candidates."),
    ]
    tr = "".join(f"<tr><td><strong>{esc(a)}</strong><div class='muted' style='font-size:0.8rem'>{esc(b)}</div></td>"
                 f"<td>{esc(c)}</td></tr>" for a, b, c in rows)
    return f"""
<section id="sources">
  <p class="kicker">Beyond the public TESS light curves</p>
  <h2>Other telescopes, and the newest TESS images</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>Two things limit a search of public TESS light curves: the newest full-frame images have no light curves
  yet, and TESS's 21″ pixels cannot say which star in a crowded field dims. Both were attacked with other data,
  all public. ISRO's AstroSat, the one Indian space observatory with public data, observes in the ultraviolet and
  X-rays and has nothing to add on optical transits; ESA's CHEOPS and PLATO archives hold nothing on these stars
  yet.</p>
  </div>
  <div class="tablewrap"><table class="data"><thead><tr><th>Data</th><th>What it gave</th></tr></thead><tbody>{tr}</tbody></table></div>
</section>"""


def section_calibration(sv, s2, comp):
    q = sv[sv["product"] == "qlp"]
    sp = sv[sv["product"] == "spoc"]
    funnel = [
        ("Light curves searched", f"{int(sv.n_real.sum()):,}", f"{len(q)} full-frame (QLP) sectors and {len(sp)} two-minute sectors"),
        ("Dips at 6σ or more", f"{int(sv.n_events.sum()):,}", "box filter, 10 durations from 0.75 to 18 h, red-noise-aware threshold"),
        ("Clean dips", f"{int(sv.n_clean.sum()):,}", "after removing spacecraft and CCD systematics, blends, edges, gaps and background artefacts"),
        ("Single-transit candidates (Track B)", f"{int(sv.n_trackb.sum()):,}", "≤ 2 clean dips per sector, dwarf host, planet-sized"),
    ]
    if s2:
        funnel.append(("Planet-shaped after re-fitting", f"{s2['n_pl_real']:,}",
                       f"of {s2['n_vetted_real']:,} re-fitted: flat-bottomed, symmetric, unique, in raw flux"))
    frows = "".join(f"<tr><td>{a}</td><td class='r m'>{b}</td><td class='muted'>{c}</td></tr>" for a, b, c in funnel)
    nrows = ""
    if s2:
        for prod, name in (("spoc", "Two-minute (SPOC)"), ("qlp", "Full-frame (QLP)")):
            b = s2["by_product"][prod]
            rr = 1e3 * b["real"] / max(b["n_real_lc"], 1)
            fr = 1e3 * b["flip"] / max(b["n_flip_lc"], 1)
            nrows += (f"<tr><td>{name}</td><td class='r m'>{b['n_real_lc']:,}</td><td class='r m'>{b['real']}</td>"
                      f"<td class='r m'>{rr:.2f}</td><td class='r m'>{b['n_flip_lc']:,}</td><td class='r m'>{b['flip']}</td>"
                      f"<td class='r m'>{fr:.2f}</td></tr>")
        rates = [(1e3 * x["flip"] / x["n_flip_lc"], x["sector"]) for x in s2["by_sector"] if x["n_flip_lc"] > 0]
        med_rate = float(np.median([r for r, _ in rates])) if rates else float("nan")
        max_rate, max_sec = max(rates) if rates else (float("nan"), "–")
    crow = ""
    if comp:
        crow = "".join(f"<tr><td>{r['bin']}</td><td class='r m'>{r['n']}</td><td class='r m'>{r['found']:.0%}</td></tr>"
                       for r in comp["rows"])
    return f"""
<section id="calibration">
  <p class="kicker">Calibration</p>
  <h2>How much to trust a dip</h2>
  <div class="two">
    <div>
      <h3>From {sv.n_real.sum() / 1e6:.1f} million light curves to a handful</h3>
      <div class="tablewrap"><table><tbody>{frows}</tbody></table></div>
    </div>
    <div>
      <h3>False alarms, measured</h3>
      <p>One light curve in ten was also searched upside down, through the identical pipeline. Anything found
      there is a false alarm by construction. Counted here: planet-shaped dips longer than 2 hours on stars hotter
      than 4,000 K, per thousand light curves.</p>
      <div class="tablewrap"><table><thead><tr><th>Data</th><th class="r">LCs</th><th class="r">Found</th>
      <th class="r">/1000</th><th class="r">Flipped LCs</th><th class="r">Found</th><th class="r">/1000</th></tr></thead>
      <tbody>{nrows}</tbody></table></div>
      <p>In two-minute data the flipped light curves produced nothing, so planet-shaped dips there are essentially
      all astrophysical (transits or eclipses). In full-frame data the flipped rate varies by sector (median
      {med_rate:.1f}, highest {max_rate:.1f} per 1000 in sector {max_sec}), driven by scattered light, which adds
      upward bumps that turn into dips when flipped. It is an upper bound: systematics are vetoed when they hit many
      stars at once, and flipping only 10% of light curves makes that veto ten times weaker for the flipped set. A
      single full-frame dip therefore needs the pixel and period checks before it counts; all three headline
      candidates are in two-minute data.</p>
    </div>
  </div>
  <div class="two">
    <div>
      <h3>Recovering known planets blind</h3>
      <p>For TOIs with periods over 15 days and depths over 1.5 ppt, we predicted every transit that fell in a
      surveyed full-frame sector ({comp['n'] if comp else 0} transits) and asked whether the survey listed it as a
      single-transit candidate. Shallow ones are missed more often, as expected; the rest are lost to data gaps,
      ephemeris drift since the TOI was found, and the deliberately strict cuts.</p>
      <div class="tablewrap"><table><thead><tr><th>Depth</th><th class="r">Transits</th><th class="r">Recovered</th></tr></thead>
      <tbody>{crow}</tbody></table></div>
    </div>
    <div>
      <h3>Known objects found again</h3>
      <ul class="checks">
        <li class="check check-pass"><span class="ci" aria-hidden="true">✓</span><div><strong>TIC 168789840.</strong>
        The sextuply eclipsing sextuple star (Powell et al. 2021): all three eclipse periods, 1.570, 8.22 and 1.306 d, recovered in every sector tested.</div></li>
        <li class="check check-pass"><span class="ci" aria-hidden="true">✓</span><div><strong>TOI-2583.01.</strong>
        A confirmed planet whose two transits the survey found and passed as planet-shaped, blind.</div></li>
        <li class="check check-pass"><span class="ci" aria-hidden="true">✓</span><div><strong>TOI-7714 (TIC 423785115).</strong>
        Found as a new candidate from sector 107, then matched to a TOI created in April 2026, after the older TOI list.</div></li>
        <li class="check check-pass"><span class="ci" aria-hidden="true">✓</span><div><strong>HD 100546.</strong>
        The Herbig Be star's irregular dimming was flagged by the anomaly track.</div></li>
        <li class="check check-pass"><span class="ci" aria-hidden="true">✓</span><div><strong>WASP-18, WASP-126.</strong>
        Recovered by the periodic search, with centroid offsets of 0.04 pixel in the difference images.</div></li>
        <li class="check check-pass"><span class="ci" aria-hidden="true">✓</span><div><strong>Confirmed planets that never became TOIs.</strong>
        Found blind by the single-transit track and identified only afterwards, by position, in NASA's confirmed-planet table:
        {"; ".join(f"{n} ({d})" for n, d in RECOVERED)}.</div></li>
        <li class="check check-warn"><span class="ci" aria-hidden="true">!</span><div><strong>One lost.</strong>
        {"; ".join(f"{n}, {d}" for n, d in MISSED)}. The deeper dips that triggered the veto elsewhere in its data
        were not checked by eye; the veto may be throwing out other planets too.</div></li>
      </ul>
    </div>
  </div>
</section>"""


def section_anomalies(td):
    rows = ""
    if td is not None and len(td):
        x = td[td.verdict == "anomaly candidate"]
        rows = f"{len(x)} dips passed the anomaly cuts and the all-sector check. "
    return f"""
<section id="trackd">
  <p class="kicker">Track D · anomalies</p>
  <h2>Dippers and other odd dips</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>This track keeps dips that no planet makes: strongly asymmetric ones, and ones that last longer than
  15 hours. {rows}Looked at by eye, they fall into familiar groups: the known dipper HD 100546 (recovered), eclipses
  of a nearby binary leaking into a full-frame aperture, fast-rotating M dwarfs whose spots the detrending
  could not follow, very short-period binaries sampled into stroboscopic patterns, and long-period eclipsing
  binaries caught in one or two deep eclipses (TIC 304213412: 18-hour, 17%-deep eclipses repeating over six years;
  TIC 320329350: a 32%-deep, 17-hour eclipse). None is yet a convincing new kind of object, so none is claimed
  here.</p>
  </div>
</section>"""


def section_tracka(ta_n, ta_total, triage):
    sig = ""
    if triage is not None and len(triage):
        clean = triage[(triage.kind == "signals") & (triage.hard_flags.fillna("") == "")]
        flip = triage[(triage.kind == "signals_inverted") & (triage.hard_flags.fillna("") == "")]
        if len(clean):
            sig = (f" Of {triage.tic.nunique():,} stars with a signal passing the search thresholds, {clean.tic.nunique()} have one "
                   f"with no red flag; the strongest has SDE {clean.sde.max():.1f}, against {flip.sde.max():.1f} for the strongest "
                   f"flag-free signal in the same stars' flipped light curves." if len(flip) else "")
    return f"""
<section id="tracka">
  <p class="kicker">Track A · small planets</p>
  <h2>A deeper search of nearby red dwarfs</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>{ta_total:,} K and M dwarfs within 150 parsecs, with years of TESS coverage but few two-minute sectors, are
  being searched with every sector stitched together, for planets down to about 1.5 Earth radii. {ta_n:,} have been
  searched so far.{sig} No candidate has cleared that bar yet. The search continues and will be reported in the
  same way, against its own flipped-data baseline.</p>
  </div>
</section>"""


def section_methods():
    return f"""
<section id="methods">
  <p class="kicker">Method</p>
  <h2>How the search works</h2>
  <div class="two">
    <div><h3>Data</h3><p>Everything comes from the public MAST archive on Amazon S3: QLP full-frame light curves
    (sectors 1–98), SPOC two-minute light curves (to sector 107) with their data-validation reports, full-frame image
    cubes for pixel work, and the TESS Input Catalog and Gaia DR3. Light curves are read straight from the bucket;
    nothing is pre-selected by other pipelines.</p></div>
    <div><h3>Detrending and dips</h3><p>Each light curve is binned to 10 minutes and flattened with a robust
    sliding biweight. A box filter at ten durations finds dips; each threshold uses the measured scatter of the
    filtered series itself, so slow correlated noise is counted. Upward outliers (flares) are clipped; downward
    points never are.</p></div>
    <div><h3>Systematics</h3><p>Dips that strike many stars on the same CCD, or across the spacecraft, at the same
    time are vetoed. So are dips at data edges, in gaps, during background spikes, or absent from the raw
    aperture flux, and dips shared by nearby stars (only the deepest is kept).</p></div>
    <div><h3>Vetting</h3><p>Survivors are re-fitted with trapezoid, asymmetric and Gaussian models; checked against
    every other sector of the star; matched against TOIs, SPOC detections and the literature; and examined in
    the pixels with difference images and Gaia neighbours.</p></div>
  </div>
  <p class="muted col">Full details, thresholds and their sources: <a href="{REPO}/docs/METHODS.md">docs/METHODS.md</a>.
  Code: <a href="{REPO}/pipeline">pipeline/</a>.</p>
</section>"""


def section_caveats(lit=None):
    import exofop_check
    named = exofop_check.named_stars()
    n_lit = len(named)
    n_known = sum(1 for t in named if prior.prior(t)["kind"] is not None)
    return f"""
<section id="next">
  <p class="kicker">Limits and next steps</p>
  <h2>What would make these discoveries</h2>
  <div class="two">
    <div>
      <h3>What this report cannot rule out</h3>
      <ul class="checks">
        <li class="check check-info"><span class="ci" aria-hidden="true">i</span><div><strong>Unresolved companions.</strong>
        A faint eclipsing binary within an arcsecond or two of the target mimics a planet. Gaia excludes wider ones; high-resolution imaging and radial velocities exclude the rest.</div></li>
        <li class="check check-info"><span class="ci" aria-hidden="true">i</span><div><strong>What was checked, and what was not.</strong>
        On 11–12 October 2026 every star this report names was looked up on ExoFOP (TOIs, community TOIs, the periods
        and notes people have posted, follow-up), in SIMBAD, in VizieR (TESS and Gaia DR3 eclipsing binaries, VSX,
        published quadruple and multiply eclipsing lists), in NASA's confirmed-planet table, including the transit
        times of planets on stars up to 2′ away, and in the literature (web search, and ADS full text by every
        identifier). Of the {n_lit} stars, {n_known} turned out to be on record in some form, most of them as
        community candidates on ExoFOP; each table says what was found. Not matched row by row: candidate tables
        published only inside papers, such as the T16 hunt's 10,000 candidates from TESS's first year.</div></li>
        <li class="check check-info"><span class="ci" aria-hidden="true">i</span><div><strong>Candidates, not discoveries.</strong>
        Every object here is a candidate until follow-up observations rule out the impostors above; the search is complete
        for TESS sectors 1–107, and the vetting is deepest for stars brighter than T = 10.5.</div></li>
      </ul>
    </div>
    <div>
      <h3>Turning a candidate into a planet</h3>
      <p><strong>Report it.</strong> HIP 56966 and HD 11506 already have community TOIs on ExoFOP, so the new
      period and the second transit belong there, as an update to the existing entries, the way another user
      posted HIP 71639's period in May. The bright single transits with nothing on record would be new community
      TOIs. A short Research Note of the AAS is the usual citable record for an amateur find.</p>
      <p><strong>Rule out impostors.</strong> HIP 71639 and HIP 56966 already have Gemini speckle imaging on ExoFOP
      showing no close companion. HD 11506 has none yet; it and the new single transits need speckle or adaptive-optics
      imaging and a few high-precision radial velocities.</p>
      <p><strong>Weigh it.</strong> HIP 71639's expected signal is about 3 m/s, well within reach of HARPS or
      ESPRESSO for a V = 6 star. A transit from CHEOPS, or the full 7-hour transit on 20 April or 14 May 2027 from
      Chile, would confirm the period.</p>
    </div>
  </div>
</section>"""


# ------------------------------------------------------------------------------------------- build

PLATES = [
    ("TIC60479865_duo.png", "TIC 60479865",
     "K dwarf, T = 11.2. Two 1.1%-deep, 6-hour transits, in sector 33 (2021) and sector 107 (2026): a Saturn-sized "
     "duo with 34 allowed periods from 32 days up."),
    ("TIC206385646_single.png", "TIC 206385646",
     "F star, T = 10.8. One 0.9%-deep, 13-hour transit in sector 97. Jupiter-sized; the long duration points to "
     "a period of a year or more."),
    ("TIC358050122_single.png", "TIC 358050122",
     "F star, T = 11.3. One 1.0%-deep, 12-hour transit in sector 97, flat-bottomed, in nine sectors of otherwise "
     "flat data."),
    ("TIC54002406_single.png", "TIC 54002406",
     "F star, T = 11.5. One 0.6%-deep, 13-hour transit in sector 97. Jupiter-sized."),
    ("TIC341179881_B.png", "TIC 341179881",
     "F star, T = 11.9. 1.7%-deep, 8.7-hour transits that repeat; a large planet or a small star."),
    ("TIC265527315_B.png", "TIC 265527315",
     "K dwarf, T = 11.4. 3.6%-deep, 7.3-hour transits that repeat; 1.6 Jupiter radii, so more likely a brown "
     "dwarf or M dwarf. A secondary-eclipse search is still to do."),
]


def build():
    sv = survey_stats()
    s2 = stage2_stats(sv)
    comp = completeness()
    fit = json.load(open(OUT / "transit_fits.json"))
    ex = json.load(open(OUT / "extras.json"))
    tb = pd.read_csv(CAND / "trackb_final.csv") if (CAND / "trackb_final.csv").exists() else None
    tc = pd.read_csv(CAND / "trackc_candidates.csv") if (CAND / "trackc_candidates.csv").exists() else None
    td = pd.read_csv(CAND / "trackd_final.csv") if (CAND / "trackd_final.csv").exists() else None
    singles = json.load(open(CAND / "trackb_singles.json")) if (CAND / "trackb_singles.json").exists() else {}
    dv = json.load(open(CAND / "deep_vet.json")) if (CAND / "deep_vet.json").exists() else {}
    litp = RESULTS / "candidates" / "literature_checks.json"
    lit = json.load(open(litp)) if litp.exists() else None
    ta_dir = RESULTS / "track_a"
    ta_n = len(glob.glob(str(ta_dir / "stars" / "*.json")))
    ta_total = len(pd.read_csv(ta_dir / "targets_search.csv"))
    triage = pd.read_csv(ta_dir / "triage.csv") if (ta_dir / "triage.csv").exists() else None
    n_q = int((sv["product"] == "qlp").sum())
    n_s = int((sv["product"] == "spoc").sum())
    today = datetime.now(timezone.utc).strftime("%-d %B %Y")
    f71 = fit["HIP71639"]
    rp71 = f71["rp_rs"][0] * 1.596 * R_EARTH_PER_RSUN
    rp115 = fit["HD11506"]["rp_rs"][0] * 1.35 * R_EARTH_PER_RSUN
    f56 = ""
    if "HIP56966" in fit:
        rp56 = fit["HIP56966"]["rp_rs"][0] * 1.069 * R_EARTH_PER_RSUN
        f56 = (f'<li>{chip("New period", "new")}<div><strong>HIP 56966 (HD 101544): a {rp56:.1f} {RE} planet candidate '
               f'every {fit["HIP56966"]["period_d"][0]:.3f} days around a V = 8.8 Sun-like star.</strong> Its 2019 transit '
               f'was posted to ExoFOP in 2025 as a single-transit community candidate (Salinas et al.), with no period. '
               f'This survey found three more transits across seven years and the period they share.</div></li>')
    fnew = ""
    if "TIC402710280" in fit and "TOI2493" in fit:
        fa, fb = fit["TIC402710280"], fit["TOI2493"]
        exa, exb = ex.get("TIC402710280", {}), ex.get("TOI2493", {})
        rpa = fa["rp_rs"][0] * 0.817 * R_EARTH_PER_RSUN
        rpb = fb["rp_rs"][0] * 0.652 * R_EARTH_PER_RSUN
        fnew = (f'<li>{chip("New candidate", "new")}<div><strong>TIC 402710280: a {rpa:.1f} {RE} planet candidate every '
                f'{fa["period_d"][0]:.3f} days around a K dwarf 67 parsecs away (T = 9.2).</strong> Not a TOI, not a '
                f'community TOI, nothing on record anywhere. Five transits from 2021 to 2026, two of them only in our own '
                f'light curves of the newest TESS full-frame images; the dip is on the star and has a planet\'s shape.</div></li>'
                f'<li>{chip("Period new", "new")}<div><strong>TOI-2493.01: a {rpb:.1f} {RE} candidate every '
                f'{fb["period_d"][0]:.2f} days, receiving {exb.get("s_earth", float("nan")):.1f}× Earth\'s insolation.</strong> '
                f'A single-transit TOI since 2021 with no period; the period hunt found two more transits in 2025 and 2026. '
                f'TESS should see the next one in mid-November 2026.</div></li>')
    fph = section_periodhunt(summary_only=True)
    passed = [int(k) for k, r in dv.items() if "error" not in r and deep_verdict(r)[0] == "candidate"]
    n_dv = sum(1 for t in passed if prior.prior(t)["kind"] is None)
    n_dv_known = len(passed) - n_dv
    n_dv_all = sum(1 for r in dv.values() if "error" not in r)
    tbf = trackb_shortlist(tb) if tb is not None else None
    n_tb = len(tbf) if tbf is not None else 0
    n_bright = int((tbf.Tmag < 10).sum()) if tbf is not None else 0
    n_tc = int((tc.P2_recurs >= 2).sum()) if tc is not None else 0
    def plate_note(t):
        note = prior.summary(int(t.split()[-1]))
        return " Nothing on record." if note == "nothing on record" else f" On record: {note}."
    plates = section_plates([(f"img/{p}", t, esc(x + plate_note(t))) for p, t, x in PLATES])
    body = f"""<title>Ad Astra TESS Survey</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wdth,wght@100..125,400..800&family=IBM+Plex+Mono:wght@400;500;600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400..700;1,8..60,400&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header class="mast">
  <p class="kicker">Amateur search of NASA TESS data · survey report · {today}</p>
  <h1>Ad Astra</h1>
  <p class="dek">A blind search of {int(sv.n_real.sum()) / 1e6:.1f} million TESS light curves for transits and
  eclipses that no catalogue lists yet, extended with our own light curves of the newest TESS images and with
  Gaia, ZTF and KELT data. Checked against ExoFOP, SIMBAD, VizieR, NASA's archive and the literature, most of what
  it found is already known to someone. What is left: a sub-Neptune candidate every 36.3 days around the nearby K
  dwarf TIC 402710280, which is on no list; the 79-day period of TOI-2493.01, a temperate sub-Neptune candidate;
  the period of a sub-Neptune around HIP 56966; a second transit of a small planet around HD 11506; and {n_dv}
  single transits on bright stars that nobody has posted.</p>
  <div class="strip">
    <div><span class="label">Sectors searched</span><span class="v">{n_q + n_s} of 107</span></div>
    <div><span class="label">Light curves</span><span class="v">{int(sv.n_real.sum()):,}</span></div>
    <div><span class="label">Red-dwarf deep search</span><span class="v">{ta_n:,} of {ta_total:,}</span></div>
  </div>
</header>
<nav class="toc" aria-label="Sections"><a href="#findings">Findings</a><a href="#tic402710280">TIC 402710280</a>
<a href="#toi2493">TOI-2493</a><a href="#hip56966">HIP 56966</a><a href="#hd11506">HD 11506</a><a href="#hip71639">HIP 71639</a>
<a href="#periodhunt">Periods for TOIs</a><a href="#bright">Bright candidates</a><a href="#trackb">All transit candidates</a>
<a href="#trackc">Multiple stars</a><a href="#sources">Other data</a><a href="#trackd">Anomalies</a><a href="#calibration">Calibration</a>
<a href="#methods">Method</a><a href="#next">Next steps</a></nav>

<section id="findings">
  <p class="kicker">Summary</p>
  <h2>What the search found</h2>
  <ul class="findings">
    {fnew}
    {f56}
    <li>{chip("Second transit", "new")}<div><strong>HD 11506: a {rp115:.1f} {RE} transiting candidate in a known
    three-planet system.</strong> The 2020 transit was posted as a single-transit community candidate in 2024, but
    with a depth and duration ten times too large. This survey measures it properly and finds a second transit in
    August 2026, 2,149 days later, so the period is 25.9 days or longer. The star's three giant planets, found by
    radial velocities, do not fit.</div></li>
    <li>{chip("Confirmed independently", "known")}<div><strong>HIP 71639: a {rp71:.1f} {RE} planet candidate every
    {f71['period_d'][0]:.3f} days around a V = 6 F star.</strong> Planet Hunters TESS posted its first transit in
    2020, and another ExoFOP user posted this period from the same four transits in May 2026. The survey found it
    independently, five months later, and agrees within 0.1σ. It is still not a TOI.</div></li>
    <li>{chip(f"{n_dv} bright", "new")}<div><strong>{n_dv} transit candidates on bright stars that nobody has
    posted.</strong> Of {n_dv_all} shortlisted stars brighter than T = 10.5, {n_dv + n_dv_known} pass the pixel,
    period and binarity tests; {n_dv_known} of those were already on record, as community candidates, confirmed
    planets (HD 3167, HD 278555, K2-290), binaries, or once a transit of HD 80606 b on the neighbouring star. The
    {n_dv} left are mostly single transits of giant planets, among them the bright stars HD 41179 and HD 43475, plus
    a grazing companion every 24.05 days around CD−34 9313 and two smaller ones (4.2 and 5.8 {RE}).</div></li>
    <li>{chip("None new", "known")}<div><strong>Triple and quadruple stars: every unpublished one is a blend.</strong>
    Of {n_tc} stars with two independent eclipse sets, the six whose eclipses all come from the target are already
    published quadruples or a sextuple. In every other one, a pixel test traces one eclipse set to a neighbouring
    star, and ZTF's sharper images see those eclipses on the named neighbour wherever they reach.</div></li>
    {fph}
    <li>{chip("Measured", "ok")}<div><strong>The false-alarm floor, and what was found again.</strong> Flipped light
    curves produced no planet-shaped dips in two-minute data. Blind, the survey found nine confirmed planetary
    systems that never became TOIs, 83 community candidates, published quadruple candidates, a sextuple star and a dipper.</div></li>
  </ul>
</section>
{section_tic402710280(fit, ex)}
{section_toi2493(fit, ex)}
{section_hip56966(fit, ex)}
{section_hd11506(fit, ex)}
{section_hip71639(fit, ex)}
{section_periodhunt()}
{section_deepvet(dv, lit, tb)}
{section_trackb(tb, singles, lit)}
<section id="plates" style="padding-top:8px">
  <h3>Selected dossiers</h3>
  <p class="muted col" style="font-size:var(--step--1)">Automatic dossiers: the dip with its fit (top left), the
  difference image with Gaia neighbours (top right), and every TESS sector of the star. Select one to enlarge.</p>
  {plates}
</section>
{section_trackc(tc, lit)}
{section_sources()}
{section_anomalies(td)}
{section_tracka(ta_n, ta_total, triage)}
{section_calibration(sv, s2, comp)}
{section_methods()}
{section_caveats(lit)}
<footer>
  <p>Data: TESS mission (NASA), served by MAST from its public Amazon S3 bucket; light curves from the MIT Quick-Look
  Pipeline (QLP), TESS-SPOC and SPOC; TESS Input Catalog v8.2; Gaia DR3 (ESA). TOI tables: ExoFOP (February 2026)
  and the NASA Exoplanet Archive (5 October 2026). BTJD = BJD − 2457000.</p>
  <p>Generated {today} from <a href="{REPO}/results">the results files</a> by
  <a href="{REPO}/pipeline/report.py">pipeline/report.py</a>.</p>
</footer>
</div>
<script>{JS}</script>
"""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(body)
    print("wrote", OUT / "index.html", len(body) // 1024, "KB")



# ----------------------------------------------------------------------------------- deep-vetted list

NOTES = {
    70803960: "K2-290 c, a known 48.4-day warm Jupiter (the system also has a 9.2-day planet, K2-290 b).",
    162305348: "Short (1.9 h) for a 35-day orbit around an F star: a grazing crossing, so the companion may be larger than it looks.",
    457134361: "One 11.6-hour, 1%-deep transit; at least 17.6 days, and long for a Sun-like star, so likely hundreds of days.",
    4011632: "Gaia reports an astrometric binary (RUWE 5.7) with a variable radial velocity: the dip may belong to the companion.",
    221567884: "HD 148797, with two confirmed transiting planets (42.1 and 68.2 d): the dips of two depths are the two planets.",
    459801869: "Two transits 20.4 days apart in sector 23 and every shorter alias excluded: a warm giant on an F star, seen only in full-frame data.",
    137966491: "HD 278555 b, a confirmed 28.84-day planet that never became a TOI; the survey found its period blind.",
    441010857: "Three transits on a 24.054-day period; the 1.7-hour duration means a grazing crossing.",
    318707987: "HD 3167 (= HIP 2736), a known K2 system: the single transit is HD 3167 c (29.85 d), and the 0.96-day SPOC detection is HD 3167 b.",
    430746529: "One 16-hour transit of a bright A star (HIP 109818): a giant on an orbit of at least 65 days.",
    402710280: "Two 4.2 Earth-radius transits on a K dwarf; six periods of 1454/n days remain.",
    64038538: "Two 4.4 Earth-radius transits on a K dwarf; six periods of 1116/n days remain.",
}


_CV = None


def cv(tic):
    """Difference-image centroid result for a transit candidate (centroid_vet.py), or {}."""
    global _CV
    if _CV is None:
        p = CAND / "centroid_vet.json"
        _CV = json.load(open(p)) if p.exists() else {}
    return _CV.get(str(int(tic)), {})


def cv_text(r):
    """The centroid result in words, with its numbers."""
    v = r.get("verdict")
    c = r.get("centroid") or {}
    if not v or v == "error":
        return "not tested"
    if v == "no usable centroids":
        return "no usable centroid (dip too shallow in the pixels)"
    near = (r.get("closest_neighbours") or [None])[0]
    s = f"{r['offset_arcsec']:.1f}″ ± {c['err']:.1f}″ from the star"
    if near:
        s += f"; nearest star able to fake it ({near['name']}) {near['sigma']:.0f}σ away"
    return f"{v}: {s}"


def deep_verdict(r):
    g = r.get("gaia") or {}
    c = cv(r["tic"]) if "tic" in r else {}
    v = c.get("verdict", "")
    if v.startswith("on G"):
        return "off-target", v[3:]
    if g.get("ruwe", 0) > 1.4 or g.get("nss", 0) > 0:
        return "binary host", None
    p = r.get("period") or {}
    n_ev = len([e for e in r.get("events", []) if e.get("sector") is not None])
    if ("aliases" in p and not any(a["status"] == "allowed" for a in p["aliases"])) or (n_ev >= 2 and not p):
        return "no consistent period", None
    if v.startswith("target or") or v.startswith("offset"):
        return "ambiguous pixels", v
    return "candidate", None


def section_deepvet(dv, lit, tb):
    if not dv:
        return ""
    litm = {int(x["tic"]): x for x in (lit or [])}
    rows = []
    for k, r in dv.items():
        if "error" in r:
            continue
        tic = int(k)
        st = r["star"]
        b = tb[tb.tic == tic]
        rp = float(b.rp_rearth.iloc[0]) if len(b) else float("nan")
        ev = [e for e in r["events"] if e["sector"] is not None]
        p = r.get("period") or {}
        if "aliases" in p:
            ok = [a["P"] for a in p["aliases"] if a["status"] == "allowed"]
            per = (f"{ok[0]:.3f} d" if len(ok) == 1 else
                   f"{p['P_longest']:.0f} d / n, {len(ok)} of {len(p['aliases'])} allowed" if ok else "none fits")
        elif p.get("single"):
            per = f"≥ {p['p_min_allowed']:.0f} d" if p.get("p_min_allowed") else "–"
        else:
            per = "–"
        pix = cv_text(cv(tic))
        g = r.get("gaia") or {}
        verdict, _ = deep_verdict(r)
        pr = prior.prior(tic)
        lt = prior.summary(tic)
        if verdict == "candidate":
            verdict = {"blend": "neighbour's planet", "confirmed": "published", "ctoi": "published",
                       "multiple": "published", "binary": "binary host"}.get(pr["kind"], verdict)
        note = NOTES.get(tic, "")
        rows.append({"tic": tic, "hip": st.get("HIP"), "tmag": st.get("Tmag"), "rp": rp, "n_ev": len(ev),
                     "sectors": ", ".join(str(e["sector"]) for e in ev), "per": per, "pix": pix,
                     "ruwe": g.get("ruwe"), "nss": g.get("nss"), "tce": len((r.get("known") or {}).get("tces", [])),
                     "toi": len((r.get("known") or {}).get("toi", [])), "lit": lt, "verdict": verdict,
                     "note": note})
    df = pd.DataFrame(rows).sort_values(["verdict", "tmag"], key=lambda s: s.map({"candidate": 0}).fillna(1) if s.name == "verdict" else s)
    kinds = {"candidate": "new", "off-target": "eb", "binary host": "eb", "no consistent period": "known",
             "ambiguous pixels": "known", "published": "known", "neighbour's planet": "eb"}

    def star(r):
        h = f" · HIP {int(r.hip)}" if r.hip == r.hip and r.hip else ""
        return f"<span class='m'>TIC {r.tic}</span>{h}"
    cols = [("Star", star, "l"), ("T", lambda r: fmt(r.tmag, 1), "r"), (f"Rp ({RE})", lambda r: fmt(r.rp, 1), "r"),
            ("Transits", lambda r: f"{r.n_ev} <span class='muted'>(S{esc(r.sectors)})</span>", "l"),
            ("Period", lambda r: esc(r.per), "l"), ("Pixels", lambda r: esc(r.pix), "l"),
            ("RUWE", lambda r: fmt(r.ruwe, 2), "r"), ("TCE", lambda r: str(r.tce), "r"),
            ("On record", lambda r: f"<span class='muted' style='font-size:0.8rem'>{esc(r.lit)}</span>", "l"),
            ("Verdict", lambda r: chip(r.verdict, kinds.get(r.verdict, "known")) +
             (f"<div class='muted' style='font-size:0.78rem;margin-top:4px;max-width:22rem'>{esc(r.note)}</div>" if r.note else ""), "l")]
    n_c = int((df.verdict == "candidate").sum())
    counts = df.verdict.value_counts().to_dict()
    n_pass = n_c + counts.get("published", 0) + counts.get("neighbour's planet", 0)
    n_pub = counts.get("published", 0)
    rejected = ", ".join(f"{v} {k}" for k, v in counts.items() if k not in ("candidate", "published", "neighbour's planet"))
    df = df[df.verdict.isin(["candidate", "published", "neighbour's planet"])].sort_values(
        ["verdict", "rp"], key=lambda s: s.map({"candidate": 0, "published": 1}).fillna(2) if s.name == "verdict" else s)
    return f"""
<section id="bright">
  <p class="kicker">Track B · deep vetting</p>
  <h2>The brightest remaining candidates, checked one by one</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>Every shortlisted star brighter than T = 10.5 went through the same checks as the headline candidates:
  Gaia binarity, TOI and SPOC lookups, the depth of each dip in every sector, the period aliases, and a
  pixel-level test of where the dip comes from (difference-image centroids, compared with every Gaia star within
  150″ bright enough to fake the dip). {n_pass} of {sum(counts.values())} pass ({rejected} do not). Each was then looked up on ExoFOP
  (TOIs, community TOIs and the periods people have posted), in SIMBAD, in VizieR's eclipsing-binary and
  multiple-star catalogues, in NASA's confirmed-planet table (including the transit times of planets on
  neighbouring stars) and in the literature. {n_pub} were already on record, mostly as single-transit community
  candidates, and one dip is a transit of HD 80606 b on the neighbouring star. <strong>{n_c} have nothing on record
  anywhere</strong>; they head the table, smallest first. HIP 56966 above came out of this list.</p>
  </div>
  {table(df, cols)}
  <p class="muted" style="font-size:var(--step--1)">Pixels: where the light that disappears in transit comes
  from, by difference-image centroids (2-minute pixels where available, full-frame cutouts otherwise), measured
  from the star's own image on each frame and combined over sectors with a 2″ systematic floor per sector; the
  σ is how far the nearest Gaia star able to fake the dip lies from that position. RUWE above 1.4 hints at an
  unresolved binary.</p>
</section>"""


if __name__ == "__main__":
    build()
