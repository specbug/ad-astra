"""Build the survey report page (results/report/index.html) from the results files.

  python report.py

Numbers come from results/survey/*/summary.json, results/candidates/*, results/report/transit_fits.json and
results/report/extras.json; the narrative for the headline candidates lives in CURATED below. Figures are
drawn as inline SVG styled by the page's colour tokens, so they follow the viewer's light/dark theme.
"""
import glob
import html
import json
import math
from datetime import datetime, timezone

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import RESULTS

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
    <p class="kicker">Candidate 1 · Track B (long-period transits) · SPOC 2-minute data</p>
    <h2>HIP 71639: a warm Neptune every 23.9 days</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("Planet Hunters single-transit candidate", "known")}{chip("SPOC TCE with a 31× wrong period", "tce")}</div>
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
        nearby companion; the ground-based data were left for a later paper. With one transit there was no
        period, and web searches on 11 October 2026 found no paper giving one.</p>
        <p>NASA's SPOC pipeline also detected the signal. Its multi-sector searches of sectors 1–39 and then 1–65
        found the transits in sectors 11, 38 and 65 and reported a threshold-crossing event (multiple-event
        statistic 26.6). One transit fell in each sector, two years apart, so the pipeline took the gap between
        them as the period: <span class="num">740.42 d</span>. At that period the transit is far too short for the
        star (its fit implied a stellar density of 8.6 ± 6.5 {RHOSUN}), and the centroid test was flagged. In the
        TOI lists available to us it never became a TOI.</p>
        <p>What this survey adds is the period. The fourth transit, found in sector 102 (April 2026), arrives
        {abs(e['s102_offset_min']):.0f} minutes from where 106 cycles of a {f['period_d'][0]:.3f}-day period put
        it. 740.42 days is exactly 31 of those cycles. ExoFOP itself, where the candidate's follow-up notes live,
        could not be read from this environment, so a period recorded only there would have been missed.</p>
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
    ])
    return f"""
<section id="hip56966" class="cand">
  <div class="cand-head">
    <p class="kicker">Candidate 3 · Track B, found during vetting · SPOC 2-minute data</p>
    <h2>HIP 56966: a sub-Neptune every 19.8 days</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("No SPOC TCE", "new")}</div>
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
        ("info", "Not yet ruled out:", "an eclipsing companion closer than ~1″, or two different planets each "
         "transiting once."),
    ])
    return f"""
<section id="hd11506" class="cand">
  <div class="cand-head">
    <p class="kicker">Candidate 2 · Track B (duo transit) · SPOC 2-minute data</p>
    <h2>HD 11506: a fourth planet in a known planetary system?</h2>
    <div class="chips">{chip("Not a TOI", "new")}{chip("No SPOC TCE", "new")}{chip("Host has 3 RV planets", "known")}</div>
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
                 453147896: "NGTS-39 b", 65612701: "TIC 65612701 b"}
# Confirmed planets without TOIs that Track B found blind (checked against the published periods), and one it lost.
RECOVERED = [("NGTS-38 b", "180.5 d, one transit"), ("TIC 279401253 b", "76.8 d, one transit"),
             ("TIC 4672985 b", "69.0 d, one transit"), ("NGTS-39 b", "58.2 d"),
             ("HD 148797 b and c", "42.1 and 68.2 d"), ("HD 3167 c", "29.85 d, one transit"),
             ("HD 278555 b", "28.84 d, period recovered"), ("TIC 65612701 b", "10.22 d, period recovered")]
MISSED = [("TIC 245076932 b", "21.6 d: all nine transits found, then rejected as an eclipsing binary")]


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
        kp = r.get("known_planet")
        return kp.split("; ")[0] if isinstance(kp, str) and kp else CONFIRMED_LIT.get(int(r.tic), f"TIC {int(r.tic)}")
    ps_txt = (f", and all {len(show)} were matched by position against NASA's confirmed-planet table of 5 October 2026"
              if "known_planet" in show else "")
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
        kp, st = published(r)
        mark = (" " + chip("published", "known") if kp or st == "known in this role" else
                " <span class='muted'>(named elsewhere)</span>" if st == "mentioned elsewhere" else "")
        if kp:
            mark += f"<div class='muted' style='font-size:0.78rem;margin-top:2px'>{esc(kp)}</div>"
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
  "published" marks a star already known in that role.{pub_txt}</p>
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


def section_trackc(tc, lit=None):
    if tc is None or not len(tc):
        return ""
    litm = {int(x["tic"]): x["status"] for x in (lit or [])}
    x = tc[(tc.P2_recurs >= 2)].copy()
    x["known"] = x.is_toi.astype(bool)
    x = x.sort_values(["P2_recurs", "snr2_max"], ascending=False).head(30)
    cols = [("Star", lambda r: f"<span class='m'>TIC {int(r.tic)}</span>", "l"),
            ("T mag", lambda r: fmt(r.Tmag, 1), "r"), ("Sectors", lambda r: esc(r.sectors), "l"),
            ("P<sub>1</sub> (d)", lambda r: fmt(r.P1, 4), "r"), ("Depth 1", lambda r: f"{r.depth1 * 100:.1f}%", "r"),
            ("P<sub>2</sub> (d)", lambda r: fmt(r.P2, 4), "r"), ("Depth 2", lambda r: f"{r.depth2 * 100:.2f}%", "r"),
            ("Recurs", lambda r: f"{int(r.P2_recurs)} sectors", "r"),
            ("Eclipse-like", lambda r: fmt(r.box_over_sine, 1), "r"),
            ("Known?", lambda r: (chip("TOI", "known") if r.known else
                                 chip("published", "known") if litm.get(int(r.tic)) == "known in this role" else
                                 "<span class='muted'>named, other role</span>" if litm.get(int(r.tic)) == "mentioned elsewhere" else
                                 "<span class='muted'>no mention found</span>" if int(r.tic) in litm else
                                 "<span class='muted'>not checked</span>") +
             (f"<div class='muted' style='font-size:0.78rem;margin-top:4px;max-width:20rem'>{esc(TRACKC_NOTES[int(r.tic)])}</div>"
              if int(r.tic) in TRACKC_NOTES else ""), "l")]
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
  period recurs in two or more sectors, the strongest evidence that it is real.</p>
  <p>Two caveats apply to every row. Most such pairs are two unrelated binaries blended in the same TESS pixels,
  which only a pixel-level analysis per eclipse set can tell apart; and the published catalogues of TESS quadruples
  and triples (Kostov et al. 2022, 2024 and others) could not be downloaded from this environment, so novelty is
  not yet checked against them. Web searches for the top sixteen found three already published: the sextuple
  TIC 168789840 (all three of its eclipse periods recovered blind), a catalogued quadruple, and a pulsating binary
  whose pulsation, folded at twice its period, passed our eclipse-shape test. That test compares the signal with a
  single sine wave and will be extended to harmonics.</p>
  </div>
  {table(x, cols)}
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
    n_lit = len(lit or [])
    n_known = sum(1 for x in (lit or []) if x["status"] == "known in this role")
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
        <li class="check check-info"><span class="ci" aria-hidden="true">i</span><div><strong>Lists we could not reach.</strong>
        ExoFOP (with its community TOIs and follow-up notes), SIMBAD, VizieR and arXiv were blocked from the compute environment. Known-object checks used the NASA TOI table of 5 October 2026, ExoFOP's TOI list of February 2026, SPOC's own detections, and web searches by TIC, Gaia and HIP number, which is how HIP 71639 turned up as a Planet Hunters TESS candidate. Of the {n_lit} stars this report names, {n_known} turned out to be published already in the same role, and are marked; a position match against NASA's confirmed-planet table caught several that web searches missed because the literature uses another name (HD 3167 for HIP 2736, HD 278555, HD 148797). Large candidate lists published as tables (the T16 hunt's 10,000 candidates from TESS's first year, the TESS Triple-9 catalogues, the eclipsing-quadruple lists) could not be matched row by row, and any candidate here could still be a community TOI that no indexed page names.</div></li>
        <li class="check check-info"><span class="ci" aria-hidden="true">i</span><div><strong>Candidates, not discoveries.</strong>
        Every object here is a candidate until follow-up observations rule out the impostors above; the search is complete
        for TESS sectors 1–107, and the vetting is deepest for stars brighter than T = 10.5.</div></li>
      </ul>
    </div>
    <div>
      <h3>Turning a candidate into a planet</h3>
      <p><strong>Report it.</strong> For HIP 71639, send the period to the Planet Hunters TESS team, whose
      unpublished radial velocities may already constrain the planet, and file it on ExoFOP next to their community
      TOI. For the others, ExoFOP accepts community TOIs (CTOIs) from registered users; a short Research Note of the
      AAS is the usual citable record for an amateur find.</p>
      <p><strong>Rule out impostors.</strong> Speckle or adaptive-optics imaging (both stars are bright enough for
      any facility) and a few high-precision radial velocities.</p>
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
        f56 = (f'<li>{chip("Not a TOI", "new")}<div><strong>HIP 56966: a {rp56:.1f} {RE} planet candidate every '
               f'{fit["HIP56966"]["period_d"][0]:.3f} days around a V = 8.8 Sun-like star.</strong> Four transits across '
               f'seven years, found while vetting the survey\'s single-transit list; SPOC never saw two in one search.</div></li>')
    known_tics = {int(x["tic"]) for x in (lit or []) if x["status"] == "known in this role"}
    if tb is not None and "known_planet" in tb:
        known_tics |= set(tb.tic[tb.known_planet.fillna("").astype(str) != ""].astype(int))
    n_dv = sum(1 for k, r in dv.items() if "error" not in r and deep_verdict(r)[0] == "candidate" and int(k) not in known_tics)
    n_dv_all = sum(1 for r in dv.values() if "error" not in r)
    tbf = trackb_shortlist(tb) if tb is not None else None
    n_tb = len(tbf) if tbf is not None else 0
    n_bright = int((tbf.Tmag < 10).sum()) if tbf is not None else 0
    n_tc = int((tc.P2_recurs >= 2).sum()) if tc is not None else 0
    plates = section_plates([(f"img/{p}", t, esc(x)) for p, t, x in PLATES])
    body = f"""<title>Ad Astra TESS Survey</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wdth,wght@100..125,400..800&family=IBM+Plex+Mono:wght@400;500;600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400..700;1,8..60,400&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header class="mast">
  <p class="kicker">Amateur search of NASA TESS data · survey report · {today}</p>
  <h1>Ad Astra</h1>
  <p class="dek">A blind search of {int(sv.n_real.sum()) / 1e6:.1f} million TESS light curves for transits and
  eclipses that no catalogue lists yet. Three bright, nearby stars show transits of Neptune-sized planets that
  are not TESS Objects of Interest; for one of them, a citizen-science single-transit candidate since 2021, the
  search supplies the missing period.</p>
  <div class="strip">
    <div><span class="label">Sectors searched</span><span class="v">{n_q + n_s} of 107</span></div>
    <div><span class="label">Light curves</span><span class="v">{int(sv.n_real.sum()):,}</span></div>
    <div><span class="label">Red-dwarf deep search</span><span class="v">{ta_n:,} of {ta_total:,}</span></div>
  </div>
</header>
<nav class="toc" aria-label="Sections"><a href="#findings">Findings</a><a href="#hip71639">HIP 71639</a>
<a href="#hd11506">HD 11506</a><a href="#hip56966">HIP 56966</a><a href="#bright">Bright candidates</a><a href="#trackb">All transit candidates</a><a href="#trackc">Multiple stars</a>
<a href="#trackd">Anomalies</a><a href="#calibration">Calibration</a><a href="#methods">Method</a><a href="#next">Next steps</a></nav>

<section id="findings">
  <p class="kicker">Summary</p>
  <h2>What the search found</h2>
  <ul class="findings">
    <li>{chip("New period", "new")}<div><strong>HIP 71639: a {rp71:.1f} {RE} planet candidate every
    {f71['period_d'][0]:.3f} days around a V = 6 F star.</strong> Planet Hunters TESS volunteers flagged one of its
    transits (a single-transit candidate since 2021, with no period). NASA's pipeline saw three by 2023 but filed
    them as one transit every 740 days, a period the transit's short duration does not fit. A fourth transit in
    April 2026 pins the true period: 740 days is exactly 31 orbits.</div></li>
    <li>{chip("Not a TOI", "new")}<div><strong>HD 11506: a {rp115:.1f} {RE} transiting candidate in a known
    three-planet system.</strong> Two matching transits, 2,149 days apart, around a star whose three giant planets
    were found by radial velocities. The period is 25.9 days or longer.</div></li>
    {f56}
    <li>{chip(f"{n_dv} bright", "new")}<div><strong>{n_dv} more candidates on bright stars survive one-by-one
    vetting.</strong> Of {n_dv_all} shortlisted stars brighter than T = 10.5, these pass the pixel source test, the
    period-alias test and Gaia's binarity checks, and no catalogue we could reach lists them. They include a warm
    giant on a 20.4-day orbit (HIP 63987), a grazing companion every 24.05 days (TIC 441010857) and a single
    16-hour transit of the A star HIP 109818. Two stars an earlier version of this report counted here are known
    planetary systems that never became TOIs, HD 3167 (HIP 2736) and HD 278555 (TIC 137966491). The full
    shortlist holds {n_tb} single and paired transits with no TOI.</div></li>
    <li>{chip(f"{n_tc} stars", "known")}<div><strong>Possible triple and quadruple stars.</strong> Eclipsing
    binaries with a second, independent set of eclipses that recurs across sectors. Most will be blends of two
    binaries; the catalogues to check novelty against were out of reach.</div></li>
    <li>{chip("Measured", "ok")}<div><strong>The false-alarm floor.</strong> Flipped light curves produced no
    planet-shaped dips in two-minute data. The survey recovers known planets, a sextuple star and a dipper blind.</div></li>
  </ul>
</section>
{section_hip71639(fit, ex)}
{section_hd11506(fit, ex)}
{section_hip56966(fit, ex)}
{section_deepvet(dv, lit, tb)}
{section_trackb(tb, singles, lit)}
<section id="plates" style="padding-top:8px">
  <h3>Selected dossiers</h3>
  <p class="muted col" style="font-size:var(--step--1)">Automatic dossiers: the dip with its fit (top left), the
  difference image with Gaia neighbours (top right), and every TESS sector of the star. Select one to enlarge.</p>
  {plates}
</section>
{section_trackc(tc, lit)}
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


def deep_verdict(r):
    g = r.get("gaia") or {}
    comb = (r.get("pixels") or {}).get("combined", {})
    worst = min(comb.items(), key=lambda kv: kv[1]["dchi2"]) if comb else None
    if worst and worst[1]["dchi2"] < -9:
        return "off-target", worst[0]
    if g.get("ruwe", 0) > 1.4 or g.get("nss", 0) > 0:
        return "binary host", None
    p = r.get("period") or {}
    n_ev = len([e for e in r.get("events", []) if e.get("sector") is not None])
    if ("aliases" in p and not any(a["status"] == "allowed" for a in p["aliases"])) or (n_ev >= 2 and not p):
        return "no consistent period", None
    if worst and worst[1]["dchi2"] < 9:
        return "ambiguous pixels", worst[0]
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
        comb = (r.get("pixels") or {}).get("combined", {})
        if not comb:
            pix = "no neighbour bright enough"
        else:
            w = min(comb.items(), key=lambda kv: kv[1]["dchi2"])
            pix = (f"on target (next best Δχ² +{w[1]['dchi2']:.0f})" if w[1]["dchi2"] >= 9 else
                   f"favours {w[0]} (Δχ² {w[1]['dchi2']:.0f})" if w[1]["dchi2"] <= -9 else
                   f"inconclusive ({w[0]}, Δχ² {w[1]['dchi2']:+.0f})")
        g = r.get("gaia") or {}
        verdict, _ = deep_verdict(r)
        lt = litm.get(tic, {}).get("status", "not checked")
        kp = b.known_planet.iloc[0] if len(b) and "known_planet" in b else None
        kp = kp if isinstance(kp, str) and kp else None
        if kp:
            lt = "confirmed planet"
        if verdict == "candidate" and lt in ("known in this role", "confirmed planet"):
            verdict = "published"
        note = NOTES.get(tic, "")
        if kp:
            note = f"NASA lists {kp}. " + note
        rows.append({"tic": tic, "hip": st.get("HIP"), "tmag": st.get("Tmag"), "rp": rp, "n_ev": len(ev),
                     "sectors": ", ".join(str(e["sector"]) for e in ev), "per": per, "pix": pix,
                     "ruwe": g.get("ruwe"), "nss": g.get("nss"), "tce": len((r.get("known") or {}).get("tces", [])),
                     "toi": len((r.get("known") or {}).get("toi", [])), "lit": lt, "verdict": verdict,
                     "note": note})
    df = pd.DataFrame(rows).sort_values(["verdict", "tmag"], key=lambda s: s.map({"candidate": 0}).fillna(1) if s.name == "verdict" else s)
    kinds = {"candidate": "new", "off-target": "eb", "binary host": "eb", "no consistent period": "known",
             "ambiguous pixels": "known", "published": "known"}

    def star(r):
        h = f" · HIP {int(r.hip)}" if r.hip == r.hip and r.hip else ""
        return f"<span class='m'>TIC {r.tic}</span>{h}"
    cols = [("Star", star, "l"), ("T", lambda r: fmt(r.tmag, 1), "r"), (f"Rp ({RE})", lambda r: fmt(r.rp, 1), "r"),
            ("Transits", lambda r: f"{r.n_ev} <span class='muted'>(S{esc(r.sectors)})</span>", "l"),
            ("Period", lambda r: esc(r.per), "l"), ("Pixels", lambda r: esc(r.pix), "l"),
            ("RUWE", lambda r: fmt(r.ruwe, 2), "r"), ("TCE", lambda r: str(r.tce), "r"),
            ("Literature", lambda r: esc(r.lit), "l"),
            ("Verdict", lambda r: chip(r.verdict, kinds.get(r.verdict, "known")) +
             (f"<div class='muted' style='font-size:0.78rem;margin-top:4px;max-width:22rem'>{esc(r.note)}</div>" if r.note else ""), "l")]
    n_c = int((df.verdict == "candidate").sum())
    counts = df.verdict.value_counts().to_dict()
    rejected = ", ".join(f"{v} {k}" for k, v in counts.items() if k != "candidate")
    df = df[df.verdict.isin(["candidate", "published"])].sort_values("rp")
    return f"""
<section id="bright">
  <p class="kicker">Track B · deep vetting</p>
  <h2>The brightest remaining candidates, checked one by one</h2>
  <div class="col" style="display:grid;gap:12px">
  <p>Every shortlisted star brighter than T = 10.5 went through the same checks as the headline candidates:
  Gaia binarity, TOI and SPOC lookups, the depth of each dip in every sector, the period aliases, and a
  pixel-level source test against every Gaia star within 150″ that could fake the dip with an eclipse of at
  most 50%. {n_c} of {sum(counts.values())} survive and are listed below, smallest first ({rejected}).
  HIP 56966 above came out of this list.</p>
  </div>
  {table(df, cols)}
  <p class="muted" style="font-size:var(--step--1)">Pixels: Δχ² of the best alternative source relative to the
  target, summed over transits (2-minute pixels where available, full-frame cutouts otherwise). Positive means
  the target fits better. RUWE above 1.4 hints at an unresolved binary.</p>
</section>"""


if __name__ == "__main__":
    build()
