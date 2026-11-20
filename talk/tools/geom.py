import sys, csv, json, math, numpy as np
D = sys.argv[1]; IMG = sys.argv[2]
K, B, ARS, P, U1, U2 = 0.04209, 0.2243, 53.63, 36.33986, 0.39272, 0.20682

def flux(z):
    """Quadratic limb-darkened transit flux at sky separation z (stellar radii), by direct integration."""
    if z >= 1 + K: return 1.0
    n = 60
    r = (np.arange(n) + 0.5) / n * K
    th = (np.arange(2 * n) + 0.5) / (2 * n) * 2 * np.pi
    R, TH = np.meshgrid(r, th)
    x = z + R * np.cos(TH); y = R * np.sin(TH)
    d2 = x * x + y * y
    ins = d2 < 1
    mu = np.sqrt(np.clip(1 - d2, 0, 1))
    I = (1 - U1 * (1 - mu) - U2 * (1 - mu) ** 2) * ins
    blocked = np.sum(I * R) * (K / n) * (2 * np.pi / (2 * n))
    total = np.pi * (1 - U1 / 3 - U2 / 6)
    return 1 - blocked / total

def model_ppm(t_h):
    zx = ARS * 2 * np.pi * (t_h / 24) / P
    return (flux(math.hypot(B, zx)) - 1) * 1e6

def load(s):
    rows = list(csv.DictReader(open(f"{D}/tic402710280_s{s:03d}_window.csv")))
    t = np.array([float(r["bjd_tdb"]) for r in rows]); f = np.array([float(r["flux_normalised"]) for r in rows])
    tc = float(rows[0]["centre_bjd_tdb"])
    return (t - tc) * 24, (f - 1) * 1e6

def binned(t, y, w, lo=-8, hi=8):
    edges = np.arange(lo, hi + 1e-9, w); out = []
    for a in edges[:-1]:
        m = (t >= a) & (t < a + w)
        if m.sum() >= 1: out.append((a + w / 2, float(np.mean(y[m]))))
    return out

W, H, YT, YB = 400, 260, 1400.0, -3400.0
X = lambda t: (t + 8) / 16 * W
Y = lambda p: (YT - p) / (YT - YB) * H
def dots(pts):
    return "".join(f"M{X(t):.1f} {Y(p):.1f}h0" for t, p in pts if YB < p < YT)

out = {}
allpts = []
for s in (39, 93, 100, 102):
    t, y = load(s)
    out[f"faint{s}"] = dots(binned(t, y, 1 / 6))
    out[f"bins{s}"] = dots(binned(t, y, 0.5))
    allpts.append((t, y))
    print(s, len(t), "depth@0 10-min mean", round(float(np.mean(y[np.abs(t) < 1.5]))))
tm = np.arange(-8, 8.001, 0.125)
out["model"] = "M" + "L".join(f"{X(a):.1f} {Y(model_ppm(a)):.1f}" for a in tm)
out["y0"] = round(Y(0), 1); out["ydep"] = round(Y(model_ppm(0)), 1); print("model depth ppm", model_ppm(0))

# folded curve for the cover: 15-minute bins of the four transits, +-9 h (window data span +-19 h)
T = np.concatenate([a for a, _ in allpts]); YY = np.concatenate([b for _, b in allpts])
CW, CH = 1664, 300
fb = binned(T, YY, 0.25, -9, 9)
cx = lambda t: (t + 9) / 18 * CW
cy = lambda p: (700 - p) / (700 + 2700) * CH
out["cover"] = "".join(f"M{cx(t):.1f} {cy(p):.1f}h0" for t, p in fb)
out["cover_model"] = "M" + "L".join(f"{cx(a):.1f} {cy(model_ppm(a)):.1f}" for a in np.arange(-9, 9.001, 0.125))

# transit diagram: star radius 220 px; planet x = zx * R; curve below with the dip exaggerated to 80 px
R = 220
xs = np.arange(-400, 401, 8)
def dip(xpx):
    z = math.hypot(B, xpx / R); return (1 - flux(z)) / (1 - flux(B))
out["diagram_curve"] = "M" + "L".join(f"{x:.0f} {80 * dip(x):.1f}" for x in xs)
out["planet_r"] = round(K * R, 1); out["chord_y"] = round(B * R, 1)
print("ingress x px", R * math.sqrt((1 + K) ** 2 - B ** 2))

# pixels: 13x13 median out-of-transit image, asinh stretch, 14 levels, one path per level (row 0 at the bottom)
img = np.load(IMG)
v = np.arcsinh((img - np.percentile(img, 5)) / 150.0); v = np.clip(v / v.max(), 0, 1)
lv = np.minimum((v * 14).astype(int), 13)
c0, c1 = np.array([0x1b, 0x19, 0x18]), np.array([0xf7, 0xf0, 0xe3])
px = 48; paths = []
for L in range(14):
    d = "".join(f"M{c * px} {(12 - r) * px}h{px}v{px}h-{px}z" for r in range(13) for c in range(13) if lv[r, c] == L)
    if d:
        col = c0 + (c1 - c0) * (L / 13)
        paths.append(f'<path fill="#{int(col[0]):02x}{int(col[1]):02x}{int(col[2]):02x}" d="{d}"/>')
out["pixels"] = "".join(paths)
json.dump(out, open(sys.argv[3], "w"), indent=0)
for k, val in out.items(): print(k, len(str(val)))
