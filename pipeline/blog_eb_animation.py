# Animation for the blog post: why a diluted eclipse of the 3.5" neighbour by an equal-sized star gives a V-shaped dip,
# while the planet gives sharp edges and a flat-ish bottom. Planet: the preprint's posterior medians. Eclipse:
# illustrative k = 1 pair whose eclipse removes 74% of the neighbour's light (b solved for that depth) and whose
# half-depth width matches the planet's (a/R solved), diluted by 1/354. Writes docs/blog/eb-vs-planet.{mp4,gif}.
# Usage: python blog_eb_animation.py [font_dir]   (font_dir holding DM Sans / JetBrains Mono .ttf; optional)
import sys, glob, math, subprocess, numpy as np, matplotlib, batman
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.patches import Circle
from scipy.optimize import brentq

R = '/home/user/ad-astra/'
OUT = R + 'docs/blog/eb-vs-planet'
P = 36.33986
DIL = 354.0

if len(sys.argv) > 1:
    for f in glob.glob(sys.argv[1] + '/*.ttf'):
        font_manager.fontManager.addfont(f)
names = {f.name for f in font_manager.fontManager.ttflist}
SANS = 'DM Sans' if 'DM Sans' in names else 'DejaVu Sans'
MONO = 'JetBrains Mono' if 'JetBrains Mono' in names else 'DejaVu Sans Mono'


def oklch(L, C, h):
    a, b = C * math.cos(math.radians(h)), C * math.sin(math.radians(h))
    l_, m_, s_ = L + 0.3963377774 * a + 0.2158037573 * b, L - 0.1055613458 * a - 0.0638541728 * b, L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    rgb = (4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
           -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
           -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)
    enc = lambda c: 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055
    return tuple(min(1, max(0, enc(c))) for c in rgb)


PAPER, INK, INK2, RULE, ACCENT = oklch(.97, .012, 85), oklch(.22, .02, 270), oklch(.50, .015, 270), oklch(.86, .01, 90), oklch(.62, .18, 35)
STAR_C, STAR_L = np.array(oklch(.95, .05, 85)), np.array(oklch(.80, .11, 62))
COMP_C, COMP_L = np.array(oklch(.86, .025, 250)), np.array(oklch(.70, .04, 250))


def model(k, b, a, u, t_h):
    p = batman.TransitParams(); p.t0 = 0; p.per = P; p.rp = k; p.a = a; p.inc = math.degrees(math.acos(b / a))
    p.ecc = 0; p.w = 90; p.u = list(u); p.limb_dark = 'quadratic'
    return 1 - batman.TransitModel(p, np.atleast_1d(t_h) / 24).light_curve(p)


# planet (preprint posterior medians)
PK, PB, PA, PU = 0.04209, 0.224, 53.63, (0.393, 0.207)
# eclipse: equal-sized pair, 74% of the neighbour's light at mid-eclipse, same half-depth width as the planet
EK, EU = 1.0, (0.4, 0.25)
th = np.linspace(-6, 6, 4801)
hw = lambda d: np.ptp(th[d > d.max() / 2])
EB = brentq(lambda b: model(EK, b, 60, EU, 0.0)[0] - 0.74, 0.2, 1.5)
EA = brentq(lambda a: hw(model(EK, EB, a, EU, th)) - hw(model(PK, PB, PA, PU, th)), 20, 200)
lp, le = model(PK, PB, PA, PU, th), model(EK, EB, EA, EU, th)
pos = lambda a, t: a * math.sin(2 * math.pi * t / 24 / P)


def disc(c0, c1, u, n=420):
    y, x = np.mgrid[-1:1:n * 1j, -1:1:n * 1j]; r2 = x * x + y * y
    mu = np.sqrt(np.clip(1 - r2, 0, 1)); I = 1 - u[0] * (1 - mu) - u[1] * (1 - mu) ** 2
    w = ((I - I[r2 <= 1].min()) / (1 - I[r2 <= 1].min()))[..., None]
    img = np.concatenate([c1 + (c0 - c1) * w, np.ones_like(w)], axis=2)
    edge = np.clip((1 - np.sqrt(r2)) * n / 2, 0, 1)  # one-pixel antialiased rim
    img[..., 3] = edge
    return img


plt.rcParams.update({'font.family': SANS, 'text.color': INK, 'axes.facecolor': PAPER, 'figure.facecolor': PAPER, 'savefig.facecolor': PAPER})
W, H = 12.0, 6.4
fig = plt.figure(figsize=(W, H), dpi=120)
star_img, comp_img = disc(STAR_C, STAR_L, PU), disc(COMP_C, COMP_L, EU)
XL, YL = 3.4, 1.62


def sky(y0, label, col):
    ax = fig.add_axes([0.045, y0, 0.40, 0.33]); ax.set_xlim(-XL, XL); ax.set_ylim(-YL, YL); ax.set_aspect('equal'); ax.axis('off')
    ax.imshow(star_img, extent=(-1, 1, -1, 1), zorder=1, interpolation='bilinear')
    fig.lines.append(plt.Line2D([0.05, 0.068], [y0 + 0.355] * 2, transform=fig.transFigure, color=col, lw=2.2, solid_capstyle='butt'))
    fig.text(0.076, y0 + 0.355, label, family=MONO, fontsize=8.6, color=INK2, va='center')
    return ax


ax1 = sky(0.555, 'PLANET  ·  CROSSING MY STAR', INK)
ax2 = sky(0.105, 'EQUAL-SIZED STAR  ·  ECLIPSING THE NEIGHBOUR', ACCENT)
planet = Circle((pos(PA, -6), -PB), PK, color=INK, zorder=3); ax1.add_patch(planet)
comp = ax2.imshow(comp_img, extent=(-1, 1, -1, 1), zorder=3, interpolation='bilinear')
c1 = fig.text(0.445, 0.86, '', family=MONO, fontsize=9.5, color=INK, ha='right', va='center')
c2 = fig.text(0.445, 0.41, '', family=MONO, fontsize=9.5, color=ACCENT, ha='right', va='center')
fig.text(0.375, 0.86, 'light lost', fontsize=8.6, color=INK2, ha='right', va='center')
fig.text(0.375, 0.41, 'light lost', fontsize=8.6, color=INK2, ha='right', va='center')

ax = fig.add_axes([0.53, 0.25, 0.43, 0.62])
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlim(-6, 6); ax.set_ylim(-2550, 350)
ax.axhline(0, color=RULE, lw=0.9, zorder=0)
for y in (-1000, -2000): ax.axhline(y, color=RULE, lw=0.6, ls=(0, (1, 3)), zorder=0)
ax.set_yticks([0, -1000, -2000]); ax.set_yticklabels(['0', '−1000', '−2000'])
ax.set_xticks(range(-6, 7, 2)); ax.set_xticklabels([f'{v:+d} h'.replace('-', '−') if v else '0' for v in range(-6, 7, 2)])
ax.tick_params(length=0, labelsize=8.6, colors=INK2, pad=6)
for t in ax.get_xticklabels() + ax.get_yticklabels(): t.set_family(MONO)
fig.text(0.53, 0.905, 'BRIGHTNESS CHANGE IN THE COMBINED LIGHT  ·  PPM', family=MONO, fontsize=8.6, color=INK2)
yp, ye = -lp * 1e6, -le / DIL * 1e6
ax.plot(th, yp, color=INK, lw=1, alpha=0.10, zorder=1); ax.plot(th, ye, color=ACCENT, lw=1, alpha=0.16, zorder=1)
LP, = ax.plot([], [], color=INK, lw=2.0, zorder=3, solid_capstyle='round')
LE, = ax.plot([], [], color=ACCENT, lw=2.0, zorder=2, solid_capstyle='round')
DP, = ax.plot([], [], 'o', ms=5.5, color=INK, zorder=4)
DE, = ax.plot([], [], 'o', ms=5.5, color=ACCENT, zorder=4)
clock = fig.text(0.96, 0.905, '', family=MONO, fontsize=8.6, color=INK2, ha='right')
fig.text(0.53, 0.115, 'The neighbour supplies 1/354 of the light, so its 74% eclipse\nshrinks to the same 0.21% dip. Only the shape gives it away.',
         fontsize=10, color=INK2, linespacing=1.5, va='top')
n1 = ax.text(-2.95, -1350, 'sharp edge:\n14 min to ingress', family=MONO, fontsize=8.2, color=INK, ha='right', va='center', alpha=0)
n2 = ax.text(3.25, -1150, 'deepens for 5 h,\nthen recovers', family=MONO, fontsize=8.2, color=ACCENT, ha='left', va='center', alpha=0)

FPS, MOVE, HOLD0, HOLD1, FADE = 30, 7.0, 0.6, 2.4, 0.5
n0, nm, n1f = int(HOLD0 * FPS), int(MOVE * FPS), int(HOLD1 * FPS)
NF = n0 + nm + n1f


def frame(i):
    t = -6 + 12 * min(max(i - n0, 0), nm) / nm
    m = th <= t
    LP.set_data(th[m], yp[m]); LE.set_data(th[m], ye[m])
    dp, de = np.interp(t, th, yp), np.interp(t, th, ye)
    DP.set_data([t], [dp]); DE.set_data([t], [de])
    planet.center = (pos(PA, t), -PB)
    x = pos(EA, t); comp.set_extent((x - 1, x + 1, -EB - 1, -EB + 1))
    c1.set_text(f'{-dp / 1e4:5.2f} %'); c2.set_text(f'{-de * DIL / 1e4:5.0f} %')
    clock.set_text(f't = {t:+5.1f} h'.replace('-', '−').replace('+', '+'))
    k = (i - n0 - nm) / (FADE * FPS)
    ease = 0 if k <= 0 else 1 - (1 - min(k, 1)) ** 3
    n1.set_alpha(ease); n2.set_alpha(ease)
    return ()


anim = FuncAnimation(fig, frame, frames=NF, blit=False)
anim.save(OUT + '.mp4', writer=FFMpegWriter(fps=FPS, codec='libx264', bitrate=3000, extra_args=['-pix_fmt', 'yuv420p', '-preset', 'slow']))
frame(NF - 1); fig.savefig(OUT + '-final.png', dpi=120)
subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', OUT + '.mp4', '-vf',
                'fps=20,scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=96:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle',
                OUT + '.gif'], check=True)
print(f'b={EB:.3f} a/R={EA:.2f} planet depth={lp.max()*1e6:.0f} ppm eclipse={le.max():.3f} diluted={le.max()/DIL*1e6:.0f} ppm')
