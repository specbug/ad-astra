# Animation for the blog post: TIC 402710280.01 crossing its star, with the light curve traced underneath.
# Planet from the preprint's posterior medians, drawn to scale. Same visual language as blog_eb_animation.py.
# Writes docs/blog/transit.{mp4,gif} and transit-final.png.
# Usage: python blog_transit_animation.py [font_dir]   (font_dir holding DM Sans / JetBrains Mono .ttf; optional)
import sys, glob, math, subprocess, numpy as np, matplotlib, batman
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.patches import Circle

R = '/home/user/ad-astra/'
OUT = R + 'docs/blog/transit'
P = 36.33986

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


PAPER, INK, INK2, RULE = oklch(.97, .012, 85), oklch(.22, .02, 270), oklch(.50, .015, 270), oklch(.86, .01, 90)
STAR_C, STAR_L = np.array(oklch(.95, .05, 85)), np.array(oklch(.80, .11, 62))

# planet (preprint posterior medians)
PK, PB, PA, PU = 0.04209, 0.224, 53.63, (0.393, 0.207)
TW = 4.5  # hours either side of mid-transit
th = np.linspace(-TW, TW, 3601)
bp = batman.TransitParams(); bp.t0 = 0; bp.per = P; bp.rp = PK; bp.a = PA; bp.inc = math.degrees(math.acos(PB / PA))
bp.ecc = 0; bp.w = 90; bp.u = list(PU); bp.limb_dark = 'quadratic'
yp = -(1 - batman.TransitModel(bp, th / 24).light_curve(bp)) * 1e6
pos = lambda t: PA * math.sin(2 * math.pi * t / 24 / P)
t14 = np.ptp(th[yp < 0])  # first to last contact, hours


def disc(c0, c1, u, n=600):
    y, x = np.mgrid[-1:1:n * 1j, -1:1:n * 1j]; r2 = x * x + y * y
    mu = np.sqrt(np.clip(1 - r2, 0, 1)); I = 1 - u[0] * (1 - mu) - u[1] * (1 - mu) ** 2
    w = ((I - I[r2 <= 1].min()) / (1 - I[r2 <= 1].min()))[..., None]
    img = np.concatenate([c1 + (c0 - c1) * w, np.ones_like(w)], axis=2)
    img[..., 3] = np.clip((1 - np.sqrt(r2)) * n / 2, 0, 1)
    return img


plt.rcParams.update({'font.family': SANS, 'text.color': INK, 'axes.facecolor': PAPER, 'figure.facecolor': PAPER, 'savefig.facecolor': PAPER})
fig = plt.figure(figsize=(8.0, 8.0), dpi=120)
L, Wd = 0.09, 0.84

# sky
XL = 2.15
sky = fig.add_axes([L, 0.50, Wd, 0.40]); sky.set_xlim(-XL, XL); sky.set_ylim(-1.15, 1.15); sky.set_aspect('equal'); sky.axis('off')
sky.imshow(disc(STAR_C, STAR_L, PU), extent=(-1, 1, -1, 1), zorder=1, interpolation='bilinear')
planet = Circle((pos(-TW), -PB), PK, color=INK, zorder=3); sky.add_patch(planet)
fig.text(L, 0.935, 'TIC 402710280.01  ·  TO SCALE', family=MONO, fontsize=9, color=INK2, va='center')
fig.text(L + Wd - 0.085, 0.935, 'light lost', fontsize=9, color=INK2, ha='right', va='center')
lost = fig.text(L + Wd, 0.935, '', family=MONO, fontsize=10, color=INK, ha='right', va='center')

# light curve
ax = fig.add_axes([L, 0.10, Wd, 0.32])
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlim(-TW, TW); ax.set_ylim(-2650, 300)
ax.axhline(0, color=RULE, lw=0.9, zorder=0)
for y in (-1000, -2000): ax.axhline(y, color=RULE, lw=0.6, ls=(0, (1, 3)), zorder=0)
ax.set_yticks([0, -1000, -2000]); ax.set_yticklabels(['0', '−1000', '−2000'])
ax.set_xticks(range(-4, 5, 2)); ax.set_xticklabels([f'{v:+d} h'.replace('-', '−') if v else '0' for v in range(-4, 5, 2)])
ax.tick_params(length=0, labelsize=9, colors=INK2, pad=6)
for t in ax.get_xticklabels() + ax.get_yticklabels(): t.set_family(MONO)
fig.text(L, 0.455, 'BRIGHTNESS  ·  PPM', family=MONO, fontsize=9, color=INK2, va='center')
clock = fig.text(L + Wd, 0.455, '', family=MONO, fontsize=9, color=INK2, ha='right', va='center')
ax.plot(th, yp, color=INK, lw=1, alpha=0.10, zorder=1)
LP, = ax.plot([], [], color=INK, lw=2.0, zorder=3, solid_capstyle='round')
DP, = ax.plot([], [], 'o', ms=5.5, color=INK, zorder=4, clip_on=False)
# closing annotation: the duration, first to last contact
h = t14 / 2
span = [ax.plot([-h, h], [-2450, -2450], color=INK2, lw=0.9, alpha=0)[0],
        ax.plot([-h, -h], [-2510, -2390], color=INK2, lw=0.9, alpha=0)[0],
        ax.plot([h, h], [-2510, -2390], color=INK2, lw=0.9, alpha=0)[0]]
span_t = ax.text(0, -2380, f'{t14:.1f} hours', family=MONO, fontsize=9, color=INK2, ha='center', va='bottom', alpha=0)

FPS, MOVE, HOLD0, HOLD1, FADE = 30, 6.5, 0.6, 2.4, 0.5
n0, nm, n1 = int(HOLD0 * FPS), int(MOVE * FPS), int(HOLD1 * FPS)
NF = n0 + nm + n1


def frame(i):
    t = -TW + 2 * TW * min(max(i - n0, 0), nm) / nm
    m = th <= t
    LP.set_data(th[m], yp[m])
    d = np.interp(t, th, yp)
    DP.set_data([t], [d])
    planet.center = (pos(t), -PB)
    lost.set_text(f'{-d / 1e4:4.2f} %')
    clock.set_text(f't = {t:+4.1f} h'.replace('-', '−'))
    k = (i - n0 - nm) / (FADE * FPS)
    a = 0 if k <= 0 else 1 - (1 - min(k, 1)) ** 3
    for ln in span: ln.set_alpha(a)
    span_t.set_alpha(a)
    return ()


anim = FuncAnimation(fig, frame, frames=NF, blit=False)
anim.save(OUT + '.mp4', writer=FFMpegWriter(fps=FPS, codec='libx264', bitrate=2500, extra_args=['-pix_fmt', 'yuv420p', '-preset', 'slow']))
frame(NF - 1); fig.savefig(OUT + '-final.png', dpi=120)
subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', OUT + '.mp4', '-vf',
                'fps=20,scale=800:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=64:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle',
                OUT + '.gif'], check=True)
print(f'depth={-yp.min():.0f} ppm  T14={t14:.2f} h')
