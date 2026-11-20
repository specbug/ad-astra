import sys, numpy as np
a = np.load(sys.argv[1])[0:18, 0:32]            # rows 0-17, cols 0-31; star at row 9, col 21
lo = 201.0
v = np.arcsinh(np.clip(a - lo, 0, None) / 150.0)
v = v / v.max()
c0 = np.array([0x13, 0x11, 0x10]); c1 = np.array([0xf7, 0xf0, 0xe3])
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1920 1080" width="1920" height="1080" shape-rendering="crispEdges">',
       '<path fill="#131110" d="M0 0h1920v1080H0z"/>']
by = {}
for r in range(18):
    for c in range(32):
        q = round(float(v[r, c]) * 40) / 40          # 41 levels
        if q == 0: continue
        rgb = np.round(c0 + (c1 - c0) * q).astype(int)
        by.setdefault('#%02x%02x%02x' % tuple(rgb), []).append(f'M{c*60} {(17-r)*60}h60v60h-60z')
for k, d in sorted(by.items()):
    out.append(f'<path fill="{k}" d="{"".join(d)}"/>')
out.append('</svg>')
open(sys.argv[2], 'w').write(''.join(out))
print(len(''.join(out)), len(by))
