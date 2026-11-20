import sys, json, numpy as np
from astropy.io import fits
f = fits.open(sys.argv[1])
d = f[1].data
t, flux, q = d["TIME"], d["FLUX"], d["QUALITY"]
ok = (q == 0) & np.isfinite(t)
tc = 2459384.6094 - 2457000
oot = ok & (np.abs(t - tc) > 0.25) & (np.abs(t - tc) < 1.5)
img = np.nanmedian(flux[oot], axis=0)
print(img.shape, float(img.min()), float(img.max()), oot.sum())
np.save(sys.argv[2], img)
h = f[1].header
print({k: h.get(k) for k in ("1CRV4P","2CRV4P","1CRPX4","2CRPX4")}, f[0].header.get("SECTOR"), f[0].header.get("CAMERA"), f[0].header.get("CCD"))
print(np.round(img).astype(int))
