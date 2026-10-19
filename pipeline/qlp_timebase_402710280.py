"""TIC 402710280.01: why the QLP time stamps of Sectors 100 and 102 are offset (qlp_recheck_402710280.py measures it).

For each QLP light curve, TIME is fitted with a straight line in CADENCENO. A time stamp carrying the target's own
barycentric correction departs from a straight line by as much as that correction changes across the sector (minutes);
one carrying a single correction for the whole sector does not depart at all. The rate at which the target's
correction changes at each transit (astropy, builtin ephemeris, geocentric observer) is compared with the drift of
the measured offset.

  python qlp_timebase_402710280.py
Output: results/candidates/qlp_timebase_402710280.json
"""
import json

import env  # noqa: F401
import numpy as np

from env import RESULTS
from qlp_recheck_402710280 import qlp_raw

RA, DEC = 248.5236, -63.5497


def bary_rate(t_btjd):
    """d(light-travel time to the barycentre)/dt for the target, in s per day, at BTJD t."""
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    c = SkyCoord(RA, DEC, unit="deg")
    loc = EarthLocation.from_geocentric(0, 0, 0, unit="m")
    ltt = [Time(t + 2457000, format="jd", scale="tdb", location=loc).light_travel_time(c, "barycentric").to(u.s).value
           for t in (t_btjd - 0.5, t_btjd + 0.5)]
    return float(ltt[1] - ltt[0])


def main():
    rec = json.load(open(RESULTS / "candidates" / "qlp_recheck_402710280.json"))["retime_qlp"]
    out = {}
    for s in (93, 100, 102):
        cad, t = qlp_raw(s)
        ok = np.isfinite(t)
        cad, t = cad[ok], t[ok]
        a = np.polyfit(cad.astype(float), t, 1)
        r = (t - np.polyval(a, cad.astype(float))) * 86400
        tc = rec[str(s)]["tc"]
        out[s] = {"n": int(t.size), "span_d": float(t[-1] - t[0]), "step_s": float(a[0] * 86400),
                  "departure_from_line_s": {"rms": float(r.std()), "range": float(np.ptp(r))},
                  "offset_drift_s_per_d": rec[str(s)]["time_offset"]["slope_s_per_d"],
                  "bary_rate_s_per_d": bary_rate(tc)}
        print(s, out[s])
    json.dump(out, open(RESULTS / "candidates" / "qlp_timebase_402710280.json", "w"), indent=1)


if __name__ == "__main__":
    main()
