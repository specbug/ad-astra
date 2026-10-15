# Request for association review: three community TOIs merged into TOIs

Draft note to ExoFOP staff (exofop-support@ipac.caltech.edu), 15 October 2026. Evidence: pipeline/merged_ctoi_test.py and results/candidates/merged_ctoi_test.json; light-curve strips in the report's Corrections section. Not to be sent as parameter uploads: the TOIs' own values are not in question here, only whether the community TOIs belong to them.

---

Dear ExoFOP team,

Three single-transit community TOIs are listed as promoted to TOIs on the same stars, but in the TESS light curves their dips do not belong to those TOIs. For each star I measured every dip the same way (trapezoid fit for mid-time, depth and duration, quadratic local baseline, red-noise-scaled errors), and tested every period that could put all the dips on one ephemeris against the rest of the data. Could you review the associations? I am not asking for any TOI parameter to change.

**CTOI 156716001.01 (TIC 156716001) → TOI-7602.01.** The community TOI's dips (S6 BTJD 1485.415, S33 BTJD 2224.516, S88 BTJD 3702.773): 2454 ± 110 ppm, duration 10.7 h. TOI-7602.01's transits on its catalogue ephemeris (12 measured): 805 ± 48 ppm, median duration 5.9 h (individual fits 3.2–9.5 h). No single period puts all the dips on one ephemeris without predicting transits that the data do not show (no period fits all the dip times).

**CTOI 369532319.01 (TIC 369532319) → TOI-5967.01.** The community TOI's dips (S16 BTJD 1755.777, S56 BTJD 2831.744): 1209 ± 105 ppm, median duration 11.3 h (individual fits 8.1–14.5 h). TOI-5967.01's transits on its catalogue ephemeris (4 measured): 1367 ± 81 ppm, median duration 4.7 h (individual fits 3.9–4.7 h). No single period puts all the dips on one ephemeris without predicting transits that the data do not show (no period fits all the dip times). TOI-5967.01's own transit is present in sector 16 at its catalogue ephemeris (BTJD 1753.108 ± 0.009, 914 ± 217 ppm), 2.7 d before the community TOI's dip.

**CTOI 67711670.01 (TIC 67711670) → TOI-6033.01.** The community TOI's dips (S18 BTJD 1799.034, S85 BTJD 3622.805): 4408 ± 282 ppm, median duration 8.6 h (individual fits 6.4–10.9 h). TOI-6033.01's transits on its catalogue ephemeris (2 measured): 5497 ± 107 ppm, duration 6.0 h. No single period puts all the dips on one ephemeris without predicting transits that the data do not show (the one period that fits all four times, 0.564 d, predicts 107 other transits, none of which is there). Separately, TOI-6033.01's 18.0401-d period (two transits in sector 58) predicts a transit in sector 18 near BTJD 1803.57 that is not seen (≤ ~500 ppm within ±12 h, where ~5500 ppm is expected), so that TOI's period may itself need a look.

All three community-TOI signals are on the target in difference-image centroids (each within its errors), though that does not exclude a source closer than a few arcseconds. I can send the per-event measurements and light-curve plots.

Thank you,
Rishit Vora
