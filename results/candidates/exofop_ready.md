# ExoFOP-ready values (12 October 2026)

Prepared from this survey's fits; BJD_TDB = BTJD + 2457000. Posting needs an ExoFOP account (the CTOI form or
"Upload planet parameters" on the target page). Nothing here has been posted.

## TIC 402710280.01 — new community TOI (CTOI)
| Field | Value |
|---|---|
| TIC | 402710280 (TYC 9041-03677-1, Gaia DR3 5828441289812582400) |
| Period (d) | 36.33972 ± 0.00007 |
| Transit epoch (BJD) | 2460474.8014 ± 0.0016 |
| Depth (ppm) | 1449 ± 229 |
| Duration (h) | 5.19 ± 0.18 |
| Rp/R* | 0.0381 ± 0.0029 |
| Impact parameter | 0.46 ± 0.32 |
| Planet radius (R_Earth) | 3.4 (R* = 0.817 R_Sun, TIC v8.2) |
| Transits used | S39 2384.607, S66 3111.411, S93 3838.214, S100 4092.595, S102 4128.936 (BTJD) |
| Disposition | PC |
Notes: five transits; S39 and S66 from TESS-SPOC, S93 from QLP, S100 and S102 from our own photometry of the
sector 100/102 full-frame images (no public FFI light curves of those sectors yet). Every predicted transit inside
the data is present. Difference-image centroids (five transits) 1.4" ± 1.2" from the target; the only unresolved
Gaia neighbour (G = 16.2 at 2.2", parallax 0.69 mas, background) would need a ~73% eclipse, excluded by the
transit shape (diluted-EB fits worse by dchi2 >= 65 for k >= 0.5). Stellar density from the transit shape
consistent with the K dwarf (S102 transit, on a ramp 1.4 d after an orbit start, excluded from the shape fit).
Gaia DR3 RUWE 1.00. No secondary eclipse; odd/even depths agree (2.3 sigma). Next transits: 2026-11-02 10:54, 2026-12-08 19:03, 2027-01-14 03:12, 2027-02-19 11:22 UTC.

## TOI-2493.01 — period for a single-transit TOI
| Field | Value |
|---|---|
| TIC | 123664207 |
| Period (d) | 79.0302 ± 0.0004 |
| Transit epoch (BJD) | 2461043.5093 ± 0.0042 |
| Depth (ppm) | 1467 ± 209 |
| Duration (h) | 6.10 ± 0.31 |
| Rp/R* | 0.0383 ± 0.0026 |
| Planet radius (R_Earth) | 2.7 (R* = 0.652 R_Sun) |
| Transits used | S33 2225.82, S98 4043.53, S107 4280.63 (BTJD), SPOC 2-min PDCSAP |
Notes: the S33 single transit (TOI comment "single at TBJD~2226") plus two more in sectors 98 and 107; every
shorter alias predicts transits inside sectors 6, 33, 87, 98 or 107 that are absent. Centroids (2-min pixels)
3.4" ± 2.0" from the target; fitted stellar density consistent with the K5 dwarf. Each transit is a 7-9 sigma
event. Insolation ~1.3 S_Earth. Next transits: 2026-11-15 03:07, 2027-02-02 03:50, 2027-04-22 04:34 UTC
(the first likely in TESS sector 110).

# Corrections to existing ExoFOP entries (12 October 2026)

Found by comparing the period hunt with periods already posted (`posted_check.py`, `posted_cases.py`,
`corrections_fit.py`); every transit was measured in the light curves and the transits were tested at the pixel
level. These are corrections to other people's entries: post them as planet-parameter updates on the target page
with a short note; for the TOI ones the TESS TOI team can also be told through ExoFOP.

## CTOI 438225366.01 — the period is 20.8615 d, not 41.7219 d
| Field | Value |
|---|---|
| TIC | 438225366 |
| Posted | P = 41.72193 ± 0.00222 d (gsun, 2022-04-30; "update period (Sectors 06,33,43,44,45)") |
| Period (d) | 20.861509 ± 0.000019 |
| Transit epoch (BJD) | 2459889.07987 ± 0.00048 |
| Depth (ppm) | 8253 ± 100 (fitted Rp/R*² ; box depth 7.7–8.5 ppt per transit) |
| Duration (h) | 7.04 ± 0.07 |
| Impact parameter | 0.75 ± 0.02 |
| Transits | BTJD 1470.513 (S6), 2221.520 (S33), 2492.718 (S43), 2534.440 (S45), 3243.725 (S71), 3264.586 (S72), 3285.448 (S72), 3681.813 (S87) |
Notes: 8 transits of equal depth; 4 of them (S43, S45, S71, S72 at 3285.45) fall at odd cycles of the 41.72-d period,
where it allows none. Odd/even depths 7773 vs 8026 ppm (2.0 sigma); no eclipse at phase 0.5 (58 ± 70 ppm).
Centroids on target (1.2"). Stellar density from the shape 0.178 ± 0.014 solar vs TIC 0.199. The host is an A star
(8260 K, 2.17 R_Sun, 650 pc): Rp/R* = 0.091 means about 1.9 R_Jup, so a low-mass star is more likely than a planet.

## TOI-2633.01 — the listed ephemeris misses the transits
| Field | Value |
|---|---|
| TIC | 386131764 (CTOI 386131764.01, steuer) |
| Listed | T0 = BJD 2459229.361575, P = 18.1960491 ± 0.0032 d (TESS project, 2022-07-06) |
| Period (d) | 18.214811 ± 0.000030 |
| Transit epoch (BJD) | 2459964.00532 ± 0.0014 |
| Depth (ppm) | 3740 ± 140 |
| Duration (h) | 4.62 ± 0.09 |
| Transits | BTJD 1506.771 (S7, QLP), 2235.367 and 2253.582 (S34), 3674.344 (S87), 3692.559 and 3710.774 (S88) |
Notes: the light curve is flat at the listed epoch and at 5 of the 6 other times the listed ephemeris predicts
(S7, S34, S87, S88, S108). The transits follow the community TOI's own epoch (BJD 2459253.5816). Pixel test:
target or a G = 17.6 star 7" away (unresolved; that star would need a ~90% eclipse). The transit implies a stellar
density of 1.2 vs the TIC's 0.40 solar: eccentric orbit, a smaller star, or a blend.

## CTOI 156716001.01 is not TOI-7602.01
| Field | Value |
|---|---|
| TIC | 156716001 |
| ExoFOP | CTOI (hawthorn 2023, epoch BJD 2458485.406; essack 2024, epoch 2459224.4578) merged into TOI-7602.01 (QLP faint search, P = 10.29779 d, depth 794 ppm, 3.8 h) |
| Its transits (BJD) | 2458485.406 (S6), 2459224.525 (S33), 2460702.770 (S88, new) |
| Depth (ppm) | 2303 (+124/-78) |
| Duration (h) | 8.9 ± 0.15 |
| Periods allowed | 739.12/n d, n = 1-10 (>= 67 d for a circular orbit); stellar density favours n = 6-7 (123.2 or 105.6 d) |
Notes: the S6 transit is 24 h (96 sigma) off TOI-7602.01's ephemeris; TOI-7602.01's own shallow transits are in the
data at other times (one of them a day later in S6). All three transits on target (1.3"). A second transiting
candidate on the star.

## CTOI 369532319.01 is not TOI-5967.01
| Field | Value |
|---|---|
| TIC | 369532319 |
| ExoFOP | PHT CTOI (eisner 2020, single transit BJD 2458755.81) merged into TOI-5967.01 (P = 15.343585 d) |
| Its transits (BJD) | 2458755.810 (S16), 2459831.560 (S56, new) |
| Depth (ppm) | 1228 (+194/-120) |
| Duration (h) | 5.8 ± 0.4 |
| Periods allowed | 1075.75/n d, about 30 values >= 20 d; density too uncertain to choose |
Notes: the S16 transit is 2.7 d after a TOI-5967.01 transit (65 h, 306 sigma off its ephemeris), and the S56 one
4.4 d after another; TOI-5967.01's own transits are all present. Both CTOI transits on target. Each is a 9-10 sigma
event.

## CTOI 67711670.01 is not TOI-6033.01
| Field | Value |
|---|---|
| TIC | 67711670 |
| ExoFOP | TSTPC single transit (dragomir 2022, BJD 2458798.9854) merged into TOI-6033.01 (SPOC S58, P = 18.04015 ± 0.0017 d) |
| Its transits (BJD) | 2458798.985 (S18), 2460622.800 (S85, new) |
| Depth (ppm) | 6136 ± 230 |
| Duration (h) | 6.4 ± 0.25 |
| Periods allowed | 1823.80/n d (35 values >= 12 d); stellar density favours about 57-63 d (62.89 d, n = 29, is allowed) |
Notes: the S18 transit is 110 h (45 sigma) off TOI-6033.01's ephemeris and the S85 one 2.8 d before a predicted
TOI-6033.01 transit; no single period fits all four transits. Both on target (4.5").

## TOI ephemerides that no longer find their transits (toi_audit.py, toi_audit_review.py)
All planet-candidate TOIs (PC/APC) were checked against every TESS sector. The corrected ephemerides are in
`toi_audit_corrections.csv` (one row per flagged TOI): columns `P`, `P_err`, `T0` (BTJD; BJD = BTJD + 2457000),
`T0_err`, `status` ("corrected": the listed ephemeris is off on 1 Dec 2026 by more than a transit duration and more
than 3 sigma of its own error; "refined": off by more than a duration but within its error), `solid` (>= 5 transits
seen, >= 70% of those covered, chi2_red < 5), `already_posted` (a matching period already uploaded to ExoFOP by
someone), `offset_h` (listed minus corrected transit time on 1 Dec 2026). Rows with `odd_even_sigma` > 3 may be
eclipsing binaries at twice the period.
