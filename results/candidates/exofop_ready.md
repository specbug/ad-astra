# ExoFOP-ready values (12 October 2026)

Prepared from this survey's fits; BJD_TDB = BTJD + 2457000. Posting needs an ExoFOP account (the CTOI form or
"Upload planet parameters" on the target page). Nothing here has been posted.

## TIC 402710280.01 — new candidate (NOT uploadable yet: ExoFOP requires prior publication)
Revised 15 October 2026 after the independent audit (pipeline/audit_402710280.py, blend_shape_402710280.py,
fit_new.py). Since 19 August 2026 ExoFOP accepts new community candidates only after peer-reviewed publication
and an approved Published Candidate Upload Request; these values are for the paper and a later upload.

| Field | Value |
|---|---|
| TIC | 402710280 (TYC 9041-03677-1, Gaia DR3 5828441289812582400) |
| Period (d) | 36.33990 ± 0.00053 (errors inflated by sqrt(chi2_red): the mid-times scatter, see below) |
| Transit epoch (BJD_TDB) | 2460838.2019 ± 0.0122 (the sector-93 transit) |
| Depth (ppm) | 2097 ± 46 (mean of the four transits' fitted depths; Rp/R*^2 = 1773 +60/-45) |
| Duration (h) | 5.27 ± 0.05 |
| Rp/R* | 0.0421 +0.0007/-0.0005 |
| Impact parameter | 0.23 ± 0.20 |
| Stellar density from the shape | 1.56 +0.12/-0.33 solar (TIC: 1.61) |
| Planet radius (R_Earth) | 3.75 ± 0.24 (R* = 0.817 ± 0.05 R_Sun, TIC v8.2) |
| Disposition | PC |

| Transit | Product | Mid-time (BTJD) | O-C (min) | Depth (ppm) |
|---|---|---|---|---|
| S39 | TESS-SPOC | 2384.6094 ± 0.0013 | +5 | 2203 ± 85 |
| S93 | QLP | 3838.1847 ± 0.0012 | -25 | 2050 ± 82 |
| S100 | own FFI photometry | 4092.5851 ± 0.0024 | +6 | 2075 ± 150 |
| S102 | own FFI photometry | 4128.9449 ± 0.0015 | +34 | 2045 ± 86 |

Notes: four transits; the predicted S66 transit (BTJD ~3111.40) falls on a stray-light ramp (TESS-SPOC good data
start 0.6 h before ingress, QLP flags the whole window) and is not used. Mid-times scatter about a straight line
(chi2 = 461 for 2 dof): timing systematics (S102 near an orbit start) or TTVs. Shape fitted with each transit
centred on its own mid-time. Difference-image centroids (four transits) 1.5" ± 1.3" from the target; the only
unresolved Gaia neighbour (G = 16.15 at 2.2", parallax 0.69 mas, background) would need a ~74% eclipse; diluted
eclipse fits are worse than a planet by dchi2 >= 142 (k = 0.3-1.2, 2227 points). Gaia DR3 RUWE 1.00. No secondary
eclipse. Odd/even not testable (one odd transit). TRICERATOPS 1.1 on these inputs (10 runs, no contrast curve):
FPP = 0.0026 ± 0.0037 (worst run 0.012), NFPP = 4e-24 (planet on target 85%, diluted by an unseen companion 14%,
planet around an unseen companion 0.3%); default 5x5-pixel apertures: FPP 0.0010, NFPP 6e-7. Below the numerical
thresholds, but not a validation without high-resolution imaging. The 12 October run (FPP 0.015) used the S66
window and a 1449 ppm depth and is withdrawn. Next transits (±20 min): 2027-02-19 11:31 UTC (TESS S114),
2027-05-03 03:50 UTC (full transit at night from Chile).

## TOI-2493.01 — period for a single-transit TOI
| Field | Value |
|---|---|
| TIC | 123664207 |
| Period (d) | 79.0302 ± 0.0004 |
| Transit epoch (BJD) | 2461043.5093 ± 0.0042 |
| Depth (ppm) | 1678 (fitted transit depth; Rp/R*^2 = 1467 ± 209) |
| Duration (h) | 6.10 ± 0.31 |
| Rp/R* | 0.0383 ± 0.0026 |
| Planet radius (R_Earth) | 2.7 (R* = 0.652 R_Sun) |
| Transits used | S33 2225.82, S98 4043.53, S107 4280.63 (BTJD), SPOC 2-min PDCSAP |
Notes: the S33 single transit (TOI comment "single at TBJD~2226") plus two more in sectors 98 and 107; every
shorter alias predicts transits inside sectors 6, 33, 87, 98 or 107 that are absent, and our own sector-108
photometry is flat at the P/3 and P/4 transit times (-620 ± 200 and -240 ± 180 ppm vs ~1500). The TOI table's
2.09 h duration (and depth 1704 ppm, epoch 2459225.8233 ± 0.0039) likely comes from SPOC's 15.1 d TCE fit; every
transit here lasts 6.1 h, and this ephemeris puts the S33 transit 2.2 sigma earlier (2459225.815). Upload as an
additional parameter set on the existing TOI. Centroids (2-min pixels)
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
