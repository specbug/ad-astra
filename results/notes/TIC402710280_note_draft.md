# DRAFT — A sub-Neptune candidate on a 36.34-day orbit around the K dwarf TIC 402710280

*Draft research note (RNAAS-style, under 1000 words). No author list; not submitted. Numbers from this
repository's results (12 October 2026).*

## Summary
We report TIC 402710280.01, a transiting sub-Neptune candidate on a 36.3397-day orbit around a
nearby K dwarf (TIC 402710280 = TYC 9041-03677-1; T = 9.24, G = 9.78, 66.8 pc). Five transits were found in TESS
data from 2021 to 2026, two of them only in our own photometry of the sector 100 and 102 full-frame images, for
which no public light curves exist yet. The star is not a TOI or community TOI, and we found no prior mention of
it in the literature.

## Data and detection
A blind single-transit search of public TESS light curves found two transits, in sector 39 (TESS-SPOC, BTJD
2384.607) and sector 93 (QLP, BTJD 3838.214), 1453.6 d apart. Light curves we made from the sector 99–108
full-frame image cubes (astrocut cutouts, background from faint pixels, aperture chosen for minimum 1-hour
scatter) show two more transits, in sectors 100 (BTJD 4092.595) and 102 (BTJD 4128.936). The only period
consistent with all four and with the rest of the data, P = 36.340 d, predicts a fifth transit in sector 66,
which is present (TESS-SPOC, BTJD 3111.411). Every shorter alias predicts transits inside the data that are
absent, and the only other transit predicted during TESS observations falls in a data gap (sector 103).

## Transit model
A joint fit (batman + emcee, exposure-integrated for each cadence) of the sectors 39, 66, 93 and 100 transits
gives P = 36.33972 ± 0.00007 d, T0 = BTJD 3474.8014 ± 0.0016,
Rp/R★ = 0.0381 ± 0.0029, b = 0.46 ± 0.32, T14 = 5.19 ± 0.18 h, and a stellar density of
1.27 ± 0.81 ρ☉, consistent with the TIC value for the K dwarf (1.61 ρ☉; fixing it changes χ² by 0.2).
The sector 102 transit lies on a ramp 1.4 d after the start of an orbit and is left out of the shape fit.
With R★ = 0.82 R☉ the planet radius is 3.4 R⊕; the orbit (a = 0.206 AU) receives
10 times Earth's insolation (T_eq ≈ 498 K).

## Vetting
- Pixels: difference-image centroids of all five transits, measured relative to the star's own image on each
  frame (Gaia DR3 positions propagated to the epoch), lie 1.4″ ± 1.2″
  from the target; every Gaia star bright enough to cause the dip lies ≥ 22σ away except a G = 16.2 background
  star 2.2″ away (parallax 0.69 mas vs. 14.96 mas).
- That star would need a ~73%-deep eclipse; fitting the transits as a diluted eclipse with a true radius ratio
  k ≥ 0.5 is worse than the planet model by Δχ² ≥ 65 (Δχ² = 27 even at k = 0.3).
- No secondary eclipse at phase 0.5 (≤ 150 ppm); odd and even transit depths agree (2.3σ).
- Gaia DR3: RUWE 1.00, no non-single-star solution, a single image in all scans.
- Our pixel test, run blind on 40 confirmed TESS planets, never placed a planet's transits on a neighbour.

## Predictions
Next transits (UTC): 2026-11-02 10:54, 2026-12-08 19:03, 2027-01-14 03:12, 2027-02-19 11:22. TESS is expected to observe the star
again in sectors 113–114. The expected radial-velocity semi-amplitude is ≈ 2.4 m/s for a
11 M⊕ planet, within reach of HARPS or ESPRESSO for a V ≈ 10 star; high-resolution imaging would
exclude companions closer than Gaia resolves.

## Data availability
Code and results: https://github.com/specbug/ad-astra (pipeline/ffi_new.py, fit_new.py, centroid_vet.py;
results/candidates/new_candidates.json).
