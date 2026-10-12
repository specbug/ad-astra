# Native 120-second SPOC check around TIC 201308024

A separately declared S1/S2 native-SPOC follow-up downloaded four exact-identity public products for the bright neighbor TIC 201308023: two light curves and two target-pixel files. All original FITS, URI/HTTP/SHA256 provenance, supplied checksum verification, timing metadata and failures are retained in `data/native-spoc-201308023` and `runs/native-spoc-201308023*`. No later-sector flux was accessed.

The primary TIC 201308024 falls inside both actual 11×11 target-pixel stamps: S1 zero-based coordinates (8.752, 1.822), S2 (9.219, 4.064). Full 2×2 and 3×3 primary apertures fit within the stamps and do not overlap the bright-neighbor 3×3 aperture. Native time ingestion separately checks TICID, SECTOR, TIMESYS=TDB, TIMEUNIT=d, BJDREFI/BJDREFF and the TIME-column reference. The original strict reader rejected the publisher's column label `BJD - 2457000, days`; its source, declaration and all four failures were preserved. A second declared run accepts that exact native form only when its embedded reference agrees with BJDREFI+BJDREFF. Six timing/reference/identity tests pass.

The analysis fixes the two previously selected periods, 2.03270208338902 and 4.194743460281676 days; it does not optimize periods. Quality-zero, finite samples with positive native errors are used without clipping. Weighted sin/cos fits include separate offsets and slopes across gaps over 0.5 days. Pixel errors use quadrature sums, without a covariance model. Pixel flux already has the SPOC pipeline background removed; a remaining border-median covariate excludes primary/neighbor 3×3 masks. Optional models include POS_CORR1/2 and the summed neighbor aperture. Native LC control models use SAP_BKG, MOM_CENTR1/2 and POS_CORR1/2. There are 18,277 usable samples in S1 and 18,299 in S2; all 48 declared fits complete.

| Sector and fixed period | Primary 2×2 base / joint controls (ppm) | Primary 3×3 base / joint controls (ppm) | Neighbor native SAP base / controls (ppm) | Neighbor PDCSAP base / controls (ppm) |
|---|---:|---:|---:|---:|
| S1, 2.032702 d | 1486 / 1755 | 1756 / 1763 | 38.6 / 35.0 | 9.8 / 8.7 |
| S1, 4.194743 d | 1312 / 899 | 970 / 857 | 101.5 / 104.8 | 5.3 / 4.7 |
| S2, 2.032702 d | 679 / 403 | 498 / 389 | 19.4 / 15.6 | 15.0 / 14.8 |
| S2, 4.194743 d | 2412 / 2552 | 2392 / 2445 | 5.3 / 13.1 | 12.0 / 12.8 |

Values are sinusoidal semiamplitudes, normalized separately by each aperture/LC median flux; unequal aperture losses and contaminating flux prevent direct conversion of these fractional amplitudes into intrinsic stellar variability. `joint controls` means neighbor aperture, border median and position corrections. Formal residual-scaled white-noise amplitude errors for the S2 4.194743-day joint fits are 34.0 ppm (2×2) and 28.9 ppm (3×3); correlated errors, pointing-dependent losses and postselection are omitted, so these are not calibrated discovery significances.

The approximately 4.19-day S2 signal survives native 120-second primary-aperture extraction and nuisance controls, while the neighbor's standard LC has much smaller modulation at that period. This is useful processing/cadence and spatial contrast evidence. It does not establish intrinsic source attribution or novelty. The same detector and underlying observations are used, no PRF deblend or aperture-loss correction has been applied, and other nearby sources remain possible. Phase plots show substantial outliers and scatter; robust fitting and day-block sensitivity checks remain appropriate before promoting the result. S1 dominance near 2.03 days versus S2 near 4.19 days remains descriptive, rather than a stable physical-period determination.

Machine-readable results and four fixed-period figures are in `runs/native-spoc-201308023-analysis`. Its declaration binds the retrieval and both source hashes. No holdout or new period search is included. A next useful comparison is the Gaia-informed TGLC extraction for the exact primary source, with the same fixed periods and explicit product quality/timing checks.

## Declared robust and day-block sensitivity addendum

A subsequent separately declared postselection check fixes only the S2 4.194743-day period and retains the same joint primary-aperture nuisance models. Its source, declaration, original-input hashes, full fits, individual deletions and consistency verification are retained in `runs/native-spoc-201308023-sensitivity`. Huber IRLS uses a fixed 2.5 threshold, standardizes residuals by each sample's own photon error, and estimates residual scale by MAD. Each available integer `floor(BTJD)` day is omitted in turn, with all nuisance coefficients and gap groups refitted. Native neighbor SAP/PDCSAP use their instrument-control model. There is no new search or calibrated false-alarm calculation.

| S2 series | Full WLS / Huber amplitude (ppm) | Huber leave-one-day amplitude range (ppm) | Huber phase change relative full fit (cycles) |
|---|---:|---:|---:|
| Primary 2×2 | 2551.5 / 2553.5 | 2499.3–2657.9 | −0.00431 to +0.00873 |
| Primary 3×3 | 2445.2 / 2447.1 | 2373.4–2544.7 | −0.00489 to +0.00958 |
| Neighbor SAP | 13.12 / 13.37 | 3.35–33.26 | −0.343 to +0.437 |
| Neighbor PDCSAP | 12.81 / 12.40 | 3.73–27.76 | −0.0843 to +0.1290 |

The primary-aperture signal survives downweighting outliers and removal of any one observed day, with small changes in fitted phase. All Huber fits converge. This reduces the concern that a few outlying samples or one observed day account for the signal. It does not guarantee an intrinsic primary source or establish a physical period; correlated pointing/aperture effects or another nearby source may remain. The neighbor's small fitted amplitudes have unstable phase under day deletion, so their phase is not interpreted. All native input hashes remain unchanged, and the new ordinary-WLS full fits reproduce the original four matching amplitudes to within 0.000001 ppm.
