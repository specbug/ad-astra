# TIC 201308024: alternative FFI extraction and localization

This exploratory check finds target-centered approximately 4.19-day modulation in sector 2, while sector-1 source attribution remains uncertain. It supports further investigation of stellar variability; it does not establish a new variable star, a rotation period, or a planet. TESSCut and QLP use the same TESS full-frame images, so this is an alternative reduction rather than an independent instrument confirmation.

The initial method was declared in `runs/qlp-sap-pixels-201308024/declaration.json` before downloading any pixels: only sectors 1 and 2, a 7×7 cutout, quality zero, no clipping, fixed periods 2.032702083 and 4.194743460 days in both sectors, three harmonics, per-gap offset and slope, and baseline-only versus outer-border median nuisance models. Target-centered nominal 3×3, 5×5 and 7×7 square extractions compare raw, border-subtracted and border-nuisance models. Gaia positions include proper motion from the catalog epoch to each sector. Source files, provenance, original ZIPs, extracted FITS and checksums are preserved under `data/qlp-sap-pixels-201308024`; derived JSON and four figures are in the corresponding `runs` directory.

Both products have the requested RA/Dec in RA_OBJ/DEC_OBJ, explicit TDB time, BJD reference 2457000, day units and flux/error units e-/s. All three HDUs in each FITS supplied CHECKSUM and DATASUM cards and passed the original-byte integrity validator. There are 1,267 and 1,196 quality-zero samples in sectors 1 and 2, covering BTJD 1325.3262–1353.1598 and 1354.1390–1381.4928. Only these two sectors were retrieved. The 90-arcsecond Gaia proper-motion cone contains 12 rows, matching the independently counted original selection cone.

After inspecting the initial maps, a separate sensitivity declaration was written before its calculation: same fixed periods and baselines, Huber weighting at 3 residual MAD sigma for 12 iterations, and a 180-arcsecond Gaia cone with an independent count. This is explicitly post-inspection sensitivity, not a pre-specified validation or a new period search. It uses the already downloaded pixels. `runs/qlp-sap-pixels-201308024-sensitivity` contains the declaration/results, and `data/qlp-sap-pixels-201308024-sensitivity` preserves the complete 34-row cone and count.

The expanded cone reveals Gaia DR3 4906017863775668992, G=9.70135, 92.29 arcseconds from the requested target position, with proper-motion-adjusted pixels (−0.25,6.06) in sector 1 and (−1.49,3.70) in sector 2. Its wings dominate the left cutout edge. Another G=13.9463 source, Gaia DR3 4906017863775668864, lies nearby at 93.88 arcseconds from the target. Both fall just outside the original 90-arcsecond catalog filter. Their Gaia phot_variable_flag values are NOT_AVAILABLE, which does not mean constant. Neither source is identified here as a physically associated companion. The target Gaia DR3 4906011953900670976 has G=12.2878 and pixels (3.72,3.82)/(3.22,3.06).

## Quantitative sensitivity

The table gives fundamental semiamplitude in ppm from robust target-centered square extractions. Three harmonics were simultaneously fit. Values refer to each extraction's median total flux, so dilution and background change the denominator.

| Sector | Fixed period (d) | Aperture | Raw | Border-subtracted | Border nuisance |
|---|---:|---|---:|---:|---:|
| 1 | 2.032702083 | 3×3 | 916 | 1,960 | 1,241 |
| 1 | 2.032702083 | 5×5 | 2,797 | 1,820 | 952 |
| 1 | 4.194743460 | 3×3 | 2,076 | 1,122 | 917 |
| 1 | 4.194743460 | 5×5 | 1,801 | 685 | 415 |
| 2 | 2.032702083 | 3×3 | 641 | 509 | 361 |
| 2 | 2.032702083 | 5×5 | 1,031 | 684 | 361 |
| 2 | 4.194743460 | 3×3 | 2,371 | 2,469 | 1,821 |
| 2 | 4.194743460 | 5×5 | 2,309 | 2,219 | 1,273 |

In sector 2 at 4.194743460 days, the pixel with greatest fundamental amplitude is (3,3), the target-nearest pixel. Its amplitude is 1.7396 e-/s in the robust baseline model and 1.5409 e-/s with the border nuisance. Robust weighting does not weaken this conclusion: the 3×3 raw extraction downweights zero samples, and the nuisance extraction downweights only 11 of 1,196. The modulation footprint resembles the target's central flux distribution, unlike the 2.0327-day footprint concentrated near the left bright-source edge. This is useful descriptive evidence for target-localized sector-2 variability.

In sector 1, large outliers near BTJD 1348 materially change the unrobust fits. After robust weighting both fixed-period pixel maxima are at (0,6), close to the bright G=9.70 neighbor. The target aperture still has approximately 2-day modulation, but its amplitude changes with aperture/background treatment, preventing secure attribution. Rotation spots can change harmonic dominance and phase, so this sector disagreement alone does not rule out target rotation. It also does not prove that the two selected periods reflect the same astrophysical process.

The nominal 7×7 target-centered mask is truncated in sector 1 because the target rounds to pixel (4,4) inside a 7×7 cutout: it contains only 36 pixels. It contains all 49 pixels in sector 2. JSON records actual pixel counts. These largest masks include strong neighbor wings and are not clean target photometry; they are primarily sensitivity diagnostics. A 9×9 or larger cutout with more complete background treatment would be a useful next reduction.

The border median shows strong orbital stray-light trends. Positive-flux centroid diagnostics are also stored, but neither their formal errors nor pixel amplitude maxima yield a source-localization confidence interval. Border subtraction can remove stellar wings and introduce their modulation into an aperture; nuisance regression can absorb real stellar variability correlated with background. Pixel models omit temporal/interpixel covariance, PRF blending and a calibrated position-dependent aperture-loss model. Nearby unresolved sources remain possible. No false-alarm probability, novelty clearance, independent-instrument confirmation, or holdout test is claimed.

A productive next step is a larger-cutout reduction with explicit neighbor treatment, followed by literature and independent-observatory checks of the target and bright neighbor. The two original selected periods and all sensitivity limitations must remain visible in that work.

## Declared 11×11 sensitivity

The additional `runs/qlp-sap-pixels-201308024-11x11/declaration.json` was written before downloading new 11×11 cutouts. Original 7×7 source/results were copied and hashed in its `source-snapshot`. Separate data and source paths preserve those earlier results. The new method retains both fixed periods, Huber3/12 iterations, three harmonics, per-gap offsets/slopes, and quality-zero samples. It compares two median backgrounds: all outer-border pixels, and the border excluding pixels within 2.5 pixels of every Gaia G<17 source in the preserved 180″ cone. Twenty-two and seventeen clean border pixels remain in S1/S2. Target and neighbor 3×3/5×5 masks exclude pixels within 1.5 pixels of the other star. All actual masks and counts are stored. The neighbor's S2 5×5 mask is truncated to 20 pixels at the cutout edge; its 3×3 is complete. The target masks contain 9/25 pixels in both sectors.

The larger cutout shows that the 7×7 sector-2 global pixel maximum depended on the frame boundary: the bright neighbor now also has visible approximately 4.2-day pixel modulation. Thus the earlier target-nearest maximum must not be treated as exclusive source localization. The target nevertheless retains a distinct approximately 4.2-day modulation footprint, whereas the approximately 2-day footprint is much stronger at the neighbor.

With the Gaia-masked border nuisance model, the complete 3×3 aperture fundamental amplitudes are:

| Sector | Period (d) | Target (e-/s) | Target (ppm) | Neighbor (e-/s) | Neighbor (ppm) |
|---|---:|---:|---:|---:|---:|
| 1 | 2.032702083 | 4.19 | 1,279 | 4.69 | 165 |
| 1 | 4.194743460 | 3.08 | 940 | 7.33 | 257 |
| 2 | 2.032702083 | 1.34 | 427 | 7.25 | 277 |
| 2 | 4.194743460 | 5.86 | 1,871 | 1.66 | 64 |

Sector-2 target 4.19-day amplitude also persists in the 5×5 extraction at 1,357 ppm with the masked-background nuisance, or 2,155 ppm after subtracting that background. Changing from all-border to bright-masked-border changes the 3×3 target nuisance estimate only from 1,881 to 1,871 ppm. Source-normalized ppm values are diluted by each aperture's mean flux; the target/neighbor comparison above therefore also gives absolute modulation in e-/s.

Joint regressions include the fixed harmonics alongside background and background-subtracted neighbor flux. An additional declared mode includes whole-cutout positive-flux centroid x/y as nuisances. In S2, target 3×3 4.19-day amplitude remains 1,903 ppm after neighbor/background regression and 1,754 ppm with centroids added; target 5×5 estimates are 1,366 and 1,033 ppm. The corresponding approximately 2-day amplitudes fall to 310/200 ppm (3×3) and 251/158 ppm (5×5). All joint designs are full rank. Error-weighted condition numbers range approximately 9.7–45.3; the 3×3 S2 models have 13.4 without centroids and 32.9 with centroids. Target/neighbor raw-versus-background-subtracted Pearson correlations are 0.295 (3×3 S2) and 0.265 (5×5 S2), and are descriptive rather than causal attribution.

At S2 approximately 4.2 days, the complete 3×3 target and neighbor fundamental vectors differ in phase by −0.445 cycles and the target absolute amplitude is 3.53 times the neighbor's. At approximately 2 days, the target amplitude is only 0.184 times the neighbor's, with phase difference −0.106 cycles. These comparisons follow from the same fitted fixed-period coefficient vectors; they are not separate experiments. Opposite-phase modulation can also arise from aperture loss, motion or shared processing, so this strengthens the motivation for a PRF/position model rather than eliminating that alternative.

The 11×11 sensitivity therefore keeps target-centered S2 approximately 4.2-day variability as a worthwhile lead, while weakening any claim based on a unique global pixel maximum. The data do not establish a single coherent S1/S2 rotation period, clean source attribution, or novelty. The bright neighbor has an unresolved nearby Gaia G=13.95 source, and no PRF deblending or calibrated position-dependent aperture-loss model has been performed. Both cutouts again passed original-byte CHECKSUM/DATASUM validation in all three HDUs, with exact requested coordinates, explicit timing/flux units and literal S1/S2 guards. No download or analysis failures occurred.
