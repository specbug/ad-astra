# Second-cohort exploratory vetting

Independent inspection of the four period-match priorities and all eight contact sheets finds no convincing transiting-planet candidate. All 24 targets were viewed in original SAP, upstream-corrected KSPSAP, and fixed-S1 transit folds. This is postselection inspection of two already searched sectors; it is not untouched validation. No holdout flux was opened. The machine-readable audit is `runs/qlp-overlapfield-vetting/audit.json`, with dossier and image hashes, full local-event summaries, inverse controls, and decisions.

| TIC | Matched periods S1 / S2 (days) | S1 fixed-window median depth / S2 (ppm) | Positive measurable windows S1 / S2 | Decision |
|---|---:|---:|---:|---|
| 425933704 | 0.629802 / 0.629377 | 216 / −1.9 | 24/29 / 15/30 | Reject promotion: phase disagreement |
| 231910457 | 0.656955 / 1.320185 | 446 / 113 | 21/28 / 14/24 | Reject promotion: factor-two match lacks persuasive repeating shape |
| 234490600 | 0.994329 / 2.984700 | 325 / −0.5 | 17/24 / 12/25 | Reject promotion: fixed S1 signal absent in S2 |
| 257655479 | 2.392358 / 0.795502 | 183 / −83 | 10/11 / 5/11 | Reject promotion: fixed S1 S2 depth reverses sign |

Depths in the table use the existing fixed running-median flattening. Positive windows are signs of local sample medians, not statistically significant transit events.

TIC 425933704 is the only direct period match. Its independently fitted S2 centers lie a median 0.17135 days from the nearest S1 prediction, compared with a 0.052-day geometric overlap tolerance; zero centers overlap. Its short-duration box has a 6.35% duty fraction. Applying the actual S1 ephemeris to S2 gives a median depth of −1.9 ppm, and unflattened KSPSAP/SAP give −21.6/−32.5 ppm. S1 local dips persist in SAP (median 218 ppm) and unflattened KSPSAP (203 ppm), and no single event dominates the descriptive positive-depth sum (largest share 8.5%). Thus the S1 feature is not explained solely by one outlier or this running median, but it does not repeat at the same phase in S2. A similar numerical period is insufficient to establish a planet. The plot shows a separate positive S1 outlier and broad SAP baseline changes; these are retained in inspection rather than clipped.

TIC 231910457 has a factor-two match, not a direct period match. Its fixed S1 box has 6.09% duty fraction. S2 medians are 113 ppm after flattening, 116 ppm unflattened KSPSAP, and 107 ppm SAP, far below S1's 446/434/425 ppm. Local depths span roughly −1035 to +1713 ppm in S1 and −726 to +1115 ppm in S2. The folded S2 central dip is weak amid other phase structure and scatter. An eclipsing binary can produce harmonic matches, but these data do not justify that classification either; neither significance nor independent timing coherence is established.

TIC 234490600 has a factor-three match. Its 6.44% duty box gives mixed-sign S2 depths; SAP and unflattened KSPSAP medians of 42/50 ppm remain much smaller than S1's 322/316 ppm. The plotted S2 time series contains low-frequency structure but no convincing narrow dip at the fixed S1 phase. This is not a persistent ephemeris.

TIC 257655479's matched S1 rank-three peak was searched after 34.1% of samples were already masked; its duty fraction is 8.36%. The S1 feature survives the alternate SAP view (median 267 ppm), but S2 medians are negative in all three views: −83 ppm flattened, −103 ppm unflattened KSPSAP, and −115 ppm SAP. Largest positive-event depth share is 22% in S1 and 55% in S2; the latter is a descriptive sum fraction, not integrated signal or statistical leverage. Its S2 fold lacks a corresponding dip.

None of these four targets has a qualifying direct or harmonic match in the sign-inverted controls. That does not rehabilitate the original matches: inversion is an asymmetric-noise diagnostic, not a calibrated false-positive probability. None has strongest Lomb–Scargle period agreement between sectors.

The all-target inspection did reveal SAP behavior worth a separate processing audit. TIC 201308024 shows repeated swings of several thousand ppm on a few-day scale in both sectors; TIC 355767420 shows broad, roughly two-week troughs approaching 6000 ppm; TIC 281593180 shows substantial broad modulation. KSPSAP strongly suppresses all three. These are followup priorities for SAP period fitting, comparison-star/common-mode analysis, background and centroid comparison, and independent aperture extraction before interpreting them as rotation or other stellar variability. Other targets have similar edge/gap trends, which makes spacecraft or reduction effects a serious competing explanation. No unambiguous missed transit or independently coherent variable-star discovery was seen in the contact sheets.

The prior catalog review's zero nearby known-variable matches and Gaia 90-arcsecond prefilter do not prove novelty or exclude uncataloged contamination. A negative transit vetting decision here also does not prove the stars are constant or lack planets below the pilot sensitivity.
