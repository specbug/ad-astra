# Independent ASAS-SN ground photometry for TIC 201308024

Public **optical measurements**, rather than a catalog nonmatch, are now retrieved for the best unresolved second-pilot lead and its bright neighbor. The legacy [ASAS-SN Photometry Database](https://asas-sn.osu.edu/photometry) provides existing light curves without an account or computation request. Its landing page says the database was updated 07/02/2021; these are historical V-band measurements, not newly processed/current Sky Patrol photometry.

## Identity and retrieval

| Object | Public source | ASAS-SN coordinates (degrees) | Separation from supplied Gaia/SIMBAD coordinate | CSV measurements |
|---|---|---|---|---|
| TIC 201308024 | [AP37233362](https://asas-sn.osu.edu/photometry/ed536415-2c44-5c67-adc2-01042a699136) | 3.80301, −59.70528 | 0.157 arcsec | 258 |
| TIC 201308023 / CD-60 23 | [AP37233381](https://asas-sn.osu.edu/photometry/2b960b56-6180-5f5a-9e9f-f7e6596e1653) | 3.84931, −59.69516 | 1.603 arcsec | 258 |
| Unrelated field control | AP37233363 | 3.82641, −59.69489 | Not assigned to a TIC | 260 |

Target identity was checked with a 0.1-arcmin cone. The first 2-arcmin target cone returned the target and AP37233363; the latter was provisionally labeled “bright neighbor” before coordinate validation. This was corrected: it is approximately 40 arcsec from the actual bright neighbor, and has V≈15.63 rather than V≈9.95. The exact bright-neighbor cone returned AP37233381. Retrieved field-control files were renamed and retained; `availability-summary.json` records the correction.

HTML source pages, full CSV exports, the first JSON response page, query responses and UTC/URL/SHA256 provenance are in `data/ground-validation-201308024/`. JSON exports are paginated: those saved responses are **not** the complete curves. For target and bright neighbor, the full CSV exactly matches all 258 observations embedded in each original HTML page. Their page-level summary says 252 and 251 epochs respectively; these stale/inconsistent summaries were not substituted for the actual exported records. The field-control CSV has two additional rows at HJD 2457067.51515 with magnitude-error 99.99 (one magnitude 0); the HTML embedded array omits them. The discrepancy is preserved, and that faint source is unsuitable as a precision-matched control.

## Time, errors and image limitations

Both main curves cover HJD **2456791.89379–2458385.65923**, a 1,593.765-day baseline, with 73 `bf` and 185 `bh` camera observations. Source-page observation tables label the filter V; the CSV itself lacks a filter column. The website explicitly states that HJD is computed for **field centers**, with approximately 200-second timing precision, and recommends a separate reprocessing pipeline for more accurate timing. The saved page does not declare the time scale or exposure midpoint convention. We therefore do not treat these times as BJD_TDB or extrapolate a TESS phase reference over years.

The CSV has `hjd,camera,mag,mag_err,flux,flux_err` but no quality/saturation flag. Target and bright-neighbor median magnitude errors are both 0.020 mag. This appears to be a floor and is far above an individual 1,754-ppm sinusoid (~0.00190 mag). `Blend=false` is a database annotation, not proof of isolation. Typical ASAS-SN spatial resolution is much coarser than Gaia; the primary/bright-neighbor separation is 92.29 arcsec, but unresolved close companions and PSF wings remain possible. We have not reconstructed the survey PSF/aperture or quantified contamination. The bright neighbor's V≈9.95 and large scatter warrant saturation/systematics review; no per-observation saturation status is supplied and saturation is **not established** here.

## Declared fixed-period check

`scripts/ground_fixed_periods_201308024.py` writes its source/input hashes and experiment declaration **before fitting**. It uses the previously identified TESS periods exactly: **2.03270208338902 and 4.194743460281676 days**, with weighted sine/cosine magnitude fits and additive camera/season offsets. Seasons split at gaps longer than 90 days. All finite HJD/magnitude rows with their own positive finite magnitude error are included; there is no clipping, period optimization, TESS phase extrapolation or calibrated significance estimate. Results and plots are in `runs/ground-validation-201308024-v2/`.

| Curve | Fixed period (d) | Magnitude semiamplitude | Formal white-noise amplitude error (mag) | Residual-scaled amplitude error (mag) |
|---|---|---|---|---|
| Target | 2.032702 | 0.001089 | 0.001831 | 0.001345 |
| Target | 4.194743 | 0.002098 (~1,932 ppm) | 0.001741 | 0.001274 |
| Bright neighbor | 2.032702 | 0.003423 | 0.001783 | 0.007767 |
| Bright neighbor | 4.194743 | 0.011644 | 0.001780 | 0.007718 |

The target's offset-only empirical RMS is **0.01448 mag**, and its quoted-error reduced chi-square is 0.537. Its 4.19-day fitted amplitude is only about 1.2 times the formal white-noise amplitude error. The bright neighbor has offset-only RMS **0.08577 mag** and reduced chi-square **18.83**, so its tiny nominal quoted-error amplitude error is misleading. Multiplying by residual scatter does not model red noise, but makes the inadequacy clear. These fits neither confirm the TESS signal nor establish that the primary or neighbor is constant.

The unweighted sampling-window modulus at one cycle/day is **0.771** and at one cycle/year **0.434**, warning of strong ground-survey daily/seasonal structure. No alias search or false-alarm probability was computed.

Four predeclared injected 1,754-ppm fractional-flux sinusoids at the exact 4.19-day period and phases 0, 0.25, 0.5 and 0.75 cycles were transformed to magnitude exactly, using the observed ground times/errors and camera/season design. Paired fitted coefficient increments recover ~**0.001904 mag** in each case, verifying response to an added sinusoid. The total fitted amplitudes after injection range from **0.00153 to 0.00370 mag**, of order the quoted-error uncertainty. This is an algebraic injection-response check, **not a detection-completeness or significance calibration**; upstream ASAS-SN reduction was not replayed. V and TESS passbands and time-variable stellar amplitudes also limit direct comparison.

The scientific conclusion is **independent data obtained, current ground check inconclusive**. The curves are useful for future controlled modeling, but do not provide discovery or novelty clearance.

## Failures and API availability

The first fitting attempt completed its calculations/figures but failed while serializing a NumPy integer rank-derived degree of freedom. Its original source, declaration and failure record remain in `runs/ground-validation-201308024/`; the corrected integer serialization and separate completed output are in `...-v2/`. The scientific model was unchanged.

The official modern `pyasassn` client source and README were retrieved from the ASAS-SN GitHub repository. Its documented API host `asassn-lb01.ifa.hawaii.edu:9006` returned proxy **403 “Domain forbidden”** on a schema GET. This is a network prerequisite, not evidence that authentication is required or target data are absent. No bypass or new photometry-computation submission was attempted. The legacy database independently supplied the measurements above. OGLE root/OCVS availability probes are retained but no OGLE light curve was retrieved in this task.
