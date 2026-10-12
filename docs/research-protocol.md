# Reconstructed QLP pilot protocol

The aim is to find defensible leads in public TESS photometry while retaining failed controls, prior discoveries and ambiguous explanations. This is an exploratory pilot, not a calibrated occurrence-rate or completeness study. ADS is excluded at the user's request.

## Selection before science flux

The fixed sky field is centered at RA 90°, Dec −66.5°, radius 5°. TIC cuts are Tmag 8–12, radius 0.3–1.5 solar radii, contamination ratio below 0.1 and absolute Galactic latitude above 15°. Seed 20261009 fixes the target order. The complete retrieved TIC population contains 1,155 rows. Frozen metadata exclusions include NASA TOIs and confirmed planets, ExoFOP CTOIs, SIMBAD variable/binary types and specified previously discussed targets. Publication-based amendments are retained separately before opening science data. An absent catalog match is not a novelty claim.

Discovery requires QLP sectors 1 and 2, with at least two metadata-confirmed later products in a separated season. Product selection uses exact TIC, sector, URI, observation identifier and publication provenance. The science runner rejects all sectors other than 1/2, independent of the named holdout list. This is necessary because the previous runtime's complete protected-target manifest was lost. No later science flux may be accessed without first resolving that conflict.

## Declared screening

Use quality-zero rows with finite times and positive finite flux and reported errors. Stable-sort times and keep the first qualifying row at each duplicate time. Do not sigma-clip. Normalize by the median. For transit screening, divide by a 0.5-day running median separately across gaps exceeding 0.5 day. Search periods 0.5–min(13, baseline/2) days, durations 0.04/0.065/0.1/0.15/0.2 days, and a dense uniform-frequency grid whose worst nearest-grid phase drift is at most a quarter of the shortest duration over the baseline. Stop rather than silently coarsen a grid above 200,000 points.

Inspect up to three box-least-squares peaks, masking 1.5 durations around previously selected peaks. Require at least three sampled windows and three positive local median depths for an inspection priority. These counts do **not** count statistically significant events. Compare all up-to-nine cross-sector peak pairs, recording direct and factor-two/three coincidences within 1%. Conditional depth errors omit correlated noise and search selection.

For stellar variability, use the strongest weighted Lomb–Scargle sinusoid in the median-normalized KSPSAP series without the additional median division. QLP has already processed this flux, so it is not raw detector photometry. Search periods 0.05–min(20, baseline/2) days, with five frequency samples per peak. Record amplitude and cross-sector period agreement; assign no survey false-alarm probability or stellar classification.

For every discovery curve, reflect its cleaned median-normalized flux as `2 - flux` before additional detrending and apply the identical search. Preserve errors and input rows. Fail that control explicitly if reflection creates nonpositive samples. These sign controls expose symmetric noise and variability; they are neither independent null samples nor an estimator of the false discovery rate.

## Controls and subsequent interpretation

The unchanged screening code recovered 18/18 injections in both period and phase: periods 1/3/7 days, depths 1,000/3,000/10,000 ppm, in π Mensae sectors 1 and 27. Injection occurs into an already processed QLP product before our additional median filter, not before the upstream QLP extraction. The real π Mensae c signal is recovered in sector 1 but missed in sector 27 in all three selected peaks. A strong localized artifact dominates that sector. Retain this failure; do not tune the pipeline to make it pass.

After screening, retain all target results. Inspect SAP versus KSPSAP, detrending effects, individual events, odd/even behavior, cross-sector phase prediction and sign controls. Follow up promising stellar behavior with independent public metadata, nearby-source information and published period catalogs. Post-screen diagnostics are exploratory and must be labeled accordingly. Any later validation hypothesis needs its own frozen prediction before untouched data are opened, with the old-holdout protection resolved first.

Scripts bind configurations and source hashes before downloading science files. They reject overwrites, source changes, identity mismatches and cache tampering. Keep every failed query or target; no hidden retry, sample replacement or discarded result is permitted. Source review approval means readiness for exploration, not scientific validation.

## Preservation

Retain source, exact dependency inventory, public metadata responses, original FITS bytes, provenance, declarations and reports. Git ignores `data/` and `runs/`; a code commit alone is insufficient. A local archive is a recovery aid, not a guarantee across destruction of the entire environment. Saved environment configuration likewise does not publish or snapshot files until the platform's publication step succeeds.
