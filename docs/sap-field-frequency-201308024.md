# Fixed-frequency field comparison

A separately declared postselection check evaluated all 24 previously opened QLP SAP curves in both sectors at **4.194743460281676 days**. Each curve uses the same per-gap offset/slope sinusoidal model, with ordinary and Huber-2.5 fits and no frequency optimization. Original FITS input hashes and the helper source hash are bound in `runs/sap-field-frequency-201308024/declaration.json`.

TIC 201308024 has the largest robust fitted semiamplitude of all 24 in each sector. S1 gives **770 ppm**, compared with **309 ppm** for the next largest curve. S2 gives **2,505 ppm**, compared with **264 ppm** for the next largest curve. Both runners-up are TIC 355767420, on the same camera/CCD as the primary in their respective sectors. The primary's S2 amplitude is approximately 9.5 times that runner-up.

This makes a broadly shared detector-wide signal less compelling as the whole explanation. These stars are the already selected and examined exploratory sample, not a matched null population; amplitudes and ranks have no calibrated false-alarm probability. Position-dependent effects and unresolved contamination remain possible. The new comparison is not independent validation or a physical-period measurement.

All 48 comparison rows are retained in `runs/sap-field-frequency-201308024/results.json`; the helper is `scripts/sap_field_frequency_201308024.py`. No additional flux was retrieved or protected season opened.
