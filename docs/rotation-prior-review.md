# Scoped published rotation-period review

This review checks five exploratory leads against published TESS rotation work. It is not an exhaustive literature clearance. The ADS topic search was truncated (145 papers, newest 100 returned); absence from that list cannot establish novelty.

The authoritative CDS ReadMe for Claytor et al. (2024), `J/ApJ/962/47`, identifies table 5 as 7,245 detected rotation periods and table 6 as 35,664 rotation nondetections. Exact TIC queries for 201308024, 355767420, 281593180, 149392671, and 150104580 returned no detected-period rows. TIC 149392671 appears in **table 6, the nondetection sample**, with photometric range 2,142 ppm and CDPP 210 ppm. The other four have no row in either table. A nondetection does not establish constant flux, or invalidate the nearby known RR Lyrae contamination found in our analysis.

Initial strict parsing rejected VizieR metadata units `(10-6)` and `log(Sun)`. These were reviewed against the ReadMe: the former describes ppm and the latter solar logarithmic abundances. The local parser accepts only these specific additional unit warnings and generated XML-ID warnings; the shared catalogue parser was not changed. Original failed parsing records and raw bytes remain preserved in `data/rotation-prior-review`.

Boyle, Bouma, and Mann (2026), *The TESS All-Sky Rotation Survey*, reports 1,046,317 stellar variability periods, estimating that 93% are rotation periods. The candidate CDS path `J/ApJS/284/75` is **unavailable**, not a zero-match catalogue: the ReadMe returned HTTP 404 and the VizieR response explicitly says the catalogue was not found. Its aggregate catalogue and code are advertised on Zenodo (current DOI 10.5281/zenodo.19917941; arXiv v1 references an earlier DOI). Zenodo access failed with a network-proxy 403. Aggregate-catalogue membership therefore remains unchecked.

The official MAST TARS HLSP page and ReadMe are accessible. They describe published sector-periodogram statistics in FITS headers. We queried exact TIC metadata restricted to sectors 1 and 2, downloaded those products only, and read header statistics without inspecting flux arrays. These are alternate reductions of the same TESS exposures, rather than observations by another telescope. Highest periodogram peaks are not automatically accepted stellar rotation periods.

| TIC | Published TARS sector-1/2 metadata | Implication |
|---|---|---|
| 281593180 | S1 PER3 2.3599487 d; S2 PER3 2.3708787 d. Dominant peaks 18.1322385 and 19.0262182 d. | A roughly 2.37-day periodicity already appears in publicly released TARS measurements; do not describe its discovery as novel. |
| 150104580 | S2 PER4 1.1989171 d; S1 PER5 1.142793 d. Dominant peaks 10.2997045 and 7.3480909 d. | The approximately 1.2-day signal has a published sector-2 counterpart; this does not identify its physical source. |
| 149392671 | Highest five header peaks differ between sectors and exclude 0.477 d. | This limited peak list cannot exclude known RR Lyrae contamination or establish novelty. |
| 201308024 | No exact-TIC indexed TARS S1/S2 observation returned. | Aggregate membership and alternate target-ID matching remain unresolved. |
| 355767420 | No exact-TIC indexed TARS S1/S2 observation returned. | Aggregate membership and alternate target-ID matching remain unresolved. |

The TARS header for TIC 149392671 stores GAIADR3 as a floating-point value, losing integer precision. It must not replace the independently verified exact Gaia identifier. TIC 281593180 uses Gaia DR3 4906669771092140928 in TARS; its TIC metadata lists the distinct DR2 identifier 4906669771091724928. A change of Gaia release is not a positional mismatch by itself.

Raw ReadMes, HTML, MAST request/response pages, VOTables, S1/S2 FITS files, and URL/UTC/byte-count/SHA-256 provenance are retained in `data/rotation-prior-review`. `claytor-results.json`, `tars-observation-summary.json`, and `tars-header-periods.json` contain the scoped results. No held-out sector flux was opened.
