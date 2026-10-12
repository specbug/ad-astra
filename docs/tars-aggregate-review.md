# TARS aggregate catalogue review

The official Zenodo record for Boyle, Bouma and Mann (2026), DOI 10.5281/zenodo.19917941, became accessible during continued research. Its raw JSON and URL/UTC/SHA-256 provenance are preserved in `data/tars-aggregate-review`. The record describes table 2 as the default catalogue of 1,046,317 adopted rotation periods, one row per star. Table 4 contains 39 million sector-level measurements for 7.5 million targets and is a different membership test.

This check downloads only the published default table 2 CSV ZIP and the official `HOW_TO_USE.md`. Publisher MD5 checksums and file sizes are verified before parsing. No new light curves or held-out flux are retrieved or opened. Table 2 can incorporate published measurements from sectors 1–96; consulting a published catalogue is a prior-art check, not an independent validation observation.

The exact TIC queries cover 201308024, its bright neighboring TIC 201308023, and exploratory leads 355767420, 281593180, 149392671 and 150104580. The complete scan found exactly 1,046,317 rows, matching the publisher's advertised count. The archive is 136,182,994 bytes and passes published MD5 `49e0b7538840dea83ecbde3a48387fd6`.

| TIC | Default adopted-period catalogue result |
|---|---|
| 201308024 | No exact TIC row. Independently verified Gaia ID 4906011953900670976 also absent from both DR2 and DR3 identifier columns. |
| 201308023 (bright neighbor) | No exact TIC row. Independently verified Gaia ID 4906017863775668992 also absent from both identifier columns. |
| 355767420 | No exact TIC row. |
| 281593180 | Adopted 17.1732 ± 1.1848 days; `flag_multiple_periods=True`, `flag_doubled_period=True`, `n_secs=2`, contributing sectors 68 and 95. Their raw catalogue periods are 8.68 and 8.50 days; the adopted result applies harmonic doubling. `flag_possible_binary=False`, `final_n_contams=0`, `flag_one_cvz_sector=False`. These flags are published classifications, not guarantees. |
| 149392671 | No exact TIC row. |
| 150104580 | No exact TIC row. |

TIC 281593180 therefore has a published long-period variability interpretation. Our sector-1/2 13-day search-boundary behavior cannot be claimed as a new rotation discovery. The published adopted period and inconsistency flag also caution against treating 17.1732 days as an unquestioned physical rotation period. The shorter approximately 2.37-day header peaks discussed in the earlier review remain a distinct published measurement.

`results.json` retains every matched table-2 column, exact queries and the query-script SHA-256. `gaia-identity-results.json` records the independent full-table exact-Gaia scan. `scripts/review_tars_aggregate.py` reproduces the checksum verification and TIC/Gaia checks. The table-4 parent catalogue has not been downloaded or checked here, and exact Gaia checks are not a complete coordinate crossmatch. Absence from the default rotation catalogue cannot establish absence from the full parent measurements, constant flux, or novelty. The primary TIC 201308024 has independently verified Gaia DR2 parallax 1.67589 ± 0.21422 mas and fails the official TARS DR2 parallax >2 mas selection, despite its larger DR3 parallax.
