# Alternate TESS product availability for TIC 201308024

A metadata-only MAST audit on 2026-10-12 found useful alternate reductions in the already inspected sectors 1 and 2. Requests explicitly restricted `sequence_number` to 1/2. Exact TIC and `TIC <ID>` names were queried together. Every raw request, response, HTTP provenance, response SHA256, paging count and completion record is retained in `data/alternate-tess-201308024`; no later-sector products were queried or downloaded by this audit.

| Star | Provider | S1 observations | S2 observations |
|---|---|---:|---:|
| TIC 201308024 | QLP | 1 | 1 |
| TIC 201308024 | TASOC | 2 | 2 |
| TIC 201308024 | GSFC-ELEANOR-LITE | 1 | 1 |
| TIC 201308024 | TGLC | 1 | 1 |
| TIC 201308024 | T16 | 1 | 1 |
| TIC 201308023 | Native SPOC | 1 | 1 |
| TIC 201308023 | TESS-SPOC | 1 | 1 |
| TIC 201308023 | QLP | 1 | 1 |
| TIC 201308023 | TASOC | 2 | 2 |
| TIC 201308023 | GSFC-ELEANOR-LITE | 1 | 1 |
| TIC 201308023 | TGLC | 1 | 1 |
| TIC 201308023 | T16 | 1 | 1 |
| TIC 201308023 | TARS | 1 | 1 |

The complete all-provider exact-name queries contain 12 observations/16 products for the primary and 18 observations/28 products for the bright neighbor. A separate native `obs_collection=TESS` query contains zero observations for the primary and two observations/four products for the neighbor. This establishes absence only under these exact names and this collection filter. It does not establish global archive absence.

The native neighbor SPOC observations are 120-second cadence, public, and each offers a light curve and target-pixel file: S1 obsid 60843726, S2 obsid 60879525. TPF sizes are 49,184,640 and 48,355,200 bytes. The neighbor is approximately 92 arcseconds from the primary; actual stamp inclusion must be verified, rather than assumed from typical stamp sizes. A separately declared [native SPOC follow-up](native-spoc-201308023.md) performs that check and fixed-period fits.

Both attempted 5-arcsecond `Mast.Caom.Filtered.Position` supplement queries returned archive `ERROR`, not zero matches. Their original requests/responses and explicit failed states are retained in `supplement-results.json`; no narrow-cone completeness claim is made. Exact-name queries and product paging completed without truncation. Alternate reductions use the same TESS instrument and largely the same underlying images, so agreement provides processing/cadence robustness rather than independent instrumental validation.

TGLC's Gaia-source-identified primary products (Gaia 4906011953900670976) are public S1/S2 light curves at obsids 108443670 and 113356767. These are a promising additional source-deblending check. Availability is not a prior variability or novelty clearance, and no period was optimized in this metadata audit.
