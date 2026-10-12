# Second QLP cohort: scoped prior-catalog review

The fixed second cohort of 24 targets has no returned entries within 10 arcseconds in ASAS-SN, VSX, or GCVS. All three requested tables were present for every target. The Gaia DR3 queries returned 27 unique sources within the target-centered 5-arcsecond union cones, with no published variability class or rotation-period values. All 27 carry `phot_variable_flag=NOT_AVAILABLE`, which does not mean constant. These are scoped catalog non-matches, not evidence of discovery.

The input is `data/qlp-overlapfield-frozen/manifest.json`, SHA-256 `d5d6e0bc825adc700cac0fd28c5f4b2d0b9fe0a532a645dcf1b1eb4ef70d6a4b`. The review completed on 2026-10-12 at 13:29:28 UTC. It ran independently alongside the light-curve analysis and was not used to revise cohort membership.

| Catalog | Query scope | Coverage and result |
| --- | --- | --- |
| ASAS-SN `II/366/catv2021` | 10-arcsecond cones; dated catalog snapshot | 24 tables present; zero rows |
| AAVSO VSX `B/vsx/vsx` | 10-arcsecond cones | 24 tables present; zero rows |
| GCVS `B/gcvs/gcvs_cat` | 10-arcsecond cones | 24 tables present; zero rows |
| Gaia DR3 classifier and rotation tables | Four six-target batches, 5-arcsecond cones | 7, 7, 7, and 6 sources returned; no classes or rotation periods |
| Independent Gaia source counts | Identical union cones, querying the source table independently | Counts 7, 7, 7, and 6 match returned rows in every batch |

Raw service responses and URL/hash provenance are retained under `data/qlp-overlapfield-prior/`. The main result is `catalog-review.json`; the count amendment was declared before its queries in `coverage-amendment.json`, and results are in `gaia-coverage-check.json`. No catalog queries failed or returned missing requested tables. The parser recorded 168 reviewed metadata warnings about generated XML IDs; unreviewed warnings, service errors, overflows, and VizieR row caps would fail the query rather than become non-matches.

`gaia-positional-review.json` records spherical separations and every source within 5 arcseconds of each target. All 24 positions have at least one match; the largest nearest-source separation is about 2.204 arcseconds. These are positional candidates only: no exact catalog-identity assignment was made, and coordinates were not propagated between epochs.

The prior selection also excluded Gaia-flagged variables within 90 arcseconds. That wider filter has incomplete sensitivity: unclassified variables and uncataloged faint contaminants remain possible. TESS pixels span roughly 21 arcseconds, so target-centered narrow catalog cones do not cover an entire photometric aperture. Catalog snapshots and indexing are incomplete; no ADS queries were performed in this particular review. Any interesting signal still requires source attribution, full literature and cross-catalog checking, and independent observations before a novelty claim.
