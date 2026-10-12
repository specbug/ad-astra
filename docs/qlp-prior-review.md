# QLP pilot prior-catalog review

The reviewed 24-target QLP pilot has no entries within 10 arcseconds in the returned VSX, GCVS, or ASAS-SN variable-star tables, and no classification or published rotation period in the returned Gaia DR3 5-arcsecond cones. **This is a scoped catalog result, not evidence that the aperture signals are new.** A wider follow-up found a known RR Lyrae 29.25 arcseconds from TIC 149392671 with a period matching the strongest short-period pilot lead; see [the stellar metadata audit](qlp-stellar-metadata.md).

The final machine-readable target review is `data/qlp-prior-review-v2/catalog-review-v3-final.json`, completed at 2026-10-11 11:45:19 UTC. Its input is the final provisional manifest, SHA-256 `e9e71df0ad45c05499a112272976aa2d5755058818afd4d2b11907582f7ef900`, containing the same 24 targets as the frozen pilot. This review followed light-curve inspection and is exploratory prior checking, not a pre-registered selection filter.

| Check | Scope and stored version | Result |
| --- | --- | --- |
| AAVSO VSX via VizieR `B/vsx/vsx` | 10-arcsecond cones; response labeled 2026-10-05 | All 24 tables present; no returned rows |
| GCVS `B/gcvs/gcvs_cat` | 10-arcsecond cones; GCVS 5.1, October 2020 | All 24 tables present; no returned rows |
| ASAS-SN `II/366/catv2021` | 10-arcsecond cones; response describes the 2021 version and a 2024 database extraction | All 24 tables present; no returned rows |
| Gaia DR3 classifier and rotation tables | 5-arcsecond cones around manifest positions | Four successful batches; 28 Gaia sources total; no returned class or rotation-period values |
| Gaia identity review | Exact TIC/Gaia aliases in the saved SIMBAD audits plus spherical separations | 21 targets with unique exact alias evidence; three position-only matches |

The three position-only cases are TIC **150187010**, **150357844**, and **176984190**, with nearest Gaia sources 5476898800729228032, 5478079126460256640, and 5281974626279231104 at approximately 0.113, 0.213, and 0.160 arcseconds. Their proximity supports investigation but is not an exact-identity proof. The identity review is saved in `data/qlp-prior-review-v2/gaia-identity-review-v3.json`. Coordinates were not propagated between catalog epochs, and additional Gaia rows are neighbors rather than extra target detections.

Earlier review outputs are retained. The VizieR parser initially rejected legacy metadata warnings. The corrected reader now records generated-XML-ID warnings and the known percent-unit warning, without altering science values, and rejects unreviewed warnings, service errors, overflows, or row caps. Missing requested tables are reported as unchecked. The final review contains all three requested tables for every target. The final v3 target list supersedes the earlier v2 review, whose membership differs; old summaries must not be used as final-cohort evidence.

The selected targets also passed the catalog and SIMBAD checks described in [the selection report](qlp-selection.md), including the preserved exclusion lists and fetched reference titles. SIMBAD bibliography is incomplete; article titles alone are not a comprehensive period search. ADS was omitted at the user's request. ASAS-SN and GCVS snapshots are dated, and Gaia does not classify every variable star. TESS pixels are about 21 arcseconds wide, so checking only target-centered 5- or 10-arcsecond cones leaves much of a photometric aperture unexamined. This review cannot establish a first discovery, rule out previously reported periods in unsearched literature, or attribute an aperture signal to the nominal target.

The wider metadata audit gives a concrete example of that limitation: VSX OID **4318237**, Gaia DR3 **4757558810204649856**, also known as **OGLE LMC-RRLYR-34351**, lies outside the original target cone and has published period **0.477114 days**. A separate known RR Lyrae lies 45.99 arcseconds from TIC 150104580. These objects are documented as contamination checks; catalog and pixel/aperture evidence must be assessed together before promoting a target.
