# QLP stellar and neighbor metadata audit

The strongest repeatable short-period QLP lead, the roughly 0.477-day signal in the TIC 149392671 aperture, matches a **known neighboring RR Lyrae**. Gaia DR3 gives that neighbor a period of **0.4771135961 ± 0.0000108556 days**; VSX independently lists its identity and a rounded period of **0.477114 days**. The neighbor is only **29.2515 arcseconds (1.393 TESS pixels)** from the target. This supplies a specific, plausible contamination explanation. It does not support a claim of newly discovered variability in TIC 149392671.

This follow-up was performed after inspection of the pilot. It is exploratory, and uses public metadata only. No new light curves, images, or protected-season data were retrieved. Metadata retrieval dates and exact queries are preserved in the sidecars under `data/qlp-stellar-metadata/`; this report describes the 2026-10-11 retrievals. ADS was not used.

## Identities and nearby known variables

| QLP target | Exact Gaia DR3 target identity | Known neighboring variable | Separation | Published Gaia period |
| --- | --- | --- | --- | --- |
| TIC 149392671 / TYC 8883-946-1 | 4757558780144090368 | Gaia 4757558810204649856 / OGLE LMC-RRLYR-34351 / OGLE LMC598.06.000413 | 29.2515 arcsec | 0.4771135961 days, RRab |
| TIC 150104580 / TYC 8901-833-1 | 5476797645662534144 | Gaia 5476797641361944064 / OGLE LMC-RRLYR-37553 | 45.9918 arcsec | 0.5442901668 days, RRab |

Target identities are linked by exact TIC and Gaia DR3 aliases in SIMBAD, rather than by a nearest-coordinate assignment. Neighbor OGLE aliases likewise share the exact Gaia source ID in SIMBAD. VSX entries are within 0.03 arcseconds of these Gaia neighbors: OID **4318237** for the 0.477-day source and OID **4437813** for the 0.544-day source. VSX can incorporate Gaia information, so its agreement is a catalog cross-reference, not an independent photometric experiment.

The 90-arcsecond Gaia cones contain 95 sources around TIC 149392671 and 89 around TIC 150104580; independent `COUNT(*)` queries agree with the returned row counts. Neither target has a row in the queried Gaia DR3 variability summary, classifier, rotation, short-timescale, or eclipsing-binary tables. Their `NOT_AVAILABLE` flags do **not** establish constancy. The wider cones were essential: the earlier target-centered 5- or 10-arcsecond catalog checks could not find either RR Lyrae.

For TIC 150104580, the nearby RR Lyrae period is not a direct match to the roughly 1.20-day QLP peak. The arithmetic difference-frequency period is `1 / (1 / 0.5442901668 - 1) = 1.1943788068 days`, close to the sector-1 peak. This alone does not establish an alias: TESS observations are nearly continuous, and actual sampling and detrending must be tested. The target remains an unresolved low-amplitude aperture signal, with no new-variable classification justified by these metadata.

## Blend amplitudes and the known white-dwarf system

Let `r = 10**(-0.4 * (G_neighbor - G_target))`. If only the target and one neighbor contribute, both have equal aperture throughput, and the Gaia-G flux ratio equals their TESS-band flux ratio, the neighbor semiamplitude needed to produce observed fractional semiamplitude `a` is `a * (1 + r) / r`.

| Aperture / possible source | Gaia G | Conditional intrinsic semiamplitude needed in sectors 1 / 2 |
| --- | --- | --- |
| TIC 149392671 / known RR Lyrae at 29.25 arcsec | 19.255 | 24.3% / 28.9% |
| TIC 149392671 / white dwarf at 6.20 arcsec | 18.887 | 17.3% / 20.6% |
| TIC 149392671 / unrelated brighter neighbor at 24.45 arcsec, Gaia 4757555859566328320 | 15.939 | 1.17% / 1.39% |
| TIC 150104580 / known RR Lyrae at 45.99 arcsec | 19.247 | 13.2% / 12.5% |

The 0.477-day RR Lyrae has a Gaia fitted peak-to-peak G variation of **1.2013 mag** and RP variation of **1.0333 mag**. Its variability is sufficiently large that dilution is a plausible explanation for the few-hundred-ppm target-aperture signal. The other RR Lyrae has peak-to-peak G/RP variations of **0.5478/0.3273 mag**. These are different-band, different-epoch measurements; they do not predict exact TESS amplitudes. Unknown aperture throughput and color can change the estimates in either direction. Other contaminating sources can also contribute. These calculations do not spatially localize a signal or eliminate a neighbor.

The close white dwarf is **Gaia DR3 4757558814501499392**, distinct from the RR Lyrae. TIC 149392671 is the bright companion of this already published common-proper-motion system: [Raddi et al. (2022), A&A 658 A22](https://doi.org/10.26093/cds/vizier.36580022), [VizieR table J/A+A/658/A22/table2](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/658/A22/table2), system 2655. The table assigns the white dwarf to `SourceID1` and the bright TIC target to `SourceID2`. The Gaia-2016 separation is 6.20487 arcseconds; their similar proper motions support that published association. This is not a newly found white dwarf or stellar system.

## Spectroscopic checks

Gaia has no `vbroad` value for either target. Its radial-velocity constancy statistics do not indicate variability in these records: `rv_chisq_pvalue = 0.9313` for TIC 149392671 and `0.5079` for TIC 150104580. RUWE values are 0.9269 and 0.7793. These diagnostics neither classify the photometric signal nor rule out a companion or variability.

Exact Gaia identities were also recovered in [GALAH DR3](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/506/150/stars) and [GALAH DR4](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/other/PASA/42.51/allstar). DR3 records additionally match the SIMBAD GALAH identifiers.

| Target | GALAH DR3 broadening | GALAH DR4 broadening | Equatorial speed if QLP period were rotation and TIC radius correct |
| --- | --- | --- | --- |
| TIC 149392671 | 5.42 ± 2.11 km/s | 5.19 ± 1.66 km/s | 96.0 km/s for 0.477 days and 0.9051 solar radii |
| TIC 150104580 | 6.25 ± 2.34 km/s | 4.57 ± 1.59 km/s | 49.6 km/s for 1.20 days and 1.1757 solar radii |

Both DR3 records have `Flagsp=0` and `redflag=0`. Both DR4 records have `flagSp=0` and `flagSpFit=0`; TIC 150104580 also has `flagRed=0`, while TIC 149392671 has **`flagRed=3072`**, which has not been interpreted here. Its DR4 value must retain that caveat. Catalog definitions say the fitted broadening uses `sme.vsini` with `sme.vmac=0`, so it includes broadening beyond a clean projected-rotation measurement. The much smaller measured broadening is uncomfortable for a simple rapid-rotation interpretation at ordinary inclination, but cannot rule out a pole-on rotation axis or an incorrect photometric interpretation. For TIC 149392671, the matching known RR Lyrae offers the more direct explanation.

A RAVE DR6 row for the exact TIC 150104580 alias `J060928.0-641857` gives heliocentric RV 28.10 ± 1.09 km/s and primary morphology `n`; the Gaia RV is 27.151 ± 0.557 km/s. This is a limited consistency check, not a binary test. RAVE returns no row for TIC 149392671 in this 5-arcsecond query. The initial attempted GALAH table identifier `III/350/galah3` was absent from the response and was treated as unchecked; the correct GALAH catalogs above were then queried. No absent-table result was converted into a non-match.

SIMBAD also returns APOGEE/SDSS DR17 projected-velocity measurements of 1.511 and 2.501 km/s, respectively, citing `2022ApJS..259...35A`. Their source flags and uncertainties have not been examined, so these values are recorded as a literature lead only and are not used in the numerical interpretation above.

## Reproduction and limits

- `scripts/qlp_stellar_metadata.py` and `data/qlp-stellar-metadata/summary.json`: exact identities, 90-arcsecond source lists, Gaia auxiliary-table results, conditional blend calculations.
- `scripts/qlp_neighbor_metadata.py` and `data/qlp-stellar-metadata/completion/summary.json`: known-neighbor periods and exact-ID GALAH checks. Running the script reuses source- and URL-verified public metadata caches.
- `data/qlp-stellar-metadata/neighbor-rrlyrae.json`: original Gaia RR Lyrae response; schema column definitions are in `completion/rrlyrae-column-definitions.json`.
- `completion/rrlyrae-catalog-identities.json` and its raw XML sidecars: VSX identity and period checks at the neighbor positions.
- `data/qlp-prior-review-v2/white-dwarf-kinematics-149392671.tsv` and `raddi-ReadMe.txt`: published white-dwarf pair and column definitions.

FITS light curves were not accessed by this metadata audit. VizieR warning exceptions are narrow, recorded metadata defects: generated XML IDs, the non-primary `obs.field` UCD, and legacy percent/log-gravity unit strings. Other warnings, error responses, and row-cap results fail parsing. Native catalog values and raw responses remain preserved. Gaia's `epoch_g` is saved but its time origin was not established in this audit, so no cross-mission phase claim is made from that number. Catalog coverage remains incomplete, and non-detections do not prove novelty.
