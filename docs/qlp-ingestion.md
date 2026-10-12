# Reconstructed QLP ingestion and real control

This is a fresh reconstruction on 2026-10-10 UTC after the earlier runtime's
uncommitted research files were lost. It does not restore or reproduce the old
88-star survey. This experiment downloads only the known pi Men control
(TIC 261136679) in QLP sectors 1 and 27.

The new reader and archive client successfully ingest both real files. The
new, fixed search recovers the known period as its strongest peak in sector 1,
but **fails the same test in sector 27**. The latter failure is retained; the
reader passing its tests does not establish scientific sensitivity.

## Fresh observations

The NASA Exoplanet Archive `pscomppars` table gives pi Men c a period of
6.2678399 days, a transit midpoint of BJD TDB 2458425.789204, and a duration of
2.952024 hours. The corresponding `ps` default row attributes these to
Kunovac Hodžić et al. (2021). The raw catalogue responses are saved alongside
the results. These are known-planet parameters, not a new finding.

The search uses `KSPSAP_FLUX` with its own errors, quality-zero finite positive
samples, a 0.5-day running median in segments split by gaps longer than
0.5 days, and no outlier clipping. It searches periods 0.5–13 days, capped at
half the measured baseline, with five box durations. The frequency spacing
bounds the nearest-grid accumulated timing drift by one quarter of the
shortest tested duration. The strongest peak must be within 1% of the known
period to satisfy this limited recovery check.

| Fresh control | Accepted samples | Median cadence | Strongest period | Within 1% of known period? |
|---|---:|---:|---:|---|
| Sector 1 | 993 | 1800 seconds | 6.264816 days | Yes |
| Sector 27 | 3302 | 600 seconds | 6.892805 days | No |

At the fixed published period **and phase**, descriptive local median depths
are approximately 297, 248, 319 and 339 ppm in sector 1, and 328, 83, 234 and
54 ppm in sector 27. These are unweighted descriptive values, not independent
significances. Sector 27 contains much larger localized excursions. Its
strongest box-search peak and its phase-optimized fit at the known period
must not be mistaken for recovery at the published transit times.

The first search output and its exact script are preserved as
`original-search-results.json` and `original-search-script.py`. The subsequent
run adds the published-ephemeris diagnostics and writable plotting caches;
its period search is unchanged and produces identical strongest periods.

## Reader scope and errors

`src/adastra/qlp.py` supports only the legacy QLP schema in sectors 1–73 and
80–93. It explicitly checks the sector, `TIMESYS=TDB`, day-valued `TIMEUNIT`,
and the BJD reference. Actual control headers use `BJDREFR`; the reader also
accepts `BJDREFF`, but rejects conflicting values. It requires an explicit
integer and fractional reference instead of guessing a time origin.

Both real files contain `SAP_FLUX`, `KSPSAP_FLUX`, `KSPSAP_FLUX_ERR`, `TIME`,
and `QUALITY`. Neither contains `SAP_FLUX_ERR`. A SAP read therefore has
`flux_error=None`; it never borrows corrected-flux uncertainties. The reader
preserves every sample and quality value. Filtering is a separate explicit
`good_mask()` operation. `read_qlp_lightcurve()` provides a mapping interface
whose `time` is explicitly BTJD, or BJD TDB minus 2457000.

Neither real control has `CHECKSUM` or `DATASUM` FITS cards. The archive client
verifies FITS structure and saves a SHA-256 digest with each download. The
digest identifies subsequent cached bytes; **it is not a publisher-supplied
checksum**. The client verifies checksum cards when a file supplies them.

MAST exposes these products as `obs_collection=HLSP`,
`provenance_name=QLP`, with numeric `target_name` strings. Restricting the
query to the TESS collection would miss them. Metadata for other sectors
may be returned by the query, but the script downloads only sectors 1 and 27.

## Release boundaries

The [MAST QLP page](https://archive.stsci.edu/hlsp/qlp), retrieved afresh and
archived with this experiment, says:

- Sectors 74–79 had barycentric timing errors of up to about ±3 minutes;
  corrected products were redelivered as version 02. This initial reader
  conservatively excludes that range even when a corrected version exists.
- QLP switched to the TESS-Gaia Light Curve methodology in sector 94. Those
  files require separate schema and behavior validation.
- Sectors 99–104 have a separate timing-correction issue. The page marks
  affected sectors that have not yet been redelivered. They are excluded.

The machine-readable outputs are in
`runs/reconstruction/qlp-control/results.json`, with the declaration, raw
MAST and NASA catalogue responses, FITS download provenance, and two PNG/PDF
figures beside it. Generated data are not source-controlled by default.

## Reproduce and test

```bash
bash scripts/setup.sh
.venv/bin/python -m pytest -q tests/test_qlp.py
.venv/bin/python scripts/qlp_control.py
```

The fresh targeted test run passed 12 tests covering the time-reference
variants, conflicting references, quality preservation, absent SAP errors,
required corrected-flux errors, time scale, and excluded sector ranges.
Tests use synthetic FITS files; the two public controls provide the separate
real-file check. Additional control recovery and false-positive experiments
are needed before using this search as a discovery filter.
