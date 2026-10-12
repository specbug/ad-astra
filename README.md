# Ad astra

A reproducible amateur-astronomy research workspace for public TESS observations: exploratory transit searches and unusual stellar variability, followed by catalog and literature checks. **No new planet or novel stellar source has been established.**

This workspace was reconstructed after an environment reset. The earlier 88-star analysis survives only in conversation history; see [recovery notes](docs/recovery-notes.md). Current results must be assessed from the current saved artifacts.

Setup, using the existing isolated checkout:

```bash
cd /workspace/ad-astra
bash scripts/setup.sh
.venv/bin/python -m pytest -q
```

Python 3.12 and the direct dependency versions in `requirements-validated.txt` were tested. No persistent service or ADS credential is required. Direct dependencies are pinned; this is not a complete transitive dependency lock.

The reconstructed pipeline reads legacy MIT QLP light curves with explicit timing, quality and identity checks. It verifies supplied FITS checksums against the original bytes and records retrieval provenance and SHA-256 hashes. Some QLP products supply no FITS checksum: local hashes track the downloaded content but do not replace a publisher checksum.

Two completed pilots now contain 48 distinct stars, 96 light curves, and 96 sign-reversed controls. None of their transit priorities justified promotion. The first sample's strongest repeated variability matches a known nearby RR Lyrae; pixel attribution remains inconclusive. A weaker signal remains unresolved. The second sample also prompted a separate investigation of variation in less-processed SAP measurements. Read the [current progress](docs/research-progress.md), [first-pilot report](docs/research-report.md), [second-pilot vetting](docs/qlp-overlapfield-vetting.md), [neighbor audit](docs/qlp-stellar-metadata.md), and [pixel check](docs/qlp-rrlyrae-pixels.md).

The retained follow-up lead is TIC 201308024: roughly 4.19-day sector-2 variability persists in [native two-minute pixel measurements](docs/native-spoc-201308023.md) and [Gaia-informed PSF photometry](docs/tglc-201308024.md), survives robust/day-deletion checks and is strongest among 24 field comparisons at that frequency. Independent [ASAS-SN photometry](docs/ground-validation-201308024.md) is too noisy for confirmation. Its absence from the [million-star adopted-period catalogue](docs/tars-aggregate-review.md) does not establish novelty. Rotation, binarity and intrinsic attribution remain unconfirmed; [Gaia/GALAH context](docs/qlp-sap-stellar-201308024.md) and position-dependent aperture effects require further investigation.

Independent ESA Gaia epoch observations show strong variability in both known RR neighbors. ADS literature access works, with public query snapshots and scoped limitations preserved. The [archive assessment](docs/archive-expansion.md) records actual coverage and access failures for ESA, ISRO and other observatories; [ADS review](docs/ads-literature-review.md) explains why identifier nonmatches do not establish novelty.

Science retrieval uses only sectors 1 and 2. Later seasons are reserved and must remain unopened: the lost original manifest prevents auditing every historical holdout. Target overlap with the old cohort is unknown. Selection records, frozen configurations, raw public metadata, all failures, known-signal controls, and inverted-flux controls are retained.

See [QLP ingestion](docs/qlp-ingestion.md), [selection](docs/qlp-selection.md), and the [research protocol](docs/research-protocol.md). Run outputs live under `runs/`; downloaded observations and catalog snapshots live under `data/`. Both are ignored by Git and require separate preservation.

The search ranks possibilities; conditional signal-to-noise values and positive-depth event windows do not establish detection significance. Public catalog nonmatches do not establish novelty. Known π Mensae controls recovered 18/18 deeper injected signals but missed the real shallow planet in one of two sectors, so nondetections cannot establish absence.
