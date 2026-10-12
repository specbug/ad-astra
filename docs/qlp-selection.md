# QLP metadata pilot selection

This is a newly retrieved, metadata-only selection following a runtime reset. It
does not reconstruct the historical 88-star manifest, and the 24 selected targets
must not be described as 24 previously unexamined stars. The full overlap with the
historical survey is unknown. TIC 149975676 appears in the conversation's earlier
injection study; TIC 167419543 appears in its earlier SIMBAD audit. Those references
are reasons to investigate prior work, not replacement evidence for the missing
artifacts.

## Current result

The script selected 24 provisional targets from the first 34 ranks of a seeded,
metadata-defined pool. All 34 QLP observation queries completed successfully. No
scientific flux was downloaded or inspected by the selection script. The manifest
remains provisional until the separate data reader and known-signal controls pass
review.

The raw TIC response contains all 1,155 filtered cone rows, with MAST reporting
1,155 filtered rows and one complete page. Local checks and exclusions leave 1,138
stars. The contemporaneous catalog snapshots contain 8,148 TOI rows, 6,445 confirmed
planet rows, and 5,138 ExoFOP CTOI rows. Seventeen TIC identities are excluded from
the cone pool; individual reasons overlap.

The provisional manifest is `data/qlp-pilot/provisional-manifest.json`, SHA-256
`d5a4aed0cd90dfaa980e44114c4db8acc9e08e8c37c4b6ae0f28e436c2a010ec`.
Raw source responses, request JSON, UTC retrieval times, hashes, complete-page
checks, the shuffled pool, and all selection/exclusion decisions are retained
under `data/qlp-pilot/`. These are retrieval records from 2026-10-10 UTC, not a public
preregistration or proof that the catalogs contain all prior scientific reports.

## Rules declared before the new flux measurements

- TIC cone: ICRS RA 90°, Dec −66.5°, radius 5°.
- TESS magnitude 8–12 inclusive; stellar radius 0.3–1.5 solar radii inclusive;
  catalog contamination ratio from zero to strictly below 0.1; absolute Galactic
  latitude strictly above 15°.
- Exclude all current TOI, confirmed-planet, or CTOI matches by usable TIC identity
  or a 5-arcsecond coordinate match, irrespective of candidate disposition. This
  deliberately conservative exclusion is not a probabilistic crossmatch.
- Explicitly exclude the earlier two modulation leads, TIC 150250788 and
  167011368, and the reference repository's TIC 370612437, 252721481, and 350883383.
- Sort remaining numeric TIC IDs, shuffle with Python's deterministic random seed
  20261009, and select the first 24 passing the observation metadata rules.
- QLP observations use MAST provenance `QLP`, numeric `target_name`, and an exact
  `_llc.fits` data URI. QLP belongs to the `HLSP` collection; filtering the collection
  to `TESS` would incorrectly produce no rows.
- New discovery data are **sectors 1 and 2 only**. Both are required; later sectors
  cannot substitute for a missing early sector. This restriction preserves the
  old protocol's holdouts despite the loss of its per-target manifest: that
  protocol required at least two earlier discovery sectors before a holdout.
- Safe QLP schema scope is sectors 1–73 or 80–93. Sectors 74–79 and 94 onward are
  conservatively excluded. Available metadata outside this range remain recorded.
- Reserve the latest eligible observing season at sector 27 or later, requiring
  at least two sectors and at least 180 days between discovery and holdout.
  Consecutive metadata observations are assigned to different seasons when their
  end-to-start gap is at least 180 days. If the latest season has only one sector,
  reject the target instead of choosing an earlier season. No holdout flux is read.

The selected TIC IDs, in seed order, are:

`167169952, 364322406, 149503544, 167419543, 167420753, 277099094, 167306038,
150065109, 149121579, 167006262, 167084404, 41334745, 150187010, 150357844,
41031619, 41108527, 150358170, 41028187, 41035042, 149975676, 176957913,
31386156, 31380486, 150389127`.

## Scope and reproducibility

Run `python scripts/select_qlp_pilot.py` from the repository root. The script uses
only Python's standard library. Existing raw snapshots are reused only if their
URL and SHA-256 agree with the recorded sidecar. Selection parameters and source
hash must match the existing declaration or the script requires a new output
directory. MAST pagination must complete with the reported row count; a repeated
page, metadata cap breach, or incomplete response fails the selection.

The first probe with a 120-second timeout failed on the relatively expensive TIC
filtered cone query. A request allowing 240 seconds completed and its result was
then retrieved through the documented script. The raw probes are auxiliary files;
the authoritative snapshot and provenance are `tic-query/response-001.json` and its
sidecar. This timeout did not remove any TIC row or bias the final 1,155-row pool.

Four targeted tests pass for the fixed discovery sectors, latest-season singleton
rejection, exclusion of sector 94 onward, and coordinate parsing/RA wraparound.

The selection step does not audit SPOC data-validation products, variable-star
catalogs, papers, image contamination, or astrophysical classification. Missing
catalog matches do not establish novelty. QLP is another reduction of TESS images,
so agreement with historical SPOC photometry would not be an independent
observation. The magnitude, radius, contamination and sector-coverage cuts also
limit completeness; no population occurrence rate is inferred.
