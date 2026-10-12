# Preserved research artifacts

These immutable archives preserve the research data ignored by Git, alongside the readable source, reports and tests in this repository.

The **combined `ad-astra-research-20261012` checkpoint** contains 3,911 files: both 24-star pilots, 96 original science FITS inputs and 96 sign controls, raw catalog responses, controls, plots, declarations, source snapshots, all recorded query/analysis failures, ADS public responses, ESA Gaia epoch data, and subsequent pixel and stellar-context investigations. The external inventory records every member's path, size and SHA256. Its core and figure archives contain complementary files and must both be retained. Credentials, virtual environments, Python caches, Git history and incomplete `.part` files are excluded.

The **original `ad-astra-qlp-pilot-20261012` checkpoint** preserves the earlier first-pilot state. Keep its restoration erratum: the original archived setup example incorrectly changes to the canonical workspace path. Its bytes and original checksums are unchanged.

Each archive is below GitHub's 100 MiB individual-file limit. Archives are committed as ordinary Git objects so preservation does not depend on an external LFS store or expiring attachment. Use GitHub's **Download raw file** for binary archives, or clone this branch. Do not extract them over an existing checkout.

From this directory, verify the combined checkpoint:

```bash
sha256sum --check ad-astra-research-20261012-SHA256SUMS.txt
```

For portable restoration, see [the restoration guide](../docs/research-checkpoint.md). Historical provenance paths and timestamps remain unchanged. The earlier lost 88-star artifacts were never recovered and are not included. No new planet, star or binary has been established; see [current scientific outcomes](../docs/research-progress.md).

The archives predate this preservation PR, so their embedded source/reports retain that earlier state. The working source and PR documentation are also committed separately. Do not overwrite checkpoint files; create a new prefix for future preservation.
