# Completed QLP pilot checkpoint

The checkpoint contains the reconstructed first pilot: 24 targets, 48 sector-1/2
science light curves, their 48 sign-reversed analyses, controls, all-target plots,
catalog responses, pixel cutouts, source, tests, dependency pins and reports.
The frozen source and all 48 science-input hashes are checked before archiving.
The new clean-field cohort is excluded, including its selection script and data.
The earlier lost 88-star artifacts cannot be restored from this checkpoint.

Download and keep the four original checkpoint files and the restoration erratum
together:

- `/workspace/shared/downloads/ad-astra-qlp-pilot-20261012-core.tar.gz`
- `/workspace/shared/downloads/ad-astra-qlp-pilot-20261012-figures.tar.gz`
- `/workspace/shared/downloads/ad-astra-qlp-pilot-20261012-inventory.json`
- `/workspace/shared/downloads/ad-astra-qlp-pilot-20261012-SHA256SUMS.txt`
- `/workspace/shared/downloads/ad-astra-qlp-pilot-20261012-RESTORE-ERRATUM.md`

Core contains source, documentation, numerical outputs and public raw scientific
data; figures contains PNG/PDF and other image files. All archive member paths
start with `ad-astra/`. Git history, virtual environments, Python caches,
credentials and environment configuration secrets are excluded. The inventory
records every included relative archive path, original source path, size,
modification date and SHA256, plus archive sizes and hashes. Every member is read
back and compared with this inventory after archiving. Archive SHA256 values live
in the external inventory and checksum file: embedding an archive's own hash in
this archived document would create a circular dependency.

To restore into a new directory after downloading:

```bash
sha256sum --check ad-astra-qlp-pilot-20261012-SHA256SUMS.txt
mkdir restored-pilot
tar -xzf ad-astra-qlp-pilot-20261012-core.tar.gz -C restored-pilot
tar -xzf ad-astra-qlp-pilot-20261012-figures.tar.gz -C restored-pilot
cd restored-pilot/ad-astra
python3.12 -m venv .venv
.venv/bin/python -m pip install --disable-pip-version-check --no-cache-dir -r requirements-validated.txt
.venv/bin/python -m pip install --disable-pip-version-check --no-deps -e .
.venv/bin/python -m pytest -q
```

Use a fresh destination to avoid overwriting another working tree. These commands
require Python 3.12 and create a virtual environment in the restored checkout,
install the recorded dependencies and run tests; package downloads require
internet access. Public science artifacts are already
included, so restoration does not require opening protected later-sector data or
an ADS token. Existing provenance records retain their original paths and UTC
timestamps; do not rewrite those historical records merely because the restore
directory differs.

The original archived copy of this document incorrectly used `bash
scripts/setup.sh` after changing to the restored directory. That script explicitly
changes to `/workspace/ad-astra`, so it would install and test that canonical
checkout instead. Use the explicit commands above for restoration elsewhere.
The original archives, inventory and checksum file have been preserved unchanged;
the external `RESTORE-ERRATUM.md` records this correction. This live document now
differs from its inventoried archived version for that reason. The setup script
itself has not been modified.

The helper is `scripts/create_checkpoint.py`. Repeating it requires a new
`--prefix` because it refuses to overwrite an existing checkpoint. It verifies
integrity and archives files; it does not repeat the science or test suite.

**A local archive is not guaranteed to survive replacement of the entire cloud
environment.** Retain downloaded copies outside this environment. Saved setup
configuration alone does not preserve research data. No external publication or
repository push is performed by checkpoint creation.
