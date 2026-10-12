# Combined research checkpoint

The combined checkpoint preserves both completed 24-target TESS QLP pilots:
48 distinct targets, 96 sector-1/2 science inputs and 96 sign-reversed analyses.
It includes controls, selection metadata (including the failed coverage proposal),
public ADS responses, ESA/Gaia archive assessments and available epoch products,
SAP follow-up diagnostics, source, tests, dependency pins and current reports.
The earlier lost 88-star artifacts are absent. Named later-sector QLP holdouts
remain protected; their metadata are included but their light curves are unopened.

The helper `scripts/create_research_checkpoint.py` audits both run declarations,
frozen manifests/configurations, control-result hashes and current plus snapshotted
screening source. It verifies all 96 original science-input hashes and requires
48 distinct targets. It refuses a cached protected QLP URI. Scientific FITS inputs
in `data/cache` are retained; runtime/Python caches, virtual environments, Git
history, credentials and partial files are excluded. Every archive member is read
back and compared with its inventory; all included source files are rehashed after
archive creation to detect concurrent changes. SHA256 proves file identity, not
scientific validity, significance or novelty.

The output prefix is `ad-astra-research-20261012` under
`/workspace/shared/downloads`, with `-core.tar.gz`, `-figures.tar.gz`,
`-inventory.json` and `-SHA256SUMS.txt` suffixes. These outputs are created only
after the final reports and diagnostics are ready. The helper refuses to overwrite
existing outputs. First-pilot archives remain unchanged.

Download and retain all four combined files outside this environment. Restore to
a fresh directory using the explicit commands below; Python 3.12 and internet
access for package installation are required.

```bash
sha256sum --check ad-astra-research-20261012-SHA256SUMS.txt
mkdir restored-research
tar -xzf ad-astra-research-20261012-core.tar.gz -C restored-research
tar -xzf ad-astra-research-20261012-figures.tar.gz -C restored-research
cd restored-research/ad-astra
python3.12 -m venv .venv
.venv/bin/python -m pip install --disable-pip-version-check --no-cache-dir -r requirements-validated.txt
.venv/bin/python -m pip install --disable-pip-version-check --no-deps -e .
.venv/bin/python -m pytest -q
```

The first checkpoint's original archived restoration instructions incorrectly
used `bash scripts/setup.sh` after changing to the restored directory. That script
changes explicitly to `/workspace/ad-astra`; use the commands above when restoring
elsewhere. The original archives and their hashes remain preserved, with a separate
`ad-astra-qlp-pilot-20261012-RESTORE-ERRATUM.md` beside them.

Historical provenance records retain their original paths and UTC timestamps.
Do not rewrite them to match a restoration directory. Existing archived public
ADS responses need no token to read. Credentials are neither required for local
restoration nor included in this checkpoint. Archive creation runs no new science,
network requests or tests and performs no push or external publication. Local
archives can disappear if the cloud environment is replaced, so download copies.
