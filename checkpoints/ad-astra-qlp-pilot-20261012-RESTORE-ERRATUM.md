# Restoration correction for the 12 October 2026 QLP pilot checkpoint

The original core archive contains an incorrect command in
`ad-astra/docs/checkpoint.md`. Its restoration example extracts into
`restored-pilot/ad-astra` and then invokes `bash scripts/setup.sh`. The included
setup script explicitly changes to `/workspace/ad-astra`, so that command would
install and test the canonical checkout instead of the restored one (or fail if
the canonical checkout does not exist).

Use these explicit commands instead, from the directory containing the downloaded
checkpoint files. Python 3.12 must be installed, and package installation needs
internet access:

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

Choose a fresh restoration directory to avoid overwriting another checkout.
These commands keep environment installation and tests in the restored directory.
They do not download science data or open protected later TESS sectors.

The error was verified directly in the archive. The archived document's SHA256 is
`37e71340efeff1347fafe580e852eacd02cd6ab7b038fd2bae209c0b85914631`.
Both original archives, their inventory and `SHA256SUMS.txt` remain unchanged;
this erratum is an additional file and is not included in that original inventory.
The live `docs/checkpoint.md` was corrected, while the validated setup script was
left unchanged. No science results changed and the test suite was not rerun for
this documentation correction.

Retain this erratum alongside the original four checkpoint files. A separate
`ad-astra-qlp-pilot-20261012-RESTORE-ERRATUM.sha256` verifies this erratum's bytes.
Local archives still need downloaded copies elsewhere to survive possible cloud
environment replacement.
