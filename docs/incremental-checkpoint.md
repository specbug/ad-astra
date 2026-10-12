# Incremental research checkpoint

The original checkpoints under `checkpoints/` are immutable. New work is preserved as an overlay relative to `checkpoints/ad-astra-research-20261012-inventory.json`, which inventories 3,911 files across the combined two-pilot core and figures archives. This avoids duplicating unchanged data while preserving new or changed code, reports, raw responses, scientific data and results.

`scripts/create_incremental_checkpoint.py` compares SHA256 identities in the baseline inventory with every current allowlisted file in `data`, `runs`, `scripts`, `docs`, `src`, `tests`, and the four root project files. Scientific files under `data/cache` remain eligible; unchanged files are omitted by identity. Virtual environments, Git metadata, Python caches, credential-like paths and unfinished `.part` files are excluded. Deleted baseline paths are explicitly listed; restoring the overlay does not silently delete them.

Before archiving, the helper verifies both required base archives and the unchanged frozen sources, manifests and 96 science inputs from the two pilots. It checks that named later-sector holdout light curves are absent from the QLP cache. This is an integrity check of the recorded pilot and cache, not a general audit of every possible external data source.

Each archive member is read back and checked against its size and SHA256. Included source files are rehashed after writing. Every existing checkpoint is checked for unchanged bytes. The full overlay stays outside Git in `/workspace/shared/downloads`; GitHub copies use ordered parts of at most 48 MiB with names such as `ad-astra-round3-20261012.tar.gz.001`. These are complete chunk files, not unfinished `.part` downloads. Their ordered concatenation is checked against the full archive size and SHA256.

The external inventory records the required base archive hashes, every overlay member, every ordered chunk, the reconstructed archive identity, absent baseline paths and the pilot integrity audit. The checksum list covers the full local archive, chunks and inventory. If only GitHub chunks were downloaded, checking the full checksum list before reconstruction will report the missing full archive; reconstruct it first.

## Build after research files are finalized

Run from the repository, choosing a fresh prefix if an earlier checkpoint exists:

```bash
.venv/bin/python scripts/create_incremental_checkpoint.py --prefix ad-astra-round3-20261012
```

Do not run while another process is editing the research files. Copy the ordered chunks, inventory and checksum list into a new `checkpoints/round3/` directory for the PR; keep the oversized full archive outside Git. The original base checkpoint files must remain unchanged.

## Portable restoration

Use an empty destination. Place both base archives and the overlay chunks, inventory and checksum list in a separate downloads directory. Verify the base archive SHA256 values against `required_base_archives` in the overlay inventory, and the inventory against its checksum-list entry. Concatenate the chunks in the inventory's numeric `order`; do not rely on filesystem directory ordering. For the default prefix:

```bash
cat ad-astra-round3-20261012.tar.gz.[0-9][0-9][0-9] > ad-astra-round3-20261012.tar.gz
sha256sum -c ad-astra-round3-20261012-SHA256SUMS.txt
mkdir restored
 tar -xzf ad-astra-research-20261012-core.tar.gz -C restored
 tar -xzf ad-astra-research-20261012-figures.tar.gz -C restored
 tar -xzf ad-astra-round3-20261012.tar.gz -C restored
```

The generated chunks have a zero-padded numeric suffix. For unusually large archives exceeding 999 chunks, concatenate by the inventory's explicit numeric order instead of the example glob. Verify all restored overlay members against their recorded `archive_path`, `size_bytes` and `sha256`. Check `baseline_paths_absent_now` before deciding whether any older files should be removed.

Install dependencies from the restored directory without invoking the cloud setup script, whose working directory is fixed to `/workspace/ad-astra`:

```bash
cd restored/ad-astra
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-validated.txt
.venv/bin/python -m pip install --no-deps -e .
.venv/bin/python -m pytest -q
.venv/bin/python -m pip check
```

An overlay is not a standalone archive: both base archives are required. Byte integrity does not establish the scientific interpretation or novelty of any result. Credentials and environment configuration are not part of the research checkpoint.
