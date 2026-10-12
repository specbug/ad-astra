# Third-round research overlay

These four ordered chunk files preserve 319 new or modified research files relative to the combined two-pilot checkpoint. They include the publisher-checksummed million-star TARS default catalogue, independent ASAS-SN photometry, native two-minute SPOC light curves/pixels, TGLC deblended products, fits, figures, declarations, failures, source and reports. No confirmed discovery is established.

Each chunk is at most 48 MiB, below GitHub's individual-file limit. The oversized full archive remains outside Git; reconstruct it from the ordered `.001`–`.004` files. Its SHA256 is `8aa2da9fbb4f5f87ea3d4a0153b4124b08424bf65f63a2d71e9648621fc131bf`. Every member's size/SHA256, source bytes and chunk concatenation were verified. Original checkpoints are unchanged.

From this directory:

```bash
cat ad-astra-round3-20261012.tar.gz.001 ad-astra-round3-20261012.tar.gz.002 ad-astra-round3-20261012.tar.gz.003 ad-astra-round3-20261012.tar.gz.004 > ad-astra-round3-20261012.tar.gz
sha256sum --check ad-astra-round3-20261012-SHA256SUMS.txt
```

The inventory records both required base archive hashes, the ordered chunk identities, all overlay files and the unchanged pilot-input audit. Restore the two base archives first, then overlay this archive in a fresh destination. See [full restoration instructions](../../docs/incremental-checkpoint.md) and [scientific outcomes](../../docs/research-progress.md).

The reconstructed full archive is a local restoration output; do not commit it as a single Git blob. Retain all chunks, inventory and checksums. Credentials, Git history, runtime caches and incomplete `.part` downloads are excluded. The scientific data and outputs remain under ignored `data/` and `runs/` paths when restored.
