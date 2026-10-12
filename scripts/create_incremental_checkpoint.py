#!/usr/bin/env python3
"""Create a verified research overlay and GitHub-sized multipart copies, without networking."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
from create_research_checkpoint import audit_pilot, digest, now, selected_files, verify_archive


def split_archive(archive: Path, chunk_bytes: int = 48 * 1024 * 1024) -> list[dict]:
    """Split losslessly; refuse to replace previous parts and verify concatenation."""
    if chunk_bytes < 1:
        raise ValueError("chunk_bytes must be positive")
    if list(archive.parent.glob(archive.name + ".[0-9]*")):
        raise FileExistsError("Multipart output already exists")
    chunks = []
    with archive.open("rb") as source:
        while content := source.read(chunk_bytes):
            destination = archive.with_name(archive.name + f".{len(chunks) + 1:03d}")
            with destination.open("xb") as output:
                output.write(content)
            chunks.append({"order": len(chunks) + 1, "filename": destination.name,
                           "size_bytes": len(content), "sha256": digest(destination)})
    combined = hashlib.sha256()
    total = 0
    for chunk in chunks:
        path = archive.parent / chunk["filename"]
        if digest(path) != chunk["sha256"]:
            raise ValueError("Chunk changed during verification")
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                combined.update(block)
                total += len(block)
    if combined.hexdigest() != digest(archive) or total != archive.stat().st_size:
        raise ValueError("Concatenated multipart bytes differ from archive")
    return chunks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--baseline", type=Path, default=Path("checkpoints/ad-astra-research-20261012-inventory.json"))
    parser.add_argument("--output", type=Path, default=Path("/workspace/shared/downloads"))
    parser.add_argument("--prefix", default="ad-astra-round3-20261012")
    args = parser.parse_args()
    if Path(args.prefix).name != args.prefix:
        parser.error("prefix must be a filename component")
    root = args.root.resolve()
    baseline_path = args.baseline if args.baseline.is_absolute() else root / args.baseline
    baseline = json.loads(baseline_path.read_text())
    prior = {entry["archive_path"]: entry for entry in baseline["files"]}
    if len(prior) != len(baseline["files"]):
        raise ValueError("Duplicate baseline paths")
    preserved = {str(p.relative_to(root)): digest(p) for p in (root / "checkpoints").rglob("*") if p.is_file()}
    required = []
    for entry in baseline["archives"]:
        path = root / "checkpoints" / Path(entry["path"]).name
        if digest(path) != entry["sha256"] or path.stat().st_size != entry["size_bytes"]:
            raise ValueError(f"Baseline archive mismatch: {path}")
        required.append({"filename": path.name, "sha256": entry["sha256"], "size_bytes": entry["size_bytes"]})
    audit = audit_pilot(root)
    entries = []
    current = set()
    for path in selected_files(root):
        name = "ad-astra/" + path.relative_to(root).as_posix()
        current.add(name)
        sha = digest(path)
        if name in prior and prior[name]["sha256"] == sha:
            continue
        entries.append({"archive_path": name, "source_path": str(path), "size_bytes": path.stat().st_size,
                        "sha256": sha, "change": "modified" if name in prior else "new"})
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    archive = output / (args.prefix + ".tar.gz")
    inventory_path = output / (args.prefix + "-inventory.json")
    sums = output / (args.prefix + "-SHA256SUMS.txt")
    partial = archive.with_name(archive.name + ".part")
    if any(p.exists() for p in (archive, partial, inventory_path, sums)):
        raise FileExistsError("Checkpoint already exists; choose a new prefix")
    with tarfile.open(partial, "w:gz", compresslevel=6) as tar:
        for entry in entries:
            path = Path(entry["source_path"])
            content = path.read_bytes()
            if len(content) != entry["size_bytes"] or hashlib.sha256(content).hexdigest() != entry["sha256"]:
                raise ValueError(f"Source changed while archiving: {path}")
            info = tarfile.TarInfo(entry["archive_path"])
            info.size = len(content)
            info.mode = 0o755 if path.stat().st_mode & 0o111 else 0o644
            info.mtime = int(path.stat().st_mtime)
            tar.addfile(info, io.BytesIO(content))
    verified = verify_archive(partial, entries)
    partial.rename(archive)
    chunks = split_archive(archive)
    for entry in entries:
        if digest(Path(entry["source_path"])) != entry["sha256"]:
            raise ValueError("Source changed before checkpoint completion")
    for relative, sha in preserved.items():
        if digest(root / relative) != sha:
            raise ValueError(f"Immutable checkpoint changed: {relative}")
    inventory = {"created_utc": now(), "scope": "New and SHA256-modified allowlisted research files relative to the combined two-pilot checkpoint",
                 "baseline_inventory": {"filename": baseline_path.name, "sha256": digest(baseline_path), "file_count": len(prior)},
                 "required_base_archives": required, "preserved_checkpoint_sha256": preserved,
                 "archive": {"filename": archive.name, "sha256": digest(archive), "size_bytes": archive.stat().st_size},
                 "chunks": chunks, "chunk_limit_bytes": 48 * 1024 * 1024,
                 "verified_members": verified, "files": entries,
                 "baseline_paths_absent_now": sorted(set(prior) - current),
                 "restore_semantics": "Extract both base archives, then concatenate chunks in recorded order, verify archive SHA256, and extract overlay. Absent baseline paths are recorded, not deleted automatically.",
                 "pilot_integrity_audit": audit,
                 "verification": "Every member size/SHA256 and included source checked; concatenated chunks match full archive size/SHA256; preexisting checkpoints unchanged.",
                 "limitations": ["No secrets, virtual environment or Git history included. Cache scientific artifacts remain eligible.", "SHA256 verifies byte identity, not scientific novelty."]}
    inventory_path.write_text(json.dumps(inventory, indent=2) + "\n")
    paths = [archive] + [output / chunk["filename"] for chunk in chunks] + [inventory_path]
    sums.write_text("".join(f"{digest(path)}  {path.name}\n" for path in paths))
    print(json.dumps({"archive": str(archive), "inventory": str(inventory_path), "checksums": str(sums),
                      "files": verified, "chunks": chunks}, indent=2))


if __name__ == "__main__":
    main()
