#!/usr/bin/env python3
"""Export an explicitly selected working-tree snapshot without development history."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile


PROJECT_NAME = "IONP-MR-Simulation"
PROJECT_VERSION = "0.0.1"
MANIFEST_NAME = "PUBLICATION_MANIFEST.json"


def publication_paths(source: Path) -> list[str]:
    """Reject missing files, symlinks, and private/generated paths before copying."""
    paths = []
    seen = set()
    forbidden = {".git", ".codex", ".agents", ".vscode", "testing", "data",
                 "build", ".venv", "__pycache__"}
    for line in (source / "publication-files.txt").read_text(encoding="utf-8").splitlines():
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        relative = PurePosixPath(name)
        if (relative.is_absolute() or ".." in relative.parts or "\\" in name
                or ":" in name or str(relative) != name):
            raise ValueError(f"Invalid publication path: {name}")
        if name.casefold() in seen:
            raise ValueError(f"Duplicate publication path: {name}")
        seen.add(name.casefold())
        if any(part.casefold() in forbidden or part.casefold().startswith("build-")
               for part in relative.parts):
            raise ValueError(f"Excluded directory in publication path: {name}")
        if relative.parts[0].casefold() == "ini" and name != "ini/README.md":
            raise ValueError(f"Study inputs cannot be published: {name}")
        if relative.suffix.casefold() in {".bin", ".png", ".pyc"}:
            raise ValueError(f"Generated artifact cannot be published: {name}")
        current = source
        for part in relative.parts:
            current = current / part
            if current.is_symlink():
                raise ValueError(f"Symlink cannot be published: {name}")
        if not current.is_file():
            raise ValueError(f"Required publication file is missing: {name}")
        paths.append(name)
    if not paths:
        raise ValueError("The publication inclusion list is empty")
    return sorted(paths)


def write_archive(tree: Path, archive: Path) -> None:
    """Normalize timestamps, ownership and modes for reproducible archives."""
    with archive.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as tar:
                for path in sorted(tree.rglob("*")):
                    if not path.is_file():
                        continue
                    data = path.read_bytes()
                    info = tarfile.TarInfo(f"{PROJECT_NAME}/{path.relative_to(tree).as_posix()}")
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = 0
                    tar.addfile(info, io.BytesIO(data))


def export(source: Path, output: Path) -> tuple[Path, int]:
    paths = publication_paths(source)
    if output == source:
        raise ValueError("Choose a separate publication output directory")
    archive_name = f"{PROJECT_NAME}-{PROJECT_VERSION}.tar.gz"
    artifact_names = (PROJECT_NAME, archive_name, "SHA256SUMS")
    for name in artifact_names:
        target = output / name
        if target.exists() or target.is_symlink():
            raise ValueError(f"Publication output already exists: {target}; choose a new output directory")
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".publication-", dir=output) as scratch:
        staging = Path(scratch)
        tree = staging / PROJECT_NAME
        tree.mkdir()
        records = []
        for name in paths:
            data = (source / name).read_bytes()
            target = tree / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            target.chmod(0o644)
            records.append({"path": name, "bytes": len(data),
                            "sha256": hashlib.sha256(data).hexdigest()})
        manifest = {"format_version": 1, "project": PROJECT_NAME,
                    "version": PROJECT_VERSION, "files": records}
        (tree / MANIFEST_NAME).write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        archive = staging / archive_name
        write_archive(tree, archive)
        sums = []
        for relative in (archive_name, f"{PROJECT_NAME}/{MANIFEST_NAME}"):
            digest = hashlib.sha256((staging / relative).read_bytes()).hexdigest()
            sums.append(f"{digest}  {relative}\n")
        (staging / "SHA256SUMS").write_text("".join(sums), encoding="utf-8")
        for name in artifact_names:
            # Refuse to replace files even if a second exporter finished meanwhile.
            target = output / name
            if target.exists() or target.is_symlink():
                raise ValueError(f"Publication output appeared during export: {target}")
            (staging / name).rename(target)
    return output / PROJECT_NAME, len(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("build/publication"),
                        help="Destination for the source tree, archive and checksums")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    try:
        tree, count = export(source, args.output_dir.resolve())
    except (OSError, ValueError) as error:
        parser.exit(1, f"Publication export failed: {error}\n")
    print(f"Exported {count} files to {tree}")
    print(f"Archive and SHA256SUMS: {tree.parent}")


if __name__ == "__main__":
    main()
