#!/usr/bin/env python3
"""Check publication isolation, current-file copying and reproducible archives."""

import importlib.util
import json
from pathlib import Path
import tarfile
import tempfile
import unittest


TOOL = Path(__file__).resolve().parents[1] / "tools/prepare_publication.py"
SPEC = importlib.util.spec_from_file_location("prepare_publication", TOOL)
publication = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publication)


class PublicationExportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "source"
        self.source.mkdir()
        self.write("README.md", "working-tree content\n")
        self.write("ini/README.md", "local inputs\n")
        self.write("ini/private.txt", "private study input\n")
        self.write(".git/config", "private history metadata\n")
        self.write("tools/result.bin", "generated data\n")
        self.write("publication-files.txt", "README.md\nini/README.md\n")

    def write(self, name, content):
        path = self.source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def test_curated_archive_copies_current_contents_without_private_files(self):
        self.write("README.md", "updated, uncommitted content\n")
        tree, count = publication.export(self.source, self.root / "first")
        self.assertEqual(count, 2)
        self.assertEqual((tree / "README.md").read_text(), "updated, uncommitted content\n")
        manifest = json.loads((tree / publication.MANIFEST_NAME).read_text())
        self.assertEqual({item["path"] for item in manifest["files"]},
                         {"README.md", "ini/README.md"})
        archive = next(tree.parent.glob("*.tar.gz"))
        with tarfile.open(archive) as tar:
            self.assertEqual(set(tar.getnames()), {
                f"{publication.PROJECT_NAME}/{name}" for name in
                ("README.md", "ini/README.md", publication.MANIFEST_NAME)})
        second_tree, _ = publication.export(self.source, self.root / "second")
        self.assertEqual(archive.read_bytes(), next(second_tree.parent.glob("*.tar.gz")).read_bytes())
        self.assertEqual((tree.parent / "SHA256SUMS").read_bytes(),
                         (second_tree.parent / "SHA256SUMS").read_bytes())
        with self.assertRaises(ValueError):
            publication.export(self.source, self.root / "first")
        self.assertEqual((tree / "README.md").read_text(), "updated, uncommitted content\n")

    def test_invalid_inclusions_fail_before_creating_output(self):
        for bad in ("../outside.txt", "/absolute.txt", "ini/private.txt", ".git/config",
                    "tools/result.bin", "missing.txt", "README.md\nREADME.md"):
            with self.subTest(path=bad):
                self.write("publication-files.txt", bad + "\n")
                output = self.root / "rejected"
                with self.assertRaises(ValueError):
                    publication.export(self.source, output)
                self.assertFalse(output.exists())

    def test_symlinked_source_is_rejected(self):
        target = self.source / "linked.md"
        try:
            target.symlink_to(self.source / "README.md")
        except OSError:
            self.skipTest("Symlinks are unavailable for this user/platform")
        self.write("publication-files.txt", "linked.md\n")
        with self.assertRaises(ValueError):
            publication.export(self.source, self.root / "rejected")


if __name__ == "__main__":
    unittest.main()
