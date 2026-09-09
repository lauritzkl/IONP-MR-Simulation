# Public source release

Version 0.0.1 is distributed under the [MIT license](../LICENSE), copyright
2026 Lauritz Klünder. The [third-party notices](../THIRD_PARTY_NOTICES.md)
identify the retained vendored code. The source package contains synthetic
examples and preserves the simulator's implemented behavior.

## Create the snapshot

Use Python 3.10+ from the development repository root:

```bash
python3 test/publication_export.py
python3 tools/prepare_publication.py --output-dir build/publication
```

The exporter reads [publication-files.txt](../publication-files.txt), an
explicit list of permitted files, and copies their current working-tree
contents. This includes uncommitted source and documentation changes. It does
not read Git history or automatically include unlisted files. Review the list
when adding public source files.

The output directory contains:

- `IONP-MR-Simulation/`: the source tree, ready for a fresh initial commit.
- `IONP-MR-Simulation-0.0.1.tar.gz`: the same tree in a source archive.
- `SHA256SUMS`: checksums for the archive and its inclusion manifest.

The tree's `PUBLICATION_MANIFEST.json` lists every included source file with
its byte count and SHA-256 hash. It does not list itself. Archive timestamps,
permissions and ownership are normalized, so repeat exports of identical files
on the same Python/zlib installation produce identical archives. Existing
release outputs are never replaced: move the previous release or choose a new
`--output-dir` when exporting again.

Development history, study inputs, generated results, binary grids, plots and
local settings are excluded. Their originals remain in the development
repository. The public tree retains `ini/README.md` for input-root detection;
the actual examples are under `examples/`. Publish the fresh tree, rather than
the development repository whose history still contains excluded material.

## Validate the exported source

Run the commands from inside the exported tree, or extract the archive into
a separate directory and use that tree. No development-tree files are needed.

```bash
cmake --preset linux-gcc-release
cmake --build build/linux-gcc-release --parallel
ctest --test-dir build/linux-gcc-release --output-on-failure -LE benchmark
python3 test/publication_examples.py \
  --executable build/linux-gcc-release/src/ionp/IONP \
  --source-dir . --work-dir build/publication-example-check
```

Repeat the build and CTest commands with `linux-gcc-debug`. The example check
runs each synthetic input twice, validates numerical outputs and reconstructs
the cell signal from complex MIS contributions. All test checks remain active
in release builds. It also saves grid and trajectory files for exercising the
[plotting tools](../README.md#optional-plotting).

From the directory containing `SHA256SUMS`, Linux users can run
`sha256sum -c SHA256SUMS`. Compare the manifest with the source tree before
publication; build artifacts created during validation are ignored by Git and
are not part of the manifest or archive.

## Validation status and limitations

Linux/GCC is the locally validated release platform. Native Windows/MSVC jobs
are included as advisory CI checks; they do not block the Linux release.
They request LLVM OpenMP for the existing unsigned loop indices. The runtime
source also contains implicit filesystem-path-to-string conversions in
histogram loading that require a native MSVC portability review. Windows builds
have not been executed locally; use the Linux instructions under WSL for the
validated workflow. No Windows binaries or runtime DLLs are bundled.

The separate manual benchmark workflow retains its existing 100× query-speedup
threshold. The initial audit measured about 96× in RelWithDebInfo while all
six correctness tests passed; Debug passed all seven tests. This benchmark is
sensitive to machine load and is not a total-simulation performance promise.

Existing CLI error handling, output naming, numeric parsing and MSE indexing
limitations are documented in the
[simulation reference](simulations.md#reproducibility-and-limitations). This
release preparation does not change them. Use fresh output directories and
check completion logs, stderr and the actual output time columns.

See [release-validation.md](release-validation.md) for the completed local
checks associated with this preparation. Those checks cover software behavior;
research use still requires convergence and model validation for the study.

Remote publication is a separate step. Initialize any new Git history inside
the exported source tree; do not copy the development repository's `.git`
directory into it.
