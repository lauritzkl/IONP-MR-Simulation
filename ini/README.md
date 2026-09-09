# Local simulation inputs

The public source release contains synthetic inputs in [examples](../examples/README.md).
Run them with explicit `--ini-dir`, `--ini-root`, and `--config` paths as shown
there. This directory is retained because the simulator uses `ini/` together
with `src/ionp/` to locate the repository root.

You can place your own inputs here. With no arguments, the simulator scans
only `.txt` files directly in this directory; it does not scan subdirectories.
Local inputs are ignored by Git and excluded from the publication exporter.
Study inputs in the development repository are retained locally and are not
part of the curated public source.
