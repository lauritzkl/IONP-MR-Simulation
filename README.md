# IONP MR Simulation

A C++ simulator for magnetic resonance signal formation around iron-oxide
nanoparticles (IONPs). It places particles, optional aggregates and spherical
cells in a simulation volume, propagates diffusing free and optional bound
protons, and calculates FID, spin-echo (SE), or multi-spin-echo (MSE) signals.
Outputs include signal magnitude, particle positions, and optional trajectories
and magnetic-field grids.

Start with the synthetic [examples](examples/README.md). The
[simulation reference](docs/simulations.md) explains configuration keys, geometry
formats, field solvers, importance sampling and output interpretation. The
examples demonstrate operation; their small sample counts and short runs are
not evidence of convergence or experimental validation.

For reading or extending the implementation, the [source guide](docs/source-guide.md)
maps the simulation pipeline, units, ownership, random streams and weighting rules
to the source files.

## Build

The simulator requires a C++17 compiler and OpenMP. The preset commands below
require CMake 3.21 or newer and Ninja. Configure and build do not download
third-party source code. Python is needed for the publication checks and the
optional plotting tools; it is not required to build the simulator.

From the repository root on Linux with GCC:

```bash
cmake --preset linux-gcc-release
cmake --build build/linux-gcc-release --parallel
```

This builds the simulator and test executables. Use `linux-gcc-debug` in both
commands for a debug build. Clang presets are also provided, with a compatible
OpenMP installation required.

For a manual Makefile build, CMake 3.16 or newer is supported:

```bash
cmake -S . -B build-make -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build-make --parallel
```

On Windows, the locally validated route is GCC inside WSL using the Linux
commands above. Native MSVC builds are experimental and have advisory CI jobs;
existing histogram path conversions still need a native portability review.

For native build evaluation, use Visual Studio 2022 C++ tools, a Windows SDK,
and CMake 3.30+. The existing unsigned OpenMP loops require the LLVM OpenMP
runtime. Configure with:

```text
cmake -S . -B build-windows -G "Visual Studio 17 2022" -A x64 -DOpenMP_RUNTIME_MSVC=llvm
cmake --build build-windows --config RelWithDebInfo --parallel
```

If using a Windows MSVC Ninja preset from an x64 developer prompt, also pass
`-DOpenMP_RUNTIME_MSVC=llvm`. Running the result requires the matching installed
LLVM OpenMP DLL directory on `PATH`; the Windows CI workflow locates it under
Visual Studio's `VC/Redist/MSVC` directory. These native commands have not been
validated locally. See [Microsoft's OpenMP documentation](https://learn.microsoft.com/en-us/cpp/build/reference/openmp-enable-openmp-2-0-support?view=msvc-170)
and [CMake's runtime option](https://cmake.org/cmake/help/latest/module/FindOpenMP.html#variable:OpenMP_RUNTIME_MSVC)
for prerequisites. Platform presets are in [CMakePresets.json](CMakePresets.json).

## First Run

Run a single synthetic FID simulation from the repository root:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples \
  --ini-root examples \
  --config examples/fid.txt \
  --output-dir data/fid \
  --log-dir data/log/fid \
  --seed 123
```

The simulator creates the output directories. Expect one `*_signal_FID_*.dat`
file and one `*_IONPs_*.dat` file under `data/fid`, plus a log under
`data/log/fid`. Signal columns are **magnitude, time in seconds**; particle
columns are **x, y, z in meters**. Check the log and stderr as well as the output
files: the current executable can report a failed simulation and still return
exit status zero.

| Config | Demonstrates |
| --- | --- |
| [fid.txt](examples/fid.txt) | FID with three fixed synthetic IONP positions. |
| [se.txt](examples/se.txt) | A single spin echo with the same geometry. |
| [mse.txt](examples/mse.txt) | Multiple spin echoes, with full and echo-sampled signals. |
| [cell_importance.txt](examples/cell_importance.txt) | One cell, membrane crossings, and overlapping importance proposals with automatic weights. |

Replace `fid` in the command with the desired example name. Each example uses
the included [positions file](examples/IONP_pos/positions.txt); study inputs
are not needed. Keep each run in its own output and log directories.

## Run Options

`--config` paths are relative to the current working directory. Geometry paths
are controlled separately: `--ini-dir` supplies `IONP_pos` and `Agg_pos`, while
`--ini-root` supplies `Cell_pos`, `Cell_load`, `Histogram`, and
`Radius_Histogram`. See the [CLI reference](docs/simulations.md#cli-reference)
for root inference and all supported options.

Frequently used options are:

- `--no-log`: print stdout diagnostics to the terminal.
- `--record-positions`: save proton trajectories and reflection events.
- `--save-region-signals`: save importance-sampling diagnostics. For overlapping
  proposals, the complex `mis_contribution` files sum to the free-proton signal;
  proposal `mean` files describe their individual samples.
- `--save-bgrid <file>`: export the interpolation grid when grid calculation is
  enabled. Use a separate filename for each run.
- `--seed <integer>`: seed placement and diffusion. For repeatable runs, also
  preserve the inputs, build environment, and thread settings.

Historical (legacy) IONP placement is the default. You can omit
`--legacy-placement`; existing commands that include it behave the same way.

Without `--config`, the program runs all `.txt` configs directly inside
`--ini-dir`, sorted by path. With no options, it looks for `<repo>/ini`; the
public source package keeps that directory as an empty input location. The
included examples are selected explicitly as above. The retained `--example`
option generates a separate MSE execution example under `Testing/example`;
its meter-scale geometry is a software demonstration.

Custom file selectors such as `IONP_position_file` choose one `.txt` or `.dat`
leaf filename from its input folder. Their activation rules, cell-position/load
pairing requirements and legacy folder scans are described under
[custom geometry files](docs/simulations.md#custom-geometry-files).

## Optional Plotting

Use Python 3.10 or newer with NumPy and Matplotlib. On Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-plotting.txt
```

On Windows, create the environment with `py -3 -m venv .venv`, activate
`.venv\Scripts\activate` in Command Prompt, and use the same pip command.

Export and plot a field grid:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples \
  --ini-root examples \
  --config examples/fid.txt \
  --output-dir data/fid-grid \
  --log-dir data/log/fid-grid \
  --seed 123 \
  --save-bgrid data/fid-grid/bgrid.bin
python tools/plot_bgrid.py data/fid-grid/bgrid.bin --output data/fid-grid/bgrid.png
```

For trajectories, run `examples/cell_importance.txt` with
`--record-positions`. Then pass one actual trajectory filename from the log:

```bash
python tools/plot_trajectories.py \
  --file "path/to/the/generated_ProtonTrajectories_file.dat" \
  --radius 1e-6 --center 0 0 0 --proton 0
```

Replace the placeholder with the complete filename, including its numeric
prefix. The trajectory tool accepts one file and opens an interactive figure;
it requires a graphical display. For grid plots on a machine without a display,
set `MPLBACKEND=Agg` and use `--output`.

## Tests and Reproducibility

After building all targets, run the correctness suite:

```bash
ctest --test-dir build/linux-gcc-release --output-on-failure -LE benchmark
```

Run the separate performance check on a release build:

```bash
ctest --test-dir build/linux-gcc-release --output-on-failure -L benchmark
```

For a Visual Studio build, use `--test-dir build-windows -C RelWithDebInfo`. Tests
cover parsing, geometry selection, cell-load radii, reflection, importance
sampling and field evaluation. The benchmark has a hardware-sensitive speedup
threshold; see [verification](docs/simulations.md#verification) for its limits.

For research runs, save the config, auxiliary inputs, seed, source revision,
build options and environment alongside the results. Use one config per
process when comparing seeded runs: a batch shares a random stream, so earlier
configs affect later ones. Check convergence with sample count, time step and
grid resolution, and use independent seeds to assess sampling variability.
See [reproducibility and limitations](docs/simulations.md#reproducibility-and-limitations)
for operational and numerical details.

## Preparing a Public Source Release

Follow [the publication guide](docs/publication.md) to export the source and
curated examples while keeping study inputs and generated results local. The
export is a working-tree snapshot, so review its manifest and run its checks
before publishing it.

## License

[MIT](LICENSE), Copyright (c) 2026 Lauritz Klünder. Retained third-party
components are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## AI Assistance Disclosure

The original simulation design and core functionality were developed without
generative AI. OpenAI Codex was subsequently used to assist with substantial
modifications to the codebase. The project authors retain full responsibility
for the final implementation and simulation results.
