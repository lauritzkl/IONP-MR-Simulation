# Synthetic examples

These small inputs demonstrate the simulator without requiring study data.
The particle coordinates are synthetic and also used by the runtime regression
fixtures. Each example uses three fixed particles, a 12-micrometre cubic
volume, a 12 × 12 × 12 magnetic-field grid, one thread, and a
100-microsecond acquisition. The MSE example uses a 10-microsecond time step;
the other examples use 1 microsecond. Lengths in the inputs are metres;
times are seconds.

| Input | Sequence and purpose | Free protons |
| --- | --- | --- |
| `fid.txt` | Free induction decay | 24 |
| `se.txt` | Spin echo with a pulse at 50 microseconds | 24 |
| `mse.txt` | Multi-spin echo with 20-microsecond echo spacing | 24 |
| `cell_importance.txt` | Single-cell FID with overlapping importance proposals | 128 |

These are quick execution examples. Their small grids, proton counts and run
times do not establish numerical convergence or validate a physical model.

## Build and run

From the repository root, use CMake 3.21 or newer, Ninja and a C++17 compiler
with OpenMP:

```bash
cmake --preset linux-gcc-release
cmake --build build/linux-gcc-release --target IONP
```

Run each example separately so its output is easy to identify:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/fid.txt \
  --output-dir data/fid --log-dir data/log/fid --seed 123

./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/se.txt \
  --output-dir data/se --log-dir data/log/se --seed 123

./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/mse.txt \
  --output-dir data/mse --log-dir data/log/mse --seed 123

./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/cell_importance.txt \
  --output-dir data/cell_importance --log-dir data/log/cell_importance \
  --save-region-signals --seed 123
```

For other build systems or Windows, use the executable path from the
[main build guide](../README.md#build). Repeating an example adds numbered files
to its output directory; use a fresh directory to compare runs. A fixed seed
and the same executable and settings reproduce these single-thread examples;
identical floating-point results across different toolchains are not promised.

## Read the outputs

Signal files contain magnitude followed by time in seconds. Particle files
contain the three `x y z` coordinates from `IONP_pos/positions.txt`. FID and SE
write one full signal curve. MSE writes both the full curve and a shorter
sampled curve: for this example, the sampled times are 0, 20, 40, 60 and
80 microseconds. Both MSE files use the same naming pattern, so distinguish
them by their row count. The full curve reaches 100 microseconds. The MSE
timestep matches the retained signal-record interval; see the
[MSE output limitation](../docs/simulations.md#signal-files) before using finer
timesteps in short MSE acquisitions.

The cell example places a cell of radius 1 micrometre at the origin. One
particle lies fully inside and two lie outside. Its proposals draw 64 protons
uniformly from accessible water, 32 inside the cell, and 32 outside the cell
within 3 micrometres of its centre. Proposal volumes and per-proton MIS weights
are calculated automatically. No manual weights are supplied.

With `--save-region-signals`, each active proposal writes a complex `mean` file
and a `mis_contribution` file, with columns `real imaginary time`. Sum the
three complex contributions and take the magnitude to reconstruct the main
signal. This example has no bound protons. The unweighted proposal means are
diagnostic outputs and cannot be added to reconstruct the signal.

To generate field-grid and trajectory data for the plotting tools, run:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/cell_importance.txt \
  --output-dir data/cell_diagnostics --log-dir data/log/cell_diagnostics \
  --record-positions --save-region-signals \
  --save-bgrid data/cell_diagnostics/bgrid.bin --seed 123
```

See the [plotting instructions](../README.md) for Python dependencies and the
[simulation reference](../docs/simulations.md) for configuration and output
details. Output filenames truncate radii, timesteps and echo times to integer
display units; use the input and signal time column for precise values.
For file-loaded geometry, signal and particle filenames retain the configured
`IONP_num = 0` label even though three particles are loaded. Read the actual
particle count from the output rows and logs.

## Validate the examples

The publication check uses only the Python standard library. Choose an absent
or empty work directory; the script refuses to reuse nonempty output folders.

```bash
python3 test/publication_examples.py \
  --executable build/linux-gcc-release/src/ionp/IONP \
  --source-dir . --work-dir build/publication-example-check
```

It runs all four examples twice with seed 123, checks completed simulations,
finite signals and expected times, particle geometry, field grids and
trajectories, reconstructs the cell signal from MIS contributions, and compares
the repeated numerical output files. Logs and outputs remain in the work
directory for inspection.
