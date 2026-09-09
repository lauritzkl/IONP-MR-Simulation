# Simulation Workflow And Reference

This document explains how the IONP MR simulator works, how to prepare
configuration files, how to run simulations, and how to interpret generated
data. It is written for research use first; code internals are included only
where they affect simulation behavior.

## Contents

- [Simulation pipeline](#simulation-pipeline)
- [Running simulations](#running-simulations) and [CLI reference](#cli-reference)
- [Config files](#config-files), [geometry files](#custom-geometry-files), and [importance sampling](#importance-sampling)
- [Synthetic examples](#example-configs)
- [Output files](#output-files) and [logs](#logs)
- [Reproducibility and limitations](#reproducibility-and-limitations)
- [Verification](#verification)

## Simulation Pipeline

Each simulation run follows this pipeline:

1. **Collect configs.** The executable reads either one `--config` file or all
   `.txt` files directly inside `--ini-dir`. Directory scanning is not
   recursive.
2. **Parse parameters.** `Data::readData()` parses `key = value` lines. Empty
   lines and lines starting with `#` are ignored. Keys are case-insensitive and
   ignore spaces and underscores, so `IONP_radius`, `IONP radius`, and
   `ionpradius` normalize to the same key if that alias exists.
3. **Create the simulation matrix.** The matrix length/origin/node settings
   define the physical simulation volume and optional B-field interpolation
   grid. If no volume length is given, the code derives a cubic volume from the
   particle count and the built-in volume fraction.
4. **Place IONPs, aggregates, and cells.** IONPs can be placed randomly from
   concentration/cell settings, placed in aggregates, loaded from position
   files, or placed in/out of cells. Historical (legacy) placement is the
   default. `--legacy-placement` remains accepted for existing commands and
   selects the same behavior; it is no longer needed.
5. **Configure importance sampling.** If enabled, free protons use either
   overlapping single-region proposals with per-proton MIS weights or the
   legacy exclusive cell/joint strata.
6. **Build the B-field grid.** If matrix node counts are nonzero, the simulator
   precomputes a B-field grid and interpolates it during proton propagation.
   Every grid node always includes every magnetic field source. With
   `Surrounding_voxels = Yes`, those sources are the realized central IONPs plus
   either exact translated copies or iron-conserving superparticle templates in
   the 26 immediately adjacent voxels. The retained `--legacy-bgrid` option is
   a compatibility no-op because the full source calculation is the default
   again.
7. **Propagate protons.** Free protons diffuse through the volume with periodic
   outside-volume wrapping. Cell membranes can be given a symmetric crossing
   probability; denied entry and exit attempts reflect specularly. Free protons
   currently inside cells use `Cell_D` when configured. IONP cores remain
   impermeable. The coating/shell can likewise be given a crossing probability;
   denied coating crossings reflect at the coating boundary. If a coating
   diffusion coefficient is configured, free protons currently in the shell use
   that coefficient for their diffusion step. Bound protons are sampled in
   positive-thickness IONP shells using either `Bound_Prot_n` or the
   population-specific bound controls.
8. **Apply MR sequence timing.** The simulator accumulates phase from the local
   B-field and optional T2 noise. SE and MSE runs apply refocusing at the
   configured echo spacing.
9. **Aggregate signal.** Standard runs average `exp(i phase)` over all free and
   bound proton samples. Overlapping-proposal runs apply an initial-position
   weight to each free proton; legacy importance sampling combines exclusive
   region means by normalized weights.
10. **Write outputs.** Signal, IONP, aggregate, optional cube signal, optional
    trajectory, optional B-field grid, and logs are written to the configured
    output directories.

## Running Simulations

Build the simulator as described in the [README](../README.md#build). From the
repository root, select a synthetic example with explicit paths:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/fid.txt \
  --output-dir data/fid --log-dir data/log/fid --seed 123
```

Use fresh output/log directories for each independent run. Add `--no-log` to
print diagnostics to the terminal. To run the entire synthetic set in a batch,
omit `--config` and use a separate output directory:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples \
  --output-dir data/example-batch --log-dir data/log/example-batch --seed 123
```

A batch shares a random stream, so its later results differ from individually
seeded runs. The scan includes only `.txt` files directly in the selected
directory and sorts their paths; it is not recursive.

The executable detects a repository root by searching ancestors of the current
directory and then the executable path for both `ini/` and `src/ionp/`. With no
arguments, it scans `<repo>/ini` and uses `<repo>/data/tmp` and `<repo>/data/log`.
The public release includes only a README in `ini/`, so there is no default
batch until you add your own configs. The curated examples live in `examples/`.

The retained generated-example option also works from the repository root:

```bash
./build/linux-gcc-release/src/ionp/IONP --example --seed 123
```

It writes `Testing/example/ini/example_config.txt` and outputs under
`Testing/example/data`. This generated config has meter-scale geometry and is
an execution demonstration, not a realistic research configuration.

## CLI Reference

| Option | Effect |
| --- | --- |
| `--ini-dir <dir>`, `-i <dir>` | Directory scanned for `.txt` files when `--config` is not set. Also used as the base for manual `IONP_pos` and `Agg_pos` loading. |
| `--ini-root <dir>` | Base directory for auxiliary input folders: `Cell_pos`, `Cell_load`, `Histogram`, and `Radius_Histogram`. Defaults to `<repo>/ini`; an explicit `--ini-dir` triggers the root inference described below. |
| `--config <file>` | Run exactly one config file, relative to the current working directory. Does not change auxiliary input roots. |
| `--output-dir <dir>` | Directory for signal, IONP, aggregate, cube-signal, and trajectory outputs. |
| `--log-dir <dir>` | Directory for redirected stdout logs. |
| `--no-log` | Do not redirect stdout to a log file. |
| `--seed <integer>` | Set the random seed for reproducible placement and diffusion. |
| `--example` | Create and run a small generated MSE example. |
| `--record-positions` | Write full proton trajectory samples and boundary-hit events. |
| `--save-region-signals` | In legacy importance sampling, write each active exclusive-region mean. In overlapping-proposal mode, write both the unweighted source-proposal mean and its additive MIS-weighted contribution. |
| `--legacy-bgrid` | Compatibility no-op retained for existing commands. Every grid point already includes every IONP by default. |
| `--legacy-placement` | Explicitly select the default historical placement; retained for command compatibility. |
| `--legacy-ionp-check` | Disable spatial-index IONP overlap checks and loop through all IONPs for proton/IONP checks. |
| `--save-bgrid <file>` | Save the calculated B-field grid to a binary file. |

When `--ini-root` is omitted and `--ini-dir` is supplied, the root is the
selected directory if it is named `ini` or contains a recognized auxiliary
folder. Otherwise its parent is used if that parent satisfies the same rule;
otherwise the selected directory itself is used. Recognized folder names are
`Cell_pos`, `Cell_load`, `Histogram`, `Radius_Histogram`, `IONP_pos`, and `Agg_pos`.
Using both options explicitly avoids dependence on the surrounding folders.

CLI options are processed in order. `--example` resets all input/output paths
when encountered, so any intended path overrides must follow it. Unknown
options, including `--help`, and options missing a value are currently ignored.
Use the documented commands; this release preserves that parser behavior.

## Config Files

Configuration files are plain text:

```text
key = value
# comments are allowed
```

Lengths are in meters; timing and T2 constants are in seconds. Diffusion
coefficients are in m²/s and crossing probabilities are dimensionless values
in [0, 1]. Bulk `Concentration` and `Concentration_extern` are interpreted as
µmol/L of iron; cell iron-load values and `Cell_conc` are in fg/cell. A nonzero
`Magnetic Moment` supplies a dipole moment in A·m² (equivalently J/T).

The radius-based field model uses compiled constants from
[include/mptmacros.h](../include/mptmacros.h): outside-water diffusion
`D_D_CONST = 3e-9 m²/s`, IONP magnetization `D_M_IONP = 3.8e5 A/m`, proton
gyromagnetic ratio `D_GI = 2.6752218744e8 rad/(s·T)`, `D_PI = 3.14159265`, and
permeability `D_MO = 4 * D_PI * 1e-7`. The equatorial field factor is
`B_fixed = D_MO * D_M_IONP / 3`. These are the implementation's retained model
values, not claims of newly verified physical constants. Changing them requires
rebuilding and can change the simulation results.

### Geometry And Grid

| Key | Aliases | Meaning |
| --- | --- | --- |
| `Matrix_x0`, `Matrix_y0`, `Matrix_z0` | | Matrix origin components. |
| `Origin` | | Sets all origin components to the same value. |
| `Matrix_xlength`, `Matrix_ylength`, `Matrix_zlength` | `xLength`, `yLength`, `zLength` | Simulation volume dimensions. |
| `Length` | | Sets all three dimensions to the same value. |
| `Matrix_xnodes`, `Matrix_ynodes`, `Matrix_znodes` | `xNodes`, `yNodes`, `zNodes` | B-field grid node counts. Use zero nodes to disable the interpolation grid. |
| `Nodes` | | Sets all three node counts to the same value. |
| `BField_solver` | | Proton-position field evaluator: `current` (default) or the opt-in `corrected_hybrid`. The latter requires at least two nodes per axis and stationary IONPs. |
| `Surrounding_voxels` | | `Yes` or `No`. When enabled, add translated magnetic sources in the 26 immediately adjacent voxels. These are exact copies unless superparticle compression is configured. Defaults to `No` when omitted. |
| `Surrounding_voxel_superparticles` | | Non-negative maximum number of magnetic-only particles used per nonempty cell or extracellular group in each surrounding voxel. `0` (the default) keeps exact copies. A positive value requires `Surrounding_voxels = Yes` and `Magnetic Moment = 0`. |
| `Voxel position` | `Voxel_pos` | Voxel position used by cell/placement logic. |
| `Cube Origin`, `Cube Size` | | Defines a cube region for extra cube-signal output. |

### Particle And Aggregate Settings

| Key | Aliases | Meaning |
| --- | --- | --- |
| `IONP_radius` | `IONP_r`, `r_ionp`, `radius_ionp` | Primary IONP radius. |
| `IONP_r1` | `IONP_radius1`, `r_ionp1` | Secondary radius for bimodal radius distributions. |
| `IONP_per` | `IONP_probability`, `IONP_distribution_probability` | Probability of drawing from the primary radius distribution when `IONP_r1` is set. Accepts `0.21` or legacy percent form `21`; defaults to `0.21`. |
| `IONP_r_std` | `IONP_radius_std`, `r_ionp_std` | Primary radius distribution width. |
| `IONP_r_std1` | `IONP_radius_std1`, `r_ionp_std1` | Secondary radius distribution width. |
| `IONP_num` | `IONP_number`, `IONP_n` | Requested number of IONPs. |
| `IONP_position_file` | | Optional leaf filename selecting one custom IONP geometry from `<ini-dir>/IONP_pos`. |
| `Agg_radius` | `Agg_r`, `r_agg` | Aggregate radius. |
| `Agg_r_std` | `Agg_radius_std`, `r_agg_std` | Aggregate radius distribution width. |
| `Agg_num` | `Agg_number`, `Agg_n` | Number of aggregates. |
| `Agg_position_file` | | Optional leaf filename selecting one custom aggregate geometry from `<ini-dir>/Agg_pos`. |
| `Coating` | `shell`, `coating_thickness` | Intended total outer coating radius. Each realized IONP stores `max(0, Coating - core radius)` as its shell thickness, so a core larger than `Coating` has no shell and remains its own outer radius. |
| `Coating_D` | `D_coating`, `coating_diffusion`, `coating_diffusion_coefficient` | Diffusion coefficient used by free protons while inside the IONP coating/shell, in m^2/s. If omitted, the outside-water value `D_D_CONST` is used. |
| `Coating_permeability` | `Coating_P`, `coating_probability`, `coating_permeation_probability` | Probability in `[0, 1]` that a free proton crossing the outer coating boundary is accepted. If omitted, crossings are always accepted and no coating-boundary check is performed. |
| `Concentration` | `conc` | IONP concentration inside cells or globally, depending on cell settings. |
| `Concentration_extern` | `conc_ext` | IONP concentration outside cells. |
| `Magnetic Moment` | | Magnetic moment override. If zero, the radius-based magnetization model is used. |
| `Radius_hist` | `Radius_histogram` | Use radius histogram data from `<ini-root>/Radius_Histogram`. |

If `IONP_num = 0`, the simulator loads custom IONP positions from
`<ini-dir>/IONP_pos`. If `IONP_num > 0`, `Agg_num = 0`, and `Agg_radius > 0`,
it loads custom aggregate positions from `<ini-dir>/Agg_pos`. The optional
filename selectors and compatibility behavior are detailed under
[Custom Geometry Files](#custom-geometry-files).

When `IONP_num > 0`, `Agg_num > 0`, `Agg_radius > 0`, and `Conc_ext > 0`,
the configured aggregate count and radius apply only to the extracellular
population. Normal intracellular placement runs first when cells are present.
Aggregate centers are then placed inside the voxel, outside all cells, and
without aggregate overlap. `Conc_ext` determines the total iron mass distributed
equally across the aggregates; each aggregate is filled with the configured
IONP radius distribution and coating. IONP centers are sampled uniformly by
volume throughout the available aggregate sphere and rejected and resampled if
they overlap another IONP, matching the legacy intracellular-aggregate
placement method. Very dense coated-particle configurations can exhaust the
random-placement attempts; increase `Agg_radius` or `Agg_num` in that case.
Legacy intracellular aggregate configurations with `Conc_ext = 0` retain their
existing behavior.

### Cell Settings

| Key | Aliases | Meaning |
| --- | --- | --- |
| `Cell Origin` | `Cell_orig` | Cell center for single-cell cases or initial reference point. |
| `Cell Radius` | `Cell_rad` | Spherical cell radius. |
| `Cell Number` | `Cell_num` | Number of cells to place. |
| `Cell_permeability` | `Cell_P`, `cell_probability`, `cell_permeation_probability`, `Cellular_permeability` | Probability in `[0, 1]` that a free proton crossing a cell membrane is accepted, for both entry and exit. Denied crossings reflect specularly. The default is `0`, making cells impermeable. Percent values are not accepted. |
| `Cell_D` | `D_cell`, `cell_diffusion`, `cell_diffusion_coefficient` | Diffusion coefficient used by free protons while inside cells, in m^2/s. If omitted, the outside-water value `D_D_CONST` is used. Coating diffusion still takes precedence while inside an IONP coating. |
| `Cell Concentration` | `Cell_conc` | Median cell iron load. On the non-aggregate intracellular placement path, a nonzero value represents each cell load with `IONP_num` equal iron-equivalent particles; their derived core radius overrides the configured radius distribution inside cells only. Extracellular particles still use the configured radius distribution. |
| `Cell Concentration Std` | `Cell_conc_std` | Cell concentration spread. |
| `Cell Concentration Distribution` | `Cell_conc_dist` | `Normal` or `Lognormal`. |
| `Cell_hist` | `Cell_histogram` | Use cell iron-load histogram data from `<ini-root>/Histogram`. |
| `Cell_position_file` | | Optional leaf filename selecting one cell-center file from `<ini-root>/Cell_pos`. |
| `Cell_load_file` | | Optional leaf filename selecting one cell-load file from `<ini-root>/Cell_load` on the non-aggregate path. |

When `Cell Number = 0` and `Cell Radius > 0`, some placement paths load cell
positions from `<ini-root>/Cell_pos` and cell loads from
`<ini-root>/Cell_load`.

### Custom Geometry Files

The four optional selectors let configs scanned from the same `--ini-dir`
choose different custom geometries:

| Canonical key | Numeric activation rule | File location |
| --- | --- | --- |
| `IONP_position_file` | `IONP_num = 0` | `<ini-dir>/IONP_pos/<filename>` |
| `Agg_position_file` | `IONP_num > 0`, `Agg_num = 0`, and `Agg_radius > 0` | `<ini-dir>/Agg_pos/<filename>` |
| `Cell_position_file` | `Cell_num = 0` and `Cell_radius > 0`, when the selected placement path supports file-based cells | `<ini-root>/Cell_pos/<filename>` |
| `Cell_load_file` | The existing non-aggregate file-based cell-load path | `<ini-root>/Cell_load/<filename>` |

Key parsing remains case-insensitive and ignores spaces and underscores. For
example, `Cell_position_file`, `cell position file`, and
`CELLPOSITIONFILE` identify the same setting.

Selector values have these rules:

- Use only a leaf filename, such as `cells cohort 3.txt`. Absolute paths and
  names containing directory components are rejected.
- The extension must be `.txt` or `.dat`, and the resolved target must be a
  regular file. Parsing it must produce at least one input record.
- An active selector loads exactly the named file. It does not combine that
  file with other files in the folder.
- If an active selector is omitted, compatibility behavior scans and
  concatenates every regular `.txt` and `.dat` file in the corresponding
  folder.
- A selector outside its activation branch is not validated or loaded. The log
  warns that it was ignored. Active loads log whether an explicit file or the
  legacy folder scan was used, including resolved paths.

Custom IONP placement has the same precedence as before: `IONP_num = 0` enters
the manually positioned IONP branch before cell or aggregate placement is
considered. A custom aggregate selector is used only with a positive requested
IONP count, zero requested aggregates, and a positive aggregate radius.

For non-aggregate file-based cells with `Cell_hist = Yes`,
`Cell_position_file` and `Cell_load_file` are a pair. Set both to select one
matched pair of files, or set neither to retain the legacy scans of both
`Cell_pos` and `Cell_load`; setting only one is an error. After loading, the
number of cell loads must equal the number of cell centers.

Aggregate placement may use `Cell_position_file` for its cell geometry, but it
does not consume per-cell data from `Cell_load`. Its existing concentration and
`<ini-root>/Histogram` behavior remains unchanged, so a configured
`Cell_load_file` is ignored with a warning on that path.

Select one IONP file for a manual-geometry simulation:

```ini
Length = 1.0E-4
IONP_num = 0
IONP_position_file = ionps batch 07.dat
```

Select one custom aggregate file:

```ini
IONP_num = 250
Agg_num = 0
Agg_radius = 2.5E-7
Agg_position_file = aggregates run 12.txt
```

Select matching cell-center and load files on the non-aggregate path:

```ini
IONP_num = 100
Agg_num = 0
Agg_radius = 0
Cell_num = 0
Cell_radius = 1.0E-5
Cell_hist = Yes
Cell_position_file = cells cohort 3.txt
Cell_load_file = loads cohort 3.dat
```

`Cell_permeability` is a dimensionless Monte Carlo probability evaluated once
per membrane-crossing attempt, not a physical permeability coefficient with
units. Consequently, results for values strictly between `0` and `1` can depend
on the diffusion time step: changing the step size can change how often a
trajectory attempts to cross a membrane.

For reduced intracellular free-water diffusion, add for example:

```ini
Cell_D = 1.5E-9
```

### Auxiliary Data Formats

Geometry and histogram files contain whitespace-delimited numbers without
headers or comments. Put each complete record on its own line:

| File type | Record | Units |
| --- | --- | --- |
| IONP, aggregate or cell positions | `x y z` | meters |
| Cell loads | `iron_load` | fg/cell, in the same order as cell centers |
| Cell-load histogram | `bin_value frequency` | fg/cell and integer count |
| Radius histogram | `bin_value frequency` | meters and integer count |

Frequencies are read numerically and converted to integers. Use nonnegative
integer frequencies and at least one positive bin. These readers stop at the
first nonnumeric or incomplete record and may retain a partial dataset without
an error; verify the loaded record counts in the logs. Config numeric parsing
also accepts some numeric prefixes rather than validating the entire token.
Supply finite numeric values without units appended to the value.

Keep exactly one intended histogram file in each active `Histogram` or
`Radius_Histogram` folder. Existing loaders visit files in filesystem order and
replace the histogram sample vector for each file. Multiple files therefore
make selection order-dependent; geometry selectors do not select histograms.

### MR Sequence And Runtime

| Key | Aliases | Meaning |
| --- | --- | --- |
| `Proton_num` | `Prot_num`, `Proton_number`, `Prot_number` | Number of free proton samples in standard mode. |
| `Bound_Prot_n` | `Bound_Proton_number` | Bound proton samples per IONP. |
| `T_E` | `echo_time`, `TE` | Echo time / simulation duration. |
| `time_step` | `dt`, `delta_t` | Base diffusion time step. |
| `variable_Time` | `variable_timestep` | Explicitly enable or disable variable diffusion substeps. Accepts `Yes`, `No`, `true`, `false`, `ja`, `nein`. |
| `echo_spacing` | `SE_time` | Echo spacing for SE/MSE refocusing. |
| `Sequence` | `MRSequence` | `FID`, `SE` / `SpinEcho`, or `MSE` / `MultiSpinEcho`. |
| `Threads` | `num_threads` | OpenMP thread count. Use `0` for max available threads. |
| `T2out` | `background_T2` | Background/outside-cell T2 relaxation term. |
| `T2in` | `cell_T2` | Inside-cell T2 relaxation term. |
| `v_IONP` | `velocity_IONP` | IONP velocity for moving-particle simulations. |
| `Direction` | `move_direction` | IONP movement direction vector. |

If `variable_Time` is omitted, variable stepping is enabled only when
`time_step = 0`. If `variable_Time = Yes`, a positive supplied `time_step`
remains the base step and can be subdivided near IONPs; `time_step = 0` uses
the automatic radius/coating-derived base step.

With the default `BField_solver = current`, a variable substep evaluates the
field directly from every IONP; normal steps use the grid wherever the existing
region logic permits it. With `BField_solver = corrected_hybrid`, every normal
and variable step instead uses the same corrected evaluator, so crossing the
variable-step threshold does not switch field models.

The corrected evaluator starts with an interpolation of the complete grid. It
then uses a sparse particle index to find only nearby IONPs and replaces each
near particle's interpolated contribution with its exact dipole field. For grid
spacing `Delta = max(dx, dy, dz)`, particle `i` is fully corrected through
`R0 = max(8 * outer_radius_i, 3 * Delta)`. A cubic smoothstep tapers the
correction to zero at `R1 = R0 + Delta`. Distant IONPs remain represented by
the complete grid and are not double-counted. The magnetic amplitude still
uses the core radius; the outer radius controls only the spatial correction
range.

Enable it in a selected simulation with:

```ini
BField_solver = corrected_hybrid
```

#### Surrounding magnetic voxels

Enable one shell of neighboring magnetic sources with:

```ini
Surrounding_voxels = Yes
```

The simulator forms a `3 x 3 x 3` magnetic environment from the central voxel
and its 26 immediate neighbors. By default, each neighbor receives an exact
translated copy of the complete realized central IONP collection; it does not
draw an independent random distribution. Core radii, coating thicknesses, and
relative particle positions are therefore identical in every exact copy. For
`N` physical IONPs, the exact mode presents `27 * N` magnetic sources.

Large surrounding collections can instead be compressed into iron-conserving
superparticles:

```ini
Surrounding_voxels = Yes
Surrounding_voxel_superparticles = 10
```

The positive count is a maximum per nonempty group, not a global maximum. The
realized central IONPs are partitioned into one group for every containing cell
and one extracellular group. A group containing at most the configured count
is copied exactly. A larger group is represented by that many equal, uncoated,
magnetic-only particles. For a group with physical core radii `r_i` and proxy
count `K`, every proxy has core radius
`cbrt(sum(r_i^3) / K)`. Thus each cell's iron and the extracellular iron are
conserved separately, including any overshoot produced while realizing a
concentration target.

Proxy centers are drawn uniformly and without overlap inside their source cell,
or inside the voxel and outside every cell for the extracellular group.
Particles that straddle a cell membrane or belong ambiguously to multiple cells
are rejected rather than assigned to the wrong iron pool. If the requested
proxies cannot fit, the simulation fails with a diagnostic recommending a
larger count. One random template is built and translated into all 26
neighbors, so neighbors do not receive independent draws. Supplying `--seed`
makes that template reproducible across runs. If the resulting template has
`T` sources, the field calculation sees
`N + 26 * T` sources; a one-cell-only case with more than 10 physical IONPs and
a cap of 10 therefore adds only 260 surrounding sources.

The central `N` particles always remain fully explicit and are the sole input
to physical collision and coating-boundary handling, cell and population
classification, importance sampling, bound-proton generation, concentration
reporting, and IONP output. Protons also remain confined to the central
simulation coordinates by the existing periodic wrap: a proton that crosses
one voxel face re-enters through the opposite face rather than entering a
neighboring voxel.

All magnetic evaluation paths use the expanded source collection: B-field grid
construction, direct field evaluation (including the variable-step fallback),
corrected-hybrid near-source corrections, and the magnetic proximity checks
that select direct evaluation or smaller time steps. In moving direct-field
simulations, surrounding magnetic sources participate in the same movement
updates.

Omitting the key, or setting `Surrounding_voxels = No`, retains the previous
central-only behavior and cost. Omitting
`Surrounding_voxel_superparticles`, or setting it to `0`, retains exact
surrounding replicas. A positive count requires `Surrounding_voxels = Yes` and
`Magnetic Moment = 0`, because proxy magnetization is derived from its
iron-conserving radius. The original physical collection and a magnetic spatial
index are retained separately, so memory and run-time growth depend on the
template size, selected evaluator, and rest of the simulation. All three voxel
dimensions must be positive when surrounding voxels are enabled.

This option is a finite one-voxel shell, not an infinite periodic or Ewald
dipole summation. Truncating the environment to 26 neighbors therefore does not
guarantee that the calculated B-field is exactly equal on opposing voxel faces.

### Importance Sampling

Importance sampling reduces the number of simulated free protons by placing
more starting points in regions of interest and correcting their contributions
to the signal. The implementation represents the two choices explicitly as
`LegacyJoint` and `OverlappingProposal` modes:

- **Overlapping proposals (recommended):** configure simple marginal regions,
  such as `cell_inside` or `intracellular_ionp_near`. These proposals may
  overlap; per-proton multiple-importance-sampling (MIS) weights prevent double
  counting.
- **Legacy exclusive strata:** retain the original cell-only or 27-way joint
  region configuration and fixed region weights.

#### Overlapping proposals

The ten available proposal-count keys are:

| Key | Starting positions drawn uniformly from |
| --- | --- |
| `IS_protons_uniform` | All accessible free-water volume. This count must be positive. |
| `IS_protons_cell_inside` | Inside any cell. |
| `IS_protons_cell_near` | Outside cells but within `IS_cell_near_cutoff` of the nearest cell center. |
| `IS_protons_cell_far` | The complementary far-cell region. |
| `IS_protons_intracellular_ionp_near` | Near an IONP located completely inside a cell. |
| `IS_protons_intracellular_ionp_intermediate` | At intermediate distance from an intracellular IONP. |
| `IS_protons_intracellular_ionp_far` | Far from every intracellular IONP. |
| `IS_protons_extracellular_ionp_near` | Near an IONP outside all cells. |
| `IS_protons_extracellular_ionp_intermediate` | At intermediate distance from an extracellular IONP. |
| `IS_protons_extracellular_ionp_far` | Far from every extracellular IONP. |

Each value is a requested number of trajectories, not a region weight or a
fraction of the total. The free-proton count is the sum of all ten values and
replaces `Proton_num`. Omitted proposal keys have count zero. At least one of
the new keys selects overlapping-proposal mode; that mode cannot be mixed with
legacy `IS_protons_inside/near/far`, joint free-region keys, or any manual free
weight key. Bound-proton weights remain allowed.

For example:

```text
ImportanceSampling = Yes
IS_cell_near_cutoff = 0.000012
IS_intracellular_ionp_near_cutoff = 0.000001
IS_intracellular_ionp_far_cutoff = 0.000004
IS_extracellular_ionp_near_cutoff = 0.000001
IS_extracellular_ionp_far_cutoff = 0.000004

IS_protons_uniform = 3000
IS_protons_cell_inside = 2000
IS_protons_intracellular_ionp_near = 2500
IS_protons_extracellular_ionp_near = 2500
```

This requests 10,000 free trajectories. A starting point may simultaneously be
inside a cell, near an intracellular IONP, and near an extracellular IONP. No
special cross-product count is required.

IONP distances are measured from the nearest realized outer surface:

```text
surface distance = center distance - (core radius + stored shell thickness)
```

Near covers `0 <= distance <= X`, intermediate covers `X < distance <= Y`,
and far covers `distance > Y`. Any positive proposal count for an IONP
population, including its far proposal, requires that population's cutoffs to
satisfy `0 < X < Y`. Free protons are never initialized at negative surface
distance. `cell_inside` and `cell_near` require a positive `Cell Radius`, and
`cell_near` requires `IS_cell_near_cutoff > Cell Radius`.

The positive uniform proposal guarantees that every accessible free-water
starting position remains represented. If a requested nonuniform proposal has
zero estimated support in the realized random geometry, or cannot be
initialized within the retry limit, its complete count is moved to the uniform
proposal with a warning. For example, near/intermediate proposals have no
support when their IONP population is absent. A far proposal remains valid in
that case because every free-water point is far from the absent population. The
total free-proton count does not change, and a missing random overlap can no
longer prevent startup.

##### Exact overlap weights

Let `D` be the accessible free-water volume. For proposal `k`, let `A_k`
be its region, `V_k` its accessible volume, `n_k` its realized count after
any redistribution, and `N = sum_k n_k`. Its density and the deterministic
mixture density are

```text
q_k(x)   = 1[x is in A_k] / V_k
q_mix(x) = sum_k (n_k / N) q_k(x).
```

For every simulated free proton starting at `x_i`, the simulator first forms

```text
h_i = 1 / sum_k (n_k * 1[x_i is in A_k] / V_k)
```

and then normalizes the free-proton coefficients:

```text
w_i = h_i / sum_j h_j.
```

The free signal is therefore

```text
S_free(t) = sum_i w_i * exp(i * phase_i(t)).
```

For example, suppose the active counts are uniform `n_U`, intracellular-near
`n_I`, and cell-inside `n_C`, with volumes `V_D`, `V_I`, and `V_C`.
A proton in the intracellular-near/cell-inside overlap has

```text
h = 1 / (n_U / V_D + n_I / V_I + n_C / V_C),
```

whereas an intracellular-near proton outside cells has

```text
h = 1 / (n_U / V_D + n_I / V_I).
```

The weight depends on every active proposal containing the initial point, not
on which proposal happened to generate it. Region membership and the MIS
weight are fixed from the initial position; later diffusion across region
boundaries does not require reclassification. The required volumes are
estimated with a deterministic-seed 200,000-point pilot before trajectory
propagation. The displayed MIS formula is exact conditional on those estimated
fractions; the pilot ratios and final self-normalization retain a small
finite-sample ratio-estimator bias.

When bound samples are present, let `f_D = V_D / V_box` and let `b_in` and
`b_out` be the configured intracellular/extracellular bound weights. The final
complex signal before taking its magnitude is

```text
S(t) = (f_D * S_free(t) + b_in * S_bound_in(t) + b_out * S_bound_out(t))
       / (f_D + b_in + b_out).
```

With the common `IS_weight_bound` interface, that weight is split between the
two bound populations in proportion to their realized sample counts, matching
legacy behavior.

##### Computational load and timing

MIS weights are evaluated during initialization and aggregation rather than
inside the diffusion/B-field loop. The 200,000-point volume pilot adds a fixed
startup cost, which can dominate short simulations. Initialization cost depends
on geometry and rejection frequency. Aggregation adds a weight multiplication
for each free proton and recorded time point.

Each free proton stores a `double` coefficient and a proposal label, plus vector
bookkeeping. Logs report volume-pilot, proposal-initialization, propagation and
aggregation timings so comparisons can be measured on the actual configuration
and hardware. No fixed percentage overhead is guaranteed.

Runtime depends on the total requested samples and their trajectories. Adding
uniform samples without reducing another proposal increases the work. Starts
near IONPs can spend more time in the 64-way variable-substep path; overlap with
both IONP populations does not apply that subdivision twice.

#### Bound sampling in either mode

The bound controls are shared by `LegacyJoint` and `OverlappingProposal`:

| Key | Meaning |
| --- | --- |
| `IS_weight_bound` | Signal weight for bound protons. Required when `Bound_Prot_n > 0`. |
| `IS_bound_protons_intracellular_ionp` | Bound coating samples per coated intracellular IONP. The old `..._ionp_inside` name remains accepted. |
| `IS_bound_protons_extracellular_ionp` | Bound coating samples per coated extracellular IONP. The old `..._ionp_outside` name remains accepted. |
| `IS_weight_bound_intracellular_ionp` | Signal weight for intracellular-IONP bound samples. The old `..._ionp_inside` name remains accepted. |
| `IS_weight_bound_extracellular_ionp` | Signal weight for extracellular-IONP bound samples. The old `..._ionp_outside` name remains accepted. |

#### Legacy exclusive strata

Legacy particle-aware configurations classify each starting point along three
axes:

```text
cell location × distance to intracellular IONPs × distance to extracellular IONPs
```

Each axis has three states, producing up to 27 exclusive free-proton strata. A
point that is simultaneously near a cell, an intracellular IONP, and an
extracellular IONP belongs to one joint tuple and contributes to the signal
once.

The words have one fixed meaning throughout the C++ code, INI keys, logs, and
output filenames:

| Axis | States | What it describes |
| --- | --- | --- |
| `cell` | `inside`, `near`, `far` | The proton's location relative to the nearest cell. |
| `intracellular_ionp` | `near`, `intermediate`, `far` | The proton's surface distance to the nearest IONP located completely inside a cell. |
| `extracellular_ionp` | `near`, `intermediate`, `far` | The proton's surface distance to the nearest IONP outside all cells. |

Double underscores separate the axes in a joint region name. For example,
`cell_inside__intracellular_ionp_near__extracellular_ionp_far` means that the
proton is inside a cell, near an intracellular IONP, and far from every
extracellular IONP.

| Key | Meaning |
| --- | --- |
| `ImportanceSampling` | `Yes` or `No`. Enables importance sampling. |
| `IS_cell_near_cutoff` | Distance from cell center separating near and far cell regions. The legacy name `IS_near_cutoff` remains accepted. |
| `IS_intracellular_ionp_near_cutoff` | X surface distance for IONPs whose complete outer sphere is inside a cell. The legacy name `IS_ionp_inside_near_cutoff` remains accepted. |
| `IS_intracellular_ionp_far_cutoff` | Y surface distance for intracellular IONPs. The legacy name `IS_ionp_inside_far_cutoff` remains accepted. |
| `IS_extracellular_ionp_near_cutoff` | X surface distance for IONPs outside all cells. The legacy name `IS_ionp_outside_near_cutoff` remains accepted. |
| `IS_extracellular_ionp_far_cutoff` | Y surface distance for extracellular IONPs. The legacy name `IS_ionp_outside_far_cutoff` remains accepted. |
| `IS_protons_inside` | Number of free protons sampled inside cells. |
| `IS_protons_near` | Number of free protons sampled outside the cell but inside `IS_near_cutoff`. |
| `IS_protons_far` | Number of free protons sampled outside the near region. |
| `IS_weight_inside` | Signal weight for the inside-cell region. |
| `IS_weight_near` | Signal weight for the near-cell region. |
| `IS_weight_far` | Signal weight for the far-cell region. |

Joint free counts and optional explicit weights use named keys:

```text
IS_protons_cell_<C>__intracellular_ionp_<I>__extracellular_ionp_<E>
IS_weight_cell_<C>__intracellular_ionp_<I>__extracellular_ionp_<E>
```

`<C>` is `inside`, `near`, or `far`; `<I>` and `<E>` are `near`,
`intermediate`, or `far`. For example:

```text
IS_protons_cell_inside__intracellular_ionp_near__extracellular_ionp_far = 500
IS_weight_cell_inside__intracellular_ionp_near__extracellular_ionp_far = 0.015
```

Existing joint keys using `ionp_inside`/`ionp_outside` and `middle` remain valid,
so old simulation configs do not need to be rewritten.

IONP distances are measured from the nearest realized surface:

```text
surface distance = center distance - (core radius + stored shell thickness)
```

Near covers `0 <= distance <= X`, intermediate covers `X < distance <= Y`,
and far covers `distance > Y`. Free protons are never initialized at negative
surface distance, including for an IONP population whose distance axis is not
actively stratified.

Validation rules:

- Inside/near-cell strata require a positive `Cell Radius`.
- `IS_cell_near_cutoff` may be omitted when no near-cell stratum is sampled.
  When configured, or when a near-cell stratum is requested, it must be larger
  than `Cell Radius`.
- Each configured IONP cutoff pair must satisfy `0 < X < Y`.
- A requested non-far intracellular/extracellular IONP stratum requires at least one
  realized particle in that population. An absent population fixes its axis to
  `far`.
- An IONP is intracellular only when its complete realized outer sphere fits inside
  a cell. Particles intersecting a membrane are rejected as ambiguous.
- Region proton counts must be non-negative.
- At least one free proton sample is required.
- Weights must be non-negative.
- A positive weight requires samples in the same joint stratum.
- Population-specific bound samples require coated particles and positive
  population-specific bound weights.
- Legacy and joint free keys cannot be mixed. Common and population-specific
  bound controls likewise cannot be mixed.

If the free-region weights are omitted or left at zero, the simulator estimates
joint volume weights automatically for the active free regions. The estimator
mixes uniform-volume proposals with shell-centered proposals and accounts for
points covered by multiple shells/populations in the proposal density.
Joint configurations with explicit weights use the same estimator as a
pre-propagation support check, so requested zero-volume tuples fail before any
trajectory is advanced.

```text
Automatic active joint importance-sampling weights:
  cell_inside__intracellular_ionp_near__extracellular_ionp_far: 0.015
```

In importance-sampling mode, `Proton_num` is replaced by the sum of active joint
free samples. Output names include those samples plus active bound samples.
All free strata may dynamically use variable substeps near IONPs; bound samples
remain on the fixed-step path.

## Example Configs

The [examples guide](../examples/README.md) provides four complete standalone
configs: [FID](../examples/fid.txt), [SE](../examples/se.txt),
[MSE](../examples/mse.txt), and [cell importance sampling](../examples/cell_importance.txt).
All use the included three-particle synthetic geometry, a single thread, and
short acquisitions. The cell example uses overlapping uniform, inside-cell and
near-cell proposals with automatically calculated MIS weights. No study inputs
or histograms are required.

These inputs demonstrate execution and file formats. Increase sampling and
check convergence before interpreting a run scientifically.

## Output Files

Output file names start with an incrementing file index based on the current
number of regular files in the output directory. A single simulation can write
multiple consecutive files. The log file name includes the next log index and
the first output index seen when the run starts.

### Signal Files

Signal files are named like:

```text
164_signal_FID_184nm_10IONPs_60Protons_5000ns_dt_100ms_te.dat
```

Format:

```text
signal_magnitude    time_seconds
```

Example:

```text
1                   0
0.99945368653399    1e-05
0.999205270917944   2e-05
```

For MSE runs, the simulator writes the full sampled signal and then a second
MSE signal selected from that full curve. Both use the `signal_MSE` naming
pattern, so use log order, prefix order, or row count to distinguish them.
The sampled curve retains time zero and excludes the final nominal echo.
The supplied MSE example produces five sampled rows at 0, 20, 40, 60 and
80 microseconds, while its full curve also includes 100 microseconds.

**Retained MSE indexing limitation:** the full curve normally records every
10 microseconds when the timestep is smaller. For acquisitions shorter than
0.1 seconds, the MSE selector still calculates indices using the simulation
timestep. This can select wrong times or omit rows. For example, `T_E = 1e-4`,
`time_step = 1e-6`, and `echo_spacing = 2e-5` produce only the zero-time sampled
row. The curated MSE example uses `time_step = 1e-5` to align those indices.
Always inspect the actual time column; the full signal remains available.
This publication cleanup does not change the output-selection algorithm.

With `--save-region-signals`, legacy exclusive sampling writes one unweighted
complex mean per active region, for example:

```text
165_signal_complex_FID_region_cell_near__intracellular_ionp_far__extracellular_ionp_far_184nm_10IONPs_60Protons_5000ns_dt_100ms_te.dat
```

Their three-column format is:

```text
signal_real    signal_imaginary    time_seconds
```

The normal signal file remains the magnitude of the weighted complex
combination of all active regions. Therefore, applying the normalized legacy
region weights to these complex values, summing them, and then taking the
magnitude reconstructs the combined signal. Legacy cell-only importance
sampling uses the same output mechanism; its filenames show the selected cell
state and `far` for both unstratified IONP axes.

Overlapping-proposal mode instead writes two complex files for each active
source proposal:

```text
*_proposal_<proposal-name>_mean_*.dat
*_proposal_<proposal-name>_mis_contribution_*.dat
```

The `mean` file is the unweighted mean of trajectories drawn from that proposal
and is intended for diagnostics. Because proposals overlap, these means cannot
be combined with one fixed weight per file. The `mis_contribution` file already
contains that proposal's additive per-proton MIS-weighted contribution; summing
all such complex contribution files and then taking the magnitude reconstructs
the free-proton part of the combined signal.

### IONP And Aggregate Files

IONP files are named like:

```text
165_IONPs_FID_184nm_10IONPs_60Protons_5000ns_dt_100ms_te.dat
```

Format:

```text
x    y    z
```

Aggregate files use the same coordinate format and have `Aggs` in the file
name.

### Cube Signal Files

If `Cube Size` is nonzero, standard runs can write `SignalCube` files for the
signal contribution from protons inside that cube. Cube signal output is
disabled in importance-sampling mode because the weighted sampling scheme does
not represent uniform cube occupancy directly.

### Proton Trajectory Files

When `--record-positions` is used, trajectory files are named like:

```text
87_ProtonTrajectories_FID_1nm_1IONPs_1Protons_2000000ns_dt_10ms_te.dat
```

Format:

```text
# time[s]    proton_index    kind    x    y    z
```

`kind` values:

- `0`: regular trajectory sample.
- `1`: boundary-hit point for a denied crossing/specular reflection. This event
  kind identifies the reflection event generically rather than a specific type
  of boundary.
- `2`: post-reflection position, supported by plotting tools if emitted.

Plot a trajectory:

```bash
python tools/plot_trajectories.py \
  --file data/tmp/87_ProtonTrajectories_FID_1nm_1IONPs_1Protons_2000000ns_dt_10ms_te.dat \
  --radius 1e-5 \
  --center 0 0 0 \
  --proton 0
```

### B-Field Grid Files

Use `--save-bgrid <file>` to write the B-field grid:

```bash
./build/linux-gcc-release/src/ionp/IONP \
  --ini-dir examples --ini-root examples --config examples/fid.txt \
  --output-dir data/fid-grid --log-dir data/log/fid-grid --seed 123 \
  --save-bgrid data/fid-grid/bgrid.bin
```

The binary format is:

1. `uint32 nx`
2. `uint32 ny`
3. `uint32 nz`
4. `nx * ny * nz` contiguous `double` values

Values use index `z * nx * ny + x * ny + y` (z-major order). Integers and
doubles are written in the host's native encoding; there is no byte-order,
geometry or units metadata. The plotting helper assumes a little-endian
32-bit dimension header and reads native 64-bit floating-point values, as on
the tested x86-64 platform. Preserve the config alongside the binary grid.
Reusing `--save-bgrid` across runs overwrites that file.

With the [optional plotting dependencies](../README.md#optional-plotting), plot a grid
or compare two previously saved grids:

```bash
python tools/plot_bgrid.py data/fid-grid/bgrid.bin --output data/fid-grid/bgrid.png
python tools/plot_bgrid.py data/tmp/bgrid_old.bin data/tmp/bgrid_new.bin --output data/tmp/bgrid_compare.png
```

## Logs

Logs are written to `data/log` unless `--no-log` is used. A log records:

- Start and finish time.
- Config path.
- Parsed IONP/cell/sequence/runtime settings.
- Realized core-radius, stored shell-thickness, and outer-radius ranges,
  including zero-shell and internal/external particle counts.
- Derived matrix dimensions if no dimensions were supplied.
- Effective importance-sampling proton count.
- Auto-estimated importance-sampling weights.
- Step count and real time step.
- Whether variable stepping and/or B-field grid interpolation are active.
- Whether surrounding magnetic voxels and superparticle compression are
  enabled, together with the per-group cap, template/replica/total source
  counts, and iron-volume conservation.
- OpenMP thread count.
- Output file paths.
- Elapsed runtime.

The synthetic cell-importance example reports 128 requested free samples,
three active proposals, and IONP populations of one intracellular and two
extracellular particles. Check the `finished computation` line and inspect
stderr as well as the output files. Timings depend on hardware, build, inputs
and thread count; record those details when reporting performance.

## Compatibility Notes

The following retained features matter when comparing configurations or results
from different development versions:

- **Explicit variable-time control.** `variable_Time` can now force variable
  stepping on or off independently from `time_step = 0`.
- **Importance sampling.** The recommended mode uses overlapping uniform,
  cell, and population-specific IONP proposals with automatic per-proton MIS
  weights. Legacy exclusive cell and joint configurations remain supported.
- **Effective proton naming.** Output file names use the sum of all requested
  free-proposal (or legacy free-region) and active bound sample counts.
- **MSE aggregation cleanup.** MSE echo-sampled output is derived from the
  aggregated full signal, reducing duplicated aggregation logic.
- **Trajectory recording.** `--record-positions` writes regular samples and
  reflection events for debugging diffusion/cell-boundary behavior.
- **Full B-field grid restored.** Every grid node is calculated from every
  IONP. `--legacy-bgrid` remains accepted but no longer changes the result.
- **Opt-in corrected hybrid.** `BField_solver = corrected_hybrid` combines the
  complete grid with sparse exact near-source corrections at proton positions.
- **Opt-in surrounding magnetic voxels.** `Surrounding_voxels = Yes` adds
  translated magnetic sources in the 26 adjacent voxels to every magnetic
  field evaluation while leaving physical-particle behavior in the central
  voxel unchanged. `Surrounding_voxel_superparticles` can cap each cell and
  extracellular replica group while conserving its iron; exact replication
  remains the default for backward compatibility.
- **B-field grid export.** `--save-bgrid` writes the calculated grid for
  plotting or comparison.

## Reproducibility And Limitations

Run one config per process with an explicit seed, thread count and separate
empty output/log directories. Save the complete config and auxiliary files,
source revision or publication manifest, compiler and standard-library versions,
CMake options, operating system and OpenMP settings. The process seeds its RNG
once before a sorted batch; earlier configs consume random draws. Identical
results across compilers, standard libraries or platforms are not guaranteed.

Individual simulation exceptions print `Simulation failed for` to stderr but
the program continues and can exit zero. Check stderr, the completion log and
expected nonempty outputs when scripting. Unknown CLI options are ignored.
Input readers can accept partial numeric data; verify logged geometry counts.

Output indices use the number of existing files plus one. Removing selected
old files, mixing unrelated files, or running processes concurrently in one
output directory can cause filename collisions and truncation. Use a fresh
output and log directory for each run and a unique grid-export path.

Grid storage alone needs `8 * nx * ny * nz` bytes, before other allocations.
Signal/phase storage grows with proton count and recorded time samples;
trajectory recording adds further memory and disk use. Pilot calculations and
geometry placement can dominate short runs. There is no universal runtime bound.

Assess convergence using proton count, time step, grid resolution and
independent seeds. The tests check implemented behavior and selected numerical
properties; they do not establish experimental agreement for a particular
material or configuration. Check the [MSE output limitation](#signal-files)
before using the separately sampled MSE file.

## Verification

Build all targets and run correctness tests in the selected build directory:

```bash
cmake --preset linux-gcc-release
cmake --build build/linux-gcc-release --parallel
ctest --test-dir build/linux-gcc-release --output-on-failure -LE benchmark
python3 test/publication_examples.py \
  --executable build/linux-gcc-release/src/ionp/IONP \
  --source-dir . --work-dir build/publication-example-check
```

Use Python 3.10+ for the example check; it requires only the standard library.
It checks completed runs and outputs, repeatability, geometry and MIS signal
reconstruction. Choose a new or empty work directory. CTest covers parsing,
cell-load radii, reflection, importance sampling, custom geometry selection,
surrounding-voxel replication and field evaluation. Test assertions are active
in both Debug and release configurations.

Run the hardware-sensitive benchmark separately:

```bash
ctest --test-dir build/linux-gcc-release -L benchmark --verbose
```

It uses 75,000 IONPs in a 171-micrometer volume with 250 grid nodes per axis
and requires near-field queries to be at least 100 times faster than direct
all-particle evaluation. It uses a zero-valued grid and excludes grid
construction; it measures query cost, not total simulation speed or numerical
accuracy. Shared-machine timing can fail the threshold even when correctness
tests pass. The manual benchmark workflow retains this threshold.
