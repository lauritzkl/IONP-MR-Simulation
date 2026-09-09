# Reading the simulation source

Start with the [working examples](../examples/README.md) and
[simulation reference](simulations.md) for the user interface. This guide maps
that interface to the implementation. Source comments describe the current
behavior, including legacy details; the readability cleanup does not change
the numerical model, configuration, random draws or output formats.

## Follow one run

1. [main.cpp](../src/ionp/main.cpp) resolves the CLI and input/output paths,
   selects a single configuration or a sorted batch, and calls
   `run_single_simulation`. A requested seed is set once for the invocation.
2. [Data.cpp](../src/ionp/Data.cpp) normalizes configuration keys, converts
   values, checks supported combinations and resolves defaults. Auxiliary
   readers load particle positions, cell positions/loads and histograms.
3. [SimSpace.cpp](../src/ionp/SimSpace.cpp) constructs the requested physical
   geometry. Placement retries reject overlaps; cell and extracellular
   populations retain their own input and distribution rules.
4. `main.cpp` resolves the integration step and material parameters, constructs
   [MRsequence](../src/ionp/MRsequence.hpp), and calls `Sequence`.
5. [MRsequence.cpp](../src/ionp/MRsequence.cpp) prepares magnetic sources and
   any field grid, assigns free/bound samples, propagates proton positions and
   phases, and combines the recorded phases into complex signals.
6. `Data` writes the signal and geometry files. Optional trajectory and grid
   writers supply their respective diagnostics. Consult the reference for
   filenames and columns; the usual main signal file stores magnitude and time.

Historical (legacy) placement is the default in the executable and in the
`simulationSpace_inCell`/`simulationSpace_outCell` API defaults. The
`--legacy-placement` flag remains accepted for existing commands. The overlap
grid implementation remains available to C++ callers through an explicit
`use_overlap_grid=true` argument.

## Module map

| Files | Responsibility and useful starting points |
| --- | --- |
| [Point3D.hpp](../src/ionp/Point3D.hpp), [Particle.hpp](../src/ionp/Particle.hpp) | Coordinate arithmetic, initial/current positions, placement and displacement primitives. |
| [Ionp.hpp](../src/ionp/Ionp.hpp), [Aggregate.hpp](../src/ionp/Aggregate.hpp), [Proton.hpp](../src/ionp/Proton.hpp) | Core/coating geometry, aggregate placement regions, and spin phase/cell state. Their small `.cpp` files retain compatibility translation units. |
| [Data.hpp](../src/ionp/Data.hpp) | Parsed configuration, explicit-option flags, auxiliary input readers and output overloads. |
| [SimSpace.hpp](../src/ionp/SimSpace.hpp) | Cell/extracellular placement, aggregate construction and surrounding-voxel magnetic templates. |
| [SimMatrix.hpp](../src/ionp/SimMatrix.hpp), [SimMatrix.cpp](../src/ionp/SimMatrix.cpp) | Scalar dipole field, all-source grid construction, interpolation stencils and binary grid output. |
| [BFieldSolver.hpp](../src/ionp/BFieldSolver.hpp), [CorrectedHybridField.cpp](../src/ionp/CorrectedHybridField.cpp) | Solver selection and a grid with smooth exact-field corrections near sources. |
| [MRsequence.hpp](../src/ionp/MRsequence.hpp), [MRsequence.cpp](../src/ionp/MRsequence.cpp) | Spatial candidate indices, material boundaries, sequence timing, propagation and aggregation. |
| [ImportanceSampling.hpp](../src/ionp/ImportanceSampling.hpp), [ImportanceSampling.cpp](../src/ionp/ImportanceSampling.cpp) | Region/proposal labels, accessible volume estimates, allocations and estimator weights. |
| [RandNumberBetween.hpp](../src/ionp/RandNumberBetween.hpp), [RandomList.cpp](../src/ionp/RandomList.cpp) | Run-wide seed allocation and per-object random engines/distributions. |
| [OutputUtils.hpp](../src/ionp/OutputUtils.hpp) | Existing count-based output numbering. |
| [mptmacros.h](../include/mptmacros.h), [config.hpp.in](../configured_files/config.hpp.in) | Compiled model constants/casts and the CMake project-metadata template. |
| [Vector3D.hpp](../src/ionp/Vector3D.hpp), [profiler.h](../src/ionp/profiler.h) | Legacy endpoint-vector and timer helpers; neither defines the main diffusion algorithm. |

The three wrapper headers in `include/` retain historical external-library
includes. They do not add argument-parser or Eigen requirements to the default
simulator build. Vendored CSV source retains its original comments and license.

## Geometry, units and ownership

Positions, radii and displacements use metres; times use seconds; fields use
tesla; phase uses radians; diffusion coefficients use m²/s. Membrane permeability
parameters are crossing probabilities. Concentration inputs have their own units
listed in the [configuration reference](simulations.md#config-files).

`Particle::x0()` stores a reference/initial position and `xt()` the current
position. The displacement methods update `xt`; time bookkeeping and boundary
enforcement belong to the sequence. `Proton::inCell()` and `Cellnum()` together
describe current cell membership. The historical zero-argument `setPhase()`
method reads phase; `setPhase(value)` writes it.

`Ionp::radius()` is the magnetic core radius. `coating()` stores shell thickness,
while the configuration's coating value specifies the total outer radius.
`ionp_shell_thickness()` performs that conversion and `outer_radius()` returns
core radius plus thickness. These quantities have different roles in field
amplitude, collision tests and material diffusion.

`main.cpp` owns the run's configuration and geometry vectors. `MRsequence`
stores configuration and a copy of matrix metadata, while `Sequence` receives
mutable state and output vectors by reference. Field-grid storage is supplied
explicitly. `CorrectedHybridField` borrows the matrix and grid and copies source
coordinates/amplitudes: its borrowed objects must outlive it and remain stable
during queries. Spatial-query scratch buffers belong to individual callers.

Physical IONPs in the central voxel determine obstacles and bound-water jobs.
Surrounding-voxel copies, including optional compressed magnetic proxies,
contribute to the field separately. Compression preserves the represented sum
of core radii cubed; its proxies do not replace central collision geometry.

## Reading the simulation kernel

`MRsequence.cpp` is organized as helpers followed by the large `Sequence`
orchestrator. Read these groups in order:

- **Random streams and timing:** `ProtonRandomStream`, `proton_random_seed`,
  `build_refocusing_schedule`, `update_phase` and `phase_record_interval`.
  Refocusing changes the sign of phase at rounded base-step pulse indices.
  Integration steps and recorded-output intervals are distinct.
- **Candidate queries and materials:** `IonpSpatialIndex`, `CellSpatialIndex`,
  coating-state helpers and diffusion-coefficient selection. Candidate lists
  are conservative; exact sphere tests determine contact or membership.
- **Boundary walking:** sphere intersections and `walk_cell_boundaries`.
  A proposed displacement can encounter cell membranes, IONP cores/coatings
  and periodic voxel faces. The walker processes events in travelled-path
  order, reflects or crosses as appropriate, and retains remaining path length.
  Periodic cell images preserve the physical cell index. Failed walks request
  a retry from the saved state.
- **Initial sampling:** legacy joint-region and overlapping-proposal helpers
  reject inaccessible starts. Overlapping shells require coverage correction
  so their intersections do not receive excess probability.
- **Execution and aggregation:** `Sequence` prepares sources and sampling
  jobs, runs free and bound spins, checks the common recorded time axis, then
  combines `exp(i*phase)` in the existing order. Histories use
  `[proton][recorded time]` indexing. Optional diagnostic buffers share the
  corresponding signal's time axis.

The dipole calculation returns the z component for z-aligned moments.
`SimMatrix` stores grid values at `z*nx*ny + x*ny + y`, so y varies fastest.
Grid coordinates are centered on zero; the stored `x0` does not shift them.
Interpolation clamps grid indices rather than wrapping a point periodically.
The hybrid evaluator starts with full-grid interpolation and adds, for each
nearby source, a smooth weight times its exact field minus its interpolated
contribution. This subtraction avoids counting that source twice. Coincident
source/query centres remain singular in the dipole helpers.

## Sampling weights and reproducibility

Legacy importance sampling uses exclusive joint strata: cell region × distance
to intracellular IONPs × distance to extracellular IONPs. It averages within
each active stratum and then combines the means with region weights.

Overlapping proposals label the distribution that generated each initial
position. Several proposals can cover the same point. For a covering proposal
with sample count `n_k` and accessible volume fraction `f_k`, the mixture-density
term is `n_k/f_k`; the raw sample weight is the reciprocal of their sum. These
weights are normalized across free samples. Bound-water contributions are
combined separately. The complex `mis_contribution` outputs sum to the
corresponding combined free-water signal; ordinary proposal means do not.
Initial sampling labels and estimator weights remain fixed as spins move.

Volume estimation uses separate deterministic pilot generators. Propagation
streams derive from the run seed, proton index and fixed stream tags. Some
geometry helpers instead allocate seeds from a shared counter. In particular,
the legacy `RandomList` constructor deliberately requests two seeds even though
it uses only one. Preserve constructor order, stream tags, random draw counts,
batch order and reduction order when investigating reproducibility. Record
thread count, compiler, standard library, build configuration and complete
inputs as described in the [reproducibility notes](simulations.md#reproducibility-and-limitations).

## Retained implementation details

The comments call out details that should not be silently changed during a
readability edit:

- Several `Particle::placeRandom` overloads restore `x0` for `init=false`,
  sometimes after consuming random draws; other overloads draw a fresh `xt`.
  `Particle::move` normalizes its direction argument in place.
- `Vector3D` is a legacy two-endpoint helper. Its `atan` angle formulas omit
  quadrant correction, `stretch` changes only the final endpoint, and its
  equality/inequality operators currently use the same expression.
- `SimMatrix` copies geometry/save settings but does not copy its legacy public
  auxiliary vectors. Active field arrays are passed explicitly.
- The optional profiler relies on declarations supplied by its including
  context, keeps header-local static storage and has no thread synchronization.
- Unknown CLI options, numeric-prefix parsing, caught per-run failures,
  output-name reuse and historical MSE extraction formulas retain the
  [documented runtime limitations](simulations.md#reproducibility-and-limitations).

Use the [release validation report](release-validation.md) for the checks
performed on this snapshot. Execution and regression checks establish software
behavior; scientific use also requires convergence and model validation.
