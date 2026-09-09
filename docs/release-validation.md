# Release validation: 0.0.1

Validated on 2026-09-06 using Ubuntu 22.04.5 LTS, x86-64, an Intel Core i7-12700,
GCC 11.4.0, CMake 3.22.1, Ninja 1.13.0 (Kitware build), and GCC OpenMP 4.5.
Python checks used Python 3.10.12. Plotting used NumPy 2.2.6 and Matplotlib
3.10.8 with the headless Agg backend.

## Initial publication checks

These checks describe the initial packaging pass, before the subsequent source
readability work recorded below.

| Check | Result |
| --- | --- |
| Fresh exported-source Debug build | Passed without compiler warnings |
| Fresh exported-source RelWithDebInfo build | Passed without compiler warnings |
| Correctness CTest suite | 6/6 passed in each build |
| Release assertion settings | Verified `-UNDEBUG` on test compilation only; simulation compilation retains its existing release settings |
| Required OpenMP dependency | A deliberately disabled OpenMP discovery fails during configuration |
| Simulation source and headers at initial packaging | All 36 files byte-identical to that working-tree baseline |
| Existing source, inputs, vendored files, plotting tools and fixtures | 259-file preservation check passed |
| Before/after runtime comparisons | All 14 fixtures matched byte-for-byte in Debug and RelWithDebInfo: 28 comparisons |
| Synthetic examples | FID, SE, MSE and cell importance sampling passed twice per build, with identical repeated outputs |
| Cell MIS reconstruction | Sum of complex proposal contributions reconstructed the magnitude within `1e-10` relative / `1e-12` absolute tolerance |
| Plotting | Single-grid and grid-comparison CLI outputs rendered; identical-grid difference was zero; trajectory CLI rendered through a headless `show` replacement |
| Exporter tests | Curated membership, current-file copying, deterministic archives, existing-output protection, invalid/missing input and symlink rejection passed |
| Documentation | Local links and heading anchors checked; documented examples executed |
| Source manifest/archive | Included-file hashes and exact archive membership checked; private/generated paths excluded |

Runtime comparisons used the preserved pre-cleanup binaries and freshly built
exported binaries, with the same compiler, build configuration, seed 123 and
fixture-specific thread count. They compared all generated `.dat` files,
including optional trajectories and importance-region diagnostics. The fixtures
exercise stationary and moving IONPs, coatings, bound protons, cells,
permeability, variable steps, importance sampling, cube signals and MSE, with
one- and four-thread cases. Timestamped logs were excluded from comparisons.

The example checks also validated three synthetic IONP coordinates, finite
11-row full signals, grid dimensions and values, trajectory proton IDs, and
the five sampled MSE times. The examples are execution checks, not converged
research datasets.

## Source readability cleanup

A subsequent pass on the same date added file/API and algorithm comments,
removed obsolete commented-out implementations, clarified legacy behavior,
and tidied whitespace and long declarations. It also added the
[source guide](source-guide.md). Production edits were limited to comments
and whitespace; no executable C++ tokens were changed.

| Check | Result |
| --- | --- |
| Source/header lexical comparison | All 92,280 non-comment tokens identical across 37 files, including the CMake header template and legacy wrappers |
| Exported-source GCC builds | Debug and RelWithDebInfo built without warnings, independently of the development tree |
| Correctness CTest suite (`-LE benchmark`) | 6/6 passed in each build |
| Seeded runtime comparisons | All 14 fixtures matched the preserved pre-readability outputs in both builds: 28 comparisons, 118 byte-identical `.dat` files |
| Synthetic examples | All four passed twice per build, including completion diagnostics, finite signals, expected geometry, grids, trajectories and repeated-output equality |
| Cell MIS reconstruction | Passed within the example check's existing output-precision tolerance |
| Exporter and documentation | All three exporter tests passed; local Markdown paths/anchors checked, including the new source guide |

The lexical comparison retained literals, identifiers, numbers and operators
and discarded only comments and whitespace. Runtime comparisons again used
seed 123 and the unchanged fixture-specific thread counts with the same GCC
and build configurations. The source-readability package contained the same
source/header bytes as that independently built copy. Study inputs, vendored
code and the previous publication snapshot were preserved.

Plotting tools and instructions were unchanged by this pass; their completed
checks are recorded above. The timing benchmark was not rerun for these
comment/whitespace edits; the earlier measurements below remain historical
results, not new performance claims.

## Restored default placement

At the maintainer's request, historical placement (previously selected with
`--legacy-placement`) was restored as the default after the readability pass.
The executable defaults to legacy placement, and both `SimSpace` placement
functions default to `use_overlap_grid=false`. The flag remains accepted for
existing commands. C++ callers can still request the overlap-grid path
explicitly. The placement implementations themselves were not changed.

This is an intentional change to default selection, separate from the
comment-only cleanup above. Both local GCC builds were updated, and all six
correctness tests passed in Debug and RelWithDebInfo.

For each build, 19 seeded cases compared the preserved pre-change executable
with `--legacy-placement`, the updated executable with no placement flag, and
the updated executable with `--legacy-placement`. All three produced identical
data: 74 files per build, including trajectories and optional regional signals.
The cases comprised the 14 runtime fixtures plus five short synthetic generated
geometries covering intracellular/extracellular count and concentration paths
and a mixed population. Runs used seed 123, matched build configurations and
unchanged case-specific thread counts. Completion diagnostics and finite signals
were checked in addition to process exit status.

All four documented examples also passed twice per build, including seeded
repeat equality and MIS reconstruction. No timing benchmark was rerun for this
default-selection change.

## Separate performance result

One initial-packaging benchmark invocation per configuration passed the unchanged
100× query-speedup threshold:

| Build | Direct queries | Corrected-hybrid queries | Ratio |
| --- | --- | --- | --- |
| Debug | 0.268431 s | 0.00108941 s | 246.401× |
| RelWithDebInfo | 0.0527408 s | 0.000381767 s | 138.149× |

The initial RelWithDebInfo audit had measured 95.9637× and failed the same
threshold. These individual timing measurements show host variability; they
do not demonstrate a performance improvement from this cleanup. The benchmark
times near-field queries with a zero-valued grid and excludes grid construction.
It remains separate from routine correctness CI.

## Limits of this validation

Native Windows/MSVC and Clang builds were not executed locally. Windows CI is
advisory pending the documented source portability review and runtime setup.
GitHub workflows were inspected locally but have not been run remotely during
this preparation. Interactive desktop plotting was not exercised.

The existing CLI, output-collision, input-parsing and MSE sampled-output
limitations remain documented in [simulations.md](simulations.md#reproducibility-and-limitations).
Apart from the requested placement-default selection described above, no
numerical or runtime fixes were introduced. Identical results across other
toolchains and scientific validity for a particular study are not established
by these checks.
