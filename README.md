# Tension Distribution for a Planar Four-Cable Robot

MATLAB implementations and reproducible comparisons for computing cable tensions subject to force equilibrium and lower/upper tension bounds.

The proposed **Direct Method (DM)** minimizes the tension L2 norm for a feasible force request by evaluating at most two scalar candidates. Its geometric derivation assumes four counterclockwise-ordered anchors and a platform position strictly inside their convex frame.

## Quick start

The small example needs **MATLAB only**. The package was tested with MATLAB R2024b Update 4; no Python, compilation, or additional downloads are needed for this example.

1. Download and extract the repository ZIP, or clone the repository.
2. Set MATLAB's **Current Folder** to the folder containing this README and run:

```matlab
addpath('examples');
result = run_example;
```

The example displays the requested force, four computed cable tensions, achieved force, equilibrium error, and whether the tension bounds are satisfied. It also returns these values in `result`.

**Expected result:** all example checks pass, tension bounds are satisfied, and equilibrium error is below `1e-8` N. The example uses a known feasible request, positions in **mm**, forces and tensions in **N**, and equal cable bounds of **3.8–25 N**. It records no performance timings.

### Compare all six methods

With **Optimization Toolbox** installed, run:

```matlab
comparison = run_example('all');
```

| Method | Implementation |
| --- | --- |
| Direct Method (DM) | [`wt_direct.m`](src/wt_direct.m) |
| Improved closed-form method (ICFM) | [`wt_pott.m`](src/wt_pott.m) |
| Bounded quadratic programming (QP) | [`wt_qp_bounded.m`](src/wt_qp_bounded.m) |
| Vertex-based VTDA-L2 | [`wt_gouttefarde.m`](src/wt_gouttefarde.m) |
| Analytic-centre ACTA | [`wt_acta.m`](src/wt_acta.m) |
| Extended X-ACTA | [`wt_xacta.m`](src/wt_xacta.m) |

The first five methods receive a shared wrench-exertion-capability (WEC) target; X-ACTA receives the original force. WEC leaves this example's feasible request unchanged. ACTA and X-ACTA are equation-based MATLAB implementations, not verified copies of the original authors' software. Their analytic-centre objectives differ from DM's tension L2 objective. See the [method sources](docs/METHOD_SOURCES.txt) for attribution.

## Requirements

| Task | Dependencies |
| --- | --- |
| Direct example | MATLAB |
| All-method example | MATLAB and Optimization Toolbox |
| New benchmark capture | MATLAB, Optimization Toolbox, Parallel Computing Toolbox, and Java |
| Saved-data reproduction | Base MATLAB and Java; original saved data supplied separately |

## Package guide

| Location | Contents |
| --- | --- |
| [`src/`](src/) | Seven numerical functions, with input/output contracts and returned statuses in their help text |
| [`examples/`](examples/) | Small example with independent physical checks |
| [Function reference](docs/FUNCTION_REFERENCE.txt) | Dimensions, units, cable ordering, bounds, and method-specific statuses |
| [Method sources](docs/METHOD_SOURCES.txt) | Comparison methods and their references |
| [Validation](docs/VALIDATION.txt) | Tested commands, results, and scope of verification |
| [`benchmark/`](benchmark/README.md) | Capture and saved-data reproduction instructions, settings, and costs |
| [`provenance/`](provenance/) | Original/published source identities and function-name mapping |
| [Changes](CHANGES.txt) · [Notices](NOTICES.txt) | Presentation changes, source attribution, and license information |

Returned tensions follow the original anchor-column order. A finite result or a legacy status flag alone does not establish success; independently check equilibrium and tension bounds. WEC can change an infeasible request's magnitude or signed direction, some force lines have no target, and zero requests have a separate fallback. Read the [function reference](docs/FUNCTION_REFERENCE.txt) before supplying different inputs.

## Benchmark reproduction

### Small capture and saved-data replay

Start with the **192-request smoke test**, which uses two process workers. Then reconstruct its tables from the saved outputs:

```matlab
addpath('benchmark');
capture = run_benchmark(fullfile(pwd,'smoke_results'));
tables = reproduce_saved_results(fullfile(pwd,'smoke_results'), ...
    fullfile(pwd,'smoke_tables'));
```

Both output directories must be new. Worker startup can dominate this small capture. Saved-data reproduction copies its inputs and makes no new solver or performance-timing calls.

### Full experiment: substantial computational cost

The full experiment has **30,965,760 requests** and eight recorded timing series. Its original capture loop took **5,843.6 seconds (about 97.4 minutes)** on eight workers of the reported i9-12900K system, excluding startup and later analysis.

- Original raw data: approximately **10.52 GB**.
- New full capture: conservatively allow **40 GB** of free disk space.
- Independent saved-data replay: allow a further **12 GB** and several GB of working memory.

These are resource estimates, not completion-time guarantees. Read the [full benchmark instructions](benchmark/README.md) before starting a full capture or reproducing the published results from the separately supplied data archive. The example and default smoke command never start a full run. Large raw data and generated results are excluded from this repository and package.

## Validation and numerical preservation

The numerical statements are preserved after the documented function-name mapping. Public function names have no revision suffixes; solver settings and timing boundaries are unchanged. Setup and validation wrappers operate outside the measured routines.

The clean-extraction examples, 192-request smoke capture, and saved-data replay passed. See [validation details](docs/VALIDATION.txt) for the historical full-data checks and the limit on full-run testing, and [changes](CHANGES.txt) for the documented presentation updates.
