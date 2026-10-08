TENSION DISTRIBUTION FOR A PLANAR FOUR-CABLE ROBOT

Start here

The proposed Direct Method (DM) computes cable tensions for a feasible force
request by evaluating at most two scalar candidates. It minimizes the actual
tension L2 norm subject to equilibrium and lower/upper tension bounds. Its
geometric derivation assumes four counterclockwise-ordered anchors and a
platform position strictly inside their convex frame.

The simplest example needs MATLAB only. This package was tested with MATLAB
R2024b Update 4. No Python, compilation, download or path editing is required.

1. Extract the ZIP and set MATLAB's Current Folder to the extracted folder
   containing this README.
2. Run these two lines in the Command Window:

   addpath('examples');
   result = run_example;

The example prints the requested force, four computed cable tensions, achieved
force, equilibrium error and whether the tension bounds are satisfied. It uses
a known feasible request, positions in mm, forces/tensions in N, and equal
3.8--25 N cable bounds. Expected: bounds satisfied = yes; equilibrium error
below 1e-8 N; all example checks passed = yes. The returned result contains the
same values. This example records no performance timings.

Optional comparison example

With Optimization Toolbox installed, use:

   comparison = run_example('all');

This prints the same checks for DM, improved closed-form method (ICFM), bounded
quadratic programming (QP), vertex-based VTDA-L2, analytic-centre ACTA and
extended X-ACTA. The first five use a shared wrench-exertion-capability (WEC)
target; X-ACTA receives the original force. The example's feasible request is
unchanged by WEC. ACTA and X-ACTA are equation-based MATLAB implementations,
not verified copies of the original authors' software. They use analytic-centre
objectives, not DM's tension L2 objective. Sources and exact APIs are in docs/.

Where to look

  src/         Seven numerical functions; help text lists inputs and statuses.
  examples/    The small example above, including independent physical checks.
  docs/        Function contracts, method references and test summary.
  benchmark/   Capture/reproduction entry points, settings and cost information.
  provenance/  Original/published source identities and executable-line checks.
  CHANGES.txt  Presentation changes and new setup/validation wrappers.
  NOTICES.txt  Source attribution and license information.

The returned tensions follow the original anchor-column order. Independent
checks are necessary: a finite result or a legacy status flag alone is not a
success certificate. WEC can change an infeasible request's magnitude or signed
direction; some force lines have no target, and zero requests have a separate
fallback. See docs/FUNCTION_REFERENCE.txt before supplying different inputs.

Benchmark reproduction -- optional and substantially more expensive

New benchmark capture requires Optimization Toolbox AND Parallel Computing
Toolbox. The captured harness uses process workers. Start with the 192-request
smoke run, then reconstruct its saved tables without rerunning the methods:

   addpath('benchmark');
   capture = run_benchmark(fullfile(pwd,'smoke_results'));
   tables = reproduce_saved_results(fullfile(pwd,'smoke_results'), ...
       fullfile(pwd,'smoke_tables'));

Both output directories must be new. The smoke run uses two process workers;
startup can dominate this small run. Saved-data reproduction uses base MATLAB
and Java, copies its inputs, and makes no new solver or performance-timing calls.

The FULL experiment has 30,965,760 requests and eight recorded timing series.
The original capture loop took 5,843.6 seconds (about 97.4 minutes) on eight
workers of the reported i9-12900K system, excluding startup and later analysis.
The original raw files occupy about 10.52 GB; conservatively allow 40 GB for
a new full capture and a further 12 GB for an independent saved-data replay,
plus several GB of working memory. Runtime depends on the
machine; these are resource estimates, not completion-time guarantees.

Read benchmark/README.txt before requesting a full run or reproducing the
published results from the separately supplied original data archive. A full
run is never started by the example or by the default smoke command. Large
raw data and generated results are intentionally not included in this ZIP.

Validation and preservation

The published numerical functions differ from the measured sources only in
full-line comments and blank lines. Their executable lines, solver settings
and timing boundaries are preserved. New wrappers handle paths, fresh output
directories and validation outside the measured routines. See CHANGES.txt and
docs/VALIDATION.txt for the tested commands and the full-run testing limit.
