BENCHMARK CAPTURE AND SAVED-DATA REPRODUCTION

This folder provides two distinct entry points. run_benchmark makes new timing
observations; reproduce_saved_results only analyzes an existing saved capture.
The default capture is a small 192-input smoke check, not the full experiment.
The measured harness is supplied in captured/ with revision suffixes removed
from function names and their references. Generated log/default-directory names
were also updated. The summary helper omits an unused plot-source dependency;
all numerical calculations, settings and timing boundaries are unchanged. The
wrapper stages the harness beside copies of the seven package src/ functions
so its source-resolution checks continue to work. Source identities and the
original-to-current name mapping are in provenance/SOURCE_MANIFEST.csv.

Quick test from the extracted package root

Use fresh directory names for each attempt. These paths are portable MATLAB
examples; neither wrapper contains a machine-specific input or output path.

addpath(fullfile(pwd,'benchmark'));
captureDir = fullfile(tempdir,'reviewer_smoke_capture_01');
analysisDir = fullfile(tempdir,'reviewer_smoke_analysis_01');
run_benchmark(captureDir,struct('profile','smoke','workers',2,'pin_pcores',false));
reproduce_saved_results(captureDir,analysisDir);

The smoke profile selects 4 by 4 poses, 3 magnitudes and 4 directions: 192
input rows. It exercises the frame edge as well as interior poses. Expected
algorithm failures are recorded and do not make the capture incomplete. The
final capture manifest must have status raw_complete; the saved reproduction
receipt must have status passed and input_files_unchanged=true. Startup and
worker warmup dominate this small test; allow a few minutes. It does not
reproduce the full-run counts or timing distributions.

MATLAB requirements

Validated release: R2024b Update 4 on 64-bit Windows. Both wrappers need Java
for SHA-256 and path handling. run_benchmark needs Parallel Computing Toolbox
for process workers and Optimization Toolbox for the bounded-QP comparator.
It closes a pool it created. An existing process pool must have the requested
worker count and remains open. The wrappers restore the caller's path and
working directory. The captured harness restores worker warning states and
any affinity settings it changes.

Saved-data analysis alone needs base MATLAB R2024b, including tables, quantile
and v7.3 MAT-file support; it does not call Optimization or Parallel Computing
Toolbox. The standalone Direct example elsewhere in this package also avoids
those benchmark-only dependencies. Linux/macOS skip Windows affinity changes;
their timing measurements constitute a new hardware/software experiment.

Reproduce existing saved results without new timings

Call with an extracted saved-run directory and a new, disjoint output tree:

reproduce_saved_results(savedRunDir,newAnalysisDir);

The input requires manifest.mat, geometry.mat, wec.mat, dm_postwec.mat,
icfm_postwec.mat, qp_postwec.mat, vtda_postwec.mat, acta_postwec.mat,
xacta_raw.mat and wecdm_raw.mat. The measured series must be complete.
The wrapper hashes every required file before copying, verifies each physical
copy, then checks the original again after analysis. It places the copies
and captured-helper outputs under newAnalysisDir/captured_replay.
Fail-fast stubs shadow all seven public solver functions during replay.
No timing observation is recomputed or replaced, and attempted failures and
observed long tails stay in their original timing populations.

Nine selected CSVs are published at the new output root:
T_domain, T_timing, T_quality, T_postwec_components, T_iterative_domains,
T_wec_behavior, T_pipeline_parity, T_force_approximation and
T_complete_paths_matched. Their names end in .csv. saved_reproduction_receipt.json
records input identities, source hashes, options and successful completion.
These hashes establish which inputs were replayed; without the original
inventory they do not independently authenticate a downloaded historical run.

The captured summary helper also writes T_L2_matched.csv and a descriptive
report inside captured_replay. This is the historical QP-relative analysis
with different filtering, NOT the final DM-relative norm comparison. It is
not published at the output root. The wrapper never generates the old captured
diagnostic plots, and their plotting source is not shipped. The summary's
unconditional provenance-source list was shortened by one filename so it can
run with plotting disabled; its numeric calculations are unchanged. The
source identities and metadata edits are in ANALYSIS_METADATA_CHANGE.txt and
provenance/SOURCE_MANIFEST.csv at the package root.
No generated historical figures, run logs or rebuttal text are shipped here.

Two optional analyses

For the paired interior-force subsets, set options.interior_table=true.
This additionally requires the original independent_readback_audit.mat and
publishes T_complete_paths_matched_interiors.csv after checking its all-seven
method physical-audit completion. This does not rerun that independent audit.

For the final actual-tension L2 difference, set options.norm_statistics=true
and options.inventory_file to the original FULL_DATA_INVENTORY.csv. This is
supported only for the original complete 30,965,760-row saved run:

opts = struct('norm_statistics',true,'inventory_file',inventoryFile);
reproduce_saved_results(savedFullRunDir,newAnalysisDir,opts);

The original archive is named ACTA_full_raw_data.zip (about 10.6 GB), supplied
separately by the authors. Its recorded SHA-256 is
3e99e0ca71b8096354d51749400054dee03af5bb78a6a2ca37af415023afe933.
Extract it and use its review_v11_full directory as savedFullRunDir.
The matching original inventory is included as
provenance/FULL_DATA_INVENTORY.csv at the package root. For example, set
inventoryFile = fullfile(pwd,'provenance','FULL_DATA_INVENTORY.csv') before
running the optional command. No public data-download URL is available here.

The optional norm postprocessor writes newAnalysisDir/norm_analysis and uses
the original 24,433,774-row common unchanged-force starting cohort. It
uses delta = norm(T_method,2) - norm(T_DM,2), in N, from the actual returned
tension vectors. It recomputes norms and physical checks, preserves signed
differences and saves exact per-method valid/excluded masks. It needs the ten
raw files above plus original masks.mat, T_quality.csv, completion_receipt.json,
and these original source_snapshot files: Compare_Review_v11.m,
wt4_2024_minmax_v4.m, wt_pott_v2.m, wt_qp_bounded_v11.m and wt_gouttefarde_v2.m.
The inventory validates those original captured identities; these historical
filenames must remain unchanged in the saved-run directory. They are not the
renamed, documented copies delivered in src/. No optimizer or tension
distribution function is invoked.

It can also run separately without the general timing replay:

compute_norm_differences_saved(savedFullRunDir,newNormDir,inventoryFile);

Full-run dimensions, resources and cost

The measured full grid contains 64 x 64 poses in [0,315]^2 mm, 21 force
magnitudes in [0,10] N and 360 directions: 30,965,760 original input rows.
There are 30,005,640 spatial-interior rows and 960,120 frame-edge rows.
Zero x/y axes are interior symmetry axes; only x=315 or y=315 is a frame edge.
The six-method common unchanged-force timing population has 24,433,774 rows.

The original host used Intel Core i9-12900K hardware, 16 physical/24 logical
processors, MATLAB R2024b Update 4 and eight process workers. Each worker
completed 192 warmup boundary calls (eight repetitions of three representative
requests across WEC plus seven methods/pipelines). Warmups are discarded;
their internal solver iterations are not repetitions of the measured grid.
The host-specific affinity mapping was applied and restored. See
SETTINGS_AND_TIMING.txt for the recorded policy and numerical acceptance.

The captured loop took 5,843.591001 seconds (97.393 minutes), including
parallel scheduling, independent post-checks and raw storage inside that
loop, excluding pool startup, warmup and later analysis. This is the observed
historical cost, not a runtime guarantee. The ten raw files occupy about
10.52 GB in their saved compressed representation. Saved replay makes real
copies: allow about 12 GB additional free disk plus several GB of working RAM
(roughly 5-7 GB for the full summary's arrays and quantile temporaries; this is
an estimate, not a recorded peak). Optional norm analysis adds working memory
and hashing/reading time but produces compact masks, not a full delta array.
Small smoke runs use a tiny fraction of these storage/memory requirements.

A new full capture is expensive and is never started by either default:

run_benchmark(newFullDir,struct('profile','full','workers',8,'pin_pcores',true));

Use sufficient headroom for potentially less-compressible v7.3 output (about
40 GB is a conservative full-capture disk allowance); actual storage depends
on the platform and outputs. The explicit representative profile is 436,968
rows. Its cost and distributions differ from both smoke and full. Fresh
outputs are required; no resume or partial-file summarization is attempted.

Important interpretation

All per-method time fields are raw elapsed seconds; published table times are
microseconds. Ratios of table medians are ratios of marginal medians, not
medians of paired ratios. Observed maxima are not WCET guarantees. A single
pass removes repeated-grid best-of selection; it does not eliminate internal
Newton iterations or pre-measurement warmup. Shared-WEC tension distribution timing and
complete raw-request pipeline timing have different boundaries and must not
be interchanged. See SETTINGS_AND_TIMING.txt before using a timing claim.

Function reference

run_benchmark(outputDir,options): outputDir is a fresh path; options is a scalar
struct. Supported fields are profile (smoke/representative/full), positive
integer workers/chunk_size, logical pin_pcores, and integer seed. The returned
manifest gives ninputs, nposes, nforces, completed_chunks, measured_attempts,
source identities and status='raw_complete' only after successful capture.
MAT raw arrays use N rows: tensions are N-by-4, elapsed/error/norm arrays are
N-by-1, diagnostic matrices are N-by-20 where applicable. Native failures are
data, separate from a failed or incomplete benchmark execution.

reproduce_saved_results(inputDir,outputDir,options): paths are disjoint and the
output is fresh. options supports logical interior_table/norm_statistics and
inventory_file (path). Its result has status='passed', ninputs, published_tables,
input_files_unchanged and before/copy/after hashes. It throws on incomplete
captures, missing inputs, altered hashes or helper failures. It does not prove
mathematical correctness independently of the saved checks; optional norm
analysis does independently reconstruct geometry and per-pair acceptance.

compute_norm_differences_saved(inputDir,outputDir,inventoryFile): optional
full-run DM-relative analysis; explicit paths are recommended. Its returned
receipt has status='passed' only after checking the original inventory and
all physical mask/readback invariants. It writes three method rows with
median/q2.5/q97.5/min/max in N, counts and exact masks. In its mask MAT,
candidate_mask is N-by-1, pair_valid_mask and exclusion_reason_bits are N-by-3
in ICFM/bounded-QP/VTDA order. Native positive QP status, bound <=1e-9 N and
raw-force 2-norm <=1e-8 N are required; signed differences never filter rows.
Quantiles use Hyndman-Fan type 5, h=N*p+0.5, with interpolation and clamping.

captured/ is a private provenance dependency, not an additional public API.
compare_methods(config) is the renamed measured harness delegated to by the
capture wrapper. summarize_results(folder,false) consumes the complete N-row
capture and returns tables/summary, recording analysis_complete separately.
analyze_interior_infeasible(folder) returns two paired-error table rows
after checking the existing seven-method independent audit. The summary's
plotting branch is unused by the wrapper and its plot helper is not shipped;
call the documented wrapper, which always passes false for plotting.
