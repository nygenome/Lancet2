---
description: Run the full Lancet2 pipeline on both test profiles back-to-back. Germline (NA12878 / chr1) then somatic (HCC1395 tumor + HCC1395BL normal / chr4). Confirms exit-0 and a reasonable variant count for each; truth comparison is out of scope.
allowed-tools: Bash
---

# /e2e-pipeline-test — end-to-end smoke test

Run the full Lancet2 pipeline twice in succession against the
project's test data. The germline stage exercises the single-sample
(case-only) code path on `NA12878.final.cram` against GRCh38 chr1.
The somatic stage exercises the primary tumor/normal case-control
path on the `HCC1395` / `HCC1395BL` pair against GRCh38 chr4. Both
stages confirm a clean exit and report the variant count;
truth-set comparison is a separate harness and is not done here.

Use this command after completing a change that touches the
pipeline-level code in `src/lancet/core/` or any layer below, to
confirm the binary still produces variant calls. Use it before
invoking `fresh-reviewer` on a substantive change as a final
sanity check. Use the `LANCET_TEST_*_REGION_SMALL` regions (1 Mb
subsets) via the `sanitizer-build-analysis` skill for fast iteration;
`/e2e-pipeline-test` always runs the full-chromosome regions.

## Prerequisites

The test data must be downloaded once: `bash data/download_test_data.sh`
from the repository root pulls everything from
`gs://lancet2-test-datasets/test_harness_data/` into `data/`. The
required env vars (`LANCET_TEST_GERMLINE_*`, `LANCET_TEST_SOMATIC_*`)
are configured in `.claude/settings.local.json` (gitignored); copy
`.claude/settings.local.json.example` and adjust if your checkout
location differs.

## CLI shape assumption

The single-sample germline invocation passes the lone sample via
`--normal` (the CLI's `--normal` flag is `required=true`). The
case-control mode flag in `pipeline_runner.cpp` is computed as
`has_label(CASE) && has_label(CTRL)`, so a control-only run is
supported and the SHARED/CTRL/CASE INFO headers will not be
emitted. The equivalent advanced form is `--sample <path>:ctrl`;
both work. The script verifies the expected flags against
`Lancet2 pipeline --help` before running; if the help text shows a
different flag (e.g., `--single-sample`), the script reports the
mismatch and stops rather than failing in an opaque way deep in
the pipeline.

## What this command does

The full execution lives in `.claude/scripts/e2e-pipeline-test.sh` and is
invoked here as a thin wrapper:

```bash
bash "$CLAUDE_PROJECT_DIR/.claude/scripts/e2e-pipeline-test.sh"
```

The script's outline:

1. Build the release binary if missing.
2. Verify the CLI shape against `Lancet2 pipeline --help`.
3. Validate per-profile env vars and data-file existence.
4. Run the germline stage; capture variant counts.
5. Run the somatic stage; capture variant counts.
6. Print a summary table with both stages' result, wall time, and
   raw variant count. PASS/FAIL filtering is downstream-only — Lancet2
   emits FILTER="." on every record by design — so no PASS column.

Stages skip cleanly with a `skipped` rather than `fail` status
if their required env vars are unset or data files are absent.

## What the user sees

A successful run prints two stage banners with file paths and
elapsed wall times, then a final summary table. A skipped stage
shows `skipped` so the developer can tell apart missing-data
situations from real failures.

## Maintenance

The script's CLI-shape verification (`for required_flag in
--tumor --normal --reference --region --out-vcfgz`) tracks the
flag surface of `src/lancet/cli/cli_interface.cpp`. If a flag is
renamed or removed there, update the verification list in
`.claude/scripts/e2e-pipeline-test.sh`. The `/audit-bundle` slash command
catches drift between this script and `cli_interface.cpp`.

The variant-count parsing depends on `bcftools view`'s output
format from the `hts-tools` pixi env. If the env is renamed or
bcftools is replaced, update the script's `pixi run -e
hts-tools bcftools view` invocations.
