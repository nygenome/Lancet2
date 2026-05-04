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
7. Print an "ARTIFACT MANIFEST:" line followed by the per-stage VCF
   output paths, their `.tbi` indexes, and the per-stage tee'd log
   files. The manifest path itself is the file at
   `/tmp/lancet2_e2e_artifacts_<run_ts>.list`.

Stages skip cleanly with a `skipped` rather than `fail` status
if their required env vars are unset or data files are absent.

**Post-run cleanup prompt (mandatory).** After the script returns,
the assistant MUST parse the `ARTIFACT MANIFEST:` line from the
script's output, locate the manifest file, and ask the user via
the `AskUserQuestion` tool whether to delete the listed artifacts.
This step exists because:

- The script does NOT auto-delete its `/tmp/` outputs. Multi-GB VCFs
  or tens-of-MBs log files would otherwise vanish on the next OS
  /tmp-reaper sweep without the developer being asked.
- A subset of post-run analyses (truth-set comparison, profile
  attribution, debug archaeology) explicitly need the artifacts to
  remain. The default answer is "cleanup" because the common case
  is a smoke-test run, but the prompt MUST be presented so the
  developer can opt out.

The `AskUserQuestion` shape:

- `question`: "Clean up the e2e artifacts listed above?"
- `multiSelect`: false
- `options`:
  - `Yes, delete them all` — `rm -f` each path in the manifest, then
    `rm -f` the manifest file itself. Skip paths that no longer exist.
  - `Keep, I need them for analysis` — leave every file in place and
    print the manifest path so the developer can act later.

The deletions go through the project's `block_dangerous_bash` hook,
so each `rm` must use an absolute path. `/tmp/` is in the hook's
approved deletion roots; pass each path on its own command, not via
glob expansion.

## What the user sees

A successful run prints two stage banners with file paths and
elapsed wall times, live per-window EtaTimer progress (the script
uses `tee` so progress is visible while the stage runs, not just
at completion), a final summary table, and an artifact manifest
listing the VCFs / indexes / per-stage logs the run produced.

The run ends with a single `AskUserQuestion` prompt asking the
developer whether to clean up the manifest's contents. A skipped
stage shows `skipped` so the developer can tell apart missing-data
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
