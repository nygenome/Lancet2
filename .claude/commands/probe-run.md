# /probe-run — Step 2 of probe tracking: forensic Lancet2 invocation

Re-run Lancet2 with `--probe-variants` and `--probe-results` to collect forensic data on every missed truth variant. This is step 2 in the probe tracking workflow; step 1 (`/probe-concordance`) produced the missed variants list this step needs as input.

The probe system records raw facts per `(probe_id, window, comp_id, k)` attempt — no attribution logic, just a TSV of survival counts and structural flags. Step 3 (`/probe-analyze`) does the attribution. Both `--probe-variants` and `--probe-results` flags are bidirectionally required; without either, the entire probe system is inactive (zero overhead in production).

## When to use

Run this when `/probe-concordance` has produced `missed_variants.txt` and you want to know what happened to those missed variants inside Lancet2's pipeline. Run with `--verbose` so step 3's `--log` flag can sub-classify the `not_processed` category.

## Procedure

### Phase 1 — Verify inputs

The user typically supplies an output directory from step 1; default to the most recent `notes/probe-debug-<YYYY-MM-DD>/` directory if not specified. Check:

```bash
ls <outdir>/missed_variants.txt
```

If the file is missing, surface the error and suggest running `/probe-concordance` first. Do not proceed with synthetic input — that pattern silently produces wrong results.

If the user did not name `--normal`, `--tumor`, `--reference`, or output paths, gather them via `AskUserQuestion`. Suggest defaults from the `test-data-locations` skill. For single-sample germline, the lone CRAM goes as `--normal` (the CLI's `--normal` flag is `required=true`); the case-control mode flag in `pipeline_runner.cpp` only fires when both case and control labels are present. Do not silently use `--tumor` for a germline run — the CLI rejects that combination.

### Phase 2 — Build staleness check

Step 2 takes 10-30+ minutes on a full chromosome. Running it against a stale binary is expensive. The check has two parts:

```bash
# (a) Working tree must be clean.
if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
  echo "Working tree has staged/unstaged/untracked changes — rebuild required."
  needs_rebuild=true
fi

# (b) Embedded SHA must match HEAD.
if [ "$needs_rebuild" != "true" ] && [ -x cmake-build-release/Lancet2 ]; then
  embedded_sha=$(./cmake-build-release/Lancet2 --version 2>&1 | awk -F'-' '{print $NF}' | tr -d '[:space:]')
  head_sha=$(git rev-parse --short=10 HEAD)
  if [ "$embedded_sha" != "$head_sha" ]; then
    echo "Binary embedded SHA ($embedded_sha) != HEAD SHA ($head_sha) — rebuild required."
    needs_rebuild=true
  fi
fi

# (c) Binary missing.
if [ ! -x cmake-build-release/Lancet2 ]; then
  echo "Release binary missing — building first."
  needs_rebuild=true
fi
```

The version string format is `<tag>-<branch>-<short10-sha>`. The SHA is set at CMake **configure** time from `src/lancet/version.h.inc`. The pixi `build-release` task `depends-on` `configure-release`, so a single `pixi run build-release` covers both:

```bash
if [ "$needs_rebuild" = "true" ]; then
  pixi run build-release
  if [ $? -ne 0 ]; then
    echo "Build failed — aborting /probe-run."
    exit 1
  fi
fi
```

Tell the user when rebuild is needed and why, before invoking the build. Do not silently rebuild — the user may want to abort.

### Phase 3 — Verify CLI shape

Same pattern as `/e2e-pipeline-test`:

```bash
help_out=$(./cmake-build-release/Lancet2 pipeline --help 2>&1)
for required_flag in "--probe-variants" "--probe-results" "--normal" "--reference" "--region" "--out-vcfgz" "--verbose"; do
  if ! echo "$help_out" | grep -q -- "$required_flag"; then
    echo "❌ Lancet2 pipeline --help does not mention $required_flag."
    echo "   The CLI shape this command assumes has changed; update /probe-run."
    exit 1
  fi
done
```

This catches the case where flag names change in cli_interface.cpp without this command being updated.

### Phase 4 — Run the pipeline

```bash
out_vcf=<outdir>/probe_run_$(date +%s).vcf.gz
debug_log=<outdir>/lancet2_debug.log

./cmake-build-release/Lancet2 pipeline \
    --normal "$NORMAL_PATH" \
    [--tumor "$TUMOR_PATH"] \
    --reference "$REFERENCE_PATH" \
    --region "$REGION" \
    --probe-variants "<outdir>/missed_variants.txt" \
    --probe-results "<outdir>/probe_results.tsv" \
    --num-threads "$(nproc)" \
    --out-vcfgz "$out_vcf" \
    --verbose 2> "$debug_log"
```

`--verbose` is critical for step 3's sub-classification of `not_processed`. Even if the user said they don't need it, include it; the disk cost of the log is small compared to the analytical value it unlocks.

Stream stderr to the debug log file; do not silently swallow it. Tail a summary to chat as the run progresses (the script tails won't work directly through the slash command, so just surface key markers — start time, completion, exit code).

### Phase 5 — Confirm and hand off

After the run completes, verify:

```bash
ls -la <outdir>/probe_results.tsv <outdir>/lancet2_debug.log
wc -l <outdir>/probe_results.tsv
```

If `probe_results.tsv` doesn't exist or is empty, something went wrong inside the pipeline — surface the tail of the debug log.

End with:

> Step 2 complete. Probe results: `<outdir>/probe_results.tsv` (N rows). Debug log: `<outdir>/lancet2_debug.log`. Next: `/probe-analyze` to attribute each missed variant to a pipeline stage.

Do not auto-invoke step 3.

## When NOT to use this command

Do not use this command without `missed_variants.txt` from step 1. A hand-curated probe-variants file is technically possible but not the documented workflow.

Do not use this command without `--verbose`. The cost of skipping it is a 20-50% opaque chunk of the funnel attribution in step 3.

Do not use this command if you cannot afford 10-30+ minutes of run time. For quick iteration on smaller scope, restrict `--region` to a 1 Mb subregion (the `LANCET_TEST_*_REGION_SMALL` env vars).

## Maintenance

The required-flag list in phase 3 must match what step 1's outputs and step 3's inputs assume. If `cli_interface.cpp` adds, removes, or renames a probe-related flag, update this command. The `/audit-probe-pipeline` slash command catches drift.

The version-string parsing in phase 2 assumes the format `<tag>-<branch>-<short10-sha>`. If `src/lancet/version.h.inc` changes that format (e.g., switches to a different separator), update the awk command.
