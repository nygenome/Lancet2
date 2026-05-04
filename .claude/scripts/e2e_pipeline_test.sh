#!/usr/bin/env bash

# Fixture presence check
required_vars=(
    LANCET_TEST_GERMLINE_CRAM
    LANCET_TEST_GERMLINE_REFERENCE LANCET_TEST_GERMLINE_REGION
    LANCET_TEST_SOMATIC_TUMOR LANCET_TEST_SOMATIC_NORMAL
    LANCET_TEST_SOMATIC_REFERENCE LANCET_TEST_SOMATIC_REGION
)
for var in "${required_vars[@]}"; do
    if [ -z "${!var:-}" ]; then
        echo "X required env var $var is not set" >&2
        echo "  Source the project's test fixture script before running this command." >&2
        exit 1
    fi
done
# Lancet2 /e2e-pipeline-test end-to-end smoke test
#
# Run the full pipeline against both test profiles back-to-back:
#   germline (NA12878 / chr1, single-sample) then somatic (HCC1395 +
#   HCC1395BL / chr4, tumor/normal). Confirms exit-0 and a reasonable
#   variant count on each; truth-set comparison is out of scope.
#
# This script is invoked by the /e2e-pipeline-test slash command (.claude/commands/e2e-pipeline-test.md)
# but is also runnable standalone as `bash .claude/scripts/e2e-pipeline-test.sh`.
#
# Required env vars (typically loaded from .claude/settings.local.json):
#   LANCET_TEST_GERMLINE_CRAM, _REFERENCE, _REGION
#   LANCET_TEST_SOMATIC_TUMOR, _NORMAL, _REFERENCE, _REGION
#
# Stages skip cleanly if their env vars are unset or data files are
# absent — a missing-data run reports `skipped`, not `fail`.

set +e

# ── Artifact manifest ─────────────────────────────────────────────────────
# Every path the script writes to /tmp during this run gets recorded here.
# Reasons (a) so the developer can clean up explicitly rather than waiting
# for the OS's /tmp reaper to silently delete a multi-GB VCF after a few
# days, (b) so the post-run AskUserQuestion step (handled by Claude per
# `.claude/commands/e2e-pipeline-test.md`) has a single source of truth
# for which files are eligible for cleanup, and (c) so a developer who
# wants to keep the artifacts for further analysis (truth-set comparison,
# pprof attribution, etc.) gets an explicit list to copy-paste.
LANCET_E2E_RUN_TS=$(date +%s)
LANCET_E2E_MANIFEST="/tmp/lancet2_e2e_artifacts_${LANCET_E2E_RUN_TS}.list"
: > "$LANCET_E2E_MANIFEST"
echo "Artifact manifest: $LANCET_E2E_MANIFEST (paths appended as the run progresses)"

record_artifact() {
  # Args: one path per call. Filters out empty strings and already-recorded
  # entries so the manifest stays a clean deduplicated list.
  local path="$1"
  [ -z "$path" ] && return 0
  if ! grep -Fxq -- "$path" "$LANCET_E2E_MANIFEST" 2>/dev/null; then
    printf '%s\n' "$path" >> "$LANCET_E2E_MANIFEST"
  fi
}

# ── Thread budget ─────────────────────────────────────────────────────────
# Use 2/3 of the machine's logical cores rather than all of them. Saturating
# every core makes the workstation unresponsive while the test runs and gives
# Lancet2 nothing in return — the per-window dispatcher is bottlenecked on
# I/O at the high end of the thread count, not CPU. 2/3 is a heuristic: it
# leaves headroom for system tasks, the developer's editor, and the bcftools
# subprocess at the end of each stage. Floor at 1 for tiny machines (single-
# core CI runners).
LANCET_E2E_TOTAL_CORES=$(nproc)
LANCET_E2E_NUM_THREADS=$(( (LANCET_E2E_TOTAL_CORES * 2) / 3 ))
[ "$LANCET_E2E_NUM_THREADS" -lt 1 ] && LANCET_E2E_NUM_THREADS=1
echo "Thread budget: ${LANCET_E2E_NUM_THREADS} of ${LANCET_E2E_TOTAL_CORES} cores (2/3)"

# ── Validate prerequisites ────────────────────────────────────────────────

if [ ! -x cmake-build-release/Lancet2 ]; then
  echo "Release binary missing; building first ..."
  pixi run --quiet build-release
  if [ $? -ne 0 ]; then
    echo "❌ pixi run build-release failed; aborting /e2e-pipeline-test"
    exit 1
  fi
fi

# ── Verify the CLI shape we expect ────────────────────────────────────────

help_out=$(./cmake-build-release/Lancet2 pipeline --help 2>&1)
for required_flag in "--tumor" "--normal" "--reference" "--region" "--out-vcfgz"; do
  if ! echo "$help_out" | grep -q -- "$required_flag"; then
    echo "❌ Lancet2 pipeline --help does not mention $required_flag."
    echo "   The CLI shape this command assumes has changed; update /e2e-pipeline-test or"
    echo "   .claude/commands/e2e-pipeline-test.md to match the current flags."
    exit 1
  fi
done

# ── Validate env vars and data files ──────────────────────────────────────

check_profile() {
  # Args: profile_name, var_name1, var_name2, ...
  # Exits with skip code 0 (and a message) if any is missing or any path file is absent.
  local profile="$1"; shift
  local missing=""
  for var in "$@"; do
    local val="${!var}"
    if [ -z "$val" ]; then
      missing="$missing $var"
    elif [[ "$var" == *_REGION* ]]; then
      :  # regions are strings, not files
    elif [ ! -f "$val" ]; then
      echo "❌ $profile: $var=$val does not exist"
      return 2
    fi
  done
  if [ -n "$missing" ]; then
    echo "ℹ $profile skipped: env vars not set:$missing"
    echo "  Populate .claude/settings.local.json from settings.local.json.example."
    return 1
  fi
  return 0
}

# ── Run the germline stage ────────────────────────────────────────────────

run_germline() {
  check_profile "germline" \
    LANCET_TEST_GERMLINE_CRAM \
    LANCET_TEST_GERMLINE_REFERENCE \
    LANCET_TEST_GERMLINE_REGION
  local rc=$?
  if [ $rc -ne 0 ]; then return $rc; fi

  echo "─── /e2e-pipeline-test germline (chr1, single-sample NA12878) ─────────────"
  echo "Sample:    $LANCET_TEST_GERMLINE_CRAM"
  echo "Reference: $LANCET_TEST_GERMLINE_REFERENCE"
  echo "Region:    $LANCET_TEST_GERMLINE_REGION"

  local out_vcf stage_log
  out_vcf=/tmp/lancet2_e2e_germline_${LANCET_E2E_RUN_TS}.vcf.gz
  stage_log=/tmp/lancet2_e2e_germline_${LANCET_E2E_RUN_TS}.log
  record_artifact "$out_vcf"
  record_artifact "${out_vcf}.tbi"
  record_artifact "$stage_log"
  echo "Output:    $out_vcf"
  echo "Log:       $stage_log"
  echo ""

  # `tee` streams Lancet2's progress lines (per-window EtaTimer ticks,
  # SKIPPED/FOUND breakdown, PeakRSS) to stdout in real time AND mirrors
  # the same content into $stage_log for postmortem inspection. The
  # previous `| tail -10` buffered everything until the stage ended,
  # losing all live visibility.
  local start end elapsed
  start=$(date +%s)
  ./cmake-build-release/Lancet2 pipeline \
    --normal "$LANCET_TEST_GERMLINE_CRAM" \
    --reference "$LANCET_TEST_GERMLINE_REFERENCE" \
    --region "$LANCET_TEST_GERMLINE_REGION" \
    --num-threads "$LANCET_E2E_NUM_THREADS" \
    --out-vcfgz "$out_vcf" 2>&1 | tee "$stage_log"
  rc=${PIPESTATUS[0]}
  end=$(date +%s); elapsed=$((end - start))

  if [ $rc -ne 0 ]; then
    echo "❌ germline pipeline exited with code $rc (wall: ${elapsed}s)"
    GERMLINE_RC=$rc; GERMLINE_ELAPSED=$elapsed; GERMLINE_TOTAL=0
    return $rc
  fi

  # Lancet2 itself does not emit a PASS/FAIL FILTER value — the FILTER
  # column on every output record is "." (untestable) by design. PASS-vs-
  # FAIL filtering is the responsibility of the downstream Python ML
  # module that consumes Lancet2's VCFs. We therefore report only the
  # total raw-variant count, not a PASS count.
  local n_total=0
  if [ -f "$out_vcf" ]; then
    n_total=$(pixi run --quiet -e hts-tools bcftools view "$out_vcf" 2>/dev/null | grep -v '^#' | wc -l)
  fi
  echo "✓ germline complete in ${elapsed}s — ${n_total} raw variants"
  GERMLINE_RC=0; GERMLINE_ELAPSED=$elapsed; GERMLINE_TOTAL=$n_total
}

# ── Run the somatic stage ─────────────────────────────────────────────────

run_somatic() {
  check_profile "somatic" \
    LANCET_TEST_SOMATIC_TUMOR \
    LANCET_TEST_SOMATIC_NORMAL \
    LANCET_TEST_SOMATIC_REFERENCE \
    LANCET_TEST_SOMATIC_REGION
  local rc=$?
  if [ $rc -ne 0 ]; then return $rc; fi

  echo ""
  echo "─── /e2e-pipeline-test somatic (chr4, HCC1395 tumor / HCC1395BL normal) ────"
  echo "Tumor:     $LANCET_TEST_SOMATIC_TUMOR"
  echo "Normal:    $LANCET_TEST_SOMATIC_NORMAL"
  echo "Reference: $LANCET_TEST_SOMATIC_REFERENCE"
  echo "Region:    $LANCET_TEST_SOMATIC_REGION"

  local out_vcf stage_log
  out_vcf=/tmp/lancet2_e2e_somatic_${LANCET_E2E_RUN_TS}.vcf.gz
  stage_log=/tmp/lancet2_e2e_somatic_${LANCET_E2E_RUN_TS}.log
  record_artifact "$out_vcf"
  record_artifact "${out_vcf}.tbi"
  record_artifact "$stage_log"
  echo "Output:    $out_vcf"
  echo "Log:       $stage_log"
  echo ""

  # See germline path for the rationale on `tee`-vs-`tail` for live
  # progress capture.
  local start end elapsed
  start=$(date +%s)
  ./cmake-build-release/Lancet2 pipeline \
    --tumor "$LANCET_TEST_SOMATIC_TUMOR" \
    --normal "$LANCET_TEST_SOMATIC_NORMAL" \
    --reference "$LANCET_TEST_SOMATIC_REFERENCE" \
    --region "$LANCET_TEST_SOMATIC_REGION" \
    --num-threads "$LANCET_E2E_NUM_THREADS" \
    --out-vcfgz "$out_vcf" 2>&1 | tee "$stage_log"
  rc=${PIPESTATUS[0]}
  end=$(date +%s); elapsed=$((end - start))

  if [ $rc -ne 0 ]; then
    echo "❌ somatic pipeline exited with code $rc (wall: ${elapsed}s)"
    SOMATIC_RC=$rc; SOMATIC_ELAPSED=$elapsed; SOMATIC_TOTAL=0
    return $rc
  fi

  # See germline path for the rationale: Lancet2 emits FILTER="." on
  # every record; PASS filtering is downstream-only.
  local n_total=0
  if [ -f "$out_vcf" ]; then
    n_total=$(pixi run --quiet -e hts-tools bcftools view "$out_vcf" 2>/dev/null | grep -v '^#' | wc -l)
  fi
  echo "✓ somatic complete in ${elapsed}s — ${n_total} raw variants"
  SOMATIC_RC=0; SOMATIC_ELAPSED=$elapsed; SOMATIC_TOTAL=$n_total
}

# ── Drive both stages and summarize ───────────────────────────────────────

GERMLINE_RC=-1; SOMATIC_RC=-1
run_germline
run_somatic

echo ""
echo "─── /e2e-pipeline-test summary ─────────────────────────────────────────────"
printf "  %-10s %-12s %12s %12s\n" "stage" "result" "wall(s)" "raw_variants"
case $GERMLINE_RC in
  -1) printf "  %-10s %-12s %12s %12s\n" "germline" "skipped" "-" "-" ;;
   0) printf "  %-10s %-12s %12s %12s\n" "germline" "ok" "$GERMLINE_ELAPSED" "$GERMLINE_TOTAL" ;;
   *) printf "  %-10s %-12s %12s %12s\n" "germline" "fail($GERMLINE_RC)" "$GERMLINE_ELAPSED" "$GERMLINE_TOTAL" ;;
esac
case $SOMATIC_RC in
  -1) printf "  %-10s %-12s %12s %12s\n" "somatic" "skipped" "-" "-" ;;
   0) printf "  %-10s %-12s %12s %12s\n" "somatic" "ok" "$SOMATIC_ELAPSED" "$SOMATIC_TOTAL" ;;
   *) printf "  %-10s %-12s %12s %12s\n" "somatic" "fail($SOMATIC_RC)" "$SOMATIC_ELAPSED" "$SOMATIC_TOTAL" ;;
esac
echo "──────────────────────────────────────────────────────────────"

# Surface the artifact manifest in a clearly-marked block so the post-run
# AskUserQuestion step (driven by the Claude command spec at
# `.claude/commands/e2e-pipeline-test.md`) has a stable anchor to parse.
# The literal "ARTIFACT MANIFEST:" prefix is matched by that step — do
# not change without updating the command spec.
echo ""
echo "─── /e2e-pipeline-test artifacts ───────────────────────────────"
echo "ARTIFACT MANIFEST: $LANCET_E2E_MANIFEST"
if [ -s "$LANCET_E2E_MANIFEST" ]; then
  while IFS= read -r path; do
    if [ -e "$path" ]; then
      printf "  %-7s  %s\n" "exists" "$path"
    else
      printf "  %-7s  %s\n" "missing" "$path"
    fi
  done < "$LANCET_E2E_MANIFEST"
else
  echo "  (no artifacts recorded)"
fi
echo "──────────────────────────────────────────────────────────────"
echo ""
echo "These files are NOT auto-deleted. The Claude post-run hook will ask"
echo "via AskUserQuestion whether to clean them up. Default answer: yes."
echo "Reply 'keep' if a follow-up analysis (truth-set comparison, profile"
echo "attribution, debug session) needs them."

if [ $GERMLINE_RC -gt 0 ] || [ $SOMATIC_RC -gt 0 ]; then
  exit 1
fi
exit 0
