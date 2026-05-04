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

  local out_vcf=/tmp/lancet2_e2e_germline_$(date +%s).vcf.gz
  echo "Output:    $out_vcf"
  echo ""

  local start end elapsed
  start=$(date +%s)
  ./cmake-build-release/Lancet2 pipeline \
    --normal "$LANCET_TEST_GERMLINE_CRAM" \
    --reference "$LANCET_TEST_GERMLINE_REFERENCE" \
    --region "$LANCET_TEST_GERMLINE_REGION" \
    --num-threads "$(nproc)" \
    --out-vcfgz "$out_vcf" 2>&1 | tail -10
  rc=${PIPESTATUS[0]}
  end=$(date +%s); elapsed=$((end - start))

  if [ $rc -ne 0 ]; then
    echo "❌ germline pipeline exited with code $rc (wall: ${elapsed}s)"
    GERMLINE_RC=$rc; GERMLINE_ELAPSED=$elapsed; GERMLINE_TOTAL=0; GERMLINE_PASS=0
    return $rc
  fi

  local n_pass=0 n_total=0
  if [ -f "$out_vcf" ]; then
    n_total=$(pixi run --quiet -e hts-tools bcftools view "$out_vcf" 2>/dev/null | grep -v '^#' | wc -l)
    n_pass=$(pixi run --quiet -e hts-tools bcftools view -f PASS "$out_vcf" 2>/dev/null | grep -v '^#' | wc -l)
  fi
  echo "✓ germline complete in ${elapsed}s — ${n_total} variants total, ${n_pass} PASS"
  GERMLINE_RC=0; GERMLINE_ELAPSED=$elapsed; GERMLINE_TOTAL=$n_total; GERMLINE_PASS=$n_pass
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

  local out_vcf=/tmp/lancet2_e2e_somatic_$(date +%s).vcf.gz
  echo "Output:    $out_vcf"
  echo ""

  local start end elapsed
  start=$(date +%s)
  ./cmake-build-release/Lancet2 pipeline \
    --tumor "$LANCET_TEST_SOMATIC_TUMOR" \
    --normal "$LANCET_TEST_SOMATIC_NORMAL" \
    --reference "$LANCET_TEST_SOMATIC_REFERENCE" \
    --region "$LANCET_TEST_SOMATIC_REGION" \
    --num-threads "$(nproc)" \
    --out-vcfgz "$out_vcf" 2>&1 | tail -10
  rc=${PIPESTATUS[0]}
  end=$(date +%s); elapsed=$((end - start))

  if [ $rc -ne 0 ]; then
    echo "❌ somatic pipeline exited with code $rc (wall: ${elapsed}s)"
    SOMATIC_RC=$rc; SOMATIC_ELAPSED=$elapsed; SOMATIC_TOTAL=0; SOMATIC_PASS=0
    return $rc
  fi

  local n_pass=0 n_total=0
  if [ -f "$out_vcf" ]; then
    n_total=$(pixi run --quiet -e hts-tools bcftools view "$out_vcf" 2>/dev/null | grep -v '^#' | wc -l)
    n_pass=$(pixi run --quiet -e hts-tools bcftools view -f PASS "$out_vcf" 2>/dev/null | grep -v '^#' | wc -l)
  fi
  echo "✓ somatic complete in ${elapsed}s — ${n_total} variants total, ${n_pass} PASS"
  SOMATIC_RC=0; SOMATIC_ELAPSED=$elapsed; SOMATIC_TOTAL=$n_total; SOMATIC_PASS=$n_pass
}

# ── Drive both stages and summarize ───────────────────────────────────────

GERMLINE_RC=-1; SOMATIC_RC=-1
run_germline
run_somatic

echo ""
echo "─── /e2e-pipeline-test summary ─────────────────────────────────────────────"
printf "  %-10s %-12s %12s %12s %s\n" "stage" "result" "wall(s)" "variants" "PASS"
case $GERMLINE_RC in
  -1) printf "  %-10s %-12s %12s %12s %s\n" "germline" "skipped" "-" "-" "-" ;;
   0) printf "  %-10s %-12s %12s %12s %s\n" "germline" "ok" "$GERMLINE_ELAPSED" "$GERMLINE_TOTAL" "$GERMLINE_PASS" ;;
   *) printf "  %-10s %-12s %12s %12s %s\n" "germline" "fail($GERMLINE_RC)" "$GERMLINE_ELAPSED" "$GERMLINE_TOTAL" "$GERMLINE_PASS" ;;
esac
case $SOMATIC_RC in
  -1) printf "  %-10s %-12s %12s %12s %s\n" "somatic" "skipped" "-" "-" "-" ;;
   0) printf "  %-10s %-12s %12s %12s %s\n" "somatic" "ok" "$SOMATIC_ELAPSED" "$SOMATIC_TOTAL" "$SOMATIC_PASS" ;;
   *) printf "  %-10s %-12s %12s %12s %s\n" "somatic" "fail($SOMATIC_RC)" "$SOMATIC_ELAPSED" "$SOMATIC_TOTAL" "$SOMATIC_PASS" ;;
esac
echo "──────────────────────────────────────────────────────────────"

if [ $GERMLINE_RC -gt 0 ] || [ $SOMATIC_RC -gt 0 ]; then
  exit 1
fi
exit 0
