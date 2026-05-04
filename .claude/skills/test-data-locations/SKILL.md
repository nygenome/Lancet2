---
name: test-data-locations
description: Use whenever the user asks "which test data should I use", "what's in the dataset", "which BAM has X", or when planning a test/benchmark/sanitizer run and you need to pick the right fixture. Documents every file at gs://lancet2-test-datasets/test_harness_data/, the env vars that map to them, and which fixture fits which workflow (somatic, germline, truth comparison, T2T, longdust calibration). Read-only reference.
allowed-tools: Read, Glob, Grep
---

# Lancet2 test-data reference

This skill is the source of truth for what test data Lancet2 has and which fixture fits which job. The dataset lives at `gs://lancet2-test-datasets/test_harness_data/` and is downloaded by `bash data/download_test_data.sh` into `data/`. Per-developer paths are configured in `.claude/settings.local.json` (gitignored); the example file shows the default layout.

## Datasets at a glance

The dataset has five logical groups:

| Purpose | Sample(s) | Reference | Region |
|:--------|:----------|:----------|:-------|
| Germline (single-sample) | NA12878 (CRAM) | GRCh38 | chr1 |
| Somatic (tumor/normal) | HCC1395 / HCC1395BL (BAM) | GRCh38 | chr4 |
| Truth comparison (small variants) | GIAB | GRCh38 | chr1 |
| Truth comparison (large variants) | Manta | GRCh38 | chr1 |
| T2T reference (optional) | n/a (reference only) | CHM13v2.0 | full |
| Longdust calibration | (calibration data) | n/a | n/a |

## Decision tree: which fixture for which job

**"I want to run the pipeline end-to-end on a small representative input."**
- Tumor/normal somatic mode → `LANCET_TEST_SOMATIC_TUMOR` + `LANCET_TEST_SOMATIC_NORMAL`, region `LANCET_TEST_SOMATIC_REGION_SMALL` (a 1 Mb subregion of chr4) for fast iteration; `LANCET_TEST_SOMATIC_REGION` (full chr4) for `/e2e-pipeline-test`.
- Germline single-sample mode → `LANCET_TEST_GERMLINE_CRAM` only, region `LANCET_TEST_GERMLINE_REGION_SMALL` (1 Mb of chr1) for fast iteration; `LANCET_TEST_GERMLINE_REGION` (full chr1) for `/e2e-pipeline-test`. Note that single-sample germline runs by passing the lone CRAM as `--normal` (the CLI's `--normal` flag is `required=true`); the case-control mode flag in `pipeline_runner.cpp` only fires when both case and control samples are present.

**"I'm validating a correctness change against truth and need a high-confidence callset."**
- Small variants (SNVs and small indels) → use the GIAB truth VCF for chr1: `expected_small_variants_giab.chr1.vcf.gz`. This pairs naturally with the germline NA12878 fixture above.
- Large variants (SVs) → use the Manta truth VCF for chr1: `expected_large_variants_manta.chr1.vcf.gz`. Same pairing.
- Truth comparison is a separate harness from `/e2e-pipeline-test` (which only checks exit-0 and variant count). The bundle exposes the probe variant forensic pipeline as the truth-comparison workflow: `/probe-concordance` produces the missed-variants TSV, `/probe-run` collects forensic data, `/probe-analyze` attributes each missed variant to a pipeline stage. The `probe-tracking` skill is the operational playbook; the `probe-interpreter` subagent does the interpretive work. For ad-hoc comparison without the forensic pipeline, `bcftools isec` works fine; the probe pipeline is for "why was this variant missed" rather than just "was it missed."
- Somatic case (HCC1395 / chr4) currently has no truth VCFs in the dataset. Step 1 of the probe pipeline cannot be run for the somatic case until truth VCFs are added. Step 2 and step 3 work fine if the user has a hand-curated `missed_variants.txt`.

**"I'm running sanitizers and need a fast-iterating fixture."**
- Use the `_REGION_SMALL` variants (1 Mb subregions). Sanitizer builds run ~3-5x slower than Release, so the small region keeps the iteration loop tight.
- Both somatic and germline have small-region variants; pick somatic if your code touches the tumor/normal-specific paths (combined_scorer, sample_mask, FORMAT field semantics for SHARED/CTRL/CASE), germline otherwise.

**"I need to profile and the data should be representative."**
- Use the full-region variants (`_REGION` not `_REGION_SMALL`). Profiles on the 1 Mb subregions miss whole-pipeline characteristics (window batching, threading saturation, allocation patterns at scale). The somatic full-chr4 run is the standard profiling fixture; the germline full-chr1 run is a good second.

**"I need T2T-CHM13 for stratified validation."**
- `LANCET_TEST_REFERENCE_T2T` points to `chm13v2.0.fa.gz` if present. This is for cross-reference work — comparing how Lancet2 calls the same biology against GRCh38 versus CHM13. Not required for normal development.

**"I need longdust calibration data."**
- Longdust calibration data ships in the dataset. It is consumed by tests/benchmarks that exercise the polar-features path. If your work is in that area, the data is already on disk after `download_test_data.sh`; the test/benchmark itself knows the path.

## Env-var reference

The skills and the `/e2e-pipeline-test` command consume these env vars. Copy `.claude/settings.local.json.example` to `.claude/settings.local.json` and adjust paths if your checkout location differs from `data/`.

```
LANCET_TEST_GERMLINE_CRAM         data/NA12878.final.cram
LANCET_TEST_GERMLINE_REFERENCE    data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
LANCET_TEST_GERMLINE_REGION       chr1
LANCET_TEST_GERMLINE_REGION_SMALL chr1:1000000-2000000

LANCET_TEST_SOMATIC_TUMOR         data/chr4_with_pairs.HCC1395_SAMN10102573_SRR7890893.bam
LANCET_TEST_SOMATIC_NORMAL        data/chr4_with_pairs.HCC1395BL_SAMN10102574_SRR7890943.bam
LANCET_TEST_SOMATIC_REFERENCE     data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
LANCET_TEST_SOMATIC_REGION        chr4
LANCET_TEST_SOMATIC_REGION_SMALL  chr4:60000000-61000000

LANCET_TEST_REFERENCE_T2T         data/chm13v2.0.fa.gz   (optional)
```

Skills that consume these vars check for presence and skip cleanly if unset rather than failing in an opaque way. If you find a skill failing because a var is unset, copy the example file and re-run.

## What's in the dataset and what isn't

The dataset includes BAM/CRAM reads, the GRCh38 reference and FASTA index, the optional CHM13 reference, and pre-computed truth VCFs for chr1 (GIAB + Manta).

The dataset does NOT include other regions of NA12878, the rest of the HCC1395 BAMs (the project has subset to chr4 with pairs for fast iteration), arbitrary chromosome 22 fixtures (ignore any older docs that reference `chr22:25M-26M` — the standard fixtures are chr1 germline and chr4 somatic), or any custom truth sets. If your workflow needs data not in the dataset, that is a discussion with the maintainer rather than a bundle problem.

## Maintenance

When the dataset changes upstream (new files added, files reorganized at the GCS bucket, regions changed, new truth sources added), update three places in coordinated fashion: this skill (the decision tree and the env-var table), `.claude/settings.local.json.example` (the env vars themselves), and `AGENTS.md` (the brief paragraph about test data). All three should agree; drift between them produces confusing failures the first time someone reaches for a fixture that was renamed.

The decision tree is the most useful section to keep current. If you find yourself asking "which fixture should I use" and the answer is not obvious from the table, add a row.
