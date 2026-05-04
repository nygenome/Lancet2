---
name: probe-tracking
description: Use whenever running, planning, or troubleshooting Lancet2's probe variant forensic pipeline — the three-step workflow of truth_concordance.py → Lancet2 with --probe-variants/--probe-results → analyze_probe_results.py. Covers operational mechanics: which inputs go where, the flag dance between steps, the data layout, --verbose requirements, and how to derive paths for the somatic case (currently lacks truth VCFs). Use proactively when the user mentions probe tracking, missed variants, truth concordance, lost_at_stage, sensitivity debugging, or specificity debugging. For interpretation of an existing report — "what does this funnel mean, where do I dig in next?" — delegate to the `probe-interpreter` subagent instead.
allowed-tools: Read, Glob, Grep, Bash
---

# Probe-tracking on Lancet2

This skill is the operational playbook for the probe variant forensic pipeline. Its job is to get the workflow run correctly: right inputs, right flags, right output paths, in the right order. The interpretive layer — what to do when stage X dominates the funnel, what code to read for `geno_reads_reassigned`, how to act on the report — lives in the `probe-interpreter` subagent.

The discipline these are split: when the user is *running* the pipeline, this skill applies. When the user is *interpreting* a report that has already been produced, the subagent applies.

## What the pipeline answers

For every truth variant Lancet2 missed, exactly where in the pipeline did it get lost? The architecture has two halves: a C++ data emitter that records raw facts per `(probe_id, window, comp_id, k)` attempt with no attribution logic, and a Python attribution engine (`scripts/analyze_probe_results.py`) that derives `lost_at_stage` from those facts using a 27-level bottom-up cascade. The full conceptual model lives in `docs_dev/subsystems/probe_tracking.md` — read that document any time the workflow's structure feels unclear.

## The three steps

The workflow has three steps, composable in any combination. Step 1's outputs are step 2's and step 3's inputs, but you can run any step alone if its prerequisites already exist on disk.

### Step 1: Truth concordance

`scripts/truth_concordance.py` compares a truth VCF to a Lancet2 output VCF. Outputs that downstream steps depend on are `missed_variants.txt` (consumed by step 2 as `--probe-variants`) and `concordance_details.txt` (optional enrichment for step 3 §1).

```
pixi run -e hts-tools python3 scripts/truth_concordance.py \
    --truth-small <giab_chr1.vcf.gz> \
    --truth-large <manta_chr1.vcf.gz> \
    --lancet <lancet_output.vcf.gz> \
    --ref <reference.fa.gz> \
    --samples <sample.bam_or.cram> [<sample2.bam_or.cram> ...] \
    --output-dir <outdir> \
    --mode {small,large,all} \
    --workers <N>
```

Required: `--lancet`, `--ref`, `--samples`. At least one of `--truth-small` or `--truth-large`. Mode defaults to `all`; workers defaults to 16 (raise for read-evidence pass on large samples; the bwa-mem2 realignment in tier-2 sensitivity is the slow part).

The script invokes the `hts-tools` pixi env, which provides Polars, Rich, pysam, edlib, tqdm, and bwa-mem2.

### Step 2: Lancet2 probe run

Re-run Lancet2 with the missed variants to collect forensic data:

```
./cmake-build-release/Lancet2 pipeline \
    -n <normal.bam_or.cram> -t <tumor.bam_or.cram> \
    -r <reference.fa> -o <out.vcf.gz> \
    --probe-variants <outdir>/missed_variants.txt \
    --probe-results <outdir>/probe_results.tsv \
    --verbose
```

`--probe-variants` and `--probe-results` are bidirectionally required (one without the other is a CLI error). Without either, the entire probe system is inactive — zero overhead, which is the production default.

`--verbose` is strongly recommended even though it is technically optional. Without it, step 3's `--log` flag can't sub-classify the `not_processed` category, which is typically 20-50% of probes. Always include `--verbose` and capture the debug log:

```
./cmake-build-release/Lancet2 pipeline ... --verbose 2> <outdir>/lancet2_debug.log
```

For single-sample germline runs, pass the lone sample as `--normal` (the CLI flag is `required=true`); the case-control mode flag in `pipeline_runner.cpp` is `has_label(CASE) && has_label(CTRL)`, so a control-only run is supported.

### Step 3: Probe analysis

```
pixi run -e hts-tools python3 scripts/analyze_probe_results.py \
    --probe-results <outdir>/probe_results.tsv \
    --missed-variants <outdir>/missed_variants.txt \
    --concordance-details <outdir>/concordance_details.txt \
    --log <outdir>/lancet2_debug.log \
    --output-dir <outdir> \
    --view {scorecard,funnel,survival,breakdown,genotyper,targets,deepdive,windows,all}
```

Required: `--probe-results`, `--missed-variants`. Optional but high-value: `--concordance-details` (enriches §1 scorecard), `--log` (enables the 6 sub-stages of `not_processed`). `--view` defaults to `all`; use individual views during iterative debugging to skip rendering the rest.

Outputs that land in `--output-dir`:

- `probe_analysis_report.txt` — full rich-text report (8 sections)
- `probe_stage_attribution.txt` — TSV: one row per probe with final `lost_at_stage`, type, size
- `probe_survival_matrix.txt` — TSV: one row per `(probe, window, comp, k)` with 6 survival counts

## Output directory convention

Default to `notes/probe-debug-<YYYY-MM-DD>/` for every run. The bundle's protected-paths hook explicitly does not block writes there, the path is gitignored, and dating the directory preserves history across runs. The `data/` directory is the script defaults but pollutes the source tree and overwrites previous runs — don't use it except for one-off ad-hoc work.

When using `notes/probe-debug-<date>/`, pass it as `--output-dir` to all three steps so the missed-variants TSV, the probe results TSV, the debug log, and the analysis outputs all live in one directory. This makes the run self-contained and easy to delete, archive, or share.

## Test-data inputs

The `test-data-locations` skill catalogs the available fixtures; consult it before picking inputs. The relevant pairings for the probe pipeline:

- **Germline (NA12878 / chr1)**: full small + large truth available (`expected_small_variants_giab.chr1.vcf.gz`, `expected_large_variants_manta.chr1.vcf.gz`). Use `--mode all` with both `--truth-small` and `--truth-large`.
- **Somatic (HCC1395 / chr4)**: no truth VCFs currently in the dataset. Step 1 cannot be run as documented for the somatic case until truth VCFs are added. Step 2 and step 3 can be run if a hand-curated `missed_variants.txt` is available, but the standard workflow is blocked here.

When the user asks to debug the somatic case, surface the truth-VCF gap explicitly rather than failing inside step 1. The pipeline can still be exercised on synthetic missed variants for development, but real specificity/sensitivity work needs the truth set.

## Build staleness check (for step 2)

Step 2 takes long enough that running it against a stale binary is expensive. Before invoking Lancet2 with `--probe-variants`, verify the binary is current:

1. Is the working tree clean? Run `git status --porcelain`. Non-empty output (staged, unstaged, OR untracked changes) means you cannot trust the binary's embedded SHA — always rebuild.
2. If clean, parse the SHA from `./cmake-build-release/Lancet2 --version`. The output format is `<tag>-<branch>-<short10-sha>` (set at CMake **configure** time from the lancet_version.h template). Compare the SHA to `git rev-parse --short=10 HEAD`. If they differ, rebuild: `pixi run build-release` (it `depends-on` `configure-release`, which regenerates the version header, so a single command is enough).
3. If the binary doesn't exist at all, build first.

The check is cheap (two git commands) and catches the common "I forgot to rebuild after pulling" mistake before sinking 10-30 minutes into step 2.

## Common operational pitfalls

**Forgetting `--verbose` in step 2.** Step 3's `--log` flag is what unlocks the 6 sub-stages of `not_processed` (ref_all_n, ref_repeat, inactive, low_coverage, no_alt_haplotype, other_variant_called). Without `--verbose`, step 3 still works but the `not_processed` category is opaque — usually the largest single category, so opacity here hurts.

**Mismatched paths between steps.** If step 1 wrote to `data/missed_variants.txt` but step 2 was given `notes/probe-debug-2026-04-28/missed_variants.txt`, you will get a file-not-found error or worse, the wrong file. Always pass `--output-dir` consistently across all three steps and reuse the same directory.

**Wrong `--mode` for the available truth.** `--mode small` with only `--truth-large` is a silent no-op. The script does not error; it produces a report with empty sections. Always pass `--mode` matching what truth files you provided.

**Using single-sample germline with `--tumor` instead of `--normal`.** The CLI rejects this. The germline path uses `--normal <cram>` (control-only mode); the equivalent advanced form is `--sample <path>:ctrl`.

**Probe results from a stale binary getting analyzed by the latest script.** The C++ schema and the Python cascade are tightly coupled. If the binary is older than `analyze_probe_results.py`, the column set may not match. Run `/audit-probe-pipeline` if you suspect drift; in normal operation this is rare because both ship together in the repo.

**Output directory not gitignored.** If you put outputs in a path not covered by `.gitignore`, your `git status` will be polluted. The `notes/probe-debug-*/` convention is gitignored; if you use a different layout, add it to `.gitignore` first.

## Composing the steps

The three commands are independent but commonly chain. The bundle exposes `/probe-concordance`, `/probe-run`, and `/probe-analyze` as slash commands; each command's "Next step" section points at the next one in the chain. There is no monolithic `/probe-debug` command on purpose — the user often wants to run only one or two of them, and forcing all three would obscure what each step contributes.

For a fresh end-to-end run:

```
/probe-concordance   # produces missed_variants.txt + concordance_details.txt
/probe-run           # produces probe_results.tsv + lancet2_debug.log
/probe-analyze       # produces probe_analysis_report.txt and hands off to probe-interpreter
```

For iterative debugging on already-collected data:

```
/probe-analyze --view funnel    # render just one section
/probe-analyze --view targets   # render another section
```

For re-running step 2 after a Lancet2 source change:

```
/probe-run                       # check builds, capture fresh probe_results.tsv
/probe-analyze                   # re-attribute and compare to previous report
```

## When NOT to use this skill

Do not use this skill for interpreting a probe analysis report — that is the `probe-interpreter` subagent's job. The split between operational and interpretive is intentional: this skill is the playbook for getting clean data; the subagent is the consultant for what to do with it.

Do not use this skill for general Lancet2 sensitivity/specificity questions divorced from the probe workflow. If the user asks "why is Lancet2's specificity dropping" without referencing missed variants or probe tracking, that's a `fresh-reviewer` or `assembly-and-calling-expert` question.

Do not use this skill on the somatic fixture (HCC1395 / chr4) for step 1 until truth VCFs are added to the dataset. Surface the gap to the user explicitly.

## Maintenance

The flag tables and command shapes in this skill drift if the upstream scripts change. The `/audit-probe-pipeline` slash command catches that drift by walking `docs_dev/subsystems/probe_tracking.md` against the actual source files; run it after any change to `scripts/truth_concordance.py`, `scripts/analyze_probe_results.py`, `src/lancet/cbdg/probe_*.cpp`, or `src/lancet/core/probe_diagnostics.cpp`.

The "common operational pitfalls" section is the part most likely to grow over time as new failure modes surface. Add to it; do not delete entries from it without justification — pitfalls that occurred once usually recur.

The somatic-truth-gap note will become stale when truth VCFs are added for the somatic case. Remove the gap note then and add the relevant fixture rows to the test-data pairings.
