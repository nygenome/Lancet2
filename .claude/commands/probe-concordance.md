# /probe-concordance — Step 1 of probe tracking: truth concordance

Run `scripts/truth_concordance.py` to compare a truth VCF against a Lancet2 output VCF, classify every truth variant by match level, and produce the inputs that step 2 (`/probe-run`) and step 3 (`/probe-analyze`) consume.

This is one of three composable commands for the probe tracking workflow. The full workflow is `/probe-concordance` → `/probe-run` → `/probe-analyze`. Each command can also stand alone if its inputs already exist on disk. The operational mechanics — what each script does, what flag goes where, how to handle the somatic case — are documented in the `probe-tracking` skill; reach for it when the operational details get unclear.

## When to use

Run this when you have a Lancet2 output VCF and a corresponding truth VCF, and you want to know which truth variants were called correctly, partially, or missed. The output `missed_variants.txt` is the prerequisite for step 2 — without it, you have nothing to probe.

For the germline NA12878 / chr1 fixture, both small (GIAB) and large (Manta) truth VCFs are available; default to `--mode all`. For the somatic HCC1395 / chr4 fixture, no truth VCFs ship with the test data — surface that gap to the user before running and ask whether they have a private truth set.

## Procedure

### Phase 1 — Gather inputs

If the user did not supply explicit paths, use `AskUserQuestion` to gather:

- Truth VCF(s) — at least one of `--truth-small` (GIAB-style SNVs+indels) or `--truth-large` (Manta-style SVs)
- Lancet2 output VCF — `--lancet`
- Reference FASTA — `--ref` (required even for non-CRAM samples; the read-evidence engine needs it)
- Sample BAM/CRAM file(s) — `--samples`
- Output directory — defaults to `notes/probe-debug-<YYYY-MM-DD>/`

Suggest defaults from the `test-data-locations` skill. If the user names "germline", suggest the GIAB+Manta+NA12878 set; if they name "somatic" before truth VCFs are added, surface the gap explicitly.

If the user runs the command bare with no args at all, walk through `AskUserQuestion` to gather every input. Do not silently fail.

### Phase 2 — Confirm output directory

Compute the default output directory: `notes/probe-debug-$(date +%Y-%m-%d)/`. If it already exists, ask the user whether to (a) reuse it (additional outputs land alongside existing ones), (b) timestamp it (`-HHMMSS` suffix), or (c) use a different path. Default to (b) if the existing directory has files in it, (a) if it is empty.

`mkdir -p` the directory before running.

### Phase 3 — Run the script

```bash
pixi run -e hts-tools python3 scripts/truth_concordance.py \
    --truth-small <path>          # if applicable
    --truth-large <path>          # if applicable
    --lancet <path> \
    --ref <path> \
    --samples <path> [<path> ...] \
    --output-dir <outdir> \
    --mode {small,large,all} \
    --workers <N>
```

Choose `--mode` based on which truth VCFs were provided. `--workers` defaults to 16; on machines with many cores, raise it (the bwa-mem2 realignment in tier-2 sensitivity analysis is the slow part). Capture both stdout and stderr; the script writes a console report to stdout and progress to stderr.

The script may take several minutes on a full chromosome with many variants. Stream output to the user as it runs; do not silently wait.

### Phase 4 — Surface what was produced

After the script completes, list the files written to the output directory:

```bash
ls -la <outdir>
```

Highlight the three downstream-relevant outputs:

- `missed_variants.txt` — required input for step 2 (`/probe-run`).
- `concordance_details.txt` — optional enrichment for step 3 (`/probe-analyze`).
- `truth_concordance_report.txt` — rich diagnostic report for human review (sections §1–§4: concordance funnel, sensitivity, specificity).

Read the rich report's §1 funnel summary (concordance level distribution: L0/LD/L1/L2/L3/MISS) and surface the top-line numbers in chat. The user wants to know the missed-variant count before deciding whether to proceed to step 2.

### Phase 5 — Next step

End the response with:

> Step 1 complete. Next: `/probe-run` to collect forensic data for the missed variants. The `--probe-variants` input for that command is `<outdir>/missed_variants.txt`.

Do not auto-invoke step 2; the user often wants to inspect the concordance report first.

## When NOT to use this command

Do not use this command if you only have a Lancet2 VCF and no truth VCF. The whole point is the comparison; without truth, there is nothing for the script to compare against.

Do not use this command on the somatic HCC1395 fixture until truth VCFs are added to the test data. The germline fixture works end-to-end.

Do not use this command for general VCF comparison or Mendelian-error analysis. The script is specifically tuned for the truth-vs-Lancet2 case and emits Lancet2-flavored outputs.

## Maintenance

This command's flag list and defaults track `scripts/truth_concordance.py`. If the script's argparse changes, update this command. The `/audit-probe-pipeline` slash command catches this drift.
