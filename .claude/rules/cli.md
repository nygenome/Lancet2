---
description: Lancet2 cli/ layer rules — CLI11-based flag surface with AddOpt template, BGZF stream open-flush-close lifecycle that produces tail-readable output, BuildVcfHeader as the canonical FORMAT/INFO field declaration site, case-control mode computed once at startup, cloud-credential pre-flight validation. Load when editing src/lancet/cli/.
paths:
  - "src/lancet/cli/**"
---

# cli/ layer rules

`cli/` is the orchestration entry point: parse flags, validate
inputs, open the output VCF, write the header, configure the worker
pool, hand control to `core::PipelineExecutor`, then close cleanly.
The layer is small (~880 lines) but every line of `pipeline_runner.
cpp::Run` is ordering-sensitive. The CLI surface is also the
strictest part of the API — once a flag is shipped, removing it
breaks every downstream user.

## Required flags: --normal, --reference, --out-vcfgz

`cli_interface.cpp` declares only three CLI options as
`required(true)` via the `AddOpt` template's last parameter:
`-n,--normal`, `-r,--reference`, `-o,--out-vcfgz`. **`-t,--tumor` is
optional** — Lancet2 supports control-only mode (single-sample
germline). The `--sample <path>:<role>` advanced form is the
forward-facing unified surface; `--normal`/`--tumor` are kept as
ergonomic shorthands that map to it.

When adding a flag, use the `AddOpt` template (or `AddFlag` for
booleans) and assign it to the appropriate `GRP_*` group constant
(`GRP_REQUIRED`, `GRP_DATASETS`, `GRP_REGIONS`, `GRP_PARAMETERS`,
`GRP_FLAGS`, `GRP_OPTIONAL`). Don't add new group constants without
understanding how `--help` output groups flags.

## Removing or renaming a flag is a breaking change

Every flag in `cli_interface.cpp` is part of the public API. Renaming
`--num-threads` to `--threads`, even if the docs are updated,
breaks every wrapper script and every CI pipeline that runs Lancet2.
Flag removal/renaming requires:

1. Explicit user approval in the change request
2. A deprecation cycle: keep the old name as an alias for one release
3. CHANGELOG entry under "Breaking changes"
4. `/audit-bundle` pairing 3 (CLI flag surface) catches drift between
   `cli_interface.cpp` and the `e2e_pipeline_test.sh` script's flag verification

Use `vcf-validator`-style review for flag changes — they're
analogous to schema changes in their downstream impact.

## case-control mode is computed ONCE at startup, immutable thereafter

`mIsCaseCtrlMode = has_label(CASE) && has_label(CTRL)` in
`ValidateAndPopulateParams` runs **before** any worker starts.
Subsequent code (VCF header builder, read collector, active region
detector) reads the flag without re-deriving it. **Do not** compute
case/control state from sample lists at multiple sites — divergent
derivations are a classic source of drift bugs.

The flag specifically gates SHARED/CTRL/CASE INFO header emission
in `vcf_header_builder.cpp`. Single-sample (CTRL-only) runs emit no
case-vs-control INFO fields; mixed runs emit all three.

## VCF header is committed BEFORE workers start

The lifecycle in `PipelineRunner::Run`:

```
output_vcf.Open(...)                  // BgzfOstream::Open
output_vcf << BuildVcfHeader(*params) // header text
output_vcf.flush()                    // header committed to disk
                                      // (--- workers may now start ---)
window_builder.SortInputRegions()     // deterministic genomic ordering
PipelineExecutor::Execute(output_vcf) // workers + flush loop
output_vcf.Close()                    // BGZF EOF marker block
```

Two ordering invariants:

1. **`output_vcf.flush()` after the header** — without this, the
   header lives in BGZF's internal buffer; if a worker crashes
   before the buffer fills, the on-disk file has zero content. The
   explicit flush forces a BGZF block boundary and makes the file
   tail-readable from the moment workers start.

2. **`output_vcf.Close()` after `Execute` returns** — finalizes the
   BGZF stream's EOF marker. A path that returns from `Run` without
   calling `Close()` produces a silently truncated file that
   `bcftools view` reports as "premature EOF". The `std::exit(EXIT_
   SUCCESS)` at the end of `Run` runs after `Close`; preserve this.

## SortInputRegions before batch emission, not after

`window_builder.SortInputRegions()` runs **before**
`PipelineExecutor::Execute`. Sorting input → sorted output windows →
workers see windows in genomic order → the chunked sorted flush in
`core::VariantStore` works correctly. Sort after batching and the
out-of-order completion logic in `FlushCompletedVariants` fights
against pre-shuffled inputs. The window builder's `SortInputRegions`
is the single source of input ordering.

## BuildVcfHeader is the canonical FORMAT/INFO declaration site

`vcf_header_builder.cpp::FORMAT_STR_HEADER` is a constexpr raw string
holding the entire VCF v4.5 header template with named placeholders
(`{RUN_TIMESTAMP}`, `{FULL_VERSION_TAG}`, `{FULL_COMMAND_USED}`,
`{REFERENCE_PATH}`, `{CONTIG_HDR_LINES}`, `{CONDITIONAL_INFO_LINES}`,
`{ANNOTATION_INFO_LINES}`). **The unconditional FORMAT and INFO field
declarations live here as one block of `##FORMAT`/`##INFO` lines.**
The case/control-specific INFO fields (SHARED, CTRL, CASE) are
declared in a separate constexpr raw string `CASE_CTRL_INFO_HDR_LINES`
in the same file, gated by `mIsCaseCtrlMode` at header-emission
time. Together those two raw strings are the canonical declaration
sites for VCF schema fields.

The list is the single source of truth for which fields the VCF
emits. Code in `caller/variant_call.cpp` (and downstream) populates
these fields; if the field is in `FORMAT_STR_HEADER` but no code
populates it, every record gets `.` for it. If code emits a field
not declared here, every record fails VCF v4.5 validation.

When adding a FORMAT or INFO field:

1. Add the `##FORMAT=<...>` or `##INFO=<...>` line to the appropriate
   block in `vcf_header_builder.cpp`.
2. Update the value-emission code in `caller/variant_call.cpp` (or
   wherever the field is computed).
3. Invoke `vcf-validator` to verify schema correctness — it's the
   canonical home for the project's schema invariants.
4. The `/audit-vcf-schema` slash command catches drift between
   header declarations and emission code.

The CASE_CTRL_INFO_HDR_LINES block (SHARED/CTRL/CASE) is gated by
`mIsCaseCtrlMode`. Adding a new INFO field that's also case-control-
specific should follow the same gating pattern.

## CONTIGS_BUFFER_SIZE = 524'288 is sized for human references

The `contig_hdr_lines.reserve(524'288)` pre-allocation is sized for
typical human reference contig counts (a few hundred). Mouse and
plant references with thousands of contigs may benefit from a
larger buffer; the reservation is a hint, so wrong size is just a
realloc, not a correctness issue.

## Cloud credentials validated upfront

The cloud pre-flight in `OpenOutputVcf` calls
`hts::ValidateCloudAccess` for `gs://` and `s3://` URIs **before**
any work starts. Without this, libcurl's 5MB multipart threshold
defers the HTTP handshake to `BgzfOstream::Close()` — a 40-hour
pipeline that silently authenticates only at the end is the worst
failure mode possible. Adding new cloud-output destinations follows
the same pattern: validate at startup, fail fast.

## Pipeline exit is via std::exit, not return

`PipelineRunner::Run` is `[[noreturn]]` and ends in `std::exit(EXIT_
SUCCESS)` (and `EXIT_FAILURE` on cloud auth, missing file, invalid
flags, etc.). Do not refactor to return a status code from `Run` —
the exit ordering ensures `output_vcf.Close()` and
`MergePerWorkerGraphShards()` complete before destructors run on the
stack. Changing this risks UAF on the BGZF stream during stack
unwind.

## Probe tracking setup happens after VCF setup, not before

`SetupProbeTracking` runs after `OpenOutputVcf` because probe-index
construction is expensive (~minutes for chrosome-scale fixtures)
and the user wants to see the VCF file appear before the probe
indexing begins — visible progress matters. The function is a no-op
when `--probe-variants` is unset (zero overhead in production).

## SetupPerWorkerGraphShards: `.tar.gz` extension is enforced

The `--out-graphs-tgz` value must end in `.tar.gz` — the merger uses
this extension for the final archive. The check is a hard exit
because a typo (`.tar` or `.tgz`) silently writes nothing useful. If
the path doesn't pass, exit with `LOG_CRITICAL` rather than
"helpfully" appending `.tar.gz` (the user wanted something specific;
guess wrong and they overwrite the wrong file).
