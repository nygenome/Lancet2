---
description: Compose a Lancet2-conformant commit message from the staged diff. Reads .claude/commit-style.json (grounded in .chglog/config.yml and docs_dev/style/cpp_style.md § Git commit messages) for rules, drafts subject and body in the project's two-section pattern, asks for confirmation, then commits.
allowed-tools: Bash, Read
---

# /commit — compose a conformant commit message

Read the staged diff, derive the conventional-commit type from what changed, draft a subject line and (if the change is substantial) a two-section body in the shape `docs_dev/style/cpp_style.md` § Git commit messages prescribes, present the proposal, and on confirmation run `git commit`. The validate_commit_message hook will independently double-check the message at the moment of commit.

## When to use this command

Use this command instead of writing `git commit -m "..."` by hand. It gets the type, subject length, body shape, and bullet format right on the first try in nearly every case, which avoids the round-trip of writing a message, having the validator block it, and rewriting.

Do not use this command for merge commits, revert commits, or fixup/squash commits; those have their own conventions and the validator skips them anyway.

## Procedure

First, inspect what is staged. Run `git status --short` and `git diff --cached --stat` to see the file inventory. If nothing is staged, stop and tell the user to `git add` first; do not assume which files should be staged. If a mix of staged and unstaged changes exists, ask the user whether they intended to stage everything or only what is currently in the index.

Second, read the staged diff. Run `git diff --cached` and read the output. For diffs over five hundred lines, fall back to `git diff --cached --stat` plus targeted reads of the most-changed files. The goal is to understand what changed semantically and to know which files should appear in the file-list bullets.

Third, load the style rules. Read `.claude/commit-style.json` to confirm the allowed types, the trivial-commit patterns, and the body-shape expectations. The rules there are derived directly from `.chglog/config.yml` and `docs_dev/style/cpp_style.md` § Git commit messages; treat them as authoritative.

Fourth, derive the type. The Lancet2 chglog config filters on exactly four types, and each maps to a section in the generated CHANGELOG.md:

- `feat` → New Features. Use when the diff adds a user-visible capability: a new CLI flag, a new VCF FORMAT or INFO or FILTER field, a new subcommand, a new pipeline option.
- `fix` → Bug Fixes. Use when the diff corrects incorrect behavior in existing code. Includes correctness fixes, parser fixes, off-by-one fixes, threading-correctness fixes.
- `perf` → Performance Improvements. Use when the diff is a measurable performance change with no behavior change. Must be backed by a benchmark or wall-clock measurement; without measurement evidence, classify it as `chore` instead.
- `chore` → Refactoring. The catch-all. Use for refactors, internal cleanup, build changes, dependency bumps, test additions, documentation updates, CI changes, and anything else that should appear in the changelog but does not fit the three above.

The Lancet2 config does NOT use `refactor`, `docs`, `test`, `build`, `ci`, `style`, or `revert`. A message starting with any of those prefixes parses but is silently dropped by chglog. Map all of them to `chore`. The validator will reject them too.

Prefer the more specific type when more than one fits, in this order: `feat` > `fix` > `perf` > `chore`. A diff that adds a new feature and incidentally fixes a bug is still a `feat` (mention the bug fix in the body's context paragraph).

Fifth, do NOT add a scope. The Lancet2 chglog header pattern is `^(\w*)\:\s(.*)$` — type and subject only, no scope group. Writing `fix(caller): ...` does not parse. If the change is layer-specific, mention the layer in the subject text or in the body, e.g., `fix: preserve missing-value emission for caller's ASMD field on indels`.

Sixth, draft the subject line. The subject states what the change does in imperative mood, lowercase first letter, no trailing period, fitting within `max_subject_length` (seventy-two characters). Be specific about what changed, not what file was edited. "fix: preserve missing-value emission for ASMD on indels" is good; "fix: update sample_format_data.cpp" is not. If the subject would exceed the length limit, the right move is almost always to split the change into two commits, not to compress the subject into jargon.

Seventh, decide trivial versus substantial. `docs_dev/style/cpp_style.md` § Git commit messages exempts trivial commits from carrying a body. A change is trivial when it is one of:

- A typo fix in a comment, a help string, or documentation.
- An executable-bit change (`chmod +x`).
- A formatting-only change (clang-format, whitespace normalization, line-ending fixes).
- A pinned-dependency bump (`pixi.lock`, `pixi.toml` version bump, conda recipe version).
- A `.gitignore` addition or similar repository-hygiene change.

For trivial commits, write only the subject. Match one of the trivial patterns the validator recognizes: `chore: bump <thing> to <version>`, `chore: update <thing> to v<n>`, `chore: pin <thing>`, `chore: format ...`, `chore: clang-format ...`, `fix: typo ...`, `chore: typo ...`, `chore: exec bit ...`, `chore: gitignore ...`. The validator passes these silently regardless of commit type. If the subject matches one of these patterns but the staged diff is unexpectedly large (default threshold thirty lines insertions+deletions), the validator soft-warns rather than passes silently, so an accidental "trivial" subject on a five-hundred-line refactor gets caught.

For everything else, the change is substantial and the body shape depends on the diff size. The validator's threshold is thirty lines (insertions plus deletions); below that, a body is encouraged but not required (a one-paragraph context note is enough); above that, the validator hard-blocks a missing body. When in doubt, draft the body — the few seconds spent on a context paragraph save reviewers and future-you many minutes when the change is being read months later. Body-required is keyed on the diff size, not the commit type, so a thirty-line `feat:` and a thirty-line `chore:` are treated the same way.

Eighth, draft the substantial body in the two-section pattern. `docs_dev/style/cpp_style.md` § Git commit messages is precise about the shape. The body has two sections separated by a blank line.

The first section is one to three short context paragraphs explaining the *why*: the problem being solved, the user-visible effect, the rationale behind the chosen approach. Reference symbols, flags, and file paths in backticks (`` `--probe-variants` ``, `` `BuildHaplotypes` ``, `` `path.{h,cpp}` ``). Keep paragraphs short — the goal is review-ability, not exhaustive documentation. Lines wrap at the configured `max_body_line_length` (one hundred characters).

The second section is a file-list. Each file (or paired-file group) gets one bullet in the form `- path: change`. Use the brace shorthand `name.{h,cpp}` for header/source pairs that change together. When the diff spans multiple subsystems, group bullets under language sub-headers (`C++ changes:`, `Python changes:`, `Shell changes:`, etc.) with a blank line before each sub-header. The shape looks like:

```
fix: preserve missing-value emission for ASMD on indels

ASMD was emitting `0.0` for indel paths where the alt-allele soft-clip
subtraction underflowed past the reference length, producing a value
that downstream filters interpreted as a real (poor) alignment quality
score. The correct emission is `.` (untestable) per `vcf_header_builder.cpp`.

The fix moves the underflow check to `compute_asmd` and treats any negative
intermediate as untestable rather than clamping to zero.

- src/lancet/caller/variant_call.{h,cpp}: route untestable ASMD to `.`
- src/lancet/caller/sample_format_data.cpp: skip ASMD on indel paths
- tests/caller/asmd_test.cpp: regression case for indel underflow
```

For a multi-subsystem diff:

```
feat: add --probe-variants for known-truth diagnostics

The probe-diagnostics infrastructure in `cbdg/probe_index.{h,cpp}` could
not be exercised from the CLI; this exposes it via `--probe-variants <vcf>`,
loading the provided VCF, hashing each variant's k-mer set, and tracking
graph paths that recover those k-mers. The output goes to `<out>.probe.tsv`.

C++ changes:
- src/lancet/cli/cli_interface.{h,cpp}: register `--probe-variants` flag
- src/lancet/cli/cli_params.h: add `probe_variants_path` field
- src/lancet/core/probe_diagnostics.cpp: wire CLI input to probe runner

Python changes:
- python/score/probe_summary.py: parse the new `<out>.probe.tsv` format
```

Read the staged file inventory and produce one bullet per changed file (or per paired group). Match each bullet to the right language sub-header by file extension: `.cpp`, `.cc`, `.h`, `.hpp` go under "C++ changes:"; `.py` under "Python changes:"; `.sh`, `.bash` under "Shell changes:"; `cmake/` and `CMakeLists.txt` under "CMake changes:"; and so on. Skip language sub-headers entirely when the diff is single-language; they help only when the diff spans subsystems.

Ninth, present the proposal. Show the full message in a fenced block, ready to be committed. Above the block, state in one sentence which type was chosen and whether the body is trivial-omitted or substantial-included. Ask for confirmation before committing.

Tenth, on confirmation, commit. Use one `-m` for the subject and one `-m` per body section: `git commit -m "<subject>" -m "<context paragraphs as one block with blank lines between>" -m "<file-list bullets, with sub-headers if used>"`. Each `-m` becomes a paragraph separated by a blank line in the final commit message, which is what the chglog parser and the validator both expect. The validate_commit_message hook will run; if it blocks, read the error and fix the message. If the user declines or wants edits, do not commit; ask what to change and produce a revised proposal.

## Heuristics for hard cases

A diff that touches `src/lancet/caller/` and adds a new VCF FORMAT field is a `feat` whose context paragraph must mention the field's `Number=`, `Type=`, and `Description=` from the header builder. Run `vcf-validator` against the change before committing if you have not already.

A diff that touches `src/lancet/cbdg/` plus `tests/cbdg/` is `fix` or `feat`, not `chore`; the test is in service of the source change. The test file goes in the file-list section under the same C++ sub-header as the source it exercises.

A diff that bumps `pixi.lock` along with code changes is two commits, not one. Stage and commit the `pixi.lock` change separately as `chore: bump pixi.lock to <version>` first (trivial, no body needed), then the code change (substantial, two-section body).

When in doubt about whether a change is `fix` or `chore`, ask: would a user notice this change in their VCF output, in their CLI behavior, or in their wall-clock time? If yes, it is `fix` (or `feat` or `perf`); if no, it is `chore`.

When a body would need a fourth or fifth context paragraph, the change is too large for one commit. Split it.

## What the command does NOT do

It does not stage files. Staging is the developer's decision; conflating staging with committing produces too-large commits.

It does not run the test suite, the linter, or the validator before committing. The Stop hook and `/fix-and-validate` are for that.

It does not push. After committing, the user reviews the result and pushes manually.

It does not amend. If the user wants to amend the previous commit, that is `git commit --amend`; this command is for new commits.
