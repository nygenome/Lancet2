# Skills

Twelve skills live here. Each one is a procedure that Claude applies with judgment — not a checklist, not a rulebook. Project-wide rules go in `AGENTS.md`; deterministic enforcement goes in hooks; everything that benefits from being written down once and applied consistently across multiple invocations goes here.

## What's here

- **`add-cpp-test`** — adding a unit or integration test, scaffolding a Catch2 fixture, doing TDD work. Enforces the project's test layout and TDD discipline. Has a `references/` directory with the Catch2 idioms cheat sheet and the bench-binary invocation reference.
- **`profile-and-optimize`** — full performance-improvement workflow: build the profiling tree, capture against real data, analyze via pprof (optionally with the `perf-analyst` subagent), validate with the google/benchmark suite, validate VCF identity. Composes with `perf-analyst`. Has a `references/` directory with four reference files (google/benchmark idioms, bench invocation, gperftools internals, pprof and analyze_profile.py).
- **`sanitizer-build-analysis`** — full procedure from "I have a repro" to "merged with regression test." Drives pixi-managed sanitizer build trees (cmake-build-{asan,msan,tsan,ubsan}). Has a `references/sanitizer_matrix.md` with per-sanitizer detail. Different from the `sanitizer-expert` agent, which does analysis-only of an existing report.
- **`systematic-debugging`** — methodology for any bug, test failure, or unexpected behavior. The Iron Law: no fixes without root-cause investigation first. Has two reference files: `root-cause-tracing.md` walks the trace-to-cause procedure; `defense-in-depth.md` covers tiered prevention so the same bug class doesn't recur.
- **`doc-sync`** — keeps in-source comments, developer docs in `docs_dev/`, and user docs in `docs/` in sync with each other and with code. Two modes: Mode A is targeted sync after a code change (find every cross-layer reference, propose a coordinated edit set, gate on a consistency check); Mode B is a periodic semantic audit (file-by-file walk of a subdirectory applying a five-pass methodology — structural compliance, comment quality, language compliance, header-vs-implementation split, cross-file synchronization — producing a tiered backlog of findings). Excludes `docs_dev/investigations/` (immutable by design), agent-memory files, `.claude/` READMEs, `AGENTS.md`, and `CHANGELOG.md`.
- **`test-data-locations`** — documents the GCS dataset (gs://lancet2-test-datasets/test_harness_data/), the env vars that map to each fixture, and a decision tree for which fixture fits which workflow.
- **`probe-tracking`** — operational playbook for the probe variant forensic pipeline (truth_concordance → Lancet2 with --probe-variants → analyze_probe_results). Covers inputs, the flag dance between steps, the data layout, --verbose requirements, and the somatic-no-truth-VCF gap. For interpreting an existing report, delegates to the `probe-interpreter` subagent.
- **`external-interface-changes`** — interview-driven workflow for VCF FORMAT/INFO/FILTER changes and CLI flag changes. Walks the five-operation matrix (add / rename / cardinality-change / remove / silent-semantic-change), proposes coordinated edits across vcf_header_builder.cpp + caller code + tests + docs, and invokes the `vcf-validator` subagent before merge.
- **`python-script-sync`** — maintains the contract between the C++ pipeline and the Python tooling under `scripts/` (analyze_*, run_clang_*, build_*, bump_version, update_changelog, etc.). Catches drift when the C++ side changes a TSV column or log format that a script consumes.
- **`cmake-sync`** — covers the full chain of changes for adding source files, layers, dependencies, build options, link-line changes, and compile-flag changes. Aware of the protected-paths boundary around `cmake/` and `pixi.toml`.
- **`clang-tidy-discipline`** — procedural how-to for resolving clang-tidy violations and shaping NOLINT suppressions. Covers diagnose-before-acting (most violations are real and should be fixed, not suppressed), the canonical scoped suppression forms (`NOLINTNEXTLINE` for one line; `NOLINTBEGIN`/`NOLINTEND` for a block), the cognitive-complexity ceiling that gates further function growth, and the project-specific reason `clang-tidy --fix` is forbidden. Defers to `docs_dev/style/cpp_style.md` for the canonical rule statements.
- **`release-notes`** — distills a CHANGELOG range into a categorized impact summary (breaking / schema / feature / fix / perf / internal). Complements but does not replace the chglog-generated CHANGELOG.md; produces the announcement-style narrative for users.

## Why skills and not AGENTS.md content?

AGENTS.md is paid in input tokens on every message of every session (Claude Code reads it through the CLAUDE.md wrapper). A multi-step procedure that's used once a week shouldn't pay that ongoing cost. Skills load on demand — Cem Karaca's writeup ([medium.com/@cem.karaca](https://medium.com/@cem.karaca/my-claude-md-was-eating-42-000-tokens-per-conversation-heres-how-i-fixed-it-85ffba809bd4)) shows the cost differential is real: moving a 1,200-line CLAUDE.md to skills cut always-loaded tokens by 94%.[^cem] The same principle applies to AGENTS.md.

A useful test: would Claude consistently get this wrong without it being a standing instruction? If yes, AGENTS.md (it must shape every conversation). If no — if Claude only needs the procedure when the user is doing the thing — skill.

[^cem]: <https://medium.com/@cem.karaca/my-claude-md-was-eating-42-000-tokens-per-conversation-heres-how-i-fixed-it-85ffba809bd4>

## Why skills and not slash commands?

Slash commands are user-typed shortcuts. Skills are procedures Claude invokes based on context — when you say "this is slow, can you profile it" Claude should reach for `profile-and-optimize` without you typing `/profile`. Slash-command form would force the explicit invocation and bury the procedure behind keystrokes.

The exception is workflows that benefit from explicit user intent — `/fix-and-validate`, `/e2e-pipeline-test`, `/commit`, `/brainstorm`, `/spec`, `/execute-spec`, `/wrap-branch`, `/sync-cost-model`. Each of these earns slash-command form because the user wants deliberate invocation rather than context-driven inference.

## Description-as-trigger discipline

The official Claude Code skills documentation says descriptions should be triggers that tell Claude when to fire, not labels that describe what's inside.[^skills-docs] "Use when the user asks to..." or "Trigger on..." or "When investigating..." are the right shapes. If the description reads as a noun phrase ("Performance optimization guide"), the skill won't fire when it should.

When you find a skill isn't firing when it should, the first fix is the description. Read several days of session transcripts and note the words used when this skill's help was wanted; weave those into the trigger.

[^skills-docs]: <https://code.claude.com/docs/en/skills>

## Skill structure conventions

- Skills are folders, not files. A skill named `foo` lives at `.claude/skills/foo/SKILL.md` (and may include supporting files in `references/`, `scripts/`, `examples/` subdirectories).
- The YAML frontmatter at the top of `SKILL.md` is required. It must include `name` (matching the directory) and `description` (as a trigger). It may include `allowed-tools` (allowlist).
- The body should give goals and constraints rather than railroad Claude into prescriptive steps. The community-curated guidance is consistent on this: prescriptive checklists make skills brittle.[^community]
- Include a "When NOT to use this skill" section at the end. The negative case sharpens the trigger and prevents the skill from firing on adjacent-but-different problems.

[^community]: shanraisshan/claude-code-best-practice (GitHub) and similar curated repos converge on this principle.

## Maintenance lifecycle

### Adding a new skill

A skill earns its place when:

1. The procedure is repeatable enough that re-deriving it each time wastes effort.
2. The procedure has judgment calls (otherwise it's a hook or a script, not a skill).
3. The audience is Claude doing work, not a human reading docs (otherwise it's a doc in `docs/`).

Workflow for adding: write the SKILL.md with frontmatter (name, description-as-trigger, allowed-tools). Body: a short rationale up front, then the procedure as goal-driven steps with constraints, then a "When NOT to use" section. Add the new skill to the "What's here" section of this README. If the skill replaces or overlaps with existing content (in AGENTS.md, in another skill, in a slash command), explicitly note the overlap and decide whether the older content should be deleted or kept narrower.

Cost reasoning: each skill description pays per-session baseline cost (Claude needs to know what skills exist to decide when to invoke); the body is paid only when invoked. See `../cost-model.md`.

### Deleting a skill

Delete a skill when:

- It hasn't fired in the last quarter and you have no specific upcoming workflow that needs it.
- Its content has migrated to a more appropriate home (AGENTS.md, a hook, a doc in `docs/`).
- Its trigger overlaps with another skill and the other one fires more reliably.

The bias should be toward deletion. Every skill description pays per-session baseline cost. Skills you don't use accumulate noise in Claude's "what skills are available" awareness, which dilutes the trigger signal for skills that do matter.

Procedure: delete the directory, update this README's "What's here" section, add a brief micro-changelog entry at the bottom. Git history preserves the file content.

### Refactoring a skill

Common patterns:

**Trigger refinement.** When a skill fires at the wrong times or doesn't fire at the right ones, rewrite the description. The body usually doesn't need to change — the issue is routing.

**Body pruning.** Skill bodies accumulate edge cases. When a skill body crosses ~150-200 lines, walk every section and ask: is this used? If a section hasn't been hit in two quarters, delete.

**Merging.** Two skills with overlapping triggers may merge into one. The signal is inconsistent routing — Claude fires skill A on a query that skill B should have handled, and vice versa, because the descriptions overlap on the same trigger phrases. Merging into a single skill with the union of expertise resolves the routing ambiguity at the cost of a longer body.

**Promotion or demotion.** A skill that's grown to cover something that should be enforced rather than suggested may need to become a hook. A skill whose content is reference rather than procedure may need to become a doc in `docs/`.

### Reviewing the skill set

Quarterly. Walk all twelve skills in turn:

- Has this skill fired? When? Was the result useful?
- Does the description still reflect what the skill does? Has the codebase changed in ways that make the procedure stale?
- Are the cited file paths and tool invocations still correct? (Skills that mention `pixi run X` need updating when X is renamed.)
- Could this skill be merged with another, or split into two?

Re-grounding is the most likely maintenance need. Skills that say specific things about the codebase (paths, function names, test fixtures) drift as the code drifts. Periodic re-grounding catches this before a skill produces wrong instructions.

### Retiring a skill

Retire a skill when it has been on the deletion candidate list for two quarters in a row. Articulating "this should retire" twice without acting on it is a signal that the skill is providing some marginal value worth examining; either commit to keeping it (and write down why) or commit to retiring it.

## Cost model

See `../cost-model.md` for the full breakdown of how skill descriptions and bodies cost context tokens. Briefly: each skill's description pays per-session baseline cost; the body costs nothing until invocation, at which point it loads into the main session's context (not an isolated one like a subagent body). This makes skills cheaper than subagents at the margin, and the right choice when context isolation isn't needed.
