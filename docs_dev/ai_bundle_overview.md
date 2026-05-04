# About this bundle (`.claude/` and `docs_dev/`)

This document explains *what* the Claude Code bundle in this repo is, *why* it exists, and *how the pieces fit together*. It's for new contributors (human or AI) who walk into the repo and wonder why there's a `.claude/` directory and a `docs_dev/` directory and what to do with them.

For the operational details — how each subdirectory works, how to add or retire something, what each mechanism costs in context tokens — read the README inside the relevant subdirectory and `.claude/cost-model.md`.

## What this bundle is

Configuration that runs alongside Lancet2 development. It exists to make Claude consistent and correct on this specific codebase — not as a general-purpose Claude Code starter. Two trees:

- `.claude/` — agents, skills, slash commands, hooks, settings, scripts, and the cost model. This is what Claude Code reads.
- `docs_dev/` — developer documentation: style conventions, subsystem deep-dives, workflow guides, architectural decision records, and investigation writeups. Excluded from Claude's default load via `.claudeignore`; opt-in for both humans and Claude.

Top-level files (`AGENTS.md`, `CLAUDE.md`) sit at the repo root. `AGENTS.md` is the canonical project-memory document, cross-tool portable. `CLAUDE.md` is a thin wrapper that imports it via `@AGENTS.md` so Claude Code finds the content through its conventional path.

## Design philosophy

The bundle implements the standard Claude Code mechanism layering. Each mechanism is used for what it's best at; using one mechanism to do another's job is the most common bundle anti-pattern.

**Hooks** for rules that must be enforced deterministically. Layer direction, naming, dangerous bash, protected paths, commit format, pre-commit summary. Hooks fire every time, run outside the model's context window, and cannot be persuaded out of by any context signal. Every hook is technical debt, so the set is kept minimal — eight hooks, collectively encoding the project's hardest correctness rules.

**Subagents** for delegation patterns where context isolation pays for itself. Pre-merge review, VCF schema validation, performance analysis, sanitizer report analysis, variant-discovery code review, probe-tracking interpretation. Six of them. The discriminator is whether the work would pollute the main session's context if done inline; if yes, it earns a subagent slot.

**Skills** for procedures that benefit from being written down once and applied consistently. Twelve of them. Skills load on demand and cost nothing baseline; their cost is paid only when invoked, into the main session's context.

**Slash commands** for one-keystroke invocations of repeatable workflows. The sixteen commands group by purpose: validation (`/fix-and-validate`, `/e2e-pipeline-test`), commit hygiene (`/commit`), feature kickoff (`/brainstorm`, `/spec`, `/execute-spec`), end-of-branch (`/wrap-branch`), document scaffolding (`/arch-decision-record`, `/investigate`), probe-tracking forensics (`/probe-concordance`, `/probe-run`, `/probe-analyze`), and bundle hygiene (`/sync-cost-model`, `/audit-bundle`, `/audit-vcf-schema`, `/audit-probe-pipeline`). These are user-initiated and explicit by design — workflows the user wants to start deliberately rather than have inferred from context.

**`AGENTS.md`** for rules that should shape every conversation. Layer architecture, naming conventions, commit-message rules, the chr1/chr4 (not chr22) test data correction. Kept under 200 lines per Anthropic's costs guide; per-message cost is the binding constraint.

## What this bundle does NOT do

It does not duplicate the project's style conventions. Those live in `docs_dev/style/` (entry point: `docs_dev/style/README.md`), and the bundle references them rather than restating them.

It does not encode general C++ best practice. Claude knows that already; restating it just adds noise that dilutes the project-specific signal.

It does not have one agent per source-tree layer. Taxonomy is not delegation — when Claude's natural action when asked about a layer is to read the actual files, an agent that says "consult me about this layer" provides little value over what reading provides. The bundle has one expert that spans the two algorithmically-densest layers (`assembly-and-calling-expert` covering `cbdg/` + `caller/`) plus task-shaped specialists for review, schema, perf, sanitizer reports, and probe-tracking interpretation.

It does not push every workflow into a slash command. Most workflows should be inferred from context (skill descriptions trigger them) rather than invoked by hotkey. Slash commands are reserved for workflows where explicit user intent matters.

It does not include MCP servers. MCP servers add per-process maintenance burden and per-call token cost, and the bundle has not yet identified a workflow where that overhead pays for itself. Re-evaluate if the workflow shifts.

## How contributors use the bundle

Open Claude Code in the repo root and start working. The hooks fire automatically on relevant tool invocations. The skills and agents activate when Claude judges them relevant (their descriptions are written as triggers, not as labels — "Use when..." rather than "This is..."). `AGENTS.md` applies always. Slash commands run when you type them.

For a substantive feature, use `/spec` first — it interviews you, drafts a spec in `notes/<FEATURE>/`, and produces a kickoff prompt for a fresh execution session. For day-to-day work, just talk to Claude; the bundle's defaults will route correctly most of the time.

Each subdirectory under `.claude/` and `docs_dev/` has its own README that explains what's there and how to keep it healthy. Read those when you need to add, modify, or retire something in that directory.

## Cost model

The reasoning behind why there are six subagents and not thirteen, why `AGENTS.md` is small, why `/spec` is a slash command rather than an agent — all of this comes from the cost-model document, which explains how each mechanism is paid for in context tokens. See `.claude/cost-model.md`. Every directory README links there for the cost reasoning rather than restating it.

## Maintenance principles

The bundle is developer infrastructure. Developer infrastructure goes wrong in two specific ways: it grows unbounded as you add things you might need someday, or it stagnates as you stop tending it because it works "well enough." The principles below are aimed at both failure modes.

### When to add to the bundle

Add when you find yourself frustrated by Claude doing something wrong on this codebase, AND the wrongness is repeatable enough that fixing it once will save you from fixing it again. Single-occurrence frustrations are not bundle-worthy; pattern-of-frustration is.

The right tier for the addition follows from the four-mechanism framework above. Ask in this order: can a hook enforce this? (deterministic, blocks). Can a skill encode this? (procedure with judgment). Does this need a subagent? (delegated context-isolated work). Should this just be an `AGENTS.md` note? (rule that shapes everything). When you can answer "yes" to two tiers, prefer the lower-cost one — skills over `AGENTS.md`, hooks over skills, slash commands over hooks for user-initiated rituals.

### When to delete from the bundle

Delete when a thing in the bundle has not pulled its weight in the last quarter. Specific signals: a subagent you have not invoked, a skill that has not fired, a slash command you have not typed, a hook that has not blocked anything (or has blocked things you wished it hadn't).

The bias should be toward deletion. Adding feels productive; deletion feels like loss. But every kept-but-unused element pays context-token cost on every session, and the cumulative noise dilutes the signal of what's actually project-specific. When in doubt, delete; if you find you needed the deleted thing, the git history makes it cheap to bring back.

### When to refactor

Refactor when a thing in the bundle is being used but the structure has decayed. Specific signals: an `AGENTS.md` section that's grown to half a page, a skill body that's accumulated edge cases until it reads like a checklist, a subagent description that no longer matches what the agent actually does, a hook with cumulative caveats that have made it hard to reason about.

Refactoring usually means demoting (`AGENTS.md` content → path-scoped rule, `AGENTS.md` content → skill, skill → reference doc), splitting (one mechanism doing two jobs → two mechanisms each doing one), or merging (two mechanisms with overlapping triggers → one with a clear trigger). The most common refactor is "this `AGENTS.md` section has gotten big enough that it should be a skill or a path-scoped rule" — Cem Karaca's writeup ([medium.com/@cem.karaca](https://medium.com/@cem.karaca/my-claude-md-was-eating-42-000-tokens-per-conversation-heres-how-i-fixed-it-85ffba809bd4)) walks through the `CLAUDE.md` → skills variant of this transformation in detail and is worth re-reading when `AGENTS.md` feels heavy.

### When to review the bundle as a whole

Quarterly. Sit with the bundle for an hour, walk every file, ask of each: is this still pulling its weight? Has the codebase changed in ways that make this stale? Has Claude Code shipped a new capability that supersedes this? Are descriptions still working as triggers, or have they decayed back into labels?

The quarterly review is also when to run `/sync-cost-model` to re-verify the cost reasoning against current Anthropic guidance.

### When to retire something entirely

Retire when a thing has been on the deletion candidate list for two quarters in a row. Some things sit in "I might need this" purgatory for years; the two-quarters rule forces a decision.

If a thing is providing genuine value and you can articulate when, write that articulation into the relevant README and keep it. If you can't articulate it, retire.
