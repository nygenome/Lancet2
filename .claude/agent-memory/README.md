# Agent memory

This directory carries per-agent memory files that persist across
sessions. Each subagent has its own file (`<agent-name>.md`) where it
records the kinds of knowledge that compound over time — bug patterns
encountered, struct-layout decisions made, past-PR resolutions worth
remembering, architectural understanding evolution, and a REJECTED
decisions log so old anti-patterns don't get re-proposed.

## Scope: project vs user

Five of the six memory files are **project-scoped**. They contain
knowledge that's specific to Lancet2 and useful to any developer who
opens the project — bug patterns in cbdg/ graph code, architectural
decisions about VCF schema evolution, etc. These files are
git-tracked and shared across the team.

One file is **user-scoped**: `perf-analyst.md`. Performance
observations are tied to specific hardware, kernel versions, BIOS
settings, and even ambient room temperature on the day a benchmark
ran. A finding like "this loop is 14% faster with prefetch on chr4"
generalizes poorly across the team's heterogeneous laptops/desktops/
clusters; storing such findings as if they were universal pollutes
the shared memory. The user-scoped file lives in the same directory
for discoverability but is excluded from git via a
`.claude/agent-memory/perf-analyst.md` line in `.gitignore`.

## Memory file format

Each agent's memory file is a single markdown document with three
top-level sections, in this order:

1. **Active knowledge** — facts and patterns the agent should pull
   into context whenever invoked. Kept short. New entries displace
   stale ones; the agent prunes during `/audit-bundle` review or on
   demand.

2. **Decision log** — chronological record of significant decisions
   the agent participated in (e.g., "decided to use Welford's algorithm
   for OnlineStats merge over a naive Σx²/n on 2026-03-15 because of
   Q-score precision loss"). Append-only during normal operation;
   compacted during quarterly review.

3. **REJECTED decisions** — patterns/approaches that were considered
   and explicitly rejected, with the reason. This section prevents
   the agent from re-proposing the same anti-pattern in a future
   session ("we already tried X, here's why it didn't work"). Without
   this, the same conversation re-runs every quarter.

## Knowledge kinds per agent

Each agent's memory specializes by what knowledge accumulates fastest:

- `assembly-and-calling-expert.md` — bug patterns in graph
  construction, struct-layout decisions in cbdg/, raw_variant.cpp
  and variant_bubble.cpp identity-redesign history, SPOA
  configuration tradeoffs encountered.
- `fresh-reviewer.md` — review patterns that recur across PRs (e.g.,
  "third time this quarter someone forgot to update vcf-validator
  delegation"), false positives the reviewer learned not to flag.
- `vcf-validator.md` — schema invariants that have been re-derived
  in past reviews, edge cases (multi-ALT, missing-value), downstream
  consumers that broke when fields changed.
- `sanitizer-expert.md` — known-benign sanitizer warnings (with
  rationale), patterns of new warnings that turned out to be real
  bugs.
- `probe-interpreter.md` — common probe-tracking interpretive
  patterns ("when probes are missing in graph but present in raw
  reads, suspect anchor failure first").
- `perf-analyst.md` (user-scoped) — hardware-specific observations,
  per-laptop benchmark baselines, profile-and-optimize cycles that
  worked or didn't.

## Git tracking and chglog exclusion

Project-scoped memory files are git-tracked and bundled into the
next feat/fix/perf/chore commit (no separate "update agent memory"
commits — they're invisible auxiliary state, not user-visible work).
The `pre_commit_summary.sh` hook surfaces the count of bundled memory
files in the commit summary so the user knows when memory is
implicitly traveling with their work.

Project-scoped memory updates always bundle into typed commits
(`feat:`, `fix:`, `perf:`, `chore:`) rather than landing in standalone
commits. Git-chglog filters on commit type, so memory updates appear
under their host commit's primary type rather than as their own
CHANGELOG entry. No additional `.chglog/config.yml` filter is needed
— the design relies on Claude never producing memory-only commits.
The `release-notes` skill operates on commit ranges (not individual
file paths within a commit), so memory file changes contribute to
the release narrative only through the host commit's user-impact
classification, not as standalone entries.

## Audit and accuracy review

`/audit-bundle` Pairing 12 reviews each project-scoped memory file
during PR review for accuracy. The audit specifically looks for:

- **Stale active knowledge** that contradicts current source code
  (e.g., a memory entry referring to a struct field that's been
  renamed).
- **Drift between memory and source** — if memory says "the SPOA
  Match score is +1" but msa_builder.h says `MSA_MATCH_SCORE = 0`,
  that's a finding.
- **Bloat** — memory files over 200 lines that have not been
  compacted.
- **Decision log entries** older than 1 year that haven't been
  promoted to AGENTS.md, retired, or moved to a layer rule.

The user-scoped `perf-analyst.md` is exempt from the team-PR audit
(it's per-developer), but the agent itself is encouraged to audit
its own file periodically.

## How to add a new agent

If a new subagent is added to `.claude/agents/`, create a matching
`<agent-name>.md` here with the three-section template:

```markdown
# <agent-name> agent memory

## Active knowledge

(empty)

## Decision log

(empty)

## REJECTED decisions

(empty)
```

If the new agent's knowledge is hardware/setup-specific (like
`perf-analyst`), add the file path to `.gitignore` to make
it user-scoped (gitignored, not shared with the team).
