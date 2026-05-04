# Agents

Six subagents live here. Each one earns its place by satisfying two conditions that distinguish a subagent from a skill or an AGENTS.md note: the work is large enough that doing it in the main session would push out useful state, and only the conclusion needs to come back, not the search trail.

## What's here

- **`fresh-reviewer`** — pre-merge review of a diff against main. Read-only, opinionated, no memory of how the change was implemented. Returns a prioritized findings list with file:line references. Flags undisclosed NOLINT suppressions and bare-NOLINT additions.
- **`assembly-and-calling-expert`** — variant-discovery pipeline expert covering `cbdg/` (de Bruijn graph assembly) and `caller/` (MSA, scoring, genotyping, VCF construction). Read-only. Use for review or reasoning about correctness, scoring math, FORMAT/INFO/FILTER semantics, or graph algorithm invariants.
- **`vcf-validator`** — VCF v4.5 schema validation when adding, removing, or changing FORMAT/INFO/FILTER fields. Read-only.
- **`perf-analyst`** — interprets gperftools/pprof CPU profiles, identifies hotspots, proposes optimizations with trade-off analysis.
- **`sanitizer-expert`** — analyzes existing sanitizer reports (ASan/TSan/UBSan/MSan output) when you want focused interpretation without a full reproduce-fix workflow. The full workflow lives in the `sanitizer-build-analysis` skill (the agent does analysis, the skill does the procedure).
- **`probe-interpreter`** — reads a probe tracking analysis report and produces focused recommendations grounded in C++ source. Maps each `lost_at_stage` to the responsible code, proposes specific files/functions/parameters to investigate, frames in terms of sensitivity vs specificity. May write a written-up analysis next to the raw probe outputs in `notes/probe-debug-<date>/`. The operational mechanics of running the workflow live in the `probe-tracking` skill.

## Why subagents and not skills?

Skills don't have isolated context windows. The six agents here all benefit specifically from context isolation: `fresh-reviewer` needs to read the diff cold, `vcf-validator` walks the v4.5 spec without polluting the main session, `perf-analyst` digs into profile state that would otherwise dominate the conversation, `assembly-and-calling-expert` reasons about deep correctness without bringing all of that reasoning back into the main session, `sanitizer-expert` analyzes multi-thousand-line traces, `probe-interpreter` walks a 27-stage attribution cascade and traces each stage to its C++ owner.

For information-only "expert" content where context isolation doesn't pay for itself, prefer a skill — skills are easier to maintain, easier to compose, and cheaper at the margin (no per-session description cost). The choice between subagent and skill for "expert" content is whether the context-isolation benefit matters for the use case.

## Why few agents?

The bar is delegation, not taxonomy: an agent earns its slot when context isolation pays for the per-session description cost (see `.claude/cost-model.md`). One agent per source-tree layer or one per kind of question doesn't clear that bar — when Claude's natural action when asked about a layer is to read the actual files, the agent provides no delegation value over what reading provides.

The six kept here all pass the delegation test. Pre-merge review legitimately benefits from a fresh context (otherwise you're reviewing your own writing through your own writing). VCF schema validation is bulky read-only work. Performance analysis benefits from context isolation while trade-offs are reasoned through. Sanitizer-report analysis is a focused task with structured output. The variant-discovery expert covers the two most algorithmically-dense layers (`cbdg/` and `caller/`) as a delegation pattern (when reasoning gets deep, hand it off). Probe-interpretation maps a structured stage-attribution report onto specific C++ source locations and is bulky enough that doing it inline would crowd out the calling session's other state.

## Description-as-trigger discipline

Each agent's description field is what Claude uses to decide when to delegate to that agent. The description should read as a routing rule, not as a label. "Use proactively before a commit, after writing a feature, or when the user says..." routes well; "Pre-merge code reviewer" reads as a label and doesn't fire reliably. Anthropic's official subagents documentation is explicit on this: "the description field is what Claude uses to decide when to delegate. Be specific about the trigger conditions, not just the capability."[^subagents-docs]

When you find an agent isn't firing when it should, the first fix is the description, not the body. Tighten the description to include the user phrases you actually say when you want that agent's help.

[^subagents-docs]: <https://code.claude.com/docs/en/sub-agents>

## Maintenance lifecycle

### Adding a new agent

Adding an agent should pass two filters before you write a single line of YAML.

**Filter 1 — does it need context isolation?** If the work could be done in the main session without polluting context, it should be a skill, not an agent. Ask: does this work involve reading hundreds of lines of state (a profile, a sanitizer trace, a long spec) that would otherwise sit in the main session for the rest of the work? If no, it's a skill.

**Filter 2 — would I actually delegate to it?** Imagine the workflow that would invoke this agent. If your imagined invocation is "Claude, please consult the X expert about Y" and you cannot articulate why you wouldn't just have Claude read X's files directly, the agent provides no value.

If both filters pass, write the agent. The frontmatter must include `name`, `description` (as a trigger, not a label), `tools` (allowlist; restrict to read-only when the agent is for analysis or review), and `model` (typically `opus` for reasoning-heavy work; `sonnet` for routine tasks; `haiku` for simple sweeps). Update this README's "What's here" section.

Cost reasoning for the addition: see `../cost-model.md`. The added per-session baseline cost is roughly 100-200 tokens for the description; the agent body is paid only on invocation.

### Deleting an agent

Delete an agent when it has not been invoked in the last quarter. Concretely: if you sit with the bundle for the quarterly review and you cannot remember the last time the agent fired (or the last time you typed `@agent-name`), it has not been pulling its weight.

The deletion procedure: remove the markdown file, update this README's "What's here" section, add a brief note to the per-directory micro-changelog at the bottom of this file. Git history preserves the file content if you ever need it back.

Common deletion patterns: agents that overlap with another agent's trigger (the lower-priority one goes); agents whose description has decayed to a label (delete; the agent wasn't firing anyway); per-layer experts that turn out to be reference content rather than delegation work (downgrade to skill, then delete the agent).

### Refactoring an agent

Three refactor patterns are common.

**Tightening the description.** When an agent isn't firing reliably, rewrite the description to include the actual phrases the user says. Read several days of session transcripts and note the words used when this agent's help was wanted; weave those into the trigger.

**Pruning the body.** Agent bodies accumulate "watch for" items over time. When a body crosses ~250 lines, walk every item and ask: have I seen this issue in the last quarter? If no, delete. If you need it later, git history will surface it.

**Splitting or merging.** If an agent's description has grown to cover three or four genuinely different triggers, splitting into two agents (each with a tighter trigger) often lets each fire more reliably. Conversely, when two agents have overlapping triggers and Claude routes inconsistently between them, merging into one with the union of expertise can be a win.

### Reviewing the agent set

Quarterly. Walk all six agents in turn:

- Has this agent been invoked? When? Was the result useful?
- Does the description still reflect what the agent does? Have I been saying things that should fire it but don't?
- Does the body still match the codebase? Have any of the cited file:line invariants moved?
- Are the tool restrictions still right? Should it be more restrictive (e.g., a review agent that's been editing files)?

The grounding for what the agent claims about source should match the actual source. If `assembly-and-calling-expert` cites `genotype_likelihood.cpp:45` for CMLOD computation and that line has moved, the agent will give wrong file:line citations until the body is updated. Periodic re-grounding is the price of having agents that cite specifically.

### Retiring an agent

Retire an agent when it has been on the deletion candidate list for two quarters in a row. The two-quarters rule prevents the "I might need this" purgatory.

If you've decided to retire and then realize there's a workflow that legitimately uses the agent, articulate the workflow in the agent's body before retiring is reversed. The articulation forces you to confront whether the use is real.

## Cost model

See `../cost-model.md` for the full breakdown of how subagent descriptions and bodies cost context tokens. Briefly: each agent's description pays per-session baseline cost (loaded into Claude's awareness at session start so the main agent knows what subagents exist); the body costs nothing until invocation, at which point it loads into the subagent's isolated context window (not the main session's). This makes subagents the right choice when you want bulky reasoning kept out of the main conversation.
