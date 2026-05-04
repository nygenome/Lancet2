# Bundle cost model

This document explains how each Claude Code mechanism (CLAUDE.md, hooks, subagents, skills, slash commands) costs context tokens, and the reasoning behind the bundle's structure. Every directory README links here. The model evolves as Claude Code changes; the `/sync-cost-model` slash command refreshes this file against the current state of Anthropic's documentation.

This document is the single point that the directory READMEs cite for cost reasoning. Keep cost discussions out of the directory READMEs — they should reference this and stay narrow.

## What costs what

The five mechanisms have meaningfully different cost shapes. Treating them the same is the most common bundle-design error.

**CLAUDE.md / AGENTS.md** is paid per message in every session. The full text is loaded into Claude's context as a system-message prefix at every interaction. In Lancet2's setup, CLAUDE.md is a thin wrapper that imports AGENTS.md (`@AGENTS.md`); AGENTS.md carries the canonical content. Anthropic's costs guide ([code.claude.com/docs/en/costs](https://code.claude.com/docs/en/costs)) explicitly recommends keeping CLAUDE.md under 200 lines, citing the per-message tax as the reason — that line target applies to AGENTS.md after the migration since the import resolves to its full content. This is the most expensive mechanism per character of content; reserve it for rules that genuinely should shape every conversation.

**Subagent description fields** are paid per session. When a session starts, Claude Code scans `.claude/agents/` and loads each subagent's frontmatter into Claude's awareness so the main agent knows what subagents exist and what they do. A subagent body is NOT loaded baseline — it only enters context when the subagent is actually invoked (and even then, in its own isolated context window). The cost per subagent is on the order of 50-200 tokens of always-loaded description, paid every session. The current six subagents × ~100 tokens each = ~600 tokens of session-baseline cost. The choice of how many subagents to keep is partly a cost decision.

**Skill description fields** are similarly paid per session — Claude needs to know what skills exist to decide when to fire them. A skill body is NOT loaded baseline; it loads only when invoked, into the main session's context (not an isolated one like a subagent's). The trade-off: skills are cheaper baseline than subagents (no isolated-context overhead), but invoking a skill pollutes the main session's context with the skill body (whereas a subagent invocation does not).

**Slash commands** are essentially free baseline. Claude Code lists slash command names so the user can type them, but the command bodies are not loaded into Claude's context until the user explicitly invokes one. A slash command costs nothing until typed; once typed, the body is added to the main session's context.

**Hooks** cost zero context tokens. They are shell scripts that run outside the model's context entirely. They communicate with the model via stdout/stderr and exit codes. The tradeoff is that hooks cannot do anything that requires LLM judgment — they are deterministic.

## How the bundle's choices map to this

The bundle currently has 6 subagents, 12 skills, 16 slash commands, 8 hooks, AGENTS.md targeting ~150-180 lines after pruning, and a CLAUDE.md that is a 3-line wrapper importing AGENTS.md. The reasoning, mechanism by mechanism:

**Why few subagents.** Subagent descriptions are loaded every session, and the marginal value of adding another subagent has to exceed the per-session cost. The bar is whether the work would pollute the main session's context if done inline; if yes, it earns a subagent slot. The six current subagents all clear that bar: pre-merge review (fresh context), VCF schema validation (bulky read-only spec walk), performance analysis (profile state would dominate the main session), sanitizer report analysis (multi-thousand-line traces), variant-discovery code review (deep correctness reasoning over `cbdg/`+`caller/`), and probe-tracking interpretation (a 27-stage attribution cascade traced to specific C++ source). Taxonomy-driven agents — one per source-tree layer or one per kind of question — provide little delegation value when Claude's natural action is to read the actual files anyway.

**Why some skills are reference skills.** `assembly-and-calling-expert` is a subagent because the work it does (reasoning about caller/cbdg correctness) benefits from context isolation when the analysis is deep. Reference content that's just information could equivalently be a skill — the choice between subagent and skill for "expert" content is whether the context-isolation benefit matters for the use case.

**Why slash commands rather than skills for /fix-and-validate, /e2e-pipeline-test, /commit, /brainstorm, /spec, /execute-spec, /wrap-branch.** These are user-initiated workflows where the user wants explicit control over when they fire. Skills can be invoked without the user typing a command, which is good for "consult on demand" patterns but bad for "run this specific procedure right now." For `/sync-cost-model` specifically, the user-initiated discipline matters even more — the maintenance ritual should be deliberate, not implicit.

**Why hooks for layer/naming/commit-message rules.** These are rules that must be enforced even when Claude has been told to ignore them. Hooks are deterministic and cannot be persuaded. AGENTS.md cannot enforce; it can only request.

**Why AGENTS.md is small.** Per-message cost. Rules that shape every conversation belong here; rules that apply only in specific workflows belong in skills or path-scoped rules under `.claude/rules/`.

## Anthropic-published guidance and citations

Sources current as of 2026-04. Re-verify when running `/sync-cost-model`.

- **CLAUDE.md ≤ 200 lines:** [code.claude.com/docs/en/costs](https://code.claude.com/docs/en/costs) — "Aim to keep CLAUDE.md under 200 lines by including only essentials."
- **Subagent description as trigger, not summary:** [code.claude.com/docs/en/sub-agents](https://code.claude.com/docs/en/sub-agents) — "the description field is what Claude uses to decide when to delegate. Be specific about the trigger conditions."
- **Skills load on demand:** [code.claude.com/docs/en/skills](https://code.claude.com/docs/en/skills) — "Auto-compaction carries invoked skills forward within a token budget. When the conversation is summarized to free context, Claude Code re-attaches the most recent invocation of each skill."
- **Hook lifecycle (25 events, exit code 2 blocks):** [code.claude.com/docs/en/hooks](https://code.claude.com/docs/en/hooks) and various community writeups.
- **Cost reduction case study:** Cem Karaca, "My CLAUDE.md Was Eating 42,000 Tokens Per Conversation" (Medium, 2026-02) — demonstrated 94% always-loaded token reduction by moving CLAUDE.md content to skills.

## Concrete numbers (rough, for planning only)

The numbers below are rough approximations. They are useful for planning relative costs (is this 100 tokens or 10,000 tokens?), not for budgeting precise quotas. Re-measure if you need accuracy.

| Mechanism | Where the cost lands | Rough size |
|:----------|:---------------------|:-----------|
| AGENTS.md (150 lines) via @AGENTS.md import | Every message, every session | ~1,500-2,500 tokens always loaded |
| Subagent description (each) | Session start, every session | ~100-200 tokens always loaded |
| Skill description (each) | Session start, every session | ~100-200 tokens always loaded |
| Subagent invocation | When invoked, isolated context | Subagent body + reasoning, kept out of main |
| Skill invocation | When invoked, main context | Skill body added to main session |
| Slash command name | Always available, not loaded | ~0 tokens until typed |
| Slash command body | When user types | Added to main session |
| Hook execution | Outside model context | ~0 tokens |

For the current bundle (6 subagents + 12 skills + ~150-180 line AGENTS.md + 6 path-scoped rules whose descriptions are always loaded): always-loaded baseline is roughly 4,000-6,500 tokens per session start, depending on description sizes and AGENTS.md final length. That's reasonable — well under the 200K context window, and small enough that a session has plenty of room for actual work. Path-scoped rule bodies (122-180 lines each) load only on glob match, so a session that touches one or two layers loads an additional 1,500-3,000 tokens of rule content. Total typical-session baseline (with rules loaded for an active layer): 5,500-9,500 tokens, comfortably below the 200K threshold where Opus pricing changes.

## Maintenance: what triggers an update to this document

Update this document when any of the following happens. The `/sync-cost-model` slash command formalizes the refresh procedure.

**Anthropic publishes new guidance.** New Claude Code features, changes to costs guidance, new lifecycle events. The Anthropic-published guidance section above lists the docs to re-check.

**Claude Code releases major changes.** New mechanism types (e.g., when plugins were introduced, the cost model needed an update). Changes to how skills load (e.g., the auto-compaction rules around skills). Changes to subagent inheritance.

**The bundle's structure changes meaningfully.** When any mechanism count changes substantively (e.g., the agent set grows or shrinks by half, or the skill set crosses ±5), update the "How the bundle's choices map to this" section to reflect the current numbers and reasoning.

**Community converges on a new pattern that affects cost.** New community-curated guidance on bundle structure that contradicts what this document says — e.g., if the practitioner consensus shifts on whether subagents or skills are the right home for "expert" content.

The cost numbers are the most volatile piece. Approximate token sizes change as Claude Code optimizes its bundling and as model context windows evolve. If the numbers in the table feel stale, run `/sync-cost-model` to verify against current docs and re-measure your own bundle.
