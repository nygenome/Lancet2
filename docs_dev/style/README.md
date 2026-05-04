# Style Conventions

This directory holds the project's writing and code conventions. Five documents, split by audience:

| Document | Audience | When to read |
|:---|:---|:---|
| `website_docs.md` | Lancet2 user reading the rendered MkDocs site | Authoring or editing a page in `docs/` |
| `code_comments.md` | Developer with the source open | Writing or reviewing inline comments in `src/` |
| `cpp_style.md` | Anyone editing C++ in `src/` | Adding code, hitting a clang-tidy finding, structuring a struct |
| `test_style.md` | Anyone editing C++ in `tests/` | Adding a TEST_CASE, picking a tag, writing a property test, generating a fixture |
| `sync_and_verification.md` | Anyone committing | Cross-document sync rule and the pre-commit verification checklist |

Voice, tone, and audience expectations differ between website docs and code comments — see the focused documents for each. The principles below apply across all four.

## Core Principles

These five rules apply equally to website documentation, inline code comments, C++ source, and dev docs. They are the shared foundation underneath the audience-specific guidance in each focused document.

1. **Source code is the single source of truth.** Every claim in any document must be verifiable against the current codebase. Never document behaviour from memory or assumption. Read the implementation first.

2. **High density, no filler.** Every sentence should teach the reader something they cannot trivially infer. Omit preambles like "In this section, we will discuss..." and start with the content. Filler adverbs ("strictly", "fundamentally", "essentially") rarely earn their place — strike them.

3. **Explain the *why*, not just the *what*.** Stating that `Match = 0` is useless without explaining SIMD overflow prevention. Parameters, thresholds, and design decisions always need rationale. The code shows the *what*; words exist to convey the *why*.

4. **Quantify wherever possible.** Prefer "reduces WGS runtime by ~80%" over "significantly improves performance." Prefer "≥2 reads at the same position" over "multiple reads." When a claim has a number behind it, use the number.

5. **Address the reader's real question.** For every feature, the reader asks: *what does it do, when should I use it, what's the trade-off?* Answer all three. A document that explains only what without when or trade-off is incomplete.

## The synchronization rule

Every claim in code comments, website docs, and dev docs must stay in sync with the actual logic in the codebase. Stale comments and stale docs are worse than no documentation. The full audit procedure for keeping these synchronized lives in `sync_and_verification.md` — read it before any substantive refactor that touches a name, a metric, or a behavioural semantic.

The one exception to the synchronization rule is `docs_dev/investigations/` — investigations are immutable snapshots of a moment in time and are deliberately *not* maintained against the evolving codebase.

## When the rules conflict

If a rule in a focused document contradicts a Core Principle here, the focused document wins for its specific context. The hub principles are the floor; the focused documents add audience-specific rules on top. For example: the "explain *why*, not just *what*" principle holds everywhere, but how that explanation is phrased differs between website docs (where the audience may be a clinician) and code comments (where the audience is a developer with the source open).

## How to use this directory

If you are writing user-facing content, read `website_docs.md` end-to-end before drafting. The voice and structure rules there are tighter than they appear at first glance, and a draft written without reading them tends to need substantial editing before merge.

If you are writing or editing code comments, `code_comments.md` is the reference. The six comment-type templates there (block headers, ASCII diagrams, pipeline annotations, formula comments, data-structure tables, design-decision comments) cover ~90% of comment-writing situations in the project.

If you are touching C++ source, `cpp_style.md` is the rulebook. clang-format and clang-tidy enforce most of it automatically, but several rules — quote-vs-angle-bracket includes, struct member layout, the prefer-`<algorithm>`-and-Abseil rule — are enforced by code review.

If you are adding or editing tests, `test_style.md` is the rulebook for `tests/`. It covers the per-symbol PascalCase tag convention, the `namespace lancet::<layer>::tests` wrapper, the determinism discipline (no `random_device`, no wall clock, no own-filesystem dependencies), the property-test rules, and the test-fixture-script layout (`tests/scripts/<name>.py` produces `tests/data/<layer>/<name>.tsv`).

Before every commit, run through `sync_and_verification.md`'s pre-commit verification checklist. The list is short by design.
