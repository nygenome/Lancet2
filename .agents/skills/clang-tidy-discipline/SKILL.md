---
name: clang-tidy-discipline
description: Use when clang-tidy reports a violation that needs to be addressed, when adding a NOLINT suppression of any kind, when a function approaches the cognitive-complexity ceiling, or when invoking clang-tidy directly. Trigger on "clang-tidy says...", "I need to suppress this warning", "NOLINT", "NOLINTBEGIN", "this function is too complex", "should I refactor or NOLINT", "is it OK to use clang-tidy --fix", "what's the right NOLINT form for this case". Walks the procedure for resolving a clang-tidy violation, the discipline for justifying and shaping NOLINT suppressions, the complexity ceiling that gates further function growth, and the project-specific reason `clang-tidy --fix` is forbidden. Defers to docs_dev/style/cpp_style.md for the canonical rule statements.
allowed-tools: Read, Glob, Grep, Edit
---

# Clang-tidy discipline on Lancet2

`docs_dev/style/cpp_style.md` is the source of truth for clang-tidy rules. This skill is the **procedural how-to**: when a clang-tidy violation appears, how to decide whether to fix or suppress; when adding a NOLINT, what shape the suppression takes and where the rationale lives; when a function approaches the complexity ceiling, what the alternatives are.

The first answer to "should I suppress this?" is always **no — fix the underlying issue**. Reach for a suppression only when the check is genuinely wrong about this specific case, the diagnostic cannot be fixed without harming clarity or correctness, and you can articulate why in writing.

## Step 1 — Diagnose before acting

Read the full clang-tidy output for the violation. Note the check name (e.g., `cppcoreguidelines-pro-type-reinterpret-cast`), the file:line, and the explanation. Before touching code:

1. **Is the check correct about this code?** Read the cited line and surrounding context. Most violations are real.
2. **Can the underlying code be restructured to satisfy the check?** Most can. The check name often suggests the fix (rename to satisfy `readability-identifier-naming`, restructure to satisfy `readability-cognitive-complexity`, etc.).
3. **Is the check genuinely wrong here?** Rare but real: htslib C FFI requiring `reinterpret_cast`, struct members that legitimately must be `const&`, etc. Only at this stage is suppression a candidate.

If you can fix the underlying issue, fix it. Suppression is a last resort.

## Step 2 — When suppression IS the right call: the form

Suppressions must be **scoped**. The forms permitted by Lancet2:

- `// NOLINTNEXTLINE(check-name)` for a single line.
- `// NOLINTBEGIN(check-name)` ... `// NOLINTEND(check-name)` for a block of related suppressions (e.g., a struct's members that must all be `const&`).

The forms **forbidden** by Lancet2:

- Bare `// NOLINT` (no check name) — pollutes the statement line, hides what is being suppressed.
- Inline `// NOLINT(check-name)` (same line as the code) — same scannability problem.
- Inline `// NOLINTNEXTLINE(check-name) -- rationale here` (rationale on the same line as the directive) — rationale must be ABOVE the directive, on a `//` comment line of its own.

The `validate_cpp_identifiers.py` PreToolUse hook hard-blocks bare-NOLINT and missing-rationale additions inside `src/lancet/`.

The structural reason matters: bare inline NOLINT (`int x; // NOLINT(check)`) is fragile under clang-format. When clang-format wraps long lines, it can move the comment away from the statement it suppresses, silently invalidating the suppression. Scoped forms anchor to surrounding lines, not the wrapped statement, so they survive any future formatting pass. This is not stylistic preference — it is the only form robust to the project's auto-format policy.

## Step 3 — When suppression IS the right call: the rationale

Every NOLINT directive must carry a `//` rationale comment on the line(s) immediately above it. The rationale answers two questions: **why is the check wrong about this case**, and **what would change to make the suppression unnecessary**. One sentence is usually enough; longer when the case is subtle.

The canonical form:

```cpp
// htslib's bam1_t accessor requires raw integer arithmetic on a packed
// field; the cppcoreguidelines check is correct in general but the
// pattern is forced by the C library's ABI.
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
auto const* aptr = reinterpret_cast<unsigned long long const*>(first.data());
```

For block suppressions:

```cpp
// These struct members must be const& because the type is non-copyable
// and the consuming SPOA layer requires reference semantics; the
// cppcoreguidelines check would force pointer-or-value semantics that
// breaks the SPOA contract.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
absl::Span<u8 const> const mQuery;         // 16B
absl::Span<u8 const> const mTarget;        // 16B
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
```

`NOLINTEND` carries no rationale (it's a closing token; the rationale lives on the matching `NOLINTBEGIN`).

## Step 4 — Disclose the suppression

Every newly-added NOLINT suppression must be mentioned in the change author's summary of work (the description that accompanies the diff). This is the disclosure mechanism: it prevents suppressions from being quiet workarounds for real violations. The `fresh-reviewer` agent enforces this — a NOLINT in the diff that isn't disclosed in the summary is a finding regardless of whether the suppression itself is justified.

When refactoring, do not carry old suppressions forward without re-evaluation. Remove all NOLINTs in the touched code and re-add only those that survive a fresh review. The accumulation of suppressions across a codebase is a signal worth surfacing in `/audit-bundle`.

## Function complexity ceilings

`docs_dev/style/cpp_style.md` § "Function complexity thresholds" specifies the cognitive-complexity ceiling (35), max parameters (6), and max nesting depth (4). These are clang-tidy-enforced. A function exceeding the ceiling has three resolution paths, in order of preference:

1. **Decompose** — split into smaller helper functions with clear single responsibilities.
2. **Refactor structure** — collapse branching into a table-driven dispatch, lift invariants out of inner loops, replace nested conditions with early returns.
3. **NOLINT, last resort** — only if the function is at a true complexity floor (e.g., a finite-state machine that genuinely has 40 states; a parser dispatch over a fixed grammar) and decomposition would obscure the structure. Rationale must explain why the structure is irreducible.

## Invoking clang-tidy directly

The project provides `pixi run lint-check` (read-only). `pixi run lint-fix` is **deliberately not provided**. Clang-tidy's auto-fix mode (`clang-tidy --fix` or the removed `--fix` flag on `scripts/run_clang_tidy.py`) is **forbidden** for Lancet2 work, including via direct binary invocation.

Why: the auto-fix path has historically broken compilation in this project. `modernize-use-trailing-return-type` rewrote lambdas wrong; `modernize-use-ranges` broke sort calls without `<=>`; `misc-include-cleaner` reordered includes destructively. Silent miscompilations have shipped via `--fix`-edited code. The risk-vs-reward is wrong: clang-tidy violations need human judgment about whether to fix the structure or suppress the diagnostic; mechanical rewriting trades that judgment for speed.

The right invocation pattern is `pixi run lint-check`, read the output, edit by hand, re-run. `pixi run iwyu-fix` (which auto-fixes includes and re-formats) is fine — IWYU's rewriter is conservative and the project has no failure history with it. Clang-tidy auto-fix is the exception.

## When NOT to use this skill

Do not use this skill for clang-tidy violations in vendored code under `cmake-build-*/_deps/`. That code is read-only by protected-paths; the only fix path is at the Lancet2 callsite. Do not use this skill to argue for raising the cognitive-complexity ceiling project-wide; that's a bundle-wide change discussed in `docs_dev/architecture/`.

## References

- `docs_dev/style/cpp_style.md` § "NOLINT usage rules" — canonical rule statements with examples.
- `docs_dev/style/cpp_style.md` § "Function complexity thresholds" — ceiling values and rationale.
- `.clang-tidy` — the actual check configuration.
