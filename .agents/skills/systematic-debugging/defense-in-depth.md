# Defense-in-Depth Validation

## Overview

When you fix a bug caused by invalid data, adding validation at one place feels sufficient. But that single check can be bypassed by different code paths, refactoring, or mocks.

**Core principle:** validate at EVERY tier the data passes through. Make the bug structurally impossible.

> A note on terminology: this file uses **"tier"** (not "layer") for the four defensive checkpoints. The codebase already uses "layer" for the architectural stack (`base/` → `hts/` → `cbdg/` → `caller/` → `core/` → `cli/`), and conflating the two would be confusing.

## Why Multiple Tiers

Single validation: "We fixed the bug."
Multiple tiers: "We made the bug impossible."

Different tiers catch different cases:
- **Entry-point validation** catches most bugs at the API boundary.
- **Business-logic validation** catches edge cases the entry-point check missed.
- **Environment guards** prevent dangerous operations in specific contexts (debug-only assertions, sanitizer-only checks).
- **Debug instrumentation** captures forensic context when the other tiers somehow fail.

## The Four Tiers

### Tier 1: Entry-Point Validation

**Purpose:** reject obviously invalid input at the public API boundary.

In Lancet2 idiom, this typically takes the form of:
- Parameter validation in the constructor or the function that the rest of the codebase calls into.
- Returning `absl::Status` / `absl::StatusOr<T>` for recoverable errors with a precise diagnostic.
- Throwing for programmer errors that should never reach the entry point in correctly-written calling code.

The check should be specific enough that the caller can fix the issue from the error alone, without re-reading the callee's source.

### Tier 2: Business-Logic Validation

**Purpose:** ensure the data makes sense for *this specific operation*, not just for the type system.

Even if Tier 1 lets a value through (e.g., a non-empty path), Tier 2 may know constraints Tier 1 doesn't (e.g., this operation requires a directory, not a file; this region must lie within the reference; this k must be odd and within the configured range).

The check at this tier should reference the operation's invariants explicitly, ideally as a `LANCET_ASSERT` (defined in `src/lancet/base/assert.h` — compiles to a no-op in Release while still preserving the documentation effect of the assertion in source).

### Tier 3: Environment Guards

**Purpose:** prevent dangerous operations in contexts where they shouldn't run.

Examples in Lancet2:
- Sanitizer-build-only checks gated on the `LANCET_SANITIZE_BUILD` macro that confirm invariants which would be too expensive to check in Release.
- Test-only assertions that catch test-environment misuse (e.g., a test fixture passing an empty path that would silently degrade the operation in production).
- Concurrency contracts annotated with `ABSL_GUARDED_BY(mutex)` so clang's thread-safety analyzer rejects unguarded access at compile time. This is a static-analysis tier that doesn't even reach runtime — the bug is caught before the binary builds.

### Tier 4: Debug Instrumentation

**Purpose:** capture context when the other three tiers somehow fail.

Use `LOG_DEBUG` to record the operation's inputs and the relevant state immediately before the dangerous step. In a debug-build crash, the last log line tells the next developer where to start. Strip or guard heavy instrumentation behind a debug-only macro before merging.

## Applying the Pattern

When you find a bug:

1. **Trace the data flow.** Where does the bad value originate? Where is it used?
2. **Map all checkpoints.** List every tier the value passes through between origin and use.
3. **Add validation at each tier.** Pick the right kind for the tier (entry-point status return, business-logic assertion, environment guard, debug log).
4. **Test each tier in isolation.** Try to bypass Tier 1 — does Tier 2 catch it? Try to misuse the test fixture — does Tier 3 catch it? Without per-tier tests, you don't actually know which tier is doing the work.

## Key Insight

All four tiers serve different cases. During refactoring or when adding new code paths, each tier catches bugs the others miss:

- Different code paths bypass entry-point validation — business-logic asserts catch them.
- Test mocks bypass business-logic checks — environment guards catch them.
- Edge cases on different platforms or build modes need environment-specific guards.
- Debug logging identifies structural misuse that no static check anticipated.

**Don't stop at one validation point.** Add checks at every tier.

## A note on cost

Tiers 1, 2, and 3 are essentially free in production: `absl::Status` returns are cheap, `LANCET_ASSERT` compiles out, `ABSL_GUARDED_BY` is compile-time only. Tier 4 (debug logging) has a cost and should be guarded. The discipline is not "add expensive checks everywhere" — it's "add the right kind of check at each tier where the data is observable."
