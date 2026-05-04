# Hooks

Hooks are deterministic. They fire every time, can block actions before they execute (`PreToolUse` with exit code 2), and run outside Claude's context window so they cost zero tokens. Skills and `AGENTS.md` are probabilistic — Claude can ignore or forget them. Hooks cannot be ignored.

This is the entire reason hooks exist in this bundle: there are rules that must be enforced even when Claude has been told to ignore them, even when the rule is buried under context, even when Claude is focused on something else. The eight hooks here are the rules that meet that bar.

## What's here

### Blocking hooks (PreToolUse, exit 2 on violation)

- **`block_protected_paths.py`** — refuses edits to `pixi.toml`, `cmake/`, `.github/`, `data/`, `pixi.lock`, `Dockerfile`, `cmake-build-*/`, `_deps/`, `.pixi/`, `CHANGELOG.md`. Substring match. The blocklist is conservative; `notes/scratch/` and `notes/<feature>/` are NOT blocked.
- **`block_dangerous_bash.py`** — validates bash commands before execution. Three layers: (1) **blanket-banned patterns** (`git push`, `git reset --hard origin/...`, `dd if=`, `mkfs.<fs>`, the fork-bomb literal, `chmod -R 777 /`, `chown -R`, `xargs rm`) — these have no safe variant; the maintainer handles them by hand outside Claude Code. `xargs rm` is in this list because its paths come from stdin and cannot be statically validated. (2) **path-validated deletions** (`rm`, `rmdir`, `shred`, `find … -delete`, `find … -exec rm`) — the hook parses out the path arguments and verifies each is absolute, contains no shell metacharacters that would re-expand at execution time (`$`, backtick, `*`, `?`, `[`, `]`, `{`, `}`), resolves cleanly (no `..` escape, no symlink to outside the allowlist), is not the bare allowlist root, and lands under `/tmp/`, `/var/tmp/`, `/scratch/`, OR inside `$CLAUDE_PROJECT_DIR` AND is gitignored (verified via `git check-ignore`). Any path that fails validation blocks the whole command. (3) **user-approval prompt** (downstream of the hook) — `settings.json` does NOT list `Bash(rm:*)` / `Bash(rmdir:*)` / `Bash(shred:*)` in either allow or deny, so Claude Code prompts the user per invocation. The hook prints a verdict line ("safe: under /tmp/", "safe: gitignored under project root") to stderr so the user-approval prompt surfaces the categorization the hook applied. Together this forms a two-key gate: the hook validates path safety AND the user approves each delete by hand. Word-boundary regex matching on the blanket-banned layer avoids false positives on words containing `rm` (firmware, storm, warm).
- **`validate_layer_direction.py`** — enforces the six-layer dependency rule from CMakeLists.txt: base → hts → cbdg → caller → core → cli. Blocks any include that goes against the chain.
- **`validate_cpp_identifiers.py`** — catches `mPascalCase` violations on members, `using namespace std` in headers, bare `assert()`, `std::format` and `std::print` (the project uses fmtlib via spdlog), and bare `// NOLINT` (any inline same-line form, with or without a check name; only scoped `NOLINTNEXTLINE` and `NOLINTBEGIN`/`NOLINTEND` with rationale on the line(s) above are permitted).
- **`validate_commit_message.py`** — block-or-pass validator for `git commit -m` messages, grounded in `.chglog/config.yml` and `docs_dev/style/cpp_style.md` § Git commit messages. Hard-blocks invalid types, malformed subjects, and substantive changes (>30 staged lines) without a body. Expands `$(cat <<'EOF' ... EOF)` HEREDOCs so the parsed text reflects what git actually receives. Never emits soft warnings.
- **`pre_commit_gate.sh`** — verifies that `/fix-and-validate` has passed against the EXACT working-tree state about to be committed. Does NOT re-run validation. Compares `git stash create` (an atomic commit-object hash representing the full working-tree state) against a marker at `.claude/cache/validation-state.txt` written by `/fix-and-validate` on full success. Match → silent exit 0 (~50ms). Mismatch or missing marker → exit 2 with remediation pointing at `/fix-and-validate`. Skips silently when the staged diff has no entries matching C++ source, CMake files, build-helper scripts (`scripts/*.{py,sh}`), or `.clang-tidy`/`.clang-format` — those commits don't need `/fix-and-validate` and don't trigger the marker check. The actual heavy validation (Release build + lint + IWYU + test, 3-7 min warm) runs only inside `/fix-and-validate`, on demand; the gate itself is hash-comparison fast.

### Informational hooks (always exit 0)

- **`pre_commit_summary.sh`** — prints a one-paragraph summary of the staged change (file count, lines, layers touched, kinds of files) before a commit executes. Heuristic warnings for multi-layer-without-docs, source-without-tests, and VCF schema changes. Never blocks.
- **`stop_reminder.sh`** — Stop-event hook that prints one line if the working tree has unstaged/staged C++/CMake changes (and is silent otherwise). The Stop event fires after every assistant turn (not just at session end), so the hook is intentionally lightweight: a reminder, not a validation gauntlet. Heavy validation lives in `/fix-and-validate`, invoked deliberately by the user. Runs in ~50ms regardless of state.

The Lancet2 statusline (`.claude/ui/statusline.sh`) is wired through the separate `statusLine` settings slot, not a hook event; it is configured alongside the hooks for organizational simplicity but is not part of the hook lifecycle.

## Why hooks vs skills vs AGENTS.md

A rule is hookable when (a) it has zero false positives — every match is a real violation — and (b) you want it enforced regardless of what Claude is currently focused on.

- **Layer-direction violations:** hookable (the dependency chain is in CMakeLists, every upward include is wrong).
- **Naming violations** like bare `assert()` or `using namespace std` in headers: hookable (the rules are crisp, the source enforces them via clang-tidy).
- **Comment language quality:** not hookable (judgment calls about whether a word is filler in this specific context).
- **Choice of which algorithm to use:** not hookable (judgment about clarity vs. cleverness).

A rule that has any judgment call belongs in a skill or in AGENTS.md, not a hook. Hooks that try to encode judgment produce false positives that train Claude (and you) to ignore the hook output.

## Why use hooks at all when AGENTS.md could state the rule

AGENTS.md is read at session start (loaded by Claude Code through the CLAUDE.md wrapper) and at the start of each message, but it's not enforced — Claude can be persuaded out of any rule there by a strong enough context signal ("but in this case I should ignore the rule because…"). Hooks override the persuasion.

For the layer rule specifically: AGENTS.md explains the architecture so Claude understands the design; the hook prevents accidental violations even when Claude has correctly understood the design but slipped while writing.

## Why use hooks at all when CI catches violations

Feedback latency. CI catches the violation after you push; the hook catches it before the edit is even applied. The faster the feedback, the smaller the rework cost.

## Hook calibration philosophy

Every Lancet2 hook is block-or-pass. Hooks emit stderr only when blocking; they do not emit soft-warnings on commits that pass. A check that would fire as a heuristic — flagging a stylistic concern, suggesting an alternative, noting an edge case — does not belong in a hook. Heuristic checks belong in a skill or in `AGENTS.md` where Claude can apply judgment. Soft-warnings train both the user and Claude to ignore stderr output, and the cost of that training (a real block missed because the user has stopped reading hook output) is much higher than the value of any heuristic warning.

The bar for hard blocks: (a) the violation is explicit and unambiguous in the source text, (b) reviewers have flagged it repeatedly in the past, and (c) the cost of an erroneous block is trivially recoverable (the user adds an inline comment override or fixes the offending pattern and re-runs).

When you propose a new check, classify it against this bar. If you're not sure whether the pattern is FP-near-zero, don't ship it as a hook. Skills, `AGENTS.md`, and CI gates are all valid alternatives.

## The hook ceiling

The hook count is technical debt. Every hook has to be maintained as the codebase changes (the layer-direction hook needs updates when layer rules change; the naming hook needs updates when conventions evolve; the gate scripts need updates when pixi tasks evolve; etc.). Eight is a comfortable working count for this project. Adding a ninth hook should require articulating why the existing eight aren't enough — and should ideally pair with retiring one that's no longer earning its keep.

The eight break down as: three deterministic source-edit enforcement (`block_protected_paths`, `validate_layer_direction`, `validate_cpp_identifiers`), one Bash-safety enforcement (`block_dangerous_bash`), three commit-time gates (`validate_commit_message`, `pre_commit_summary`, `pre_commit_gate`), and one Stop-time reminder (`stop_reminder`).

Specific anti-patterns to avoid:

- **A hook that duplicates a skill.** The hook fires deterministically; the skill is suggestive. If both fire, the skill's effort is wasted. Pick one.
- **A hook whose rule changes frequently.** Frequent edits to a hook's regex or path list means the rule has soft edges that the regex can't capture. Convert to a skill (where Claude can apply judgment) or remove.
- **A hook with a long allowlist of exceptions.** Exceptions are how false positives leak back in. If the hook needs more than a handful of exceptions, the rule isn't actually deterministic.

## Maintenance lifecycle

### Adding a new hook

A new hook earns its place when:

1. The rule is deterministic — every match is a genuine violation.
2. The rule needs to be enforced rather than suggested. A skill or AGENTS.md note isn't enough.
3. CI doesn't catch it tightly enough. The feedback-latency improvement justifies the maintenance cost.

Procedure: write the script (Python or bash, both supported). Convention: the script reads the JSON payload from stdin, exits 0 to allow, exits 2 to block. Print clear error messages to stderr — Claude reads stderr and uses it to explain the violation to the user. Wire the hook into `settings.json` under the appropriate event matcher (`PreToolUse` for blocking, `Stop` for end-of-turn reminders). Add to this README's "What's here" section. Write a smoke test under `notes/scratch/` for non-trivial logic; the commit-message validator was multi-case smoke-tested at every behavior change.

The hook event lifecycle is documented in Anthropic's hooks reference.[^hooks-docs] The most useful events are `PreToolUse` (can block; happens before the tool runs) and `Stop` (end-of-turn; observational). `PostToolUse`, `SessionStart`, `UserPromptSubmit`, and `PreCompact` exist but are not used in this project.

[^hooks-docs]: <https://code.claude.com/docs/en/hooks>

### Deleting a hook

Delete a hook when:

- It has not blocked anything in the last quarter (it's not catching violations because the rule is internalized, or because the violation never happens).
- Its rule has been superseded by a clang-tidy check, a CMake check, or another deterministic enforcement.
- It has accumulated false positives or exceptions to the point where you regularly bypass it.

Procedure: remove the script, remove the entry from `settings.json`, update this README, grep `.claude/**` and `AGENTS.md`/`docs_dev/**`/`docs/**` for any reference to the script's name and update those too. Git history preserves the script if you ever need it back.

### Refactoring a hook

Common patterns:

**Tightening the regex.** A hook that fires false positives needs its match logic constrained. The fix is usually a more specific pattern (e.g., a path prefix rather than a substring, or a word-boundary regex rather than a substring). Resist the urge to add an allowlist of exceptions; that's a smell that the rule isn't deterministic.

**Splitting one hook into two.** When a hook's logic has grown to handle two distinct rules, splitting often makes both clearer.

**Adding fast-path skips.** A hook that fires on every Bash command but only does work in a specific subset (e.g., `git commit`) should short-circuit early — both `validate_commit_message.py` and `pre_commit_gate.sh` follow this pattern. The fast-path keeps the per-command cost near-zero for the 99% of bash invocations the hook doesn't need to inspect.

**Adding informational output.** Many hooks benefit from richer stderr output without changing their block/allow logic. The `pre_commit_summary` hook is informational-only; the validation hooks include diagnostic context (specific line numbers, the offending pattern, a suggested fix) so Claude can explain the violation to the user accurately. When a hook blocks, the stderr text is what Claude has to work with — invest in it.

### Reviewing the hook set

Quarterly. Walk all hooks:

- Has this hook fired in the last quarter? Did it catch real violations?
- Is the rule still deterministic? Have any judgment calls crept into the regex?
- Does the error message still make sense? Does Claude know how to explain it to the user?
- Could this hook be deleted (CI now catches it; clang-tidy enforces it; the rule was internalized)?

The grounding for what each hook checks should match the actual source. The layer-direction hook's chain is grounded in CMakeLists.txt; if the layer chain changes, the hook needs updating. The naming hook's rules are grounded in `.clang-tidy`; if those change, the hook needs updating. The pre-commit gate's step list is grounded in pixi.toml; if pixi tasks change, the gate needs updating. Periodic re-grounding catches drift.

### Retiring a hook

Retire a hook that has been on the deletion candidate list for two quarters running. Hooks are technical debt; the bias should be toward fewer hooks rather than more.

The exception is the `pre_commit_summary` and `stop_reminder` hooks — they don't enforce anything, so they don't carry the false-positive risk. Informational hooks are cheap to keep around if they're useful at all.

## Cost model

Hooks cost zero context tokens. They run outside the model's context window entirely. The cost they do have is wall-clock and maintenance — every hook is code that has to be kept in sync with the rules it enforces. See `../cost-model.md` for how this fits into the overall mechanism cost picture.
