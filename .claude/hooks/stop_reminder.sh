#!/bin/bash
# Lancet2 Stop hook: stop_reminder
#
# Claude Code fires the Stop hook after EVERY assistant turn, not at
# session end. Running heavy validation here (the previous behavior)
# meant a 5-15 min build+test+iwyu+lint cycle after every Q&A turn,
# regardless of whether C++ files changed. The new design moves all
# heavy work to the pre_commit_gate hook (gates `git commit` only) and
# the /check slash command (explicit user invocation). This hook is
# now a quiet reminder: silent if no relevant files dirty, one-liner
# otherwise.

set +e

relevant=$(git status --porcelain -- \
    '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' '*.hxx' \
    'CMakeLists.txt' '**/CMakeLists.txt' '*.cmake' '**/*.cmake' \
    2>/dev/null)

if [ -n "$relevant" ]; then
    count=$(printf '%s\n' "$relevant" | wc -l | tr -d ' ')
    echo "⚠ $count unstaged/staged C++/CMake change(s); run /check or commit to validate"
fi

exit 0
