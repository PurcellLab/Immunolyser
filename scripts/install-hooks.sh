#!/bin/sh
# Install the pre-commit hook from scripts/pre-commit into .git/hooks/.
# Run once after cloning: sh scripts/install-hooks.sh

set -e

REPO_ROOT="$(git rev-parse --show-toplevel)"
HOOK_SRC="$REPO_ROOT/scripts/pre-commit"
HOOK_DST="$REPO_ROOT/.git/hooks/pre-commit"

cp "$HOOK_SRC" "$HOOK_DST"
chmod +x "$HOOK_DST"
echo "Pre-commit hook installed at $HOOK_DST"
echo "The hook runs: pytest tests/ -x -q before every commit."
