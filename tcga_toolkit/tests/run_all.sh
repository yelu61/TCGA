#!/usr/bin/env bash
# Run every smoke test in this directory. Each test exits non-zero on failure.
# Usage:
#   bash tcga_toolkit/tests/run_all.sh

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.." || exit 1

status=0
echo "== test_introspection.R =="
Rscript tests/test_introspection.R || status=1
echo
echo "== test_templates.R =="
Rscript tests/test_templates.R || status=1
echo
echo "== test_dispatch.R =="
Rscript tests/test_dispatch.R || status=1
echo
echo "== test_python_adapter.py =="
if command -v python3 >/dev/null 2>&1; then
  python3 tests/test_python_adapter.py || status=1
else
  echo "  SKIP  python3 not on PATH"
fi

echo
if [ $status -eq 0 ]; then
  echo "All smoke tests passed."
else
  echo "Some smoke tests failed."
fi
exit $status
