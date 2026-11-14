#!/usr/bin/env bash
# Quad meshing and parameterization viewers. The vector-field viewer lives in
# ../vector-field-decomposition.
#
#   ./build.sh gui-build   compile the remaining viewers
set -euo pipefail

cd "$(dirname "$0")"
MODE="${1:-gui-build}"
if [[ $# -gt 0 ]]; then
  shift
fi

build() {
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
  cmake --build build --parallel "$(sysctl -n hw.logicalcpu)"
}

case "$MODE" in
  gui-build)
    build
    ;;
  *)
    echo "usage: ./build.sh gui-build" >&2
    exit 1
    ;;
esac
