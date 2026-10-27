#!/usr/bin/env bash
# Native principal-curvature field viewer, using this repo's Polyscope.
#
#   ./build.sh gui-build   compile only -> build/bin/main
#   ./build.sh gui         compile and open bunny.obj
set -euo pipefail

cd "$(dirname "$0")"
MODE="${1:-gui}"
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
  gui)
    build
    exec ./build/bin/main "$@"
    ;;
  *)
    echo "usage: ./build.sh [gui-build|gui] [mesh.obj]" >&2
    exit 1
    ;;
esac
