#!/usr/bin/env bash
# The metal waveguide.  First the horn's own radiation at five resolutions ('emit': TE11 sent
# out through the throat), which by reciprocity gives the coupling of every Gaussian beam and
# shows the convergence; then the three launches, each also in air for the incident field,
# with the turn-on movie frames.  On a 16-core desktop: ~2 h for the design runs and ~4 h for
# the launches at 120 px/wavelength (the six launch runs go side by side, NP_EACH processes each).
set -euo pipefail
source "$(dirname "$0")/env.sh"
horn() { "$MPIRUN" -np "$1" "$PYTHON" -m fiberlaunch.horn "${@:2}"; }
d="$DATA/horn"
if [ -z "${SKIP_DESIGN:-}" ]; then
  for r in 20 40 60 80 120; do
    [ -f "$d/design/emit_res$r/meta.json" ] || horn "$NP" --kind emit --res "$r" --half-width 6 --out "$d/design/emit_res$r"
  done
  "$PYTHON" scripts/horn_design.py "$d"/design/emit_res{20,40,60,80,120} | tee media/horn_design.txt
fi
res=${RES:-120}
mkdir -p "$d/res$res"
for case in ${CASES:-narrow matched wide}; do
  for kind in air horn; do
    out="$d/res$res/$case"
    [ "$kind" = air ] && out="${out}_air"
    extra=()
    [ "$kind" = horn ] && extra=(--movie)
    [ -f "$out/meta.json" ] || horn "${NP_EACH:-5}" --kind "$kind" --case "$case" --res "$res" "${extra[@]}" \
        --out "$out" > "$out.log" 2>&1 &
  done
done
wait
