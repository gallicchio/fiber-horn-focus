#!/usr/bin/env bash
# The metal waveguide.  First each horn's own radiation at several resolutions ('emit': TE11 sent
# out through the throat), which by reciprocity gives the coupling of every Gaussian beam, the
# best waist and focus, and the convergence (scripts/horn_design.py).  Then three sets of three
# launches (params.HORN_SETS), each launch also in air for the incident field, with the turn-on
# movie frames:
#   smooth        the smooth horn, the fiber's three beams focused on its mouth
#   smooth_focus  the smooth horn, the beams around its best waist, focused where that is best
#   corrugated    the corrugated horn, the same for its best waist and focus
# On a 16-core desktop: ~3 h for the design runs and ~12 h per set at 120 px/wavelength (the six
# runs of a set go side by side, NP_EACH processes each).
set -euo pipefail
source "$(dirname "$0")/env.sh"
horn() { "$MPIRUN" -np "$1" "$PYTHON" -m fiberlaunch.horn "${@:2}"; }
d="$DATA/horn"
if [ -z "${SKIP_DESIGN:-}" ]; then
  for r in 20 40 60 80 120; do
    [ -f "$d/design/emit_res$r/meta.json" ] || horn "$NP" --kind emit --res "$r" --half-width 6 --out "$d/design/emit_res$r"
  done
  for r in 40 60 80 120; do
    out="$d/design/corrugated_emit_res$r"
    # The slots' edges fall on cell boundaries, where Meep's averaging depends on where the box
    # starts; at 120 the box starts where the launches' does, so the slots are the same cells.
    zs=()
    [ "$r" = 120 ] && zs=(--z-source -8)
    [ -f "$out/meta.json" ] || horn "$NP" --kind emit --wall corrugated --res "$r" --half-width 6 "${zs[@]}" --out "$out"
  done
  "$PYTHON" scripts/horn_design.py "$d"/design/emit_res{20,40,60,80,120} | tee media/horn_design.txt
  "$PYTHON" scripts/horn_design.py "$d"/design/corrugated_emit_res{40,60,80,120} | tee media/horn_design_corrugated.txt
fi
res=${RES:-120}
for set in ${SETS:-smooth smooth_focus corrugated}; do
  mkdir -p "$d/$set/res$res"
  for case in ${CASES:-narrow matched wide}; do
    for kind in air horn; do
      out="$d/$set/res$res/$case"
      [ "$kind" = air ] && out="${out}_air"
      extra=()
      [ "$kind" = horn ] && extra=(--movie)
      [ -f "$out/meta.json" ] || horn "${NP_EACH:-5}" --set "$set" --kind "$kind" --case "$case" --res "$res" \
          "${extra[@]}" --out "$out" > "$out.log" 2>&1 &
    done
  done
  wait
done
