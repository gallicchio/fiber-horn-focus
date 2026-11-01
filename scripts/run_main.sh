#!/usr/bin/env bash
# The three launches into 120 wavelengths of fiber, with the turn-on movie frames and the
# stack of transverse planes; the source 8 wavelengths out (params.Z_SOURCE_FAR), so that the
# pictures show the beam converging.  ~3 h per launch at 16 px/wavelength on a 16-core desktop.
set -euo pipefail
source "$(dirname "$0")/env.sh"
res=${RES:-16}
d="$DATA/main/res$res"
zs=$("$PYTHON" -c "from fiberlaunch import params as P; print(P.Z_SOURCE_FAR)")
[ -f "$d/mode/meta.json" ] || fdtd --kind mode --res "$res" --out "$d/mode"
for case in ${CASES:-narrow matched wide}; do
  [ -f "$d/${case}_air/meta.json" ] || fdtd --kind air --case "$case" --res "$res" --z-source "$zs" --out "$d/${case}_air"
  [ -f "$d/${case}/meta.json" ] || fdtd --kind fiber --case "$case" --res "$res" --z-source "$zs" --movie --xy-planes \
      --out "$d/$case"
done
