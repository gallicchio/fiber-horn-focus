#!/usr/bin/env bash
# The three launches into 120 wavelengths of fiber, with the turn-on movie frames and the
# stack of transverse planes.  ~2.5 h per launch at 16 px/wavelength on a 16-core desktop.
set -euo pipefail
source "$(dirname "$0")/env.sh"
res=${RES:-16}
d="$DATA/main/res$res"
[ -f "$d/mode/meta.json" ] || fdtd --kind mode --res "$res" --out "$d/mode"
for case in ${CASES:-narrow matched wide}; do
  [ -f "$d/${case}_air/meta.json" ] || fdtd --kind air --case "$case" --res "$res" --out "$d/${case}_air"
  [ -f "$d/${case}/meta.json" ] || fdtd --kind fiber --case "$case" --res "$res" --movie --xy-planes --out "$d/$case"
done
