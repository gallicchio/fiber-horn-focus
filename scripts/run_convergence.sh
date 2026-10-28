#!/usr/bin/env bash
# Resolution convergence: the three launches into a short (12 wavelength)
# fiber at 10, 12, 16 and 20 pixels per wavelength.  ~2.5 h on a 16-core desktop.
set -euo pipefail
source "$(dirname "$0")/env.sh"
for res in ${RESOLUTIONS:-10 12 16 20}; do
  d="$DATA/convergence/res$res"
  [ -f "$d/mode/meta.json" ] || fdtd --kind mode --res "$res" --out "$d/mode"
  for case in narrow matched wide; do
    [ -f "$d/${case}_air/meta.json" ] || fdtd --kind air --case "$case" --res "$res" --length 12 --out "$d/${case}_air"
    [ -f "$d/${case}/meta.json" ] || fdtd --kind fiber --case "$case" --res "$res" --length 12 \
        --mode-planes 6 10 --out "$d/$case"
  done
done
