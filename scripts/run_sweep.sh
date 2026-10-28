#!/usr/bin/env bash
# Coupling efficiency against beam NA: nine waists from 4x too small to 4x too big, each
# launched into a short (10 wavelength) fiber: the guided power is fixed from the facet on,
# so a short fiber measures it as well as a long one.  ~2 h at 16 px/wavelength.
set -euo pipefail
source "$(dirname "$0")/env.sh"
res=${RES:-16}
hw=16
d="$DATA/sweep/res$res"
for f in 0.25 0.3333 0.5 0.6667 1 1.5 2 3 4; do
  w0=$("$PYTHON" -c "from fiberlaunch import params as P; print(f'{P.W_MATCHED*$f:.4f}')")
  [ -f "$d/w$w0/air/meta.json" ] || fdtd --kind air --w0 "$w0" --res "$res" --length 10 --half-width $hw --out "$d/w$w0/air"
  [ -f "$d/w$w0/fiber/meta.json" ] || fdtd --kind fiber --w0 "$w0" --res "$res" --length 10 --half-width $hw \
      --mode-planes 5 8 --out "$d/w$w0/fiber"
done
