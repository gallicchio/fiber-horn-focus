"""Was the field steady before it was Fourier-transformed?  From the movie frames (Ex on
the x-z plane, every 0.34 periods), fit a sinusoid over successive 8-period windows at
points along the fiber and report how much its complex amplitude changes over the last
10 periods.  Points where the field is below 10^-3 of the on-axis field are skipped: there
the frames' float32 rounding, not the physics, sets the change.

    python scripts/check_steady.py [--root data/main/res16]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def amplitude(s, t):
    """Least-squares fit of s = Re(A exp(-i 2 pi f t)); the frames are not a whole number
    per period, so a plain projection would leak."""
    M = np.stack([np.cos(2 * np.pi * P.FREQ * t), np.sin(2 * np.pi * P.FREQ * t)], 1)
    c, *_ = np.linalg.lstsq(M, s, rcond=None)
    return c[0] - 1j * c[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    a = ap.parse_args()
    for name in a.cases:
        if not os.path.exists(os.path.join(a.root, name, "meta.json")):
            continue
        run = load_case(a.root, name)
        F, t = run.frames()
        F = F[:len(t)]  # the frame file has one spare slot at the end
        side = run.side("xz")
        x, z = side["x"], side["z"]
        worst, n = 0.0, 0
        late, last = [(t > t1 - 8) & (t <= t1) for t1 in (t[-1] - 10, t[-1])]
        for zz in (10, 30, 60, 90, 115):
            iz = int(np.argmin(abs(z - zz)))
            axis = None
            for xx in (0.0, 3.0, 6.0, 10.0):
                s = np.asarray(F[:, int(np.argmin(abs(x - xx))), iz], dtype=float)
                A0, A1 = amplitude(s[late], t[late]), amplitude(s[last], t[last])
                axis = abs(A1) if axis is None else axis
                if abs(A1) >= 1e-3 * axis:
                    worst, n = max(worst, abs(A1 - A0) / abs(A1)), n + 1
        print(f"{name:8s} largest relative change of the field over the last 10 periods "
              f"({n} points, z = 10 to 115, x = 0 to 10): {worst:.1e}")


if __name__ == "__main__":
    main()
