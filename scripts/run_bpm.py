"""Follow the light far down the fiber with the BPM, starting from the 3D FDTD field.

1. validation: BPM from the FDTD field at z = 5 lambda (close to the facet, so even the
   steepest uncoupled light is still inside the FDTD window), compared with the FDTD
   itself at the end of its box (guided power, total power, field correlation);
2. the long run: the same start, through the whole cladding and its coating, to 6000
   lambda (about 3 cm of SMF-28 once the scaling is undone).

    python scripts/run_bpm.py [--root data/main/res16] [--out data/bpm] [--length 6000]
"""

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import bpm as B  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def fdtd_plane(run, z_want):
    x, y, z, mm = run.planes()
    i = int(np.argmin(abs(z - z_want)))
    pb = run.plane_budget()
    return x, y, z[i], np.asarray(mm["Ex"][i], dtype=np.complex128), pb["total"][i], pb["guided"][i], i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "bpm"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    ap.add_argument("--length", type=float, default=6000.0)
    ap.add_argument("--z-start", type=float, default=5.0)
    ap.add_argument("--snap-every", type=float, default=25.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    bpm = B.BPM()
    print(f"BPM grid {bpm.n}^2, dx {bpm.dx:.4f}, cladding radius {B.R_CLAD:.2f}, coating to {B.R_COAT:.2f}")
    path = os.path.join(a.out, "summary.json")
    summary = json.load(open(path)) if os.path.exists(path) else {}  # keep the launches not rerun
    for name in a.cases:
        if not os.path.exists(os.path.join(a.root, name, "xy_coords.npz")):
            print("no FDTD planes for", name)
            continue
        run = load_case(a.root, name)
        x, y, z0, Ex, p0, g0, i0 = fdtd_plane(run, a.z_start)
        u0 = bpm.from_fdtd(x, y, Ex, z0, p0, guided=g0)
        # ---- validation against the FDTD at the far end of its box -----------------------
        xe, ye, z1, Ex1, p1, g1, i1 = fdtd_plane(run, 115.0)
        val = bpm.run(u0.copy(), z1 - z0, z0=z0, sample_every=z1 - z0)
        u1 = bpm.from_fdtd(xe, ye, Ex1, z1, p1)
        # compare inside the FDTD's window only: beyond it the FDTD light is already in its PML
        X, Y = np.meshgrid(bpm.x, bpm.x, indexing="ij")
        win = (abs(X) <= xe.max()) & (abs(Y) <= ye.max())
        ub = val["final"] * win
        corr = abs(np.vdot(u1, ub)) ** 2 / (np.vdot(u1, u1).real * np.vdot(ub, ub).real)
        bpm_total_win = float(np.sum(abs(ub) ** 2) * bpm.dA)
        pb = run.plane_budget()
        v = dict(z_start=float(z0), z_end=float(z1),
                 fdtd_guided_start=float(pb["guided"][i0]), bpm_guided_start=float(abs(bpm.guided_amp(u0)) ** 2),
                 fdtd_guided_end=float(pb["guided"][i1]), bpm_guided_end=float(val["guided"][-1]),
                 fdtd_total_end=float(pb["total"][i1]), bpm_total_end_in_fdtd_window=bpm_total_win,
                 field_correlation_end=float(corr))
        print(name, json.dumps(v))
        # ---- the long run -------------------------------------------------------------------
        t0 = time.time()
        out = bpm.run(u0, a.length - z0, z0=z0, sample_every=5.0, snap_every=a.snap_every)
        print(f"{name}: {a.length:.0f} lambda in {time.time() - t0:.0f} s; guided {out['guided'][0]:.4f} -> "
              f"{out['guided'][-1]:.4f}, total {out['total'][0]:.4f} -> {out['total'][-1]:.4f}")
        out.pop("final")
        np.savez_compressed(os.path.join(a.out, f"{name}.npz"), R=run.R, r_clad=B.R_CLAD, r_coat=B.R_COAT,
                            length_to_smf28=B.LENGTH_TO_SMF28, **out)
        summary[name] = dict(validation=v, guided_end=float(out["guided"][-1]), total_end=float(out["total"][-1]),
                             in_glass_end=float(out["in_glass"][-1]), R=float(run.R))
    with open(path, "w") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
