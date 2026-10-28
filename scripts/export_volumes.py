"""Export the guided and rest intensity volumes of each launch for 3D rendering, all on one
brightness scale (per unit of incident power).

Runs in the simulation environment (it uses the plane stacks and the mode projection)
and writes small float16 volumes that scripts/render_3d.py (render environment) reads.

    python scripts/export_volumes.py [--root data/main/res16] [--out data/volumes]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402
from movie_flythrough import Planes, share_ref  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "volumes"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    ap.add_argument("--xlim", type=float, default=12.0)
    ap.add_argument("--stride", type=int, default=2, help="transverse downsampling")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stacks = share_ref({name: Planes(load_case(a.root, name), a.xlim) for name in a.cases})
    for name, st in stacks.items():
        k = a.stride
        n = (st.keep.sum() // k) * k
        gs, rs = [], []
        for i in range(len(st.z)):
            _, g, r = st.split(i)
            gs.append(g[:n, :n].reshape(n // k, k, n // k, k).mean((1, 3)))
            rs.append(r[:n, :n].reshape(n // k, k, n // k, k).mean((1, 3)))
        x = st.x[st.keep][:n].reshape(-1, k).mean(1)
        np.savez_compressed(os.path.join(a.out, f"{name}.npz"), x=x, z=st.z, ref=st.ref,
                            guided=(np.array(gs) / st.ref).astype(np.float16),
                            rest=(np.array(rs) / st.ref).astype(np.float16))
        print(name, np.array(gs).shape)


if __name__ == "__main__":
    main()
