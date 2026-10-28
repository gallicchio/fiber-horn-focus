"""Table of the convergence study: reflection and coupling against resolution."""

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import fibermode  # noqa: E402
from fiberlaunch import params as P  # noqa: E402

DATA = os.environ.get("DATA", os.path.join(os.path.dirname(__file__), "..", "data"))


def main():
    rows = []
    for d in sorted(glob.glob(os.path.join(DATA, "convergence", "res*")), key=lambda s: int(s.split("res")[-1])):
        res = int(d.split("res")[-1])
        for case in P.CASE_ORDER:
            try:
                air = json.load(open(os.path.join(d, f"{case}_air", "meta.json")))
                fib = json.load(open(os.path.join(d, case, "meta.json")))
            except FileNotFoundError:
                continue
            pin = air["P_net_refl_plane"]
            last = fib["mode_planes"][-1]
            rows.append(dict(res=res, case=case, R=1 - fib["P_net_refl_plane"] / pin,
                             eta=(last["P_fwd"] - last["P_bwd"]) / pin, bwd=last["P_bwd"] / last["P_fwd"],
                             neff=last["neff"], wall=fib["wall_seconds"]))
    m = fibermode.lp01()
    print(f"theory: Fresnel R(n_clad) = {1 - fibermode.fresnel_T(P.N_CLAD):.4f}, "
          f"R(n_eff) = {1 - fibermode.fresnel_T(m['neff']):.4f}, HE11 n_eff = {fibermode.he11_neff():.6f}")
    for case in P.CASE_ORDER:
        th = fibermode.gaussian_overlap(P.CASES[case]["w0"])
        print(f"{case}: scalar overlap {th:.4f}, x Fresnel(n_eff) {th * fibermode.fresnel_T(m['neff']):.4f}")
    print(f"{'res':>4} {'case':>8} {'R':>8} {'eta':>8} {'|a-|^2/|a+|^2':>14} {'n_eff':>10} {'wall s':>7}")
    for r in rows:
        print(f"{r['res']:4d} {r['case']:>8} {r['R']:8.5f} {r['eta']:8.5f} {r['bwd']:14.2e} {r['neff']:10.6f} {r['wall']:7.0f}")
    # Richardson: the facet's staircase error is second order, X(res) = X_inf + c / res^2
    out = {}
    print("extrapolated to infinite resolution (fit X = X_inf + c / res^2):")
    for case in P.CASE_ORDER:
        rr = [r for r in rows if r["case"] == case]
        if len(rr) < 2:
            continue
        h2 = np.array([1.0 / r["res"] ** 2 for r in rr])
        fit = {}
        for key in ("R", "eta"):
            c, x_inf = np.polyfit(h2, [r[key] for r in rr], 1)
            fit[key] = float(x_inf)
            fit[key + "_res16"] = next((r[key] for r in rr if r["res"] == 16), None)
        fit["resolutions"] = [r["res"] for r in rr]
        out[case] = fit
        print(f"  {case:>8}: R {fit['R']:.4f}   eta {fit['eta']:.4f}   (from res {fit['resolutions']})")
    with open(os.path.join(os.path.dirname(__file__), "..", "media", "convergence.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
