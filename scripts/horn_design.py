"""Size the horn: the coupling of every Gaussian beam into its throat, from one 'emit' run.

By reciprocity, an incoming beam g couples into the guide's TE11 with efficiency
    η = |∫(E_g × H_h - E_h × H_g)·ẑ dA|² / (16 P_g P_h),
where h is the field the horn radiates when TE11 is sent out of the throat with power P_h,
both on one plane just in front of the mouth (z = -0.5).  The beams are the x-polarized
Gaussians focused on the mouth (z = 0), built exactly (fiberlaunch/beam.py).  For the cylindrical
runs it also finds the best beam focused deeper, inside the horn.

    python scripts/horn_design.py data/horn/design/emit_res60 [...]

Runs of fiberlaunch/horn_cartesian.py (a front_xy.npz plane instead of dft_rz.npz) are read too;
those are the 3D Cartesian runs whose staircased walls do not converge.
"""

import argparse
import functools
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import j0

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import beam  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch.horn_analysis import COMPS, guide_mode, project  # noqa: E402

WAISTS = np.round(np.concatenate([np.arange(0.5, 1.4, 0.1), np.arange(1.4, 2.21, 0.05), np.arange(2.4, 3.01, 0.2)]), 3)
Z_FRONT = -0.5


@functools.lru_cache
def beam_power(w0):
    """Power of the whole beam (it is the same on every plane; the waist is the compact one)."""
    rho = np.linspace(0, 4 * w0 + 2, 1201)
    return beam.power(rho, beam.family(rho, 0.0, w0))


def aperture_theory(wall="smooth", mouth=P.MOUTH_RADIUS, flare=P.FLARE_DEG):
    """The textbook estimate: the mouth is filled with the mode of a guide of the mouth's radius,
    on the spherical wavefront of a cone of this flare (centred on its apex), and nothing else.
    The smooth wall's mode is TE11; the corrugated wall's is the balanced hybrid HE11, linearly
    polarized with E ∝ J0(2.405 ρ/a) and H = ẑ × E.  Returns η(w0, depth), the beam's waist at
    z = depth (inside the horn for depth > 0)."""
    rho = np.linspace(0, mouth, 801)
    if wall == "corrugated":
        u = j0(2.4048 * rho / mouth)
        mode = [u, -u, 0 * rho, u, u, 0 * rho]
    else:
        mode = guide_mode(rho, mouth, "TE", P.TE11_ROOT)
    apex = mouth / np.tan(np.radians(flare))  # apex to mouth, along the axis
    ph = np.exp(1j * beam.K * (np.hypot(rho, apex) - apex))  # the wave leaving the horn
    h = [ph * mode[0], ph * mode[1], 0 * rho, -ph * mode[3], -ph * mode[4], 0 * rho]  # heading out (-z)

    def eta(w0, depth=0.0):
        g = beam.family(rho, -depth, w0)
        return abs(beam.reaction(rho, g, h) - beam.reaction(rho, h, g)) ** 2 / (
            16 * beam_power(w0) * abs(beam.power(rho, h)))

    return eta


def cylindrical(run, m):
    d = np.load(os.path.join(run, "dft_rz.npz"))
    r, z = d["r"], d["z"]
    F = beam.from_meep([d[c].astype(np.complex128) for c in COMPS])
    h = [f[np.argmin(abs(z - Z_FRONT))] for f in F]
    zf = z[np.argmin(abs(z - Z_FRONT))]
    # P_h: the TE11 power heading out of the throat, between the source and the horn
    kt = np.argmin(abs(z - (m["horn_length"] + 0.3 * m["throat_length"])))
    te11 = guide_mode(r, m["throat_radius"], "TE", P.TE11_ROOT)
    _, P_h = project(r, [f[kt] for f in F], te11)

    def eta(w0, depth=0.0, spectrum=None):  # waist at z = depth, inside the horn for depth > 0
        g = beam.family(r, zf - depth, w0, spectrum=spectrum)
        if spectrum is None:
            P_g = beam_power(w0)
        else:  # another beam with this waist (beam.launched_spectrum): its power on a wide plane
            rho = np.linspace(0, 60.0, 2401)
            P_g = beam.power(rho, beam.family(rho, 0.0, w0, spectrum=spectrum))
        return abs(beam.reaction(r, g, h) - beam.reaction(r, h, g)) ** 2 / (16 * P_g * P_h)

    return eta, f"plane carries {abs(beam.power(r, h)) / P_h:.4f} of the TE11 power sent out"


def cartesian(run, m):
    h = dict(np.load(os.path.join(run, "front_xy.npz")))
    x, y = h["x"], h["y"]
    dA = (x[1] - x[0]) * (y[1] - y[0])
    z = float(np.atleast_1d(h["z"])[0])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho, phi = np.hypot(X, Y), np.arctan2(Y, X)

    def cart(f):  # family profiles -> Cartesian transverse fields on the grid
        c, s = np.cos(phi), np.sin(phi)
        fr = [np.interp(rho, rr, v.real) + 1j * np.interp(rho, rr, v.imag) for v in f]
        Er, Ep, Hr, Hp = fr[0] * c, fr[1] * s, fr[3] * s, fr[4] * c
        return dict(Ex=Er * c - Ep * s, Ey=Er * s + Ep * c, Hx=Hr * c - Hp * s, Hy=Hr * s + Hp * c)

    def cz(a, b):
        return np.sum(a["Ex"] * b["Hy"] - a["Ey"] * b["Hx"]) * dA

    P_h = -0.5 * np.real(np.sum(h["Ex"] * np.conj(h["Hy"]) - h["Ey"] * np.conj(h["Hx"]))) * dA
    rr = np.linspace(0, np.max(rho) + 0.1, 1001)

    def eta(w0):
        g = cart(beam.family(rr, z, w0))
        P_g = beam_power(w0)
        return abs(cz(g, h) - cz(h, g)) ** 2 / (16 * P_g * P_h)

    return eta, "3D Cartesian run (P_h from the plane itself)"


def extrapolate(res, eta):
    """Richardson extrapolation from the three finest grids: η(res) = η∞ - c / res^p."""
    (h1, h2, h3), (e1, e2, e3) = 1 / np.asarray(res[-3:], float), eta[-3:]
    ratio = (e2 - e1) / (e3 - e2)

    def f(p):
        return (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - ratio

    lo, hi = 0.25, 6.0
    if f(lo) * f(hi) > 0:
        return None, None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
    p = 0.5 * (lo + hi)
    return e3 + (e3 - e2) * h3 ** p / (h2 ** p - h3 ** p), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    series = []  # (res, {name: η}) of the cylindrical runs
    walls = {json.load(open(os.path.join(run, "meta.json"))).get("wall", "smooth") for run in a.runs}
    assert len(walls) == 1, f"one kind of wall at a time, not {walls}"
    wall = walls.pop()
    theory = aperture_theory(wall)
    e = np.array([theory(w) for w in WAISTS])
    j = int(np.argmax(e))
    best = minimize(lambda v: -theory(*v), [WAISTS[j], 2.0], method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-6))
    print(f"aperture theory ({'HE11' if wall == 'corrugated' else 'TE11'} on the mouth, spherical phase): "
          f"best w0 {WAISTS[j]:.2f}: {e[j]:.4f};  "
          + "  ".join(f"{k} {theory(c['w0']):.4f}" for k, c in P.CASES.items())
          + f";  best focused inside: w0 {best.x[0]:.2f}, {best.x[1]:.1f} λ in: {-best.fun:.4f}")
    for run in a.runs:
        m = json.load(open(os.path.join(run, "meta.json")))
        cyl = os.path.exists(os.path.join(run, "dft_rz.npz"))
        eta, note = (cylindrical if cyl else cartesian)(run, m)
        e = np.array([eta(w) for w in WAISTS])
        j = int(np.argmax(e))
        print(f"{os.path.basename(run)}: res {m['res']}, mouth {m['mouth_radius']} flare {m['flare_deg']} "
              f"L {m['horn_length']:.2f}  ({note})")
        print("   " + "  ".join(f"{w:.2f}:{x:.4f}" for w, x in zip(WAISTS, e)))
        print(f"   best w0 {WAISTS[j]:.2f} (NA {P.beam_na(WAISTS[j]):.3f}): {e[j]:.4f};  the fiber's three beams:  "
              + "  ".join(f"{k} (w0 {c['w0']:.3f}) {eta(c['w0']):.4f}" for k, c in P.CASES.items()))
        if cyl:
            best = minimize(lambda v: -eta(*v), [1.45, 3.5], method="Nelder-Mead",
                            options=dict(xatol=1e-3, fatol=1e-6))
            print(f"   focused inside the horn instead: best w0 {best.x[0]:.2f}, "
                  f"{best.x[1]:.1f} λ in from the mouth: {-best.fun:.4f}")
            series.append((m["res"], {**{k: eta(c["w0"]) for k, c in P.CASES.items()}, "best focus inside": -best.fun}))
    series.sort(key=lambda s: s[0])
    if len(series) >= 3:
        res = [s[0] for s in series]
        print(f"extrapolated to an infinitely fine grid from res {res[-3:]}:")
        for k in [*P.CASES, "best focus inside"]:
            inf, p = extrapolate(res, np.array([s[1][k] for s in series]))
            e = np.array([s[1][k] for s in series])
            note = (f"->  {inf:.4f}  (error ~ 1/res^{p:.1f})" if inf is not None and 1 <= p <= 4 else
                    f"(not yet smooth enough to extrapolate; the three finest differ by {np.ptp(e[-3:]):.4f})")
            print(f"   {k:17s} " + "  ".join(f"{r}:{v:.4f}" for r, v in zip(res, e)) + "   " + note)


if __name__ == "__main__":
    main()
