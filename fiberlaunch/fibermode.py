"""Analytic step-index fiber modes and the launch-efficiency theory the FDTD is checked against.

* the scalar LP01 mode (weakly guiding), its profile, and its overlap with a Gaussian beam
  focused on the facet (the textbook "mode-matching" coupling efficiency);
* the exact vector HE11 propagation constant, to check MPB/Meep's discretized fiber;
* the ray-optics prediction, which is what "the beam NA fits inside the fiber NA"
  intuition gives.
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import j0, j1, jv, jvp, k0, k1, kv, kvp

from . import params as P


def _vw(u, V):
    return np.sqrt(V**2 - u**2)


def lp01(a=P.CORE_RADIUS, n1=P.N_CORE, n2=P.N_CLAD, lam=P.WAVELENGTH):
    """Scalar LP01 mode.  Returns dict(u, w, V, b, neff)."""
    V = 2 * np.pi * a * np.sqrt(n1**2 - n2**2) / lam

    def f(u):
        w = _vw(u, V)
        return u * j1(u) / j0(u) - w * k1(w) / k0(w)

    u = brentq(f, 1e-9, min(2.404825557695773, V) - 1e-12)
    w = _vw(u, V)
    b = (w / V) ** 2
    neff = np.sqrt(n2**2 + b * (n1**2 - n2**2))
    return dict(u=u, w=w, V=V, b=b, neff=neff)


def lp01_field(r, a=P.CORE_RADIUS, **kw):
    """LP01 transverse field (any units; 1 on axis) at radius r."""
    m = lp01(a=a, **kw)
    r = np.asarray(r, dtype=float)
    out = np.empty_like(r)
    inside = r < a
    out[inside] = j0(m["u"] * r[inside] / a) / j0(m["u"])
    out[~inside] = k0(m["w"] * r[~inside] / a) / k0(m["w"])
    return out


def he11_neff(a=P.CORE_RADIUS, n1=P.N_CORE, n2=P.N_CLAD, lam=P.WAVELENGTH):
    """Exact (vector) effective index of the HE11 mode of a step-index fiber."""
    V = 2 * np.pi * a * np.sqrt(n1**2 - n2**2) / lam
    nu = 1
    r = (n2 / n1) ** 2

    def f(u):
        w = _vw(u, V)
        J = jvp(nu, u) / (u * jv(nu, u))
        K = kvp(nu, w) / (w * kv(nu, w))
        return (J + K) * (J + r * K) - nu**2 * (1 / u**2 + 1 / w**2) * (1 / u**2 + r / w**2)

    # HE11 is the lowest root, below the first zero of J0 (2.405) and of J1 (3.83)
    us = np.linspace(1e-3, min(V, 2.4048) - 1e-6, 4000)
    fs = f(us)
    i = np.flatnonzero(np.sign(fs[:-1]) != np.sign(fs[1:]))[0]
    u = brentq(f, us[i], us[i + 1])
    w = _vw(u, V)
    b = (w / V) ** 2
    return np.sqrt(n2**2 + b * (n1**2 - n2**2))


def _radial_grid(rmax=60.0, n=600001):
    r = np.linspace(0, rmax, n)
    return r, r[1] - r[0]


def gaussian_overlap(w0, a=P.CORE_RADIUS, **kw):
    """Power coupling of a Gaussian (waist w0 on the facet, on axis) into LP01.

    This is the scalar mode-matching overlap only; multiply by the facet's Fresnel
    transmission to compare with a full-wave calculation."""
    r, dr = _radial_grid(rmax=max(60.0, 12 * w0))
    psi = lp01_field(r, a=a, **kw)
    g = np.exp(-(r**2) / w0**2)
    num = (np.sum(psi * g * r) * dr) ** 2
    den = np.sum(psi**2 * r) * dr * np.sum(g**2 * r) * dr
    return num / den


def optimal_waist(**kw):
    res = minimize_scalar(lambda w: -gaussian_overlap(w, **kw), bounds=(0.3, 10), method="bounded")
    return res.x, -res.fun


def gaussian_gaussian(w0, wm):
    """Overlap of two Gaussians with waists w0 and wm (the usual approximation)."""
    return (2 * w0 * wm / (w0**2 + wm**2)) ** 2


def ray_optics(w0, a=P.CORE_RADIUS, na=P.NA_FIBER, lam=P.WAVELENGTH):
    """What ray optics predicts: the fraction of rays that hit the core *and* lie inside
    the acceptance cone.  A Gaussian beam's rays (its Wigner function) are Gaussian in
    both position (1/e^2 radius w0) and angle (1/e^2 half-angle lam/(pi w0)), and a ray in
    a step-index fiber is bound if its angle to the axis is inside the acceptance cone,
    whatever its skewness."""
    theta0 = lam / (np.pi * w0)
    return (1 - np.exp(-2 * a**2 / w0**2)) * (1 - np.exp(-2 * na**2 / theta0**2))


def fresnel_T(n):
    return 1 - ((n - 1) / (n + 1)) ** 2


if __name__ == "__main__":
    m = lp01()
    print(f"n_core = {P.N_CORE:.5f}, n_clad = {P.N_CLAD}, NA = {P.NA_FIBER}, a = {P.CORE_RADIUS}")
    print(f"V = {m['V']:.4f}, LP01 neff = {m['neff']:.6f}, HE11 neff = {he11_neff():.6f}")
    w, eta = optimal_waist()
    print(f"best Gaussian waist {w:.4f} (beam NA {P.beam_na(w):.4f}), overlap {eta:.4f}")
    print(f"beat length lambda/(neff - n_clad) = {1 / (m['neff'] - P.N_CLAD):.1f}")
    for name in P.CASE_ORDER:
        c = P.CASES[name]
        print(
            f"{name:8s} w0 {c['w0']:.3f} NA {c['na']:.4f}  overlap {gaussian_overlap(c['w0']):.4f}"
            f"  x Fresnel {gaussian_overlap(c['w0']) * fresnel_T(m['neff']):.4f}"
            f"  ray optics {ray_optics(c['w0']) * fresnel_T(m['neff']):.4f}"
        )
