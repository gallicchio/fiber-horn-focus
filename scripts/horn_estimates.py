"""Quick estimates for the metal guide, to set beside the simulation (no simulation data needed).

1. The Gaussian's power outside the mouth: what a beam focused on the mouth sends onto the face.
2. How well the horn's mouth field matches a Gaussian: TE11 of a pipe of the mouth's radius on a
   flat mouth, and on the spherical wavefront of the 8° cone.  Then the mouth fields of better
   horns of the same size, on a flat mouth: a corrugated horn's HE11 (x-polarized J0, zero at the
   wall) and a dual-mode (Potter) horn's TE11 + TM11, mixed in the proportion that matches best.
3. Ray optics: every bounce off a cone's wall keeps a ray's closest distance to the apex (each
   tangent plane of a cone contains the apex), and the ray turns back at that distance.  So a ray
   reaches the throat only if its line passes within the throat's distance of the apex.  Rays are
   drawn from each Gaussian at its waist (positions and directions both Gaussian).

    python scripts/horn_estimates.py | tee media/horn_estimates.txt
"""

import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import j0

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import beam  # noqa: E402
from fiberlaunch.horn_analysis import guide_mode  # noqa: E402
from horn_design import WAISTS, aperture_theory, beam_power  # noqa: E402

TM11_ROOT = 3.8317  # first zero of J1
HE11_ROOT = 2.4048  # first zero of J0


def heading_out(f):
    return [f[0], f[1], f[2], -f[3], -f[4], f[5]]


def best_gaussian(rho, h):
    """(η, w0) of the Gaussian focused on the mouth that couples best to the mouth field h."""
    def eta(w0):
        g = beam.family(rho, 0.0, w0)
        return abs(beam.reaction(rho, g, h) - beam.reaction(rho, h, g)) ** 2 / (
            16 * beam_power(w0) * abs(beam.power(rho, h)))
    w = minimize_scalar(lambda w0: -eta(w0), bounds=(0.5, 2.5), method="bounded").x
    return eta(w), w


def other_horns(a=P.MOUTH_RADIUS):
    """Best Gaussian coupling of the flat mouth fields of a smooth, a corrugated and a dual-mode horn."""
    rho = np.linspace(0, a, 801)
    te = guide_mode(rho, a, "TE", P.TE11_ROOT)
    tm = guide_mode(rho, a, "TM", TM11_ROOT)
    f = j0(HE11_ROOT * rho / a)
    he = [f, -f, 0 * rho, f, f, 0 * rho]  # x-polarized E = J0, H = ẑ × E: a balanced hybrid mode
    mix = lambda c: heading_out([u + c * v for u, v in zip(te, tm)])
    c = minimize_scalar(lambda c: -best_gaussian(rho, mix(c))[0], bounds=(-2, 2), method="bounded").x
    p_tm = c ** 2 * beam.power(rho, tm) / beam.power(rho, [u + c * v for u, v in zip(te, tm)])
    return [("smooth wall, TE11", *best_gaussian(rho, heading_out(te))),
            ("corrugated wall, HE11", *best_gaussian(rho, heading_out(he))),
            (f"dual-mode (Potter), TE11 + TM11 with {100 * p_tm:.0f}% of the power", *best_gaussian(rho, mix(c)))]


def rays(w0, n=2_000_000, seed=1):
    """Fractions of a Gaussian's rays that hit the face, turn back in the horn, reach the throat."""
    a = P.MOUTH_RADIUS
    z_apex = a / np.tan(np.radians(P.FLARE_DEG))
    r_throat = P.THROAT_RADIUS / np.sin(np.radians(P.FLARE_DEG))  # apex to the throat, along the wall
    rng = np.random.default_rng(seed)
    x, y = rng.normal(0, w0 / 2, (2, n))  # intensity exp(-2 r²/w0²)
    sx, sy = rng.normal(0, 1 / (2 * np.pi * w0), (2, n))  # direction sines, half-width λ/(π w0)
    keep = sx ** 2 + sy ** 2 < 1
    x, y, sx, sy = x[keep], y[keep], sx[keep], sy[keep]
    d = np.stack([sx, sy, np.sqrt(1 - sx ** 2 - sy ** 2)], axis=1)
    to_apex = np.stack([-x, -y, np.full_like(x, z_apex)], axis=1)
    b = np.linalg.norm(np.cross(to_apex, d), axis=1)  # closest distance of the ray's line to the apex
    inside = np.hypot(x, y) < a
    return np.mean(~inside), np.mean(inside & (b > r_throat)), np.mean(inside & (b <= r_throat)), r_throat / z_apex


def main():
    waists = {name: c["w0"] for name, c in P.CASES.items()}
    print("Gaussian power outside the mouth (%):",
          "  ".join(f"{n} {100 * np.exp(-2 * P.MOUTH_RADIUS ** 2 / w ** 2):.1f}" for n, w in waists.items()))
    for label, flare in (("flat mouth", 1e-4), (f"spherical wavefront ({P.FLARE_DEG:g}° cone)", P.FLARE_DEG)):
        eta = aperture_theory(flare=flare)
        e = np.array([eta(w) for w in WAISTS])
        j = int(np.argmax(e))
        print(f"TE11 on the mouth, {label}: best Gaussian w0 {WAISTS[j]:.2f} couples {100 * e[j]:.1f}%;  "
              + "  ".join(f"{n} {100 * eta(w):.1f}" for n, w in waists.items()))
    print("other horns with the same mouth, flat (best Gaussian):")
    for name, eta, w in other_horns():
        print(f"  {name}: w0 {w / P.MOUTH_RADIUS:.2f} × the mouth radius ({w:.2f}), {100 * eta:.1f}%")
    print("ray optics (%):                 off the face  turn back in the horn  reach the throat")
    for n, w in waists.items():
        face, back, through, sin_acc = rays(w)
        print(f"  {n:8s}                       {100 * face:5.1f}        {100 * back:5.1f}               {100 * through:5.1f}")
    print(f"a ray through the centre of the mouth reaches the throat if sin(angle) < {sin_acc:.3f} "
          f"({np.degrees(np.arcsin(sin_acc)):.1f}°)")


if __name__ == "__main__":
    main()
