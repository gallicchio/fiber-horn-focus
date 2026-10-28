"""The exact Gaussian beam and the metal guide's modes (no simulation data needed)."""

import numpy as np
import pytest

from fiberlaunch import beam
from fiberlaunch.horn_analysis import MODES, cross, guide_mode, project, reverse

RHO = np.linspace(0, 0.36, 2001)  # across a guide of the throat's radius
ROOTS = {name: (kind, root) for name, kind, root in MODES}


def test_beam_is_the_gaussian_at_its_waist():
    """Ex = e_r cos²φ - e_φ sin²φ is exp(-ρ²/w0²) on the waist, for a beam wide enough that
    leaving out the evanescent waves changes nothing."""
    rho = np.linspace(0, 9, 301)
    f = beam.family(rho, 0.0, 3.0)
    assert np.allclose(f[1], -f[0])
    assert np.max(abs(f[0] - np.exp(-rho ** 2 / 9.0))) < 1e-5


def test_beam_power_is_the_paraxial_power():
    rho = np.linspace(0, 14, 1401)
    assert beam.power(rho, beam.family(rho, 0.0, 3.0)) == pytest.approx(np.pi * 9 / 4, rel=1e-3)


@pytest.mark.parametrize("w0", [1.0, 1.66])
def test_beam_power_is_the_same_on_every_plane(w0):
    rho = np.linspace(0, 4 * w0 + 6, 2401)  # (a tighter focus sends some light past this at grazing angles)
    p = [beam.power(rho, beam.family(rho, z, w0)) for z in (0.0, -2.0)]
    assert p[1] == pytest.approx(p[0], rel=2e-3)


def test_a_beam_couples_fully_into_its_time_reverse():
    rho = np.linspace(0, 10, 1001)
    f = beam.family(rho, -1.0, 1.66)
    back = tuple(np.conj(e) for e in f[:3]) + tuple(-np.conj(h) for h in f[3:])
    assert beam.coupling(rho, f, back) == pytest.approx(1.0, abs=1e-9)


def test_meep_round_trip():
    rho = np.linspace(0, 3, 31)
    f = beam.family(rho, 0.5, 1.0)
    for a, b in zip(beam.from_meep(beam.to_meep(f)), f):
        assert np.allclose(a, b)


@pytest.mark.parametrize("name", ["TE11", "TM11", "TE12"])
def test_project_finds_a_pure_mode(name):
    rho = np.linspace(0, 1.0, 2001)  # a guide of radius 1 λ, where all three propagate
    m = guide_mode(rho, 1.0, *ROOTS[name])
    P = 0.5 * np.real(cross(rho, m, m))
    p_in, p_out = project(rho, m, m)
    assert p_in == pytest.approx(P, rel=1e-9) and p_out < 1e-12 * P
    out = m[:3] + tuple(-h for h in m[3:])
    p_in, p_out = project(rho, out, m)
    assert p_out == pytest.approx(P, rel=1e-9) and p_in < 1e-12 * P


@pytest.mark.parametrize("name", ["TE11", "TM12"])
def test_project_separates_the_directions_on_a_curved_wavefront(name):
    rho = np.linspace(0, 1.5, 3001)
    ph = np.exp(-2j * np.pi * (np.hypot(rho, 4.0) - 4.0))  # converging on a point 4 λ ahead
    m = tuple(c * ph for c in guide_mode(rho, 1.5, *ROOTS[name]))
    P = 0.5 * np.real(cross(rho, m, m))
    field = tuple(0.6 * a + 0.8j * b for a, b in zip(m, reverse(m)))
    p_in, p_out = project(rho, field, m)
    assert p_in == pytest.approx(0.36 * P, rel=1e-9) and p_out == pytest.approx(0.64 * P, rel=1e-9)


def test_modes_are_orthogonal():
    rho = np.linspace(0, 1.0, 4001)
    te11, tm11 = (guide_mode(rho, 1.0, *ROOTS[n]) for n in ("TE11", "TM11"))
    p_in, p_out = project(rho, tm11, te11)
    assert p_in + p_out < 1e-6 * 0.5 * np.real(cross(rho, tm11, tm11))


def test_throat_is_single_mode():
    """Only TE11 propagates in the throat (radius 0.36 λ)."""
    alive = [name for name, kind, root in MODES if guide_mode(RHO, 0.36, kind, root) is not None]
    assert alive == ["TE11"]
