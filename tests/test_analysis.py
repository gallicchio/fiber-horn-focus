"""The non-paraxial coupling estimate, on synthetic fields (no simulation data needed)."""

import numpy as np
import pytest

from fiberlaunch import fibermode as fm
from fiberlaunch.analysis import Mode, exact_beam_coupling, flux_z


def lp01_mode(x):
    """An x-polarized LP01 'Mode' on the grid x (H = n_eff E in units with Z0 = 1)."""
    X, Y = np.meshgrid(x, x, indexing="ij")
    psi = fm.lp01_field(np.hypot(X, Y)).astype(complex)
    n = fm.lp01()["neff"]
    z = np.zeros_like(psi)
    f = dict(Ex=psi, Ey=z, Ez=z, Hx=z, Hy=n * psi)
    m = Mode.__new__(Mode)
    m.x = m.y = x
    m.dA = (x[1] - x[0]) ** 2
    p = flux_z(f, m.dA)
    m.f = {k: v / np.sqrt(p) for k, v in f.items()}
    m.norm = 2 * flux_z(m.f, m.dA)
    return m, n


def gaussian(x, w0, pad=0):
    xx = np.pad(x, pad, mode="linear_ramp", end_values=(x[0] - pad * (x[1] - x[0]), x[-1] + pad * (x[1] - x[0])))
    X, Y = np.meshgrid(xx, xx, indexing="ij")
    E = np.exp(-(X**2 + Y**2) / w0**2)
    return dict(x=xx, y=xx, Ex=E, Ey=0 * E)


@pytest.fixture(scope="module")
def grid():
    dx = 1 / 16
    return (np.arange(512) - 255.5) * dx  # voxel centres, like the FDTD's planes


@pytest.mark.parametrize("w0", [3.0, 5.0])
def test_paraxial_beams_reduce_to_overlap_times_fresnel(grid, w0):
    mode, n = lp01_mode(grid)
    eta, T = exact_beam_coupling(gaussian(grid, w0), mode, n)
    assert T == pytest.approx(fm.fresnel_T(n), abs=2e-4)
    assert eta == pytest.approx(fm.gaussian_overlap(w0) * fm.fresnel_T(n), rel=2e-3)


def test_no_glass_means_no_reflection(grid):
    mode, _ = lp01_mode(grid)
    _, T = exact_beam_coupling(gaussian(grid, 1.0), mode, 1.0)
    assert T == pytest.approx(1.0, abs=1e-9)


def test_a_wider_incident_plane_than_the_mode_is_padded(grid):
    mode, n = lp01_mode(grid)
    a = exact_beam_coupling(gaussian(grid, 2.0), mode, n)
    b = exact_beam_coupling(gaussian(grid, 2.0, pad=32), mode, n)
    assert b == pytest.approx(a, rel=1e-6)
