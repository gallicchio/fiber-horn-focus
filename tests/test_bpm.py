"""The long-distance BPM: it must conserve the guided mode and reproduce the launch overlap."""

import numpy as np
import pytest

from fiberlaunch import bpm as B
from fiberlaunch import fibermode as fm
from fiberlaunch import params as P


@pytest.fixture(scope="module")
def small():
    return B.BPM(n=256, half_width=16.0, coating=False)


def test_mode_matches_analytic_lp01(small):
    psi = fm.lp01_field(small.r)
    psi /= np.sqrt(np.sum(psi**2) * small.dA)
    assert np.sum(psi * small.mode) * small.dA > 0.999


def test_guided_mode_is_conserved(small):
    out = small.run(small.mode.astype(complex), 200.0, sample_every=50.0)
    assert out["total"][-1] == pytest.approx(1.0, abs=2e-4)
    assert out["guided"][-1] == pytest.approx(1.0, abs=2e-4)


@pytest.mark.parametrize("name", P.CASE_ORDER)
def test_launch_overlap(small, name):
    w0 = P.CASES[name]["w0"]
    u = small.gaussian(w0)
    expect = fm.gaussian_overlap(w0) * fm.fresnel_T(P.N_CLAD)
    assert abs(small.guided_amp(u)) ** 2 == pytest.approx(expect, abs=5e-3)


def test_guided_power_constant_while_cladding_light_moves(small):
    u = small.gaussian(P.CASES["narrow"]["w0"])
    out = small.run(u, 300.0, sample_every=100.0)
    assert np.ptp(out["guided"]) < 1e-3
