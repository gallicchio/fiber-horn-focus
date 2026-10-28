"""The analytic fiber theory the simulations are compared against."""

import numpy as np
import pytest

from fiberlaunch import fibermode as fm
from fiberlaunch import params as P


def test_single_mode():
    assert P.V_NUMBER < 2.405


def test_lp01_against_rudolph_neumann_fit():
    m = fm.lp01()
    b_fit = (1.1428 - 0.9960 / m["V"]) ** 2  # good to ~0.2 % for 1.5 < V < 2.5
    assert m["b"] == pytest.approx(b_fit, rel=5e-3)
    assert P.N_CLAD < m["neff"] < P.N_CORE


def test_vector_and_scalar_modes_agree_for_weak_guidance():
    assert fm.he11_neff() == pytest.approx(fm.lp01()["neff"], abs=2e-4)


def test_matched_waist_is_the_optimum_and_close_to_marcuse():
    w, eta = fm.optimal_waist()
    assert P.W_MATCHED == pytest.approx(w, rel=1e-3)
    assert eta > 0.99
    V = P.V_NUMBER
    marcuse = P.CORE_RADIUS * (0.65 + 1.619 / V**1.5 + 2.879 / V**6)
    assert w == pytest.approx(marcuse, rel=0.01)


def test_overlap_is_a_fraction_and_falls_off_both_ways():
    w = P.W_MATCHED
    etas = [fm.gaussian_overlap(f * w) for f in (1 / 3, 1 / 2, 1, 2, 3)]
    assert all(0 < e <= 1 for e in etas)
    assert etas[2] == max(etas)
    # the Gaussian approximation is symmetric in w0 / wm -> wm / w0; the true mode nearly so
    assert etas[0] == pytest.approx(etas[4], abs=0.04)
    assert fm.gaussian_gaussian(1.0, 1.0) == 1.0


def test_ray_optics_misses_the_too_narrow_launch():
    w0 = P.CASES["narrow"]["w0"]
    assert P.CASES["narrow"]["na"] < P.NA_FIBER  # inside the acceptance cone ...
    assert fm.ray_optics(w0) < 0.5 * fm.gaussian_overlap(w0)  # ... yet ray optics is way off


def test_fresnel():
    assert fm.fresnel_T(1.45) == pytest.approx(1 - (0.45 / 2.45) ** 2)
    assert np.isclose(fm.fresnel_T(1.0), 1.0)
