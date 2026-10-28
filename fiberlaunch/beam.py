"""The exact (non-paraxial) x-polarized Gaussian beam in cylindrical components.

The beam is the sum of plane waves whose x-component at the waist (z = 0) is the Gaussian
Ex = exp(-ρ²/w0²); evanescent waves are left out, so it is a solution of Maxwell's equations
in air for every w0.  Linear polarization along x on a beam with a round angular spectrum
has only cos φ and sin φ dependence:

    E_r = e_r(ρ) cos φ    E_φ = e_φ(ρ) sin φ    E_z = e_z(ρ) cos φ
    H_r = h_r(ρ) sin φ    H_φ = h_φ(ρ) cos φ    H_z = h_z(ρ) sin φ

so it is exactly half angular order m = +1 and half m = -1, and any body of revolution on
the same axis keeps it that way.  These six 'family' profiles are what this module returns
(units with Z0 = 1 and the wavelength 1, Meep's units).
"""

import numpy as np
from scipy.special import j0, j1

K = 2 * np.pi


def family(rho, z, w0, n_theta=4000):
    """(e_r, e_phi, e_z, h_r, h_phi, h_z) at radii rho on the plane z (waist at z = 0).

    The angular-spectrum integrals are done over the plane-wave angle θ (κ = k sin θ), which
    removes the 1/k_z of the longitudinal components at grazing angles."""
    rho = np.asarray(rho, float)[:, None]
    th = (np.arange(n_theta) + 0.5) * (0.5 * np.pi / n_theta)
    dth = 0.5 * np.pi / n_theta
    s, c = np.sin(th), np.cos(th)
    kap, kz = K * s, K * c
    A = np.pi * w0 ** 2 * np.exp(-(kap * w0) ** 2 / 4) * np.exp(1j * kz * z)
    # κ dκ = k² s c dθ;  κ dκ / k_z = k s dθ
    w_e = A * K * K * s * c * dth / (2 * np.pi)  # weight of ∫ κ dκ (...) / 2π
    w_z = A * K * s * dth / (2 * np.pi)  # weight of ∫ κ dκ (...) / (2π k_z)
    x = kap * rho
    J0, J1 = j0(x), j1(x)
    J2 = np.where(x > 1e-3, 2 * J1 / np.maximum(x, 1e-3) - J0, x * x / 8)
    i0 = (w_e * J0).sum(1)
    e_z = -1j * (w_z * kap * J1).sum(1)
    b = (w_z * (kz ** 2 + 0.5 * kap ** 2) * J0).sum(1) / K
    a = (w_z * 0.5 * kap ** 2 * J2).sum(1) / K
    h_z = -1j * (w_e * kap * J1).sum(1) / K
    return i0, -i0, e_z, b + a, b - a, h_z


def power(rho, f):
    """Time-averaged power through a z-plane of the x-polarized beam described by the family
    profiles f on the radii rho: (1/2) Re ∫(E × H*)·ẑ dA = (π/2) Re ∫(e_r h_φ* - e_φ h_r*) ρ dρ."""
    e_r, e_p, _, h_r, h_p, _ = f
    return 0.5 * np.pi * np.real(np.trapezoid((e_r * np.conj(h_p) - e_p * np.conj(h_r)) * rho, rho))


def reaction(rho, f, g):
    """∫(E_f × H_g)·ẑ dA for two x-polarized fields given by their family profiles."""
    return np.pi * np.trapezoid((f[0] * g[4] - f[1] * g[3]) * rho, rho)


def coupling(rho, f, g):
    """Power coupling between two fields travelling in opposite directions (reciprocity):
    |∫(E_f × H_g - E_g × H_f)·ẑ dA|² / (16 P_f P_g)."""
    num = abs(reaction(rho, f, g) - reaction(rho, g, f)) ** 2
    return num / (16 * abs(power(rho, f)) * abs(power(rho, g)))


def to_meep(f):
    """The m = +1 half of the field, as Meep's cylindrical (m = 1) components
    (Er, Ep, Ez, Hr, Hp, Hz): a cos φ profile c contributes c/2, a sin φ profile s, s/2i."""
    e_r, e_p, e_z, h_r, h_p, h_z = f
    return e_r / 2, e_p / 2j, e_z / 2, h_r / 2j, h_p / 2, h_z / 2j


def from_meep(F):
    """Inverse of to_meep: family profiles of the linearly polarized field whose m = +1 half
    is the Meep field F = (Er, Ep, Ez, Hr, Hp, Hz)."""
    Er, Ep, Ez, Hr, Hp, Hz = F
    return 2 * Er, 2j * Ep, 2 * Ez, 2j * Hr, 2 * Hp, 2j * Hz
