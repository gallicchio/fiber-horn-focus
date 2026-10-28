"""Long-distance continuation: a scalar beam-propagation method (BPM) with the whole fiber
cross-section (core, cladding, coating).

FDTD can afford ~120 wavelengths of fiber.  What happens to the light that did not
couple plays out over millimetres, so this follows it with a split-step Fourier BPM:

* free propagation in the cladding index uses the exact angular-spectrum phase
  exp(i (sqrt(k^2 - kt^2) - k) dz) (no paraxial approximation for the cladding light);
* the core, coating and absorber are phase/amplitude screens exp(i k0 (n^2 - n_ref^2) dz / 2 n_ref);
* the coating has a higher index than the glass, as real acrylate coatings do, so
  light reaching the cladding surface is not totally reflected; the coating and an
  absorbing jacket beyond it remove what crosses into them.

The glass dimensions are those of SMF-28 (125 um cladding, ~190 um primary coating)
scaled by the same factor as the core (see params), so the long-distance behaviour maps
onto SMF-28 by stretching lengths by (0.25/0.14)^2 ~ 3.2.

The scalar field u stands for the dominant transverse component (Ex) of the FDTD field,
with the e^{i n_ref k0 z} carrier removed.  It is checked against the 3D FDTD in
``scripts/run_bpm.py``.
"""

import numpy as np
from scipy import fft as sfft

from . import params as P

# SMF-28: core radius 4.1 um, NA 0.14, 125 um cladding, primary coating to 190 um, at 1.55 um.
SCALE = (4.1 / 1.55) / P.CORE_RADIUS  # SMF-28 core in wavelengths / our core
R_CLAD = 62.5 / 1.55 / SCALE  # 22.9 wavelengths
R_COAT = 95.0 / 1.55 / SCALE  # 34.8 wavelengths
LENGTH_TO_SMF28 = SCALE**2  # longitudinal scaling for the same diffraction physics
# The coating: n 1.48 against silica's 1.444 at 1.55 um, with weak absorption.  Shrinking
# transverse lengths by SCALE and stretching z by SCALE^2 keeps the (paraxial) physics
# only if every index step n^2 - n_clad^2, real and imaginary, grows by SCALE^2, as the
# core's does (NA 0.14 -> 0.25).  Then reflection at the cladding surface per bounce and
# the absorption per SMF-28 length are both those of the real fiber.
N_CLAD_SMF, N_COAT_SMF, K_COAT_SMF = 1.444, 1.48, 1e-3
_dn2 = SCALE**2 * (N_COAT_SMF**2 - N_CLAD_SMF**2 + 2j * N_COAT_SMF * K_COAT_SMF)
_n_coat = np.sqrt(P.N_CLAD**2 + _dn2)
N_COAT, K_COAT = float(_n_coat.real), float(_n_coat.imag)  # 1.559, 0.0030


class BPM:
    def __init__(self, n=768, half_width=None, dz=0.5, coating=True, workers=-1):
        hw = half_width or (R_COAT + 6.0)
        self.n, self.dx = n, 2 * hw / n
        self.x = (np.arange(n) - n // 2) * self.dx
        X, Y = np.meshgrid(self.x, self.x, indexing="ij")
        self.r = np.hypot(X, Y)
        self.k0 = 2 * np.pi / P.WAVELENGTH
        self.nref = P.N_CLAD
        self.dz = dz
        self.workers = workers
        self.coating = coating
        idx = np.full((n, n), P.N_CLAD, dtype=complex)
        idx[self.r < P.CORE_RADIUS] = P.N_CORE
        if coating:
            idx[self.r >= R_CLAD] = N_COAT + 1j * K_COAT
            # jacket: smoothly ramped absorber out to the window edge
            ramp = np.clip((self.r - R_COAT) / (hw - R_COAT), 0, 1)
            idx = idx + 1j * 0.05 * ramp**2 * (self.r >= R_COAT)
        else:
            ramp = np.clip((self.r - (hw - 6)) / 6, 0, 1)
            idx = idx + 1j * 0.05 * ramp**2
        self.index = idx
        kx = 2 * np.pi * sfft.fftfreq(n, self.dx)
        KX, KY = np.meshgrid(kx, kx, indexing="ij")
        # Exact phase for every propagating plane wave.  The sharp core boundary keeps
        # regenerating a little evanescent spectrum (kt > n k0) from the guided mode; damping
        # it, as the true propagator would, drains the mode by ~1 % per 100 wavelengths, so
        # those components just keep the phase they have at kt = n k0 (unitary step).
        kz = np.sqrt(np.maximum((self.k0 * self.nref) ** 2 - KX**2 - KY**2, 0.0))
        self._prop = lambda d: np.exp(1j * (kz - self.k0 * self.nref) * d)
        self.H = self._prop(dz)
        self.Hhalf = self._prop(dz / 2)
        self.screen = np.exp(1j * self.k0 * (idx**2 - self.nref**2) / (2 * self.nref) * dz)
        self.mode = self._fundamental_mode()
        self.dA = self.dx**2

    def _fft(self, u):
        return sfft.fft2(u, workers=self.workers)

    def _ifft(self, u):
        return sfft.ifft2(u, workers=self.workers)

    def _fundamental_mode(self, iters=600):
        """The BPM's own LP01 (core in an unbounded cladding), by imaginary-distance
        propagation, so the guided/unguided split has no discretization leak."""
        idx = np.where(self.r < P.CORE_RADIUS, P.N_CORE, P.N_CLAD)
        d = 0.5
        kx = 2 * np.pi * sfft.fftfreq(self.n, self.dx)
        KX, KY = np.meshgrid(kx, kx, indexing="ij")
        kz = np.sqrt((self.k0 * self.nref) ** 2 - KX**2 - KY**2 + 0j)
        Hd = np.exp((np.real(kz) - self.k0 * self.nref) * d / 2) * (np.imag(kz) == 0)
        sc = np.exp(self.k0 * (idx**2 - self.nref**2) / (2 * self.nref) * d)
        u = np.exp(-(self.r**2) / 2.0**2).astype(complex)
        for _ in range(iters):
            u = self._ifft(Hd * self._fft(u))
            u = u * sc
            u = self._ifft(Hd * self._fft(u))
            u /= np.sqrt(np.sum(abs(u) ** 2))
        u = np.real(u * np.conj(u[self.n // 2, self.n // 2]) / abs(u[self.n // 2, self.n // 2]))
        return u / np.sqrt(np.sum(u**2) * self.dx**2)

    # ---- inputs ---------------------------------------------------------------------
    def gaussian(self, w0):
        """A Gaussian focused on the facet, after a normal-incidence Fresnel loss."""
        t2 = 1 - ((P.N_CLAD - 1) / (P.N_CLAD + 1)) ** 2
        u = np.exp(-(self.r**2) / w0**2).astype(complex)
        return u * np.sqrt(t2 / (np.sum(abs(u) ** 2) * self.dA))

    def from_fdtd(self, x, y, Ex, z0, power, guided=None):
        """Interpolate an FDTD Ex plane (complex) onto the BPM grid, remove the carrier,
        and scale so that sum |u|^2 dA is that plane's power (as a fraction of the
        incident power).  With ``guided`` (the FDTD's own guided fraction there), the
        guided and the rest are scaled separately, so the BPM starts from exactly the
        FDTD's split: |E|^2 is not proportional to flux for both alike."""
        from scipy.interpolate import RegularGridInterpolator

        it_r = RegularGridInterpolator((x, y), Ex.real, bounds_error=False, fill_value=0)
        it_i = RegularGridInterpolator((x, y), Ex.imag, bounds_error=False, fill_value=0)
        X, Y = np.meshgrid(self.x, self.x, indexing="ij")
        pts = np.stack([X.ravel(), Y.ravel()], -1)
        u = (it_r(pts) + 1j * it_i(pts)).reshape(X.shape)
        u *= np.exp(-1j * self.k0 * self.nref * z0)
        u *= np.sqrt(power / (np.sum(abs(u) ** 2) * self.dA))
        if guided is None:
            return u
        a, _, rest = self.split(u)
        return (np.sqrt(guided) * a / abs(a) * self.mode
                + rest * np.sqrt((power - guided) / (np.sum(abs(rest) ** 2) * self.dA)))

    # ---- diagnostics ----------------------------------------------------------------
    def guided_amp(self, u):
        return np.sum(u * self.mode) * self.dA

    def budget(self, u):
        p = abs(u) ** 2 * self.dA
        a = self.guided_amp(u)
        return dict(total=p.sum(), guided=abs(a) ** 2, in_glass=p[self.r < R_CLAD].sum(),
                    near_core=p[self.r < 3 * P.CORE_RADIUS].sum())

    # ---- propagation ----------------------------------------------------------------
    def split(self, u):
        """Guided part a * mode and the rest u - a * mode."""
        a = self.guided_amp(u)
        return a, a * self.mode, u - a * self.mode

    def run(self, u, length, z0=0.0, sample_every=10.0, snap_every=None, snap_stride=2, callback=None):
        """Propagate u from z0 to z0 + length (Strang splitting).

        Every sample_every: the power budget and the x-z cut through the axis, split into
        guided and rest.  Every snap_every: the transverse guided and rest intensities,
        downsampled by snap_stride."""
        nsteps = int(round(length / self.dz))
        every = max(1, int(round(sample_every / self.dz)))
        snap = int(round(snap_every / self.dz)) if snap_every else None
        c = self.n // 2
        rec = dict(z=[], total=[], guided=[], in_glass=[], near_core=[], side_guided=[], side_rest=[])
        snaps = dict(z=[], guided=[], rest=[])

        def record(v, z):
            b = self.budget(v)
            for k in ("total", "guided", "in_glass", "near_core"):
                rec[k].append(b[k])
            rec["z"].append(z)
            _, g, r = self.split(v)
            rec["side_guided"].append((abs(g[:, c]) ** 2).astype(np.float32))
            rec["side_rest"].append((abs(r[:, c]) ** 2).astype(np.float32))

        def take_snap(v, z):
            _, g, r = self.split(v)
            k = snap_stride
            snaps["z"].append(z)
            snaps["guided"].append((abs(g) ** 2).reshape(self.n // k, k, self.n // k, k).mean((1, 3)).astype(np.float32))
            snaps["rest"].append((abs(r) ** 2).reshape(self.n // k, k, self.n // k, k).mean((1, 3)).astype(np.float32))

        record(u, z0)
        if snap:
            take_snap(u, z0)
        U = self._fft(u) * self.Hhalf  # leading half step
        for i in range(1, nsteps + 1):
            u = self._ifft(U) * self.screen
            U = self._fft(u)
            want_rec = i % every == 0 or i == nsteps
            want_snap = bool(snap) and i % snap == 0
            if want_rec or want_snap:
                v = self._ifft(U * self.Hhalf)  # the field at z exactly
                z = z0 + i * self.dz
                if want_rec:
                    record(v, z)
                if want_snap:
                    take_snap(v, z)
                if callback:
                    callback(i, nsteps, v, z)
            U = U * self.H
        out = {k: np.array(v) for k, v in rec.items()}
        out["final"] = self._ifft(U * self._prop(-self.dz / 2))
        if snap:
            out.update({f"snap_{k}": np.array(v) for k, v in snaps.items()})
            out["snap_x"] = self.x.reshape(-1, snap_stride).mean(1)
        out["x"] = self.x
        return out
