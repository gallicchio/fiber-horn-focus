"""Turn raw FDTD output into the power budget and the "guided" / "rest" field split.

Every steady-state field here is a complex DFT amplitude with the e^{-i omega t}
convention, so the time-averaged power through a plane is Re[integral (E* x H) . z dA]
(Meep's flux convention: the same constant multiplies everything, and only ratios to the
incident power are reported).

The guided part of the field at a plane is found by projecting onto HE11 with the
mode-orthogonality relation

    a_fwd = [ int (E x h*) . z  +  int (e* x H) . z ] / [ int (e x h*) . z  +  int (e* x h) . z ]

where (e, h) is the HE11 mode *as the FDTD itself propagates it* (recorded in the "mode"
run on the same grid), so discretization does not leak into the split.  The rest of the
field, E - a e, is everything that is not in the guided mode: cladding and radiation
modes, which is the light this project is about.
"""

import json
import os

import numpy as np

from . import params as P

TRANSVERSE = ("Ex", "Ey", "Hx", "Hy")


def load_meta(run_dir):
    with open(os.path.join(run_dir, "meta.json")) as f:
        return json.load(f)


def flux_z(f, dA, mask=None):
    """Re int (E* x H) . z dA for a dict (or npz) of transverse components, optionally
    only where ``mask`` (x, y) is true."""
    sz = np.conj(f["Ex"]) * f["Hy"] - np.conj(f["Ey"]) * f["Hx"]
    if mask is not None:
        sz = sz * mask
    return float(np.real(sz.sum(axis=(-2, -1))) * dA) if np.ndim(sz) == 2 else np.real(sz.sum(axis=(-2, -1))) * dA


class Mode:
    """HE11 (x-polarized) on the DFT grid, scaled to unit power."""

    def __init__(self, mode_dir):
        d = np.load(os.path.join(mode_dir, "mode_xy.npz"))
        self.x, self.y = d["x"], d["y"]
        self.dA = (self.x[1] - self.x[0]) * (self.y[1] - self.y[0])
        f = {k: d[k].astype(np.complex128) for k in ("Ex", "Ey", "Ez", "Hx", "Hy")}
        # unit power, and real-positive Ex on axis
        p = flux_z(f, self.dA)
        i0, j0 = np.argmin(abs(self.x)), np.argmin(abs(self.y))
        ph = f["Ex"][i0, j0] / abs(f["Ex"][i0, j0])
        self.f = {k: v / np.sqrt(p) / ph for k, v in f.items()}
        e, h = self.f, self.f
        n1 = np.sum(e["Ex"] * np.conj(h["Hy"]) - e["Ey"] * np.conj(h["Hx"])) * self.dA
        self.norm = n1 + np.conj(n1)  # = 2 Re N

    def project(self, F):
        """Forward HE11 amplitude of the field F (dict of arrays, last two axes x, y).
        |a|^2 is the guided power in flux units."""
        e = self.f
        t1 = np.sum(F["Ex"] * np.conj(e["Hy"]) - F["Ey"] * np.conj(e["Hx"]), axis=(-2, -1))
        t2 = np.sum(np.conj(e["Ex"]) * F["Hy"] - np.conj(e["Ey"]) * F["Hx"], axis=(-2, -1))
        return (t1 + t2) * self.dA / self.norm

    def project_backward(self, F):
        e = self.f
        t1 = np.sum(F["Ex"] * np.conj(e["Hy"]) - F["Ey"] * np.conj(e["Hx"]), axis=(-2, -1))
        t2 = np.sum(np.conj(e["Ex"]) * F["Hy"] - np.conj(e["Ey"]) * F["Hx"], axis=(-2, -1))
        return (t1 - t2) * self.dA / self.norm


class Run:
    """One launch: the fiber run, its air (incident) run, and the reference mode."""

    def __init__(self, fiber_dir, air_dir, mode_dir):
        self.dir = fiber_dir
        self.meta = load_meta(fiber_dir)
        self.air = load_meta(air_dir)
        self.P_inc = self.air["P_net_refl_plane"]
        self.mode = Mode(mode_dir)
        self.incident_facet = np.load(os.path.join(air_dir, "incident_facet.npz"))
        self.name = os.path.basename(os.path.normpath(fiber_dir))

    # ---- scalar budget --------------------------------------------------------------
    @property
    def R(self):
        """Facet reflection.  In the homogeneous air between source and facet the
        incident/reflected cross terms of the time-averaged flux vanish, so the net
        flux there is exactly P_inc - P_refl."""
        return 1 - self.meta["P_net_refl_plane"] / self.P_inc

    def mpb_table(self):
        """Per mode plane: z, guided (MPB), total forward flux in the window, all / P_inc."""
        rows = self.meta["mode_planes"]
        return dict(z=np.array([r["z"] for r in rows]),
                    guided=np.array([r["P_fwd"] for r in rows]) / self.P_inc,
                    backward=np.array([r["P_bwd"] for r in rows]) / self.P_inc,
                    total=np.array([r["flux"] for r in rows]) / self.P_inc)

    # ---- the stack of transverse planes -----------------------------------------------
    def has_planes(self):
        return os.path.exists(os.path.join(self.dir, "xy_coords.npz"))

    def planes(self):
        c = np.load(os.path.join(self.dir, "xy_coords.npz"))
        mm = {n: np.load(os.path.join(self.dir, f"xy_{n}.npy"), mmap_mode="r")
              for n in ("Ex", "Ey", "Ez", "Hx", "Hy")}
        return c["x"], c["y"], c["z"], mm

    def plane_budget(self, chunk=16):
        """Guided amplitude and power, and total flux, at every transverse plane."""
        if getattr(self, "_pb", None) is not None:
            return self._pb
        x, y, z, mm = self.planes()
        dA = (x[1] - x[0]) * (y[1] - y[0])
        hw = self.meta["half_width"] + 1e-9  # the samples beyond are inside the side PML
        inside = (abs(x)[:, None] <= hw) & (abs(y)[None, :] <= hw)
        nz = len(z)
        a = np.zeros(nz, complex)
        tot = np.zeros(nz)
        for s in range(0, nz, chunk):
            F = {k: np.asarray(mm[k][s:s + chunk], dtype=np.complex128) for k in TRANSVERSE}
            a[s:s + chunk] = self.mode.project(F)
            tot[s:s + chunk] = flux_z(F, dA, inside)
        self._pb = dict(z=z, a=a, guided=abs(a) ** 2 / self.P_inc, total=tot / self.P_inc)
        return self._pb

    def plane_fields(self, i, a=None):
        """Total and rest |E|^2 at plane i (rest = minus the guided part)."""
        x, y, z, mm = self.planes()
        F = {k: np.asarray(mm[k][i], dtype=np.complex128) for k in ("Ex", "Ey", "Ez", "Hx", "Hy")}
        if a is None:
            a = self.mode.project(F)
        e = self.mode.f
        tot = sum(abs(F[k]) ** 2 for k in ("Ex", "Ey", "Ez"))
        rest = sum(abs(F[k] - a * e[k]) ** 2 for k in ("Ex", "Ey", "Ez"))
        return tot, rest

    def guided_fit(self, zmin=3.0):
        """In a straight fiber the guided amplitude is exactly A exp(i beta z): fit it to the
        projections of the plane stack.  Returns A, beta and the rms relative scatter of |a|
        about |A| (a check that the guided power really is constant along the fiber)."""
        if getattr(self, "_fit", None) is None:
            b = self.plane_budget()
            z, a = b["z"], b["a"]
            sel = z > zmin
            neff = self.meta["mode_planes"][0]["neff"] if self.meta.get("mode_planes") else 1.4611
            beta0 = 2 * np.pi * neff / P.WAVELENGTH
            ph = np.unwrap(np.angle(a[sel] * np.exp(-1j * beta0 * z[sel])))
            slope = np.polyfit(z[sel], ph, 1)[0]
            beta = beta0 + slope
            A = np.mean(a[sel] * np.exp(-1j * beta * z[sel]))
            scatter = np.std(abs(a[sel]) / abs(A) - 1)
            self._fit = (A, beta, scatter)
        return self._fit

    def side(self, which="xz"):
        return np.load(os.path.join(self.dir, f"dft_side_{which}.npz"))

    def side_split(self, which="xz"):
        """|E|^2 on the side plane through the axis: total, guided part and the rest.
        The guided part is A exp(i beta z) times the mode's profile on that plane."""
        d = self.side(which)
        x = d["x"] if which == "xz" else d["y"]
        z = d["z"]
        A, beta, _ = self.guided_fit()
        m = self.mode
        if which == "xz":  # the plane y = 0 is midway between two voxel-centre rows
            j = np.argsort(abs(m.y))[:2]
            prof = {k: m.f[k][:, j].mean(axis=1) for k in ("Ex", "Ey", "Ez")}
        else:
            i = np.argsort(abs(m.x))[:2]
            prof = {k: m.f[k][i, :].mean(axis=0) for k in ("Ex", "Ey", "Ez")}
        assert len(prof["Ex"]) == len(x), "side plane and mode plane grids differ"
        az = np.where(z > 0, A * np.exp(1j * beta * z), 0)
        tot = np.zeros((len(x), len(z)))
        gui = np.zeros_like(tot)
        rest = np.zeros_like(tot)
        for k in ("Ex", "Ey", "Ez"):
            E = d[k].astype(np.complex128)
            G = prof[k][:, None] * az[None, :]
            tot += abs(E) ** 2
            gui += abs(G) ** 2
            rest += abs(E - G) ** 2
        return dict(x=x, z=z, total=tot, guided=gui, rest=rest, Ex=d["Ex"])

    def frames(self):
        return (np.load(os.path.join(self.dir, "frames_xz_Ex.npy"), mmap_mode="r"),
                np.load(os.path.join(self.dir, "frames_t.npy")))

    def exact_beam_theory(self, n2=None):
        """exact_beam_coupling for this launch's incident beam, into glass of index n_eff."""
        return exact_beam_coupling(self.incident_facet, self.mode, n2 or self.meta["mode_planes"][0]["neff"])


def load_case(root, case):
    """The Run for one named launch in a campaign directory (data/main/res16)."""
    return Run(os.path.join(root, case), os.path.join(root, f"{case}_air"), os.path.join(root, "mode"))


def exact_beam_coupling(incident, mode, n2):
    """Coupling predicted from the incident beam itself, without the paraxial Gaussian:
    the incident E on the facet (from an air run) is split into plane waves, each is carried
    into glass of index n2 with its own s/p Fresnel coefficient, and the transmitted field
    just inside is projected onto HE11.  Only the incident E is used, with continuum
    plane-wave relations for every H (on the FDTD's centred grid the E/H ratio is off by
    cos(k dx/2), about 2% at 16 points per wavelength).  The incident plane may be wider
    than the mode's; the mode is padded with zeros.
    Returns (guided fraction, transmitted fraction) of the incident power."""
    k0 = 2 * np.pi * P.FREQ
    x = incident["x"]
    dx = x[1] - x[0]
    pad = (len(x) - len(mode.x)) // 2
    assert pad >= 0 and abs(x[pad] - mode.x[0]) < 1e-6 * dx, "incident and mode grids differ"
    m = Mode.__new__(Mode)
    m.f = {c: np.pad(v, pad) for c, v in mode.f.items()}
    m.norm, m.dA = mode.norm, mode.dA
    k = 2 * np.pi * np.fft.fftfreq(len(x), dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    kt = np.hypot(KX, KY)
    prop = kt < k0  # evanescent components carry no power into the glass
    safe = np.where(kt > 0, kt, 1.0)
    ux, uy = np.where(kt > 0, KX / safe, 1.0), np.where(kt > 0, KY / safe, 0.0)
    ex, ey = np.fft.fft2(incident["Ex"].astype(complex)), np.fft.fft2(incident["Ey"].astype(complex))
    es, ep = -uy * ex + ux * ey, ux * ex + uy * ey  # s, and the transverse part of p

    def fields(n, ts, tp):
        """Transverse (E, H) in a medium of index n, from its s and transverse-p parts."""
        kz = np.sqrt((n * k0) ** 2 - kt ** 2 + 0j)
        Ex, Ey = np.where(prop, ux * tp - uy * ts, 0), np.where(prop, uy * tp + ux * ts, 0)
        Ez = -(KX * Ex + KY * Ey) / np.where(prop, kz, 1)
        Hx, Hy = (KY * Ez - kz * Ey) / k0, (kz * Ex - KX * Ez) / k0
        return {c: np.fft.ifft2(v) for c, v in dict(Ex=Ex, Ey=Ey, Hx=Hx, Hy=Hy).items()}

    c1 = np.sqrt(np.maximum(k0 ** 2 - kt ** 2, 0)) / k0  # cos(angle) in air, in glass
    c2 = np.sqrt((n2 * k0) ** 2 - kt ** 2 + 0j) / (n2 * k0)
    c1s = np.where(prop, c1, 1)  # grazing and evanescent waves are masked out anyway
    t_s = np.where(prop, 2 * c1 / (c1s + n2 * c2), 0)
    t_p = np.where(prop, 2 * c2 / (n2 * c1s + c2), 0)  # on the transverse p part: t_p cos(t2) / cos(t1)
    dA = dx * dx
    p_in = flux_z(fields(1.0, es, ep), dA)
    glass = fields(n2, t_s * es, t_p * ep)
    return float(abs(m.project(glass)) ** 2 / p_in), flux_z(glass, dA) / p_in
