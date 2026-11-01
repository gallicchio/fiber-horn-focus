"""Post-processing for the metal waveguide (horn) runs: the power budget, and the split of the
field into the guide's own mode (TE11, heading in) and the rest.

The runs store the m = +1 half of the field on the (r, z) plane; beam.from_meep() turns it
into the cos φ / sin φ profiles of the linearly polarized field, and every integral over a
surface of revolution reduces to π ∫ (...) w ds (w ds dφ the area element).

The modes of a conical horn are spherical waves about the cone's apex: on a spherical cap
centred on the apex, a mode heading in or out has the same transverse profile and a uniform
phase, and for a narrow cone (8° here) that profile is the TE1n / TM1n profile of a straight
guide whose radius is the cap's arc radius.  So inside the horn the field is taken on such
caps, one for every depth at which a cap meets the wall; in the straight throat, on planes.
On each, every mode's amplitude heading in and heading out comes from Lorentz reciprocity
against the same mode heading the other way,
    a_in = R(F, m_out) / R(m_in, m_out),   R(f, g) = ∫(E_f × H_g - E_g × H_f)·n̂ dA,
which separates the directions and the modes exactly where the profiles are exact (the throat)
and to about a percent in the horn (checked on TE11 sent out from the throat).  In front of
the block the split is incident (the field of the same beam with no metal) and reflected
(total - incident).

A corrugated horn's modes are not those of a smooth guide, so there the part of the field that
ends up in the guide is found another way, which works for any horn: it is the time reverse of
the field the horn radiates when TE11 is sent out of the throat (an 'emit' run on the same grid),
scaled to the TE11 heading in in the throat.  The reflected light is split by where the incident
beam lands: on the metal face around the mouth, or inside the mouth.
"""

import json
import os

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.special import jn_zeros, jnp_zeros, jv, jvp

from . import beam
from . import params as H

COMPS = ("Er", "Ep", "Ez", "Hr", "Hp", "Hz")
K = 2 * np.pi


def radius_at(z, mouth=H.MOUTH_RADIUS, length=H.HORN_LENGTH):
    """Inner radius of the horn / throat at depth z (mouth at z = 0)."""
    z = np.asarray(z, float)
    r = mouth - (mouth - H.THROAT_RADIUS) * np.clip(z / length, 0, 1)
    return np.where(z < 0, np.nan, r)


def mode_list(n=6):
    """The m = 1 modes of a circular guide in order of cut-off: (name, kind, Bessel root)."""
    te = [("TE1%d" % (i + 1), "TE", x) for i, x in enumerate(jnp_zeros(1, n))]
    tm = [("TM1%d" % (i + 1), "TM", x) for i, x in enumerate(jn_zeros(1, n))]
    return sorted(te + tm, key=lambda m: m[2])


MODES = mode_list()


def cutoff_radius(root):
    return root / K


def guide_mode(rho, radius, kind, root):
    """Family profiles (e_r, e_φ, e_z, h_r, h_φ, h_z) of the x-polarized TE1n / TM1n mode heading
    in (+z), transverse parts only, or None below cut-off."""
    kc = root / radius
    beta = np.sqrt(K * K - kc * kc + 0j)
    if beta.imag != 0:
        return None
    beta = beta.real
    u = np.maximum(kc * rho, 1e-9)
    inside = rho < radius
    if kind == "TE":  # E_t from H_z ∝ J1 sin φ;  H_t = (β/k) ẑ × E_t
        e_r, e_p = jv(1, u) / u, -jvp(1, u)
        y = beta / K
    else:  # E_t ∝ ∇_t E_z,  E_z ∝ J1 cos φ;  H_t = (k/β) ẑ × E_t
        e_r, e_p = jvp(1, u), -jv(1, u) / u
        y = K / beta
    e_r, e_p = np.where(inside, e_r, 0.0), np.where(inside, e_p, 0.0)
    z0 = np.zeros_like(e_r)
    return e_r, e_p, z0, -y * e_p, y * e_r, z0  # ẑ × (e_r r̂ + e_φ φ̂) = e_r φ̂ - e_φ r̂


def _integral(s, x, w=None):
    """π ∫ x w ds; on a plane s = ρ and w = ρ (the default)."""
    return np.pi * np.trapezoid(x * (s if w is None else w), s)


def cross(s, f, g, w=None):
    """∫(E_f × H_g*)·n̂ dA for family profiles."""
    return _integral(s, f[0] * np.conj(g[4]) - f[1] * np.conj(g[3]), w)


def reaction(s, f, g, w=None):
    """∫(E_f × H_g)·n̂ dA for family profiles."""
    return _integral(s, f[0] * g[4] - f[1] * g[3], w)


def reverse(m):
    """The same mode heading the other way: its time reverse (E*, -H*)."""
    return tuple(np.conj(c) for c in m[:3]) + tuple(-np.conj(c) for c in m[3:])


def amplitudes(s, F, m, w=None):
    """Amplitudes (a_in, a_out) of the mode m (heading in) and of its reverse in the field F,
    and the power the mode carries at unit amplitude."""
    def R(f, g):
        return reaction(s, f, g, w) - reaction(s, g, f, w)

    back = reverse(m)
    return R(F, back) / R(m, back), R(F, m) / R(back, m), 0.5 * np.real(cross(s, m, m, w))


def project(s, F, m, w=None):
    """Powers (in, out) of the mode m (family profiles, heading in) in the field F."""
    a_in, a_out, P = amplitudes(s, F, m, w)
    return P * abs(a_in) ** 2, P * abs(a_out) ** 2


class HornRun:
    def __init__(self, root, case, emit=None):
        """The launch root/case and its beam alone, root/case_air; emit: the horn's 'emit' run at the
        same resolution, for the part of the field that reaches the guide (needed for a corrugated
        horn)."""
        self.dir = os.path.join(root, case)
        self.air_dir = os.path.join(root, f"{case}_air")
        self.meta = json.load(open(os.path.join(self.dir, "meta.json")))
        self.air = json.load(open(os.path.join(self.air_dir, "meta.json")))
        self.P_inc = self.air["P_net_refl_plane"]
        self.name = case
        self.mouth, self.length = self.meta["mouth_radius"], self.meta["horn_length"]
        self.theta0 = np.radians(self.meta["flare_deg"])
        self.z_apex = self.mouth / np.tan(self.theta0)  # where the cone would close
        self.wall = self.meta.get("wall", "smooth")
        self.focus = self.meta.get("focus", 0.0)
        self.z_source = self.meta["z_source"]
        self.emit_dir = emit
        if emit is not None:
            e = json.load(open(os.path.join(emit, "meta.json")))
            assert e["res"] == self.meta["res"] and e.get("wall", "smooth") == self.wall, (emit, self.dir)

    def radius(self, z):
        return radius_at(z, self.mouth, self.length)

    def fields(self, air=False):
        """r, z and the family profiles of the steady field, each indexed (z, r).  The beam alone
        (air) is put on the launch's grid; its box may stop just behind the mouth, and is zero
        beyond."""
        if not air and hasattr(self, "_f"):
            return self._f
        if air and hasattr(self, "_fa"):
            return self._fa
        d = np.load(os.path.join(self.air_dir if air else self.dir, "dft_rz.npz"))
        f = d["r"], d["z"], beam.from_meep([d[c].astype(np.complex128) for c in COMPS])
        if not air:
            self._f = f
            return f
        r, z, F = self.fields()
        k = _offset(z, f[1])
        assert len(f[0]) == len(r) and np.allclose(f[0], r), "the beam alone is on another grid"
        full = [np.zeros_like(c) for c in F]
        for a, g in zip(full, f[2]):
            a[k:k + len(f[1])] = g
        self._fa = r, z, full
        return self._fa

    # ---- where the power goes --------------------------------------------------------
    def landing(self):
        """Where the incident beam lands: its power through the plane just in front of the face,
        inside the mouth and on the metal around it, as fractions of the incident power."""
        r, z, G = self.fields(air=True)
        k = int(np.max(np.where(z < 0)[0]))
        g = [c[k] for c in G]
        inside = r < self.mouth
        p_in = beam.power(r[inside], [c[inside] for c in g]) / self.P_inc
        p_all = beam.power(r, g) / self.P_inc
        return dict(z=float(z[k]), inside=float(p_in), face=float(p_all - p_in), total=float(p_all))

    def budget(self):
        m = self.meta
        throat = [r["flux"] for r in m["flux_planes"] if r["z"] > self.length]
        guided = float(np.mean(throat)) / self.P_inc
        back = 1 - m["P_net_refl_plane"] / self.P_inc
        face = m["P_out_side_wall"] / self.P_inc
        horn = [r["flux"] / self.P_inc for r in m["flux_planes"]]
        return dict(guided=guided, reflected=back + face, back=back, along_face=face, total=guided + back + face,
                    horn_flux=horn, horn_flux_z=[r["z"] for r in m["flux_planes"]])

    # ---- the field on the surfaces the modes are defined on ------------------------------
    def cap_radius(self, z_edge):
        """Distance from the apex of the cap that meets the wall at depth z_edge."""
        return (self.z_apex - z_edge) / np.cos(self.theta0)

    def surface(self, z_edge, n=600):
        """The steady field on the cap about the apex that meets the wall at depth z_edge (a
        plane in the throat).  Returns the arc length s from the axis, the area weight w
        (dA = w ds dφ), the surface's radius in s, and the family profiles with the r and z
        components replaced by the ones along the surface (θ) and along its normal (inward)."""
        r, z, F = self.fields()
        if z_edge >= self.length:
            k = int(np.argmin(abs(z - z_edge)))
            return r, r, float(self.radius(z[k])), [f[k] for f in F]
        Ra = self.cap_radius(z_edge)
        th = np.linspace(0, self.theta0, n)
        rho, zz = Ra * np.sin(th), self.z_apex - Ra * np.cos(th)
        idx = np.array([(zz - z[0]) / (z[1] - z[0]), (rho - r[0]) / (r[1] - r[0])])
        g = [map_coordinates(f.real, idx, order=1, mode="nearest")
             + 1j * map_coordinates(f.imag, idx, order=1, mode="nearest") for f in F]
        c, s = np.cos(th), np.sin(th)  # θ̂ = cos θ ρ̂ + sin θ ẑ, inward normal -sin θ ρ̂ + cos θ ẑ
        prof = [g[0] * c + g[2] * s, g[1], g[2] * c - g[0] * s, g[3] * c + g[5] * s, g[4], g[5] * c - g[3] * s]
        return Ra * th, rho, Ra * self.theta0, prof

    # ---- the mode split inside ----------------------------------------------------------
    def depths(self, step=None):
        """The depths of the surfaces: one per grid row from just inside the mouth on."""
        _, z, _ = self.fields()
        keep = z[z > 0.05]
        return keep if step is None else keep[::step]

    def mode_powers(self, modes=MODES):
        """For every depth inside the metal: the power heading in and heading out in each mode,
        as fractions of the incident power.  Returns depths, {name: (in, out)}."""
        if hasattr(self, "_mp"):
            return self._mp
        zs = self.depths()
        out = {name: (np.full(len(zs), np.nan), np.full(len(zs), np.nan)) for name, _, _ in modes}
        for j, zk in enumerate(zs):
            s, w, rad, row = self.surface(zk)
            for name, kind, root in modes:
                m = guide_mode(s, rad, kind, root)
                if m is not None:
                    p_in, p_out = project(s, row, m, w)
                    out[name][0][j], out[name][1][j] = p_in / self.P_inc, p_out / self.P_inc
        self._mp = (zs, out)
        return self._mp

    def entry(self, z_at=1.0):
        """On one surface just inside the mouth (where the same modes propagate as on the mouth
        itself): the power that entered, summed over the modes heading in, and the power in
        the modes heading back out; fractions of the incident power."""
        z, mp_ = self.mode_powers()
        k = int(np.argmin(abs(z - z_at)))
        return dict(z=float(z[k]), entered=float(np.nansum([v[0][k] for v in mp_.values()])),
                    back=float(np.nansum([v[1][k] for v in mp_.values()])))

    def te11_in(self):
        """TE11's amplitude heading in on every surface, as a function of the distance from the
        apex (caps) and of depth (throat planes)."""
        zs = self.depths()
        a = np.zeros(len(zs), complex)
        for j, zk in enumerate(zs):
            s, w, rad, row = self.surface(zk)
            m = guide_mode(s, rad, "TE", H.TE11_ROOT)
            a[j] = amplitudes(s, row, m, w)[0]
        return zs, a

    def from_emit(self):
        """The part of the steady field that reaches the guide, from the time reverse (E*, -H*)
        of the emit run's field, scaled so that its TE11 heading in in the throat is the launch's.
        Family profiles on the launch's grid (zero outside the emit run's box)."""
        r, z, F = self.fields()
        d = np.load(os.path.join(self.emit_dir, "dft_rz.npz"))
        E = beam.from_meep([d[c].astype(np.complex128) for c in COMPS])
        k0, nr = _offset(z, d["z"]), len(d["r"])
        assert np.allclose(d["r"], r[:nr]), "the emit run is on another grid"
        T = [np.zeros_like(c) for c in F]
        for t, e, sign in zip(T, E, (1, 1, 1, -1, -1, -1)):
            t[k0:k0 + len(d["z"]), :nr] = sign * np.conj(e)
        k = int(np.argmin(abs(z - (self.length + 0.3 * self.meta["throat_length"]))))
        m = guide_mode(r, float(self.radius(z[k])), "TE", H.TE11_ROOT)
        alpha = amplitudes(r, [f[k] for f in F], m)[0] / amplitudes(r, [t[k] for t in T], m)[0]
        return [alpha * t for t in T]

    def side_split(self, method=None):
        """|E|² on the plane of polarization (y = 0, x = ±r): incident and reflected in front of
        the mouth; the part that reaches the guide (TE11 heading in, in the guide), and the rest,
        behind it.  Inside the horn that part is TE11 heading in on the caps (method 'caps', the
        smooth horn) or the time-reversed emit field ('emit', any horn; the default for a
        corrugated one).  Returns r, z and (z, r) arrays."""
        method = method or ("emit" if self.wall == "corrugated" else "caps")
        if getattr(self, "_side_method", None) == method:
            return self._side
        r, z, F = self.fields()
        _, _, G = self.fields(air=True)

        def e2(f):  # |E|² at φ = 0, where only the cos φ profiles are nonzero
            return abs(f[0]) ** 2 + abs(f[2]) ** 2

        tot, inc = e2(F), e2(G)
        refl = e2([a - b for a, b in zip(F, G)])
        front = z < 0
        Z, Rr = np.meshgrid(z, r, indexing="ij")
        g_r, g_z = np.zeros_like(F[0]), np.zeros_like(F[0])  # the part that reaches the guide, ρ and z parts
        if method == "emit":
            T = self.from_emit()
            horn = (Z > 0) & (Z < self.length)
            g_r[horn], g_z[horn] = T[0][horn], T[2][horn]
        else:
            # each grid point lies on the cap through it; TE11's amplitude, smooth once its
            # exp(-ikR) is taken out, is interpolated in the cap's distance R from the apex
            zs, a = self.te11_in()
            cone = zs < self.length
            Rc = self.cap_radius(zs[cone])
            smooth = a[cone] * np.exp(1j * K * Rc)
            Ra = np.hypot(Rr, self.z_apex - Z)
            th = np.arctan2(Rr, self.z_apex - Z)
            inside = (Z > 0) & (Z < self.length) & (th <= self.theta0)
            Rp = Ra[inside]
            A = (np.interp(Rp, Rc[::-1], smooth.real[::-1]) + 1j * np.interp(Rp, Rc[::-1], smooth.imag[::-1])) \
                * np.exp(-1j * K * Rp)
            u = np.maximum(H.TE11_ROOT * th[inside] / self.theta0, 1e-9)  # TE11's e_θ on each point's cap
            e_t = jv(1, u) / u
            g_r[inside] = A * e_t * np.cos(th[inside])
            g_z[inside] = A * e_t * np.sin(th[inside])
        # in the throat: TE11 heading in, on each plane
        for k in np.where(z >= self.length)[0]:
            m = guide_mode(r, float(self.radius(z[k])), "TE", H.TE11_ROOT)
            g_r[k], g_z[k] = amplitudes(r, [f[k] for f in F], m)[0] * m[0], 0.0
        gui = abs(g_r) ** 2 + abs(g_z) ** 2
        behind = Z > 0
        rest = np.where(behind, abs(F[0] - g_r) ** 2 + abs(F[2] - g_z) ** 2, 0.0)
        gui = np.where(behind, gui, 0.0)
        self._side = dict(r=r, z=z, total=tot, incident=np.where(front[:, None], inc, 0.0),
                          reflected=np.where(front[:, None], refl, 0.0), guided=gui, rest=rest)
        self._side_method = method
        return self._side


def _offset(z, zs):
    """Index in z of the first of zs, a run of the same grid."""
    k = int(np.argmin(abs(z - zs[0])))
    n = min(len(zs), len(z) - k)
    assert np.allclose(z[k:k + n], zs[:n], atol=1e-6), "not the same grid"
    return k
