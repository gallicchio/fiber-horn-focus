"""FDTD (Meep) of the same three Gaussian beams launched into a single-mode metal waveguide.

A perfectly conducting circular waveguide is single mode only while it is less than 0.77 λ
across (TE11 cuts off at radius 0.293 λ, TM01 propagates above 0.383 λ), far smaller than
the fiber's mode, so no free-space beam can be focused down to match it.  Microwave practice
feeds it through a horn, which is what this models: a smooth conical horn, flush with the face
of a thick metal block, flaring from the single-mode throat to a mouth about the size of the
fiber's core, so that the matched beam of the fiber runs (NA 0.19) is also the matched beam here.

The horn, the guide and the x-polarized beam on their axis are all symmetric about that axis,
so the fields are exactly cos φ / sin φ in angle: half angular order m = +1, half m = -1, and the
two halves are mirror images (fiberlaunch/beam.py).  Meep solves the full 3D Maxwell equations
for the m = +1 half on the (r, z) plane, which is exact for this geometry and costs a 2D run.
That matters here because Meep cannot smooth a perfect conductor below the pixel: the cone's
staircased wall converts TE11 into TM11 at any resolution a 3D Cartesian grid can afford
(fiberlaunch/horn_cartesian.py: 41, 55, 69 % coupling of the matched beam at 16, 20, 28 px/λ),
while on the (r, z) plane the staircase can be made as fine as needed.

  air    the beam alone (the same box, no metal): incident power, and the incident field, so
         that the reflected field is (total - incident)
  horn   the launch into the horn
  emit   the reverse: TE11 sent out of the throat, the field it radiates recorded in front of
         the mouth.  By reciprocity its overlap with any incoming beam is that beam's coupling
         into the throat (scripts/horn_design.py)

The beam enters through a plane of equivalent electric and magnetic currents (J = ẑ × H,
M = -ẑ × E of the exact beam, fiberlaunch/beam.py), which radiates it forward only.

Run with MPI, e.g.
    mpirun -np 8 python -m fiberlaunch.horn --kind horn --case narrow --movie
"""

import argparse
import json
import math
import os
import time

import meep as mp
import numpy as np

from . import beam
from . import params as P
from .params import FLARE_DEG, HORN_RESOLUTION, MOUTH_RADIUS, THROAT_LENGTH, THROAT_RADIUS, te11_beta

MASTER = mp.am_master()
COMPS = [mp.Er, mp.Ep, mp.Ez, mp.Hr, mp.Hp, mp.Hz]
COMP_NAMES = ["Er", "Ep", "Ez", "Hr", "Hp", "Hz"]


def log(*a):
    if MASTER:
        print(*a, flush=True)


class Box:
    """The (r, z) cell: air in front (z < 0), the metal block with the horn and throat behind."""

    def __init__(self, kind, res, half_width, mouth=MOUTH_RADIUS, flare=FLARE_DEG):
        self.kind, self.res, self.hw = kind, res, half_width
        self.A, self.flare = mouth, flare
        self.L = (mouth - THROAT_RADIUS) / math.tan(math.radians(flare))
        self.sr = half_width + P.DPML
        self.z_lo = P.Z_SOURCE - P.GAP_BEHIND_SOURCE - P.DPML
        # an even number of pixels in z, so that every edge below lies between pixels
        n = 2 * math.ceil(0.5 * res * (self.L + THROAT_LENGTH + P.DPML_END - self.z_lo))
        self.sz = n / res
        self.z_hi = self.z_lo + self.sz
        self.zc = 0.5 * (self.z_lo + self.z_hi)
        self.z_in_lo = self.z_lo + P.DPML
        self.z_in_hi = self.z_hi - P.DPML_END

    def snap(self, z):
        d = 1.0 / self.res
        return (np.floor((z - self.zc) / d) + 0.5) * d + self.zc

    def v(self, r=0.0, z=0.0):
        return mp.Vector3(r, 0, z - self.zc)

    @property
    def cell(self):
        return mp.Vector3(self.sr, 0, self.sz)

    def pml(self):
        return [mp.PML(P.DPML, direction=mp.R),
                mp.PML(P.DPML, direction=mp.Z, side=mp.Low),
                mp.PML(P.DPML_END, direction=mp.Z, side=mp.High)]

    def geometry(self):
        if self.kind == "air":
            return []
        depth = self.z_hi + 1.0
        rt = THROAT_RADIUS
        # the cone as a prism over its (r, z) outline, one unit thick in the unused φ direction
        horn = [self.v(r, z) - mp.Vector3(0, 0.5, 0)
                for r, z in ((-self.A, 0.0), (self.A, 0.0), (rt, self.L), (-rt, self.L))]
        return [
            mp.Block(center=self.v(0, 0.5 * depth), size=mp.Vector3(2 * self.sr + 2, mp.inf, depth),
                     material=mp.metal),
            mp.Prism(horn, height=1.0, axis=mp.Vector3(0, 1, 0), material=mp.air),
            mp.Block(center=self.v(0, self.L + 0.5 * (depth - self.L)), size=mp.Vector3(2 * rt, mp.inf, depth - self.L),
                     material=mp.air),
        ]

    def z_plane(self, z, half=None):
        """The disc r < half at depth z (as a center/size pair)."""
        h = self.hw if half is None else half
        return dict(center=self.v(0.5 * h, self.snap(z)), size=mp.Vector3(h, 0, 0))

    def rz_plane(self):
        zlen = self.z_in_hi - self.z_in_lo
        return dict(center=self.v(0.5 * self.hw, 0.5 * (self.z_in_lo + self.z_in_hi)),
                    size=mp.Vector3(self.hw, 0, zlen))

    def rz_grid(self, sim):
        """Coordinates of the arrays Meep returns for rz_plane(), which are indexed (z, r) at
        pixel centres, one pixel past each edge (get_array_metadata is not available in
        cylindrical coordinates)."""
        d = sim.get_array_slice_dimensions(vol=mp.Volume(**self.rz_plane()))
        nz, nr = d.dim_sizes
        lo, hi = d.min_corner, d.max_corner
        return np.linspace(lo.x, hi.x, nr), np.linspace(lo.z, hi.z, nz) + self.zc

    def side_wall(self):
        """The cylinder r = half-width between the reflection plane and the metal face: light
        scattered along the face leaves through it without crossing that plane."""
        zr = self.snap(P.Z_REFL)
        return mp.FluxRegion(center=self.v(self.hw, 0.5 * zr), size=mp.Vector3(0, 0, -zr), direction=mp.R)


def steady_time(box):
    """The light reaches the end of the throat at about 0.6 c (TE11 near cut-off), and the
    higher modes that turn back inside the horn linger near their turning points.  The beam in
    air is only used in front of the mouth, which it crosses once."""
    ramp = P.TURN_ON * (3.0 + 5.0)
    if box.kind == "air":
        return ramp - P.Z_SOURCE + 20.0
    return ramp + 2.5 * (box.z_hi - P.Z_SOURCE) + 20.0


def beam_sources(box, w0, src_time):
    """Equivalent currents on the plane z = Z_SOURCE that radiate the m = +1 half of the beam
    focused on the mouth (z = 0) forward: J = ẑ × H, M = -ẑ × E."""
    rho = np.linspace(0, box.hw, int(box.hw * 400) + 1)
    Er, Ep, _, Hr, Hp, _ = beam.to_meep(beam.family(rho, P.Z_SOURCE, w0))
    currents = {mp.Er: -Hp, mp.Ep: Hr, mp.Hr: Ep, mp.Hp: -Er}
    c = box.z_plane(P.Z_SOURCE)

    def amp(prof):
        re, im = np.real(prof), np.imag(prof)
        return lambda p: complex(np.interp(p.x + 0.5 * box.hw, rho, re), np.interp(p.x + 0.5 * box.hw, rho, im))

    return [mp.Source(src_time, component=comp, amp_func=amp(prof), **c) for comp, prof in currents.items()]


def run(kind, w0, res, half_width, out, movie=False, extra_periods=0.0, mouth=MOUTH_RADIUS, flare=FLARE_DEG,
        frame_res=32):
    t_start = time.time()
    box = Box(kind, res, half_width, mouth, flare)
    src_time = mp.ContinuousSource(frequency=P.FREQ, width=P.TURN_ON)
    if kind == "emit":  # a small x-polarized current sheet in the throat: TE11 a few λ away
        z = box.L + 0.6 * THROAT_LENGTH
        sources = [mp.Source(src_time, component=mp.Er, center=box.v(0.075, z), size=mp.Vector3(0.15, 0, 0)),
                   mp.Source(src_time, component=mp.Ep, center=box.v(0.075, z), size=mp.Vector3(0.15, 0, 0),
                             amplitude=1j)]
    else:
        sources = beam_sources(box, w0, src_time)
    sim = mp.Simulation(cell_size=box.cell, resolution=res, geometry=box.geometry(), sources=sources,
                        boundary_layers=box.pml(), dimensions=mp.CYLINDRICAL, m=1, default_material=mp.air)
    if MASTER:
        os.makedirs(out, exist_ok=True)
    T_ss = steady_time(box) + extra_periods
    dt = sim.Courant / res
    meta = dict(kind=kind, w0=w0, beam_na=P.beam_na(w0) if w0 else None, res=res, half_width=half_width,
                throat_radius=THROAT_RADIUS, mouth_radius=box.A, flare_deg=box.flare,
                horn_length=box.L, throat_length=THROAT_LENGTH, te11_beta_throat=te11_beta(THROAT_RADIUS),
                z_lo=box.z_lo, z_hi=box.z_hi, z_source=P.Z_SOURCE, z_refl=P.Z_REFL,
                nproc=mp.count_processors(), T_steady=T_ss, dt=dt, dft_periods=P.DFT_PERIODS,
                fields="Meep m = +1 half; the linearly polarized field is beam.from_meep() of it")
    log(f"[{kind}] w0={w0} res={res} cell={box.cell} horn L={box.L:.2f} T_ss={T_ss:.1f}")

    step_funcs = []
    sim.init_sim()
    vol = mp.Volume(**box.rz_plane())
    r, z = box.rz_grid(sim)
    # where the metal is, for drawing (Meep reports PEC as a large negative permittivity)
    eps = sim.get_array(component=mp.Dielectric, vol=vol)
    assert eps.shape == (len(z), len(r)), (eps.shape, len(z), len(r))
    if MASTER:
        np.savez(os.path.join(out, "eps_rz.npz"), eps=eps, r=r, z=z)
    if movie:  # Re Er on the (r, z) plane, which is Ex on the plane of polarization (x = ±r)
        stride = max(1, res // frame_res)
        steps_per_period = round(1 / dt)
        every = max(1, round(0.34 * steps_per_period))
        n_frames = int((T_ss + P.DFT_PERIODS) / (every * dt)) + 2
        shape = (len(z[::stride]), len(r[::stride]))
        mm = (np.lib.format.open_memmap(os.path.join(out, "frames_rz_Er.npy"), mode="w+",
                                        dtype=np.float16, shape=(n_frames,) + shape) if MASTER else None)
        times = []

        def grab(sim):
            if sim.timestep() % every == 0 and len(times) < n_frames:
                a = sim.get_array(component=mp.Er, vol=vol)
                if MASTER:
                    mm[len(times)] = np.real(a[::stride, ::stride]).astype(np.float16)
                times.append(sim.meep_time())

        step_funcs.append(grab)
        meta.update(frame_every_steps=every, frame_dt=every * dt, frame_stride=stride)
        if MASTER:
            np.savez(os.path.join(out, "frames_grid.npz"), r=r[::stride], z=z[::stride])

    sim.run(*step_funcs, until=T_ss)
    log(f"[{kind}] steady at t={sim.meep_time():.2f}, wall {time.time() - t_start:.0f} s")

    dfa = dict(decimation_factor=1)
    refl = sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.z_plane(P.Z_REFL)), **dfa)
    wall = sim.add_flux(P.FREQ, 0, 1, box.side_wall(), **dfa)
    if kind == "air":
        flux_z = []
    elif kind == "emit":
        flux_z = [0.5, box.L + 0.3 * THROAT_LENGTH]
    else:
        flux_z = [0.5, 0.25 * box.L, 0.5 * box.L, 0.75 * box.L, box.L + 1.0, box.L + 0.5 * THROAT_LENGTH]
    flux_planes = {z: sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.z_plane(z, box.A + 0.5)), **dfa)
                   for z in flux_z}
    plane = sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False, **box.rz_plane(), **dfa)
    sim.run(*step_funcs, until=P.DFT_PERIODS)
    log(f"[{kind}] DFT done at t={sim.meep_time():.2f}, wall {time.time() - t_start:.0f} s")

    meta["P_net_refl_plane"] = mp.get_fluxes(refl)[0]
    meta["P_out_side_wall"] = mp.get_fluxes(wall)[0]
    meta["flux_planes"] = [dict(z=float(box.snap(z)), flux=mp.get_fluxes(f)[0]) for z, f in flux_planes.items()]
    for row in meta["flux_planes"]:
        log(f"  z={row['z']:6.2f}  net flux {row['flux']:.5f}")
    log(f"  reflection plane {meta['P_net_refl_plane']:.5f}  side wall {meta['P_out_side_wall']:.3e}")
    arrs = {}
    for c, n in zip(COMPS, COMP_NAMES):
        a = sim.get_dft_array(plane, c, 0)
        if MASTER:
            arrs[n] = a.astype(np.complex64)
    if MASTER:
        assert arrs["Er"].shape == (len(z), len(r)), arrs["Er"].shape
        np.savez(os.path.join(out, "dft_rz.npz"), r=r, z=z, **arrs)
    if movie and MASTER:
        np.save(os.path.join(out, "frames_t.npy"), np.array(times))
        mm.flush()
    meta["wall_seconds"] = time.time() - t_start
    if MASTER:
        with open(os.path.join(out, "meta.json"), "w") as f:
            json.dump(meta, f, indent=1)
    log(f"[{kind}] finished, wall {meta['wall_seconds']:.0f} s")
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=["horn", "air", "emit"], default="horn")
    ap.add_argument("--case", choices=list(P.CASES), help="one of the named launch conditions")
    ap.add_argument("--w0", type=float, help="beam waist on the mouth (overrides --case)")
    ap.add_argument("--res", type=int, default=HORN_RESOLUTION)
    ap.add_argument("--half-width", type=float, default=P.HALF_WIDTH)
    ap.add_argument("--movie", action="store_true", help="record Re Er(r, z, t) frames")
    ap.add_argument("--extra-periods", type=float, default=0.0, help="run this much longer before the DFT")
    ap.add_argument("--mouth", type=float, default=MOUTH_RADIUS, help="horn mouth radius, in wavelengths")
    ap.add_argument("--flare", type=float, default=FLARE_DEG, help="horn half-angle, in degrees")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    w0 = a.w0 if a.w0 is not None else (P.CASES[a.case]["w0"] if a.case else None)
    run(a.kind, w0, a.res, a.half_width, a.out, movie=a.movie, extra_periods=a.extra_periods,
        mouth=a.mouth, flare=a.flare)


if __name__ == "__main__":
    main()
