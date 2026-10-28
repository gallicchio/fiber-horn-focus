"""The first, superseded version of fiberlaunch/horn.py: the same horn on a 3D Cartesian grid.

Kept as the record of why the horn is simulated on the (r, z) plane instead.  Meep cannot
smooth a perfect conductor below the pixel, so on a Cartesian grid the cone's wall is a
staircase, and the staircase scatters TE11 into TM11 on its way down the taper.  The coupling of
the matched beam this code finds (by reciprocity, scripts/horn_design.py on its 'emit' runs)
is 41, 55 and 69 % at 16, 20 and 28 px/λ, still climbing; the body-of-revolution version
converges to 81 %.  A 3D grid fine enough to converge is far out of reach.

  air    the beam alone (the same box, no metal): incident power, and the incident field on
         the side planes, so that the reflected field is (total - incident)
  horn   the launch into the horn
  emit   the reverse: TE11 sent out of the throat, the field it radiates recorded in front of
         the mouth (front_xy.npz)

Run with MPI, e.g.
    mpirun -np 16 python -m fiberlaunch.horn_cartesian --kind emit --res 20 --out data/horn/cartesian/emit_res20
"""

import argparse
import math
import json
import os
import time

import meep as mp
import numpy as np

from . import params as P
from .fdtd import COMP_NAMES, COMPS, MASTER, SYMMETRIES, log, save_dft_plane
from .params import FLARE_DEG, MOUTH_RADIUS, THROAT_LENGTH, THROAT_RADIUS, te11_beta

RESOLUTION = 20  # px/λ; the 3D grid cannot go much finer


class Box:
    """The cell: air in front (z < 0), the metal block with the horn and throat behind."""

    def __init__(self, kind, res, half_width, mouth=MOUTH_RADIUS, flare=FLARE_DEG):
        self.kind, self.res, self.hw = kind, res, half_width
        self.A, self.flare = mouth, flare
        self.L = (mouth - THROAT_RADIUS) / math.tan(math.radians(flare))
        self.sxy = 2 * (half_width + P.DPML)
        self.z_lo = P.Z_SOURCE - P.GAP_BEHIND_SOURCE - P.DPML
        self.z_hi = self.L + THROAT_LENGTH + P.DPML_END
        self.sz = self.z_hi - self.z_lo
        self.zc = 0.5 * (self.z_lo + self.z_hi)
        self.z_in_lo = self.z_lo + P.DPML
        self.z_in_hi = self.z_hi - P.DPML_END

    def snap(self, z):
        d = 1.0 / self.res
        return (np.floor((z - self.zc) / d) + 0.5) * d + self.zc

    def v(self, x=0.0, y=0.0, z=0.0):
        return mp.Vector3(x, y, z - self.zc)

    @property
    def cell(self):
        return mp.Vector3(self.sxy, self.sxy, self.sz)

    def pml(self):
        return [mp.PML(P.DPML, direction=mp.X), mp.PML(P.DPML, direction=mp.Y),
                mp.PML(P.DPML, direction=mp.Z, side=mp.Low),
                mp.PML(P.DPML_END, direction=mp.Z, side=mp.High)]

    def geometry(self):
        if self.kind == "air":
            return []
        depth = self.z_hi
        throat = self.z_hi - self.L
        return [
            mp.Block(center=self.v(z=0.5 * depth), size=mp.Vector3(mp.inf, mp.inf, depth),
                     material=mp.metal),
            # Cone: 'radius' at center - height/2 (the mouth), 'radius2' at center + height/2
            mp.Cone(center=self.v(z=0.5 * self.L), height=self.L, axis=mp.Vector3(0, 0, 1),
                    radius=self.A, radius2=THROAT_RADIUS, material=mp.air),
            mp.Cylinder(center=self.v(z=self.L + 0.5 * throat), height=throat, radius=THROAT_RADIUS,
                        axis=mp.Vector3(0, 0, 1), material=mp.air),
        ]

    def xy_plane(self, z, half=None):
        h = self.hw if half is None else half
        return dict(center=self.v(z=self.snap(z)), size=mp.Vector3(2 * h, 2 * h, 0))

    def side_plane(self, normal):
        zlen = self.z_in_hi - self.z_in_lo
        zmid = 0.5 * (self.z_in_lo + self.z_in_hi)
        size = mp.Vector3(2 * self.hw, 0, zlen) if normal == "y" else mp.Vector3(0, 2 * self.hw, zlen)
        return dict(center=self.v(z=zmid), size=size)

    def side_walls(self):
        """The four sides of the air gap between the reflection plane and the metal face:
        light scattered along the face leaves through these without crossing that plane."""
        zr = self.snap(P.Z_REFL)
        h = 0.5 * (0.0 - zr)
        c = 0.5 * (zr + 0.0)
        out = []
        for sgn in (1, -1):
            out.append(mp.FluxRegion(center=self.v(x=sgn * self.hw, z=c), size=mp.Vector3(0, 2 * self.hw, 2 * h),
                                     direction=mp.X, weight=sgn))
            out.append(mp.FluxRegion(center=self.v(y=sgn * self.hw, z=c), size=mp.Vector3(2 * self.hw, 0, 2 * h),
                                     direction=mp.Y, weight=sgn))
        return out


def steady_time(box):
    """The light reaches the end of the throat at about 0.6 c (TE11 near cut-off), and the
    higher modes that turn back inside the horn linger near their turning points."""
    ramp = P.TURN_ON * (3.0 + 5.0)
    return ramp + 2.5 * (box.z_hi - P.Z_SOURCE) + 20.0


def run(kind, w0, res, half_width, out, movie=False, xy_planes=False, extra_periods=0.0,
        mouth=MOUTH_RADIUS, flare=FLARE_DEG):
    t_start = time.time()
    box = Box(kind, res, half_width, mouth, flare)
    if kind == "emit":  # a small x-polarized current sheet in the throat: TE11 a few λ away
        src = mp.Source(mp.ContinuousSource(frequency=P.FREQ, width=P.TURN_ON), component=mp.Ex,
                        center=box.v(z=box.L + 0.6 * THROAT_LENGTH), size=mp.Vector3(0.3, 0.3, 0))
    else:
        src = mp.GaussianBeam3DSource(
            mp.ContinuousSource(frequency=P.FREQ, width=P.TURN_ON),
            center=box.v(z=P.Z_SOURCE), size=mp.Vector3(2 * box.hw, 2 * box.hw, 0),
            beam_x0=mp.Vector3(0, 0, 0.0 - P.Z_SOURCE),  # focus on the mouth
            beam_kdir=mp.Vector3(0, 0, 1), beam_w0=w0, beam_E0=mp.Vector3(1, 0, 0))
    sim = mp.Simulation(cell_size=box.cell, resolution=res, geometry=box.geometry(), sources=[src],
                        boundary_layers=box.pml(), symmetries=SYMMETRIES, default_material=mp.air)
    if MASTER:
        os.makedirs(out, exist_ok=True)
    T_ss = steady_time(box) + extra_periods
    dt = sim.Courant / res
    meta = dict(kind=kind, w0=w0, beam_na=P.beam_na(w0) if w0 else None, res=res, half_width=half_width,
                throat_radius=THROAT_RADIUS, mouth_radius=box.A, flare_deg=box.flare,
                horn_length=box.L, throat_length=THROAT_LENGTH, te11_beta_throat=te11_beta(THROAT_RADIUS),
                z_lo=box.z_lo, z_hi=box.z_hi, z_source=P.Z_SOURCE, z_refl=P.Z_REFL,
                nproc=mp.count_processors(), T_steady=T_ss, dt=dt, dft_periods=P.DFT_PERIODS)
    log(f"[{kind}] w0={w0} res={res} cell={box.cell} horn L={box.L:.2f} T_ss={T_ss:.1f}")

    step_funcs = []
    sim.init_sim()
    vol = mp.Volume(**box.side_plane("y"))
    if movie:
        steps_per_period = round(1 / dt)
        every = max(1, round(0.34 * steps_per_period))
        n_frames = int((T_ss + P.DFT_PERIODS) / (every * dt)) + 2
        probe = sim.get_array(component=mp.Ex, vol=vol)
        mm = (np.lib.format.open_memmap(os.path.join(out, "frames_xz_Ex.npy"), mode="w+",
                                        dtype=np.float16, shape=(n_frames,) + probe.shape) if MASTER else None)
        times = []

        def grab(sim):
            if sim.timestep() % every == 0 and len(times) < n_frames:
                a = sim.get_array(component=mp.Ex, vol=vol)
                if MASTER:
                    mm[len(times)] = a.astype(np.float16)
                times.append(sim.meep_time())

        step_funcs.append(grab)
        meta.update(frame_every_steps=every, frame_dt=every * dt)
    # where the metal is, for drawing (Meep reports PEC as a large negative permittivity)
    eps = sim.get_array(component=mp.Dielectric, vol=vol)
    x, _, z, _ = sim.get_array_metadata(vol=vol)
    if MASTER:
        np.savez(os.path.join(out, "eps_xz.npz"), eps=eps, x=np.asarray(x), z=np.asarray(z) + box.zc)

    sim.run(*step_funcs, until=T_ss)
    log(f"[{kind}] steady at t={sim.meep_time():.2f}, wall {time.time() - t_start:.0f} s")

    dfa = dict(decimation_factor=1)
    refl = sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.xy_plane(P.Z_REFL)), **dfa)
    walls = sim.add_flux(P.FREQ, 0, 1, *box.side_walls(), **dfa)
    flux_z = [] if kind == "air" else [z for z in (0.5, 0.25 * box.L, 0.5 * box.L,
                                                   0.75 * box.L, box.L + 1.0,
                                                   box.L + 0.5 * THROAT_LENGTH)]
    if kind == "emit":
        flux_z = [0.5, box.L + 0.3 * THROAT_LENGTH]
        front = sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False, **box.xy_plane(-0.5), **dfa)
    flux_planes = {z: sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.xy_plane(z, box.A + 0.5)), **dfa)
                   for z in flux_z}
    sides = {k: sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False, **box.side_plane(k), **dfa)
             for k in ("y", "x")}
    mouth = sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False, **box.xy_plane(0.0), **dfa)
    inside = None
    if xy_planes and kind == "horn":  # every voxel of the horn and throat, for the mode split
        h = box.A + 0.25
        inside = sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False,
                                    center=box.v(z=0.5 * box.z_in_hi), size=mp.Vector3(2 * h, 2 * h, box.z_in_hi),
                                    **dfa)
    sim.run(*step_funcs, until=P.DFT_PERIODS)
    log(f"[{kind}] DFT done at t={sim.meep_time():.2f}, wall {time.time() - t_start:.0f} s")

    meta["P_net_refl_plane"] = mp.get_fluxes(refl)[0]
    meta["P_out_side_walls"] = mp.get_fluxes(walls)[0]
    meta["flux_planes"] = [dict(z=float(box.snap(z)), flux=mp.get_fluxes(f)[0]) for z, f in flux_planes.items()]
    for row in meta["flux_planes"]:
        log(f"  z={row['z']:6.2f}  net flux {row['flux']:.5f}")
    log(f"  reflection plane {meta['P_net_refl_plane']:.5f}  side walls {meta['P_out_side_walls']:.3e}")
    for k, d in sides.items():
        save_dft_plane(sim, d, os.path.join(out, f"dft_side_{'xz' if k == 'y' else 'yz'}.npz"), box.zc)
    save_dft_plane(sim, mouth, os.path.join(out, "mouth_xy.npz"), box.zc)
    if kind == "emit":
        save_dft_plane(sim, front, os.path.join(out, "front_xy.npz"), box.zc)
    if inside is not None:
        x, y, z, _ = sim.get_array_metadata(dft_cell=inside)
        arrs = {}
        for c, n in zip(COMPS, COMP_NAMES):
            a = sim.get_dft_array(inside, c, 0)
            if MASTER:
                arrs[n] = a.astype(np.complex64)
        if MASTER:
            np.savez(os.path.join(out, "inside_xyz.npz"), x=np.asarray(x), y=np.asarray(y),
                     z=np.asarray(z) + box.zc, **arrs)

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
    ap.add_argument("--res", type=int, default=RESOLUTION)
    ap.add_argument("--half-width", type=float, default=P.HALF_WIDTH)
    ap.add_argument("--movie", action="store_true", help="record Ex(x, z, t) frames")
    ap.add_argument("--xy-planes", action="store_true", help="record the field in every voxel of the horn and throat")
    ap.add_argument("--extra-periods", type=float, default=0.0, help="run this much longer before the DFT")
    ap.add_argument("--mouth", type=float, default=MOUTH_RADIUS, help="horn mouth radius, in wavelengths")
    ap.add_argument("--flare", type=float, default=FLARE_DEG, help="horn half-angle, in degrees")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    w0 = a.w0 if a.w0 is not None else (P.CASES[a.case]["w0"] if a.case else None)
    run(a.kind, w0, a.res, a.half_width, a.out, movie=a.movie, xy_planes=a.xy_planes,
        extra_periods=a.extra_periods, mouth=a.mouth, flare=a.flare)


if __name__ == "__main__":
    main()
