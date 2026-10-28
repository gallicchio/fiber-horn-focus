"""3D FDTD (Meep) of a Gaussian beam focused onto the cleaved end of a single-mode fiber.

Three kinds of run share one box (``params``), so their DFT planes land on one grid:

  air    the beam alone, no fiber: incident power and the incident field on the facet
  mode   a fiber with an eigenmode source: the HE11 mode on the same grid, for projecting
  fiber  the launch itself

A run is CW: the source ramps on, the run continues until the light has reached the end
of the fiber and everything is steady, and then DFT monitors are added and accumulate
over a whole number of periods.  The DFT of a steady sinusoid over whole periods is
exactly its complex amplitude (times a constant that cancels in every ratio we report).

Run with MPI, e.g.
    mpirun -np 16 python -m fiberlaunch.fdtd --case narrow --movie --xy-planes
"""

import argparse
import json
import os
import time

import meep as mp
import numpy as np

from . import fibermode
from . import params as P

MASTER = mp.am_master()


def log(*a):
    if MASTER:
        print(*a, flush=True)


class Box:
    """Geometry of the computational cell; converts physical z (facet at 0) to Meep z."""

    def __init__(self, kind, res, fiber_length, half_width):
        self.kind, self.res, self.L, self.hw = kind, res, fiber_length, half_width
        self.sxy = 2 * (half_width + P.DPML)
        if kind == "air":
            self.z_lo = P.Z_SOURCE - P.GAP_BEHIND_SOURCE - P.DPML
            self.z_hi = 1.0 + P.DPML
        elif kind == "mode":
            self.z_lo = -P.DPML_END
            self.z_hi = 12.0 + P.DPML_END
        else:
            self.z_lo = P.Z_SOURCE - P.GAP_BEHIND_SOURCE - P.DPML
            self.z_hi = fiber_length + P.DPML_END
        self.sz = self.z_hi - self.z_lo
        self.zc = 0.5 * (self.z_lo + self.z_hi)
        # first and last physical z outside the PML
        self.z_in_lo = self.z_lo + (P.DPML if kind != "mode" else P.DPML_END)
        self.z_in_hi = self.z_hi - (P.DPML_END if kind != "air" else P.DPML)

    def snap(self, z):
        """Nearest voxel-centre z.  A DFT plane there gives E and H consistent with Meep's
        own flux to 1e-5; at a voxel boundary one of them is interpolated across the plane
        and the hand-computed flux is off by cos(k dz / 2)."""
        d = 1.0 / self.res
        zm = z - self.zc
        return (np.floor(zm / d) + 0.5) * d + self.zc

    def v(self, x=0.0, y=0.0, z=0.0):
        return mp.Vector3(x, y, z - self.zc)

    @property
    def cell(self):
        return mp.Vector3(self.sxy, self.sxy, self.sz)

    def pml(self):
        if self.kind == "mode":
            return [mp.PML(P.DPML, direction=mp.X), mp.PML(P.DPML, direction=mp.Y),
                    mp.PML(P.DPML_END, direction=mp.Z)]
        hi = P.DPML if self.kind == "air" else P.DPML_END
        return [mp.PML(P.DPML, direction=mp.X), mp.PML(P.DPML, direction=mp.Y),
                mp.PML(P.DPML, direction=mp.Z, side=mp.Low),
                mp.PML(hi, direction=mp.Z, side=mp.High)]

    def geometry(self):
        if self.kind == "air":
            return []
        z0 = self.z_lo if self.kind == "mode" else 0.0
        length = self.z_hi - z0
        zmid = z0 + 0.5 * length
        return [
            mp.Block(center=self.v(z=zmid), size=mp.Vector3(mp.inf, mp.inf, length),
                     material=mp.Medium(index=P.N_CLAD)),
            mp.Cylinder(radius=P.CORE_RADIUS, height=length, axis=mp.Vector3(0, 0, 1),
                        center=self.v(z=zmid), material=mp.Medium(index=P.N_CORE)),
        ]

    def xy_plane(self, z):
        return dict(center=self.v(z=self.snap(z)), size=mp.Vector3(2 * self.hw, 2 * self.hw, 0))

    def mode_window(self, z):
        """Small plane around the core: holds all of HE11 and keeps MPB cheap."""
        return dict(center=self.v(z=self.snap(z)), size=mp.Vector3(2 * P.MODE_HALF_WIDTH, 2 * P.MODE_HALF_WIDTH, 0))

    def side_plane(self, normal):
        """x-z plane (normal 'y') or y-z plane (normal 'x') through the axis, PML excluded."""
        zlen = self.z_in_hi - self.z_in_lo
        zmid = 0.5 * (self.z_in_lo + self.z_in_hi)
        if normal == "y":
            return dict(center=self.v(z=zmid), size=mp.Vector3(2 * self.hw, 0, zlen))
        return dict(center=self.v(z=zmid), size=mp.Vector3(0, 2 * self.hw, zlen))


SYMMETRIES = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=+1)]  # x-polarized, on axis
COMPS = [mp.Ex, mp.Ey, mp.Ez, mp.Hx, mp.Hy]
COMP_NAMES = ["Ex", "Ey", "Ez", "Hx", "Hy"]


def steady_time(box):
    """Time by which the CW field is steady everywhere in the box."""
    ramp = P.TURN_ON * (3.0 + 5.0)  # tanh ramp of ContinuousSource(width, slowness=3)
    return ramp + 1.5 * (box.z_hi - P.Z_SOURCE) + 10.0


def make_sim(box, w0):
    if box.kind == "mode":
        sources = [mp.EigenModeSource(
            mp.ContinuousSource(frequency=P.FREQ, width=P.TURN_ON),
            center=box.v(z=box.z_in_lo + 0.5), size=mp.Vector3(2 * box.hw, 2 * box.hw, 0),
            direction=mp.Z, eig_band=1, eig_parity=mp.EVEN_Y,
            eig_kpoint=mp.Vector3(0, 0, fibermode.he11_neff() * P.FREQ))]
    else:
        sources = [mp.GaussianBeam3DSource(
            mp.ContinuousSource(frequency=P.FREQ, width=P.TURN_ON),
            center=box.v(z=P.Z_SOURCE), size=mp.Vector3(2 * box.hw, 2 * box.hw, 0),
            beam_x0=mp.Vector3(0, 0, 0.0 - P.Z_SOURCE),  # focus on the facet
            beam_kdir=mp.Vector3(0, 0, 1), beam_w0=w0, beam_E0=mp.Vector3(1, 0, 0))]
    return mp.Simulation(cell_size=box.cell, resolution=box.res, geometry=box.geometry(),
                         sources=sources, boundary_layers=box.pml(), symmetries=SYMMETRIES,
                         default_material=mp.air)


def save_dft_plane(sim, dft, path, zc):
    """Gather a DFT plane (collective) and save it on the master as complex64, with
    coordinates in physical z (facet at 0)."""
    arrs = {}
    for c, n in zip(COMPS, COMP_NAMES):
        a = sim.get_dft_array(dft, c, 0)
        if MASTER:
            arrs[n] = a.astype(np.complex64)
    x, y, z, _ = sim.get_array_metadata(dft_cell=dft)
    if MASTER:
        np.savez(path, x=np.asarray(x), y=np.asarray(y), z=np.asarray(z) + zc, **arrs)


def run(kind, w0, res, fiber_length, half_width, out, movie=False, xy_planes=False,
        mode_planes=None, check_steady=False, frame_every=None):
    t_start = time.time()
    box = Box(kind, res, fiber_length, half_width)
    sim = make_sim(box, w0)
    if MASTER:
        os.makedirs(out, exist_ok=True)
    meta = dict(kind=kind, w0=w0, beam_na=P.beam_na(w0) if w0 else None, res=res,
                fiber_length=fiber_length, half_width=half_width,
                n_core=P.N_CORE, n_clad=P.N_CLAD, core_radius=P.CORE_RADIUS,
                z_lo=box.z_lo, z_hi=box.z_hi, dpml=P.DPML, dpml_end=P.DPML_END,
                z_source=P.Z_SOURCE, z_refl=P.Z_REFL, nproc=mp.count_processors())
    T_ss = steady_time(box)
    dt = sim.Courant / res
    meta.update(T_steady=T_ss, dt=dt, dft_periods=P.DFT_PERIODS)
    log(f"[{kind}] w0={w0} res={res} cell={box.cell} T_ss={T_ss:.1f}")

    step_funcs = []
    if movie:
        # instantaneous Ex on the x-z plane through the axis, for the turn-on movie
        sim.init_sim()
        vol = mp.Volume(**box.side_plane("y"))
        steps_per_period = round(1 / dt)
        every = frame_every or max(1, round(0.34 * steps_per_period))
        n_frames = int((T_ss + P.DFT_PERIODS) / (every * dt)) + 2
        probe = sim.get_array(component=mp.Ex, vol=vol)
        shape = (n_frames,) + probe.shape
        mm = (np.lib.format.open_memmap(os.path.join(out, "frames_xz_Ex.npy"), mode="w+",
                                        dtype=np.float16, shape=shape) if MASTER else None)
        times = []

        def grab(sim):
            if sim.timestep() % every == 0 and len(times) < n_frames:
                a = sim.get_array(component=mp.Ex, vol=vol)
                if MASTER:
                    mm[len(times)] = a.astype(np.float16)
                times.append(sim.meep_time())

        step_funcs.append(grab)
        meta.update(frame_every_steps=every, frame_dt=every * dt)
        eps = sim.get_array(component=mp.Dielectric, vol=vol)
        x, _, z, _ = sim.get_array_metadata(vol=vol)
        if MASTER:
            np.savez(os.path.join(out, "eps_xz.npz"), eps=eps, x=np.asarray(x),
                     z=np.asarray(z) + box.zc)

    sim.run(*step_funcs, until=T_ss)
    log(f"[{kind}] steady at t={sim.meep_time():.2f}, wall {time.time() - t_start:.0f} s")

    # ---- steady-state monitors, accumulated over whole periods -------------------------
    dfa = dict(decimation_factor=1)
    refl = sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.xy_plane(P.Z_REFL)), **dfa) \
        if kind != "mode" else None
    if kind == "air":
        facet = sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False, **box.xy_plane(0.0), **dfa)
    if kind == "mode":
        zm = box.snap(0.5 * (box.z_in_lo + box.z_in_hi) + 2.0)
        mode_dft = sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False, **box.xy_plane(zm), **dfa)
        mode_flux = sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.xy_plane(zm)), **dfa)
        mode_mon = sim.add_mode_monitor(P.FREQ, 0, 1, mp.ModeRegion(**box.mode_window(zm)), **dfa)
    mons, fluxes = {}, {}
    if kind == "fiber":
        for zp in (mode_planes or []):
            mons[zp] = sim.add_mode_monitor(P.FREQ, 0, 1, mp.ModeRegion(**box.mode_window(zp)), **dfa)
            fluxes[zp] = sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.xy_plane(zp)), **dfa)
        sides = {k: sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False,
                                       **box.side_plane(k), **dfa) for k in ("y", "x")}
        planes = []
        if xy_planes:
            dz = P.XY_PLANE_DZ
            zs = np.arange(box.z_in_lo + 1.0 + 0.5 * dz, box.z_in_hi, dz)
            zs = np.array([box.snap(z) for z in zs[zs > P.Z_SOURCE + 0.2]])
            for zp in zs:
                planes.append((zp, sim.add_dft_fields(COMPS, P.FREQ, 0, 1, yee_grid=False,
                                                      **box.xy_plane(zp), **dfa)))
    sim.run(*step_funcs, until=P.DFT_PERIODS)
    log(f"[{kind}] DFT done at t={sim.meep_time():.2f}, wall {time.time() - t_start:.0f} s")

    # ---- results ---------------------------------------------------------------------
    if refl is not None:
        meta["P_net_refl_plane"] = mp.get_fluxes(refl)[0]
    if kind == "air":
        save_dft_plane(sim, facet, os.path.join(out, "incident_facet.npz"), box.zc)
    if kind == "mode":
        meta["P_mode"] = mp.get_fluxes(mode_flux)[0]
        meta["z_mode_plane"] = zm
        r = sim.get_eigenmode_coefficients(mode_mon, [1], eig_parity=mp.EVEN_Y)
        meta["mpb_check"] = dict(P_fwd=abs(r.alpha[0, 0, 0]) ** 2, P_bwd=abs(r.alpha[0, 0, 1]) ** 2,
                                 neff=r.kpoints[0].z / P.FREQ)
        log(f"[mode] flux {meta['P_mode']:.5f}  MPB |a+|^2 {meta['mpb_check']['P_fwd']:.5f}"
            f"  |a-|^2 {meta['mpb_check']['P_bwd']:.2e}")
        save_dft_plane(sim, mode_dft, os.path.join(out, "mode_xy.npz"), box.zc)
    if kind == "fiber":
        kguess = fibermode.he11_neff() * P.FREQ
        rows = []
        for zp, mon in mons.items():
            r = sim.get_eigenmode_coefficients(
                mon, [1], eig_parity=mp.EVEN_Y,
                kpoint_func=lambda f, n: mp.Vector3(0, 0, kguess))
            rows.append(dict(z=zp, flux=mp.get_fluxes(fluxes[zp])[0], flux_window=mp.get_fluxes(mon)[0],
                             a_fwd=[r.alpha[0, 0, 0].real, r.alpha[0, 0, 0].imag],
                             a_bwd=[r.alpha[0, 0, 1].real, r.alpha[0, 0, 1].imag],
                             P_fwd=abs(r.alpha[0, 0, 0]) ** 2, P_bwd=abs(r.alpha[0, 0, 1]) ** 2,
                             neff=r.kpoints[0].z / P.FREQ, vgrp=r.vgrp[0]))
            log(f"  z={zp:6.1f}  flux={rows[-1]['flux']:.5f}  |a+|^2={rows[-1]['P_fwd']:.5f}"
                f"  |a-|^2={rows[-1]['P_bwd']:.2e}  neff={rows[-1]['neff']:.6f}")
        meta["mode_planes"] = rows
        for k, d in sides.items():
            save_dft_plane(sim, d, os.path.join(out, f"dft_side_{'xz' if k == 'y' else 'yz'}.npz"), box.zc)
        if planes:
            # stream the stack plane by plane: get_dft_array is collective and big
            zs = np.array([zp for zp, _ in planes])
            probe = sim.get_dft_array(planes[0][1], mp.Ex, 0)
            mms = {}
            if MASTER:
                for n in COMP_NAMES:
                    mms[n] = np.lib.format.open_memmap(
                        os.path.join(out, f"xy_{n}.npy"), mode="w+", dtype=np.complex64,
                        shape=(len(zs),) + probe.shape)
            for i, (zp, d) in enumerate(planes):
                for c, n in zip(COMPS, COMP_NAMES):
                    a = sim.get_dft_array(d, c, 0)
                    if MASTER:
                        mms[n][i] = a
            x, y, _, _ = sim.get_array_metadata(dft_cell=planes[0][1])
            if MASTER:
                for m in mms.values():
                    m.flush()
                np.savez(os.path.join(out, "xy_coords.npz"), x=np.asarray(x), y=np.asarray(y), z=zs)

    if check_steady and kind == "fiber":
        # a second, later window: if the field was steady, nothing changes
        mons2 = {zp: sim.add_mode_monitor(P.FREQ, 0, 1, mp.ModeRegion(**box.mode_window(zp)), **dfa)
                 for zp in (mode_planes or [])}
        refl2 = sim.add_flux(P.FREQ, 0, 1, mp.FluxRegion(**box.xy_plane(P.Z_REFL)), **dfa)
        sim.run(until=P.DFT_PERIODS)
        meta["steady_check"] = dict(P_net_refl_plane=mp.get_fluxes(refl2)[0], planes=[
            dict(z=zp, flux_window=mp.get_fluxes(m)[0],
                 P_fwd=abs(sim.get_eigenmode_coefficients(m, [1], eig_parity=mp.EVEN_Y).alpha[0, 0, 0]) ** 2)
            for zp, m in mons2.items()])

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
    ap.add_argument("--kind", choices=["fiber", "air", "mode"], default="fiber")
    ap.add_argument("--case", choices=list(P.CASES), help="one of the named launch conditions")
    ap.add_argument("--w0", type=float, help="beam waist on the facet (overrides --case)")
    ap.add_argument("--res", type=int, default=P.RESOLUTION)
    ap.add_argument("--length", type=float, default=P.FIBER_LENGTH)
    ap.add_argument("--half-width", type=float, default=P.HALF_WIDTH)
    ap.add_argument("--mode-planes", type=float, nargs="*", default=None)
    ap.add_argument("--movie", action="store_true", help="record Ex(x, z, t) frames")
    ap.add_argument("--xy-planes", action="store_true", help="record the stack of transverse planes")
    ap.add_argument("--check-steady", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    w0 = a.w0 if a.w0 is not None else (P.CASES[a.case]["w0"] if a.case else None)
    planes = a.mode_planes if a.mode_planes is not None else [z for z in P.MODE_PLANES if z < a.length - 1]
    run(a.kind, w0, a.res, a.length, a.half_width, a.out, movie=a.movie, xy_planes=a.xy_planes,
        mode_planes=planes, check_steady=a.check_steady)


if __name__ == "__main__":
    main()
