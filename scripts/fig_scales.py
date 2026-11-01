"""How big each simulation is: a zoom-out series for the too-big spot (launch NA too narrow).

  1. the facet, the first few wavelengths of the 3D FDTD (the field itself, true aspect);
  2. the whole 3D FDTD region (true aspect), inside the real cladding it does not reach;
  3. the long-distance BPM: the whole cladding and coating, 6000 wavelengths;
  4. a 1 m patch cord, for scale.

Each panel is drawn as a box inside the next one, joined by zoom lines.

    python scripts/fig_scales.py [--root data/main/res16] [--bpm data/bpm] [--out media/scales.png]
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import ConnectionPatch, Rectangle  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from fiberlaunch import bpm as B  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402
from fig_results import side_image  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
LAM_UM = 1.55
MM_PER_LAMBDA = LAM_UM * 1e-3 * B.LENGTH_TO_SMF28  # along the fiber, in SMF-28
UM_PER_LAMBDA_ACROSS = LAM_UM * B.SCALE  # across the fiber, in SMF-28

FACET = (-8.0, 16.0, -6.0, 6.0)  # panel 1: z0, z1, x0, x1 (wavelengths); the source is 8 λ out
FDTD_X = 23.6  # panel 2 shows |x| up to here, just past the real cladding surface
PATCH_CM = 100.0


class Layout:
    """Axes and text placed in inches from the top of a fixed-size figure."""

    def __init__(self, fig, t):
        self.fig, self.t = fig, t
        self.W, self.H = fig.get_size_inches()

    def axes(self, left, top, w, h):
        return self.fig.add_axes([left / self.W, 1 - (top + h) / self.H, w / self.W, h / self.H])

    def text(self, left, top, s, **kw):
        """Text whose bottom sits `top` inches below the figure top, on an opaque patch
        so that zoom lines pass behind it."""
        kw = dict(dict(va="bottom", ha="left", zorder=5, color=self.t["text2"], fontsize=10.5,
                       bbox=dict(boxstyle="square,pad=0.12", fc=self.t["surface"], ec="none")), **kw)
        return self.fig.text(left / self.W, 1 - top / self.H, s, **kw)

    def heading(self, left, panel_top, head, sub, top_axis):
        """Bold title and a one-line subtitle above a panel (above its top axis, if any)."""
        y = panel_top - (0.40 if top_axis else 0.08)
        if sub:
            self.text(left, y, sub)
            y -= 0.27
        self.text(left, y, head, color=self.t["text"], fontsize=12.5, fontweight="bold")


def zoom(fig, small, big, box, t):
    """Mark `box` (z0, z1, x0, x1) on `big` and join its top corners to the bottom
    corners of `small`, the panel that shows it enlarged."""
    z0, z1, x0, x1 = box
    big.add_patch(Rectangle((z0, x0), z1 - z0, x1 - x0, fill=False, ec=t["text"], lw=1.2, zorder=5))
    for ax_frac, zz in ((0.0, z0), (1.0, z1)):
        fig.add_artist(ConnectionPatch(xyA=(ax_frac, 0.0), coordsA=small.transAxes, xyB=(zz, x1),
                                       coordsB=big.transData, color=t["text2"], lw=0.9, alpha=0.8,
                                       ls=(0, (4, 3)), zorder=1))


def top_axis(ax, scale, unit, t):
    """A secondary z axis in SMF-28 length, units on every tick (no separate label line)."""
    sec = ax.secondary_xaxis("top", functions=(lambda v: v * scale, lambda v: v / scale))
    sec.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g} {unit}"))
    sec.tick_params(labelsize=9.5, colors=t["text2"])
    return sec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--bpm", default=os.path.join(ROOT, "data", "bpm"))
    ap.add_argument("--case", default="narrow")
    ap.add_argument("--out", default=os.path.join(ROOT, "media", "scales.png"))
    a = ap.parse_args()
    run = load_case(a.root, a.case)
    d = np.load(os.path.join(a.bpm, f"{a.case}.npz"))
    t = style.use("dark")

    W, left, pw = 13.0, 0.9, 11.6  # figure width, left margin, width of the long panels
    z0, z1, x0, x1 = FACET
    w1 = 4.4
    h1 = w1 * (x1 - x0) / (z1 - z0)  # panels 1 and 2 at true aspect
    L_fdtd = run.meta["fiber_length"]
    fz0, fz1 = run.meta["z_source"] - 1.0, L_fdtd
    h2 = pw * 2 * FDTD_X / (fz1 - fz0)
    h3, h4 = 1.9, 0.32
    gap_title, gap_axis = 0.45, 0.45  # a bottom axis label, then a heading (+ a top axis)
    tops = [0.45]
    tops.append(tops[0] + h1 + gap_axis + 1.0)
    tops.append(tops[1] + h2 + gap_axis + 1.05)
    tops.append(tops[2] + h3 + gap_axis + 0.75)
    H = tops[3] + h4 + gap_title + 0.55
    fig = plt.figure(figsize=(W, H))
    L = Layout(fig, t)
    ax1 = L.axes(left, tops[0], w1, h1)
    ax2 = L.axes(left, tops[1], pw, h2)
    ax3 = L.axes(left, tops[2], pw, h3)
    ax4 = L.axes(left, tops[3], pw, h4)
    for ax in (ax1, ax2, ax3):
        ax.set_facecolor("#000000")

    # 1. the facet: the field itself, at one instant (the steady state's real part)
    s = run.side("xz")
    x, z = s["x"], s["z"]
    kx, kz = (x >= x0 - 0.1) & (x <= x1 + 0.1), (z >= z0 - 0.1) & (z <= z1 + 0.1)
    E = s["Ex"][np.ix_(kx, kz)].astype(complex)
    i0, j0 = np.argmin(abs(x[kx])), np.argmin(abs(z[kz]))
    E = np.real(E * np.exp(-1j * np.angle(E[i0, j0])))
    img = style.field_cmap()(0.5 + 0.5 * style.signed_compress(E, abs(E).max(), 0.5))
    ax1.imshow(img, origin="lower", aspect="equal", interpolation="bilinear",
               extent=[z[kz][0], z[kz][-1], x[kx][0], x[kx][-1]])
    ax1.axvline(0, color=t["text2"], lw=0.8)
    for sgn in (-1, 1):
        ax1.plot([0, z1], [sgn * P.CORE_RADIUS] * 2, color=t["text2"], lw=0.8, alpha=0.7)
    for zz, ha, s_ in ((-0.3, "right", "air"), (0.3, "left", "glass")):
        ax1.text(zz, x1 - 0.4, s_, color=t["text"], fontsize=10, ha=ha, va="top",
                 bbox=dict(boxstyle="square,pad=0.1", fc="#000000", ec="none", alpha=0.6))
    ax1.text(z1 - 0.3, P.CORE_RADIUS + 0.25, "core", color=t["text"], fontsize=9.5, ha="right", va="bottom",
             bbox=dict(boxstyle="square,pad=0.1", fc="#000000", ec="none", alpha=0.6))
    ax1.set_xlim(z0, z1)
    ax1.set_ylim(x0, x1)
    ax1.set_xlabel("z / λ")
    ax1.set_ylabel("x / λ")
    L.heading(left, tops[0], "①  The facet: the last 8 λ of the beam in air and the first 16 λ of fiber", None, False)
    L.text(left + w1 + 0.5, tops[0] + 0.05,
           "The electric field itself at one instant (red / blue = + / −),\n"
           "for the spot that is too big.  Each red-and-blue pair of stripes\n"
           "is one wavelength: λ in air, λ / 1.46 in the glass.\n\n"
           "Drawn to scale (true aspect ratio).", va="top", linespacing=1.45, bbox=None)

    # 2. the whole 3D FDTD region
    xs, zs, img2 = side_image(run, P.HALF_WIDTH)
    ax2.imshow(img2, origin="lower", aspect="auto", interpolation="bilinear",
               extent=[zs[0], zs[-1], xs[0], xs[-1]])
    hw = run.meta["half_width"]
    for sgn in (-1, 1):  # beyond the window: glass that the 3D simulation does not contain
        ax2.add_patch(Rectangle((fz0, min(sgn * hw, sgn * FDTD_X)), fz1 - fz0, FDTD_X - hw,
                                facecolor=t["glass"], edgecolor=t["muted"], hatch="////", lw=0, alpha=0.9))
        ax2.axhline(sgn * B.R_CLAD, color=t["text"], lw=1.0)
        ax2.plot([0, fz1], [sgn * P.CORE_RADIUS] * 2, color=t["muted"], lw=0.7, alpha=0.7)
    ax2.axvline(0, color=t["text2"], lw=0.8)
    box = dict(boxstyle="square,pad=0.15", fc=t["glass"], ec="none")
    ax2.text(fz1 - 1, B.R_CLAD - 0.7, "real cladding surface (62.5 µm from the axis in SMF-28).  Hatched: glass that "
             "the 3D simulation leaves out", color=t["text"], fontsize=9.5, ha="right", va="top", bbox=box)
    ax2.text(fz1 - 1, hw + 0.6, "edge of the 3D window: light that reaches it is absorbed (\"gone sideways\")",
             color=t["text"], fontsize=9.5, ha="right", va="bottom", bbox=box)
    ax2.set_xlim(fz0, fz1)
    ax2.set_ylim(-FDTD_X, FDTD_X)
    ax2.set_xlabel("z / λ")
    ax2.set_ylabel("x / λ")
    top_axis(ax2, MM_PER_LAMBDA, "mm", t)
    win_um = 2 * hw * UM_PER_LAMBDA_ACROSS
    L.heading(left, tops[1],
              f"②  The whole 3D FDTD simulation: {L_fdtd:.0f} λ of fiber, a {2 * hw:.0f} λ × {2 * hw:.0f} λ "
              f"window  (≈ {L_fdtd * MM_PER_LAMBDA:.1f} mm of SMF-28, {win_um:.0f} µm across)",
              f"Zoomed out {(fz1 - fz0) / (z1 - z0):.0f}× from ①.  Time-averaged |E|², log scale.  Blue: the "
              "core mode.  Orange: the rest.  White: the beam in air.  Drawn to scale (true aspect ratio).", True)

    # 3. the long-distance model
    zb, xb = d["z"], d["x"]
    rc, rcl = float(d["r_coat"]), float(d["r_clad"])
    keep = abs(xb) <= rc + 1.5
    g, r = d["side_guided"][:, keep].T, d["side_rest"][:, keep].T
    ref = max(g.max(), r[:, zb > 30].max())
    ax3.imshow(style.two_channel(g, r, ref, 3.0), origin="lower", aspect="auto", interpolation="bilinear",
               extent=[zb[0], zb[-1], xb[keep][0], xb[keep][-1]])
    for sgn in (-1, 1):
        ax3.axhline(sgn * rcl, color=t["text"], lw=0.9)
        ax3.axhline(sgn * rc, color=t["muted"], lw=0.7, ls=(0, (2, 2)))
    L3 = zb[-1]
    ax3.text(L3 - 30, 0.5 * (rcl + rc), "coating: absorbs what crosses the cladding surface", color=t["text2"],
             fontsize=9.5, ha="right", va="center")
    ax3.set_xlim(0, L3)
    ax3.set_ylim(-rc - 1.5, rc + 1.5)
    ax3.set_xlabel("z / λ")
    ax3.set_ylabel("x / λ")
    top_axis(ax3, MM_PER_LAMBDA / 10, "cm", t)
    stretch = (h3 / (2 * rc + 3)) / (pw / L3)
    L.heading(left, tops[2],
              f"③  The long-distance model (BPM): {L3:.0f} λ, the whole cross-section with the coating  "
              f"(≈ {L3 * MM_PER_LAMBDA / 10:.1f} cm of SMF-28)",
              f"Zoomed out {L3 / (fz1 - fz0):.0f}× from ②, and started from the 3D field 5 λ inside the facet.  "
              f"Not to scale: stretched {stretch:.0f}× across.", True)

    # 4. a patch cord
    cm = L3 * MM_PER_LAMBDA / 10
    ax4.add_patch(Rectangle((0, -0.5), PATCH_CM, 1.0, facecolor=t["glass"], edgecolor=t["muted"], lw=0.8))
    ax4.set_xlim(0, PATCH_CM)
    ax4.set_ylim(-1, 1)
    ax4.set_yticks([])
    for sp in ("left", "right", "top"):
        ax4.spines[sp].set_visible(False)
    ax4.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g} cm"))
    L.heading(left, tops[3], "④  For scale: a 1 m patch cord of SMF-28",
              f"Zoomed out {PATCH_CM / cm:.0f}× from ③: the long-distance model covers its first {cm:.1f} cm, "
              f"the 3D simulation its first {L_fdtd * MM_PER_LAMBDA:.1f} mm.  Not to scale across "
              "(the fiber is 0.25 mm thick).", False)

    zoom(fig, ax1, ax2, FACET, t)
    zoom(fig, ax2, ax3, (fz0, fz1, -hw, hw), t)
    zoom(fig, ax3, ax4, (0, cm, -0.5, 0.5), t)
    L.text(0.15, H - 0.12, "The spot that is too big (3× the mode, beam NA 0.064: too narrow).  The simulated fiber is SMF-28 "
           "scaled to a 1.5 λ core radius, which keeps the physics: 1 λ in the pictures is "
           f"{UM_PER_LAMBDA_ACROSS:.2f} µm\nof SMF-28 at 1550 nm across the fiber and {MM_PER_LAMBDA * 1e3:.2f} µm "
           "along it (see How it was computed).", color=t["muted"], fontsize=9, va="bottom", bbox=None)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=150, facecolor=t["surface"])
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
