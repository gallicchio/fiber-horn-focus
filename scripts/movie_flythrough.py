"""Fly down the fiber: the transverse cross-section, plane by plane, for each launch.
The guided-mode part glows blue and everything else orange, so the "rest" of the light
can be followed directly.

    python scripts/movie_flythrough.py [--root data/main/res16] [--out media/cross_sections.mp4]
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


class Planes:
    """Guided / rest intensities of one launch's plane stack, computed on demand."""

    def __init__(self, run, xlim):
        self.run = run
        self.x, self.y, self.z, self.mm = run.planes()
        self.keep = abs(self.x) <= xlim
        self.A, self.beta, _ = run.guided_fit()
        self.pb = run.plane_budget()
        e = run.mode.f
        self.e = {k: e[k][np.ix_(self.keep, self.keep)] for k in ("Ex", "Ey", "Ez")}
        tot_max = 0.0
        for i in np.flatnonzero(self.z >= 3.0)[::4]:
            tot_max = max(tot_max, self.split(i)[0].max())
        self.ref = tot_max  # the brightest point 3 λ or more inside the glass; see share_ref

    def split(self, i):
        k = np.ix_(self.keep, self.keep)
        a = self.A * np.exp(1j * self.beta * self.z[i]) if self.z[i] > 0 else 0.0
        tot = np.zeros(k[0].shape[0:1] + k[1].shape[1:2])
        gui, rest = np.zeros_like(tot), np.zeros_like(tot)
        for c in ("Ex", "Ey", "Ez"):
            E = np.asarray(self.mm[c][i], dtype=np.complex128)[k]
            G = a * self.e[c]
            tot += abs(E) ** 2
            gui += abs(G) ** 2
            rest += abs(E - G) ** 2
        return tot, gui, rest


def share_ref(stacks):
    """Put every launch on one brightness scale, per unit of incident power (the scale of the
    brightest of them), so equal guided power looks equally bright."""
    per_watt = max(st.ref / st.run.P_inc for st in stacks.values())
    for st in stacks.values():
        st.ref = per_watt * st.run.P_inc
    return stacks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media", "cross_sections.mp4"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--xlim", type=float, default=12.0)
    ap.add_argument("--decades", type=float, default=3.0)
    a = ap.parse_args()

    W, H = a.width, int(a.width * 9 / 16)
    fs = 1.3 * W / 1920  # text stays legible in the half-width README preview
    t = style.use("dark", font_scale=fs)
    stacks = share_ref({c: Planes(load_case(a.root, c), a.xlim) for c in a.cases})
    z = stacks[a.cases[0]].z
    nz = len(z)

    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    fig.text(0.03, 0.95, "Flying down the fiber: the cross-section at each z",
             fontsize=19 * fs, fontweight="bold", va="center")
    hw0 = stacks[a.cases[0]].run.meta["half_width"]
    fig.text(0.03, 0.905, "Time-averaged |E|², log scale.  White: the beam in air.  Blue: light in the core's "
             f"guided mode.  Orange: everything else.  Window {2 * hw0:.0f} λ × {2 * hw0:.0f} λ.",
             fontsize=12.5 * fs, color=t["text2"], va="center")
    n = len(stacks)
    pw = 0.28
    gap = (1 - n * pw) / (n + 1)
    ims, labels = [], []
    for j, (name, st) in enumerate(stacks.items()):
        ax = fig.add_axes([gap + j * (pw + gap), 0.31, pw, pw * W / H])
        tot, g, r = st.split(0)
        im = ax.imshow(np.transpose(style.two_channel(g, r, st.ref, a.decades), (1, 0, 2)),
                       origin="lower", extent=[st.x[st.keep][0], st.x[st.keep][-1]] * 2,
                       interpolation="bilinear")
        circ = plt.Circle((0, 0), P.CORE_RADIUS, fill=False, color=t["text2"], lw=0.7, alpha=0.5)
        ax.add_patch(circ)
        ax.set_xticks([-10, 0, 10])
        ax.set_yticks([-10, 0, 10])
        ax.tick_params(labelsize=9 * fs)
        ax.set_xlabel("x / λ  (polarization →)", fontsize=10 * fs)
        if j == 0:
            ax.set_ylabel("y / λ", fontsize=10 * fs)
        c = P.CASES[name]
        ax.set_title(f"{c['title']}\n", fontsize=14 * fs, fontweight="bold", loc="left")
        ax.text(0, 1.025, c["subtitle"], transform=ax.transAxes, fontsize=11 * fs, color=t["text2"])
        lab = ax.text(0, -0.15, "", transform=ax.transAxes, fontsize=11 * fs, va="top", linespacing=1.5,
                      family="DejaVu Sans Mono")
        ims.append(im)
        labels.append(lab)
    # where we are along the fiber
    bar = fig.add_axes([0.06, 0.075, 0.88, 0.03])
    z0, z1 = z[0], z[-1]
    bar.add_patch(plt.Rectangle((0, -1), z1, 2, color=t["glass"], lw=0))
    bar.plot([0, z1], [0, 0], color=t["muted"], lw=3)
    bar.set_xlim(z0 - 0.5, z1 + 0.5)
    bar.set_ylim(-1.3, 1.3)
    bar.set_yticks([])
    for s in ("left", "right", "top"):
        bar.spines[s].set_visible(False)
    bar.tick_params(labelsize=9 * fs)
    bar.set_xlabel("z / λ   (shaded: the fiber, line: its core; left of 0: air)", fontsize=10 * fs, color=t["text2"])
    marker = bar.axvline(z0, color=t["text"], lw=2)
    ztxt = fig.text(0.97, 0.95, "", fontsize=15 * fs, ha="right", va="center")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    writer = FFMpegWriter(fps=a.fps, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    hold_start, hold_end = int(1.0 * a.fps), int(1.5 * a.fps)  # pause on the beam in air, and at the end
    with writer.saving(fig, a.out, dpi=100):
        for i in [0] * hold_start + list(range(nz)) + [nz - 1] * hold_end:
            for im, lab, st in zip(ims, labels, stacks.values()):
                tot, g, r = st.split(i)
                if z[i] > 0:
                    img = style.two_channel(g, r, st.ref, a.decades)
                else:  # in air nothing is guided yet: plain white
                    img = style.two_channel(0 * tot, 0 * tot, st.ref, a.decades, Ia=tot)
                im.set_data(np.transpose(img, (1, 0, 2)))
                if z[i] > 0:
                    gp = 100 * abs(st.A) ** 2 / st.run.P_inc
                    win = 100 * st.pb["total"][i]
                    out = max(0.0, 100 * (1 - st.run.R) - win)
                    hw = st.run.meta["half_width"]
                    lab.set_text(f"in the core mode            {gp:5.1f}%\n"
                                 f"uncoupled, in the window    {win - gp:5.1f}%\n"
                                 f"uncoupled, gone sideways    {out:5.1f}%")
                else:
                    lab.set_text(f"incident beam, in air\n(reflected at the facet: {100 * st.run.R:.1f}%)")
            marker.set_xdata([z[i], z[i]])
            ztxt.set_text(f"z = {z[i]:6.1f} λ")
            writer.grab_frame()
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
