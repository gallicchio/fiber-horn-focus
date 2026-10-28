"""The headline movie: three launches side by side as the light switches on, fills the
fiber and settles, then a cross-fade to the steady-state answer (guided part blue,
everything else orange) with the power budget.

    python scripts/movie_turnon.py [--root data/main/res16] [--out media/launch_movie.mp4]
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
from fig_results import budget, common_ref, side_image  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
BUDGET_LABELS = dict(guided="in the core mode", cladding="uncoupled, near core",
                     lost="uncoupled, gone sideways", reflected="reflected")


def field_rgb(frame, emax, xkeep):
    v = style.signed_compress(frame[xkeep].astype(np.float32), emax, 0.5)
    return style.field_cmap()(0.5 + 0.5 * v)[..., :3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media", "launch_movie.mp4"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--stride", type=int, default=1, help="use every n-th time frame")
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--xlim", type=float, default=12.0)
    ap.add_argument("--hold", type=float, default=4.0, help="seconds on the steady-state picture")
    a = ap.parse_args()

    fs = 1.3 * a.width / 1920  # text stays legible in the half-width README preview
    t = style.use("dark", font_scale=fs)
    runs = {c: load_case(a.root, c) for c in a.cases}
    refs = common_ref(runs, a.xlim)
    data = {}
    for name, run in runs.items():
        frames, times = run.frames()
        n = len(times)
        side = run.side("xz")
        x, z = side["x"], side["z"]
        xkeep = abs(x) <= a.xlim
        # field scale: the strongest steady field in the fiber, so each row is fully used
        last = np.asarray(frames[n - 8:n], dtype=np.float32)
        emax = np.percentile(abs(last[:, xkeep][:, :, z > 1]), 99.9)
        xs, zs, comp = side_image(run, a.xlim, ref=refs[name])
        data[name] = dict(frames=frames, times=times[:n], x=x[xkeep], z=z, xkeep=xkeep, emax=emax,
                          comp=comp, budget=budget(run))
    n_frames = min(len(d["times"]) for d in data.values())
    idx = list(range(0, n_frames, a.stride))
    n_hold = int(a.hold * a.fps)
    n_fade = int(1.2 * a.fps)

    W, H = a.width, int(a.width * 9 / 16)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    nrow = len(runs)
    top, bottom, left, right = 0.835, 0.085, 0.045, 0.765
    gap = 0.055
    h = (top - bottom - gap * (nrow - 1)) / nrow
    axes, ims, side_txt = [], [], []
    fig.text(left, 0.952, "Focusing a Gaussian beam into a single-mode fiber (3D FDTD)",
             fontsize=19 * fs, color=t["text"], fontweight="bold", va="center")
    sub = fig.text(left, 0.905, "", fontsize=13 * fs, color=t["text2"], va="center")
    clock = fig.text(0.985, 0.952, "", fontsize=13 * fs, color=t["text2"], ha="right", va="center")
    for i, (name, d) in enumerate(data.items()):
        y0 = top - (i + 1) * h - i * gap
        ax = fig.add_axes([left, y0, right - left, h])
        rgb = field_rgb(d["frames"][0], d["emax"], d["xkeep"])
        im = ax.imshow(rgb, origin="lower", aspect="auto", interpolation="antialiased",
                       extent=[d["z"][0], d["z"][-1], d["x"][0], d["x"][-1]])
        ax.plot([0, 0], [d["x"][0], d["x"][-1]], color=t["text2"], lw=0.8, alpha=0.7)
        for s in (-1, 1):
            ax.plot([0, d["z"][-1]], [s * P.CORE_RADIUS] * 2, color=t["text2"], lw=0.6, alpha=0.45)
        ax.set_xlim(d["z"][0], d["z"][-1])
        ax.set_ylim(d["x"][0], d["x"][-1])
        ax.set_yticks([-10, 0, 10])
        ax.set_ylabel("x / λ", fontsize=10 * fs)
        ax.tick_params(labelsize=9 * fs)
        if i < nrow - 1:
            ax.set_xticklabels([])
        c = P.CASES[name]
        ax.text(0.0, 1.03, c["title"], transform=ax.transAxes, fontsize=14 * fs,
                color=t["text"], fontweight="bold", va="bottom")
        ax.text(0.16, 1.03, c["subtitle"], transform=ax.transAxes, fontsize=12 * fs,
                color=t["text2"], va="bottom")
        axes.append(ax)
        ims.append(im)
        # the budget column
        b = d["budget"]
        ys = np.linspace(y0 + h * 0.85, y0 + h * 0.12, 4)
        items = []
        for yk, k in zip(ys, ["guided", "cladding", "lost", "reflected"]):
            sq = fig.text(right + 0.015, yk, "■", color=t[k], fontsize=15 * fs, va="center")
            tx = fig.text(right + 0.035, yk, f"{100 * b[k]:4.1f}%  {BUDGET_LABELS[k]}",
                          color=t["text"], fontsize=11.5 * fs, va="center")
            items += [sq, tx]
        side_txt.append(items)
    axes[-1].set_xlabel("z / λ     (fiber facet at z = 0; core outlined)", fontsize=10 * fs)
    z_end = data[a.cases[0]]["budget"]["z"]
    note = fig.text(right + 0.015, top + 0.005, f"steady state, {z_end:.0f} λ into the fiber", color=t["text2"],
                    fontsize=10.5 * fs, va="bottom")
    foot = fig.text(left, 0.018, "", fontsize=10.5 * fs, color=t["muted"])
    for items in side_txt:
        for it in items:
            it.set_alpha(0.0)
    note.set_alpha(0.0)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    writer = FFMpegWriter(fps=a.fps, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    with writer.saving(fig, a.out, dpi=100):
        for k in idx:
            tt = data[a.cases[0]]["times"][k]
            for im, d in zip(ims, data.values()):
                im.set_data(field_rgb(d["frames"][k], d["emax"], d["xkeep"]))
            clock.set_text(f"t = {tt:5.1f} periods")
            sub.set_text("Instantaneous electric field Eₓ in the plane of polarization, as the source "
                         "switches on (red / blue = + / − field)")
            foot.set_text("Meep 1.35, 16 px/λ, fiber V = 2.36 (like SMF-28), NA 0.25, core radius 1.5 λ. "
                          "Each row scaled to its own strongest steady field.")
            writer.grab_frame()
        # cross-fade to the steady-state split
        last = [field_rgb(d["frames"][idx[-1]], d["emax"], d["xkeep"]) for d in data.values()]
        sub.set_text("Steady state, time-averaged |E|² (log, 3 decades): blue = in the core's guided mode, "
                     "orange = everything else, white = still in air")
        foot.set_text("Meep 1.35, 16 px/λ, fiber V = 2.36 (like SMF-28), NA 0.25, core radius 1.5 λ. "
                      "One brightness scale for all rows (per watt of incident light).")
        for j in range(n_fade + n_hold):
            f = min(1.0, j / n_fade)
            f = f * f * (3 - 2 * f)
            for im, l0, d in zip(ims, last, data.values()):
                im.set_data((1 - f) * l0 + f * d["comp"])
            for items in side_txt:
                for it in items:
                    it.set_alpha(f)
            note.set_alpha(f)
            writer.grab_frame()
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
