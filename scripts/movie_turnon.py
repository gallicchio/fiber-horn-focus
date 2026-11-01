"""The headline movie: the three launches side by side as the light switches on, fills the
fiber and settles, then a cross-fade to the steady-state answer (guided part blue, everything
else orange) with the power budget.  As in the horn movies, the light comes down from the top.
It shows the first 40 λ of the fiber, where the light sorts itself out, at its true shape (x and
z to one scale), until that part is steady (100 periods), slowed to 2 periods a second; the
frames between the saved ones come from fiberlaunch.frames.  The frame is as tall as that needs:
1920 × 1616 for the defaults.

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
from fiberlaunch import frames  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402
from fig_results import budget, common_ref, side_image  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
BUDGET_LABELS = dict(guided="in the core mode", cladding="uncoupled, near the core",
                     lost="uncoupled, gone sideways", reflected="reflected")


def field_rgb(frame, emax):
    """An (x, z) field as a (z, x) picture."""
    v = style.signed_compress(frame.T, emax, 0.5)
    return style.field_cmap()(0.5 + 0.5 * v)[..., :3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media", "launch_movie.mp4"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--per-second", type=float, default=2.0, help="periods shown each second")
    ap.add_argument("--t-max", type=float, default=100.0, help="periods; up to --zmax it is steady by then")
    ap.add_argument("--zmax", type=float, default=40.0)
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
        saved, times = run.frames()
        side = run.side("xz")
        x, z = side["x"], side["z"]
        xkeep, zkeep = abs(x) <= a.xlim, (z >= run.meta["z_source"]) & (z <= a.zmax)  # from the source on
        n = int(np.searchsorted(times, a.t_max)) + 10  # and the frames the fits near t_max need
        # field scale: the strongest steady field in the fiber, so each column is fully used
        last = np.asarray(saved[len(times) - 8:len(times)], np.float32)[:, xkeep][:, :, zkeep]
        emax = np.percentile(abs(last[:, :, z[zkeep] > 1]), 99.9)
        _, _, comp = side_image(run, a.xlim, ref=refs[name])
        data[name] = dict(saved=np.asarray(saved[:n], np.float32)[:, xkeep][:, :, zkeep], times=times[:n],
                          x=x[xkeep], z=z[zkeep], emax=emax, comp=comp[:, zkeep].transpose(1, 0, 2),
                          budget=budget(run))
    t_out = np.arange(0, a.t_max, a.per_second / a.fps)
    n_hold = int(a.hold * a.fps)
    n_fade = int(1.2 * a.fps)

    # the layout in pixels: three columns with one scale for x and z, and the frame as tall as
    # that makes them
    W, k = a.width, a.width / 1920
    pad_l, pad_r, gap = 115 * k, 30 * k, 60 * k
    above, below = 230 * k, 290 * k  # the titles; the ticks, the power budget and the footer
    d0 = data[a.cases[0]]
    col_w = (W - pad_l - pad_r - 2 * gap) / len(runs)
    col_h = col_w * (d0["z"][-1] - d0["z"][0]) / (d0["x"][-1] - d0["x"][0])
    H = int(np.ceil((above + col_h + below) / 16) * 16)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    ims, faders = [], []
    fig.text(pad_l / W, 1 - 55 * k / H, "Focusing a Gaussian beam into a single-mode fiber (3D FDTD)",
             fontsize=19 * fs, color=t["text"], fontweight="bold", va="center")
    sub = fig.text(pad_l / W, 1 - 108 * k / H, "", fontsize=13 * fs, color=t["text2"], va="center")
    clock = fig.text(1 - pad_r / W, 1 - 55 * k / H, "", fontsize=13 * fs, color=t["text2"], ha="right",
                     va="center")
    for i, (name, d) in enumerate(data.items()):
        x0 = pad_l + i * (col_w + gap)
        ax = fig.add_axes([x0 / W, below / H, col_w / W, col_h / H])
        im = ax.imshow(field_rgb(d["saved"][0], d["emax"]), origin="lower", aspect="auto",
                       interpolation="antialiased", extent=[d["x"][0], d["x"][-1], d["z"][0], d["z"][-1]])
        ax.plot([d["x"][0], d["x"][-1]], [0, 0], color=t["text2"], lw=0.8, alpha=0.7)
        for s in (-1, 1):
            ax.plot([s * P.CORE_RADIUS] * 2, [0, d["z"][-1]], color=t["text2"], lw=0.6, alpha=0.45)
        ax.set_xlim(d["x"][0], d["x"][-1])
        ax.set_ylim(d["z"][-1], d["z"][0])
        ax.set_xticks([-10, 0, 10])
        ax.set_yticks(np.arange(0, a.zmax + 1, 10))
        ax.set_xlabel("x / λ", fontsize=10 * fs)
        ax.tick_params(labelsize=9 * fs)
        if i == 0:
            ax.set_ylabel("z / λ     (fiber facet at z = 0; core outlined)", fontsize=10 * fs)
        else:
            ax.set_yticklabels([])
        ims.append(im)
        c = P.CASES[name]
        y_top = (below + col_h) / H
        fig.text(x0 / W, y_top + 50 * k / H, c["title"], fontsize=14 * fs, color=t["text"], fontweight="bold",
                 va="bottom")
        fig.text(x0 / W, y_top + 14 * k / H, c["subtitle"], fontsize=12 * fs, color=t["text2"], va="bottom")
        # where the power goes, under the column, faded in at the end
        b = d["budget"]
        for j, key in enumerate(["guided", "cladding", "lost", "reflected"]):
            xk, yk = x0, below - (100 + 38 * j) * k
            faders += [fig.text(xk / W, yk / H, "■", color=t[key], fontsize=13 * fs, va="center"),
                       fig.text((xk + 26 * k) / W, yk / H, f"{100 * b[key]:.1f}% {BUDGET_LABELS[key]}",
                                color=t["text"], fontsize=11.5 * fs, va="center")]
    z_end = d0["budget"]["z"]
    foot = fig.text(pad_l / W, 22 * k / H, "", fontsize=10.5 * fs, color=t["muted"], va="bottom")
    for it in faders:
        it.set_alpha(0.0)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    writer = FFMpegWriter(fps=a.fps, codec="libx264",
                          extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    with writer.saving(fig, a.out, dpi=100):
        sub.set_text("Instantaneous electric field Eₓ in the plane of polarization, as the source "
                     "switches on (red / blue = + / − field)")
        foot.set_text("Meep 1.35, 16 px/λ, fiber V = 2.36 (like SMF-28), NA 0.25, core radius 1.5 λ. "
                      f"Each column scaled to its own strongest steady field. Shown at {a.per_second:g} "
                      "periods a second.")
        for tt in t_out:
            for im, d in zip(ims, data.values()):
                im.set_data(field_rgb(frames.at(d["saved"], d["times"], tt), d["emax"]))
            clock.set_text(f"t = {tt:5.1f} periods")
            writer.grab_frame()
        # cross-fade to the steady-state split
        last = [field_rgb(frames.at(d["saved"], d["times"], t_out[-1]), d["emax"]) for d in data.values()]
        sub.set_text("Steady state, time-averaged |E|² (log, 3 decades): blue = in the core's guided mode, "
                     "orange = everything else, white = still in air")
        foot.set_text("Meep 1.35, 16 px/λ, fiber V = 2.36 (like SMF-28), NA 0.25, core radius 1.5 λ. "
                      "One brightness scale for all three (per watt of incident light).")
        clock.set_text(f"where the power is, {z_end:.0f} λ into the fiber")
        for j in range(n_fade + n_hold):
            f = min(1.0, j / n_fade)
            f = f * f * (3 - 2 * f)
            for im, l0, d in zip(ims, last, data.values()):
                im.set_data((1 - f) * l0 + f * d["comp"])
            for it in faders:
                it.set_alpha(f)
            writer.grab_frame()
    plt.close(fig)
    print(f"wrote {a.out}: {W} x {H}, {col_w / (d0['x'][-1] - d0['x'][0]):.2f} px/λ across and "
          f"{col_h / (d0['z'][-1] - d0['z'][0]):.2f} along")


if __name__ == "__main__":
    main()
