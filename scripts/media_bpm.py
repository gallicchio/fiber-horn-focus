"""Figures and movie of the long-distance (BPM) continuation.

    python scripts/media_bpm.py [--data data/bpm] [--out media]
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

ROOT = os.path.join(os.path.dirname(__file__), "..")
LAM_UM = 1.55


def smf28_cm(z, d):
    """Equivalent SMF-28 length (cm) of z wavelengths of the scaled fiber at 1.55 um."""
    return z * LAM_UM * float(d["length_to_smf28"]) * 1e-4


def stacks(d):
    """Power fractions along z.  Light that has crossed the cladding surface into the
    (higher-index, lossy) coating counts as stripped, whether or not it is absorbed yet."""
    R = float(d["R"])
    guided = d["guided"]
    clad = d["in_glass"] - guided
    stripped = np.clip(1 - R - d["in_glass"], 0, None)
    return dict(guided=guided, cladding=clad, lost=stripped, reflected=np.full_like(guided, R))


def fig_budget(data, out, theme):
    t = style.use(theme)
    keys = ["guided", "cladding", "lost", "reflected"]
    labels = dict(style.LABELS, cladding="uncoupled, in the cladding", lost=style.LOST_BPM)
    fig, axs = plt.subplots(1, len(data), figsize=(4.3 * len(data), 4.0), sharey=True)
    for ax, (name, d) in zip(np.atleast_1d(axs), data.items()):
        z = d["z"]
        s = stacks(d)
        ax.stackplot(z, *[100 * s[k] for k in keys], colors=[t[k] for k in keys],
                     edgecolor=t["surface"], linewidth=1.5)
        ax.set_xlim(0, z[-1])
        ax.set_ylim(0, 100)
        ax.set_title(P.CASES[name]["title"], loc="left")
        ax.set_xlabel("z / λ")
        sec = ax.secondary_xaxis("top", functions=(lambda v, d=d: smf28_cm(v, d),
                                                   lambda c, d=d: c / smf28_cm(1.0, d)))
        sec.set_xlabel("same physics in SMF-28 at 1550 nm (cm)", fontsize=9, color=t["text2"])
        sec.tick_params(labelsize=9)
        # direct labels on the right edge
        acc = 0.0
        for k in keys:
            v = 100 * s[k][-1]
            if v > 4:
                ax.text(z[-1] * 0.97, acc + v / 2, f"{v:.0f}%", ha="right", va="center", fontsize=10,
                        color="#ffffff", fontweight="bold")
            acc += v
    np.atleast_1d(axs)[0].set_ylabel("percent of incident power")
    handles = [plt.Rectangle((0, 0), 1, 1, color=t[k]) for k in keys]
    ncol = 4 if len(data) >= 3 else 2  # one row when the figure is wide enough
    fig.legend(handles, [labels[k] for k in keys], ncol=ncol, loc="lower center", fontsize=10, handlelength=1.2)
    fig.tight_layout(rect=(0, 0.07 * 4 / ncol, 1, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_side(data, out):
    """The long x-z cuts, one row per launch: guided blue, the rest orange, cladding surface marked."""
    t = style.use("dark")
    fig, axs = plt.subplots(len(data), 1, figsize=(13, 2.3 * len(data) + 0.9), sharex=True)
    # one brightness scale for every row (the BPM fields are per unit of incident power)
    ref = max(max(d["side_guided"].max(), d["side_rest"][d["z"] > 30].max()) for d in data.values())
    for ax, (name, d) in zip(np.atleast_1d(axs), data.items()):
        z, x = d["z"], d["x"]
        keep = abs(x) <= float(d["r_coat"]) + 2
        g, r = d["side_guided"][:, keep].T, d["side_rest"][:, keep].T
        img = style.two_channel(g, r, ref, 3.0)
        ax.imshow(img, origin="lower", aspect="auto", extent=[z[0], z[-1], x[keep][0], x[keep][-1]],
                  interpolation="bilinear")
        for s in (-1, 1):
            ax.axhline(s * float(d["r_clad"]), color=t["text2"], lw=0.8, alpha=0.7)
        ax.text(z[-1], float(d["r_clad"]) + 0.8, "cladding surface → coating", ha="right", va="bottom",
                color=t["text2"], fontsize=9)
        c = P.CASES[name]
        ax.set_title(f"{c['title']}  ·  {c['subtitle']}", loc="left", fontsize=12)
        ax.set_ylabel("x / λ")
    np.atleast_1d(axs)[-1].set_xlabel("z / λ")
    fig.text(0.01, 0.005, "Long-distance (BPM) continuation of the 3D FDTD fields.  Log intensity, 3 decades, one scale for all rows.  "
             "Blue: the core mode.  Orange: everything else.", color=t["text2"], fontsize=10)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def movie(data, out, fps=24, width=1920):
    W, H = width, int(width * 9 / 16)
    fs = 1.3 * W / 1920  # text stays legible in the half-width README preview
    t = style.use("dark", font_scale=fs)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    fig.text(0.03, 0.95, "Farther down the fiber: where the uncoupled light ends up",
             fontsize=19 * fs, fontweight="bold", va="center")
    fig.text(0.03, 0.905, "Whole fiber cross-section (glass cladding inside the circle, coating outside).  "
             "Blue: core mode.  Orange: everything else.  Log intensity.", fontsize=12.5 * fs,
             color=t["text2"], va="center")
    n = len(data)
    pw = 0.28
    gap = (1 - n * pw) / (n + 1)
    ims, labs = [], []
    names = list(data)
    ref = max(max(d["snap_guided"].max(), d["snap_rest"][1:].max()) for d in data.values())  # shared by all
    for j, name in enumerate(names):
        d = data[name]
        ax = fig.add_axes([gap + j * (pw + gap), 0.27, pw, pw * W / H])
        xs = d["snap_x"]
        im = ax.imshow(np.transpose(style.two_channel(d["snap_guided"][0], d["snap_rest"][0], ref, 3.5), (1, 0, 2)),
                       origin="lower", extent=[xs[0], xs[-1], xs[0], xs[-1]], interpolation="bilinear")
        ax.add_patch(plt.Circle((0, 0), float(d["r_clad"]), fill=False, color=t["text2"], lw=0.8, alpha=0.8))
        ax.add_patch(plt.Circle((0, 0), float(d["r_coat"]), fill=False, color=t["muted"], lw=0.6, alpha=0.6, ls=":"))
        lim = float(d["r_coat"]) + 3
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([-30, 0, 30])
        ax.set_yticks([-30, 0, 30])
        ax.tick_params(labelsize=9 * fs)
        ax.set_xlabel("x / λ", fontsize=10 * fs)
        c = P.CASES[name]
        ax.set_title(f"{c['title']}\n", fontsize=14 * fs, fontweight="bold", loc="left")
        ax.text(0, 1.025, c["subtitle"], transform=ax.transAxes, fontsize=11 * fs, color=t["text2"])
        lab = ax.text(0, -0.13, "", transform=ax.transAxes, fontsize=11 * fs, va="top", linespacing=1.5,
                      family="DejaVu Sans Mono")
        ims.append(im)
        labs.append(lab)
    ztxt = fig.text(0.97, 0.95, "", fontsize=15 * fs, ha="right", va="center")
    z2txt = fig.text(0.97, 0.86, "", fontsize=11 * fs, ha="right", va="center", color=t["text2"])  # below the subtitle
    nsnap = min(len(data[k]["snap_z"]) for k in names)
    writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with writer.saving(fig, out, dpi=100):
        for i in list(range(nsnap)) + [nsnap - 1] * int(1.5 * fps):
            for im, lab, name in zip(ims, labs, names):
                d = data[name]
                im.set_data(np.transpose(style.two_channel(d["snap_guided"][i], d["snap_rest"][i], ref, 3.5), (1, 0, 2)))
                k = int(np.argmin(abs(d["z"] - d["snap_z"][i])))
                s = stacks(d)
                lab.set_text(f"in the core mode          {100 * s['guided'][k]:5.1f}%\n"
                             f"uncoupled, in cladding    {100 * s['cladding'][k]:5.1f}%\n"
                             f"stripped by the coating   {100 * s['lost'][k]:5.1f}%")
            zz = data[names[0]]["snap_z"][i]
            ztxt.set_text(f"z = {zz:6.0f} λ")
            z2txt.set_text(f"≈ {smf28_cm(zz, data[names[0]]):.1f} cm of SMF-28 at 1550 nm")
            writer.grab_frame()
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(ROOT, "data", "bpm"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    a = ap.parse_args()
    data = {c: np.load(os.path.join(a.data, f"{c}.npz")) for c in a.cases
            if os.path.exists(os.path.join(a.data, f"{c}.npz"))}
    for theme in ("light", "dark"):
        fig_budget(data, os.path.join(a.out, f"long_distance_budget_{theme}.png"), theme)
    fig_side(data, os.path.join(a.out, "long_distance_side.png"))
    movie(data, os.path.join(a.out, "long_distance.mp4"))


if __name__ == "__main__":
    main()
