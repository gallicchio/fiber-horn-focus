"""Figures and the movie for the metal waveguide (horn) runs, with the fiber runs alongside.

    python scripts/horn_media.py [--root data/horn/res120] [--fiber-root data/main/res16]
                                 [--only side modes budget movie]

media/horn_side_view.png              steady state on the plane of polarization, the three launches
media/horn_modes_{dark,light}.png     the power in each of the horn's own modes along it, in and out
media/fiber_vs_horn_{dark,light}.png  where the light goes: the fiber and the metal guide side by side
media/horn_movie.mp4                  the light switching on, then the steady-state picture

The horn pictures stand on end: the light comes in from the top, the mouth faces up.
"""

import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402
from fiberlaunch.horn_analysis import MODES, HornRun, cutoff_radius  # noqa: E402
from fig_results import budget as fiber_budget  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
METAL = "#4a4a46"
ZLIM = (-4.0, 18.4)  # from the source plane; mouth at 0, throat from L = 12.4 to the absorber at 18.4
SUBTITLES = dict(narrow="3× the best spot · beam NA 0.064",
                 matched="the best spot · beam NA 0.19",
                 wide="⅓ of the best spot · beam NA 0.58")


def mirror(a):
    """(z, r) array -> (z, x) image on the plane of polarization, x = -r ... r."""
    return np.concatenate([a[:, ::-1], a], axis=1)


def metal_mask(run, stride=1):
    """(z, r) pixels that are metal, or partly metal (Meep reports a perfect conductor as a large
    negative permittivity, and averages it with air on the pixels the wall cuts)."""
    eps = np.load(os.path.join(run.dir, "eps_rz.npz"))["eps"][::stride, ::stride]
    return (eps < 0) | (eps > 1.01)


def common_ref(runs, z_from=0.5):
    """One brightness scale per watt of incident light: the brightest point inside the metal
    over all launches."""
    per_watt = 0.0
    for run in runs.values():
        s = run.side_split()
        per_watt = max(per_watt, s["total"][s["z"] >= z_from].max() / run.P_inc)
    return {name: per_watt * run.P_inc for name, run in runs.items()}


def horn_image(run, xlim, ref, stride=1, decades=3.0):
    """Steady state on the plane of polarization, (z, x, 3): incident light white, reflected
    light and everything inside that is not TE11 heading in yellow, TE11 heading in blue;
    metal gray."""
    s = run.side_split()
    r = s["r"][::stride]
    keep = r <= xlim

    def part(a):
        return mirror(a[::stride, ::stride][:, keep])

    img = style.composite([(part(s["incident"]), style.GLOW["air"]),
                           (part(s["reflected"] + s["rest"]), style.GLOW["reflected"]),
                           (part(s["guided"]), style.GLOW["guided"])], ref, decades)
    img[mirror(metal_mask(run, stride)[:, keep])] = matplotlib.colors.to_rgb(METAL)
    return np.concatenate([-r[keep][::-1], r[keep]]), s["z"][::stride], img


def extent(x, z):
    dx, dz = x[1] - x[0], z[1] - z[0]
    return [x[0] - dx / 2, x[-1] + dx / 2, z[0] - dz / 2, z[-1] + dz / 2]


def pct(v):
    return f"{100 * v:.0f}%" if v >= 0.0995 else f"{100 * v:.1f}%"


def split(run):
    """The horn's budget, with the reflected light split by where it turned back."""
    b = run.budget()
    e = run.entry()
    back = min(e["back"], b["reflected"])  # the mode split is good to about half a percent
    return dict(guided=b["guided"], back=back, face=b["reflected"] - back,
                reflected=b["reflected"], entered=e["entered"], back_modes=e["back"], total=b["total"])


def label_horn(ax, t, run):
    kw = dict(color=t["text2"], fontsize=8.5, ha="left", va="center")
    ax.text(-5.6, 4.0, "metal", **kw)
    ax.text(run.mouth * 0.55 + 0.3, 0.5 * run.length, "horn", **kw)
    ax.text(P.THROAT_RADIUS + 0.35, run.length + 3.0, "single-mode\nguide", **kw)


def fig_side(runs, out, xlim=6.0, zlim=ZLIM):
    """The three launches side by side, the light coming down onto the mouth."""
    t = style.use("dark")
    refs = common_ref(runs)
    fig, axs = plt.subplots(1, len(runs), figsize=(12.0, 8.4), sharey=True)
    for i, (ax, (name, run)) in enumerate(zip(np.atleast_1d(axs), runs.items())):
        x, z, img = horn_image(run, xlim, refs[name])
        ax.imshow(img, origin="lower", aspect="equal", extent=extent(x, z), interpolation="bilinear")
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(zlim[1], zlim[0])
        b = split(run)
        ax.text(0, 1.105, P.CASES[name]["title"], transform=ax.transAxes, color=t["text"], fontsize=13,
                fontweight="bold", va="bottom")
        ax.text(0, 1.062, SUBTITLES[name], transform=ax.transAxes, color=t["text2"], fontsize=10, va="bottom")
        ax.text(0, 1.012, f"into the guide {pct(b['guided'])}  ·  reflected {pct(b['reflected'])}",
                transform=ax.transAxes, color=t["text"], fontsize=10.5, va="bottom")
        ax.set_xlabel("x / λ")
        ax.set_xticks([-5, 0, 5])
        ax.set_yticks([0, 5, 10, 15])
        if i == 0:
            ax.set_ylabel("z / λ   (mouth at z = 0; the light comes down from the top)")
            label_horn(ax, t, run)
    fig.text(0.01, 0.012, "Time-averaged |E|² on the plane of polarization, log scale (3 decades), one scale for all "
             "three (per watt in).  White: the incident beam.\nYellow: reflected light, and the light inside the "
             "horn that is not the guide's mode heading in (all of it turns back).  Blue: the guide's mode (TE11) "
             "heading in.  Gray: metal.", color=t["text2"], fontsize=9.5, va="bottom")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.88, bottom=0.12, wspace=0.08)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def mode_curves(run, near_cutoff=0.2):
    """TE11 heading in, the other modes heading in, and everything heading out, against z.
    A mode is left out within the stretch just before its cut-off where β/k < near_cutoff:
    there it turns around, and splitting the standing wave into in and out is meaningless."""
    z, mp_ = run.mode_powers()
    rad = run.radius(z)
    te11 = mp_["TE11"][0]
    other_in, out = np.zeros_like(z), np.zeros_like(z)
    for name, kind, root in MODES:
        b = np.sqrt(np.clip(1 - (root / (2 * np.pi * rad)) ** 2, 0, None))
        ok = (b >= near_cutoff) & np.isfinite(mp_[name][0])
        if name != "TE11":
            other_in += np.where(ok, mp_[name][0], 0.0)
        out += np.where(ok, mp_[name][1], 0.0)
    return z, te11, other_in, out


def cutoff_z(run, root):
    """Where the horn's radius equals the cut-off radius of a mode (None if it never does)."""
    rc = cutoff_radius(root)
    if not P.THROAT_RADIUS < rc < run.mouth:
        return None
    return run.length * (run.mouth - rc) / (run.mouth - P.THROAT_RADIUS)


def fig_modes(runs, out, theme):
    """Small multiples: inside the horn, the power heading in in TE11 and in the other modes,
    and the power heading back out, against depth; each mode's cut-off marked."""
    t = style.use(theme)
    fig, axs = plt.subplots(1, len(runs), figsize=(12.5, 4.1), sharey=True)
    for i, (ax, (name, run)) in enumerate(zip(np.atleast_1d(axs), runs.items())):
        z, te11, other, back = mode_curves(run)
        ax.axvspan(run.length, z[-1], color=t["grid"], alpha=0.6, linewidth=0)
        for mname, _, root in MODES:
            zc = cutoff_z(run, root)
            if zc is not None:
                ax.axvline(zc, color=t["muted"], lw=0.8, ls=":")
                ax.text(zc, 101, mname, rotation=90, ha="center", va="bottom", fontsize=7.5, color=t["text2"])
        ax.plot(z, 100 * other, color=t["reflected"], ls=(0, (4, 2)), label="other modes, heading in")
        ax.plot(z, 100 * back, color=t["reflected"], label="all modes, heading back out")
        ax.plot(z, 100 * te11, color=t["guided"], label="TE11, heading in")
        ax.text(z[-1] - 0.3, 100 * te11[-1] + 3, f"{100 * te11[-1]:.0f}%", ha="right", va="bottom",
                color=t["text"], fontsize=10)
        ax.text(run.length + 0.3, 92, "single-mode\nguide", fontsize=8.5, color=t["text2"], va="top")
        ax.set_title(P.CASES[name]["title"], loc="left", pad=34)
        ax.set_xlim(0, z[-1])
        ax.set_ylim(0, 100)
        ax.set_xlabel("depth z / λ (mouth at 0)")
        ax.grid(True, axis="y")
        if i == 0:
            ax.set_ylabel("percent of the incident power")
    np.atleast_1d(axs)[0].legend(loc="center", bbox_to_anchor=(0.45, 0.62), fontsize=9)
    fig.text(0.01, 0.01, "Dotted: where each higher mode is cut off (the horn's radius is too small for it); "
             "a mode's light turns back there.  The shaded stretch is the straight, single-mode guide.",
             color=t["text2"], fontsize=9, va="bottom")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_budget(runs, fibers, out, theme):
    """Stacked bars, the fiber above the metal guide for each launch: where the light goes."""
    t = style.use(theme)
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    on = t["surface"] if theme == "dark" else "#ffffff"
    hatch = dict(hatch="////", edgecolor=t["surface"], linewidth=0)
    rows, ticks = [], []
    y = 0.0
    for name in runs:
        fb, hb = fiber_budget(fibers[name]), split(runs[name])
        rows.append((y, [(fb["guided"], t["guided"], {}), (fb["cladding"], t["cladding"], {}),
                         (fb["lost"], t["lost"], {}), (fb["reflected"], t["reflected"], {})]))
        rows.append((y + 0.8, [(hb["guided"], t["guided"], {}), (hb["back"], t["reflected"], hatch),
                               (hb["face"], t["reflected"], {})]))
        ticks += [(y, "fiber"), (y + 0.8, "metal guide")]
        ax.text(-0.5, y - 0.52, f"{P.CASES[name]['title']}  ·  {P.CASES[name]['na_label']}", color=t["text"],
                fontsize=11.5, fontweight="bold", ha="left", va="bottom")
        y += 2.5
    for yy, segs in rows:
        left = 0.0
        for v, c, kw in segs:
            w = 100 * v
            ax.barh(yy, max(w - 0.25, 0.0), left=left, height=0.62, color=c, **({"linewidth": 0} | kw))
            if w > 6:
                ax.text(left + w / 2, yy, f"{w:.0f}%", ha="center", va="center", color=on, fontsize=10.5,
                        fontweight="bold", bbox=dict(boxstyle="round,pad=0.15", fc=c, ec="none") if kw else None)
            left += w
    ax.set_yticks([v for v, _ in ticks], [s for _, s in ticks])
    ax.set_ylim(y - 2.5 + 0.8 + 0.5, -1.0)
    ax.set_xlim(0, 100)
    ax.set_xlabel("percent of the incident power")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=t[k]) for k in ("guided", "cladding", "lost")]
    handles += [plt.Rectangle((0, 0), 1, 1, facecolor=t["reflected"], **hatch),
                plt.Rectangle((0, 0), 1, 1, color=t["reflected"])]
    labels = ["in the guided mode (core mode / TE11)", "fiber: uncoupled, near the core",
              "fiber: uncoupled, gone sideways", "metal: reflected after entering the horn",
              "reflected without entering (facet / metal face)"]
    ax.legend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(0.45, -0.13), fontsize=9.5,
              handlelength=1.4, columnspacing=1.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def field_rgb(frame, emax, keep, metal):
    v = style.signed_compress(mirror(np.asarray(frame, np.float32)[:, keep]), emax, 0.5)
    rgb = style.field_cmap()(0.5 + 0.5 * v)[..., :3]
    rgb[metal] = matplotlib.colors.to_rgb(METAL)
    return rgb


def movie(runs, out, fps=24, width=1920, xlim=7.0, zlim=ZLIM, t_max=75.0, hold=5.0):
    """The three launches side by side as the light switches on, fills the horn and settles;
    then a cross-fade to the steady-state split, with the power budget."""
    fs = width / 1920
    t = style.use("dark", font_scale=fs)
    refs = common_ref(runs)
    data = {}
    for name, run in runs.items():
        frames = np.load(os.path.join(run.dir, "frames_rz_Er.npy"), mmap_mode="r")
        times = np.load(os.path.join(run.dir, "frames_t.npy"))
        g = np.load(os.path.join(run.dir, "frames_grid.npz"))
        r, z = g["r"], g["z"]
        stride = int(run.meta["frame_stride"])
        keep = r <= xlim
        metal = mirror(metal_mask(run, stride)[:, keep])
        # each column's field scale: its strongest steady field (in the guide, where the power is
        # squeezed into the smallest area)
        last = np.asarray(frames[len(times) - 8:len(times)], np.float32)[:, z > zlim[0]][:, :, keep]
        air = ~metal_mask(run, stride)[z > zlim[0]][:, keep]
        _, _, comp = horn_image(run, xlim, refs[name], stride=stride)
        data[name] = dict(frames=frames, times=times, x=np.concatenate([-r[keep][::-1], r[keep]]), z=z, keep=keep,
                          metal=metal, emax=np.percentile(abs(last[:, air]), 99.5), comp=comp, b=split(run))
    n = min(int(np.searchsorted(d["times"], t_max)) for d in data.values())

    W, H = width, int(width * 9 / 16)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    top, bottom = 0.80, 0.135
    h = top - bottom
    w = h * H / W * (2 * xlim) / (zlim[1] - zlim[0])
    gap, left = 0.035, 0.05
    fig.text(left, 0.952, "Focusing a beam into a single-mode metal waveguide (Maxwell's equations, FDTD)",
             fontsize=19 * fs, color=t["text"], fontweight="bold", va="center")
    sub = fig.text(left, 0.905, "", fontsize=13 * fs, color=t["text2"], va="center")
    clock = fig.text(0.985, 0.952, "", fontsize=13 * fs, color=t["text2"], ha="right", va="center")
    ims, faders = [], []
    for i, (name, d) in enumerate(data.items()):
        x0 = left + i * (w + gap)
        ax = fig.add_axes([x0, bottom, w, h])
        im = ax.imshow(field_rgb(d["frames"][0], d["emax"], d["keep"], d["metal"]), origin="lower", aspect="auto",
                       interpolation="antialiased", extent=extent(d["x"], d["z"]))
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(zlim[1], zlim[0])
        ax.set_xticks([-5, 0, 5])
        ax.set_yticks([0, 5, 10, 15])
        ax.tick_params(labelsize=9 * fs)
        if i == 0:
            ax.set_ylabel("z / λ   (mouth at 0)", fontsize=10 * fs)
        else:
            ax.set_yticklabels([])
        ax.text(0.0, 1.055, P.CASES[name]["title"], transform=ax.transAxes, fontsize=14 * fs,
                color=t["text"], fontweight="bold", va="bottom")
        ax.text(0.0, 1.015, SUBTITLES[name], transform=ax.transAxes, fontsize=11 * fs, color=t["text2"],
                va="bottom")
        ims.append(im)
        b = d["b"]
        for k, (key, label) in enumerate([("guided", "into the guide (TE11)"), ("reflected", "reflected")]):
            yk = bottom - 0.065 - 0.035 * k
            faders.append(fig.text(x0, yk, "■", color=t[key], fontsize=15 * fs, va="center"))
            faders.append(fig.text(x0 + 0.018, yk, f"{pct(b[key]):>4}  {label}", color=t["text"],
                                   fontsize=12 * fs, va="center"))
    x_note = left + 3 * w + 2 * gap + 0.03
    note = fig.text(x_note, top, "", fontsize=11 * fs, color=t["text2"], va="top", linespacing=1.5, wrap=True)
    for f in faders:
        f.set_alpha(0.0)

    field_note = ("Red / blue: + / − field.\n\nThe metal is a perfect\nconductor. Behind the\nhorn, the guide is\n"
                  "single-mode: only TE11\ncan travel down it.\n\nEach column is scaled\nto its own strongest\n"
                  "steady field.")
    steady_note = ("Steady state,\ntime-averaged |E|²\n(log, 3 decades), one\nscale for all three.\n\n"
                   "White: incident beam.\nBlue: TE11 heading in.\nYellow: everything else,\nall of which ends up\n"
                   "reflected.\nGray: metal.")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    with writer.saving(fig, out, dpi=100):
        sub.set_text("Instantaneous electric field Eₓ in the plane of polarization, as the source switches on")
        note.set_text(field_note)
        for k in range(n):
            for im, d in zip(ims, data.values()):
                im.set_data(field_rgb(d["frames"][k], d["emax"], d["keep"], d["metal"]))
            clock.set_text(f"t = {data[next(iter(data))]['times'][k]:5.1f} periods")
            writer.grab_frame()
        last = [field_rgb(d["frames"][n - 1], d["emax"], d["keep"], d["metal"]) for d in data.values()]
        sub.set_text("Where the light goes: the guide's mode heading in (blue), and everything else (yellow)")
        clock.set_text("steady state")
        note.set_text(steady_note)
        n_fade, n_hold = int(1.2 * fps), int(hold * fps)
        for j in range(n_fade + n_hold):
            f = min(1.0, j / n_fade)
            f = f * f * (3 - 2 * f)
            for im, l0, d in zip(ims, last, data.values()):
                im.set_data((1 - f) * l0 + f * d["comp"])
            for it in faders:
                it.set_alpha(f)
            writer.grab_frame()
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "horn", f"res{P.HORN_RESOLUTION}"))
    ap.add_argument("--fiber-root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media"))
    ap.add_argument("--only", nargs="*", default=["side", "modes", "budget", "movie"])
    a = ap.parse_args()
    runs = {c: HornRun(a.root, c) for c in P.CASE_ORDER}
    os.makedirs(a.out, exist_ok=True)
    results = {}
    for name, run in runs.items():
        s = split(run)
        b = run.budget()
        results[name] = dict(w0=run.meta["w0"], beam_na=run.meta["beam_na"], res=run.meta["res"], **s,
                             along_face=b["along_face"], horn_flux=b["horn_flux"], horn_flux_z=b["horn_flux_z"])
        print(f"{name:8s} guided {s['guided']:.4f}  reflected {s['reflected']:.4f} (face {s['face']:.4f}, "
              f"turned back in the horn {s['back']:.4f})  total {s['total']:.4f};  entered {s['entered']:.4f}, "
              f"entered - back {s['entered'] - s['back_modes']:.4f}")
    with open(os.path.join(a.out, "horn_results.json"), "w") as f:
        json.dump(results, f, indent=1)
    if "side" in a.only:
        fig_side(runs, os.path.join(a.out, "horn_side_view.png"))
    for theme in ("dark", "light"):
        if "modes" in a.only:
            fig_modes(runs, os.path.join(a.out, f"horn_modes_{theme}.png"), theme)
        if "budget" in a.only:
            fibers = {c: load_case(a.fiber_root, c) for c in P.CASE_ORDER}
            fig_budget(runs, fibers, os.path.join(a.out, f"fiber_vs_horn_{theme}.png"), theme)
    if "movie" in a.only:
        movie(runs, os.path.join(a.out, "horn_movie.mp4"))


if __name__ == "__main__":
    main()
