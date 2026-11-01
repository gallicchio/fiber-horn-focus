"""Figures and the movies for the metal waveguide (horn) runs, with the fiber runs alongside.

    python scripts/horn_media.py [--sets smooth smooth_focus corrugated] [--fiber-root data/main/res16]
                                 [--only side modes budget movie]

Three sets of launches (params.HORN_SETS), each in data/horn/<set>/res120:
  smooth        the smooth horn, the fiber's three beams focused on its mouth
  smooth_focus  the smooth horn, the beams around its best waist, focused where that is best
  corrugated    the corrugated horn, the same for its best waist and focus

media/horn_side_view.png, horn_focus_side_view.png, corrugated_side_view.png
                                      steady state on the plane of polarization, the three launches
media/horn_modes_{dark,light}.png    the power in each of the smooth horn's own modes along it
                                      (the beams focused on its mouth)
media/fiber_vs_horn_{dark,light}.png  where the light goes: the fiber and the three sets
media/fiber_vs_corrugated_{dark,light}.png
                                      the same for the fiber and the corrugated horn alone (the README's)
media/horn_movie.mp4, horn_focus_movie.mp4, corrugated_movie.mp4
                                      the light switching on, then the steady-state picture
media/horn_results.json               every set's power budget, and checks

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
from fiberlaunch import beam  # noqa: E402
from fiberlaunch import frames  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402
from fiberlaunch.horn_analysis import MODES, HornRun, cutoff_radius  # noqa: E402
from fig_results import budget as fiber_budget  # noqa: E402
from horn_design import cylindrical  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
METAL = "#4a4a46"
SPOT = dict(narrow="3× the best spot", matched="the best spot", wide="⅓ of the best spot")
# file names, and the emit run (data/horn/design/<emit>_res<res>) that gives each horn's own field
SETS = dict(smooth=dict(name="horn", emit="emit"), smooth_focus=dict(name="horn_focus", emit="emit"),
            corrugated=dict(name="corrugated", emit="corrugated_emit"))


def wall_name(set_name):
    return "corrugated horn" if P.HORN_SETS[set_name]["wall"] == "corrugated" else "smooth horn"


def focus_text(set_name, short=False):
    f = P.HORN_SETS[set_name]["focus"]
    if short:
        return "focus on the mouth" if f == 0 else f"focus {f:g} λ in"
    return "focused on its mouth" if f == 0 else f"focused {f:g} λ inside it"


def subtitle(run):
    return f"{SPOT[run.name]} · beam NA {P.beam_na(run.meta['w0']):.2g}"


def zlim(run):
    """From the source plane to the absorber: mouth at 0, throat from L = 12.4 to L + 6."""
    return run.z_source, run.length + run.meta["throat_length"]


def load_set(set_name, res=P.HORN_RESOLUTION):
    root = os.path.join(ROOT, "data", "horn", set_name, f"res{res}")
    emit = os.path.join(ROOT, "data", "horn", "design", f"{SETS[set_name]['emit']}_res{res}")
    return {c: HornRun(root, c, emit=emit) for c in P.CASE_ORDER}


def clipped(run):
    """The part of the exact beam's power that misses the source plane (r > the box's radius)."""
    rho = np.linspace(0, 60.0, 2401)
    f = beam.family(rho, run.z_source - run.focus, run.meta["w0"])
    k = rho <= run.meta["half_width"]
    return 1 - beam.power(rho[k], [c[k] for c in f]) / beam.power(rho, f)


def mark_focus(ax, run, t, label=True):
    """A dotted line across the horn where the beams are focused (the waist), if inside it."""
    if run.focus <= 0:
        return
    a = float(run.radius(run.focus))
    ax.plot([-a, a], [run.focus, run.focus], color="white", lw=0.9, ls=(0, (1.5, 2)), alpha=0.85)
    if label:
        ax.text(a + 0.35, run.focus, "focus", color=t["text2"], fontsize=8.5, va="center", ha="left")


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
    """The horn's budget, with the reflected light split by where the incident beam lands: on the
    metal face around the mouth (all of which bounces off), or inside the mouth (and then, but for
    what reaches the guide, turned back inside the horn)."""
    b = run.budget()
    land = run.landing()
    face = min(land["face"], b["reflected"])
    return dict(guided=b["guided"], back=b["reflected"] - face, face=face, reflected=b["reflected"],
                landed_inside=land["inside"], total=b["total"])


def label_horn(ax, t, run):
    kw = dict(color=t["text2"], fontsize=8.5, ha="left", va="center")
    ax.text(-5.6, 4.0, "metal", **kw)
    ax.text(run.mouth * 0.55 + 0.3, 0.5 * run.length, "horn", **kw)
    ax.text(P.THROAT_RADIUS + 0.35, run.length + 3.0, "single-mode\nguide", **kw)


def fig_side(runs, out, title, xlim=7.0):
    """The three launches side by side, the light coming down onto the mouth."""
    t = style.use("dark")
    refs = common_ref(runs)
    fig, axs = plt.subplots(1, len(runs), figsize=(12.0, 9.6), sharey=True)
    for i, (ax, (name, run)) in enumerate(zip(np.atleast_1d(axs), runs.items())):
        x, z, img = horn_image(run, xlim, refs[name])
        ax.imshow(img, origin="lower", aspect="equal", extent=extent(x, z), interpolation="bilinear")
        ax.set_xlim(-xlim, xlim)
        lo, hi = zlim(run)
        ax.set_ylim(hi, lo)
        mark_focus(ax, run, t, label=(i == len(runs) - 1))
        b = split(run)
        ax.text(0, 1.105, P.CASES[name]["title"], transform=ax.transAxes, color=t["text"], fontsize=13,
                fontweight="bold", va="bottom")
        ax.text(0, 1.062, subtitle(run), transform=ax.transAxes, color=t["text2"], fontsize=10, va="bottom")
        ax.text(0, 1.012, f"into the guide {pct(b['guided'])}  ·  reflected {pct(b['reflected'])}",
                transform=ax.transAxes, color=t["text"], fontsize=10.5, va="bottom")
        ax.set_xlabel("x / λ")
        ax.set_xticks([-5, 0, 5])
        ax.set_yticks([-5, 0, 5, 10, 15])
        if i == 0:
            ax.set_ylabel("z / λ   (mouth at z = 0; the light comes down from the top)")
            label_horn(ax, t, run)
    fig.text(0.06, 0.975, title, color=t["text"], fontsize=14, fontweight="bold", va="top")
    fig.text(0.01, 0.012, "Time-averaged |E|² on the plane of polarization, log scale (3 decades), one scale for all "
             "three (per watt in).  White: the incident beam.\nYellow: reflected light, and the light inside the "
             "horn that does not reach the guide (all of it turns back).  Blue: the light that reaches the guide's "
             "mode (TE11).  Gray: metal.", color=t["text2"], fontsize=9.5, va="bottom")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.10, wspace=0.08)
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


def fig_budget(sets, fibers, out, theme):
    """Stacked bars for each launch, the fiber above each set of metal-guide launches: where the
    light goes."""
    t = style.use(theme)
    fig, ax = plt.subplots(figsize=(11.5, 2.0 + 1.55 * (1 + len(sets))))

    def ink(c):  # whichever of white and near-black text reads better on the segment's colour
        lum = lambda rgb: np.dot([0.2126, 0.7152, 0.0722],
                                 [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in rgb])
        L = lum(matplotlib.colors.to_rgb(c))
        return "#ffffff" if 1.05 / (L + 0.05) >= (L + 0.05) / (lum((0.1, 0.1, 0.1)) + 0.05) else "#1a1a1a"

    hatch = dict(hatch="////", edgecolor=t["surface"], linewidth=0)
    rows, ticks = [], []
    y = 0.0
    guided = ("guided", "guided")  # (long, short) labels, written in the segment when they fit
    for name in P.CASE_ORDER:
        fb = fiber_budget(fibers[name])
        rows.append((y, [(fb["guided"], t["guided"], {}, guided),
                         (fb["cladding"], t["cladding"], {}, ("forward, beside the core", "beside the core")),
                         (fb["lost"], t["lost"], {}, ("out through the side", "out the side")),
                         (fb["reflected"], t["reflected"], {}, None)]))
        ticks.append((y, "fiber"))
        for k, (set_name, runs) in enumerate(sets.items()):
            hb = split(runs[name])
            yy = y + 0.8 * (k + 1)
            rows.append((yy, [(hb["guided"], t["guided"], {}, guided),
                              (hb["back"], t["reflected"], hatch, ("turned back inside the horn", "turned back")),
                              (hb["face"], t["reflected"], {}, ("bounced off the face", "off the face"))]))
            ticks.append((yy, f"{wall_name(set_name)}, {focus_text(set_name, short=True)}"))
        ax.text(-0.5, y - 0.52, f"{P.CASES[name]['title']}  ·  {SPOT[name]}", color=t["text"],
                fontsize=11.5, fontweight="bold", ha="left", va="bottom")
        y += 0.8 * (1 + len(sets)) + 0.9
    ax.set_xlim(0, 100)
    renderer = fig.canvas.get_renderer()
    for yy, segs in rows:
        left = 0.0
        for v, c, kw, names in segs:
            w = 100 * v
            ax.barh(yy, max(w - 0.25, 0.0), left=left, height=0.62, color=c, **({"linewidth": 0} | kw))
            room = ax.transData.transform((left + w - 1.2, yy))[0] - ax.transData.transform((left + 1.2, yy))[0]
            for s in ([f"{w:.0f}% {n}" for n in names] if names else []) + [f"{w:.0f}%"]:
                label = ax.text(left + w / 2, yy, s, ha="center", va="center", color=ink(c), fontsize=10.5,
                                fontweight="bold",
                                bbox=dict(boxstyle="round,pad=0.15", fc=c, ec="none") if kw else None)
                if label.get_window_extent(renderer).width <= room and w > 6:
                    break
                label.remove()
            left += w
    ax.set_yticks([v for v, _ in ticks], [s for _, s in ticks])
    ax.set_ylim(y - 0.9 - 0.8 + 0.5, -1.0)
    ax.set_xlim(0, 100)
    ax.set_xlabel("percent of the incident power")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=t[k]) for k in ("guided", "cladding", "lost")]
    handles += [plt.Rectangle((0, 0), 1, 1, facecolor=t["reflected"], **hatch),
                plt.Rectangle((0, 0), 1, 1, color=t["reflected"])]
    labels = ["in the guided mode (core mode / TE11)", "fiber: forward in the cladding, beside the core",
              "fiber: out through the side of the cladding", "metal: entered the horn, turned back inside it",
              "reflected at the entrance (fiber tip / metal face)"]
    # below the axis label, the same distance however many rows there are
    below = matplotlib.transforms.offset_copy(ax.transAxes, fig=fig, y=-36, units="points")
    ax.legend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(0.4, 0), bbox_transform=below,
              fontsize=9.5, handlelength=1.4, columnspacing=1.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def field_rgb(frame, emax, metal):
    v = style.signed_compress(mirror(frame), emax, 0.5)
    rgb = style.field_cmap()(0.5 + 0.5 * v)[..., :3]
    rgb[metal] = matplotlib.colors.to_rgb(METAL)
    return rgb


def movie(runs, out, title, fps=30, per_second=2.0, width=1920, xlim=8.0, t_max=70.0, hold=5.0):
    """The three launches side by side as the light switches on, fills the horn and settles (by
    t_max, to 2%), shown at per_second periods a second, with the frames between the saved ones
    from fiberlaunch.frames; then a cross-fade to the steady-state split, with the power budget."""
    fs = width / 1920
    t = style.use("dark", font_scale=fs)
    refs = common_ref(runs)
    zl = zlim(next(iter(runs.values())))
    data = {}
    for name, run in runs.items():
        saved = np.load(os.path.join(run.dir, "frames_rz_Er.npy"), mmap_mode="r")
        times = np.load(os.path.join(run.dir, "frames_t.npy"))
        n = int(np.searchsorted(times, t_max)) + 10  # and the frames the fits near t_max need
        g = np.load(os.path.join(run.dir, "frames_grid.npz"))
        r, z = g["r"], g["z"]
        stride = int(run.meta["frame_stride"])
        keep = r <= xlim
        metal = mirror(metal_mask(run, stride)[:, keep])
        # each column's field scale: its strongest steady field (in the guide, where the power is
        # squeezed into the smallest area)
        last = np.asarray(saved[len(times) - 8:len(times)], np.float32)[:, z > zl[0]][:, :, keep]
        air = ~metal_mask(run, stride)[z > zl[0]][:, keep]
        _, _, comp = horn_image(run, xlim, refs[name], stride=stride)
        data[name] = dict(saved=np.asarray(saved[:n], np.float32)[:, :, keep], times=times[:n],
                          x=np.concatenate([-r[keep][::-1], r[keep]]), z=z, metal=metal,
                          emax=np.percentile(abs(last[:, air]), 99.5), comp=comp, b=split(run))
    t_out = np.arange(0, t_max, per_second / fps)

    W, H = width, int(width * 9 / 16)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    top, bottom = 0.80, 0.135
    h = top - bottom
    w = h * H / W * (2 * xlim) / (zl[1] - zl[0])
    gap, left = 0.035, 0.05
    fig.text(left, 0.952, title, fontsize=19 * fs, color=t["text"], fontweight="bold", va="center")
    sub = fig.text(left, 0.905, "", fontsize=13 * fs, color=t["text2"], va="center")
    clock = fig.text(0.985, 0.952, "", fontsize=13 * fs, color=t["text2"], ha="right", va="center")
    ims, faders = [], []
    for i, (name, d) in enumerate(data.items()):
        x0 = left + i * (w + gap)
        ax = fig.add_axes([x0, bottom, w, h])
        im = ax.imshow(field_rgb(d["saved"][0], d["emax"], d["metal"]), origin="lower", aspect="auto",
                       interpolation="antialiased", extent=extent(d["x"], d["z"]))
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(zl[1], zl[0])
        ax.set_xticks([-5, 0, 5])
        ax.set_yticks([-5, 0, 5, 10, 15])
        ax.tick_params(labelsize=9 * fs)
        mark_focus(ax, runs[name], t, label=False)
        if i == 0:
            ax.set_ylabel("z / λ   (mouth at 0)", fontsize=10 * fs)
        else:
            ax.set_yticklabels([])
        ax.text(0.0, 1.055, P.CASES[name]["title"], transform=ax.transAxes, fontsize=14 * fs,
                color=t["text"], fontweight="bold", va="bottom")
        ax.text(0.0, 1.015, subtitle(runs[name]), transform=ax.transAxes, fontsize=11 * fs, color=t["text2"],
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

    run0 = next(iter(runs.values()))
    wall_note = ("The horn's wall is cut\ninto slots a quarter\nwave deep.\n\n" if run0.wall == "corrugated" else "")
    focus_note = (f"Dotted: the focus,\n{run0.focus:g} λ inside the horn.\n\n" if run0.focus > 0 else
                  "The beams are focused\non the mouth.\n\n")
    field_note = ("Red / blue: + / − field.\n\n" + focus_note + wall_note + "The metal is a perfect\nconductor. "
                  "Behind the\nhorn, the guide is\nsingle-mode: only TE11\ncan travel down it.\n\n"
                  "Each column is scaled\nto its own strongest\nsteady field.\n\n"
                  f"Shown at {per_second:g} periods\na second.")
    steady_note = ("Steady state,\ntime-averaged |E|²\n(log, 3 decades), one\nscale for all three.\n\n"
                   "White: incident beam.\nBlue: the light that\nreaches the guide.\nYellow: everything else,\n"
                   "all of which ends up\nreflected.\nGray: metal.")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow"])
    with writer.saving(fig, out, dpi=100):
        sub.set_text("Instantaneous electric field Eₓ in the plane of polarization, as the source switches on")
        note.set_text(field_note)
        for tt in t_out:
            for im, d in zip(ims, data.values()):
                im.set_data(field_rgb(frames.at(d["saved"], d["times"], tt), d["emax"], d["metal"]))
            clock.set_text(f"t = {tt:5.1f} periods")
            writer.grab_frame()
        last = [field_rgb(frames.at(d["saved"], d["times"], t_out[-1]), d["emax"], d["metal"])
                for d in data.values()]
        sub.set_text("Where the light goes: the light that reaches the guide (blue), and everything else (yellow)")
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
    ap.add_argument("--sets", nargs="*", default=list(P.HORN_SETS))
    ap.add_argument("--res", type=int, default=P.HORN_RESOLUTION)
    ap.add_argument("--fiber-root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media"))
    ap.add_argument("--only", nargs="*", default=["side", "modes", "budget", "movie"])
    a = ap.parse_args()
    sets = {s: load_set(s, a.res) for s in a.sets}
    os.makedirs(a.out, exist_ok=True)
    results = {}
    for set_name, runs in sets.items():
        results[set_name] = {}
        emit = next(iter(runs.values())).emit_dir
        eta = cylindrical(emit, json.load(open(os.path.join(emit, "meta.json"))))[0]
        for name, run in runs.items():
            sp, b = split(run), run.budget()
            r = dict(w0=run.meta["w0"], beam_na=run.meta["beam_na"], focus=run.focus, wall=run.wall,
                     res=run.meta["res"], **sp, along_face=b["along_face"], horn_flux=b["horn_flux"],
                     horn_flux_z=b["horn_flux_z"], clipped=clipped(run), reciprocity=eta(run.meta["w0"], run.focus),
                     # the same for the beam the finite source plane launches, which the direct run sees
                     reciprocity_launched=eta(run.meta["w0"], run.focus, beam.launched_spectrum(
                         run.meta["w0"], run.z_source - run.focus, run.meta["half_width"])))
            if run.wall == "smooth":  # the older split, by the modes heading in and out just inside the mouth
                e = run.entry()
                r.update(entered_modes=e["entered"], back_modes=e["back"])
            results[set_name][name] = r
            print(f"{set_name:12s} {name:8s} guided {r['guided']:.4f} (reciprocity: the launched beam "
                  f"{r['reciprocity_launched']:.4f}, the whole beam {r['reciprocity']:.4f}; the part of it "
                  f"beyond the source plane {r['clipped']:.4f})  reflected "
                  f"{r['reflected']:.4f} (face {r['face']:.4f}, turned back in the horn {r['back']:.4f})  "
                  f"total {r['total']:.4f}" + (f";  by the modes: turned back {r['back_modes']:.4f}"
                                               if "back_modes" in r else ""))
    with open(os.path.join(a.out, "horn_results.json"), "w") as f:
        json.dump(results, f, indent=1)
    for set_name, runs in sets.items():
        name = SETS[set_name]["name"]
        title = f"Into a single-mode metal waveguide: a {wall_name(set_name)}, the beams {focus_text(set_name)}"
        if "side" in a.only:
            fig_side(runs, os.path.join(a.out, f"{name}_side_view.png"), title)
        if "modes" in a.only and set_name == "smooth":  # the modes are a smooth pipe's; one set shows them
            for theme in ("dark", "light"):
                fig_modes(runs, os.path.join(a.out, f"{name}_modes_{theme}.png"), theme)
        if "movie" in a.only:
            movie(runs, os.path.join(a.out, f"{name}_movie.mp4"),
                  f"Focusing a beam into a single-mode metal waveguide: {wall_name(set_name)}, "
                  f"{focus_text(set_name, short=True)}")
    if "budget" in a.only:
        fibers = {c: load_case(a.fiber_root, c) for c in P.CASE_ORDER}
        for theme in ("dark", "light"):
            fig_budget(sets, fibers, os.path.join(a.out, f"fiber_vs_horn_{theme}.png"), theme)
            if "corrugated" in sets:  # the README's: the fiber, and the horn whose movie it shows
                fig_budget({"corrugated": sets["corrugated"]}, fibers,
                           os.path.join(a.out, f"fiber_vs_corrugated_{theme}.png"), theme)


if __name__ == "__main__":
    main()
