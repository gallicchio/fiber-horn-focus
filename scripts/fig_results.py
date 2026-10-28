"""Static figures from the FDTD campaign: the steady-state split of the light (side view),
the power budget, and the guided / forward power along the fiber.

    python scripts/fig_results.py [--root data/main/res16] [--out media]
"""

import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from fiberlaunch import fibermode as fm  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import load_case  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def budget(run):
    """Where the incident power is at the end of the fiber (last mode plane)."""
    tab = run.mpb_table()
    pb = run.plane_budget()
    A, _, _ = run.guided_fit()
    guided = abs(A) ** 2 / run.P_inc
    z_end = tab["z"][-1]
    total = tab["total"][-1]
    R = run.R
    return dict(z=z_end, guided=guided, cladding=total - guided, lost=max(0.0, 1 - R - total),
                reflected=R, guided_mpb=tab["guided"][-1], plane_z=pb["z"], plane_guided=pb["guided"],
                plane_total=pb["total"])


def pct(v):
    """A fraction as a percentage: whole numbers from 10%, one decimal below."""
    return f"{100 * v:.0f}%" if v >= 0.095 else f"{100 * v:.1f}%"


def common_ref(runs, xlim, z_from=3.0):
    """One brightness scale for several launches, per unit of incident power: the brightest
    point z_from or more inside the glass over all of them, so equal guided power looks equal."""
    per_watt = 0.0
    for run in runs.values():
        s = run.side_split("xz")
        keep = abs(s["x"]) <= xlim
        per_watt = max(per_watt, s["total"][keep][:, s["z"] >= z_from].max() / run.P_inc)
    return {name: per_watt * run.P_inc for name, run in runs.items()}


def side_image(run, xlim, decades=3.0, ref=None):
    s = run.side_split("xz")
    x, z = s["x"], s["z"]
    keep = (abs(x) <= xlim)
    if ref is None:
        ref = s["total"][keep][:, z > -1.0].max()  # the brightest point from the facet on
    air = z < 0  # before the facet nothing is guided yet: show it as plain light in air
    img = style.two_channel(s["guided"][keep], np.where(air, 0, s["rest"][keep]), ref, decades,
                            Ia=np.where(air, s["total"][keep], 0))
    return x[keep], z, img


def draw_fiber(ax, x0, x1, z1, color, lw=0.8):
    """Facet and core outline on an x (vertical) - z (horizontal) image."""
    ax.plot([0, 0], [x0, x1], color=color, lw=lw, alpha=0.8)
    for s in (-1, 1):
        ax.plot([0, z1], [s * P.CORE_RADIUS] * 2, color=color, lw=lw, alpha=0.6)


def fig_side(runs, out, xlim=12.0):
    """Steady state, x-z plane through the axis: guided part blue, everything else orange."""
    t = style.use("dark")
    fig_h = 2.15 * len(runs) + 1.0
    fig, axs = plt.subplots(len(runs), 1, figsize=(13, fig_h), sharex=True)
    refs = common_ref(runs, xlim)
    for ax, (name, run) in zip(np.atleast_1d(axs), runs.items()):
        x, z, img = side_image(run, xlim, ref=refs[name])
        ax.imshow(img, origin="lower", aspect="equal",
                  extent=[z[0], z[-1], x[0], x[-1]], interpolation="bilinear")
        draw_fiber(ax, x[0], x[-1], z[-1], t["muted"])
        b = budget(run)
        c = P.CASES[name]
        ax.set_title(f"{c['title']}  ·  {c['subtitle']}", loc="left", color=t["text"], fontsize=12)
        ax.text(z[-1] - 1, x[-1] - 1.2,
                f"core mode {pct(b['guided'])}   ·   uncoupled near the core {pct(b['cladding'])}"
                f"   ·   gone sideways {pct(b['lost'])}   ·   reflected {pct(b['reflected'])}",
                ha="right", va="top", color=t["text2"], fontsize=10)
        ax.set_ylabel("x / λ")
        ax.set_yticks([-10, 0, 10])
    np.atleast_1d(axs)[-1].set_xlabel("z / λ   (facet at z = 0, light travels left to right)")
    fig.text(0.01, 0.005, "Time-averaged |E|², log scale (3 decades), the same scale for every row (per watt in).  "
             "White: light in air.  Blue: the fiber's guided mode.  Orange: everything else.",
             color=t["text2"], fontsize=10)
    fig.tight_layout(rect=(0, 0.3 / fig_h, 1, 1))  # 0.3 inch for the footnote
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_budget(runs, out, theme):
    """Horizontal stacked bars: where the light ends up, per launch."""
    t = style.use(theme)
    keys = ["guided", "cladding", "lost", "reflected"]
    fig, ax = plt.subplots(figsize=(10, 3.4))
    names = list(runs)
    for i, name in enumerate(names):
        b = budget(runs[name])
        left = 0.0
        for k in keys:
            w = 100 * b[k]
            ax.barh(i, max(w - 0.25, 0.0), left=left, height=0.62, color=t[k], linewidth=0)
            if w > 6:
                ax.text(left + w / 2, i, f"{w:.0f}%", ha="center", va="center",
                        color=t["surface"] if theme == "dark" else "#ffffff", fontsize=11, fontweight="bold")
            left += w
    ax.set_yticks(range(len(names)), [f"{P.CASES[n]['title']}\n{P.CASES[n]['na_label']}" for n in names])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel(f"percent of the incident power, {budget(runs[names[0]])['z']:.0f} λ into the fiber")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=t[k]) for k in keys]
    ax.legend(handles, [style.LABELS[k] for k in keys], ncol=2, loc="upper center",
              bbox_to_anchor=(0.45, -0.32), fontsize=10, handlelength=1.2, columnspacing=1.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_along_z(runs, out, theme):
    """Small multiples: guided power and total forward power in the window against z."""
    t = style.use(theme)
    fig, axs = plt.subplots(1, len(runs), figsize=(12, 3.6), sharey=True)
    for ax, (name, run) in zip(np.atleast_1d(axs), runs.items()):
        b = budget(run)
        z = b["plane_z"]
        inside = z > 0.5
        ax.plot(z[inside], 100 * b["plane_guided"][inside], color=t["guided"], label=style.LABELS["guided"])
        ax.plot(z[inside], 100 * b["plane_total"][inside], color=t["cladding"],
                label=f"all forward light in the {2 * run.meta['half_width']:.0f} λ × "
                      f"{2 * run.meta['half_width']:.0f} λ window")
        ax.fill_between(z[inside], 100 * b["plane_guided"][inside], 100 * b["plane_total"][inside],
                        color=t["cladding"], alpha=0.15, linewidth=0)
        ax.set_title(P.CASES[name]["title"], loc="left")
        ax.set_xlabel("z / λ")
        ax.set_ylim(0, 102)
        ax.grid(True, axis="y")
        ax.text(z[inside][-1], 100 * b["plane_guided"][inside][-1] - 3, f"{100 * b['guided']:.1f}%",
                ha="right", va="top", color=t["text"], fontsize=10)
    np.atleast_1d(axs)[0].set_ylabel("percent of incident power")
    np.atleast_1d(axs)[0].legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def rest_angles(run, z_at=3.0):
    """Angular distribution (in the glass) of the uncoupled light just inside the facet:
    the transverse spectrum of E - a e, with a the guided amplitude projected at that
    plane, as a cumulative power fraction (|E_t|^2 weighting, no obliquity factor)."""
    x, y, z, mm = run.planes()
    i = int(np.argmin(abs(z - z_at)))
    a = run.plane_budget()["a"][i]
    e = run.mode.f
    dx = x[1] - x[0]
    n = 8 * len(x)  # zero-pad: fine angular resolution near the axis
    kx = 2 * np.pi * np.fft.fftfreq(n, dx)
    kt = np.hypot(*np.meshgrid(kx, kx, indexing="ij"))
    spec = 0
    for c in ("Ex", "Ey"):
        R = np.asarray(mm[c][i], dtype=np.complex128) - a * e[c]
        spec = spec + abs(np.fft.fft2(R, (n, n))) ** 2
    k = 2 * np.pi * P.N_CLAD / P.WAVELENGTH
    s = np.clip(kt / k, 0, 1)
    theta = np.degrees(np.arcsin(s)).ravel()
    w = (spec * (kt < k)).ravel()
    order = np.argsort(theta)
    return theta[order], np.cumsum(w[order]) / w.sum()


def fig_rest_angles(runs, out, theme):
    t = style.use(theme)
    from fig_theory import CASE_COLORS
    cc = CASE_COLORS[theme]
    acc = np.degrees(np.arcsin(P.NA_FIBER / P.N_CLAD))
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    ax.axvline(acc, color=t["muted"], lw=1)
    ax.text(acc + 0.6, 4, f"core acceptance angle\n(in the glass) {acc:.1f}°", color=t["text2"], fontsize=9.5)
    label_at = dict(narrow=0.8, matched=0.55, wide=0.3)  # staggered, so the labels never collide
    for name, run in runs.items():
        th, cum = rest_angles(run)
        ax.plot(th, 100 * cum, color=cc[name], lw=2)
        y = label_at.get(name, 0.5)
        j = np.searchsorted(cum, y)
        inside = 100 * np.interp(acc, th, cum)
        b = budget(run)
        ax.text(th[j] + 0.8, 100 * y, f"{P.CASES[name]['title']}: {pct(b['cladding'] + b['lost'])} uncoupled,\n"
                f"{inside:.0f}% of that inside the cone", color=t["text"], fontsize=9.5, va="center")
    ax.set_xlim(0, 45)
    ax.set_ylim(0, 101)
    ax.set_xlabel("angle to the fiber axis, inside the glass (degrees)")
    ax.set_ylabel("% of the uncoupled light\nwithin this angle")
    ax.set_title("Which way the uncoupled light travels (just inside the facet)", loc="left")
    ax.grid(True, axis="y")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = {c: load_case(a.root, c) for c in a.cases if os.path.exists(os.path.join(a.root, c, "meta.json"))}
    results = {}
    for name, run in runs.items():
        b = budget(run)
        A, beta, sc = run.guided_fit()
        w0 = P.CASES[name]["w0"]
        # Yee dispersion: the grid's beta belongs to the physical n_eff given by
        # sin(beta dz / 2) / dz = n sin(omega dt / 2) / dt
        dz, dt, om = 1.0 / run.meta["res"], run.meta["dt"], 2 * np.pi * P.FREQ
        n_eff_phys = (np.sin(beta * dz / 2) / dz) / (np.sin(om * dt / 2) / dt)
        T = fm.fresnel_T(fm.lp01()["neff"])
        exact, _ = run.exact_beam_theory()
        results[name] = dict(
            beam_na=P.CASES[name]["na"], w0=w0, z_end=b["z"], reflected=b["reflected"],
            guided=b["guided"], guided_mpb=b["guided_mpb"], guided_rms_scatter_along_z=sc,
            uncoupled_near_axis=b["cladding"], uncoupled_beyond_window=b["lost"],
            window_half_width=run.meta["half_width"], n_eff_fdtd_grid=beta / 2 / np.pi,
            n_eff_fdtd_dispersion_corrected=n_eff_phys, n_eff_he11=fm.he11_neff(),
            theory_mode_overlap_x_fresnel=fm.gaussian_overlap(w0) * T,
            theory_ray_optics_x_fresnel=fm.ray_optics(w0) * T,
            theory_exact_beam_fresnel_overlap=exact,
            P_inc=run.P_inc, wall_hours=run.meta["wall_seconds"] / 3600)
        print(f"{name:8s} R {b['reflected']:.4f}  guided {b['guided']:.4f} (MPB {b['guided_mpb']:.4f}, "
              f"scatter {sc:.1e})  cladding {b['cladding']:.4f}  beyond window {b['lost']:.4f}  "
              f"n_eff(FDTD grid) {beta / 2 / np.pi:.5f} -> corrected {n_eff_phys:.5f} (HE11 {fm.he11_neff():.5f})\n"
              f"{'':8s} theory: paraxial overlap x Fresnel {fm.gaussian_overlap(w0) * T:.4f}, "
              f"exact beam x Fresnel(s, p) x overlap {exact:.4f}")
    with open(os.path.join(a.out, "results.json"), "w") as f:
        json.dump(results, f, indent=1)
    fig_side(runs, os.path.join(a.out, "steady_side_view.png"))
    for theme in ("light", "dark"):
        fig_budget(runs, os.path.join(a.out, f"budget_{theme}.png"), theme)
        fig_along_z(runs, os.path.join(a.out, f"power_along_z_{theme}.png"), theme)
        fig_rest_angles(runs, os.path.join(a.out, f"uncoupled_angles_{theme}.png"), theme)


if __name__ == "__main__":
    main()
