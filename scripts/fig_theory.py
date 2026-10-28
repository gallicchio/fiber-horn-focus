"""Theory figures, with the FDTD points added when the sweep has been run.

* facet_profiles: the three beams on the facet (spot too big, matched, too small), each split
  into a guided part and "the rest";
* efficiency_vs_na: coupling against beam NA - wave optics (mode overlap with the paraxial
  Gaussian, and with the exact focused beam the FDTD launches), ray optics (acceptance
  cone), and 3D FDTD.

    python scripts/fig_theory.py [--sweep data/sweep/res16] [--out media]
"""

import argparse
import glob
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import fibermode as fm  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.analysis import Mode, exact_beam_coupling  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
CASE_COLORS = dict(light=dict(narrow="#e87ba4", matched="#008300", wide="#6250d6"),
                   dark=dict(narrow="#d55181", matched="#008300", wide="#9085e9"))


def unit_gaussian(x, w0):
    """Field amplitude of a 2D Gaussian beam carrying unit power (on the line y = 0)."""
    return np.sqrt(2 / (np.pi * w0**2)) * np.exp(-(x**2) / w0**2)


def unit_mode(x):
    r = np.linspace(0, 60, 600001)
    psi = fm.lp01_field(r)
    norm = np.sqrt(2 * np.pi * np.sum(psi**2 * r) * (r[1] - r[0]))
    return fm.lp01_field(np.abs(x)) / norm


def fig_profiles(out, theme):
    """Each launch on the facet, split into the part with the shape of the fiber mode (guided)
    and the rest.  Shared axes, so the spot sizes and amplitudes compare directly."""
    t = style.use(theme)
    cc = CASE_COLORS[theme]
    x = np.linspace(-10, 10, 2001)
    m = unit_mode(x)
    fig, axs = plt.subplots(1, len(P.CASE_ORDER), figsize=(13.5, 4.4), sharey=True)
    for ax, name in zip(axs, P.CASE_ORDER):
        c = P.CASES[name]
        g = unit_gaussian(x, c["w0"])
        eta = fm.gaussian_overlap(c["w0"])
        guided = np.sqrt(eta) * m  # both carry unit power: the guided amplitude is the overlap integral
        rest = g - guided
        s = 1 / g.max()
        ax.axvspan(-P.CORE_RADIUS, P.CORE_RADIUS, color=t["glass"], lw=0, zorder=0)
        ax.text(0, 0.02, "core", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                color=t["text2"], fontsize=9)
        ax.axhline(0, color=t["muted"], lw=0.8)
        ax.plot(x, s * g, color=cc[name], lw=5, alpha=0.6, solid_capstyle="butt")  # wide, so it shows under
        ax.plot(x, s * guided, color=t["guided"], lw=2)
        ax.plot(x, s * rest, color=t["cladding"], lw=2)
        ax.set_title(f"{c['title']}\n", loc="left", fontweight="bold")
        ax.text(0, 1.02, c["subtitle"], transform=ax.transAxes, color=t["text2"], fontsize=10)
        ax.text(0.97, 0.96, f"guided   {100 * eta:4.1f}%\nthe rest {100 * (1 - eta):4.1f}%",
                transform=ax.transAxes, ha="right", va="top", fontsize=10, color=t["text"],
                family="DejaVu Sans Mono", linespacing=1.4)
        ax.set_xlim(-10, 10)
        ax.set_xticks([-10, -5, 0, 5, 10])
        ax.set_xlabel("x / λ   (on the facet)")
    axs[0].set_ylim(-1.0, 2.05)
    axs[0].set_ylabel("field, incident peak = 1")
    handles = [plt.Line2D([], [], color=t["text2"], lw=5, alpha=0.6),
               plt.Line2D([], [], color=t["guided"], lw=2), plt.Line2D([], [], color=t["cladding"], lw=2)]
    fig.legend(handles, ["incident field (the launch's own colour)",
                         "guided part: the fiber mode's share of it",
                         "the rest: incident minus guided"],
               ncol=3, loc="lower center", fontsize=10, handlelength=1.8, columnspacing=2.0)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def load_sweep(root):
    pts = []
    for d in sorted(glob.glob(os.path.join(root, "w*"))):
        try:
            air = json.load(open(os.path.join(d, "air", "meta.json")))
            fib = json.load(open(os.path.join(d, "fiber", "meta.json")))
        except FileNotFoundError:
            continue
        pin = air["P_net_refl_plane"]
        last = fib["mode_planes"][-1]  # P_fwd - P_bwd: the staggered-grid backward artifact removed
        pts.append(dict(w0=fib["w0"], na=P.beam_na(fib["w0"]), eta=(last["P_fwd"] - last["P_bwd"]) / pin,
                        R=1 - fib["P_net_refl_plane"] / pin))
    return sorted(pts, key=lambda p: p["na"])


def exact_beam_pts(sweep_root, main_root):
    """The non-paraxial wave-optics estimate (analysis.exact_beam_coupling) for every
    incident beam that has an air run: the sweep's and the three main launches'."""
    dirs = sorted(glob.glob(os.path.join(sweep_root, "w*", "air"))) + \
        [os.path.join(main_root, f"{c}_air") for c in P.CASE_ORDER]
    dirs = [d for d in dirs if os.path.exists(os.path.join(d, "incident_facet.npz"))]
    if not dirs:
        return []
    mode = Mode(os.path.join(main_root, "mode"))
    pts = []
    for d in dirs:
        w0 = json.load(open(os.path.join(d, "meta.json")))["w0"]
        eta, _ = exact_beam_coupling(np.load(os.path.join(d, "incident_facet.npz")), mode, fm.he11_neff())
        pts.append(dict(na=P.beam_na(w0), w0=w0, eta=eta))
    return sorted(pts, key=lambda p: p["na"])


def fig_efficiency(out, theme, sweep, main_pts, exact=()):
    t = style.use(theme)
    cc = CASE_COLORS[theme]
    na = np.geomspace(0.03, 1.2, 400)
    w0 = P.waist_for_na(na)
    T = fm.fresnel_T(fm.lp01()["neff"])
    wave = np.array([fm.gaussian_overlap(w) for w in w0]) * T
    ray = fm.ray_optics(w0) * T
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.axvline(P.NA_FIBER, color=t["muted"], lw=1)
    ax.text(P.NA_FIBER * 1.04, 102, f"fiber NA {P.NA_FIBER}", color=t["text2"], fontsize=9.5, va="top")
    ax.plot(na, 100 * wave, color=t["text"], lw=2)
    ax.plot(na, 100 * ray, color=t["muted"], lw=2, ls=(0, (1.5, 2)))
    if len(exact) > 1:
        ax.plot([p["na"] for p in exact], [100 * p["eta"] for p in exact], color=t["text2"], lw=1.5,
                ls=(0, (5, 3)))
        j = int(np.argmax([p["na"] for p in exact]))
        ax.text(exact[j]["na"] * 1.05, 100 * exact[j]["eta"] + 1.5,
                "the same, for the exact\nfocused beam (not paraxial)", color=t["text2"], fontsize=9.5,
                va="bottom")
    i = np.argmin(abs(na - 0.42))
    ax.text(na[i] * 0.97, 100 * wave[i] - 2, "wave optics: Gaussian\nmode overlap × Fresnel", color=t["text"],
            fontsize=10, va="top", ha="right")
    ax.text(0.105, 22, "ray optics: the rays that hit the core\ninside the acceptance cone",
            color=t["text2"], fontsize=10, va="top", ha="left")
    if sweep:
        ax.plot([p["na"] for p in sweep], [100 * p["eta"] for p in sweep], "o", ms=8,
                color=t["text"], mfc=t["surface"], mew=2, label="3D FDTD (Meep), 16 px/λ")
    for name, p in main_pts.items():
        ax.plot([p["na"]], [100 * p["eta"]], "o", ms=10, color=cc[name], mec=t["surface"], mew=1.5)
        if name == "wide":  # up and to the right, clear of the curves and the sweep points
            ax.text(p["na"] * 1.07, 100 * p["eta"] + 3, P.CASES[name]["title"], color=t["text"], fontsize=10,
                    ha="left", va="bottom")
            continue
        dy = -9 if name != "matched" else 4
        ax.text(p["na"], 100 * p["eta"] + dy, P.CASES[name]["title"], color=t["text"], fontsize=10,
                ha="center", va="bottom" if dy > 0 else "top")
    ax.set_xscale("log")
    ax.set_xlim(0.03, 1.2)
    ax.set_ylim(0, 103)
    ax.set_xticks([0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0], ["0.03", "0.05", "0.1", "0.2", "0.3", "0.5", "1"])
    ax.set_xlabel("beam NA  (λ / π w₀, the 1/e² half-angle of the focused beam)")
    ax.set_ylabel("power coupled into the core mode (%)")
    ax.grid(True, axis="y")
    if sweep:
        ax.legend(loc="upper left", fontsize=9.5)
    sec = ax.secondary_xaxis("top", functions=(lambda n: P.NA_MODE / np.maximum(np.asarray(n), 1e-9),
                                               lambda r: P.NA_MODE / np.maximum(np.asarray(r), 1e-9)))
    sec.set_xlabel("spot size on the facet / mode size", color=t["text2"], fontsize=10)
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    sec.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4]))
    sec.xaxis.set_minor_locator(NullLocator())
    sec.xaxis.set_minor_formatter(NullFormatter())
    sec.set_xticks([0.25, 0.5, 1, 2, 4], ["¼", "½", "1", "2", "4"])
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", default=os.path.join(ROOT, "data", "sweep", f"res{P.RESOLUTION}"))
    ap.add_argument("--main", default=os.path.join(ROOT, "data", "main", f"res{P.RESOLUTION}"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sweep = load_sweep(a.sweep)
    exact = exact_beam_pts(a.sweep, a.main)
    main_pts = {}
    for name in P.CASE_ORDER:
        try:  # the main runs only: a missing one is left out, never replaced by theory
            air = json.load(open(os.path.join(a.main, f"{name}_air", "meta.json")))
            fib = json.load(open(os.path.join(a.main, name, "meta.json")))
        except FileNotFoundError:
            continue
        last = fib["mode_planes"][-1]
        main_pts[name] = dict(na=P.CASES[name]["na"], eta=(last["P_fwd"] - last["P_bwd"]) / air["P_net_refl_plane"])
    ex = {round(p["w0"], 4): p["eta"] for p in exact}
    for p in sweep:
        print(f"NA {p['na']:.4f}  w0 {p['w0']:.3f}  eta {p['eta']:.4f}  R {p['R']:.4f}  "
              f"theory {fm.gaussian_overlap(p['w0']) * fm.fresnel_T(fm.lp01()['neff']):.4f}  "
              f"exact beam {ex.get(round(p['w0'], 4), float('nan')):.4f}")
    for theme in ("light", "dark"):
        fig_profiles(os.path.join(a.out, f"facet_profiles_{theme}.png"), theme)
        fig_efficiency(os.path.join(a.out, f"efficiency_vs_na_{theme}.png"), theme, sweep, main_pts, exact)


if __name__ == "__main__":
    main()
