"""The incoming beams.
* incoming_beams_{dark,light}: the fiber's three beams from far away, the 1/e² radius against
  distance in front of the focus, with the stretch the simulations show (the source is 8 λ out).
  A beam converges over about its Rayleigh range z_R = π w₀² / λ, so the spot 3× too big
  (z_R = 78 λ) is still almost parallel over the last 8 λ, while the spot 3× too small
  (z_R ≈ 1 λ) is a cone.
* incoming_simulated.png: the nine beams the horn simulations launch (the three sets of
  params.HORN_SETS; the first is also the fiber's), as the FDTD makes them with no metal in
  the way: the field itself, and the 1/e² radius measured from it against Gaussian-beam theory.
  Drawn when the horns' air runs exist.

    python scripts/fig_incoming.py [--out media]
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
from fiberlaunch import beam  # noqa: E402
from fiberlaunch import params as P  # noqa: E402
from fiberlaunch import style  # noqa: E402
from fiberlaunch.horn_analysis import COMPS  # noqa: E402
from fig_theory import CASE_COLORS  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
Z_FAR = 100.0
SET_ROWS = dict(smooth="Focused on the mouth: the fiber's own three beams",
                smooth_focus="For the smooth horn, focused 3.7 λ inside it",
                corrugated="For the corrugated horn, focused 1.7 λ inside it")


def radius(z, w0):
    """Paraxial 1/e² radius at distance z from the waist."""
    return w0 * np.sqrt(1 + (z / (np.pi * w0 ** 2)) ** 2)


def fig(out, theme):
    t = style.use(theme)
    cc = CASE_COLORS[theme]
    z = np.linspace(-Z_FAR, 0, 2001)
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    ax.axvspan(P.Z_SOURCE_FAR, 0, color=t["glass"], lw=0, zorder=0)
    ax.text(0.5 * P.Z_SOURCE_FAR, 0.97, f"in the\nsimulations\n({-P.Z_SOURCE_FAR:.0f} λ)",
            transform=ax.get_xaxis_transform(), ha="center", va="top", color=t["text2"], fontsize=9)
    ax.grid(True, axis="y", which="major")
    for name in P.CASE_ORDER:
        c = P.CASES[name]
        w0 = c["w0"]
        zr = np.pi * w0 ** 2
        ax.plot(z, radius(z, w0), color=cc[name], lw=2, label=f"{c['title']} (w₀ = {w0:.2f} λ)")
        w8 = radius(P.Z_SOURCE_FAR, w0)
        ax.text(-Z_FAR + 1.5, radius(-Z_FAR, w0) * 1.12,
                f"{c['title']}: z_R = {zr:.2g} λ;  over the last {-P.Z_SOURCE_FAR:.0f} λ, {w8:.2f} → {w0:.2f} λ",
                color=t["text"], fontsize=9.5, ha="left", va="bottom")
    ax.set_yscale("log")
    ax.set_ylim(0.3, 120)
    ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50, 100])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"])
    ax.set_xlim(-Z_FAR, 0)
    ax.set_xlabel("z / λ   (the light travels to the right; focused at z = 0)")
    ax.set_ylabel("beam radius (1/e²) / λ")
    ax.set_title("The three incoming beams, from 100 λ out", loc="left", fontweight="bold")
    ax.legend(loc="lower left", fontsize=9.5, handlelength=1.8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def air_run(set_name, case, res, xlim):
    """A horn launch's beam alone: z, x (both signs) and Ex on the plane of polarization, the
    1/e² radius of |Ex|² along it, and the run's meta.  Ex there is the family's first profile."""
    d = os.path.join(ROOT, "data", "horn", set_name, f"res{res}", f"{case}_air")
    f = np.load(os.path.join(d, "dft_rz.npz"))
    r, z = f["r"], f["z"]
    ex = beam.from_meep([f[c].astype(np.complex128) for c in COMPS])[0]
    I = np.abs(ex) ** 2
    pk = np.argmax(I, axis=1)  # out from the brightest point, in case it is off the axis
    w = np.array([r[k + np.argmax(row[k:] < row[k] * np.exp(-2))] for k, row in zip(pk, I)])
    keep = r <= xlim
    x = np.concatenate([-r[keep][:0:-1], r[keep]])
    E = np.concatenate([ex[:, keep][:, :0:-1], ex[:, keep]], axis=1)
    return z, x, E, w, json.load(open(os.path.join(d, "meta.json")))


def fig_simulated(out, res=P.HORN_RESOLUTION, xlim=10.0, z_from=0.5, dz_dots=0.5):
    t = style.use("dark")
    sets = list(P.HORN_SETS)
    fig, axs = plt.subplots(len(sets), 3, figsize=(15.0, 11.0), sharex=True, sharey=True)
    for i, set_name in enumerate(sets):
        for j, case in enumerate(P.CASE_ORDER):
            ax = axs[i, j]
            z, x, E, w, m = air_run(set_name, case, res, xlim)
            zs, foc, w0 = m["z_source"], m["focus"], m["w0"]
            band = z >= zs
            zz, Ere = z[band], E.real[band]
            img = style.signed_compress(Ere, np.abs(Ere).max(), 0.6)
            ax.imshow(img, origin="lower", aspect="equal", cmap=style.field_cmap(), vmin=-1, vmax=1,
                      extent=(x[0], x[-1], zz[0], zz[-1]), interpolation="bilinear")
            # theory as a thin line, the simulation as dots on it (off the source plane itself)
            zt = np.linspace(zs, zz[-1], 400)
            wt = radius(zt - foc, w0)
            for sgn in (1, -1):
                ax.plot(sgn * wt, zt, color=t["text2"], lw=1.0)
            zd = np.arange(zs + z_from, zz[-1], dz_dots)
            k = [np.argmin(abs(z - v)) for v in zd]
            for sgn in (1, -1):
                ax.plot(sgn * w[k], z[k], "o", ms=4.0, mfc=t["text"], mec="#000000", mew=0.6)
            ax.axhline(0, color=t["muted"], lw=0.8, ls=":")
            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(zz[-1], zs)
            ws, wm = w[np.argmin(abs(z - (zs + z_from)))], w[np.argmin(abs(z))]
            if i == 0:
                ax.text(0, 1.30, P.CASES[case]["title"], transform=ax.transAxes, color=t["text"],
                        fontsize=13, fontweight="bold", va="bottom")
            ax.text(0, 1.02, f"w₀ {w0:.2f} λ, z_R {np.pi * w0 ** 2:.2g} λ  ·  radius {ws:.2f} → {wm:.2f} λ "
                    f"from z = {zs + z_from:g} to 0", transform=ax.transAxes, color=t["text2"], fontsize=9.5,
                    va="bottom")
            ax.set_xticks([-8, -4, 0, 4, 8])
            ax.set_yticks([-8, -4, 0])
            if j == 0:
                ax.set_ylabel("z / λ")
            if i == len(sets) - 1:
                ax.set_xlabel("x / λ")
            if j == 2:
                where = "on it" if foc == 0 else f"{foc:g} λ\nbeyond it"
                ax.text(1.03, zz[-1] / (zz[-1] - zs), f"z = 0: the mouth\n(or the facet)\nthe focus: {where}",
                        transform=ax.transAxes, color=t["text2"], fontsize=9, va="center")
        axs[i, 0].text(0, 1.13, SET_ROWS[set_name], transform=axs[i, 0].transAxes, color=t["text"],
                       fontsize=11.5, va="bottom", fontweight="bold")
    fig.text(0.05, 0.975, "The incoming beams as the simulations launch them (no metal or glass in the way)",
             color=t["text"], fontsize=14, fontweight="bold", va="top")
    fig.text(0.01, 0.012, "The electric field Eₓ at one instant on the plane of polarization (red / blue = + / −), "
             "from the source plane 8 λ out to just past z = 0, where the mouth or the facet goes.\n"
             "Dots: the beam's 1/e² radius, measured from the simulated |Eₓ|².  Gray line: Gaussian-beam theory, "
             "w(z) = w₀ √(1 + (z − z_focus)² / z_R²).\nEvery beam converges.  The spot 3× too big does it over "
             "its Rayleigh range z_R, 46–78 λ, so in these 8 λ it barely narrows.",
             color=t["text2"], fontsize=9.5, va="bottom")
    fig.subplots_adjust(left=0.05, right=0.90, top=0.85, bottom=0.11, wspace=0.06, hspace=0.50)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "media"))
    a = ap.parse_args()
    for theme in ("dark", "light"):
        fig(os.path.join(a.out, f"incoming_beams_{theme}.png"), theme)
    if all(os.path.exists(os.path.join(ROOT, "data", "horn", s, f"res{P.HORN_RESOLUTION}", f"{c}_air", "dft_rz.npz"))
           for s in P.HORN_SETS for c in P.CASE_ORDER):
        fig_simulated(os.path.join(a.out, "incoming_simulated.png"))


if __name__ == "__main__":
    main()
