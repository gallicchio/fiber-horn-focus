"""Colors and helpers shared by every figure and movie.

The four destinations of the light keep one color each everywhere (bars, curves, and the
field images): guided = blue, uncoupled light near the core = orange, uncoupled light gone
sideways (out of the simulated window, or stripped by the coating) = aqua, reflected =
yellow.  In the field images the beam still in air, before the facet, glows neutral white;
where two parts overlap, the hue leans to whichever carries more light.  The values are the
first four slots of a categorical palette validated for color-vision deficiency on both
surfaces (light and dark steps differ).
"""

import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb

THEMES = {
    "light": dict(surface="#fcfcfb", text="#0b0b0b", text2="#52514e", muted="#8a8984",
                  grid="#e4e3df", glass="#e9e8e4",
                  guided="#2a78d6", cladding="#eb6834", lost="#1baf7a", reflected="#eda100"),
    "dark": dict(surface="#1a1a19", text="#f0efec", text2="#c3c2b7", muted="#8f8e88",
                 grid="#33332f", glass="#2a2a28",
                 guided="#3987e5", cladding="#d95926", lost="#199e70", reflected="#c98500"),
}

# The same four names in every figure and movie.  In the FDTD, "lost" is light that has
# left the simulated window (beyond 14 lambda of the axis, still inside a real fiber's
# cladding); in the long-distance BPM it is light absorbed in the coating.
LABELS = dict(
    guided="in the core mode",
    cladding="uncoupled, near the core",
    lost="uncoupled, gone sideways",
    reflected="reflected by the facet",
)
LOST_BPM = "uncoupled, stripped by the coating"

# Glow colors for the two-channel field images on black: brighter steps of the same hues.
GLOW = dict(guided="#5aa9ff", cladding="#ff8a3d", air="#d9d7d0", reflected="#ffc61a")


def use(theme="dark", font_scale=1.0):
    t = THEMES[theme]
    mpl.rcParams.update({
        "figure.facecolor": t["surface"], "axes.facecolor": t["surface"],
        "savefig.facecolor": t["surface"], "text.color": t["text"],
        "axes.labelcolor": t["text2"], "axes.edgecolor": t["muted"],
        "xtick.color": t["text2"], "ytick.color": t["text2"],
        "axes.grid": False, "grid.color": t["grid"], "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 11 * font_scale, "axes.titlesize": 12 * font_scale,
        "axes.labelsize": 11 * font_scale, "legend.frameon": False,
        "lines.linewidth": 2.0, "font.family": "DejaVu Sans",
    })
    return t


def field_cmap():
    """Diverging map for the instantaneous field: blue - black - red.  Berlin, with its
    dark-brown centre pulled to black so that zero field matches the black of the
    intensity images."""
    x = np.linspace(0, 1, 21)
    c = mpl.colormaps["berlin"](x)[:, :3]
    c[10] = 0.0
    return LinearSegmentedColormap.from_list("berlin_black", list(zip(x, c)))


def glow_cmap(color, name="glow"):
    """Sequential one-hue ramp from black to the hue to a whitened tip."""
    c = np.array(to_rgb(color))
    return LinearSegmentedColormap.from_list(name, [(0, 0, 0), tuple(0.55 * c), tuple(c),
                                                    tuple(1 - 0.35 * (1 - c))])


def log_norm(I, ref, decades):
    """Intensity to [0, 1]: 0 at ref * 10**-decades, 1 at ref."""
    with np.errstate(divide="ignore"):
        v = (np.log10(np.maximum(I, 1e-300) / ref) + decades) / decades
    return np.clip(v, 0, 1)


def two_channel(Ig, Ir, ref, decades=3.0, Ia=None):
    """Composite image (..., 3) on black.  Brightness: the total intensity on a log scale
    (the glow ramp).  Hue: which part carries the light there, mixed in proportion to the
    intensities: guided blue, the rest orange, and (if given) light still in air white.
    Overlapping parts mix toward the dominant one rather than adding up to white."""
    parts = [(Ig, GLOW["guided"]), (Ir, GLOW["cladding"])] + ([] if Ia is None else [(Ia, GLOW["air"])])
    return composite(parts, ref, decades)


def composite(parts, ref, decades=3.0):
    """two_channel for any list of (intensity, glow color) parts."""
    tot = sum(np.asarray(I, float) for I, _ in parts)
    c = sum(np.asarray(I, float)[..., None] * to_rgb(col) for I, col in parts) / np.maximum(tot, 1e-300)[..., None]
    v = log_norm(tot, ref, decades)[..., None]
    # the glow ramp of glow_cmap, per pixel: black, 0.55 c, c, whitened c at v = 0, 1/3, 2/3, 1
    nodes = np.stack([0 * c, 0.55 * c, c, 1 - 0.35 * (1 - c)])
    k = np.minimum((3 * v).astype(int), 2)
    w = 3 * v - k
    lo, hi = np.take_along_axis(nodes, k[None], 0)[0], np.take_along_axis(nodes, k[None] + 1, 0)[0]
    return (1 - w) * lo + w * hi


def signed_compress(E, emax, gamma=0.5):
    """Field to [-1, 1] with a power-law stretch so weak radiation stays visible."""
    return np.sign(E) * np.minimum(np.abs(E) / emax, 1) ** gamma
