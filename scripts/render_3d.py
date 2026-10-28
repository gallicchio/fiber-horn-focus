"""3D volume rendering of the three launches (render environment: PyVista, off-screen).

Each launch's steady-state intensity is drawn as a glowing volume, coloured by which part
of the light dominates at each point: the guided mode blue, everything else orange, the beam
still in air white (the same split as the 2D movies).  The half of each volume
facing the camera is cut away, so the cut face shows the inside and the back half shows
the 3D shape.  The camera moves from in front of the facet round to the side.

    python scripts/render_3d.py [--vol data/volumes] [--out media/launch_3d.mp4]
"""

import argparse
import os
import sys

import imageio.v2 as imageio
import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fiberlaunch import params as P  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
GUIDED = np.array([0x5A, 0xA9, 0xFF]) / 255
REST = np.array([0xFF, 0x8A, 0x3D]) / 255
AIR = np.array([0xD9, 0xD7, 0xD0]) / 255  # the beam before the facet, not yet split


def rgba_volume(d, decades=3.0, gamma=0.8, amax=1.0, glow=0.55, cut=True):
    g = d["guided"].astype(np.float32)  # (nz, nx, ny), 1 = the shared brightness scale
    r = d["rest"].astype(np.float32)
    x, z = d["x"], d["z"]

    def lg(v):
        return np.clip((np.log10(np.maximum(v, 1e-12)) + decades) / decades, 0, 1)

    rest_col = np.where((z < 0)[:, None, None, None], AIR, REST)  # in air nothing is guided yet
    # hue: which part carries the power here (blue where the guided mode dominates); brightness
    # and opacity: the total, on a log scale.  Colour brighter than opacity (glow < 1), so faint
    # light still reads as its hue.
    tot = g + r
    f = (g / np.maximum(tot, 1e-30))[..., None]
    v = lg(tot)
    rgb = (f * GUIDED + (1 - f) * rest_col) * v[..., None] ** glow
    a = amax * v ** gamma
    if cut:  # remove the half y > 0 (towards the camera)
        a[:, :, x > 0] = 0
    rgba = (np.concatenate([rgb, a[..., None]], -1) * 255).astype(np.uint8)
    grid = pv.ImageData(dimensions=(len(x), len(x), len(z)),
                        spacing=(x[1] - x[0], x[1] - x[0], z[1] - z[0]), origin=(x[0], x[0], z[0]))
    grid.point_data["rgba"] = np.transpose(rgba, (0, 2, 1, 3)).reshape(-1, 4)
    return grid


def glass(pl, x0, L, hw):
    """Facet and core outlines (thin lines) for one fiber centred at x = x0."""
    col = "#8f8e88"
    # the facet: a square outline at z = 0, half of it (the cut keeps y <= 0)
    pts = np.array([[x0 - hw, 0, 0], [x0 + hw, 0, 0], [x0 + hw, -hw, 0], [x0 - hw, -hw, 0], [x0 - hw, 0, 0]])
    pl.add_mesh(pv.lines_from_points(pts), color=col, line_width=1.5, opacity=0.8)
    pl.add_mesh(pv.Plane(center=(x0, -hw / 2, 0), direction=(0, 0, 1), i_size=2 * hw, j_size=hw),
                color="#c3c2b7", opacity=0.06)
    # the core: its outline on the cut plane, and its back half as a faint shell
    for s in (-1, 1):
        pl.add_mesh(pv.Line((x0 + s * P.CORE_RADIUS, 0, 0), (x0 + s * P.CORE_RADIUS, 0, L)),
                    color=col, line_width=1.0, opacity=0.6)
    th = np.linspace(np.pi, 2 * np.pi, 33)
    ring = np.stack([x0 + P.CORE_RADIUS * np.cos(th), P.CORE_RADIUS * np.sin(th), 0 * th], -1)
    pl.add_mesh(pv.lines_from_points(ring), color=col, line_width=1.2, opacity=0.8)


def camera_path(n, L, height):
    """From in front of the facets round to the side.  x is 'up' on screen (the fibers are
    stacked along x); the camera stays on the +y side, where the volumes are cut open.
    The distance keeps the whole stack (height along x, L along z) in the 30-degree view."""
    f = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, n))  # ease in and out
    out = []
    for s in f:
        phi = np.deg2rad(52 - 44 * s)  # from the +y axis towards -z (in front of the facet)
        el = np.deg2rad(12 - 6 * s)  # tilt towards +x
        dist = 2.3 * max(height, 0.55 * L) * (1.0 - 0.06 * s)
        focus = np.array([0.0, 0.0, 0.3 * L + 0.2 * L * s])
        d = np.array([np.sin(el), np.cos(el) * np.cos(phi), -np.cos(el) * np.sin(phi)])
        out.append((tuple(focus + dist * d), tuple(focus), (1, 0, 0)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vol", default=os.path.join(ROOT, "data", "volumes"))
    ap.add_argument("--out", default=os.path.join(ROOT, "media", "launch_3d.mp4"))
    ap.add_argument("--still", default=os.path.join(ROOT, "media", "launch_3d.png"))
    ap.add_argument("--cases", nargs="*", default=P.CASE_ORDER)
    ap.add_argument("--frames", type=int, default=240)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--size", type=int, nargs=2, default=(1920, 1080))
    ap.add_argument("--spacing", type=float, default=30.0)
    ap.add_argument("--amax", type=float, default=1.0, help="opacity at the brightest point")
    ap.add_argument("--gamma", type=float, default=0.8, help="opacity = amax * (log brightness)**gamma")
    ap.add_argument("--preview", action="store_true")
    a = ap.parse_args()

    pl = pv.Plotter(off_screen=True, window_size=a.size)
    pl.set_background("#000000")
    offs = {name: (1 - i) * a.spacing for i, name in enumerate(a.cases)}
    L = 0
    for name in a.cases:
        d = np.load(os.path.join(a.vol, f"{name}.npz"))
        grid = rgba_volume(d, amax=a.amax, gamma=a.gamma).translate((offs[name], 0, 0))
        # The GPU ray caster: smooth, but it ignores the opacity unit distance for RGBA data and
        # its opacity grows faster than linearly with alpha, hence the low gamma.  (The CPU
        # fixed-point caster honours the unit distance but draws steps along z.)
        vol = pl.add_volume(grid, scalars="rgba", mapper="gpu", blending="composite")
        vol.prop.interpolation_type = "linear"
        vol.mapper.SetUseJittering(True)  # no wood-grain banding along the rays
        L = float(d["z"][-1])
        hw = float(d["x"][-1])
        glass(pl, offs[name], L, hw)
        c = P.CASES[name]
        pl.add_point_labels([(offs[name] + hw + 1.5, 0, 0)], [f"{c['title']}  ({c['na_label']})"],
                            font_size=22, text_color="#f0efec", shape=None, show_points=False,
                            always_visible=True)
    pl.add_text("Gaussian beam into a single-mode fiber (3D FDTD).  White: the beam in air.  "
                "Blue: light in the core's guided mode.  Orange: the rest.", position="upper_left", font_size=12,
                color="#f0efec")
    path = camera_path(a.frames, L, a.spacing * len(a.cases))
    if a.preview:  # first, middle and last frame only
        path = [path[0], path[len(path) // 2], path[-1]]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    w = imageio.get_writer(a.out, fps=a.fps, codec="libx264", quality=None, macro_block_size=8,
                           pixelformat="yuv420p", ffmpeg_params=["-crf", "18", "-preset", "slow"])
    pl.show(auto_close=False)  # off screen, the image only updates after show() + render()
    for i, cam in enumerate(path):
        pl.camera_position = cam
        pl.render()
        img = pl.screenshot(return_img=True)
        w.append_data(img)
        if i == len(path) - 1:
            imageio.imwrite(a.still, img)
    w.close()
    pl.close()
    print("wrote", a.out)


if __name__ == "__main__":
    main()
