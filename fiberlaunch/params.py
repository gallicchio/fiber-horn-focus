"""Single source of truth for the fiber, the launch conditions and the simulation box.

All lengths are in units of the vacuum wavelength (lambda = 1).  Meep uses the same
unit, so the source frequency is 1 and one optical period is 1 time unit.  To get
microns, multiply by the wavelength in microns (e.g. x1.55 for 1550 nm).

The fiber is a step-index single-mode fiber with V = 2.36, almost exactly the V number
of SMF-28 at 1550 nm (V ~ 2.3), so its mode is the same shape relative to its core.
Its NA (0.25) is higher than SMF-28's (0.14) to keep the 3D FDTD box affordable: the
physics is identical, with transverse sizes scaled by 0.14/0.25 and angles (NAs) scaled
by 0.25/0.14.
"""

import math

WAVELENGTH = 1.0
FREQ = 1.0 / WAVELENGTH

# ---- fiber ----------------------------------------------------------------------------
N_CLAD = 1.45
NA_FIBER = 0.25
N_CORE = math.sqrt(N_CLAD**2 + NA_FIBER**2)  # 1.4714
CORE_RADIUS = 1.5
V_NUMBER = 2 * math.pi * CORE_RADIUS * NA_FIBER / WAVELENGTH  # 2.356 < 2.405: single mode

# Gaussian waist that best matches the LP01 mode (fibermode.optimal_waist(), 99.6 %).
W_MATCHED = 1.6585
# Mode field "NA": 1/e^2 far-field half-angle of the matched Gaussian, lambda/(pi w).
NA_MODE = WAVELENGTH / (math.pi * W_MATCHED)  # 0.192


def beam_na(w0):
    """Paraxial 1/e^2 divergence half-angle of a Gaussian beam with waist w0 (in air)."""
    return WAVELENGTH / (math.pi * w0)


def waist_for_na(na):
    return WAVELENGTH / (math.pi * na)


# ---- launch conditions ----------------------------------------------------------------
# An x-polarized Gaussian beam, focused (waist) exactly on the fiber facet, on axis.
# "narrow" and "wide" miss the mode by the same factor of 3 in waist, so they lose the
# same ~60 % of the light, but they lose it in completely different ways.  The keys name the
# launch NA (the question asked about); the labels name the spot on the facet, which is
# what decides the coupling: too narrow an NA is too big a spot.
CASES = {
    "narrow": dict(
        w0=3 * W_MATCHED,
        title="Spot too big",
        subtitle="3× the mode · beam NA 0.064 (too narrow)",
        na_label="beam NA 0.064",
    ),
    "matched": dict(
        w0=W_MATCHED,
        title="Spot matched",
        subtitle="the size of the mode · beam NA 0.19",
        na_label="beam NA 0.19",
    ),
    "wide": dict(
        w0=W_MATCHED / 3,
        title="Spot too small",
        subtitle="⅓ of the mode · beam NA 0.58 (too wide)",
        na_label="beam NA 0.58",
    ),
}
for _c in CASES.values():
    _c["na"] = beam_na(_c["w0"])
CASE_ORDER = ["narrow", "matched", "wide"]

# ---- simulation box (physical z: facet at z = 0, fiber along +z) -----------------------
RESOLUTION = 16  # pixels per wavelength (10.9 per wavelength inside the core)
DPML = 2.0  # PML thickness on the four sides and at the source end
DPML_END = 3.0  # PML thickness where the fiber leaves the box
HALF_WIDTH = 14.0  # half-width of the region inside the PML, x and y
Z_SOURCE = -4.0  # source plane (in air); the main runs and the horn launches use Z_SOURCE_FAR
Z_SOURCE_FAR = -8.0  # twice as far out: the pictures show twice as much of the converging beam
Z_REFL = -2.0  # flux plane in air: incident - reflected power
GAP_BEHIND_SOURCE = 1.0
FIBER_LENGTH = 120.0  # length of fiber inside the box (PML not included)

TURN_ON = 5.0  # ContinuousSource ramp width, in periods
DFT_PERIODS = 4  # steady-state DFTs are accumulated over this many whole periods

# planes (physical z) where the guided-mode content is measured with MPB
MODE_PLANES = [3.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 115.0]
MODE_HALF_WIDTH = 6.0  # MPB window around the core (HE11 power beyond it is ~1e-6)
XY_PLANE_DZ = 0.5  # spacing of the stack of transverse DFT planes

# ---- the metal waveguide (fiberlaunch/horn.py) ---------------------------------------
# A perfectly conducting circular guide fed by a smooth conical horn in the face of a thick
# metal block.  The throat is strictly single mode: TE11 propagates above radius 0.293 λ,
# TM01 only above 0.383 λ.  The mouth is about the size of the fiber's core, so that the
# matched beam of the fiber runs (NA 0.19) is also the best Gaussian for the horn.
THROAT_RADIUS = 0.36
MOUTH_RADIUS = 2.1
FLARE_DEG = 8.0  # half-angle of the cone
HORN_LENGTH = (MOUTH_RADIUS - THROAT_RADIUS) / math.tan(math.radians(FLARE_DEG))
THROAT_LENGTH = 6.0  # straight single-mode guide after the horn, before the absorber
HORN_RESOLUTION = 120  # pixels per wavelength (43 across the throat's radius); see scripts/horn_design.py
TE11_ROOT = 1.8412  # first zero of J1'
HORN_HALF_WIDTH = 18.0  # the horn launches: wide enough for the too-small beams converging from Z_SOURCE_FAR

# The corrugated horn (horn.py --wall corrugated): the same cone, its wall cut into slots
# perpendicular to the axis.  A slot a quarter wave deep (radially) is open at its mouth, so the
# wall holds neither E nor H along it and the mode becomes the balanced hybrid HE11.  Near the
# throat the slots start half a wave deep (a smooth wall, for TE11) and become a quarter wave deep
# over CORR_CONVERTER slots.  Depths follow the radial-line resonance (horn.slot_depth).
CORR_PITCH = 0.2  # slot + fin along the axis (λ/5)
CORR_SLOT = 0.15  # slot width (3/4 of the pitch)
CORR_CONVERTER = 6  # slots in the mode converter
CORR_FIRST_DEPTH = 1.0  # the slot nearest the throat: this × its half-wave depth (1: like a smooth wall)
CORR_DEPTH_SCALE = 1.0  # the other slots: this × their quarter-wave depth

# The three sets of horn launches (scripts/run_horn.sh).  Each is the three beams of CASES, too big,
# matched and too small by the same factors, around that set's best waist, all three focused at the
# depth where that best beam is focused (scripts/horn_design.py, from the reverse runs at 120 px/λ).
HORN_SETS = {
    # the fiber's own three beams, focused on the mouth
    "smooth": dict(wall="smooth", w_best=W_MATCHED, focus=0.0),
    # the smooth horn's best beam is smaller and focused inside: its mouth's wavefront is curved
    "smooth_focus": dict(wall="smooth", w_best=1.40, focus=3.7),
    # the corrugated horn: its mouth field is HE11, 98 % Gaussian, on a wavefront centred 1.7 λ inside
    "corrugated": dict(wall="corrugated", w_best=1.28, focus=1.7),
}
CASE_SCALE = {name: c["w0"] / W_MATCHED for name, c in CASES.items()}


def horn_waist(set_name, case):
    return HORN_SETS[set_name]["w_best"] * CASE_SCALE[case]


def te11_beta(radius):
    """Propagation constant (in units of k) of TE11 in a straight circular guide."""
    kc = TE11_ROOT / (2 * math.pi * radius)
    return math.sqrt(max(0.0, 1 - kc * kc))

