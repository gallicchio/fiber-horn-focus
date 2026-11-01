"""The movie frames between the saved ones (no simulation data needed)."""

import numpy as np

from fiberlaunch import frames

DT = 0.34375  # the fiber runs' frame spacing, in periods


def field(t, phase):
    """A point's field as Meep's continuous source switches on: (1 + tanh(t / 5 - 3)) / 2."""
    return 0.5 * (1 + np.tanh(t / 5 - 3)) * np.cos(2 * np.pi * t - phase)


def test_in_between_frames_follow_the_switching_on_field():
    times = np.arange(0, 60, DT)
    phases = np.linspace(0, 2 * np.pi, 7)
    saved = field(times[:, None], phases)
    for t in np.arange(1.0, 58.0, 1 / 15):
        assert np.allclose(frames.at(saved, times, t), field(t, phases), atol=2e-4)

