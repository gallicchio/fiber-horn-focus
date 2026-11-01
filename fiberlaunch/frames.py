"""Movie frames between the saved ones.

The runs save the field every 0.34 periods, about three frames a period: too few to watch.  A
standing wave blinks, and a travelling wave can seem to run backwards, like a wagon wheel in a
film.  But the source has one frequency and switches on over several periods, so at every point
the field is E(t) = Re[A(t) exp(-2πit)] (t in periods) with A changing slowly.  Fitting a
quadratic A(t) to the ten saved frames nearest t gives the field at any t in between.
Predicting each saved frame this way from the others, with it left out, misses it by 0.12% of
the steady field's peak at worst, at any time, in the fiber and the horn runs alike: about the
float16 rounding of the saved frames themselves.
"""

import numpy as np


def weights(times, t, n=10, order=2):
    """The first of the n saved frames nearest t, and weights g with E(t) = Σ g E[k0:k0 + n]:
    the least-squares fit of Re[(a0 + a1 s + a2 s² ...) exp(-2πis)], s = t_k - t, at s = 0."""
    k0 = int(np.clip(np.searchsorted(times, t) - n // 2, 0, len(times) - n))
    s = times[k0:k0 + n] - t
    M = np.stack([s ** p * f(2 * np.pi * s) for p in range(order + 1) for f in (np.cos, np.sin)], axis=1)
    return k0, np.linalg.pinv(M)[0]


def at(frames, times, t, **kw):
    """The field at time t (periods) from the saved frames (frames[k] saved at times[k])."""
    k0, g = weights(times, t, **kw)
    return np.tensordot(g, np.asarray(frames[k0:k0 + len(g)], np.float32), 1)
