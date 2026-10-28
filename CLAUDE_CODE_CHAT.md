# Chat with Claude about launching light into a single-mode fiber with too narrow an NA

Every prompt is recorded here as it was typed, with minor typos fixed, followed by
the response. The simulations, scripts, movies and documents in this directory are what
the responses produced. Times are US Pacific.

## Prompt 1

> When light entering a fiber isn't mode matched, I know that not all of it gets launched
> down the fiber. If the incoming mode is too big, I can imagine that some of it gets into
> the cladding and fails to propagate. But what happens to the "rest" of light when NA into
> fiber is too narrow? Make 3D E&M sims and show me movies. Eventually I might want to check
> this into github, so make a git repo here with a README.md that shows the video near the
> top of the page and an MIT LICENSE and everything else. Take your time to do it right.
> Install locally any electromagnetic solvers and renderers that you need.

## Prompt 2

Sent while the first response was being built:

> Keep going, but record the initial prompt and other prompts in CLAUDE_CODE_CHAT.md along
> with your responses, just like you do for the PCB designs in ~/pcbs.

## Prompt 3

Also sent while the first response was being built:

> Keep going, but could you explain the "Mellanox InfiniBand cards were causing UCX probing
> issues"? Do I have something configured incorrectly? If not, just keep going.

Answered in passing: nothing is misconfigured. The two ConnectX ports run in Ethernet (RoCE)
mode and the locked-memory limit is 8 MB. That is fine for normal use, but UCX, the transport
layer under the MPI library, probes every RDMA device it finds, tries to pin memory on it, and
prints an error each time. Within one machine MPI only needs shared memory, so the runs set
`UCX_TLS=self,sm,tcp` (in `scripts/env.sh`) and never touch the cards.

## Prompt 4

Also sent while the first response was being built:

> Keep going, but between all of the "too narrow" and "too wide" plots and pictures needs to
> be a "just right" plot or picture for comparison.

Folded in: every figure and movie shows the three launches in the order too narrow, matched
("just right"), too wide. That includes the facet profiles, which had shown only the
too-narrow split, the long-distance side view, and the resolution study, which now runs the
too-wide launch too.

## Prompt 5

Also sent while the first response was being built:

> The processors are still being used 100%. What are you waiting for? If most things are in,
> do a review of the code and your results and make sure things make sense.

Answered in passing: the processors are busy with the 3D runs themselves (three launches at
about 3 hours each, then the NA sweep and the resolution study), so nothing was waiting.
The review was done while they ran: a check of the numbers so far against theory, and an
independent read of the simulation and analysis code by a separate Claude agent. It found
no error in the physics (mode projection, normalization, symmetry, reflection, units). It
did find things that changed reported numbers or labels, all fixed:

- the facet reflection at 16 px/λ is about 12% low (3.0% against the 3.4% it converges to),
  so the resolution study now extrapolates to infinite resolution and the README quotes both;
- the 10 and 12 px/λ convergence runs predated a fix to where the monitor planes sit on the
  grid, so they are being rerun;
- the long-distance (BPM) model scaled the fiber's size but not the coating's index step,
  which made the coating strip light too fast; the coating index is now scaled too;
- the long-distance model started 20 λ in, after some steep light had left the 3D window, and
  started from a slightly different guided/unguided split than the 3D simulation; it now starts
  5 λ in from exactly the 3D split;
- the FDTD effective index printed in the results includes the grid's numerical dispersion;
  corrected for it, it is 1.46105 against 1.46107 from the exact fiber-mode equation;
- smaller label and accounting fixes: the movie said every row had the same brightness scale
  when each had its own (since then every steady-state picture, the side views, the end of the
  movies and the 3D render, shares one scale per watt of incident light, and the field rows say
  they are scaled separately), the field colours were described backwards, light that had only left
  the simulated window was called "escaped" or "radiated", the same four destinations of the
  light had different names in different figures, and the beam in air was drawn as uncoupled
  light. It is now white.

## Prompt 6

Also sent while the first response was being built:

> In the simulations, can we see enough to see that the light that doesn't make it into the
> core doesn't just go into the cladding, but that it dies away before making it to the other
> side of the simulated cladding? Or would that take a much longer simulation?

Answered in passing: not in the 3D simulation, and that is part of the answer. Nothing in the
cladding glass absorbs. The uncoupled light dies only where it crosses the cladding surface
into the higher-index coating. In 120 λ of fiber (0.6 mm of SMF-28) the too-narrow launch's
uncoupled light is still 57% of the incident power. It has spread from a radius of about 6 λ
to 12 λ (holding 90% of it) and has not reached the real cladding surface at 23 λ. The 3D
window has no cladding surface at all: the glass runs into the absorbing boundary at 14 λ, so
light that reaches it is counted as "gone sideways". The long-distance model, which has the
whole cladding and the coating, shows the stripping. 1% of the uncoupled light is gone by
165 λ, half by 1330 λ (0.6 cm of SMF-28) and three quarters by 3590 λ (1.7 cm). At 2.9 cm, 9%
of the incident power is still in the cladding, falling by 1/e every 2.7 cm. A 3D simulation
long and wide enough to show half of it stripped would cost about 500 times the 3-hour run
and several hundred GB of memory.

## Prompt 7

Also sent while the first response was being built:

> This is fine, as long as your writeup makes it clear **with pictures** what the scales for
> each of the simulations are. Maybe make a telescoping series of plots that visually expands
> or shrinks the previous one.

Folded in: `media/scales.png` (`scripts/fig_scales.py`) zooms out from the facet in four
steps, each panel drawn as a box inside the next with zoom lines between them:

1. the actual wavefronts over the first 16 λ;
2. the whole 3D simulation at its true aspect ratio, inside the real cladding it leaves out;
3. the long-distance model over 6000 λ (2.9 cm), with the coating;
4. a 1 m patch cord.

The README shows it at the start of the results, with a table of what each simulation covers.

## Prompt 8

Also sent while the first response was being built:

> When you're all done, remake all the plots. The "too wide" and "too narrow" labels seem
> backward. Try "too small" and "too big" for the beam itself as it touches the fiber tip.
> That naming is much more intuitive.

Folded in: every figure, movie and the 3D render now names the launches by the spot on the
fiber tip: "Spot too big", "Spot matched", "Spot too small", with the size relative to the
mode and the beam NA underneath (for example "3× the mode · beam NA 0.064 (too narrow)"). The
README keeps the original question as its title and says once, in section 1, that a launch NA
that is too narrow is a spot that is too big; after that it uses the spot names.

## Prompt 9

Also sent while the first response was being built:

> Then after you're totally done with the fiber, do a version of all three simulations where
> you couple into a metal microwave wave-guide tube instead of a fiber. The dimensions
> relative to the wavelength should be similar. The answer to my question, "Where does the
> light go that doesn't enter in the correct mode?" must be different for the thick metal
> wave guide. There is no cladding. So if the wave guide is truly single mode, and the metal
> is perfectly conducting, the only way not to go in is to reflect. Compare and contrast these
> situations.
>
> I don't know how you want to set up the incoming microwaves to be either matched, focused to
> a spot that's too small, or focused to a spot that's too big. You could just make these the
> initial conditions. I don't know how much you need to acknowledge that such a mismatch comes
> from a mismatched horn or parabolic dish.

To be done after the fiber work is finished; see Response 2.

## Response 1 2026-10-27 10:09 - 2026-10-28 04:05 Claude Opus 5.5

*Prompts 2 to 8 arrived while this was being built; each is folded in, as noted under it.
Prompt 9, the metal waveguide, is Response 2.*

**The answer.** The light that a too-narrow launch fails to couple goes forward, into the
cladding, and stays there for a surprisingly long way. A launch NA that is too narrow is a
focused spot that is too big for the fiber's mode. Only the part of the incoming field with the
mode's shape is guided (38% for a spot 3× too big). The rest (59%) is neither reflected nor
scattered at the facet: it enters the glass as a broad beam around the core, almost all of it
inside the core's acceptance cone, and travels beside the core without ever joining it. It leaks
out through the coating slowly. In SMF-28 terms half of it is gone after 6 mm, but a sixth is
still in the cladding after 3 cm. A spot 3× too small (NA too wide) loses about as much, but at
steep angles, and that light is gone within a few millimetres.

[![the turn-on movie](media/launch_movie.webp)](media/launch_movie.mp4)

### What was done

1. **Tools**, installed locally in two micromamba environments (`environment.yml`,
   `environment-render.yml`): Meep 1.35 (FDTD, with its MPB mode solver, MPI build) for the
   simulations, and PyVista on an off-screen VTK for the 3D render. They are separate because
   Meep's MPI build of HDF5 conflicts with VTK's.
2. **The fiber**: SMF-28's V number (2.36) in units of the wavelength, shrunk across by 1.76
   with the NA raised by the same factor. That keeps the mode's shape and every coupling number
   and makes the 3D volume 10× smaller (`fiberlaunch/params.py`, `fibermode.py`).
3. **3D FDTD** (`fiberlaunch/fdtd.py`): the three launches into 120 λ of fiber at 16 px/λ
   (about 3 h each on 16 cores), six more NAs into 10 λ of fiber, and a resolution study at
   10, 12, 16 and 20 px/λ. The beam is Meep's exact focused Gaussian, not a paraxial one.
4. **Analysis** (`fiberlaunch/analysis.py`): the field projected onto the fiber's HE₁₁ mode at
   every cross-section, the power budget, the angles of the uncoupled light, and an exact
   (non-paraxial) coupling estimate to check the simulation against.
5. **Long distances** (`fiberlaunch/bpm.py`): a beam-propagation model of the whole cladding
   and an absorbing coating, started from the 3D field 5 λ inside the facet and run to 6000 λ
   (2.9 cm of SMF-28).
6. **Media**: the turn-on movie (above), a 3D volume render, a fly-through of the
   cross-sections, the long-distance movie, and the figures in the README, among them the
   zoom-out of the scales (Prompt 7) and the facet profiles.
7. **The repository**: `README.md`, an MIT `LICENSE`, a `Makefile` that rebuilds everything,
   unit tests (`tests/`), and `.gitignore` for the 10 GB of raw data.

| | spot too big (NA 0.064) | matched (NA 0.19) | spot too small (NA 0.58) |
|---|---:|---:|---:|
| in the core mode, 3D FDTD | 38.2% | 96.4% | 37.7% |
| uncoupled, still near the core at 115 λ | 57.4% | 0.1% | 6.2% |
| uncoupled, gone sideways by 115 λ | 1.3% | 0.3% | 52.7% |
| reflected by the facet | 3.0% | 3.2% | 3.4% |
| still in the cladding after 2.9 cm (BPM) | 9.3% | 0.02% | 0.1% |
| exact focused-beam estimate of the coupling | 38.1% | 96.1% | 37.5% |
| coupling, extrapolated to an infinitely fine grid | 38.0% | 96.0% | 37.4% |

### Along the way

- **The review** asked for in Prompt 5 found no error in the physics but changed several
  numbers and labels; the list is under Prompt 5.
- **The 3D render** first came out too dim, and its opacity setting did nothing. VTK's GPU
  volume mapper ignores the opacity unit distance for colour-and-opacity volumes, and the CPU
  mapper that honours it draws visible steps. The render stays on the GPU mapper, with its
  opacity response measured and the opacity curve set to suit it, and jittered rays against
  banding.
- **The labels** changed to the spot on the fiber tip (Prompt 8), and every figure, movie and
  the 3D render were remade with them.
- **Resolution.** The resolution study finished last: from 16 to 20 px/λ the coupling moves
  by at most 0.15% of the incident power, and extrapolated it agrees with the exact
  focused-beam estimate to 0.2%. The facet's reflection converges slowly: 3.0 to 3.4% at
  16 px/λ, 3.4 to 3.8% extrapolated.

Not committed yet: the repository is initialized, and its first commit comes after Response 2.

## Response 2 2026-10-28 04:05 - 2026-10-28 09:07 Claude Opus 5.5

**The answer.** As you expected, in a single-mode, perfectly conducting guide the light that
is not coupled comes back, all of it. Guided plus reflected is 100% to within 0.01% in all
three simulations. What changes with the launch is where it turns back:

- **Spot too big** (a horn too small for its dish, which antenna engineers call *spillover*):
  70% of the power lands on the metal face around the horn's mouth and is reflected there.
  24% reaches the guide.
- **Spot too small** (a horn too big for its dish, *under-illumination*): 99% goes into
  the mouth, but as the horn's higher modes, because a small spot is a wide cone of directions.
  The horn narrows, each of those modes reaches its cut-off and turns back there, and
  75% of the power comes back out of the mouth. 24% reaches the guide.
- **Matched:** 80% reaches the guide. A smooth horn's mouth field is close to a Gaussian
  but not equal to one, so the matched beam, too, puts 16% into higher modes that turn back
  inside the horn.

So the coupling numbers are much like the fiber's (an overlap with the mode either way), but the
rest of the light goes the opposite way. The fiber's cladding is a continuum of forward-going
states that are not guided; the metal guide has only its modes, and they carry nothing below
cut-off, so the only state left is the incoming one reversed.

[![the metal-guide movie](media/horn_movie.webp)](media/horn_movie.mp4)

### What was done

1. **Why a horn.** A metal pipe is single mode only while it is less than 0.77 λ across, so it
   cannot be the size of the fiber's core; its mode is also smaller than any spot a beam can be
   focused to. A guide like that is fed through a horn. Here the horn is a smooth 8° cone cut into
   a thick metal block, 12.4 λ long, opening from the guide's radius (0.36 λ, only TE₁₁
   propagates) to 2.1 λ at the mouth. Its best Gaussian (1.60 λ waist) is within 0.1 point as
   good as the fiber's matched beam (1.66 λ), so the same three beams are matched, too big and too
   small for both. The README explains how a dish's f/D sets the spot, and how spillover and
   under-illumination map to "too big" and "too small".
2. **The beams.** Rather than starting them as initial conditions, each beam is radiated forward
   by a plane of equivalent electric and magnetic currents 4 λ in front of the mouth (the same
   exact focused Gaussian as in the fiber runs), and the run continues until it is steady. Every
   number is then a ratio of steady powers, and the movie still shows the beam arriving.
3. **The simulation.** The first version used a 3D Cartesian grid and could not converge: Meep
   cannot smooth a perfect conductor below the pixel, so the cone's wall became a staircase that
   scattered TE₁₁ into other modes (the matched coupling came out as 41, 55 and 69% at 16, 20
   and 28 points per wavelength). Since the horn, the guide and the beam share an axis, an
   x-polarized beam is exactly half angular order m = +1 and half m = −1. Meep's cylindrical mode
   solves the full 3D Maxwell equations for one angular order on the (r, z) plane, so the
   production runs are exact 3D at 2D cost: 120 points per wavelength, about 5 hours
   per launch on 5 cores (`fiberlaunch/horn.py`), plus a run of each beam with no metal for the
   incident power.
4. **Counting modes** (`fiberlaunch/horn_analysis.py`). A cone's modes are spherical waves
   about its apex, so the field is taken on spherical caps centred there, one per depth, and each
   mode's power heading in and heading back out comes from Lorentz reciprocity against the same
   mode going the other way. That separates the two directions exactly. Sent out from the guide,
   TE₁₁ comes back as 98.6 to 99.8% TE₁₁ on every cap.
5. **Checks.** The coupling also follows from a single run in reverse (TE₁₁ sent out of the
   guide, overlapped with each beam), done at 20, 40, 60, 80 and 120 points per wavelength and
   extrapolated to an infinitely fine grid (`scripts/horn_design.py`). It agrees with the direct
   runs, and the textbook aperture estimate (TE₁₁ filling the mouth on a spherical wavefront)
   agrees too.
6. **Media and README.** A movie of the three launches into the horn, the steady state on the
   plane of polarization, the power in each of the horn's modes against depth, and bars of where
   the power goes in the fiber and in the metal guide side by side; a new README section with all
   of it and a comparison table; `make horn` and `scripts/run_horn.sh` to rerun it.

| | spot too big (NA 0.064) | matched (NA 0.19) | spot too small (NA 0.58) |
|---|---:|---:|---:|
| **into the guide's mode (TE₁₁)** | **24.0%** | **80.2%** | **23.7%** |
| reflected by the metal face around the mouth | 70.3% | 4.0% | 0.8% |
| entered the horn, turned back inside it | 5.7% | 15.8% | 75.4% |
| coupling, one reverse run, extrapolated | 24.0% | 81.2% | 23.9% |
| aperture theory | 24.2% | 82.1% | 24.9% |
| *the fiber, for comparison: into the core mode* | *38.2%* | *96.4%* | *37.7%* |
| *the fiber: uncoupled, forward into the cladding* | *58.7%* | *0.4%* | *58.9%* |

### Along the way

- **The Cartesian staircase** (point 3) showed up only in the resolution series, as a coupling
  that kept climbing; that version is kept as `fiberlaunch/horn_cartesian.py`, and the README
  says why it was replaced.
- **Counting modes on flat planes** gave noisy, unphysical curves: inside a cone the wavefronts
  are curved, and a flat cut mixes the directions. Moving to caps about the apex fixed it.
- **The beam-in-air runs** first used the horn runs' long settling time; they were stopped and
  rerun with a time to suit an empty box.
- **The return loss.** The horn is well matched to the guide: TE₁₁ sent out of the guide has
  less than 0.1% of its power reflected back. So a network analyzer on the guide sees a good
  antenna, and a mis-focused beam shows up only as lost gain.
- **A dip that was the grid.** In a 40 px/λ trial the matched beam's TE₁₁ heading in seemed to
  lose 20 points halfway down the horn and win them back. At 120 px/λ it moves by only ±3
  points, near the cut-offs; the README explains what is left of it.

Committed as the repository's first commit, on `main`, not pushed: the code, the README, the
media (23 MB) and this chat log. The 13 GB of raw simulation output is left out by
`.gitignore`; `make sims` regenerates it.
