# Where the light goes when a launch into a single-mode guide is not matched: the full study

The short version, with the movies, is in the [README](README.md). This page is the full study:
what each simulation covers, how it was computed, all the numbers, the checks, and more figures and
movies.

<p align="center">
  <a href="media/launch_movie.mp4"><img src="media/launch_movie.webp" width="100%"
     alt="3D FDTD movie: a Gaussian beam switching on and filling a single-mode fiber for three spot sizes, then the steady-state split into guided (blue) and uncoupled (orange) light"></a>
  <br><sub>Full 3D Maxwell (FDTD) simulations of a focused Gaussian beam entering a cleaved single-mode
  fiber, for a spot on the fiber tip that is too big (launch NA too narrow), matched, and too small
  (launch NA too wide); the light comes down from the top. The beam switches on, fills the first
  40 λ of the fiber and settles, shown at 2 periods a second and at its true shape (x and z to
  one scale); then the picture changes to where the power ends up: <b>white</b> = the beam
  still in air, <b>blue</b> = light in the fiber's guided mode, <b>orange</b> = everything else
  (log scale, the same for all three). <a href="media/launch_movie.mp4">Full-resolution MP4</a>.</sub>
</p>

**Short answer.** It goes forward, into the cladding, and stays there for a surprisingly long
way. A launch NA that is too narrow is the same thing as a focused spot that is too big for the
fiber's mode. Only the part of the incoming field that has the mode's shape is guided: 38% for a
spot 3× too big. The rest (59%) is neither reflected nor scattered at the facet. It goes into the
glass as a broad beam around the core, almost all of it inside the core's acceptance cone, and
travels alongside the core without ever joining it. It leaks out through the fiber's coating only
slowly: in SMF-28 terms, half of it is gone after 6 mm, but a sixth of it is still in the cladding
after 3 cm, where a power meter at the end of a short pigtail will count it. A spot 3× too
*small* (launch NA too wide) loses about as much, but that light leaves at steep angles and is
gone within a few millimetres.

**In a single-mode metal waveguide** (fed through a horn, the way microwave guides are) the
answer is different, as it has to be: there is no cladding, so everything that does not end up
in the guide's one mode is reflected. A spot that is too big is mostly reflected by the metal
face around the horn's mouth. A spot that is too small goes into the mouth, but as the horn's
higher modes, and each of them turns back where the horn becomes too narrow for it
([below](#the-same-beams-into-a-single-mode-metal-waveguide)).

---

## Contents

- [The question](#the-question)
- [What the simulations show](#what-the-simulations-show)
- [The same beams into a single-mode metal waveguide](#the-same-beams-into-a-single-mode-metal-waveguide)
- [The numbers](#the-numbers)
- [More movies](#more-movies)
- [How it was computed](#how-it-was-computed)
- [How much to trust it](#how-much-to-trust-it)
- [Reproducing it](#reproducing-it)
- [Files](#files)

## The question

> When light entering a fiber isn't mode matched, not all of it gets launched down the fiber. If the
> incoming mode is too big, some of it gets into the cladding and fails to propagate. But what happens
> to the "rest" of the light when the NA into the fiber is too narrow?

## What the simulations show

### 1. "NA too narrow" and "spot too big" are the same thing

A focused Gaussian beam's waist and its divergence are locked together: w₀ = λ / (π · NA). Halve the
NA and the spot on the facet doubles. So a launch whose NA is too narrow puts a spot on the facet
that is bigger than the fiber's mode, which is exactly the "incoming mode is too big" case.

The fiber here (the same V number as SMF-28, see [below](#how-it-was-computed)) has a mode field
radius of 1.66 λ. The beam that matches it has NA 0.19, which is *less* than the fiber's NA of 0.25:
the NA on a fiber's data sheet is the acceptance cone of the core, not the divergence of its mode.
The figures and movies name each launch by its spot on the fiber tip: "spot too big" is the
launch whose NA is too narrow (NA 0.064, 3× narrower than matched, a spot 3× the mode); "spot too
small" is the launch whose NA is too wide (NA 0.58, 3× wider, a spot ⅓ of the mode).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/facet_profiles_dark.png">
  <img src="media/facet_profiles_light.png" alt="The three spots on the facet (too big, matched, too small), each split into the part shaped like the fiber mode (guided) and the rest">
</picture>

Only the part of the incoming field that has the shape of the mode is guided (blue). "The rest" is
whatever is left over: the incident field minus the guided part (orange). For the matched launch
the two shapes agree and the rest is almost nothing (0.4%). The too-big spot is broader and
flatter than the mode, so its rest has a phase-inverted dip inside the core and two broad shoulders
out in the cladding. The too-small spot is a narrow spike, so its rest is the opposite: a sharp peak
on the axis with shallow inverted wings. The power lost is about the same for both (61% and 64%);
what differs is where the rest goes, and that is what the 3D simulations show.

**All three are converging beams, but at very different rates.** Each beam narrows to its waist
over about its Rayleigh range, z_R = π w₀² / λ. That is 78 λ for the spot too big, 8.6 λ for the
matched spot and 1 λ for the spot too small. The simulations start the beam 8 λ in front of the
facet. Over those 8 λ the too-small beam closes from 4.6 λ in radius to 0.55 λ, a cone you can see.
The matched beam closes from 2.3 λ to 1.7 λ, and the too-big beam from 5.00 λ to 4.98 λ, so in the
pictures it looks parallel. It is converging, with its focus on the facet, but slowly: its radius
is 1.4 × its waist only 80 λ out.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/incoming_beams_dark.png">
  <img src="media/incoming_beams_light.png" alt="The radius of the three incoming beams against distance in front of the focus, out to 100 wavelengths; the simulations show the last 8">
</picture>

The simulated beams are these converging Gaussians. Below, each beam is launched with nothing in
its way (top row: the fiber's three beams; the other rows: the beams for the metal guide's
horns, [further down](#three-sets-of-launches-which-is-like-the-fiber)). The radius measured
from the simulated field (dots) follows Gaussian-beam theory (line) for every beam, and the
wavefronts curve around the focus. For the too-big spot neither shows much over these 8 λ: its
radius shrinks by 0.4%, and its wavefronts are flat to 0.02 λ across the beam.

![The nine incoming beams as the FDTD launches them, with the measured 1/e² radius on Gaussian-beam theory](media/incoming_simulated.png)

### What each simulation covers

Two simulations work together, and they cover very different scales. The 3D FDTD solves
Maxwell's equations in full, but only for the first 120 λ of fiber (0.6 mm of SMF-28), and only in a
28 λ window around the axis that stops short of the cladding surface. The long-distance model
(BPM) picks up the 3D field 5 λ inside the facet and carries it 6000 λ (2.9 cm of SMF-28) through
the whole cladding and its coating. Even that is the first 3% of a 1 m patch cord. The figure
zooms out from the facet in four steps, each panel drawn as a box inside the next:

![The scales: the facet, the whole 3D simulation, the long-distance model, a 1 m patch cord](media/scales.png)

| | 3D FDTD (Meep) | long-distance model (BPM) |
|---|---|---|
| length of fiber | 120 λ ≈ 0.6 mm of SMF-28 | 6000 λ ≈ 2.9 cm of SMF-28 |
| cross-section | 28 λ × 28 λ ≈ 77 µm square, inside the 125 µm cladding | the whole 125 µm cladding and an absorbing coating |
| physics | full vector Maxwell: the facet, its reflection, light at every angle | scalar waves going forward, exact for every angle |
| what it answers | how the light splits at the facet, and where the uncoupled part starts out | where the uncoupled light ends up |
| cost per launch | about 3 hours on 16 cores | 2 to 13 minutes |

### 2. The rest keeps going forward, through the cladding, beside the core

The side view through the axis, once everything has settled. In all three launches the light that
is not in the mode keeps going forward. Almost nothing comes back: the 3% reflection is the plain
Fresnel reflection of the glass, about the same for every launch. What differs is the direction
the rest takes.

- **Spot too big (NA too narrow):** the rest leaves the facet as a broad beam, about 15 λ across, that widens only
  slowly and runs almost parallel to the core. 115 λ in, 57% of the incident power is still inside
  the 28 λ window around the core, and only 1.3% has left it.
- **Spot too small (NA too wide):** the rest leaves the facet as a spike that diffracts at steep angles, up to about 30°
  in the glass. Most of it crosses the edge of the 28 λ window between 25 and 100 λ: by 115 λ, 53% of
  the incident power has gone sideways and only 6% is still near the core.

The core carries the same 38% in both. The edge of the 3D window is not the edge of the fiber:
light that leaves the window is still in the cladding glass (see
[What each simulation covers](#what-each-simulation-covers)).

![Steady state, side view of the three launches: spot too big, matched, too small](media/steady_side_view.png)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/budget_dark.png">
  <img src="media/budget_light.png" alt="Where the incident power is 115 wavelengths into the fiber, for each launch">
</picture>

### 3. It never joins the core later

The power in the core mode at every cross-section from the facet to the end of the fiber (blue),
and all the forward-going light inside the 28 λ window (orange). The guided power is flat to within
0.03% of the incident power: the uncoupled light runs alongside the core the whole way without
joining the mode, and nothing leaks out of the mode either. In a straight, uniform fiber the guided
mode and the rest are independent solutions of Maxwell's equations, and nothing moves power from
one to the other. The orange curve shows where the rest goes: for the too-big spot it hardly
drops; for the too-small spot it falls to 44% as the steep light leaves the window.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/power_along_z_dark.png">
  <img src="media/power_along_z_light.png" alt="Guided power and total forward power along the fiber">
</picture>

### 4. It is inside the acceptance cone; it is just in the wrong place

The ray picture of a fiber says that light is guided if it lands on the core inside the acceptance
cone. The uncoupled light of the too-big spot passes the direction test: decomposed into plane
waves just inside the facet, 97% of it travels within the core's acceptance angle (9.9° in the
glass), and most of it within 3°. It is not guided because it has the wrong *shape*: a broad
ring around the core with a dip inside it, which is orthogonal to the mode. And because it travels
so nearly parallel to the axis, it takes a long way to reach the outside of the cladding. The
too-small spot is the opposite: 89% of its uncoupled light is outside the cone, at angles up to
about 30°.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/uncoupled_angles_dark.png">
  <img src="media/uncoupled_angles_light.png" alt="Cumulative angular distribution of the uncoupled light inside the glass">
</picture>

### 5. Where it finally goes: out through the coating

Beyond 120 λ the long-distance model takes over. It includes the whole cladding (its surface is
22.9 λ from the axis in this scaled fiber; the 3D window ends at 14 λ) and an absorbing coating.
The uncoupled light of the too-big spot spreads slowly out to the cladding surface and is
absorbed by the coating a little at a time. Half of it is gone after 6.4 mm of SMF-28. At 2.9 cm,
the end of the run, 9% of the incident power is still in the cladding, falling by 1/e about every
2.5 cm: cladding modes, light at shallow angles that only grazes the coating now and then. The
uncoupled light of the too-small spot hits the coating at steep angles and is gone quickly: half
of it within 0.5 mm, 99% within 8 mm. In all three launches the guided power stays where the facet
put it (38%, 96%, 38%) the whole way.

<p align="center">
  <a href="media/long_distance.mp4"><img src="media/long_distance.webp" width="100%" alt="BPM movie: the whole fiber cross-section from 20 to 6000 wavelengths"></a>
</p>

![The same, far down the fiber: spot too big, matched, too small](media/long_distance_side.png)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/long_distance_budget_dark.png">
  <img src="media/long_distance_budget_light.png" alt="Where the light is, from 20 to 6000 wavelengths, for each launch">
</picture>

### 6. Why ray optics gets this wrong

Coupling efficiency against beam NA. The solid curve is the usual wave-optics estimate: the
overlap of the Gaussian with the mode, times the Fresnel transmission of the facet. The dashed curve
is the same estimate without the paraxial approximation. It sends the exact focused beam through
the facet one plane wave at a time, with its own s and p Fresnel coefficients. The dots are the 3D
FDTD: the three launches above and six more NAs from a sweep into a short (10 λ) fiber. They
follow the exact focused-beam estimate to within 0.3% of the incident power at every NA; the
paraxial estimate falls short above NA 0.3 (by 3.5% at NA 0.77).

Ray optics (dotted) counts the power that lands on the core inside the acceptance cone. On the
low-NA, big-spot side it is far too pessimistic: 16% of a spot 3× too big lands on the core, but 38% couples.
The mode reaches past the core, and the overlap adds field amplitudes coherently, so light landing
just outside the core still contributes. The reverse lesson matters more: a launch NA well inside
the fiber's NA is no guarantee of coupling. At NA 0.064, a quarter of the fiber's NA, 59% of the
light is lost, as much as at NA 0.58, more than twice the fiber's NA.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/efficiency_vs_na_dark.png">
  <img src="media/efficiency_vs_na_light.png" alt="Coupling efficiency against beam NA: wave optics, ray optics and 3D FDTD">
</picture>

### What this means at the bench

- **Spot too big (NA too narrow):** the light you lose is not reflected toward the source; that is only the usual
  3 to 4% Fresnel reflection of an uncoated facet, at any NA. It goes forward into the cladding as
  slow, shallow cladding modes and can last for centimetres. A power meter at the end of a short
  pigtail counts some of it, so the measured "coupling" reads high, then drops as the fiber gets
  longer or is bent. A cladding-mode stripper, or a few tight loops of fiber, removes it.
- **Spot too small (NA too wide):** about the same loss, but the uncoupled light leaves at steep angles and is gone
  within millimetres.
- **Matched:** the best beam NA is set by the mode field diameter, not the data-sheet NA. For
  SMF-28 at 1550 nm (mode field diameter 10.4 µm) that is λ / (π w₀) ≈ 0.095, against a quoted
  fiber NA of 0.14.

## The same beams into a single-mode metal waveguide

<p align="center">
  <a href="media/horn_movie.mp4"><img src="media/horn_movie.webp" width="100%"
     alt="FDTD movie: the same three beams switching on and filling a conical horn that feeds a single-mode metal waveguide, then the steady-state split into the guide's mode (blue) and everything else (yellow)"></a>
  <br><sub>The same three beams, now focused onto the mouth of a horn that feeds a single-mode,
  perfectly conducting circular waveguide (full Maxwell equations; the light comes down from the
  top, starting 8 λ in front of the mouth). The beam switches on, fills the horn and settles,
  shown at 2 periods a second; then the picture changes to where the
  power ends up: <b>white</b> = the incident beam, <b>blue</b> = the guide's mode heading in,
  <b>yellow</b> = everything else, all of which ends up reflected.
  <a href="media/horn_movie.mp4">Full-resolution MP4</a>.</sub>
</p>

**Short answer for the metal guide.** It comes back. A perfectly conducting guide has no
cladding. Outside its mouth is the metal face; inside, light can only travel as one of the
guide's modes, and each mode only where the guide is wide enough for it. So light that does not
end up in the guide's one mode has nowhere to go but back the way it came, and the simulations
find exactly that: guided plus reflected is 100% to within 0.01%. What differs between
the launches is *where* the light turns back.

- **Spot too big:** most of it (70%) lands on the metal face around the mouth and is
  reflected right there, as from a mirror with a hole in it. 24% reaches the guide.
- **Spot too small:** nearly all of it goes into the mouth (99.8%), but as the horn's
  *higher* modes: a small spot is a wide cone of directions, and the horn carries steep
  directions in its higher modes. Each of those turns back where the horn becomes too narrow for
  it, and 76% of the incident power comes back out of the mouth. 24% reaches the guide.
- **Matched:** 80% reaches the guide. The 20% that does not is mostly higher modes
  too: the field a smooth horn makes at its mouth is close to a Gaussian, but not equal to one.

These are the fiber's own three beams, focused on the horn's mouth as they were on the fiber's
facet. Two more sets of launches, focused at the depth inside the horn where each horn's
best beam is focused, and one of them into a better (corrugated) horn, are
[further down](#three-sets-of-launches-which-is-like-the-fiber).

### Why a metal guide needs a horn

A metal pipe is single mode only while it is less than 0.77 λ across. Its fundamental mode,
TE₁₁, needs a radius of at least 0.293 λ, and the next mode (TM₀₁) starts at 0.383 λ. The guide
here has a radius of 0.36 λ. Its mode is smaller than any spot a beam can be focused to and a
different shape, so a beam is not focused straight onto a bare single-mode guide. Microwave
practice feeds it through a horn: a gradual taper that turns the guide's mode into a beam the
size of the horn's mouth. The horn here is a smooth cone cut into the face of a thick metal
block, with a half-angle of 8°, 12.4 λ long, opening from the guide's radius to 2.1 λ at the
mouth (about the size of the fiber's core, 1.5 λ). After it, 6 λ of the single-mode guide end in
an absorber, which stands for a matched load.

That size makes the horn's matched beam the same as the fiber's. The best Gaussian focused on this
horn's mouth has a waist of 1.60 λ (0.76 × the mouth radius, the usual rule for a smooth-walled horn), and
the fiber's matched beam (1.66 λ) couples within 0.1 percentage point as well. So the same three beams, focused on
the mouth instead of the facet, are matched, too big and too small for both.

**Where a mismatched beam comes from.** In a microwave receiver the converging beam comes from
a dish or a lens, with the horn at its focus. The size of the focused spot is set by the dish's
focal ratio f/D: a shallow dish (large f/D) sends a narrow cone and makes a big spot; a deep dish
(small f/D) sends a wide cone and makes a small spot. "Spot too big" is a horn that is too small
for its dish. Antenna engineers, who usually think of the same horn transmitting, call that loss
*spillover*: the horn's beam is wider than the dish, and part of it misses. "Spot too small" is a
horn that is too big, too directive, for its dish: transmitting, it lights only the middle of the
dish (*under-illumination*, or taper loss). A feed horn is designed for the compromise between
the two.

### Where the light goes

The steady state on the plane of polarization, with the same brightness scale for all three
(per watt in). Everything yellow ends up reflected: in front of the mouth it is the reflected
light itself; inside the horn it is everything except the guide's mode heading in.

![Steady state in the metal guide and its horn: spot too big, matched, too small](media/horn_side_view.png)

- **Spot too big:** the beam is 3× wider than the best spot and much wider than the mouth. The
  part that falls on the metal around the mouth is reflected straight back toward the source,
  and stands in front of the face as a wave pattern with the incident beam. The part that falls
  into the mouth is mostly the right shape, and most of it (24% of the incident power)
  reaches the guide.
- **Spot too small:** the whole spot fits into the mouth, so almost nothing hits the face. The
  light spreads into the horn at steep angles, fills it with a pattern of higher modes, and
  comes back out, leaving the mouth as a wide, rippled beam heading back toward the source.

### Inside the horn: each wrong mode turns back at its own cut-off

The horn is a waveguide whose radius shrinks, so at each depth its modes are those of a pipe of
that radius. At the mouth seven of them can carry light of this polarization (TE₁₁, TM₁₁, TE₁₂,
TM₁₂, TE₁₃, TM₁₃, TE₁₄). Each needs a minimum radius, its cut-off. As the horn narrows they
reach their cut-offs one after another (dotted lines), and there each mode's light stops and
turns around. Only TE₁₁ fits through the throat into the guide. Counting the power mode by mode
shows it. The too-small spot's light enters mostly in the higher modes (dashed). Their power
drops in steps at the cut-offs, and the light heading back out (solid yellow) drops in the same
steps, because what turns back at a cut-off travels back out of the mouth. TE₁₁ heading in
(blue) carries on to the guide, moving by a few points on the way, mostly near the cut-offs.
That is where the horn trades a little amplitude between its modes, and where counting modes on
caps is least exact. The horn converts at most 1.4% of TE₁₁'s power (sent out from the guide, TE₁₁
stays 98.6 to 99.8% TE₁₁), but a small amplitude converted into TE₁₁ adds to or cancels part of
TE₁₁'s own, and so moves its power by more than its own share.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/horn_modes_dark.png">
  <img src="media/horn_modes_light.png" alt="Power in the horn's modes against depth, heading in and heading back out, for the three launches">
</picture>

### The fiber and the metal guide, side by side

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="media/fiber_vs_horn_dark.png">
  <img src="media/fiber_vs_horn_light.png" alt="Where the incident power goes, for a spot too big, matched and too small: the fiber, and the metal guide fed by a smooth horn (beams focused on its mouth or 3.7 wavelengths inside) or a corrugated horn">
</picture>

The chart has all three sets of launches into the metal guide
([below](#three-sets-of-launches-which-is-like-the-fiber)); the table compares the fiber with
the first, the fiber's own beams focused on the smooth horn's mouth.

| | single-mode fiber | single-mode metal guide, fed by its horn |
|---|---|---|
| what surrounds the mode | cladding: glass that carries light in any direction, just not guided | nothing that carries light: metal outside the mouth, and inside the horn only modes, each of which stops at its cut-off |
| where uncoupled light goes | forward into the cladding, then out through the coating | back the way it came: there is nowhere else |
| spot too big | 57% forward into the cladding, beside the core, lost slowly over centimetres | 70% reflected by the metal face around the mouth |
| spot too small | 59% forward at steep angles, out through the cladding within millimetres | 76% enters the horn in higher modes and turns back inside it |
| matched | 96% guided: only the Fresnel reflection is lost | 80% guided; 16% turns back inside the horn, because a smooth horn's mouth field is not quite Gaussian |
| reflected | 3 to 4% (the glass's Fresnel reflection), whatever the launch | everything that is not guided |

In both, the guided power is the overlap of the incoming field with the mode it has to match; the
coupling numbers themselves are of the same kind. What differs is what the rest can turn into. The
fiber offers the uncoupled light a continuum of forward-going states that are not guided, in the
cladding and beyond it, so it goes forward and leaks away to the side. A metal guide with
lossless walls offers nothing but its own modes, and below cut-off those carry no power, so the
only state left for the uncoupled light is the one it came in on, reversed.

Two consequences at the bench:

- **The guide never sees the mismatch.** The horn itself is well matched to the guide: TE₁₁ sent
  out of the guide has less than 0.1% of its power reflected back (a return loss better than
  30 dB). A network analyzer on the guide sees a good antenna. A badly focused beam shows up only as lost gain, with the
  missing power scattered back toward the dish and the source.
- **The fiber's loss can hide downstream; the metal guide's cannot.** A fiber's uncoupled light
  travels with the guided light for centimetres and can reach a detector. A metal guide's never
  enters the guide at all.

### Why: overfilling the core or the NA, and the horn's mouth

**Too big overfills the core, too small overfills the NA.** All three beams are focused on the
entrance and centred on the axis; only the spot's size differs. A focused spot's size and its
spread of angles are locked together (w₀ = λ / (π · NA)). The beams are diffraction limited and
so is the mode, so each covers the same area in a plot of position across the face against angle;
a mismatched beam is that area stretched one way and squeezed the other. A spot 3× too big is 3×
too wide, with a third of the mode's spread of angles: most of it lands beside the core, at angles
the fiber would hold. A spot 3× too small is 3× too narrow, with 3× the mode's spread of angles:
it all lands on the core, but much of it too steep to stay. The two lose about the same: the
coupling is 38.2% and 37.7% in the fiber, 24.0% and 24.2% in the horn, much as two Gaussians 3×
apart in size overlap by 36% whichever is bigger (the overlap is the same computed across the
face or across the angles). What differs is the light left over. The big spot's lies beside the
core and runs nearly parallel to it, so it travels on in the cladding for centimetres. The small
spot's starts on the core but is steep, so it reaches the coating within millimetres. The split
happens at the entrance. The part with the mode's shape stays mode-sized and travels on; only the
rest spreads. An ideal lens could turn either mismatched beam into the matched one without loss.

**Where the wrong angles go.** The fiber keeps only the angles within its NA. Steeper light is not
held by total internal reflection, so it leaves the core sideways and is lost through the cladding.
The horn's metal walls reflect at any angle, and a straight pipe the size of the mouth would carry
nearly all of what enters. What turns the light back is the taper. Each bounce off the converging
wall steepens a ray by 16°, twice the cone's half-angle, until the ray heads back out. Each bounce
also keeps the ray's closest distance to the cone's apex, because every tangent plane of a cone
contains the apex, and the ray turns around at that distance. So a ray reaches the guide only if
its line passes within 2.6 λ of the apex, the throat's distance from it. Through the centre of the
mouth, that means within 10° of the axis (NA 0.17), close to the matched beam's 0.19. Steeper rays
turn around nearer the mouth: the ray version of the higher modes' cut-offs. Counting rays
(`scripts/horn_estimates.py`) gets the trend, though not the numbers, as ray optics did for the
fiber:

| rays / FDTD, % of the incident power | spot too big | matched | spot too small |
|---|---:|---:|---:|
| off the face | 69.9 / 70.1 | 4.1 / 4.1 | 0.0 / 0.2 |
| turned back inside the horn | 1.1 / 5.9 | 27.2 / 15.7 | 83.5 / 75.6 |
| reach the guide | 29.0 / 24.0 | 68.7 / 80.2 | 16.5 / 24.2 |

**Why the matched beam couples better into the fiber.** The fiber's mode is almost exactly
Gaussian. Of the light that enters the glass, 99.6% is in the mode, and the only loss is the
facet's 3% Fresnel reflection. What the beam must match in the metal guide is the horn's mouth
field, not the guide's mode: sent out of the guide, TE₁₁ reaches the mouth 98.6 to 99.8% intact.
The mouth field matches a Gaussian much less well (`media/horn_estimates.txt`):

- It stops at the metal wall. The matched beam's tail outside the mouth, 4.1% of its power, lands
  on the face and is reflected.
- TE₁₁ is not round. Its width differs between the planes parallel and perpendicular to the
  polarization, and part of its field is cross-polarized. On a flat mouth it couples to the best
  Gaussian at only 86.4%.
- The 8° cone's wavefront is curved; its edge lags the centre by 0.15 λ. That brings the match to
  82.3%, against 81.2% from the FDTD extrapolated to an infinitely fine grid.

The curvature is the one part a launch can fix, and a better horn fixes most of the rest. Both are
tried [below](#three-sets-of-launches-which-is-like-the-fiber).

### Three sets of launches: which is like the fiber?

The fiber's mode is 99.6% Gaussian and has a flat wavefront on the facet. So its best beam is
focused on the facet, and the matched launch loses only the facet's reflection. The smooth horn's
mouth field is neither Gaussian nor flat, so a launch into it mixes two mismatches: the beam's
size, which this study is about, and the horn's own shape. Two more sets of launches take the
second away, one step at a time. Each set is three beams with one focus, too big, best and too
small by the same factors of 3 as the fiber's:

| set | horn | the beams' focus | waists: too big / best / too small | into the guide, % |
|---|---|---|---|---|
| focused on the mouth | smooth | on the mouth | 4.98 / 1.66 / 0.55 λ (the fiber's) | 24.0 / 80.2 / 24.2 |
| focused inside | smooth | 3.7 λ inside | 4.20 / 1.40 / 0.47 λ | 31.8 / 85.4 / 41.4 |
| corrugated | corrugated | 1.7 λ inside | 3.84 / 1.28 / 0.43 λ | 34.3 / 98.0 / 35.0 |
| *the fiber* | | on the facet | 4.98 / 1.66 / 0.55 λ | 38.2 / 96.4 / 37.7 |

The two sets focused inside need very steep too-small beams, and a source plane of finite size
misses about 5% of each (see [the box](#how-the-metal-guide-was-computed)). That steep light
carries a little of the guide's mode, out of phase with the rest. For the whole beam the
too-small spots couple 38.7% and 33.0%.

#### Focused inside the smooth horn

<p align="center">
  <a href="media/horn_focus_movie.mp4"><img src="media/horn_focus_movie.webp" width="100%"
     alt="FDTD movie: three beams focused 3.7 wavelengths inside the smooth horn, switching on, then the steady-state split into the guide's mode (blue) and everything else (yellow)"></a>
  <br><sub>The smooth horn again, with three beams around its own best beam, all focused 3.7 λ
  inside the mouth (dotted line). <a href="media/horn_focus_movie.mp4">Full-resolution MP4</a>.</sub>
</p>

The horn's wavefront at the mouth is curved, so a beam does better focused deeper, inside the
horn, where it arrives at the mouth still converging. From the reverse runs
(`media/horn_design.txt`), the best is a waist of 1.40 λ focused 3.7 λ in from the mouth. It
couples 85.4% (85.5% by reciprocity), against 80.2% for the matched beam focused on the mouth,
and converges with the grid toward roughly 87%. It reaches the mouth 1.63 λ wide, the flat
mouth's best width, with its wavefront curved on a radius of 14 λ, close to the 15 λ from the
mouth to the cone's apex. This is the "phase centre" a feed horn's data sheet gives: the point
inside the horn that its far field appears to come from.

![Steady state, the three beams focused 3.7 λ inside the smooth horn](media/horn_focus_side_view.png)

- **Spot too big:** as before, most of it (61%) lands on the face around the mouth; 32% reaches
  the guide.
- **Spot too small:** this one is no longer the too-big spot's mirror image. Focused 3.7 λ inside,
  it is still converging at the mouth, and already 2.6 λ wide there, wider than the mouth: 35% of
  it lands on the face. What enters is closer to the mode's size at the mouth (1.6 λ) than the
  too-big beam's 4.2 λ, so more of it reaches the guide: 41%.
- **Matched:** 85% reaches the guide; 4% lands on the face and 11% turns back inside the horn.
  What is left is the shape of TE₁₁ at the mouth, which no focus can fix.

#### A corrugated horn

<p align="center">
  <a href="media/corrugated_movie.mp4"><img src="media/corrugated_movie.webp" width="100%"
     alt="FDTD movie: three beams focused 1.7 wavelengths inside a corrugated horn, switching on, then the steady-state split into the guide's mode (blue) and everything else (yellow)"></a>
  <br><sub>A corrugated horn: the same cone, its wall cut into quarter-wave slots. Three beams
  around its best beam, focused 1.7 λ inside the mouth. <a href="media/corrugated_movie.mp4">Full-resolution MP4</a>.</sub>
</p>

The same 8° cone and mouth, with its wall cut into 61 slots across the axis, each 0.15 λ wide at
a pitch of 0.2 λ (`fiberlaunch/horn.py --wall corrugated`). A slot is a short-circuited radial
transmission line. A quarter wave deep, its open end makes the wall look open to the field along
the axis, just as a smooth wall is a short circuit for the field around it. So the field fades to
zero at the wall in every direction, as the fiber's mode fades into its cladding. The mode, HE₁₁,
is round, linearly polarized and 98% Gaussian (E ∝ J₀(2.405 ρ/a)). The quarter-wave depth comes
from the radial line's exact resonance: 0.26 λ at the mouth and 0.30 λ near the throat. Near the
throat, where the guide carries TE₁₁, the first 6 slots step from half a wave deep (a short
circuit at the wall, like a smooth wall) to a quarter wave, and convert TE₁₁ into HE₁₁.

Sent out of the guide, TE₁₁ leaves the mouth as a beam that couples 98.2% to the best Gaussian,
at 60, 80 and 120 points per wavelength alike: a waist of 1.30 λ focused 1.8 λ inside the mouth
(`media/horn_design_corrugated.txt`). The launches use 1.28 λ focused 1.7 λ in, the best of a first
reverse run (see the checks below), which couples within 0.04 points of it. The textbook
estimate, HE₁₁ on the cone's spherical wavefront, gives 98.2% for a waist of 1.26 λ focused 1.9 λ in. The design is not delicate. A scan at 40
points per wavelength, with converters of 3 to 10 slots and slots 0.95 to 1.1 × a quarter wave
deep, moved the coupling by at most 0.3 points.

![Steady state, the three beams focused 1.7 λ inside the corrugated horn](media/corrugated_side_view.png)

- **Spot too big:** 55% lands on the face around the mouth and 11% enters but turns back;
  34% reaches the guide.
- **Spot too small:** focused 1.7 λ inside, it reaches the mouth 1.4 λ wide, inside the 2.1 λ
  mouth, so only 6% lands on the face. What enters comes in at steep angles: 59% turns back
  inside the horn, and 35% reaches the guide.
- **Matched:** 98.0% reaches the guide. Of the 2% that does not, 0.8% is the beam's tail outside
  the mouth, reflected by the face, and 1.2% turns back inside the horn.

With the horn's own shape out of the way, the two mismatched spots lose the same again, as in
the fiber: 34% and 35% coupled. Two Gaussians 3× apart in size overlap by 36%, and the
corrugated horn's best beam couples 98%, which gives 35%.

Radio telescopes and quasi-optical systems use corrugated feeds for this reason, and they work
over a wide band. A dual-mode (Potter) horn gets close with a smooth wall: a step near the throat
turns some of TE₁₁ into TM₁₁, and if the two arrive at the mouth in phase, TM₁₁ cancels most of
TE₁₁'s field at the wall and its unround part. On a flat mouth the best mix (16% of the power in
TM₁₁) couples at 98.4% (`media/horn_estimates.txt`). But the two modes travel at different
speeds, so they stay in phase only over a narrow band.

#### Which is like the fiber?

The corrugated set. A fiber launch changes only the beam's waist. The focus stays where the
fiber's mode has its own waist, on the facet, where the mode's wavefront is flat, and the mode is
almost exactly Gaussian. For a horn, the place to focus is where the beam it radiates has its waist,
the phase centre, not the mouth. Of the three sets, only the corrugated one has both of the fiber's
properties: a mode 98% Gaussian, and every beam focused on the guide's own waist. It also gives
the fiber's pattern:

| into the guide, % | spot too big | matched | spot too small |
|---|---:|---:|---:|
| single-mode fiber | 38.2 | 96.4 | 37.7 |
| corrugated horn, focused 1.7 λ inside | 34.3 | 98.0 | 35.0 |
| smooth horn, focused 3.7 λ inside | 31.8 | 85.4 | 41.4 |
| smooth horn, focused on the mouth | 24.0 | 80.2 | 24.2 |

The two smooth sets mix in a second loss, which the fiber does not have. Focused on the mouth,
the beams are 3.7 λ away from the horn's focus as well as the wrong size. Even at the right focus,
TE₁₁'s shape costs the best beam 15%, and it splits the mismatched beams unevenly. A beam 3× too
big or too small overlaps a Gaussian by 36% either way. But the two differ from the best beam in
opposite directions: expanded in beams about the best one, their first correction has opposite
signs. So whatever part of the horn's field is not Gaussian adds to one's coupling and takes from
the other's.

The corrugated horn's mismatched spots couple a little less than the fiber's (34 to 35% against
38%). The fiber's mode is not quite Gaussian either. In the cladding it decays exponentially,
more slowly than a Gaussian, so it overlaps the spot 3× too big by 39% rather than 36%. The spot
3× too small is far from paraxial: the exact focused beam couples 37.5%, against 34.6% for a
paraxial Gaussian (`media/results.json`). The facet then reflects 3 to 4% of every launch. The
metal horn has no such reflection.

Where the uncoupled light goes is still opposite in every set: in the fiber it goes forward into
the cladding; in a metal guide it comes back.

### How the metal guide was computed

**One plane is enough.** The horn, the guide and a beam on their axis are all symmetric about
that axis. An x-polarized beam with a round spot has fields that vary around the axis only as
cos φ and sin φ, which is exactly half angular order m = +1 and half m = −1, mirror images of
each other. Meep's cylindrical mode solves the full 3D Maxwell equations for one angular order
on the (r, z) plane, so it gives the exact 3D field at the cost of a 2D run: 120 grid points
per wavelength, 43 across the guide's radius, about 7 hours per launch on 8 cores. The
first version simulated the horn on a 3D Cartesian grid and could not converge. Meep cannot
smooth a perfect conductor below the pixel, so the cone's wall became a staircase, and the
staircase scattered TE₁₁ into other modes: the matched beam's coupling came out as 41, 55 and
69% at 16, 20 and 28 points per wavelength (`fiberlaunch/horn_cartesian.py`, kept for the
record).

**The box.** 20 λ in radius (18 λ inside the absorbers), from 11 λ in front of the mouth to the
end of the guide. The beam is the same exact Gaussian as in the fiber runs, here built from
plane waves in cylindrical components (`fiberlaunch/beam.py`). It enters through a plane of
equivalent electric and magnetic currents 8 λ in front of the mouth, which radiates it forward
only (3 × 10⁻⁵ of the power goes backward). A plane of finite size cannot carry all of a very
steep beam. The too-small beams of the two sets focused inside the horn (waists 0.47 and 0.43 λ,
11.7 and 9.7 λ from their focus) have 4.8% and 5.2% of their power beyond 18 λ from the axis, and
that part is not launched (0.7% for the set focused on the mouth). It is the beam's steepest light.
Every number here is a fraction of the power that is launched, measured in a matching run with
no metal. The reverse runs give the coupling of the whole beam.

**Accounting**, as fractions of the incident power (from a matching run with no metal):

| | how it is measured |
|---|---|
| reflected | 1 − (net flux through a plane in air 2 λ in front of the mouth) / (incident power), plus the light that skims out sideways along the face (through the side of the box between that plane and the face; at most 0.74%) |
| into the guide | the flux through the guide, 1 and 3 λ behind the horn |
| reflected by the face | where the incident beam lands: its power through the plane just in front of the face, outside the mouth (from the run with no metal) |
| entered the horn, turned back inside it | reflected − (reflected by the face) |

**Counting modes in a horn.** The modes of a cone are spherical waves about its apex. On a
spherical cap centred on the apex, each mode has a uniform phase, and for a cone this narrow
its profile is that of a straight pipe whose radius is the cap's arc radius. So the field is
taken on such caps, one for every depth (on planes in the straight guide), and each mode's
power heading in and heading out comes from Lorentz reciprocity against the same mode heading
the other way, which separates the two directions exactly. Sent out from the guide, TE₁₁ is
recovered as 98.6 to 99.8% TE₁₁ on every cap at 120 points per wavelength.

**The guided part in any horn.** A corrugated horn's modes are hybrids of the slotted wall, not
those of a smooth pipe, so for it the field that ends up in the guide (blue in the pictures) is
found another way, which works for any horn. Lossless Maxwell's equations run backwards in time
as well as forwards, so the guided field of a launch is the field of the reverse run (TE₁₁ sent
out of the guide, on the same grid) reversed in time, scaled to the TE₁₁ heading in at the
throat. In the smooth horn this gives the same guided field as the caps, to 1 to 3% in |E|².

**Checks.**

- *Energy:* guided + reflected = 100% to within 0.01% for all nine launches, and the
  power heading down the horn is the same on every plane from the mouth to the guide (to 10⁻⁴).
- *Reciprocity:* the coupling also follows from one run in reverse: TE₁₁ sent out of the guide,
  its radiated field overlapped with each incoming beam in front of the mouth (Lorentz
  reciprocity; `scripts/horn_design.py`). For the beam each source plane launches, the reverse
  runs and the direct runs agree to 0.13 points or better:

  | reciprocity / direct run, % | spot too big | matched | spot too small |
  |---|---:|---:|---:|
  | smooth horn, focused on the mouth | 24.0 / 24.0 | 80.3 / 80.2 | 24.2 / 24.2 |
  | smooth horn, focused 3.7 λ inside | 31.8 / 31.8 | 85.5 / 85.4 | 41.3 / 41.4 |
  | corrugated horn, focused 1.7 λ inside | 34.3 / 34.3 | 98.1 / 98.0 | 35.0 / 35.0 |

  The whole beam's too-small spots couple 23.4%, 38.7% and 33.0%. The edge of the source plane
  changes them by 0.9, 2.6 and 2.0 points. An earlier version of the first set, with the source
  only 4 λ out on a plane 14 λ in radius, cut off less of that beam and got 23.7%. Reciprocity
  for that cut-off beam gives 23.7% too.
- *Resolution:* the reverse runs at 20, 40, 60, 80 and 120 points per wavelength
  (`media/horn_design.txt`): the matched coupling is 44, 72, 78, 79 and 80%, converging as
  1/res² toward 81.2%; the too-small spot's toward 23.9%; the too-big spot's stays at 24.0 ±
  0.2% from 60 on. Focused inside, the best beam converges more slowly: 77, 83, 84 and 85.5% at
  40 to 120, toward roughly 87%. The corrugated horn's runs at 40, 60, 80 and 120
  (`media/horn_design_corrugated.txt`) give 97.3, 98.2, 98.2 and 98.2% for its best beam, and the
  fiber's three beams move by less than half a point from 60 to 120. The slots add an
  uncertainty of their own. At 120 points per wavelength they are 18 cells wide, so their edges
  fall exactly on cell boundaries, and how Meep averages a cell there depends on rounding, and so
  on where the box starts. The first reverse run at 120 started 4 λ nearer the mouth than the
  launches. It averaged 15 of the 122 slot edges differently, by a fraction of a cell, and its
  couplings differed from the launches' by 0.4 to 1 point (33.7, 98.4 and 35.9%). Rerun on the
  launches' grid, it agrees with them to 0.13 points (the table above); the first run is kept in
  `data/horn/design/corrugated_emit_res120_zs4`. So the corrugated horn's numbers stand for the
  ideal slotted wall to about ±1 point.
- *The split of the reflected light:* by where the incident beam lands, against counting the
  horn's modes heading back out just inside the mouth (smooth horn only). For the light turned
  back inside the horn, the two give 5.9 / 5.7, 15.7 / 15.8 and 75.6 / 75.1% for the beams focused
  on the mouth, and 7.3 / 7.1, 10.8 / 11.4 and 23.9 / 24.0% for those focused inside.
- *Theory:* the textbook estimate fills the mouth with the guide's mode (TE₁₁ for a pipe of the
  mouth's radius) on a spherical wavefront centred on the cone's apex. It gives 24.2%, 82.1% and
  24.9%, close to the simulation, so the horn behaves as the standard picture says. Focused
  inside, it gives 86.6% for a waist of 1.41 λ focused 3.4 λ in. For the corrugated horn (HE₁₁ on
  the same wavefront) it gives 98.2%, the same as the simulation.

## The numbers

All as percentages of the power the beam carries in air (`media/results.json`,
`media/horn_results.json`, `data/bpm/summary.json` once the simulations are run).

| | spot too big<br>(NA too narrow) | spot matched | spot too small<br>(NA too wide) |
|---|---:|---:|---:|
| beam NA (λ / π w₀) | 0.064 | 0.19 | 0.58 |
| spot radius w₀ on the facet | 4.98 λ | 1.66 λ | 0.55 λ |
| ***3D FDTD, 115 λ into the fiber*** | | | |
| reflected by the facet | 3.0 | 3.2 | 3.4 |
| in the core mode | **38.2** | **96.4** | **37.7** |
| uncoupled, near the core (in the 28 λ window) | 57.4 | 0.1 | 6.2 |
| uncoupled, gone sideways (out of the window) | 1.3 | 0.3 | 52.7 |
| ***long-distance model, 6000 λ (2.9 cm of SMF-28)*** | | | |
| in the core mode | 38.2 | 96.4 | 37.7 |
| uncoupled, still in the cladding | 9.3 | 0.02 | 0.1 |
| uncoupled, absorbed by the coating | 49.5 | 0.4 | 58.8 |
| ***theory for the coupling*** | | | |
| mode overlap × Fresnel (paraxial) | 38.0 | 96.1 | 34.6 |
| exact focused beam, plane wave by plane wave | 38.1 | 96.1 | 37.5 |
| ray optics × Fresnel | 16.0 | 75.1 | 30.3 |
| ***metal waveguide, smooth horn, the same beams focused on its mouth (FDTD, 120 px/λ)*** | | | |
| into the guide's mode (TE₁₁) | **24.0** | **80.2** | **24.2** |
| reflected by the metal face around the mouth | 70.1 | 4.1 | 0.2 |
| entered the horn, turned back inside it | 5.9 | 15.7 | 75.6 |
| reciprocity (one reverse run), extrapolated to an infinitely fine grid | 24.0 | 81.2 | 23.9 (the whole beam) |
| aperture theory (TE₁₁ on the mouth, spherical wavefront) | 24.2 | 82.1 | 24.9 |
| ***smooth horn, waists 4.20 / 1.40 / 0.47 λ focused 3.7 λ inside*** | | | |
| into the guide's mode (TE₁₁) | **31.8** | **85.4** | **41.4** |
| reflected by the face | 60.9 | 3.9 | 34.7 |
| entered the horn, turned back inside it | 7.3 | 10.8 | 23.9 |
| ***corrugated horn, waists 3.84 / 1.28 / 0.43 λ focused 1.7 λ inside*** | | | |
| into the guide's mode (TE₁₁) | **34.3** | **98.0** | **35.0** |
| reflected by the face | 55.2 | 0.8 | 5.8 |
| entered the horn, turned back inside it | 10.5 | 1.2 | 59.2 |
| ***3D FDTD, extrapolated to an infinitely fine grid*** | | | |
| reflected by the facet | 3.4 | 3.6 | 3.8 |
| in the core mode | 38.0 | 96.0 | 37.4 |

## More movies

**In 3D.** A volume rendering of the steady state, each launch cut open along its axis, with the
camera moving from in front of the facets round to the side. Same colours as above; brightness is
log intensity on one scale for all three.

<p align="center">
  <a href="media/launch_3d.mp4"><img src="media/launch_3d.webp" width="100%" alt="3D volume rendering of the three launches: the beam in air (white), the guided mode (blue) and the uncoupled light (orange)"></a>
</p>

**Flying down the fiber.** The cross-section of the fiber at each distance from the facet, from
the beam in air to the end of the 3D simulation, with the running power budget underneath.

<p align="center">
  <a href="media/cross_sections.mp4"><img src="media/cross_sections.webp" width="100%" alt="Cross-sections of the three launches from the facet to 120 wavelengths"></a>
</p>

## How it was computed

**The fiber.** A step-index fiber with the same V number as SMF-28 (V = 2.36, single mode),
in units of the wavelength: core radius 1.5 λ, NA 0.25, cladding index 1.45. A real SMF-28 core
is 4.1 µm at 1.55 µm, i.e. 2.6 λ with NA 0.14. This fiber is that one shrunk across by 1.76 and
with the NA raised by the same factor, which keeps the guided mode's shape (in units of the core)
and the coupling numbers exactly the same while cutting the 3D volume by about 10. Diffraction
lengths along the fiber shrink by 1.76² = 3.1; the "cm of SMF-28" axes undo that.

**The launches.** A Gaussian beam focused onto the cleaved facet from air, x-polarized, at
three NAs: a spot matched to the mode (waist 1.66 λ, NA 0.19), a spot 3× too big (waist 4.98 λ,
NA 0.064: too narrow) and a spot 3× too small (waist 0.55 λ, NA 0.58: too wide). "Beam NA" is λ/(π w₀), the 1/e² half-angle. The beam is
Meep's exact focused-beam source (a complex point source), not a paraxial approximation, which
matters for the too-small spot. At NA 0.58 the exact shape of a real focus depends on the lens;
this one is an exact solution of Maxwell's equations with the Gaussian's waist.

**3D FDTD** with [Meep](https://meep.readthedocs.io) 1.35: 16 grid points per wavelength
(11 per wavelength inside the glass), a 32 λ × 32 λ × 134 λ box including PML absorbers
(28 λ × 28 λ × 120 λ of fiber, and the 8 λ of air the beam comes in through), 4-fold mirror
symmetry, about 3 hours per launch on 16 cores (the runs here shared the machine with others, 8
MPI processes each, and took 8 to 9.5 hours). The source is a continuous wave that switches on smoothly; the run continues until
every part of the field inside the fiber has reached steady state, and then the field is
Fourier-transformed over 4 more periods on:

- a stack of cross-sections every 0.5 λ (all six field components),
- the two side planes through the axis,
- flux planes for the reflection and for the transmitted power.

**Splitting the light.** The guided part at each cross-section is the projection of the
FDTD field onto the fiber's HE₁₁ mode, computed by Meep on the same grid in a separate run so
that both share the same discretization:
a = [∫(E × h*) + ∫(e* × H)]·ẑ dA / 2∫(e × h*)·ẑ dA, so that |a|² is the guided power.
Everything else is "the rest". Meep's own eigenmode coefficients (MPB) agree to 0.05%.

**Accounting for the power**, as fractions of the power the beam carries in air (from a
matching run with no fiber):

| | how it is measured |
|---|---|
| reflected | 1 − (net flux through a plane in air between the source and the facet) / (incident power) |
| in the core mode | \|a\|² from the projection, fitted as A e^{iβz} over all the cross-sections |
| uncoupled, near the core | flux through the 28 λ × 28 λ cross-section at the end of the fiber, minus the guided power |
| uncoupled, gone sideways | the remainder: light that has crossed the edge of the 28 λ window (still inside a real 125 µm cladding) |

**Long distances.** 120 λ of fiber is only 0.6 mm of SMF-28. To follow the uncoupled light
until it is gone, the FDTD field 5 λ inside the facet is handed to a beam-propagation method
(`fiberlaunch/bpm.py`): a scalar split-step Fourier method with the exact (non-paraxial)
propagator in the glass, on the whole fiber cross-section: core, the full 125 µm cladding and
a higher-index, absorbing 190 µm coating, scaled like the core. It starts from exactly the
FDTD's guided / uncoupled split and runs to 6000 λ (2.9 cm of SMF-28) in 2 to 13 minutes. Over
the 110 λ where both exist, it agrees with the 3D FDTD on the guided power to 0.2% and on the
power left inside the 3D window to 2% of the incident power; the fields themselves correlate at
0.96, 1.00 and 0.94 for the three launches (the scalar model is least exact for the steep light of
the too-small spot).

**Theory**, for comparison (`fiberlaunch/fibermode.py`): the exact HE₁₁ mode of the step-index
fiber; the overlap of the incident Gaussian with the (scalar LP₀₁) mode; the Fresnel loss at the
facet; and the ray-optics estimate, which counts the rays that land on the core inside its
acceptance cone. Because the NA 0.58 beam is far from paraxial, there is also an exact version
of the overlap estimate (`exact_beam_coupling` in `fiberlaunch/analysis.py`): the FDTD's own
incident beam in air is split into plane waves, each is sent through the facet with its own s and
p Fresnel coefficients, and the transmitted field is projected onto the HE₁₁ mode.

## How much to trust it

- **Against theory.** The coupling the 3D simulation finds is within 0.5% (relative) of the exact
  focused-beam estimate for all three launches: 38.2% against 38.1%, 96.4% against 96.1%, 37.7%
  against 37.5%. The paraxial estimate is 9% low for the too-small spot, whose NA 0.58 beam is not paraxial.
  The FDTD's guided mode travels with n_eff = 1.46105 once the grid's known numerical dispersion
  is taken out; the exact HE₁₁ value is 1.46107.
- **Resolution.** The same launches into a 12 λ fiber at 10, 12, 16 and 20 points per
  wavelength (`media/convergence.txt`): from 16 to 20 the coupling changes by at most 0.15%
  of the incident power. Extrapolated to an infinitely fine grid (the error falls as 1/res²) it is 38.0%, 96.0% and 37.4%, within 0.2% of
  the exact focused-beam estimate. The facet's reflection converges more slowly: the 3.0 to 3.4%
  found at 16 points per wavelength extrapolates to 3.4 to 3.8%, close to the 3.5% Fresnel
  reflection of a plane wave at the mode's effective index. It is the only number here that
  the grid changes noticeably.
- **Energy.** The forward power 3 λ inside the facet equals the incident minus the reflected
  power to within 10⁻⁴, for all three launches: the facet loses no power that is not counted.
- **Steady state.** Sampled from the movie frames, the field amplitude anywhere along the
  fiber changes by less than 10⁻⁴ over the last 10 periods before the Fourier transform
  (`media/steady_state.txt`), and the guided power is the same at every cross-section from the
  facet to the end of the fiber (scatter at most 3 × 10⁻⁴).
- **Movie frames.** The runs save the field every 0.34 periods, about three frames a period.
  Played back as they are, a standing wave blinks and a travelling one can seem to run
  backwards. The movies show 15 frames a period instead, at 2 periods a second. The frames in
  between are fitted: the source has one frequency and switches on over several periods, so at
  each point the field is Re[A(t) e^(−2πit)] with A changing slowly, and a quadratic A(t) fitted
  to the ten nearest saved frames gives the field in between (`fiberlaunch/frames.py`). Rebuilt
  this way from the others, each saved frame comes back to within 0.12% of the steady field's
  peak, at any time, in the fiber and the horns alike: about the float16 rounding of the saved
  frames themselves.
- **The metal guide** has its own checks (energy, reciprocity, resolution, theory): see
  [How the metal guide was computed](#how-the-metal-guide-was-computed).
- **What is idealized.** A perfectly cleaved, uncoated facet at normal incidence, a straight
  fiber, a monochromatic continuous wave, a perfectly aligned beam, and a coating modelled as a
  uniform absorbing layer (real coatings have two layers and a rough outer surface; both strip
  cladding light faster). The fiber is scaled from SMF-28 as described above, which is exact
  for the coupling and for paraxial propagation, and close otherwise.

## Reproducing it

Everything is driven by `make` (Linux, about 2½ days on a 16-core desktop, a day and a half of
it for the metal guide; about 25 GB of raw data in `data/`, which is not in git):

```sh
make envs     # two conda environments: Meep (MPI) and PyVista (off-screen)
make sims     # 3D FDTD: the three launches, the NA sweep, the resolution study; then the BPM;
              # then the metal guide: the two horns' resolution studies, and three sets of
              # three launches
make media    # figures, movies, the 3D render and the README previews
make test     # unit tests: the mode solver, the BPM, the exact-beam coupling estimate, the
              # exact beam in cylindrical components and the metal guide's modes
```

Individual pieces: `scripts/run_main.sh`, `scripts/run_sweep.sh`, `scripts/run_convergence.sh`,
`scripts/run_horn.sh`, and one run directly with
`mpirun -np 16 python -m fiberlaunch.fdtd --kind fiber --case narrow --movie --xy-planes` or
`mpirun -np 8 python -m fiberlaunch.horn --set corrugated --kind horn --case narrow --movie --out data/horn/corrugated/res120/narrow`.
`scripts/env.sh` holds the settings (number of MPI processes, paths). On a machine with
InfiniBand or RoCE adapters, it also sets `UCX_TLS=self,sm,tcp` so the MPI library uses
shared memory instead of probing the adapters, which fails noisily without an unlimited
locked-memory limit.

## Files

| | |
|---|---|
| `fiberlaunch/params.py` | the fiber, the three launches, the box |
| `fiberlaunch/fdtd.py` | the Meep simulation (mode reference, incident beam in air, launch into the fiber) |
| `fiberlaunch/analysis.py` | mode projection, power accounting, guided / uncoupled split, exact-beam coupling |
| `fiberlaunch/fibermode.py` | analytic fiber mode, overlap, ray optics, Fresnel |
| `fiberlaunch/bpm.py` | the long-distance beam propagation |
| `fiberlaunch/horn.py` | the metal guide and its horn, smooth or corrugated: Meep on the (r, z) plane (beam in air, launch into the horn, TE₁₁ sent out of the guide) |
| `fiberlaunch/beam.py` | the exact focused Gaussian beam in cylindrical components |
| `fiberlaunch/horn_analysis.py` | the metal guide's power budget, and the mode split on caps about the cone's apex |
| `fiberlaunch/horn_cartesian.py` | the first, 3D Cartesian horn, which could not converge (kept for the record) |
| `fiberlaunch/style.py` | colours shared by every figure and movie |
| `fiberlaunch/frames.py` | the movies' frames between the saved ones |
| `scripts/run_*.sh`, `run_bpm.py` | the simulation campaigns |
| `scripts/fig_*.py`, `movie_*.py`, `media_bpm.py`, `export_volumes.py`, `render_3d.py` | figures and movies (`fig_scales.py` draws the zoom-out, `fig_incoming.py` the incoming beams) |
| `scripts/check_steady.py`, `summarize_convergence.py` | the steady-state and resolution checks |
| `scripts/horn_design.py`, `horn_media.py` | the metal guide: coupling of every beam by reciprocity, the best focus inside the horn, and the resolution study; its figures and movie |
| `scripts/horn_estimates.py` | the metal guide's quick estimates: the beam outside the mouth, the mouth field's match to a Gaussian (this horn's, a corrugated horn's and a dual-mode horn's), ray optics |
| `media/` | everything shown here, plus `results.json`, `horn_results.json`, `convergence.txt`, `horn_design.txt`, `horn_design_corrugated.txt`, `horn_estimates.txt` and `steady_state.txt` |
| `tests/` | unit tests |
| `README.md`, `STUDY.md` | the short version, and this full study |
| `CLAUDE_CODE_CHAT.md` | the conversation with Claude Code that produced this |

## License

[MIT](LICENSE).
