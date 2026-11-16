# Liquid argon, after Rahman (1964)

In 1964 Aneesur Rahman followed 864 argon atoms on a CDC 3600 computer and
reported how they arrange themselves and how fast they wander: A. Rahman,
"Correlations in the Motion of Atoms in Liquid Argon", Phys. Rev. 136, A405
(1964). We ran the same liquid with our own code, written from the ground
up in Rust with no third-party crates, and compared the two.

Our diffusion coefficient is 2.57 x 10^-5 cm^2/s, 5.8% above his
2.43 x 10^-5. The structure matches the numbers others have quoted from
his paper: the first peak of g(r) is 0.03 angstrom from his, and the four
peaks of the structure factor are within 1.2% of his.

## The run

| | Rahman | ours |
|---|---|---|
| Atoms | 864 | 864 |
| Box side | 10.229 sigma | 10.229 sigma |
| Potential | Lennard-Jones, sigma 3.4 angstrom, eps/kB 120 K | the same |
| Cutoff | 2.25 sigma | 2.25 sigma, energy shifted to zero there |
| Temperature | 94.4 K | 94.79 +- 0.03 K, measured |
| Time step | 10^-14 s, predictor-corrector | 1.08 x 10^-14 s (0.005 tau), velocity Verlet |
| Measured for | about 8 ps (a secondary source) | 430 ps (200 tau) |

Shifting the energy at the cutoff leaves the forces as his, so the motion
is the same. Our atoms start on a face-centered cubic lattice with random
velocities for 94.4 K, seed 1964. For 20 tau (43 ps) the velocities are
rescaled to 94.4 K every 10 steps, and the lattice melts. Then the atoms
move freely, 10 tau to settle and 200 tau while we measure, in 10 blocks of
20 tau. These are our own choices: how Rahman prepared his liquid was not
found. Nothing holds the temperature while we measure. It stays 0.4% above
his without drifting: at the last rescaling the potential energy happened
to be a little high, and the total energy is fixed from then on. Over the
40,000 measured steps it changed by at most 3.8 x 10^-4 of itself.

## Results

| | ours | Rahman |
|---|---|---|
| Diffusion coefficient, from the MSD | 2.571 +- 0.025 x 10^-5 cm^2/s | 2.43 x 10^-5 |
| Diffusion coefficient, from the VACF | 2.559 +- 0.022 x 10^-5 cm^2/s | |
| First peak of g(r) | 3.70 angstrom | 3.67 |
| Peaks of S(k), k sigma | 6.72, 12.51, 18.41, 24.53 | 6.8, 12.5, 18.5, 24.8 |
| VACF first below zero | 0.318 ps | 0.33 |

The X-ray measurement Rahman compared with (Eisenstein and Gingrich, Phys.
Rev. 62, 261 (1942), at 91.8 K) puts the peaks of S(k) at 6.8, 12.3, 18.4
and 24.4. Each +- is one standard error of the mean over the 10 blocks.

## Diffusion

An atom in a liquid is jostled by its neighbors and wanders off in a random
walk. The diffusion coefficient D says how fast: the mean squared distance
from where the atoms started grows as 6 D t. We measure D two ways. One is
a sixth of the slope of that mean squared displacement (MSD), fitted from
2 to 10 tau. The other adds up how long each atom remembers its velocity:
D is a third of the integral of the velocity autocorrelation function
(VACF), the Green-Kubo formula. Both are averaged over every starting time
in a block. They use the same motion and agree to 0.5%.

[D from the VACF integral levels off at our D from the MSD][diffusion]

Our D is 5.8% above Rahman's, inside the 10% that PHYSICS_REFERENCE sets as
the pass line. With his atomic mass unit, which makes tau 0.4% longer, ours
converts to 2.56, 5.4% above his. Our error bar is 1%, so the gap is real
for our run. Two likely reasons, neither checked against his paper:

- His run was short: 780 steps, 7.8 ps, by a secondary source. Our blocks
  of 43 ps scatter by 3%. If the scatter grows as one over the square root
  of the length, one run of 7.8 ps scatters by about 7%. Our estimate.
- Our temperature is 0.4 K above his, and D rises with temperature. We did
  not measure by how much.

Later runs at 94.4 K and 1.374 g/cm^3 give 2.39 to 2.48 x 10^-5 (S. H. Lee
and coworkers, 1991 to 2014, 216 to 110,592 atoms; H. Heidenreich's rerun
of Rahman's setup in LAMMPS, 2.47). With today's atomic mass unit their
liquid is 1% denser than his, and a denser liquid diffuses more slowly.

A small box slows diffusion too. By Yeh and Hummer's correction, with the
viscosity Lee measured at this state, 864 atoms diffuse about 7% slower
than a large sample. That is our estimate, and it applies to his run and
ours alike. Experiment gives about 2.83 x 10^-5 at 94.4 K, and Rahman's
abstract says his D is 15% below experiment.

## Structure

g(r) is how crowded the liquid is at distance r from an atom, relative to
the average. It is zero inside an atom's core, peaks at the shell of
nearest neighbors, and settles to 1 far away.

[g(r) of liquid argon, with Rahman's first peak][gr]

S(k), the structure factor, is the Fourier transform of g(r) and what X-ray
or neutron diffraction measures. Rahman compared its peaks with X-ray data.

[S(k), with Rahman's four peaks marked on our curve][sk]

Our first peak of g(r) is at 3.70 angstrom, 0.03 angstrom past the 3.67
that Berne, Boon and Rice read from his g(r). Our peaks of S(k) are within
1.2% of his, and the last two are closer to the X-ray peaks than his are.
Our S(k) carries small ripples, about 1.2 apart in k sigma, because g(r) is
cut off at half the box, 5.1 sigma, before it has settled to 1.

The ROADMAP asks that the structure match the published curve. His g(r) is
a figure in a paper behind a paywall, and we found no copy of its numbers,
so we compared the numbers others quote from it. Drawing his whole curve
over ours needs a copy of the paper, to read the points off his figure.

## The velocity autocorrelation function

The VACF says how much of an atom's velocity is left after time t, on
average. It starts at 1 and falls as neighbors knock the atom about. It
turns negative at 0.318 ps: the atom has hit the cage of its neighbors and
bounced back. Rahman found the same at 0.33 ps. Our lowest point is -0.082
at 0.44 ps.

[The VACF, with Rahman's first zero][vacf]

## Labels

- Derived: D from the slope of the MSD (Einstein) and from the integral of
  the VACF (Green-Kubo), S(k) from g(r).
- Approximated: the error bars treat the blocks as independent, which holds
  when a block is far longer than the 2 tau the VACF lasts. The integrals
  use the trapezoid and midpoint rules. Peaks are placed at the top of a
  parabola through the three highest points.
- Assumed: the Lennard-Jones potential for argon, and Rahman's box side,
  which follows from his value of the atomic mass unit (PHYSICS_REFERENCE).

## Sources

We could read only the abstract of his paper,
https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405 . His other
numbers are quoted by: B. J. Berne, J. P. Boon and S. A. Rice, J. Chem.
Phys. 45, 1086 (1966) (D, first peak); J. C. Modla, PhD thesis, Lehigh
University (1977) (VACF zero); and H. Heidenreich's notes on the paper (D,
VACF zero, S(k) peaks, length of his run),
https://hunterheidenreich.pages.dev/notes/chemistry/molecular-simulation/classical-methods/correlations-motion-atoms-liquid-argon/ .
Later runs: S. H. Lee, J. Korean Chem. Soc. 58, 500 (2014), the papers it
cites, and https://hunterheidenreich.com/posts/rahman-1964-lammps-liquid-argon/ .
Finite size: I.-C. Yeh and G. Hummer, J. Phys. Chem. B 108, 15873 (2004).
Experiment at 94.4 K: interpolated from G. A. Cook (1961) by S. H. Lee, D.
K. Park and D. B. Kang, Bull. Korean Chem. Soc. 24, 178 (2003).

[diffusion]: argon-diffusion.svg
[gr]: argon-gr.svg
[sk]: argon-sk.svg
[vacf]: argon-vacf.svg