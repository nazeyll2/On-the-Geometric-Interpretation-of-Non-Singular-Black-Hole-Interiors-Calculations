# On-the-Geometric-Interpretation-of-Non-Singular-Black-Hole-Interiors-Calculations

Computational calculations for Lorentz boosts, torsion invariants, and diffeomorphism-related structures used in the accompanying research paper.

## Overview

This repository contains symbolic calculations related to teleparallel gravity, Lorentz covariance, torsion invariants, and geometric structures appearing in non-singular black hole interior models.

The computations are implemented using Python and SymPy.

Topics covered include:

- Lorentz boost covariance tests
- Torsion tensor constructions
- TEGR torsion scalar calculations
- Torsion invariants
- Radial shift tetrads
- Diffeomorphism-related structures
- Symbolic tensor manipulations
- Radial and transverse torsion eigenvalues


## Files

### `torsion_invarients.py`

Computes torsion invariants and verifies their behavior under radial Lorentz boosts with spin connection contributions included.

Main quantities:
- \( I_1 = T_{abc}T^{abc} \)
- \( I_2 = T_{abc}T^{cba} \)
- torsion vector \( T_a \)
- axial vector \( A^a \)
- TEGR torsion scalar


### `blackhole_lorentz_covarience.py`

Tests Lorentz covariance properties of radial shift tetrads in teleparallel gravity.

Includes:
- boosted tetrads
- generated spin connection
- transformed torsion components
- covariance consistency checks


### `black_hole_torsiyon_invaryantları.py`

Computes the torsion scalar for generic and radial shift configurations.

Includes:
- linearized torsion scalar expansion
- divergence structure of radial shifts
- analytic consistency checks
- radial torsion scalar expressions


## Requirements

Install dependencies with:

```bash
pip install sympy


## Usage

Run any script with:

```bash
python filename.py


## Notes

These calculations were developed for research purposes related to geometric and teleparallel descriptions of non-singular black hole interiors.

## License

MIT License
