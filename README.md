# LammSolvF

Fortran implementation of the SAWP (Schuck Adaptive Wave Propagation) solver for sedimentation analysis, with a Python interface.

## Reference

This implementation is based on the method described in:

> P. Schuck, "Sedimentation Analysis of Noninteracting and Self-Associating Solutes Using Numerical Solutions to the Lamm Equation," *Biophysical Journal*, 75(3):1503–1512, 1998.
> https://doi.org/10.1016/S0006-3495(98)74069-X

## Disclaimer

This library is provided as-is with no warranty. No validation against the original implementation has been attempted. Use at your own risk.

## Installation

```bash
pip install .
```

Requires a Fortran compiler (gfortran) and LAPACK/BLAS system libraries.
