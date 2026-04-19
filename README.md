# LammSolvF

Fortran implementation of a finite element solver for the Lamm equation using a moving grid, extending the method of Claverie et al. (1975). Includes a Python interface.

## Reference

This implementation is based on the method described in:

> P. Schuck, "Sedimentation Analysis of Noninteracting and Self-Associating Solutes Using Numerical Solutions to the Lamm Equation," *Biophysical Journal*, 75(3):1503–1512, 1998.
> https://doi.org/10.1016/S0006-3495(98)74069-X

## Usage

```python
import numpy as np
from lammsolvflib import sawp_solver

N     = 500
rpm   = 50000
omega = rpm * 2 * np.pi / 60  # rad/s

times, coeffs = sawp_solver(
    c0      = np.ones(N),   # uniform initial concentration, length N where N is the number of grid points used
    t0      = 0.0,          # initial time (s)
    n_steps = 200,
    omega   = omega,        # rad/s
    s       = 1e-13,        # sedimentation coefficient (5 S = 5e-13 s)
    D       = 1e-7,        # diffusion coefficient (cm²/s)
    m       = 6.0,        # meniscus radius (cm)
    b       = 7.2,        # bottom radius (cm)
    dt      = 50.0,         # time step (s)
)

# times:  shape (n_steps + 1,)   — time points in seconds
# coeffs: shape (n_steps + 1, N) — concentration at each time and grid node
```

## Disclaimer

This library is provided as-is with no warranty. No validation against the original implementation has been attempted. Use at your own risk.

## Installation

```bash
pip install .
```

Requires a Fortran compiler (gfortran) and LAPACK/BLAS system libraries.
