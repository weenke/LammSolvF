"""High-level Python wrapper for the SAWP solver."""
from __future__ import annotations

from ctypes import byref, c_double, c_int

import numpy as np

from lammsolvflib.lib import lib
from lammsolvflib.typing import Array1D, Array2D, ArrayLike1D


def sawp_solver(
    c0: ArrayLike1D,
    t0: float,
    n_steps: int,
    omega: float,
    s: float,
    D: float,
    m: float,
    b: float,
    dt: float,
) -> tuple[Array1D, Array2D]:
    """Run the SAWP time integrator.

    Args:
        c0: Initial concentration coefficients, length N. Determines the
            number of grid points N.
        t0: Reference (initial) time (s).
        n_steps: Number of time steps to integrate.
        omega: Angular velocity (rad/s). Convert from rpm via
            ``omega = rpm * 2 * np.pi / 60``.
        s: Sedimentation coefficient (s). 1 Svedberg = 1e-13 s.
        D: Translational diffusion coefficient (m²/s).
        m: Meniscus radius, lower grid boundary (m).
        b: Bottom radius, upper grid boundary (m).
        dt: Time step (s).

    Returns:
        times: Time points in seconds, shape ``(n_steps + 1,)``.
        coeffs: Concentration coefficients, shape ``(n_steps + 1, N)``.
            Row index = time point, column index = grid node.

    Raises:
        ValueError: If input validation fails or the solver returns a
            non-zero LAPACK factorisation error code.
    """
    c0 = np.ascontiguousarray(c0, dtype=np.float64)
    N = c0.shape[0]

    if c0.ndim != 1:
        raise ValueError(f"c0 must be 1-D, got shape {c0.shape}.")
    if N < 2:
        raise ValueError(f"c0 must have at least 2 elements, got {N}.")
    if n_steps < 1:
        raise ValueError(f"n_steps must be >= 1, got {n_steps}.")
    if m >= b:
        raise ValueError(
            f"Meniscus m ({m}) must be strictly less than bottom b ({b})."
        )

    # Query flat buffer lengths from the Fortran side
    times_len = c_int()
    coeffs_len = c_int()
    lib.sawp_query_sizes(N, n_steps, byref(times_len), byref(coeffs_len))

    # Allocate flat output buffers (filled in-place by the solver)
    times_out = np.zeros(times_len.value, dtype=np.float64)
    coeffs_out = np.zeros(coeffs_len.value, dtype=np.float64)

    stat = lib.sawp_solver(
        c0,
        t0,
        n_steps,
        omega,
        s,
        D,
        m,
        b,
        N,
        dt,
        times_out,
        coeffs_out,
    )

    if stat != 0:
        raise ValueError(
            f"sawp_solver returned non-zero LAPACK info code {stat}: "
            "LHS tridiagonal factorisation failed (singular or near-singular system)."
        )

    # coeffs_out is flat Fortran column-major with shape (N, n_steps+1).
    # Reshape using order='F' then transpose to row-major (n_steps+1, N).
    coeffs = coeffs_out.reshape((N, n_steps + 1), order="F").T

    return times_out, coeffs
