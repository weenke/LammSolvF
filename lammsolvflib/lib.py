"""Low-level ctypes bridge to the LammSolvF shared library."""
import ctypes
from ctypes import POINTER, c_double, c_int
import os
from pathlib import Path
import sys
from typing import Optional

import numpy as np
from numpy.ctypeslib import ndpointer

# Load shared library from the package directory (installed by CMake/scikit-build)
path = Path(__file__).parent
name = "libLammSolvF"
winmode: Optional[int] = None
if sys.platform == "win32":
    name = f"{name}.dll"
    os.add_dll_directory(str(path))
    winmode = 0
elif sys.platform.startswith("linux"):
    name = f"{name}.so"
elif sys.platform == "darwin":
    name = f"{name}.dylib"
else:
    raise ImportError(f"Unsupported platform: {sys.platform}")

lib = ctypes.CDLL(str(path / name), winmode=winmode)

# ---------------------------------------------------------------------------
#  sawp_query_sizes
#    void sawp_query_sizes(int N, int n_steps,
#                          int *times_len, int *coeffs_len);
#
#  Must be called before sawp_solver to obtain the flat buffer lengths needed
#  for times_out and coeffs_out.
# ---------------------------------------------------------------------------
lib.sawp_query_sizes.restype = None
lib.sawp_query_sizes.argtypes = [
    c_int,           # N
    c_int,           # n_steps
    POINTER(c_int),  # times_len  (out) = n_steps + 1
    POINTER(c_int),  # coeffs_len (out) = N * (n_steps + 1)
]

# ---------------------------------------------------------------------------
#  sawp_solver
#    int sawp_solver(const double *c0, double t0, int n_steps,
#                    double omega, double s, double D,
#                    double m, double b, int N, double dt,
#                    double *times_out, double *coeffs_out);
#
#  Returns 0 on success; non-zero LAPACK info code on factorisation failure.
#  coeffs_out is written in Fortran column-major order (shape N × (n_steps+1)).
#  Reshape the flat buffer with .reshape((N, n_steps+1), order='F').T to get
#  a row-major array of shape (n_steps+1, N).
# ---------------------------------------------------------------------------
lib.sawp_solver.restype = c_int
lib.sawp_solver.argtypes = [
    ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # c0[N]
    c_double,                                                    # t0
    c_int,                                                       # n_steps
    c_double,                                                    # omega
    c_double,                                                    # s
    c_double,                                                    # D
    c_double,                                                    # m
    c_double,                                                    # b
    c_int,                                                       # N
    c_double,                                                    # dt
    ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # times_out[n_steps+1]
    ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # coeffs_out[N*(n_steps+1)]
]
