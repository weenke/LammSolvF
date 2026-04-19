"""LammSolvF — Python interface to the SAWP Fortran solver."""
from importlib import metadata

from lammsolvflib.solver import sawp_solver

__all__ = ["sawp_solver"]

try:
    __version__ = metadata.version("lammsolvflib")
except metadata.PackageNotFoundError:
    __version__ = "unknown"
