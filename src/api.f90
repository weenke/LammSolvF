! =============================================================================
!  api.f90  —  Explicit C API for the SAWP solver
!
!  This module is the single translation unit that defines every symbol
!  visible to the shared library consumer (Python/ctypes, C, Julia, …).
!  All procedures carry bind(C) and use only iso_c_binding types so the
!  ABI is well-defined and stable across compilers.
!
!  Build as a shared library (link against the rest of the project):
!    gfortran -O3 -march=native -std=f2008 -fPIC \
!        schuck_matrices.f90 sawp_solver.f90 api.f90 \
!        -shared -o libsawp.so -llapack -lblas -lgfortran
!
!  With fpm, add to fpm.toml:
!    [library]
!    type = "shared"
!
!  ── Entry points ────────────────────────────────────────────────────────────
!
!  1. sawp_solver              (re-exported from sawp_solver_mod)
!       int sawp_solver(const double *c0, double t0, int n_steps,
!                       double omega, double s, double D,
!                       double m, double b, int N, double dt,
!                       double *times_out, double *coeffs_out);
!
!  2. sawp_grid_init           (C wrapper for pure grid_init)
!       void sawp_grid_init(double m, double b, int N, double *r0_out);
!
!  3. sawp_grid_propagate      (C wrapper for pure grid_propagate)
!       void sawp_grid_propagate(double t, double t0, const double *r0,
!                                double s_G, double omega,
!                                double m, double b, int N,
!                                double *rkt_out);
!
!  4. hat_basis_c              (re-exported from schuck_matrices_mod)
!       void hat_basis_c(int n_nodes, int n_r,
!                        const double *nodes, const double *r,
!                        double *P);
!
!  5. build_schuck_diags_c     (re-exported from schuck_matrices_mod)
!       void build_schuck_diags_c(int n, const double *nodes_t,
!                                 double omega, double s, double s_G,
!                                 double D, double alpha,
!                                 double *B_sub, double *B_diag,
!                                 double *B_sup,
!                                 double *J_sub, double *J_diag,
!                                 double *J_sup);
!
!  6. sawp_query_sizes         (buffer-size helper for Python callers)
!       void sawp_query_sizes(int N, int n_steps,
!                             int *times_len, int *coeffs_len);
!
!  ── ctypes usage (Python) ───────────────────────────────────────────────────
!
!    import ctypes, numpy as np
!
!    lib = ctypes.CDLL("./libsawp.so")
!
!    # --- query buffer sizes -------------------------------------------------
!    lib.sawp_query_sizes.restype  = None
!    lib.sawp_query_sizes.argtypes = [
!        ctypes.c_int, ctypes.c_int,
!        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
!    ]
!    t_len, c_len = ctypes.c_int(), ctypes.c_int()
!    lib.sawp_query_sizes(N, n_steps, ctypes.byref(t_len), ctypes.byref(c_len))
!
!    # --- allocate output buffers --------------------------------------------
!    times_out  = np.zeros(t_len.value,  dtype=np.float64)
!    coeffs_out = np.zeros(c_len.value,  dtype=np.float64)
!
!    # --- call the solver ----------------------------------------------------
!    lib.sawp_solver.restype  = ctypes.c_int
!    lib.sawp_solver.argtypes = [
!        ctypes.POINTER(ctypes.c_double),  # c0
!        ctypes.c_double,                  # t0
!        ctypes.c_int,                     # n_steps
!        ctypes.c_double,                  # omega
!        ctypes.c_double,                  # s
!        ctypes.c_double,                  # D
!        ctypes.c_double,                  # m
!        ctypes.c_double,                  # b
!        ctypes.c_int,                     # N
!        ctypes.c_double,                  # dt
!        ctypes.POINTER(ctypes.c_double),  # times_out
!        ctypes.POINTER(ctypes.c_double),  # coeffs_out
!    ]
!    info = lib.sawp_solver(
!        c0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
!        ctypes.c_double(t0), ctypes.c_int(n_steps),
!        ctypes.c_double(omega), ctypes.c_double(s),
!        ctypes.c_double(D), ctypes.c_double(m), ctypes.c_double(b),
!        ctypes.c_int(N), ctypes.c_double(dt),
!        times_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
!        coeffs_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
!    )
!    # coeffs_out is Fortran column-major (N, n_steps+1);
!    # reshape to (N, n_steps+1) then transpose to (n_steps+1, N):
!    coeffs = coeffs_out.reshape((N, n_steps + 1), order='F').T
!
! =============================================================================
module sawp_api_mod
  use iso_c_binding,       only: c_double, c_int
  use sawp_solver_mod,     only: sawp_solver
  use solvf_utilities,     only: grid_init, grid_propagate, &
                                 hat_basis_c, build_schuck_diags_c
  implicit none
  private
 
  ! Re-export the already-bind(C) routines so link-time resolution finds them
  ! in this compilation unit when the shared library is built.
  public :: sawp_solver
  public :: hat_basis_c
  public :: build_schuck_diags_c
 
  ! New C entry points defined in this module
  public :: sawp_grid_init
  public :: sawp_grid_propagate
  public :: sawp_query_sizes
 
contains
 
  ! ===========================================================================
  !  sawp_grid_init
  !
  !  C wrapper for the pure Fortran subroutine grid_init.
  !  Pure procedures cannot carry bind(C) in standard Fortran, so this thin
  !  wrapper calls grid_init and copies the result into the caller's buffer.
  !
  !  C signature:
  !    void sawp_grid_init(double m, double b, int N, double *r0_out);
  !
  !  Arguments:
  !    m      [in]  : meniscus radius (lower boundary)
  !    b      [in]  : bottom radius   (upper boundary)
  !    N      [in]  : total number of grid points
  !    r0_out [out] : caller-allocated buffer of length N; receives the grid
  ! ===========================================================================
  subroutine sawp_grid_init(m, b, N, r0_out) bind(C, name="sawp_grid_init")
    real(c_double), value, intent(in)  :: m, b
    integer(c_int), value, intent(in)  :: N
    real(c_double),        intent(out) :: r0_out(N)
 
    call grid_init(m, b, N, r0_out)
  end subroutine sawp_grid_init
 
 
  ! ===========================================================================
  !  sawp_grid_propagate
  !
  !  C wrapper for the pure Fortran subroutine grid_propagate.
  !
  !  C signature:
  !    void sawp_grid_propagate(double t, double t0,
  !                             const double *r0,
  !                             double s_G, double omega,
  !                             double m, double b, int N,
  !                             double *rkt_out);
  !
  !  Arguments:
  !    t       [in]  : current time
  !    t0      [in]  : reference (initial) time
  !    r0      [in]  : initial grid of length N (from sawp_grid_init)
  !    s_G     [in]  : shear-flow parameter
  !    omega   [in]  : angular frequency
  !    m       [in]  : meniscus boundary (pinned)
  !    b       [in]  : bottom boundary   (pinned)
  !    N       [in]  : total number of grid points
  !    rkt_out [out] : caller-allocated buffer of length N; receives r_k(t)
  ! ===========================================================================
  subroutine sawp_grid_propagate(t, t0, r0, s_G, omega, m, b, N, rkt_out) &
      bind(C, name="sawp_grid_propagate")
    real(c_double), value, intent(in)  :: t, t0
    real(c_double),        intent(in)  :: r0(N)
    real(c_double), value, intent(in)  :: s_G, omega, m, b
    integer(c_int), value, intent(in)  :: N
    real(c_double),        intent(out) :: rkt_out(N)
 
    call grid_propagate(t, t0, r0, s_G, omega, m, b, N, rkt_out)
  end subroutine sawp_grid_propagate
 
 
  ! ===========================================================================
  !  sawp_query_sizes
  !
  !  Buffer-size query helper.  The Python caller must allocate times_out and
  !  coeffs_out before calling sawp_solver.  This routine returns the required
  !  flat lengths so the caller does not embed the size formulas directly.
  !
  !  C signature:
  !    void sawp_query_sizes(int N, int n_steps,
  !                          int *times_len, int *coeffs_len);
  !
  !  Arguments:
  !    N          [in]  : number of grid points
  !    n_steps    [in]  : number of time steps
  !    times_len  [out] : required length of times_out  = n_steps + 1
  !    coeffs_len [out] : required length of coeffs_out = N * (n_steps + 1)
  !                       flat, column-major; Python reshapes with order='F'
  !                       then transposes to get row-major (n_steps+1, N)
  ! ===========================================================================
  subroutine sawp_query_sizes(N, n_steps, times_len, coeffs_len) &
      bind(C, name="sawp_query_sizes")
    integer(c_int), value, intent(in)  :: N, n_steps
    integer(c_int),        intent(out) :: times_len, coeffs_len
 
    times_len  = n_steps + 1_c_int
    coeffs_len = N * (n_steps + 1_c_int)
  end subroutine sawp_query_sizes
 
end module sawp_api_mod