! Public C API for the SAWP solver.
!
! Entry points:
!   int  sawp_solver(...)            — time integration loop
!   void sawp_grid_init(...)         — log-spaced grid initialisation
!   void sawp_grid_propagate(...)    — advance grid one time step
!   void hat_basis_c(...)            — piecewise-linear hat basis
!   void build_schuck_diags_c(...)   — tridiagonal B and J diagonals
!   void sawp_query_sizes(...)       — output buffer lengths
module sawp_api_mod
  use iso_c_binding,   only: c_double, c_int
  use sawp_solver_mod, only: sawp_solver
  use solvf_utilities, only: grid_init, grid_propagate, &
                              hat_basis_c, build_schuck_diags_c
  implicit none
  private

  public :: sawp_solver
  public :: hat_basis_c
  public :: build_schuck_diags_c
  public :: sawp_grid_init
  public :: sawp_grid_propagate
  public :: sawp_query_sizes

contains

  ! C signature:
  !   void sawp_grid_init(double m, double b, int N, double *r0_out);
  !
  ! Arguments:
  !   m      [in]  : meniscus radius (lower boundary)
  !   b      [in]  : bottom radius   (upper boundary)
  !   N      [in]  : number of grid points
  !   r0_out [out] : log-spaced grid of length N
  subroutine sawp_grid_init(m, b, N, r0_out) bind(C, name="sawp_grid_init")
    real(c_double), value, intent(in)  :: m, b
    integer(c_int), value, intent(in)  :: N
    real(c_double),        intent(out) :: r0_out(N)

    call grid_init(m, b, N, r0_out)
  end subroutine sawp_grid_init


  ! C signature:
  !   void sawp_grid_propagate(double t, double t0, const double *r0,
  !                            double s_G, double omega,
  !                            double m, double b, int N,
  !                            double *rkt_out);
  !
  ! Arguments:
  !   t       [in]  : current time
  !   t0      [in]  : reference time
  !   r0      [in]  : initial grid of length N (from sawp_grid_init)
  !   s_G     [in]  : shear-flow parameter
  !   omega   [in]  : angular frequency
  !   m       [in]  : meniscus boundary (pinned)
  !   b       [in]  : bottom boundary   (pinned)
  !   N       [in]  : number of grid points
  !   rkt_out [out] : propagated grid of length N
  subroutine sawp_grid_propagate(t, t0, r0, s_G, omega, m, b, N, rkt_out) &
      bind(C, name="sawp_grid_propagate")
    real(c_double), value, intent(in)  :: t, t0
    real(c_double),        intent(in)  :: r0(N)
    real(c_double), value, intent(in)  :: s_G, omega, m, b
    integer(c_int), value, intent(in)  :: N
    real(c_double),        intent(out) :: rkt_out(N)

    call grid_propagate(t, t0, r0, s_G, omega, m, b, N, rkt_out)
  end subroutine sawp_grid_propagate


  ! C signature:
  !   void sawp_query_sizes(int N, int n_steps,
  !                         int *times_len, int *coeffs_len);
  !
  ! Arguments:
  !   N          [in]  : number of grid points
  !   n_steps    [in]  : number of time steps
  !   times_len  [out] : n_steps + 1
  !   coeffs_len [out] : N * (n_steps + 1)
  subroutine sawp_query_sizes(N, n_steps, times_len, coeffs_len) &
      bind(C, name="sawp_query_sizes")
    integer(c_int), value, intent(in)  :: N, n_steps
    integer(c_int),        intent(out) :: times_len, coeffs_len

    times_len  = n_steps + 1_c_int
    coeffs_len = N * (n_steps + 1_c_int)
  end subroutine sawp_query_sizes

end module sawp_api_mod
