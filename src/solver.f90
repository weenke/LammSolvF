! =============================================================================
!  sawp_solver.f90
!
!  Main solver routine for the SAWP (Schuck Adaptive Wave Propagation) model.
!  Implements run_LammSim_dtOPT in Fortran with full C interoperability.
!
!  OPTIMISATION — tridiagonal storage throughout:
!    B and J are symmetric tridiagonal matrices.  B_star, lhs, and rhs_mat
!    are therefore also tridiagonal.  This implementation stores all matrices
!    as three 1-D diagonal arrays and replaces every dense LAPACK/BLAS call
!    with its O(N) tridiagonal equivalent:
!
!    Dense (old)                    Tridiagonal (new)        Cost
!    ----------------------------- ------------------------  --------
!    dgetrf  (dense LU)             dgttrf  (tridiag LU)     O(N³)->O(N)
!    dgetrs  (dense solve)          dgttrs  (tridiag solve)  O(N²)->O(N)
!    dgemv   (dense matvec)         tri_matvec (3 saxpy)     O(N²)->O(N)
!    ddot    (mass balance)         ddot  (unchanged)        O(N)
!
!    Storage: 7×N² doubles -> 7×3×N doubles  (eliminated entirely for large N)
!
!  C signature:
!    int sawp_solver(
!        const double *c0,
!        double        t0,
!        int           n_steps,
!        double        omega,
!        double        s,
!        double        D,
!        double        m,
!        double        b,
!        int           N,
!        double        dt,
!        double       *times_out,    // length  n_steps+1
!        double       *coeffs_out    // length (n_steps+1)*N, row-major
!    );
!
!  Return value: 0 on success, non-zero LAPACK info on factorisation failure.
! =============================================================================
module sawp_solver_mod
  use iso_c_binding,       only: c_double, c_int
  use solvf_utilities, only: grid_init, grid_propagate, &
                                 build_schuck_diags_c
  implicit none
  private
  public :: sawp_solver
 
  ! ---------------------------------------------------------------------------
  !  LAPACK/BLAS explicit interfaces
  !  Resolved at link time against -llapack -lblas / OpenBLAS / MKL.
  ! ---------------------------------------------------------------------------
  interface
    ! Tridiagonal LU factorisation
    ! dl(N-1), d(N), du(N-1) overwritten with LU factors in-place.
    ! du2(N-2): second super-diagonal created by partial pivoting.
    subroutine dgttrf(n, dl, d, du, du2, ipiv, info)
      import
      integer,        intent(in)    :: n
      real(c_double), intent(inout) :: dl(*)    ! sub-diagonal,   length n-1
      real(c_double), intent(inout) :: d(*)     ! main diagonal,  length n
      real(c_double), intent(inout) :: du(*)    ! super-diagonal, length n-1
      real(c_double), intent(out)   :: du2(*)   ! 2nd super-diag, length n-2
      integer,        intent(out)   :: ipiv(*)  ! pivot indices,  length n
      integer,        intent(out)   :: info
    end subroutine dgttrf
 
    ! Tridiagonal solve using factors from dgttrf
    subroutine dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
      import
      character(len=1), intent(in)    :: trans
      integer,          intent(in)    :: n, nrhs, ldb
      real(c_double),   intent(in)    :: dl(*), d(*), du(*), du2(*)
      integer,          intent(in)    :: ipiv(*)
      real(c_double),   intent(inout) :: B(ldb, *)
      integer,          intent(out)   :: info
    end subroutine dgttrs
 
    ! Dense dot product — kept for O(N) mass balance
    function ddot(n, x, incx, y, incy) result(res)
      import
      integer,        intent(in) :: n, incx, incy
      real(c_double), intent(in) :: x(*), y(*)
      real(c_double)             :: res
    end function ddot
  end interface
 
contains
 
  ! ===========================================================================
  !  SAWP SOLVER
  !
  !  All matrices stored as tridiagonal (sub, diag, super) 1-D arrays.
  !  The time loop performs only O(N) work per step.
  !
  !  Indexing (Python 0-based -> Fortran 1-based):
  !    c_next[2:-1]   = c_interp[1:-2]  ->  c_next(3:N-1) = c_interp(2:N-2)
  !    c_next[0]      = c_new[0]         ->  c_next(1)     = c_new(1)
  !    c_next[1]      = 0.5*(c[0]+c[2]) ->  c_next(2)     = 0.5*(c(1)+c(3))
  !    c_next[-1]     = mass_balance     ->  c_next(N)     = mass_balance
  ! ===========================================================================
  function sawp_solver(c0, t0, n_steps, omega, s, D, m, b, N, dt, &
                       times_out, coeffs_out)                       &
      result(info) bind(C, name="sawp_solver")
 
    ! --- Arguments -----------------------------------------------------------
    real(c_double),        intent(in)  :: c0(N)
    real(c_double), value, intent(in)  :: t0
    integer(c_int), value, intent(in)  :: n_steps
    real(c_double), value, intent(in)  :: omega, s, D, m, b
    integer(c_int), value, intent(in)  :: N
    real(c_double), value, intent(in)  :: dt
    real(c_double),        intent(out) :: times_out(n_steps + 1)
    real(c_double),        intent(out) :: coeffs_out(N, n_steps + 1)
 
    integer(c_int) :: info
 
    ! --- Local scalars -------------------------------------------------------
    real(c_double) :: s_G, alpha_dt, alpha_1
    real(c_double) :: initial_mass, current_mass_sum, inv_rad_w_last, denom
    integer        :: step, i, lapack_info
 
    ! --- Tridiagonal storage for B(t0) and J(t0) ----------------------------
    real(c_double) :: B0_sub(N-1), B0_diag(N), B0_sup(N-1)
    real(c_double) :: J0_sub(N-1), J0_diag(N), J0_sup(N-1)
 
    ! --- Tridiagonal storage for B(t0+dt) and J(t0+dt) ---------------------
    real(c_double) :: B1_sub(N-1), B1_diag(N), B1_sup(N-1)
    real(c_double) :: J1_sub(N-1), J1_diag(N), J1_sup(N-1)
 
    ! --- Time-invariant tridiagonal combinations (eq. 14) -------------------
    !  B_star  = B0 + B1
    !  lhs     = B_star - dt*J1   (factorised in-place by dgttrf)
    !  rhs_mat = B_star + dt*J0
    real(c_double) :: lhs_sub(N-1),  lhs_diag(N),  lhs_sup(N-1)
    real(c_double) :: rhs_sub(N-1),  rhs_diag(N),  rhs_sup(N-1)
 
    ! --- dgttrf output: LU factors overwrite lhs_*, du2 holds 2nd super-diag
    real(c_double) :: lhs_du2(N-2)
    integer        :: ipiv(N)
 
    ! --- Working vectors -----------------------------------------------------
    real(c_double) :: b_vec(N)      ! rhs_mat * c_curr  (tridiagonal matvec)
    real(c_double) :: c_new(N)      ! solution from dgttrs
    real(c_double) :: c_interp(N)   ! interpolated c_new onto nodes_t0 grid
 
    ! --- Interpolation precomputation ----------------------------------------
    integer        :: idx(N)
    real(c_double) :: w_L(N), w_R(N)
 
    ! --- Mass balance weights ------------------------------------------------
    real(c_double) :: radial_weights(N)
 
    ! --- Grid arrays ---------------------------------------------------------
    real(c_double) :: nodes_t0(N), nodes_tpdt(N)
 
    ! =========================================================================
    !  SETUP
    ! =========================================================================
    info     = 0_c_int
    s_G      = s
    alpha_dt = exp(s_G * omega**2 * dt)
    alpha_1  = 1.0_c_double
 
    ! --- 1. Build grids ------------------------------------------------------
    call grid_init(m, b, N, nodes_t0)
    call grid_propagate(t0 + dt, t0, nodes_t0, s_G, omega, m, b, N, nodes_tpdt)
 
    ! --- 2. Build tridiagonal diagonals for B and J at both time levels ------
    !  build_schuck_diags_c returns diagonals directly — no N×N matrix built.
    call build_schuck_diags_c(N, nodes_t0,   omega, s, s_G, D, alpha_1,  &
                               B0_sub, B0_diag, B0_sup,                    &
                               J0_sub, J0_diag, J0_sup)
    call build_schuck_diags_c(N, nodes_tpdt, omega, s, s_G, D, alpha_dt, &
                               B1_sub, B1_diag, B1_sup,                    &
                               J1_sub, J1_diag, J1_sup)
 
    ! --- 3. Form lhs and rhs_mat diagonals (eq. 14) -------------------------
    !  B_star diagonals (reused in both lhs and rhs so computed once inline):
    !  lhs  = (B0+B1) - dt*J1
    !  rhs  = (B0+B1) + dt*J0
    lhs_sub  = (B0_sub  + B1_sub)  - dt * J1_sub
    lhs_diag = (B0_diag + B1_diag) - dt * J1_diag
    lhs_sup  = (B0_sup  + B1_sup)  - dt * J1_sup
 
    rhs_sub  = (B0_sub  + B1_sub)  + dt * J0_sub
    rhs_diag = (B0_diag + B1_diag) + dt * J0_diag
    rhs_sup  = (B0_sup  + B1_sup)  + dt * J0_sup
 
    ! --- 4. Tridiagonal LU factorisation of lhs (done exactly once) ---------
    !  dgttrf overwrites lhs_sub, lhs_diag, lhs_sup with LU factors in-place.
    !  lhs_du2 receives the second super-diagonal generated by partial pivoting.
    call dgttrf(N, lhs_sub, lhs_diag, lhs_sup, lhs_du2, ipiv, lapack_info)
    if (lapack_info /= 0) then
      info = int(lapack_info, c_int)
      return
    end if
 
    ! --- 5. Precompute interpolation weights (O(N) linear scan) --------------
    block
      integer :: k
      k = 2
      do i = 1, N
        do while (k < N .and. nodes_tpdt(k) < nodes_t0(i))
          k = k + 1
        end do
        idx(i) = k
      end do
    end block
 
    idx = max(2, min(N, idx))
 
    do i = 1, N
      denom  = nodes_tpdt(idx(i)) - nodes_tpdt(idx(i) - 1)
      w_R(i) = (nodes_t0(i) - nodes_tpdt(idx(i) - 1)) / denom
      w_L(i) = 1.0_c_double - w_R(i)
      if (nodes_t0(i) < nodes_tpdt(1)) then
        w_L(i) = 1.0_c_double ;  w_R(i) = 0.0_c_double ;  idx(i) = 2
      else if (nodes_t0(i) > nodes_tpdt(N)) then
        w_L(i) = 0.0_c_double ;  w_R(i) = 1.0_c_double ;  idx(i) = N
      end if
    end do
 
    ! --- 6. Radial (mass balance) weights = column sums of B0 ---------------
    !  B0 is symmetric tridiagonal so column sums = row sums:
    !    radial_weights(1)   = B0_diag(1) + B0_sup(1)
    !    radial_weights(k)   = B0_sub(k-1) + B0_diag(k) + B0_sup(k)  k=2..N-1
    !    radial_weights(N)   = B0_sub(N-1) + B0_diag(N)
    !  No N×N matrix needed — computed directly from the diagonal arrays.
    radial_weights(1)     = B0_diag(1) + B0_sup(1)
    radial_weights(2:N-1) = B0_sub(1:N-2) + B0_diag(2:N-1) + B0_sup(2:N-1)
    radial_weights(N)     = B0_sub(N-1) + B0_diag(N)
 
    initial_mass   = ddot(N, c0, 1, radial_weights, 1)
    inv_rad_w_last = 1.0_c_double / radial_weights(N)
 
    ! --- 7. Initialise output ------------------------------------------------
    do step = 1, n_steps + 1
      times_out(step) = t0 + real(step - 1, c_double) * dt
    end do
    coeffs_out(:, 1) = c0
 
    ! =========================================================================
    !  TIME LOOP  — O(N) per step
    ! =========================================================================
    do step = 1, n_steps
 
      ! --- (a) Tridiagonal matvec: b_vec = rhs_mat * c_curr -----------------
      !  Replaces dgemv (O(N²)) with three array-section operations (O(N)).
      !  y = D*x + DL*x_shift_right + DU*x_shift_left
      b_vec = rhs_diag * coeffs_out(:, step)
      b_vec(2:N)   = b_vec(2:N)   + rhs_sub * coeffs_out(1:N-1, step)
      b_vec(1:N-1) = b_vec(1:N-1) + rhs_sup * coeffs_out(2:N,   step)
 
      ! --- (b) Tridiagonal solve: c_new = lhs^{-1} * b_vec -----------------
      !  dgttrs uses the LU factors stored in lhs_sub/diag/sup/du2 and ipiv.
      !  b_vec is copied first since dgttrs overwrites its RHS argument.
      c_new = b_vec
      call dgttrs('N', N, 1, lhs_sub, lhs_diag, lhs_sup, lhs_du2, &
                  ipiv, c_new, N, lapack_info)
      if (lapack_info /= 0) then
        info = int(lapack_info, c_int)
        return
      end if
 
      ! --- (c) Interpolate c_new onto nodes_t0 grid ------------------------
      do i = 1, N
        c_interp(i) = c_new(idx(i) - 1) * w_L(i) + c_new(idx(i)) * w_R(i)
      end do
 
      ! --- (d) Assemble c_next ----------------------------------------------
      associate(c_next => coeffs_out(:, step + 1))
        c_next(3:N-1) = c_interp(2:N-2)
        c_next(1)     = c_new(1)
        c_next(2)     = 0.5_c_double * (c_next(1) + c_next(3))
 
        ! --- (e) Mass balance ------------------------------------------------
        current_mass_sum = ddot(N - 1, c_next(1), 1, radial_weights(1), 1)
        c_next(N) = (initial_mass - current_mass_sum) * inv_rad_w_last
      end associate
 
    end do
 
  end function sawp_solver
 
end module sawp_solver_mod