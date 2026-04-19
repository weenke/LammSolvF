module solvf_utilities
    use iso_c_binding, only: c_double, c_int
    use stdlib_linalg, only: diag
    use stdlib_sorting, only: sort_index, int_index
    implicit none
    private

    public :: hat_basis_c
    public :: build_schuck_matrices_c 
    public :: build_schuck_diags_c
    public :: grid_init
    public :: grid_propagate
    public :: d_srf

contains

! ===========================================================================
!- Hat basis generation
!  C signature:
!    void hat_basis_c(int n_nodes, int n_r,
!                     const double *nodes, const double *r,
!                     double *P);
!
!  P is stored column-major, shape (n_nodes, n_r).
!  In C the caller sees a flat array of length n_nodes * n_r.
! ===========================================================================

    subroutine hat_basis_c(n_nodes, n_r, nodes, r, P) bind(C, name="hat_basis_c")
        integer(c_int), value, intent(in) :: n_nodes, n_r
        real(c_double), intent(in)        :: nodes(n_nodes)
        real(c_double), intent(in)        :: r(n_r)
        real(c_double), intent(out)       :: P(n_nodes, n_r)

        ! sort_index takes an input array, working on a copy
        real(c_double)                    :: r_sorted(n_r)
        integer(int_index)                :: order(n_r)
        integer                           :: i, orig_j, node_k
        real(c_double)                    :: rj, r_L, r_R, denom

        P = 0.0_c_double

        !Step 1: sort query points, keep track of original positions
        r_sorted = r
        call sort_index(r_sorted, order) ! r_sorted ascending; order maps back

        !Step 2: single linear scan through nodes
        node_k = 2

        do i = 1, n_r
            rj = r_sorted(i)
            orig_j = int(order(i)) ! original column in P

            !Skip pts outside the grid
            if (rj < nodes(1) .or. rj > nodes(n_nodes)) cycle

            !Advance the bracket pointer (nodes is sorted, r_sorted is sorted)

            do while (node_k < n_nodes .and. nodes(node_k) < rj)
                node_k = node_k + 1
            end do

            ! Evaluate hat weights for the interval [nodes(node_k-1), nodes(node_k)]
            r_L = nodes(node_k-1)
            r_R = nodes(node_k)
            denom = r_R - r_L
            P(node_k - 1, orig_j) = (r_R - rj) / denom
            P(node_k,     orig_j) = (rj - r_L) / denom
        end do
    end subroutine hat_basis_c

    ! ===========================================================================
    !  Building Schuck matrices (Moving-grid implementation of Schuck 1998)
    !  C signature:
    !    void build_schuck_matrices_c(int n,
    !                                 const double *nodes_t,
    !                                 double omega, double s, double s_G,
    !                                 double D,     double alpha,
    !                                 double *B_t,  double *J);
    !
    !  B_t and J are dense (n x n) arrays, col-major, length n*n each.
    ! ===========================================================================
    
    subroutine build_schuck_matrices_c(n, nodes_t, omega, s, s_G, D, alpha, &
                                        B_t, J) bind(C, name="build_schuck_matrices_c")
        integer(c_int), value, intent(in)   :: n
        real(c_double),        intent(in)   :: nodes_t(n)
        real(c_double), value, intent(in)   :: omega, s, s_G, D, alpha
        real(c_double),        intent(out)  :: B_t(n,n)
        real(c_double),        intent(out)  :: J(n,n)

        ! --- local diagonal arrays -----------------------------------------------
        real(c_double) :: B_sub(n-1),  B_diag(n),  B_sup(n-1)
        real(c_double) :: A1_sub(n-1), A1_diag(n), A1_sup(n-1)
        real(c_double) :: A2_sub(n-1), A2_diag(n), A2_sup(n-1)
        real(c_double) :: A3_sub(n-1), A3_diag(n), A3_sup(n-1)
        real(c_double) :: J_sub(n-1),  J_diag(n),  J_sup(n-1)

        real(c_double) :: r2(n)

        real(c_double) :: c1, c2, c3

        r2 = nodes_t * nodes_t ! single vectorised multiply, reused 4 times

        ! 1. Fill the four sets of tridiagonal diagonals
        call get_B_diags (nodes_t, r2, n, B_sub,  B_diag,  B_sup )
        call get_A1_diags(nodes_t, r2, n, A1_sub, A1_diag, A1_sup)
        call get_A2_diags(nodes_t, r2, n, A2_sub, A2_diag, A2_sup)
        call get_A3_diags(nodes_t, r2, n, A3_sub, A3_diag, A3_sup)

        ! 2. Scale factors
        c1 = D * alpha**(-2)
        c2 = omega**2 * s
        c3 = omega**2 * s_G

        ! 3. Combine into J diagonals (pure array arithmetic on 1-D arrays)
        J_sub  = c2*A2_sub  - c3*A3_sub  - c1*A1_sub
        J_diag = c2*A2_diag - c3*A3_diag - c1*A1_diag
        J_sup  = c2*A2_sup  - c3*A3_sup  - c1*A1_sup

        B_t = diag(B_sub, -1) + diag(B_diag) + diag(B_sup, 1)
        J   = diag(J_sub, -1) + diag(J_diag) + diag(J_sup, 1)
    end subroutine build_schuck_matrices_c

    ! ===========================================================================
    !  BUILD SCHUCK DIAGONALS  (tridiagonal-only variant)
    !
    !  Identical mathematics to build_schuck_matrices_c but returns the three
    !  diagonals of B and J directly instead of assembling N×N dense matrices.
    !  This eliminates O(N²) storage and the diag() assembly cost entirely.
    !  The solver uses these arrays with dgttrf/dgttrs and a hand-coded
    !  tridiagonal matvec, reducing every hot-path operation from O(N²) to O(N).
    !
    !  C signature:
    !    void build_schuck_diags_c(int n,
    !                              const double *nodes_t,
    !                              double omega, double s, double s_G,
    !                              double D,     double alpha,
    !                              double *B_sub,  double *B_diag, double *B_sup,
    !                              double *J_sub,  double *J_diag, double *J_sup);
    !
    !  B_sub, B_sup, J_sub, J_sup have length n-1.
    !  B_diag, J_diag have length n.
    ! ===========================================================================
    subroutine build_schuck_diags_c(n, nodes_t, omega, s, s_G, D, alpha, &
                                    B_sub, B_diag, B_sup,                &
                                    J_sub, J_diag, J_sup)                &
        bind(C, name="build_schuck_diags_c")
    integer(c_int), value, intent(in)  :: n
    real(c_double),        intent(in)  :: nodes_t(n)
    real(c_double), value, intent(in)  :: omega, s, s_G, D, alpha
    real(c_double),        intent(out) :: B_sub(n-1), B_diag(n), B_sup(n-1)
    real(c_double),        intent(out) :: J_sub(n-1), J_diag(n), J_sup(n-1)

    real(c_double) :: A1_sub(n-1), A1_diag(n), A1_sup(n-1)
    real(c_double) :: A2_sub(n-1), A2_diag(n), A2_sup(n-1)
    real(c_double) :: A3_sub(n-1), A3_diag(n), A3_sup(n-1)
    real(c_double) :: r2(n)
    real(c_double) :: c1, c2, c3

    r2 = nodes_t * nodes_t

    call get_B_diags (nodes_t, r2, n, B_sub,  B_diag,  B_sup )
    call get_A1_diags(nodes_t, r2, n, A1_sub, A1_diag, A1_sup)
    call get_A2_diags(nodes_t, r2, n, A2_sub, A2_diag, A2_sup)
    call get_A3_diags(nodes_t, r2, n, A3_sub, A3_diag, A3_sup)

    c1 = D * alpha**(-2)
    c2 = omega**2 * s
    c3 = omega**2 * s_G

    J_sub  = c2*A2_sub  - c3*A3_sub  - c1*A1_sub
    J_diag = c2*A2_diag - c3*A3_diag - c1*A1_diag
    J_sup  = c2*A2_sup  - c3*A3_sup  - c1*A1_sup

    end subroutine build_schuck_diags_c

    ! ===========================================================================
    !  t0,k log-space grid initialization
    ! Computes the SAWP starting grid of length N:
    !
    !    r(1)   = m                                          (meniscus boundary)
    !    r(k)   = m * (b/m)^((k - 3/2) / (N - 1))   k = 2..N-1  (interior)
    !    r(N)   = b                                          (bottom boundary)
    !
    !  This is a pure Fortran routine (not bind(C)) since it is internal to the
    !  utility module and is not part of the C-interoperable solver interface.
    !
    !  Arguments:
    !    m    [in]  : meniscus radius (lower boundary), r_1 = m
    !    b    [in]  : bottom radius   (upper boundary), r_N = b
    !    N    [in]  : total number of grid points (including both boundaries)
    !    r0   [out] : grid array of size N; must be allocated by the caller

    pure subroutine grid_init(m, b, N, r0)
        real(c_double), intent(in)  :: m, b
        integer(c_int), intent(in)  :: N
        real(c_double), intent(out) :: r0(N)
    
        ! k_interior holds the integer indices k = 2, 3, ..., N-1 as reals,
        ! matching np.arange(2, N) — exactly N-2 values.
        real(c_double) :: k_interior(N-2)
        integer        :: k
    
        ! --- Build the index vector [2.0, 3.0, ..., N-1.0] ----------------------
        ! Equivalent to np.arange(2, N, dtype=float64).
        ! Written as an implied-do rather than stdlib arange() so the routine
        ! can remain pure (stdlib arange is not pure as of stdlib v0.7).
        k_interior = [(real(k, c_double), k = 2, N-1)]
    
        ! --- Interior nodes: vectorised elemental exponentiation -----------------
        ! r(k) = m * (b/m)^((k - 1.5) / (N - 1)),  k = 2..N-1
        r0(2:N-1) = m * (b / m) ** ((k_interior - 1.5_c_double) / real(N-1, c_double))
    
        ! --- Fixed boundary values -----------------------------------------------
        r0(1) = m    ! r_1 = meniscus
        r0(N) = b    ! r_N = bottom
    
    end subroutine grid_init


    ! ===========================================================================
    !  r_t,k log-space grid dt propogation
    !
    !  Propagates interior grid nodes forward in time using the exponential
    !  SAWP update rule, then reassembles the full N-point grid with pinned
    !  boundaries.
    !
    !  Equation implemented:
    !    r_1(t)   = m                                         (pinned)
    !    r_k(t)   = r_k,0 * exp(S_G * omega^2 * (t - t0))   k = 2..N-1
    !    r_N(t)   = b                                         (pinned)
    !
    !  Arguments:
    !    t        [in]  : current time
    !    t0       [in]  : reference (initial) time
    !    r0       [in]  : initial grid of length N  (output of grid_init)
    !    S_G      [in]  : shear-flow parameter
    !    omega    [in]  : angular frequency
    !    m        [in]  : meniscus boundary (pinned lower value), r_1
    !    b        [in]  : bottom boundary   (pinned upper value), r_N
    !    N        [in]  : total number of grid points
    !    rk_t     [out] : propagated grid of length N; caller allocates
    ! ===========================================================================
    pure subroutine grid_propagate(t, t0, r0, S_G, omega, m, b, N, rk_t)
        real(c_double), intent(in)  :: t, t0
        real(c_double), intent(in)  :: r0(N)
        real(c_double), intent(in)  :: S_G, omega, m, b
        integer(c_int), intent(in)  :: N
        real(c_double), intent(out) :: rk_t(N)
    
        ! Scalar growth factor: computed once, broadcast over the interior section.
        real(c_double) :: alpha_t
    
        alpha_t = exp(S_G * omega**2 * (t - t0))
    
        ! --- Propagate all interior nodes in one array expression ----------------
        rk_t(2:N-1) = r0(2:N-1) * alpha_t
    
        ! --- Pin boundaries ------------------------------------------------------
        rk_t(1) = m
        rk_t(N) = b
    
    end subroutine grid_propagate

    ! ===========================================================================
    !  PRIVATE HELPERS
    !  All helpers now accept the precomputed r2(:) = r(:)**2 array to avoid
    !  redundant squaring.  The bisect_left helper has been removed entirely.
    ! ===========================================================================

    ! ---------------------------------------------------------------------------
    !  B matrix diagonals
    ! ---------------------------------------------------------------------------

    subroutine get_B_diags(r, r2, n, sub_d, diag_d, sup_d)
        integer, intent(in) :: n
        real(c_double), intent(in) :: r(n), r2(n) !r2 precomputed
        real(c_double), intent(out) :: sub_d(n-1), diag_d(n), sup_d(n-1)

        real(c_double) :: m_r, b_r
        integer :: k

        m_r = r(1); b_r = r(n)

        ! Off-diagonals: (r_{k+1}^2 - r_k^2) / 12  — use precomputed r2
        do k = 1, n-1
            sup_d(k) = (r2(k+1) - r2(k)) / 12.0_c_double
            sub_d(k) = sup_d(k)
        end do

        ! Interior diagonal
        do k = 2, n-1
            diag_d(k) = (r(k+1) - r(k-1)) * (r(k-1) + 2.0_c_double*r(k) + r(k+1)) &
                        / 12.0_c_double
        end do

        ! Boundary overrides — use r2 for the squared boundary terms
        diag_d(1) = (2.0_c_double*m_r*r(2) - 3.0_c_double*r2(1) + r2(2)) &
                    / 12.0_c_double
        diag_d(n) = (3.0_c_double*r2(n) - 2.0_c_double*b_r*r(n-1) - r2(n-1)) &
                    / 12.0_c_double
        end subroutine get_B_diags

    ! ---------------------------------------------------------------------------
    !  A1 matrix diagonals
    ! ---------------------------------------------------------------------------
    subroutine get_A1_diags(r, r2, n, sub_d, diag_d, sup_d)
        integer,        intent(in)  :: n
        real(c_double), intent(in)  :: r(n), r2(n)
        real(c_double), intent(out) :: sub_d(n-1), diag_d(n), sup_d(n-1)
    
        real(c_double) :: m_r, b_r
        integer :: k
    
        m_r = r(1) ;  b_r = r(n)
    
        ! A1 off-diagonals involve no squares; r2 accepted for interface uniformity
        do k = 1, n-1
        sup_d(k) = 0.5_c_double * (r(k+1) + r(k)) / (r(k) - r(k+1))
        sub_d(k) = sup_d(k)
        end do
    
        do k = 2, n-1
        diag_d(k) = r(k) * (r(k+1) - r(k-1)) &
                    / ((r(k) - r(k-1)) * (r(k+1) - r(k)))
        end do
    
        diag_d(1) = 0.5_c_double * (r(2) + m_r) / (r(2) - m_r)
        diag_d(n) = 0.5_c_double * (b_r + r(n-1)) / (b_r - r(n-1))
    end subroutine get_A1_diags
    
    
    ! ---------------------------------------------------------------------------
    !  A2 matrix diagonals
    ! ---------------------------------------------------------------------------
    subroutine get_A2_diags(r, r2, n, sub_d, diag_d, sup_d)
        integer,        intent(in)  :: n
        real(c_double), intent(in)  :: r(n), r2(n)
        real(c_double), intent(out) :: sub_d(n-1), diag_d(n), sup_d(n-1)
    
        real(c_double) :: m_r, b_r
        integer :: k
    
        m_r = r(1) ;  b_r = r(n)
    
        do k = 1, n-1
        sub_d(k) = (r2(k+1) + 2.0_c_double*r(k+1)*r(k) + 3.0_c_double*r2(k)) &
                    / 12.0_c_double
        sup_d(k) = -(r2(k) + 2.0_c_double*r(k)*r(k+1) + 3.0_c_double*r2(k+1)) &
                    / 12.0_c_double
        end do
    
        do k = 2, n-1
        diag_d(k) = (r(k-1) - r(k+1)) * (r(k-1) + 2.0_c_double*r(k) + r(k+1)) &
                    / 12.0_c_double
        end do
    
        diag_d(1) = (-3.0_c_double*r2(1) - 2.0_c_double*m_r*r(2) - r2(2)) &
                    / 12.0_c_double
        diag_d(n) = ( 3.0_c_double*r2(n) + 2.0_c_double*b_r*r(n-1) - r2(n-1)) &
                    / 12.0_c_double
    end subroutine get_A2_diags
    
    
    ! ---------------------------------------------------------------------------
    !  A3 matrix diagonals
    ! ---------------------------------------------------------------------------
    subroutine get_A3_diags(r, r2, n, sub_d, diag_d, sup_d)
        integer,        intent(in)  :: n
        real(c_double), intent(in)  :: r(n), r2(n)
        real(c_double), intent(out) :: sub_d(n-1), diag_d(n), sup_d(n-1)
    
        real(c_double) :: m_r, b_r
        integer :: k
    
        m_r = r(1) ;  b_r = r(n)
    
        ! General off-diagonals
        do k = 1, n-1
        sub_d(k) = ( 3.0_c_double*r2(k+1) + 2.0_c_double*r(k+1)*r(k) + r2(k)) &
                    / 12.0_c_double
        sup_d(k) = (-3.0_c_double*r2(k)   - 2.0_c_double*r(k)*r(k+1) - r2(k+1)) &
                    / 12.0_c_double
        end do
    
        ! General interior diagonal
        diag_d(1) = 0.0_c_double   ! will be overwritten below
        diag_d(n) = 0.0_c_double   ! will be overwritten below
        do k = 2, n-1
        diag_d(k) = (r(k+1) - r(k-1)) * (r(k-1) + 2.0_c_double*r(k) + r(k+1)) &
                    / 12.0_c_double
        end do
    
        ! Boundary overwrites
        diag_d(1)  =  r(2) * (m_r + r(2))              / 12.0_c_double
        sup_d(1)   = -r(2) * (m_r + r(2))              / 12.0_c_double
        sub_d(1)   =  r(2) * (m_r + 3.0_c_double*r(2)) / 12.0_c_double
    
        if (n >= 3) then
        diag_d(2)   = (-m_r*r(2) + 2.0_c_double*r(2)*r(3) + r2(3)) / 12.0_c_double
        diag_d(n-1) = ( b_r*r(n-1) - 2.0_c_double*r(n-1)*r(n-2) - r2(n-2)) &
                        / 12.0_c_double
        end if
    
        sup_d(n-1) = -r(n-1) * (b_r + 3.0_c_double*r(n-1)) / 12.0_c_double
        sub_d(n-1) =  r(n-1) * (b_r + r(n-1))              / 12.0_c_double
        diag_d(n)  = -r(n-1) * (b_r + r(n-1))              / 12.0_c_double
    
    end subroutine get_A3_diags

    ! ===========================================================================
    !  Diffusion coefficient from sedimentation coefficient and frictional ratio
    !
    !  Implements the Schuck (1998) relation:
    !    D = (sqrt(2) / (18*pi)) * kB*T * s^(-1/2) * (eta*fr)^(-3/2)
    !        * sqrt((1 - vbar*rho) / vbar)
    !
    !  Arguments:
    !    s   [in]           : sedimentation coefficient (S)
    !    fr  [in]           : frictional ratio (dimensionless)
    !    T   [in, optional] : temperature in K (default 293.15 K = 20 °C)
    !    D   [result]       : translational diffusion coefficient (m²/s)
    !
    !  Water properties at 20 °C:
    !    eta  = 1.002e-3 Pa·s   (dynamic viscosity)
    !    rho  = 998.23 kg/m³    (mass density)
    !    vbar = 0.73e-3 m³/kg   (partial specific volume, typical protein)
    ! ===========================================================================
    pure function d_srf(s, fr, T) result(D)
        use stdlib_codata,      only: BOLTZMANN_CONSTANT
        use stdlib_codata_type, only: to_real
        real(c_double), intent(in)           :: s, fr
        real(c_double), intent(in), optional :: T
        real(c_double)                       :: D

        real(c_double), parameter :: pi   = acos(-1.0_c_double)
        real(c_double), parameter :: eta  = 1.002e-3_c_double
        real(c_double), parameter :: rho  = 998.23_c_double
        real(c_double), parameter :: vbar = 0.73e-3_c_double

        real(c_double) :: kB, T_

        kB = to_real(BOLTZMANN_CONSTANT, 1.0_c_double)

        if (present(T)) then
            T_ = T
        else
            T_ = 293.15_c_double
        end if

        D = (sqrt(2.0_c_double) / (18.0_c_double * pi)) &
            * kB * T_ * s**(-0.5_c_double)              &
            * (eta * fr)**(-1.5_c_double)               &
            * sqrt((1.0_c_double - vbar * rho) / vbar)

    end function d_srf

end module solvf_utilities