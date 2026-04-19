! =============================================================================
!  test_sawp.f90
!
!  Standalone test program for the SAWP solver.
!  Calls sawp_solver() with a set of physically representative parameters,
!  then writes the full time-series to a formatted .out file.
!
!  Output file format  (sawp_results.out):
!    Line 1:   header with all input parameters
!    Line 2:   column labels
!    Line 3+:  one row per time step:
!              time   c(1)  c(2)  ...  c(N)
!
!  Build (assuming schuck_matrices.o and sawp_solver.o already compiled):
!    gfortran -O3 -march=native -std=f2008 \
!        schuck_matrices.o sawp_solver.o test_sawp.f90 \
!        -o test_sawp -llapack -lblas -lgfortran -lstdlib
! =============================================================================
program test_sawp
  use iso_c_binding,   only: c_double, c_int
  use sawp_solver_mod, only: sawp_solver
  implicit none
 
  ! ---------------------------------------------------------------------------
  !  INPUT PARAMETERS  — adjust these to change the test case
  ! ---------------------------------------------------------------------------
  !  Physical parameters
  real(c_double), parameter :: omega  = 5236.0_c_double   ! angular frequency [rad/s]
  real(c_double), parameter :: s      = 1.0e-13_c_double   ! shear parameter
  real(c_double), parameter :: D      = 1.0e-7_c_double  ! diffusion coefficient
  real(c_double), parameter :: m      = 6.0_c_double   ! meniscus radius
  real(c_double), parameter :: b      = 7.2_c_double   ! bottom radius
 
  !  Grid and time discretisation
  integer(c_int), parameter :: N       = 1000_c_int      ! number of grid points
  real(c_double), parameter :: t0      = 0.0_c_double  ! start time
  integer(c_int), parameter :: n_steps = 3000_c_int      ! number of time steps
  real(c_double), parameter :: dt      = 50.0_c_double  ! time step size

  ! --- Timing variables ------------------------------------------------------
  ! cpu_time: real scalar, seconds of CPU time
  real(c_double)  :: cpu_start, cpu_end, cpu_elapsed
  ! system_clock: integer ticks; count_rate converts to Hz
  integer(8)      :: wall_start, wall_end, wall_rate
  real(c_double)  :: wall_elapsed
 
  !  Output file
  character(len=*), parameter :: outfile = "sawp_results.out"
  ! ---------------------------------------------------------------------------
 
  ! Allocatable output arrays — sized at runtime from the parameters above
  real(c_double), allocatable :: c0(:)           ! initial coefficients, length N
  real(c_double), allocatable :: times_out(:)    ! time points, length n_steps+1
  real(c_double), allocatable :: coeffs_out(:,:) ! coefficients, shape (N, n_steps+1)
 
  integer(c_int) :: solver_info
  integer        :: step, j, iunit, alloc_stat
  character(len=256) :: iomsg
 
  ! ---------------------------------------------------------------------------
  !  1. Allocate arrays
  ! ---------------------------------------------------------------------------
  allocate(c0(N),                           stat=alloc_stat, errmsg=iomsg)
  call check_alloc(alloc_stat, iomsg, "c0")
  allocate(times_out(n_steps + 1),          stat=alloc_stat, errmsg=iomsg)
  call check_alloc(alloc_stat, iomsg, "times_out")
  allocate(coeffs_out(N, n_steps + 1),      stat=alloc_stat, errmsg=iomsg)
  call check_alloc(alloc_stat, iomsg, "coeffs_out")
 
  ! ---------------------------------------------------------------------------
  !  2. Set initial condition
  !     Gaussian bump centred in the domain: c0(k) = exp(-((r_k - r_mid)/w)^2)
  !     where r_k is the k-th Chebyshev-like grid point from grid_init and
  !     r_mid and w are the centre and width of the bump.
  !     Here we use a uniform distribution as the simplest meaningful test.
  ! ---------------------------------------------------------------------------
  c0 = 1.0_c_double / real(N, c_double)   ! uniform: total mass = 1
 
  ! ---------------------------------------------------------------------------
  !  3. Run the solver
  ! ---------------------------------------------------------------------------
  write(*, '(A)') "Running SAWP solver..."
  write(*, '(A,I0,A,I0,A, ES10.3)') "  Grid points N = ", N, &
                             ",  steps = ", n_steps, &
                             ",  dt = 0.1", dt
  
  ! Start timers immediately before the solver call
  call cpu_time(cpu_start)
  call system_clock(wall_start, wall_rate)

  solver_info = sawp_solver(c0, t0, n_steps, omega, s, D, m, b, N, dt, &
                             times_out, coeffs_out)
  

  ! Stop timers immediately after
  call cpu_time(cpu_end)
  call system_clock(wall_end)


  if (solver_info /= 0) then
    write(*, '(A,I0)') "ERROR: LAPACK returned info = ", solver_info
    stop 1
  end if
 
  ! Compute elapsed times
  cpu_elapsed  = cpu_end - cpu_start
  wall_elapsed = real(wall_end - wall_start, c_double) / real(wall_rate, c_double)

  write(*, '(A)') "  Solver completed successfully."

  ! ---------------------------------------------------------------------------
  !  4.1 Report timing
  ! ---------------------------------------------------------------------------

  write(*, '(A)') ""
  write(*, '(A)') "Timing:"
  write(*, '(A,F12.6,A)') "  Wall-clock time : ", wall_elapsed, " s"
  write(*, '(A,F12.6,A)') "  CPU time        : ", cpu_elapsed,  " s"
  write(*, '(A,F12.6,A)') "  Time per step   : ", &
      wall_elapsed / real(n_steps, c_double) * 1.0e3_c_double, " ms"
  if (cpu_elapsed > 1.0e-9_c_double) then
    write(*, '(A,F8.3)') "  CPU/Wall ratio  : ", cpu_elapsed / wall_elapsed
  end if
 
  ! ---------------------------------------------------------------------------
  !  4.2 Write output file
  ! ---------------------------------------------------------------------------
  iunit = 20
  open(unit=iunit, file=outfile, status='replace', action='write', &
       iostat=alloc_stat, iomsg=iomsg)
  if (alloc_stat /= 0) then
    write(*, '(A,A)') "ERROR opening output file: ", trim(iomsg)
    stop 1
  end if
 
  ! --- Header block ---
  write(iunit, '(A)')       "# SAWP Solver Output"
  write(iunit, '(A,F10.4)') "#   omega  = ", omega
  write(iunit, '(A,ES10.3)') "#   s      = ", s
  write(iunit, '(A,ES10.3)') "#   D      = ", D
  write(iunit, '(A,F10.4)') "#   m      = ", m
  write(iunit, '(A,F10.4)') "#   b      = ", b
  write(iunit, '(A,I10)')   "#   N      = ", N
  write(iunit, '(A,F10.4)') "#   t0     = ", t0
  write(iunit, '(A,F10.4)') "#   dt     = ", dt
  write(iunit, '(A,I10)')   "#   n_steps= ", n_steps
  write(iunit, '(A,ES10.3)') "#   t_end  = ", t0 + real(n_steps, c_double)*dt
  write(iunit, '(A,F12.6)') "#   wall_s  = ", wall_elapsed
  write(iunit, '(A,F12.6)') "#   cpu_s   = ", cpu_elapsed
  write(iunit, '(A)')       "#"
 
  ! --- Column label row ---
  ! Format: "time" then "c_001", "c_002", ... "c_NNN"
  write(iunit, '(A14)', advance='no') "# time        "
  do j = 1, N
    write(iunit, '(A2,I4.4,A6)', advance='no') "  c_", j, "      "
  end do
  write(iunit, '()')   ! newline
 
  ! --- Data rows: one row per time step ---
  ! coeffs_out is stored as (N, n_steps+1) col-major; step index is second dim.
  do step = 1, n_steps + 1
    write(iunit, '(ES14.6)', advance='no') times_out(step)
    do j = 1, N
      write(iunit, '(2X,ES14.6)', advance='no') coeffs_out(j, step)
    end do
    write(iunit, '()')   ! newline
  end do
 
  close(iunit)
 
  write(*, '(A,A,A,I0,A)') "  Results written to '", trim(outfile), &
                            "' (", n_steps + 1, " rows)."
 
  ! ---------------------------------------------------------------------------
  !  5. Quick sanity checks printed to stdout
  ! ---------------------------------------------------------------------------
  write(*, '(A)') ""
  write(*, '(A)') "Sanity checks:"
  write(*, '(A,F12.6)') "  Initial total mass  = ", &
      sum(coeffs_out(:, 1))
  write(*, '(A,F12.6)') "  Final   total mass  = ", &
      sum(coeffs_out(:, n_steps + 1))
  write(*, '(A,F12.6,A,F12.6)') "  Time range: ", &
      times_out(1), " -> ", times_out(n_steps + 1)
  write(*, '(A,F12.6)') "  Max coeff at t=0    = ", maxval(coeffs_out(:, 1))
  write(*, '(A,F12.6)') "  Max coeff at t_end  = ", &
      maxval(coeffs_out(:, n_steps + 1))
 
  ! ---------------------------------------------------------------------------
  !  6. Cleanup
  ! ---------------------------------------------------------------------------
  deallocate(c0, times_out, coeffs_out)
 
contains
 
  ! ---------------------------------------------------------------------------
  !  Internal helper: check allocation status and stop with a clear message
  ! ---------------------------------------------------------------------------
  subroutine check_alloc(stat, msg, varname)
    integer,          intent(in) :: stat
    character(len=*), intent(in) :: msg, varname
    if (stat /= 0) then
      write(*, '(A,A,A,A)') "ERROR: allocation of '", trim(varname), &
                             "' failed: ", trim(msg)
      stop 1
    end if
  end subroutine check_alloc
 
end program test_sawp