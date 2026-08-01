!-----------------------------------------------------------------------------!
! Boundary-layer DNS solver — native CUDA Fortran implementation.             !
!                                                                             !
! Incompressible Navier-Stokes, 2nd-order staggered finite differences,      !
! RK2/RK3 fractional-step projection with a cuFFT-diagonalized Poisson       !
! solve. Single GPU (sm_80), no MPI. Port of the OpenACC reference solver    !
! (A. Lozano-Duran; GPU port 2026) — see docs/CUDA_PORT_DESIGN.md.           !
!                                                                             !
! Usage:  ./boundary_layer_cuda -i case.turbb                                 !
!-----------------------------------------------------------------------------!
program boundary_layer_cuda

  use param_mod,     only: read_input_parameters, check_supported,           &
                           print_banner, param_summary, nsteps, random_init, &
                           inflow_boundary_flag, top_boundary_flag
  use grid_mod,      only: grid_generate, grid_to_device
  use field_mod,     only: fields_allocate, fields_to_device
  use ic_inflow_mod, only: generate_initial_condition,                       &
                           compute_blasius_solution_for_bc,                  &
                           compute_turbulent_solution_for_bc,                &
                           compute_blowsuction_top,                          &
                           inflow_tables_to_device
  use lund_inflow_mod, only: lund_allocate, lund_init_means
  use les_mod,       only: les_allocate, les_dbg
  use wm_bridge,     only: wm_allocate
  use poisson_mod,   only: poisson_init, poisson_finalize
  use timestep_mod,  only: timestep_init, advance_one_step, t
  use bc_kernels,    only: bc_finalize, sync_z_ghosts
  use reductions,    only: reductions_finalize
  use mpi_mod,       only: mpi_initialize, mpi_finalize_run, is_root
  use hit_inflow_mod, only: init_hit_inflow
  use io_mod,        only: io_init, read_restart,                            &
                           output_stats, output_monitor, output_snapshot,    &
                           output_boxes

  implicit none
  integer :: istep

  ! ---- setup (mirrors the reference initialize() sequence) ----
  call mpi_initialize()                   ! rank/GPU binding (no-op alone)
  if (is_root()) call print_banner()
  call read_input_parameters()            ! every rank parses the input
  block
    character(16) :: envdbg
    call get_environment_variable('BL_LES_DBG', envdbg)
    if (len_trim(envdbg) > 0) les_dbg = .true.
  end block
  call check_supported()

  call grid_generate()
  call grid_to_device()
  call fields_allocate()
  call les_allocate()                     ! LES state (no-op for LES = 0)
  call wm_allocate()                      ! vendored wall-model stack
                                          ! (no-op for WM = 0, Tauw = 0)
  if (inflow_boundary_flag == 3 .or. inflow_boundary_flag == 5) then
     call lund_allocate()                 ! Lund EMA state + plane buffers
  end if                                  ! (before read_restart: the
                                          ! companion read fills Umean_resc_To)

  if (random_init == 1) then
     call generate_initial_condition(t)   ! Blasius IC into host mirrors, t=0
  else
     call read_restart()                  ! grids + fields + t from snapshot
     call grid_to_device()                ! grids may have been overwritten
  end if
  call fields_to_device()

  if (random_init /= 1) then
     call sync_z_ghosts()                 ! rebuild the z-ghost planes the
  end if                                  ! snapshot does not store (legacy
                                          ! 'force periodicity in z' fixup)

  if (inflow_boundary_flag == 3 .or. inflow_boundary_flag == 5) then
     call compute_turbulent_solution_for_bc() ! Lund inlet/top profiles
  else
     call compute_blasius_solution_for_bc()   ! inlet/top profiles + modes
  end if
  if (top_boundary_flag == 1 .or. top_boundary_flag == 2) then
     call compute_blowsuction_top()       ! V_bs(i) table for the lid BC
  end if
  call inflow_tables_to_device()

  if (inflow_boundary_flag == 3 .or. inflow_boundary_flag == 5) then
     ! io_mod::read_restart filled Umean_resc_To from the '.mean.rescaling'
     ! companion when it exists (inflow_flag = 3 restart); the copy into
     ! the live EMA state (or the IC-mean fallback) happens here, after
     ! the fields are resident on the device.
     call lund_init_means()
  end if

  if (inflow_boundary_flag == 6) then
     call init_hit_inflow()               ! HIT plane library (reads header,
  end if                                  ! loads + uploads the first buffer)

  call poisson_init()                     ! operators, Thomas LU, cuFFT plans
  call timestep_init()                    ! RK tableaus, CFL spacing vectors
  call io_init()                          ! statistics arrays, wall clock

  if (is_root()) call param_summary()

  ! ---- time loop (mirrors reference main.f90:84-113) ----
  do istep = 1, nsteps
     call advance_one_step(istep)         ! dt + RK substeps (RHS/BC/projection)
     call output_stats(istep)             ! self-gated: nstats or istep==1
     call output_monitor(istep)           ! self-gated: nmonitor
     call output_snapshot(istep)          ! self-gated: nsave (+ .restart link)
     call output_boxes(istep)             ! self-gated: causal-campaign boxes
  end do

  ! ---- shutdown ----
  call poisson_finalize()
  call bc_finalize()
  call reductions_finalize()
  if (is_root()) write(*,*) 'Done!'
  call mpi_finalize_run()

end program boundary_layer_cuda
