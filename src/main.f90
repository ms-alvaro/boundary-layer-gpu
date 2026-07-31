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
                           print_banner, param_summary, nsteps, random_init
  use grid_mod,      only: grid_generate, grid_to_device
  use field_mod,     only: fields_allocate, fields_to_device
  use ic_inflow_mod, only: generate_initial_condition,                       &
                           compute_blasius_solution_for_bc,                  &
                           inflow_tables_to_device
  use poisson_mod,   only: poisson_init, poisson_finalize
  use timestep_mod,  only: timestep_init, advance_one_step, t
  use bc_kernels,    only: bc_finalize
  use reductions,    only: reductions_finalize
  use io_mod,        only: io_init, read_restart,                            &
                           output_stats, output_monitor, output_snapshot

  implicit none
  integer :: istep

  ! ---- setup (mirrors the reference initialize() sequence) ----
  call print_banner()
  call read_input_parameters()
  call check_supported()

  call grid_generate()
  call grid_to_device()
  call fields_allocate()

  if (random_init == 1) then
     call generate_initial_condition(t)   ! Blasius IC into host mirrors, t=0
  else
     call read_restart()                  ! grids + fields + t from snapshot
     call grid_to_device()                ! grids may have been overwritten
  end if
  call fields_to_device()

  call compute_blasius_solution_for_bc()  ! inlet/top profiles + mode tables
  call inflow_tables_to_device()

  call poisson_init()                     ! operators, Thomas LU, cuFFT plans
  call timestep_init()                    ! RK tableaus, CFL spacing vectors
  call io_init()                          ! statistics arrays, wall clock

  call param_summary()

  ! ---- time loop (mirrors reference main.f90:84-113) ----
  do istep = 1, nsteps
     call advance_one_step(istep)         ! dt + RK substeps (RHS/BC/projection)
     call output_stats(istep)             ! self-gated: nstats or istep==1
     call output_monitor(istep)           ! self-gated: nmonitor
     call output_snapshot(istep)          ! self-gated: nsave (+ .restart link)
  end do

  ! ---- shutdown ----
  call poisson_finalize()
  call bc_finalize()
  call reductions_finalize()
  write(*,*) 'Done!'

end program boundary_layer_cuda
