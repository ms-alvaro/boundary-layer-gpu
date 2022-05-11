!------------------------------------------------------!
! My boundary layer :)                                 !
!                                                      !
! Solve incompression Navier-Stokes eqs.               !
!                                                      !    
! Spatial discretization:                              !
!        - 2nd order finite differences                !
!        - Staggered mesh                              !
!                                                      !
! Temporal discretization:                             !
!        - Explicit Euler, RK2 and RK3                 !
!        - Fractional step method                      !
!                                                      !
! Boundary conditions:                                 !
!        - x inflow(t)/outflow                         !
!          or turbulent recycling inflow from Lung     !
!        - z periodic                                  !
!        - y non-slip/slip bottom, in/outflow top      !
!                                                      !
! SGS models                                           !
!        - Constant coef Smagorinsky (to be done)      !
!        - Dynamic Smagorinksy                         !
!                                                      !
! Parallelization:                                     !
!        - MPI, z-slices                               !
!                                                      !
! Required:                                            !
!       - FFTW 3.X                                     !
!       - LAPACK 3.X                                   !
!                                                      !
! Parallel, Version 0.7.1                              !
!                                                      !
! Adrian Lozano Duran                                  !
! 2016                                                 !
!------------------------------------------------------!
Program boundary_layer_FD

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use input_output
  Use initialization
  Use time_integration
  Use monitor
  Use statistics
  Use finalization
  
  ! prevent implicit typing
  Implicit None

  ! initialize everything and read input file and input flow field
  Call initialize

  ! small summary of input parameters
  Call summary
     
  ! temporal loop
  Do istep = 1, nsteps
     
     ! compute dt based on CFL
     Call compute_dt

     ! time step
     If     ( itime_step==1 ) Then
        Call compute_time_step_Euler
     Elseif ( itime_step==2 ) Then
        Call compute_time_step_RK2
     Elseif ( itime_step==3 ) Then
        Call compute_time_step_RK3
     End If

     ! compute a few statistics
     Call compute_statistics 
     !Call compute_statistics_z_modes

     ! output some key values
     Call output_monitor

     ! write snapshot if needed
     Call output_data

  End Do

  ! finalize stuff
  Call finalize

End program boundary_layer_FD
