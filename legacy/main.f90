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
!        - OpenACC GPU offloading                      !
!                                                      !
! Required:                                            !
!       - FFTW 3.X                                     !
!       - LAPACK 3.X                                   !
!       - NVIDIA HPC SDK (for GPU)                     !
!                                                      !
! Parallel, Version 0.8.0-gpu                          !
!                                                      !
! Adrian Lozano Duran                                  !
! 2016                                                 !
! GPU port 2026                                        !
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

  ! Copy main arrays to GPU
  !$acc data copy(U,V,W,P,Uo,Vo,Wo) &
  !$acc      copyin(x,xm,xg,y,ym,yg,z,zm,zg,yg_m,yg_mm) &
  !$acc      copyin(weight_y_0,weight_y_1) &
  !$acc      copyin(xg_global,yg_global,zg_global) &
  !$acc      copyin(kxx,kzz,Dyy,imode_map,kmode_map) &
  !$acc      copyin(rk_coef,rk_t,rk2_coef,rk2_t) &
  !$acc      create(term,term_1,term_2) &
  !$acc      create(rhs_uo,rhs_vo,rhs_wo,rhs_p) &
  !$acc      copy(nu_t,avg_nu_t) &
  !$acc      create(Fu1,Fu2,Fv1,Fv2,Fw1,Fw2) &
  !$acc      create(Fu3,Fv3,Fw3) &
  !$acc      copy(alpha_x,alpha_y,alpha_z) &
  !$acc      copy(alpha_xo,alpha_yo,alpha_zo) &
  !$acc      copy(V_bottom) &
  !$acc      copyin(bc_1,bc_2) &
  !$acc      copyin(U_top,V_top,W_top,U_inlet,V_inlet,W_inlet) &
  !$acc      create(plane_gpu,rhs_hat_gpu) &
  !$acc      copyin(qu_inlet,qv_inlet,qw_inlet,zmode_inlet,tmode_inlet) &
  !$acc      create(Ut_inlet,Vt_inlet,Wt_inlet) &
  !$acc      copyin(thomas_dl_fact,thomas_d_pivot,thomas_du) &
  !$acc      copyin(dct_twiddle)

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

     ! Transfer GPU -> CPU only when needed (statistics, monitor, or snapshot)
     If ( istep==1 .Or. Mod(istep,nstats)==0 .Or. Mod(istep,nmonitor)==0 .Or. Mod(istep,nsave)==0 ) Then
        !$acc update self(U,V,W,P,nu_t)
     End If

     ! compute a few statistics
     Call compute_statistics
     !Call compute_statistics_z_modes

     ! output some key values
     Call output_monitor

     ! write snapshot if needed
     Call output_data

  End Do

  !$acc end data

  ! finalize stuff
  Call finalize

End program boundary_layer_FD
