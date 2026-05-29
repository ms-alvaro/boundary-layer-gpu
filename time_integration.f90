!---------------------------------------------!
!     Module for temporal integration         !
!---------------------------------------------!
Module time_integration

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use ieee_arithmetic, Only : ieee_is_nan
  Use global
  Use equations
  Use projection
  Use boundary_conditions
  Use subgrid
  Use wallmodel

  ! prevent implicit typing
  Implicit None

Contains

  !-----------------------------------------------!
  !                Explicit Euler                 !
  !-----------------------------------------------!
  Subroutine compute_time_step_Euler

    ! equivalent to last RK step
    step_beginning = 1
    rk_step        = 3

    ! save current step
    !$acc kernels default(present)
    Uo = U
    Vo = V
    Wo = W
    !$acc end kernels
    ! Compute eddy viscosity (on CPU)
    !$acc update self(Uo,Vo,Wo)
    Call compute_eddy_viscosity(Uo,Vo,Wo,avg_nu_t,nu_t)
    ! Compute LES wall model (on CPU)
    Call compute_wall_model(Uo,Vo,Wo)
    !$acc update device(Uo,Vo,Wo,nu_t,avg_nu_t,alpha_x,alpha_y,alpha_z,V_bottom)
    ! compute rhs for U
    Call compute_rhs_u(Uo,Vo,Wo,rhs_uo)

    ! Advance U interior points
    !$acc kernels default(present)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rhs_uo
    !$acc end kernels

    ! compute rhs for V
    Call compute_rhs_v(Uo,Vo,Wo,rhs_vo)

    ! Advance V interior points
    !$acc kernels default(present)
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rhs_vo
    !$acc end kernels

    ! compute rhs for W
    Call compute_rhs_w(Uo,Vo,Wo,rhs_wo)

    ! Advance W interior points
    !$acc kernels default(present)
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rhs_wo
    !$acc end kernels

    ! Advance time
    t = t + dt

    ! boundary conditions
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

    ! projection step
    Call compute_projection_step

    ! boundary conditions
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

  End Subroutine compute_time_step_Euler

  !-----------------------------------------------!
  !                Explicit Euler                 !
  !-----------------------------------------------!
  Subroutine compute_time_step_Euler_test

    ! equivalent to last RK step
    step_beginning = 1
    rk_step        = 3

    ! Compute LES wall model
    Call compute_wall_model(U,V,W)

    ! Sanity check
    If ( Any( ieee_is_nan(U) ) ) Stop 'Error U NaNs 1!'
    If ( Any( ieee_is_nan(V) ) ) Stop 'Error V NaNs 1!'
    If ( Any( ieee_is_nan(W) ) ) Stop 'Error W NaNs 1!'

    Call apply_nonzero_Dirichlet_bc_y_bottom(V,1)

    ! Sanity check
    If ( Any( ieee_is_nan(U) ) ) Stop 'Error U NaNs 2!'
    If ( Any( ieee_is_nan(V) ) ) Stop 'Error V NaNs 2!'
    If ( Any( ieee_is_nan(W) ) ) Stop 'Error W NaNs 2!'

  End Subroutine compute_time_step_Euler_test

  !-----------------------------------------------!
  !          Explicit Runge-Kutta 2 steps         !
  !-----------------------------------------------!
  Subroutine compute_time_step_RK2

    Real(Int64) :: to

    step_beginning = 1

    ! save previous state
    to = t
    !$acc kernels default(present)
    Uo = U
    Vo = V
    Wo = W
    !$acc end kernels

    ! step 1 — eddy viscosity on CPU (1 transfer pair)
    rk_step = 1
    !$acc update self(U,V,W)
    Call compute_eddy_viscosity(U,V,W,avg_nu_t,nu_t)
    Call compute_wall_model(U,V,W)
    !$acc update device(nu_t,avg_nu_t,alpha_x,alpha_y,alpha_z,V_bottom)

    ! RHS, advance, BC, projection — all on GPU (only rhs_p transfers for FFT)
    Call compute_rhs_u(U,V,W,Fu1)
    Call compute_rhs_v(U,V,W,Fv1)
    Call compute_rhs_w(U,V,W,Fw1)

    !$acc kernels default(present)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rk2_coef(1,1)*Fu1
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rk2_coef(1,1)*Fv1
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rk2_coef(1,1)*Fw1
    !$acc end kernels
    t = to + rk2_t(rk_step)*dt

    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)
    Call compute_projection_step
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

    ! step 2 — eddy viscosity on CPU (1 transfer pair)
    rk_step = 2
    !$acc update self(U,V,W)
    Call compute_eddy_viscosity(U,V,W,avg_nu_t,nu_t)
    Call compute_wall_model(U,V,W)
    !$acc update device(nu_t,avg_nu_t,alpha_x,alpha_y,alpha_z,V_bottom)

    Call compute_rhs_u(U,V,W,Fu2)
    Call compute_rhs_v(U,V,W,Fv2)
    Call compute_rhs_w(U,V,W,Fw2)

    !$acc kernels default(present)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*( rk2_coef(2,1)*Fu1 + rk2_coef(2,2)*Fu2 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*( rk2_coef(2,1)*Fv1 + rk2_coef(2,2)*Fv2 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*( rk2_coef(2,1)*Fw1 + rk2_coef(2,2)*Fw2 )
    !$acc end kernels
    t = to + rk2_t(rk_step)*dt

    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)
    Call compute_projection_step
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

  End Subroutine compute_time_step_RK2

  !-----------------------------------------------!
  !          Explicit Runge-Kutta 3 steps         !
  !-----------------------------------------------!
  Subroutine compute_time_step_RK3

    Real(Int64) :: to

    step_beginning = 1

    ! save previous state
    to = t
    !$acc kernels default(present)
    Uo = U
    Vo = V
    Wo = W
    !$acc end kernels

    ! step 1
    rk_step = 1
    !$acc update self(U,V,W)
    Call compute_eddy_viscosity(U,V,W,avg_nu_t,nu_t)
    Call compute_wall_model(U,V,W)
    !$acc update device(U,V,W,nu_t,avg_nu_t,alpha_x,alpha_y,alpha_z,V_bottom)
    Call compute_rhs_u(U,V,W,Fu1)
    Call compute_rhs_v(U,V,W,Fv1)
    Call compute_rhs_w(U,V,W,Fw1)

    !$acc kernels default(present)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fu1
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rk_coef(1,1)*Fv1
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rk_coef(1,1)*Fw1
    !$acc end kernels
    t = to + rk_t(rk_step)*dt

    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)
    Call compute_projection_step
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

    ! step 2
    rk_step = 2
    !$acc update self(U,V,W)
    Call compute_eddy_viscosity(U,V,W,avg_nu_t,nu_t)
    Call compute_wall_model(U,V,W)
    !$acc update device(nu_t,avg_nu_t,alpha_x,alpha_y,alpha_z,V_bottom)
    Call compute_rhs_u(U,V,W,Fu2)
    Call compute_rhs_v(U,V,W,Fv2)
    Call compute_rhs_w(U,V,W,Fw2)

    !$acc kernels default(present)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*( rk_coef(2,1)*Fu1 + rk_coef(2,2)*Fu2 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*( rk_coef(2,1)*Fv1 + rk_coef(2,2)*Fv2 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*( rk_coef(2,1)*Fw1 + rk_coef(2,2)*Fw2 )
    !$acc end kernels
    t = to + rk_t(rk_step)*dt

    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)
    Call compute_projection_step
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

    ! step 3
    rk_step = 3
    !$acc update self(U,V,W)
    Call compute_eddy_viscosity(U,V,W,avg_nu_t,nu_t)
    Call compute_wall_model(U,V,W)
    !$acc update device(U,V,W,nu_t,avg_nu_t,alpha_x,alpha_y,alpha_z,V_bottom)
    Call compute_rhs_u(U,V,W,Fu3)
    Call compute_rhs_v(U,V,W,Fv3)
    Call compute_rhs_w(U,V,W,Fw3)

    !$acc kernels default(present)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fu1 + rk_coef(3,2)*Fu2 + rk_coef(3,3)*Fu3 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fv1 + rk_coef(3,2)*Fv2 + rk_coef(3,3)*Fv3 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + &
         dt*( rk_coef(3,1)*Fw1 + rk_coef(3,2)*Fw2 + rk_coef(3,3)*Fw3 )
    !$acc end kernels
    t = to + rk_t(rk_step)*dt

    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)
    Call compute_projection_step
    !$acc update self(U,V,W,Uo,Vo,Wo)
    Call apply_boundary_conditions
    !$acc update device(U,V,W)

  End Subroutine compute_time_step_RK3

  !-----------------------------------------------!
  !               test_projection                 !
  !-----------------------------------------------!
  Subroutine test_projection

    Real   (Int64) :: max_divergence
    Integer(Int32) :: i,j

    If (myid==0) Write(*,*) '----------------projection test-----------------'

    ! test velocities
    !$acc parallel loop collapse(2) default(present)
    Do i=1,nx
       Do j=1,nyg
          U(i,j,:) = i*j+j+i
       end Do
    end Do
    !$acc parallel loop collapse(2) default(present)
    Do i=1,nxg
       Do j=1,ny
          V(i,j,:) = i*j+j+i
       end Do
    end Do
    !$acc kernels default(present)
    W = 0d0
    !$acc end kernels

    ! mass correction
    !$acc update self(U,V,W)
    Call apply_global_mass_conservation(U,V,W)
    !$acc update device(U,V,W)

    ! projection step
    Call compute_projection_step

    ! mass correction
    !$acc update self(U,V,W)
    Call apply_global_mass_conservation(U,V,W)
    !$acc update device(U,V,W)

    ! divergence
    Call check_divergence(max_divergence)
    If ( myid==0 ) Write(*,*) 'max_divergence',max_divergence

    ! end
    Call sleep(3)
    Call Mpi_barrier(MPI_COMM_WORLD, ierr)
    Stop

  End Subroutine test_projection

  !-----------------------------------------------!
  !            compute dt based on CFL            !
  !-----------------------------------------------!
  ! NOTE: add rotating CFL and eddy viscosity CFL
  ! THIS IS WRONG, DEPENDS ON # PROCESSORS
  Subroutine compute_dt

    Integer(Int32) :: i, j, k
    Real   (Int64) :: lUmax, lVmax, lWmax, dt_local
    Real   (Int64) :: dt_conv_u, dt_conv_v, dt_conv_w, dt_conv
    Real   (Int64) :: dt_vis_u, dt_vis_v, dt_vis_w, dt_vis
    Real   (Int64) :: dt_max, CFLa

    CFLa = Abs( CFL )   ! absolute value of CFL ( for when CFL is dt )

    ! convective time step
    lUmax = 0d0
    lVmax = 0d0
    lWmax = 0d0
    !$acc update self(U,V,W)
    Do i=2,nxg-1
       Do j=2,nyg-1
          Do k=2,nzg-1
             lUmax = Max( lUmax,(xg(i+1)-xg(i))/Abs(U(i,j,k)) )
             lVmax = Max( lVmax,(yg(j+1)-yg(j))/Abs(V(i,j,k)) )
             lWmax = Max( lWmax,(zg(k+1)-zg(k))/Abs(W(i,j,k)) )
          End Do
       End Do
    End Do

    dt_conv_u = CFLa*lUmax
    dt_conv_v = CFLa*lVmax
    dt_conv_w = CFLa*lWmax

    dt_conv = Minval( (/dt_conv_u,dt_conv_v,dt_conv_w/) )

    ! viscous time step
    dt_vis_u = CFLa*dxmin**2d0/nu
    dt_vis_v = CFLa*dymin**2d0/nu
    dt_vis_w = CFLa*dzmin**2d0/nu

    dt_vis = Minval( (/dt_vis_u,dt_vis_v,dt_vis_w/) )

    ! time step
    dt_local = Min ( dt_conv,dt_vis )

    ! compute global minimum and communicate results to all processors
    Call MPI_Allreduce(dt_local,dt,1,MPI_real8,MPI_min,MPI_COMM_WORLD,ierr)

    ! NOTE garranz: remove alpha_mean_z since we are not using it!!!
    ! time step limiter
    !dt_max = alpha_mean_z
    !dt     = Min( dt, dt_max )

    dt_min_cfl = dt / CFLa ! save min dt for monitor output

    ! use multiple of the TS period instead
    If ( CFL<0 ) Then
       dt = -CFL
    end If

   End Subroutine compute_dt

End Module time_integration
