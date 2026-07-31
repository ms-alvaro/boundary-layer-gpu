!------------------------------------!
!   Module for boundary conditions   !
!------------------------------------!
Module boundary_conditions

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use Lund_rescaled_bc
  Use mpi
  !Use ifport  ! removed for gfortran compatibility

  ! prevent implicit typing
  Implicit None

Contains

  !--------------------------------------------------!
  !     Apply boundary conditions to velocities      !
  !              in the 3 directions                 !
  !--------------------------------------------------!
  Subroutine apply_boundary_conditions

    !---------------------------------------------------------!
    ! interior region
    Call update_ghost_interior_planes(U,1)
    Call update_ghost_interior_planes(V,2)
    Call update_ghost_interior_planes(W,3)
   
    !---------------------------------------------------------!
    ! apply velocity inflow in x
    If     ( inflow_boundary_flag == 1 .or. inflow_boundary_flag == 2 ) Then
       ! temporal modes+Blasius inflow BC in x
       Call apply_inflow_bc_x(U,V,W) 
    Elseif ( inflow_boundary_flag == 3 ) Then
       ! Lund rescaling, computed only once at the beginning of step for Uo,Vo,Wo
       Call apply_inflow_bc_x_rescaling(Uo,Vo,Wo,U,V,W)
    Elseif ( inflow_boundary_flag == 4 ) Then
       ! temporal random+Blasius inflow BC in x
       Call apply_inflow_bc_x_Blasius_random(U,V,W) 
    Elseif ( inflow_boundary_flag == 5 ) Then
       ! Lund rescaling only for fluctuations
       Call apply_inflow_bc_x_Turbulent_rescaled_flu(Uo,Vo,Wo,U,V,W) 
       !Call apply_inflow_bc_x_Turbulent_random(Uo,Vo,Wo,U,V,W) 
    End If

    !---------------------------------------------------------!
    ! apply velocity outflow in x 
    Call apply_outflow_bc_x(U,V,W) 

    !---------------------------------------------------------!
    ! apply periodicity in z 
    Call apply_periodic_bc_z(U,1)
    Call apply_periodic_bc_z(V,2)
    Call apply_periodic_bc_z(W,3)

    !---------------------------------------------------------!
    ! apply velocity boundary conditions at the top
    If ( top_boundary_flag==1 .Or. top_boundary_flag==2 ) Then
       ! V boundary condition at the top    
       Call apply_Dirichlet_bc_y_top_BlowingSuction(V,top_boundary_flag)
       ! Impose zero vorticity and dw/dy = 0
       Call zero_wz_top(U,V,W)
    Elseif (top_boundary_flag==3) Then
       Call apply_top_bc_y_Falkner_Skan(U,V,W)
    Elseif (top_boundary_flag==4) Then
       ! zero-shear tangential (du/dy=dw/dy=0) + displacement V at the lid
       Call apply_top_bc_y_zeroshear(U,V,W)
    Else
       ! U, V and W boundary condition at the top           
       Call apply_top_bc_y(U,V,W) 
    End If

    
    !---------------------------------------------------------!
    ! U boundary condition at the wall
    If     ( iwall_model==0  ) Then
       ! apply Dirichlet in y for U
       Call apply_Dirichlet_bc_y_bottom(U,2)
    Elseif ( iwall_model==11 .Or. iwall_model==15 ) Then
       ! apply non-zero Neumann in y for U
       Call apply_Nozero_Neumann_bc_y(U,alpha_x,2)
    Else
       ! apply Robin in y for U
       Call apply_Robin_bc_y_bottom(U,alpha_x,2)
       !Call apply_pseudo_Robin_bc_y_bottom(U,alpha_x,1) 
    End If

    !---------------------------------------------------------!
    ! V boundary condition at the wall
    If     ( iwall_model==0 .Or. iwall_model==11 .Or. iwall_model==15 ) Then
       ! apply Dirichlet in y for V
       Call apply_Dirichlet_bc_y_bottom(V,1)
    Elseif ( iwall_model==13 ) Then
       ! apply Dirichlet in y for V
       Call apply_nonzero_Dirichlet_bc_y_bottom(V,1)
    Else
       ! apply Robin in y for V
       !Call apply_Robin_bc_y_bottom(V,alpha_y,1) !->needs modification of poisson solver
       Call apply_pseudo_Robin_bc_y_bottom(V,alpha_y,2) 
    End If
    
    !---------------------------------------------------------!
    ! W boundary condition at the wall
    If     ( iwall_model==0 .Or. iwall_model==11 .Or. iwall_model==15 ) Then
       ! apply Dirichlet in y for W
       Call apply_Dirichlet_bc_y_bottom(W,2)
    Else
       ! apply Robin in y for W
       Call apply_Robin_bc_y_bottom(W,alpha_z,2)
       !Call apply_pseudo_Robin_bc_y_bottom(W,alpha_z,3) 
    End If

    !---------------------------------------------------------!
    ! compute boundary conditions for pseudo-pressure
    !If ( top_boundary_flag == 2 ) Then
    !   Call compute_pseudo_pressure_bc_for_top_Neumann_bc
    !End If

    !---------------------------------------------------------!
    ! enforce global mass conservation
    Call apply_global_mass_conservation(U,V,W)

  End Subroutine apply_boundary_conditions

  !---------------------------------------------------------!
  !  GPU version: all BCs on GPU using module variables     !
  !  For nprocs=1, iwall_model=0, top_flag=0, inflow_flag=1 !
  !---------------------------------------------------------!
  Subroutine apply_boundary_conditions_gpu

    Integer(Int32) :: i, j, k, n, m
    Real   (Int64) :: Uc, Qx_local, Qy_local, Q_local, Q_total, Delta_U, length_y
    Integer(Int32) :: kk

    ! Ghost cell updates: nprocs=1 -> no-op

    ! Inflow BC — entirely on GPU
    ! inflow_flag=6: HIT planes -> Ut/Vt/Wt_inlet (temporal interpolation)
    If ( inflow_boundary_flag == 6 ) Then
       Call update_hit_inlet_gpu
    End If
    ! inflow_flag=1: temporal modes -> Ut/Vt/Wt_inlet (mode synthesis)
    If ( n_modes_inlet > 0 ) Then
       !$acc parallel loop collapse(2) default(present) &
       !$acc private(n, m)
       Do j = 1, nyg
          Do k = 1, nzg
             Ut_inlet(j,k) = 0d0
             Do n = 1, n_modes_inlet
                Do m = 1, m_modes_inlet
                   Ut_inlet(j,k) = Ut_inlet(j,k) + &
                      Real( qu_inlet(j,n,m)*cdexp(dcmplx(0d0,1d0)*zmode_inlet(n)*zg(k) &
                            - dcmplx(0d0,1d0)*tmode_inlet(m)*t) , 8)
                End Do
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(2) default(present) &
       !$acc private(n, m)
       Do j = 1, ny
          Do k = 1, nzg
             Vt_inlet(j,k) = 0d0
             Do n = 1, n_modes_inlet
                Do m = 1, m_modes_inlet
                   Vt_inlet(j,k) = Vt_inlet(j,k) + &
                      Real( qv_inlet(j,n,m)*cdexp(dcmplx(0d0,1d0)*zmode_inlet(n)*zg(k) &
                            - dcmplx(0d0,1d0)*tmode_inlet(m)*t) , 8)
                End Do
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(2) default(present) &
       !$acc private(n, m)
       Do j = 1, nyg
          Do k = 1, nz
             Wt_inlet(j,k) = 0d0
             Do n = 1, n_modes_inlet
                Do m = 1, m_modes_inlet
                   Wt_inlet(j,k) = Wt_inlet(j,k) + &
                      Real( qw_inlet(j,n,m)*cdexp(dcmplx(0d0,1d0)*zmode_inlet(n)*z(k) &
                            - dcmplx(0d0,1d0)*tmode_inlet(m)*t) , 8)
                End Do
             End Do
          End Do
       End Do
    End If
    ! Apply inflow: U(1,:,:) = Blasius + perturbation (on GPU)
    !$acc parallel loop default(present)
    Do j=1,nyg
       U(1,j,:) = U_inlet(j) + Ut_inlet(j,:)
       W(1,j,:) = W_inlet(j) + Wt_inlet(j,:)
    End Do
    !$acc parallel loop default(present)
    Do j=1,ny
       V(1,j,:) = V_inlet(j) + Vt_inlet(j,:)
    End Do

    ! Outflow BC — on GPU
    Uc = U_top(nx)
    !$acc kernels default(present)
    U( nx-1,:,:) = Uo( nx-1,:,:) - rk_t(rk_step)*dt*Uc*( Uo( nx-1,:,:) - Uo( nx-2,:,:) )/(x(nx-1)-x(nx-2))
    V(nxg-1,:,:) = Vo(nxg-1,:,:) - rk_t(rk_step)*dt*Uc*( Vo(nxg-1,:,:) - Vo(nxg-2,:,:) )/(xg(nxg-1)-xg(nxg-2))
    W(nxg-1,:,:) = Wo(nxg-1,:,:) - rk_t(rk_step)*dt*Uc*( Wo(nxg-1,:,:) - Wo(nxg-2,:,:) )/(xg(nxg-1)-xg(nxg-2))
    U( nx,:,:) = U( nx-1,:,:)
    V(nxg,:,:) = V(nxg-1,:,:)
    W(nxg,:,:) = W(nxg-1,:,:)
    !$acc end kernels

    ! Periodic BC in z (nprocs=1) — on GPU
    !$acc kernels default(present)
    ! U (id=1): centers
    U(:,:,1)     = U(:,:,nzg-2)
    U(:,:,nzg-1) = U(:,:,2)
    U(:,:,nzg)   = U(:,:,3)
    ! V (id=2): centers
    V(:,:,1)     = V(:,:,nzg-2)
    V(:,:,nzg-1) = V(:,:,2)
    V(:,:,nzg)   = V(:,:,3)
    ! W (id=3): faces
    W(:,:,1)  = W(:,:,nz-1)
    W(:,:,nz) = W(:,:,2)
    !$acc end kernels

    ! Top BC — on GPU.  top_flag=4: zero-shear tangential (du/dy=dw/dy=0) +
    ! displacement V_top (lets freestream turbulence pass the lid).
    ! Else (top_flag=0): Dirichlet tangential (U=U_inf, W=0) + displacement V_top
    ! (clamps FST u'/w' to ~0 at the lid).  V_top is identical -> top pressure BC unchanged.
    If ( top_boundary_flag == 4 ) Then
       !$acc parallel loop default(present)
       Do i=1,nx
          U(i,nyg,:) = U(i,nyg-1,:)
       End Do
       !$acc parallel loop default(present)
       Do i=1,nxg
          V(i, ny,:) = V_top(i)
          W(i,nyg,:) = W(i,nyg-1,:)
       End Do
    Else
       !$acc parallel loop default(present)
       Do i=1,nx
          U(i,nyg,:) = U_top(i)
       End Do
       !$acc parallel loop default(present)
       Do i=1,nxg
          V(i, ny,:) = V_top(i)
          W(i,nyg,:) = W_top(i)
       End Do
    End If

    ! Wall BC (iwall_model=0): Dirichlet — on GPU
    !$acc kernels default(present)
    ! U: center in y -> antisymmetric
    U(:,1,:) = -U(:,2,:)
    ! V: face in y -> zero
    V(:,1,:) = 0d0
    ! W: center in y -> antisymmetric
    W(:,1,:) = -W(:,2,:)
    !$acc end kernels

    ! Global mass conservation — GPU reduction
    Qx_local  = 0d0
    length_y  = 0d0
    kk        = 1
    If ( myid==(nprocs-1) ) kk = 2
    !$acc parallel loop default(present) reduction(+:Qx_local,length_y)
    Do j = 2, nyg-1
       Qx_local = Qx_local + Sum( (U(1,j,2:nzg-kk)-U(nx-1,j,2:nzg-kk))*(y(j)-y(j-1)) )
       length_y = length_y + y(j)-y(j-1)
    End Do
    Qy_local = 0d0
    !$acc parallel loop default(present) reduction(+:Qy_local)
    Do i=2,nxg-2
       Qy_local = Qy_local + Sum( (V(i,1,2:nzg-kk)-V(i,ny,2:nzg-kk))*(x(i)-x(i-1)) )
    End Do
    Q_local = Qx_local + Qy_local
    Call MPI_Allreduce(Q_local,Q_total,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Delta_U = Q_total/length_y/Real(nzg_global-3,8)
    !$acc kernels default(present)
    U(nx-1,:,:) = U(nx-1,:,:) + Delta_U
    U(nx,:,:) = U(nx-1,:,:)
    !$acc end kernels

  End Subroutine apply_boundary_conditions_gpu
 
  Subroutine zero_wz_top( U_, V_, W_ )
     Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, W_
     Real(Int64), Dimension(:,:,:), Intent(InOut) :: V_
     Integer(Int32) :: i
     Real(Int64) :: dy

     dy = yg(nyg) - yg(nyg-1)
     Do i=1,nxg-1 
        U_(i,nyg,:) = (dy/dx)*( V_(i+1,ny,:) - V_(i,ny,:) ) + U_(i,nyg-1,:)
     End Do
     W_(:,nyg,:) = W_(:,nyg-1,:)


     Do i=1,i_rescale
        U_(i,nyg,:) = U_top(i)
        V_(i,nyg,:) = V_top(i)
        W_(i,nyg,:) = W_top(i)
     End Do

  end Subroutine zero_wz_top

  !-------------------------------------------------!
  !                Periodicity in x                 !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at x faces          !
  !            id=2-> F defined at x centers        !
  ! Output: F                                       !
  ! (Not used here)                                 !
  !-------------------------------------------------!
  Subroutine apply_periodic_bc_x(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at x faces
       F( 1,:,:) = F(nx-1,:,:)
       F(nx,:,:) = F(   2,:,:)
    Else
       ! F defined at x centers
       F(    1,:,:) = F(nxg-2,:,:)
       F(nxg-1,:,:) = F(    2,:,:) ! see note*
       F(nxg  ,:,:) = F(    3,:,:)
    End If

    ! *Note: this is done in case the initial
    ! condition is not periodic. After the first
    ! step is no longer required
  End Subroutine apply_periodic_bc_x

  !-------------------------------------------------!
  !                Neumann in x                     !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at x faces          !
  !            id=2-> F defined at x centers        !
  ! Output: F                                       !
  ! (used for sgs model)                            !
  !-------------------------------------------------!
  Subroutine apply_Neumann_bc_x(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at x faces
       F( 1,:,:) = F(    2,:,:)
       F(nx,:,:) = F( nx-1,:,:)
    Else
       ! F defined at x centers
       F(  1,:,:) = F(    2,:,:)
       F(nxg,:,:) = F(nxg-1,:,:) ! see note*
    End If

    ! *Note: this is done in case the initial
    ! condition is not periodic. After the first
    ! step is no longer required
  End Subroutine apply_Neumann_bc_x

  !-------------------------------------------------!
  !     Non-zero Neumann boundary condition in y    !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         alpha (array with derivative)           !
  !         alpha(:,1,:) -> bottom wall             !
  !         alpha(:,2,:) -> top wall                !
  !                                                 !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_Nozero_Neumann_bc_y(F,alpha,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Real   (Int64), Intent(In)    :: alpha(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       Stop 'Error: apply_Nozero_Neumann_bc_y not defined for id==1'
    Else
       ! F defined at y center
       F(:,1,:) = F(:,2,:) - alpha(:,1,:)*(yg(2)-yg(1))
    End If

  End Subroutine apply_Nozero_Neumann_bc_y

  !-------------------------------------------------!
  !             Temporal inflow BC in x             !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  U,V,W, U_inlet, V_inlet, W_inlet        !
  ! Output: U,V,W                                   !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_inflow_bc_x(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_ ,V_ ,W_

    ! local variables
    Integer(Int32) :: j

    ! compute temporal component of the inflow
    Call compute_temporal_inflow

    ! variables at centers 
    Do j=1,nyg
       U_(1,j,:) = U_inlet(j) + Ut_inlet(j,:)
       W_(1,j,:) = W_inlet(j) + Wt_inlet(j,:)
    End Do
    ! variables at faces
    Do j=1,ny
       V_(1,j,:) = V_inlet(j) + Vt_inlet(j,:)
    End Do

  End Subroutine apply_inflow_bc_x

  !-------------------------------------------------!
  !              Blasius Inflow BC in x             !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  U,V,W, U_inlet, V_inlet, W_inlet        !
  ! Output: U,V,W                                   !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_inflow_bc_x_Blasius_random(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_ ,V_ ,W_

    ! local variables
    Integer(Int32) :: j, k
    Real   (Int64) :: Amp_loc, rr

    ! Apply Blasius + random perturbations once per time step
    ! (step_beginning=1 at the start of each step, set to 0 after first BC call)
    ! Perturbations are weighted by (1-U/U_inf) to concentrate inside BL
    If ( step_beginning == 1 ) Then

       Do j=1,nyg
          Amp_loc = Amplitude_perturbations * (1d0 - U_inlet(j))
          Do k=1,nzg
             Call random_number(rr)
             U_(1,j,k) = U_inlet(j) + Amp_loc*(rr-0.5d0)
          End Do
       End Do
       Do j=1,nyg
          Amp_loc = Amplitude_perturbations * (1d0 - U_inlet(j))
          Do k=1,nz
             Call random_number(rr)
             W_(1,j,k) = W_inlet(j) + Amp_loc*(rr-0.5d0)
          End Do
       End Do
       Do j=1,ny
          Amp_loc = Amplitude_perturbations * (1d0 - U_inlet(min(j,nyg)))
          Do k=1,nzg
             Call random_number(rr)
             V_(1,j,k) = V_inlet(j) + Amp_loc*(rr-0.5d0)
          End Do
       End Do

       step_beginning = 0

    End If

  End Subroutine apply_inflow_bc_x_Blasius_random

  !-------------------------------------------------!
  !              Blasius Inflow BC in x             !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  U,V,W, U_inlet, V_inlet, W_inlet        !
  ! Output: U,V,W                                   !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_inflow_bc_x_Turbulent_random(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_ , V_,  W_

    ! local variables
    Real   (Int64) :: Amp_per, rr
    Integer(Int32) :: j, k

    If ( step_beginning == 1 ) Then

       ! compute temporal component of the inflow
       !Call compute_temporal_inflow

       ! perturbation
       Amp_per = Amplitude_perturbations

       ! variables at centers
       Do j=1,nyg
          Do k=1,nzg
             Call random_number(rr)
             U_(1,j,k) = U_inlet(j) + Amp_per*(rr-0.5d0)*(U_inlet(ny)-U_inlet(j))
          End Do
       End Do
       Do j=1,nyg
          Do k=1,nz
             Call random_number(rr)
             W_(1,j,k) = W_inlet(j) + Amp_per*(rr-0.5d0)*(U_inlet(ny)-U_inlet(j))
          End Do
       End Do
       ! variables at faces
       Do j=1,ny
          Do k=1,nzg
             Call random_number(rr)
             V_(1,j,k) = V_inlet(j) + Amp_per*(rr-0.5d0)*(U_inlet(ny)-U_inlet(j))
          End Do
       End Do

       step_beginning = 0

    End If
       
  End Subroutine apply_inflow_bc_x_Turbulent_random

  !-------------------------------------------------!
  !              Blasius Inflow BC in x             !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  U,V,W, U_inlet, V_inlet, W_inlet        !
  ! Output: U,V,W                                   !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_inflow_bc_x_Turbulent_rescaled_flu(Uo_,Vo_,Wo_,U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(In   ) :: Uo_ ,Vo_, Wo_
    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_ , V_,  W_

    ! local variables
    Integer(Int32) :: j, k, ks

    If ( step_beginning == 1 ) Then

       ! compute temporal component of the inflow
       !Call compute_temporal_inflow
       Call compute_rescaled_inflow(Uo_,Vo_,Wo_,Ut_inlet,Vt_inlet,Wt_inlet,delta_inlet,1)
       
       ! variables at centers 
       Do j=1,nyg
          Do k=1,nzg

             ! swap indeces in the spanwise direction
             ks = mod( k + 16, nzg) + 1

             U_(1,j,k) = U_inlet(j) + Ut_inlet(j,k)
          End Do
       End Do
       Do j=1,nyg
          Do k=1,nz
             ! swap indeces in the spanwise direction
             ks = mod( k + 16, nz ) + 1

             W_(1,j,k) = W_inlet(j) + Wt_inlet(j,ks)
          End Do
       End Do
       ! variables at faces
       Do j=1,ny
          Do k=1,nzg

             ! swap indeces in the spanwise direction
             ks = mod( k + 16, nzg) + 1

             V_(1,j,k) = V_inlet(j) + Vt_inlet(j,k)
          End Do
       End Do

       step_beginning = 0

    End If
       
  End Subroutine apply_inflow_bc_x_Turbulent_rescaled_flu

  !-------------------------------------------------!
  !         Lund's rescaling inflow BC in x         !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  Uo,Vo,Wo,U,V,W,delta_inlet              !
  ! Output: U,V,W                                   !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_inflow_bc_x_rescaling(Uo_,Vo_,Wo_,U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(In ) :: Uo_ ,Vo_ ,Wo_
    Real(Int64), Dimension(:,:,:), Intent(Out) :: U_  ,V_  ,W_

    ! local variables
    Integer(Int32) :: j

    If ( step_beginning == 1 ) Then ! this is only called once per step

       ! compute rescaled inflow plane
       Call compute_rescaled_inflow(Uo_,Vo_,Wo_,Ut_inlet,Vt_inlet,Wt_inlet,delta_inlet,0)
       
       ! variables at centers 
       Do j=1,nyg
          U_(1,j,:) = Ut_inlet(j,:)
          W_(1,j,:) = Wt_inlet(j,:)
       End Do
       
       ! variables at faces
       Do j=1,ny
          V_(1,j,:) = Vt_inlet(j,:)
       End Do

       step_beginning = 0
       
    End If
    
  End Subroutine apply_inflow_bc_x_rescaling

  !-------------------------------------------------!
  !                 Outflow BC in x                 !
  !          No MPI communication required          !
  !                                                 !
  ! Velocities estimated with outflow bc:           !
  !  dU/dt + Uc dU/dx = 0                           !
  !  dV/dt + Uc dV/dx = 0                           !
  !  dW/dt + Uc dW/dx = 0                           !
  !  with Uc=U_inf                                  !
  !  First order accurate                           !
  !                                                 !
  ! Input:  U,V,W,U_top                             !
  ! Output: U,V,W                                   !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_outflow_bc_x(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, V_, W_

    ! local variables
    Integer(Int32) :: j
    Real   (Int64) :: Uc

    ! convective velocity is taken as U_inf at top-outlet
    Uc = U_top(nx)

    ! explicit first order approximation in space and time
    U_( nx-1,:,:) = Uo( nx-1,:,:) - rk_t(rk_step)*dt*Uc*( Uo( nx-1,:,:) - Uo( nx-2,:,:) )/(x ( nx-1)-x ( nx-2))
    V_(nxg-1,:,:) = Vo(nxg-1,:,:) - rk_t(rk_step)*dt*Uc*( Vo(nxg-1,:,:) - Vo(nxg-2,:,:) )/(xg(nxg-1)-xg(nxg-2))
    W_(nxg-1,:,:) = Wo(nxg-1,:,:) - rk_t(rk_step)*dt*Uc*( Wo(nxg-1,:,:) - Wo(nxg-2,:,:) )/(xg(nxg-1)-xg(nxg-2))

    ! last cell is dummy
    U_( nx,:,:) = U_( nx-1,:,:)
    V_(nxg,:,:) = V_(nxg-1,:,:)
    W_(nxg,:,:) = W_(nxg-1,:,:)

  End Subroutine apply_outflow_bc_x

  !-------------------------------------------------!
  !        Dirichlet boundary condition in y        !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  ! For now only wall                               !
  !-------------------------------------------------!
  Subroutine apply_Dirichlet_bc_y_top_BlowingSuction(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: i
    
    
    If ( id==1 ) Then ! Coleman 2018
       ! F defined at y faces
       Do i=1,nxg
          ! Note that the minus sign is accounted in ( x_bs - x_g )
          F(i,ny,:) = sqrt(2d0)*Vbs_max*(x_bs-xg(i))/sigma_bs*dexp(.5d0 - ((x_bs-xg(i))/sigma_bs)**2d0 ) &
                      + phi_bs
       End Do     
    Elseif ( id==2 ) Then ! Abe 2017
       ! F defined at y faces
       Do i=1,nxg
          F(i,ny,:) = sqrt(2d0)*Vbs_max*(x_bs-xg(i))/sigma_bs*dexp(phi_bs - ((x_bs-xg(i))/sigma_bs)**2d0 )
       End Do     
    Else
       ! F defined at y centers
       Stop 'Error in apply_Dirichlet_bc_y_top_BlowingSuction'
    End If

  End Subroutine apply_Dirichlet_bc_y_top_BlowingSuction
  
  !-------------------------------------------------!
  !   Change U_outlet for global mass conservation  !
  !                                                 !
  ! To be consistent with div(vel)=0 the inflow     ! 
  ! and outflow must balance at the boundaries      !
  !                                                 !
  ! For 2nd order FD staggered and z-periodic:      !
  !   <U_inlet - dy*U_outlet>_yz +                  !
  !   <V_inlet - dx*V_outlet>_xz = 0                !
  !            and                                  !
  !   U_outlet = U*_outlet + Delta_U                !
  !                  |          |                   !
  !             estimation  correction              !
  !                                                 !
  ! Input:  U,V,W                                   !
  ! Output: U(nx-1:nx,:,:)                          !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_global_mass_conservation(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, V_, W_

    ! local variables
    Integer(Int32) :: i, j, kk
    Real   (Int64) :: Qx_local, Qy_local, Q_local, Q_total, Delta_U, length_y

    ! mass flow in x
    Qx_local  = 0d0
    length_y  = 0d0
    kk        = 1
    If ( myid==(nprocs-1) ) kk = 2 ! last periodic point in z excluded
    Do j = 2, nyg-1
       Qx_local = Qx_local + Sum( (U_(1,j,2:nzg-kk)-U_(nx-1,j,2:nzg-kk))*(y(j)-y(j-1)) )
       length_y = length_y + y(j)-y(j-1)
    End Do
    ! mass flow in y    
    Qy_local = 0d0
    Do i=2,nxg-2
       Qy_local = Qy_local + Sum( (V_(i,1,2:nzg-kk)-V_(i,ny,2:nzg-kk))*(x(i)-x(i-1)) )
    End Do

    ! total mass flow
    Q_local = Qx_local + Qy_local
    Call MPI_Allreduce(Q_local,Q_total,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)    

    ! compute velocity correction
    Delta_U      = Q_total/length_y/Real(nzg_global-3,8) ! nzg_global-3, last plane excluded
    U_(nx-1,:,:) = U_(nx-1,:,:) + Delta_U

    ! last plane is dummy, just copied
    U_(nx,:,:) = U_(nx-1,:,:)
    
  End Subroutine apply_global_mass_conservation

  !-------------------------------------------------!
  ! Periodicity in z, MPI communication required    !
  !-------------------------------------------------!
  Subroutine apply_periodic_bc_z(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id
    Integer(Int64) :: n(3)

    ! Single-processor shortcut: direct copy, no MPI needed
    If ( nprocs==1 ) Then
       If (id == 3) Then
          F(:,:,1)  = F(:,:,nz-1)
          F(:,:,nz) = F(:,:,2)
       Elseif (id == 1) Then
          F(:,:,1)     = F(:,:,nzg-2)
          F(:,:,nzg-1) = F(:,:,2)
          F(:,:,nzg)   = F(:,:,3)
       Elseif (id == 2) Then
          F(:,:,1)     = F(:,:,nzg-2)
          F(:,:,nzg-1) = F(:,:,2)
          F(:,:,nzg)   = F(:,:,3)
       Elseif (id == 4) Then
          F(:,:,1)     = F(:,:,nzg-2)
          F(:,:,nzg-1) = F(:,:,2)
          F(:,:,nzg)   = F(:,:,3)
       End If
       Return
    End If

    ! save planes
    If ( myid==0 ) Then
      ! begin planes
      If (id == 3) Then
        ! F defined at z faces
        buffer_wi(:,:)   = F(:,:,2)
      Elseif (id == 1) Then 
        ! F defined at z centers
         buffer_ui(:,:,2) = F(:,:,2) 
         buffer_ui(:,:,3) = F(:,:,3)
      Elseif (id == 2) Then
        ! F defined at z centers
         buffer_vi(:,:,2) = F(:,:,2) 
         buffer_vi(:,:,3) = F(:,:,3)
      Elseif (id == 4) Then
        ! F defined at z centers
         buffer_ci(:,:,2) = F(:,:,2) 
         buffer_ci(:,:,3) = F(:,:,3)
      End If      

    End If

    If ( myid==nprocs-1 ) Then
       ! end planes
      If (id == 3) Then
        ! F defined at z faces
        buffer_we(:,:) = F(:,:, nz-1)
      Elseif (id == 1) Then
        ! F defined at z centers
        buffer_ue(:,:) = F(:,:,nzg-2)
      Elseif (id == 2) Then
        ! F defined at z centers
        buffer_ve(:,:) = F(:,:,nzg-2)
      Elseif (id == 4) Then
        ! F defined at z centers
        buffer_ce(:,:) = F(:,:,nzg-2)
      End If
    End If
    
    ! communicate planes
    If ( myid==0 ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_wi, nxg*nyg, Mpi_real8, nprocs-1, 3,  &
        buffer_we, nxg*nyg, Mpi_real8, nprocs-1, 3, MPI_COMM_WORLD,    &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ui, nx*nyg*2, Mpi_real8, nprocs-1, 1, &
        buffer_ue, nx*nyg, Mpi_real8, nprocs-1, 1, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_vi, nxg*ny*2, Mpi_real8, nprocs-1, 2, &
        buffer_ve, nxg*ny, Mpi_real8, nprocs-1, 2, MPI_COMM_WORLD,     &
        istat, ierr)
      !Elseif (id == 4) Then
        ! Send/receive U 
      !  Call Mpi_sendrecv(buffer_ci, nxg*nyg*2, Mpi_real8, nprocs-1, nprocs-1, &
      !  buffer_ce, nxg*nyg, Mpi_real8, nprocs-1, nprocs-1, MPI_COMM_WORLD,     &
      !  istat, ierr)
      End If
    End If

    If ( myid==nprocs-1 ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_we, nxg*nyg, Mpi_real8, 0, 3, &
        buffer_wi, nxg*nyg, Mpi_real8, 0, 3, MPI_COMM_WORLD,   &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ue, nx*nyg, Mpi_real8, 0, 1,  &
        buffer_ui, nx*nyg*2, Mpi_real8, 0, 1, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive V 
        Call Mpi_sendrecv(buffer_ve, nxg*ny, Mpi_real8, 0, 2,  &
        buffer_vi, nxg*ny*2, Mpi_real8, 0, 2, MPI_COMM_WORLD,  &
        istat, ierr)
      !Elseif (id == 4) Then
        ! Send/receive U 
      !  Call Mpi_sendrecv(buffer_ce, nxg*nyg, Mpi_real8, 0, nprocs-1,  &
      !  buffer_ci, nxg*nyg*2, Mpi_real8, 0, nprocs-1, MPI_COMM_WORLD,  &
      !  istat, ierr)
      End If
    End If

    If (id == 4) Then
      If (myid == 0) Then
        Call Mpi_sendrecv(buffer_ci, nxg*nyg*2, Mpi_real8, nprocs-1, 4, &
        buffer_ce, nxg*nyg, Mpi_real8, nprocs-1, 4, MPI_COMM_WORLD,     &
        istat, ierr)
      End If
      If (myid == nprocs-1) Then
        Call Mpi_sendrecv(buffer_ce, nxg*nyg, Mpi_real8, 0, 4,  &
        buffer_ci, nxg*nyg*2, Mpi_real8, 0, 4, MPI_COMM_WORLD,  &
        istat, ierr)
      End If
   End If

    ! apply conditions
    If ( myid==0 ) Then
      If (id == 3) Then 
        F(:,:,1) = buffer_we(:,:)       ! W_global(:,:,nz_global-1)       
      Elseif (id == 1) Then
        F(:,:,1) = buffer_ue(:,:)       ! U_global(:,:,nzg_global-2)
      Elseif (id == 2) Then
        F(:,:,1) = buffer_ve(:,:)       ! U_global(:,:,nzg_global-2)
      Elseif (id == 4) Then
        F(:,:,1) = buffer_ce(:,:)       ! U_global(:,:,nzg_global-2)
      End If
    End If
    If ( myid==nprocs-1 ) Then
      If (id == 3) Then
        F(:,:,nz   ) = buffer_wi(:,:)   ! W_global(:,:,2) 
      Elseif (id == 1) Then
        F(:,:,nzg-1) = buffer_ui(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_ui(:,:,3) ! U_global(:,:,3)
      Elseif (id == 2) Then
        F(:,:,nzg-1) = buffer_vi(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_vi(:,:,3) ! U_global(:,:,3)
      Elseif (id == 4) Then
        F(:,:,nzg-1) = buffer_ci(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_ci(:,:,3) ! U_global(:,:,3)
      End If    
    End If   
    Call Mpi_barrier(MPI_COMM_WORLD, ierr)

  End Subroutine apply_periodic_bc_z

  !-------------------------------------------------!
  !      Dirichlet boundary condition at the wall   !
  !        No MPI communication required            !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_Dirichlet_bc_y_bottom(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces
       F(:,1,:) = 0d0
    Else
       ! F defined at y centers
       F(:,1,:) = -F(:,2,:)
    End If

  End Subroutine apply_Dirichlet_bc_y_bottom

  !-------------------------------------------------!
  !      Dirichlet boundary condition at the wall   !
  !        No MPI communication required            !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_nonzero_Dirichlet_bc_y_bottom(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces
       F(:,1,:) = V_bottom
    Else
       Stop 'Error! apply_nonzero_Dirichlet_bc_y_bottom'
    End If

  End Subroutine apply_nonzero_Dirichlet_bc_y_bottom

  !-------------------------------------------------!
  !       Blasius boundary condition at the top     !
  !        No MPI communication required            !
  !                                                 !
  ! Input:  U,W, U_top, V_top, W_top                !
  ! Output: U,W  at the top                         !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_top_bc_y(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, V_, W_

    ! local variables
    Integer(Int32) :: i

    ! variables at centers
    Do i=1,nx
       U_(i,nyg,:) = U_top(i)
    End Do
    ! variables at faces
    Do i=1,nxg
       V_(i, ny,:) = V_top(i)
       W_(i,nyg,:) = W_top(i)
    End Do

  End Subroutine apply_top_bc_y

  !-------------------------------------------------!
  !    Falkner-Skan boundary condition at the top   !
  !        No MPI communication required            !
  !                                                 !
  ! Input:  U,V,W                                   !
  ! Output: U,V,W at the top                        !
  !                                                 !
  ! Uinf(x) = C*(x-x_origin)^m                      !
  ! m       = betaH/(2-betaH)                       !
  ! C       = to match Uinf                         !
  !                                                 !
  ! Note: assuming Uinf=1 and x_origin = -1         !
  !-------------------------------------------------!
  Subroutine apply_top_bc_y_Falkner_Skan(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, V_, W_

    ! local variables
    Integer(Int32) :: i
    Real   (Int64) :: C_falkner, x_falkner, Uinf
        
    ! origin of the boundary layer in current frame of reference
    x_falkner = -1d0 

    ! match with Uinf at x=0
    Uinf      = 1d0
    C_falkner = Uinf / ( ( x(1)-x_falkner )**( beta_hartree/(2d0-beta_hartree) ) )

    ! variables at centers
    Do i=1,nx
       U_(i,nyg,:) = C_falkner*(xg(i)-x_falkner)**( beta_hartree/(2d0-beta_hartree) )
       U_top(i)    = U_(i,nyg,1)
    End Do    
    ! free stress in V and W
    Do i=1,nxg
       !V_(i, ny,:) = V_(i, ny-1,:) 
       !V_(i, ny,:) = Vo(i, ny-1,:) ! approx
       V_(i, ny,:) = 0d0
       W_(i,nyg,:) = W_(i,nyg-1,:) 
    End Do

  End Subroutine apply_top_bc_y_Falkner_Skan

  !-------------------------------------------------!
  !  Zero-shear (Neumann) tangential + displacement !
  !  wall-normal velocity at the top  (top_flag=4)  !
  !                                                 !
  !  du/dy = dw/dy = 0  (zero shear stress; lets    !
  !  freestream-turbulence u'/w' pass at the lid),  !
  !  V kept = V_top = U_inf*0.0160*Re_x^(-1/7)       !
  !  (same displacement/entrainment V as top_flag=0,!
  !  so the top pressure BC is unchanged).          !
  !-------------------------------------------------!
  Subroutine apply_top_bc_y_zeroshear(U_,V_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, V_, W_

    ! local variables
    Integer(Int32) :: i

    ! tangential at centers: zero shear stress (Neumann), du/dy = 0
    Do i=1,nx
       U_(i,nyg,:) = U_(i,nyg-1,:)
    End Do
    ! faces: keep displacement entrainment V; zero shear stress in W (dw/dy=0)
    Do i=1,nxg
       V_(i, ny,:) = V_top(i)
       W_(i,nyg,:) = W_(i,nyg-1,:)
    End Do

  End Subroutine apply_top_bc_y_zeroshear

  !-------------------------------------------------!
  !        Dirichlet boundary condition in y        !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  ! For now only wall                               !
  !-------------------------------------------------!
  Subroutine apply_Dirichlet_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces
       F(:, 1,:) = 0d0
       F(:,ny,:) = 0d0
    Else
       ! F defined at y centers
       F(:,  1,:) = -F(:,    2,:)
       F(:,nyg,:) = -F(:,nyg-1,:)
    End If
 
  End Subroutine apply_Dirichlet_bc_y

  !-------------------------------------------------!
  !          Neumann boundary condition in y        !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_Neumann_bc_y_bottom(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces (first order)
       F(:,1,:) = F(:,2,:)
    Else
       ! F defined at y centers (second order)
       F(:,1,:) = F(:,2,:)
    End If

  End Subroutine apply_Neumann_bc_y_bottom

  !-------------------------------------------------!
  !          Neumann boundary condition in y        !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_Neumann_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces (first order)
       F(:, 1,:) = F(:,   2,:)
       F(:,ny,:) = F(:,ny-1,:)
    Else
       ! F defined at y centers (second order)
       F(:,  1,:) = F(:,    2,:)
       F(:,nyg,:) = F(:,nyg-1,:)
    End If

  End Subroutine apply_Neumann_bc_y

  !-----------------------------------------------------------!
  !            Robin boundary condition at the wall           !
  !                                                           !
  ! bottom wall at faces (first order):                       !
  ! U(1)  = alpha*(U(2)-U(1))/(y(2)-y(1)) + beta              !
  !                                                           !
  ! bottom wall at centers (second order):                    !
  ! (U(1)+U(2))/2 = alpha (U(2)-U(1))/(yg(2)-yg(1)) + beta    !
  !                                                           !
  !                                                           !
  ! Input:  F     (array to apply boundary conditions)        !
  ! Input:  alpha (array with slip length)                    !
  !         beta  (cosntant slip velocity)                    !
  !         id id=1-> F defined at y faces   (V)              !
  !            id=2-> F defined at y centers (U,W)            !
  ! Output: F                                                 ! 
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine apply_Robin_bc_y_bottom(F,alpha,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Real   (Int64), Intent(In)    :: alpha(:,:,:)
    Integer(Int32), Intent(In)    :: id

    ! local variables
    Real   (Int64) :: beta_y, beta_y_local
    Integer(Int32) :: kk

    If ( id==1 ) Then
       ! F defined at y faces (this is first order, could move to second in future) 
       F(:,1,:) = alpha(:,1,:)*F(:,2,:) / (y(2)-y(1)) / ( alpha(:,1,:)/(y(2)-y(1)) + 1d0 )
       ! beta_y for mass correction
       kk = 1
       If ( myid==(nprocs-1) ) kk = 2 ! last periodic point in z excluded for mass conservation
       beta_y_local = Sum( F(2:nxg-2,1,2:nzg-kk) )
       Call MPI_Allreduce(beta_y_local,beta_y,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       beta_y = beta_y / Real( (nxg-3)*(nzg_global-3) , 8 )
       F(:,1,:) = F(:,1,:) - beta_y 
    Elseif ( id==2 ) Then
       ! F defined at y centers (this is second order)
       F(:,1,:) = ( 2d0*alpha(:,1,:)/(yg(2) - yg(1)) - 1d0 )*F(:,2,:) / ( 2d0*alpha(:,1,:)/(yg(2) - yg(1)) + 1d0 )
    End If
    
  End Subroutine apply_Robin_bc_y_bottom

  !-----------------------------------------------------------!
  !          pseudo-Robin boundary condition at the wall      !
  !              U^{n+1} = alpha dU^n/dy + beta               !
  !                                                           !
  ! bottom wall at faces (V, first order):                    !
  ! U(1)  = alpha*(Uo(2)-Uo(1))/(y(2)-y(1)) + beta            !
  !                                                           !
  ! bottom wall at centers (U, W, second order):              !
  ! (U(1)+U(2))/2 = alpha (Uo(2)-Uo(1))/(yg(2)-yg(1))         !
  !                                                           !
  !                                                           !
  ! Input:  F     (array to apply boundary conditions)        !
  ! Input:  alpha (array with slip length)                    !
  !         beta  (constant slip velocity)                    !
  !         id id=1,2,3-> U,V,W                               !
  ! Output: F                                                 ! 
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine apply_pseudo_Robin_bc_y_bottom(F,alpha,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Real   (Int64), Intent(In)    :: alpha(:,:,:)
    Integer(Int32), Intent(In)    :: id

    ! local variables
    Real   (Int64) :: beta_y_local
    Integer(Int32) :: kk

    If ( id==1 ) Then
       ! U defined at y centers (this is second order)
       F(:,1,:) = 2d0*alpha(:,1,:)*( Uo(:,2,:) - Uo(:,1,:) )/( yg(2) - yg(1) ) - F(:,2,:)
    Elseif ( id==2 ) Then
       ! V defined at y faces (this is first order, could move to second in future) 
       F(:,1,:) =  alpha(:,1,:)*( Vo(:,2,:)-Vo(:,1,:) )/( y(2)-y(1) )
       ! beta_y for mass correction
       kk = 1
       If ( myid==(nprocs-1) ) kk = 2 ! last periodic point in z excluded for mass conservation
       beta_y_local = Sum( F(2:nxg-2,1,2:nzg-kk) )
       Call MPI_Allreduce(beta_y_local,beta_y,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       beta_y = beta_y / Real( (nxg-3)*(nzg_global-3) , 8 )
       F(:,1,:) = F(:,1,:) - beta_y 
    Elseif ( id==3 ) Then
       ! F defined at y centers (this is second order)
       F(:,1,:) = 2d0*alpha(:,1,:)*( Wo(:,2,:) - Wo(:,1,:) )/( yg(2) - yg(1) ) - F(:,2,:)
    End If
    
  End Subroutine apply_pseudo_Robin_bc_y_bottom

  !--------------------------------------------------!
  !          Update ghost interior planes            !
  !--------------------------------------------------!
  Subroutine update_ghost_interior_planes(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)    
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    If (id == 1) Then
      !----------------------update U-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_us = F(:,:,nzg-1) ! send buffer
      Call Mpi_sendrecv(buffer_us, nx*nyg, Mpi_real8, sendto, tagto,        &
           buffer_ur, nx*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_ur ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_us = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_us, nx*nyg, Mpi_real8, sendto, tagto,        &
           buffer_ur, nx*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_ur ! received buffer

    Elseif (id == 2) Then
      !----------------------update V-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_vs = F(:,:,nzg-1) ! send buffer
      Call Mpi_sendrecv(buffer_vs, nxg*ny, Mpi_real8, sendto, tagto,        &
           buffer_vr, nxg*ny, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_vr ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_vs = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_vs, nxg*ny, Mpi_real8, sendto, tagto,        &
           buffer_vr, nxg*ny, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_vr ! received buffer
      
    Elseif (id == 3) Then
      !----------------------update W-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_ws = F(:,:,nz-1)    ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_wr ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_ws = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nz) = buffer_wr ! received buffer     

    Elseif (id == 4) Then
      !----------------------update W-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_ws = F(:,:,nzg-1)    ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_wr ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_ws = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_wr ! received buffer     
    End if  
    
  End Subroutine update_ghost_interior_planes

  !--------------------------------------------------------------!
  !                                                              !
  !          Compute top turbulent boundary conditions           !
  !                                                              !
  !  Consistency with expected turbulent BL growth:              !
  !                                                              !
  !  Vinf         = Uinf*d(delta*)/dx                            !
  !                                                              !
  !  delta*       = displacement thickness                       !
  !  delta*/x     = 0.020 *Rex^(-1/7) (version 1) used           !
  !  delta*/x     = 0.048 *Rex^(-1/5) (version 2)                !
  !  d(delta*)/dx = 0.0160*Rex^(-1/5) (version 1) used           ! 
  !  d(delta*)/dx = 0.0384*Rex^(-1/5) (version 2)                ! 
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_turbulent_solution_for_bc

    Real   (Int64) :: U_inf, Rex, d_delta
    Integer(Int32) :: i, ny_profile

    ! set U_inf
    U_inf = 1d0

    ! allocate boundary velocities (used later)
    Allocate(U_inlet (nyg_global),V_inlet ( ny_global),W_inlet (nyg_global))
    Allocate(U_outlet(nyg_global),V_outlet( ny_global),W_outlet(nyg_global))
    Allocate(U_top   (nx_global ),V_top   (nxg_global),W_top   (nxg_global))

    ! allocate auxiliary arrays (used later)
    Allocate( Ut_inlet(nyg_global,nzg) )
    Allocate( Vt_inlet( ny_global,nzg) )
    Allocate( Wt_inlet(nyg_global,nz ) )

    ! set top velocities
    U_top = U_inf
    Do i = 1, nxg_global
       Rex      = U_inf*xg_global(i)/nu
       d_delta  = 0.0160d0*Rex**(-1d0/7d0)
       V_top(i) = U_inf*d_delta
    End Do
    W_top = 0d0

    ! Impose own turbulent profile at the inlet
    If ( inflow_boundary_flag == 5) Then

       ! Need to match number of points and Rex
       If ( myid==0 ) Then 
          Write(*,*) 'Reading Turbulent profile for inlet from file'       
          ! read blasius solution
          Open(5,file=file_inflow,form='unformatted',action='read',access='stream')
          Read(5) ny_profile
          If (ny_profile/=ny_global ) Then
             Write(*,*) 'file_blasius_own',file_inflow
             Write(*,*)  'Profile doesnt match dimensions', ny_profile, ny_global
             Stop
          End If
          Read(5)  U_inlet
          Read(5)  V_inlet
          Close(5)
          W_inlet = 0d0
       End If

       ! broadcast to all processors
       Call Mpi_bcast ( U_inlet,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
       Call Mpi_bcast ( V_inlet, ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
       Call Mpi_bcast ( W_inlet,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
       
    End If       

  End Subroutine compute_turbulent_solution_for_bc

  !--------------------------------------------------!
  !                                                  !
  ! Compute Blasius solution for boundary conditions !
  !         and prepare temporal component           !
  !                                                  !
  ! Normalization: U_inf, x0 (distance leading edge) !
  !                                                  !
  ! Equations:                                       !
  !    U/U_inf = df                                  !
  !    V/U_inf = 1/sqrt(2*Rex)*(eta*df-f)            !
  !    W/U_inf = 0                                   !
  !    eta     = y*sqrt(U_inf/(2*nu*x))              !
  !    Rex     = U_inf*x/nu                          !
  !    Rex0    = U_inf*x0/nu                         !
  !                                                  !
  ! Notes:                                           !
  !    V_inf = 0.8604*U_inf*sqrt(nu/x0/U_inf)        !
  !    nu for Rex_inlet = x0*U_inf/nu (desired)      !
  !                                                  !
  ! Flags:                                           !
  !    inflow_boundary_flag = 1 -> generate Blasius  !
  !    inflow_boundary_flag = 2 -> read Blasius      !
  !                                                  !
  !--------------------------------------------------!
  ! NOTE: COULD BE IMPROVE IMPOSING EXACT BOUNDARY FOR W at top
  Subroutine compute_blasius_solution_for_bc
   
    ! solution from file (different size than the mesh)
    Real(Int64), Allocatable, Dimension(:) :: eta_source, f_source, df_source
    
    Real   (Int64) :: eta_local, Rex0, Rex_ref, w0, w1, U_inf
    Integer(Int32) :: j, jj, j0, j1, i_ref, n_source, ntotal, ny_blasius

    !----------------------------------------------------------------------!
    ! PART 1: compute Blasius

    ! set U_inf
    U_inf = 1d0

    ! allocate boundary velocities
    Allocate(U_inlet (nyg_global),V_inlet (ny_global),W_inlet (nyg_global))
    Allocate(U_outlet(nyg_global),V_outlet(ny_global),W_outlet(nyg_global))
    Allocate(U_top(nx_global),V_top(nxg_global),W_top(nxg_global))

    ! only processor 0 computes the solution
    If ( myid==0 ) Then
       
       ! read self-similar blasius solution       
       Open(4,file=file_inflow,form='formatted',action='read')
       Read(4,*) n_source
       Allocate(eta_source(n_source))
       Allocate(  f_source(n_source))
       Allocate( df_source(n_source))
       Read(4,*) eta_source
       Read(4,*)   f_source
       Read(4,*)  df_source
       Close(4)       
          
       ! inflow Rex
       Rex0 = 1d0*x(1)/nu
          
       If ( inflow_boundary_flag==1 .or. inflow_boundary_flag==3 .or. &
            inflow_boundary_flag==4 .or. inflow_boundary_flag==6 ) Then

          ! Generate own Blasius
          Write(*,*) 'Generating own Blasius for inlet'

          ! compute solution at inlet
          ! U
          U_inlet = 1d0 ! old version: U_inlet = 0d0
          i_ref   = 1
          Rex_ref = 1d0*x(i_ref)/nu
          Do j=1,nyg
             eta_local = yg(j)*(U_inf/(2d0*nu*x(i_ref)))**0.5d0
             j0 = 0
             Do jj=2,n_source
                If ( eta_source(jj)>eta_local ) Then
                   j0 = jj - 1
                   j1 = jj 
                   w1 = ( eta_local - eta_source(j0) )/( eta_source(j1) - eta_source(j0) )
                   w0 = 1d0 - w1
                   Exit
                End If
             End Do
             If ( j0>0 ) Then
                U_inlet(j) = w0*df_source(j0) + w1*df_source(j1)
             Else
                U_inlet(j) = 1d0
             End If
          End Do
          ! V
          V_inlet = 0d0
          i_ref   = 1
          Rex_ref = 1d0*xg(i_ref)/nu
          Do j=1,ny
             eta_local = y(j)*(U_inf/(2d0*nu*xg(i_ref)))**0.5d0
             j0 = 0
             Do jj=2,n_source
                If ( eta_source(jj)>eta_local ) Then
                   j0 = jj - 1
                   j1 = jj
                   w1 = ( eta_local - eta_source(j0) )/( eta_source(j1) - eta_source(j0) )
                   w0 = 1d0 - w1
                   Exit
                End If
             End Do
             If ( j0>0 ) Then 
                V_inlet(j) = w0*1d0/(2d0*Rex_ref)**0.5d0*(eta_source(j0)*df_source(j0)-f_source(j0)) + &
                     w1*1d0/(2d0*Rex_ref)**0.5d0*(eta_source(j1)*df_source(j1)-f_source(j1))
             Else
                V_inlet(j) = Maxval(V_inlet)
             End If
             
          End Do
          ! W
          W_inlet = 0d0

       Elseif ( inflow_boundary_flag==2 ) Then

          ! Read Blasius from file. 
          ! Need to match number of points and Rex
          If (inflow_boundary_flag==2) Write(*,*) 'Reading Blasius for inlet from file'

          ! read blasius solution       
          Open(5,file=file_inflow,form='unformatted',action='read',access='stream')
          Read(5) ny_blasius
          If (ny_blasius/=ny_global ) Then 
             Write(*,*) 'file_blasius_own',file_inflow
             Write(*,*)  'Blasius doesnt match dimensions', ny_blasius, ny_global
             Stop
          End If
          Read(5)  U_inlet
          Read(5)  V_inlet
          Close(5)       
          W_inlet = 0d0          

       Else
          Stop 'Error! inflow_boundary_flag unknown'
       End If          

       ! compute solution at outlet
       ! U
       U_outlet = 0d0
       i_ref    = nx-1 ! not nx
       Rex_ref  = 1d0*x(i_ref)/nu
       Do j=1,nyg
          eta_local = yg(j)*(U_inf/(2d0*nu*x(i_ref)))**0.5d0
          j0 = 0
          Do jj=2,n_source
             If ( eta_source(jj)>eta_local ) Then
                j0 = jj - 1
                j1 = jj 
                w1 = ( eta_local - eta_source(j0) )/( eta_source(j1) - eta_source(j0) )
                w0 = 1d0 - w1
                Exit
             End If
          End Do
          If ( j0>0 ) Then              
             U_outlet(j) = w0*df_source(j0) + w1*df_source(j1)
          Else
             U_outlet(j) = 1d0
          End If
       End Do
       ! V
       V_outlet = 0d0
       i_ref    = nxg-1 ! not nxg
       Rex_ref  = 1d0*xg(i_ref)/nu
       Do j=1,ny
          eta_local = y(j)*(U_inf/(2d0*nu*xg(i_ref)))**0.5d0
          j0 = 0
          Do jj=2,n_source
             If ( eta_source(jj)>eta_local ) Then
                j0 = jj - 1
                j1 = jj
                w1 = ( eta_local - eta_source(j0) )/( eta_source(j1) - eta_source(j0) )
                w0 = 1d0 - w1
                Exit
             End If
          End Do
          If ( j0>0 ) Then
             V_outlet(j) = w0*1d0/(2d0*Rex_ref)**0.5d0*(eta_source(j0)*df_source(j0)-f_source(j0)) + &
                  w1*1d0/(2d0*Rex_ref)**0.5d0*(eta_source(j1)*df_source(j1)-f_source(j1))
          Else
             V_outlet(j) = Maxval(V_outlet)
          End If
       End Do
       ! W
       W_outlet = 0d0
       
       ! compute solution at top
       ! U
       U_top = 1d0
       ! V
       Do i_ref=1,nxg
          Rex_ref      = 1d0*xg(i_ref)/nu
          j            = ny
          eta_local    = y(j)*(U_inf/(2d0*nu*xg(i_ref)))**0.5d0
          j0           = n_source
          V_top(i_ref) = 1d0/(2d0*Rex_ref)**0.5d0*(eta_source(j0)*df_source(j0)-f_source(j0))
       End Do
       ! W 
       W_top = 0d0
       
    End If

    ! broadcast to all processors
    Call Mpi_bcast ( U_inlet,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( V_inlet, ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( W_inlet,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( U_outlet,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( V_outlet, ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( W_outlet,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( U_top, nx_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( V_top,nxg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( W_top,nxg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    !----------------------------------------------------------------------!
    ! PART 2: prepare temporal component for boundary conditions

    ! allocate inlet temporal component
    Allocate( Ut_inlet(nyg_global,nzg) )
    Allocate( Vt_inlet( ny_global,nzg) )
    Allocate( Wt_inlet(nyg_global,nz ) )

    ! Skip temporal modes if no file is provided (pure Blasius, no perturbations)
    If ( Len_Trim(file_temporal_inlet)==0 .Or. &
         file_temporal_inlet(1:4)=='**TO' ) Then
       If ( myid==0 ) Write(*,*) 'No temporal inflow file -> pure Blasius (no perturbations)'
       Ut_inlet = 0d0
       Vt_inlet = 0d0
       Wt_inlet = 0d0
       n_modes_inlet = 0
       m_modes_inlet = 0
       beta_inlet    = 0d0
       omega_inlet   = 1d0
       dt_period     = 1d0
       Return
    End If

    ! read sizes and wavenumbers
    If ( myid==0 ) Then
       Open(4,file=file_temporal_inlet,form='unformatted',action='read',access='stream')
       Read(4)      ny_inlet
       Read(4) n_modes_inlet
       Read(4) m_modes_inlet
       Read(4)    beta_inlet
       Read(4)   omega_inlet
    End If

    ! broadcast to all processors
    Call Mpi_bcast (      ny_inlet,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( n_modes_inlet,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( m_modes_inlet,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    beta_inlet,1,MPI_real8,  0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   omega_inlet,1,MPI_real8,  0,MPI_COMM_WORLD,ierr )

    ! time multiple of the period
    dt_period = 2d0*pi/(abs(CFL)*omega_inlet)/Real(1000,8)

    ! allocate y-mesh, wavenumbers and coefficients
    Allocate( ymesh_inlet(     ny_inlet) ) ! wall-normal points
    Allocate( zmode_inlet(n_modes_inlet) ) ! number of spanwise modes
    Allocate( tmode_inlet(m_modes_inlet) ) ! number of temporal modes

    Allocate( qu_inlet(nyg_global,n_modes_inlet,m_modes_inlet) ) ! u modes y-interpolated
    Allocate( qv_inlet( ny_global,n_modes_inlet,m_modes_inlet) ) ! v modes y-interpolated
    Allocate( qw_inlet(nyg_global,n_modes_inlet,m_modes_inlet) ) ! w modes y-interpolated

    ! read y-mesh and modes
    If ( myid==0 ) Then

       ! temporal arrays
       ! real part
       Allocate( qu_inlet_r(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! u modes from source
       Allocate( qv_inlet_r(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! v modes from source
       Allocate( qw_inlet_r(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! w modes from source
       ! imaginary part
       Allocate( qu_inlet_i(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! u modes from source
       Allocate( qv_inlet_i(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! v modes from source
       Allocate( qw_inlet_i(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! w modes from source
       ! complex
       Allocate( qu_inlet_o(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! u modes from source
       Allocate( qv_inlet_o(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! v modes from source
       Allocate( qw_inlet_o(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! w modes from source

       ! read
       Read(4) ymesh_inlet
       Read(4) zmode_inlet
       Read(4) tmode_inlet
       Read(4) qu_inlet_r
       Read(4) qu_inlet_i
       Read(4) qv_inlet_r
       Read(4) qv_inlet_i
       Read(4) qw_inlet_r
       Read(4) qw_inlet_i
       Close(4)

       ! build wavenumber
       zmode_inlet =  beta_inlet*zmode_inlet
       tmode_inlet = omega_inlet*tmode_inlet

       ! complex modes 
       qu_inlet_o = dcmplx( qu_inlet_r, qu_inlet_i )
       qv_inlet_o = dcmplx( qv_inlet_r, qv_inlet_i )
       qw_inlet_o = dcmplx( qw_inlet_r, qw_inlet_i )

       ! interpolate qu and qw to yg_global
       qu_inlet = (0d0,0d0)
       qw_inlet = (0d0,0d0)
       Do j = 1, nyg_global
          j1 = 0          
          Do jj = 2, ny_inlet
             If ( ymesh_inlet(jj)>=yg_global(j) ) Then
                j1 = jj
                Exit
             End If
             If ( jj==ny_inlet ) j1 = ny_inlet
          End Do
          If ( j1==ny_inlet ) Then
             qu_inlet(j,:,:) = 0d0
             qw_inlet(j,:,:) = 0d0
          Else
             j0 = j1 - 1
             qu_inlet(j,:,:) = qu_inlet_o(j0,:,:) + & 
                               (yg_global(j)-ymesh_inlet(j0))*(qu_inlet_o(j1,:,:)-qu_inlet_o(j0,:,:))/(ymesh_inlet(j1)-ymesh_inlet(j0))
             qw_inlet(j,:,:) = qw_inlet_o(j0,:,:) + & 
                               (yg_global(j)-ymesh_inlet(j0))*(qw_inlet_o(j1,:,:)-qw_inlet_o(j0,:,:))/(ymesh_inlet(j1)-ymesh_inlet(j0))
          End If
       End Do

       ! interpolate qv to present y_global mesh
       qv_inlet = (0d0,0d0)
       Do j = 1, ny_global
          j1 = 0          
          Do jj = 2, ny_inlet
             If ( ymesh_inlet(jj)>=y_global(j) ) Then
                j1 = jj
                Exit
             End If
             If ( jj==ny_inlet ) j1 = ny_inlet
          End Do
          If ( j1==ny_inlet ) Then
             qv_inlet(j,:,:) = 0d0
          Else
             j0 = j1 - 1
             qv_inlet(j,:,:) = qv_inlet_o(j0,:,:) + & 
                               (y_global(j)-ymesh_inlet(j0))*(qv_inlet_o(j1,:,:)-qv_inlet_o(j0,:,:))/(ymesh_inlet(j1)-ymesh_inlet(j0))
          End If
       End Do
       
       ! deallocate
       Deallocate(qu_inlet_r,qv_inlet_r,qw_inlet_r)
       Deallocate(qu_inlet_i,qv_inlet_i,qw_inlet_i)
       Deallocate(qu_inlet_o,qv_inlet_o,qw_inlet_o)

    End If

    ! broadcast to all processors
    Call Mpi_bcast ( ymesh_inlet,     ny_inlet,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( zmode_inlet,n_modes_inlet,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( tmode_inlet,m_modes_inlet,MPI_real8,0,MPI_COMM_WORLD,ierr )

    ntotal = 2*nyg_global*n_modes_inlet*m_modes_inlet
    Call Mpi_bcast ( qu_inlet,ntotal,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( qw_inlet,ntotal,MPI_real8,0,MPI_COMM_WORLD,ierr )

    ntotal = 2*ny_global*n_modes_inlet*m_modes_inlet
    Call Mpi_bcast ( qv_inlet,ntotal,MPI_real8,0,MPI_COMM_WORLD,ierr )
    
  End Subroutine compute_blasius_solution_for_bc

  !-------------------------------------------------------------------!
  !  HIT plane inflow (inflow_flag=6)                                 !
  !                                                                   !
  !  Time-resolved y-z planes from a precursor HIT simulation,        !
  !  preprocessed (rescaled+blended) by preprocess_planes.py (v2).    !
  !  Planes are z-PERIODIC with period Lz_box_hit (must equal Lz).    !
  !                                                                   !
  !  - init_hit_inflow:      read header, allocate, load 1st buffer   !
  !  - load_hit_buffer:      read N_buffer planes, interp to          !
  !                          staggered grids (CPU)                    !
  !  - update_hit_inlet_gpu: temporal interpolation -> Ut/Vt/Wt_inlet !
  !                          (GPU kernel), reload buffer if exhausted !
  !-------------------------------------------------------------------!
  Subroutine init_hit_inflow

    Integer(Int32) :: magic
    Real   (Int64) :: L11r_f, alphaL_f, alphat_f, Uinf_f

    If ( nprocs > 1 ) Stop 'inflow_flag=6 GPU: only nprocs=1 supported'

    If ( Len_Trim(file_hit_planes)==0 ) &
         Stop 'inflow_flag=6 requires hit_file in input'

    If (myid==0) Write(*,*) 'Initializing HIT plane inflow from: ', &
                            Trim(file_hit_planes)

    Open(47, file=Trim(file_hit_planes), form='unformatted', &
         action='read', access='stream', status='old')
    Read(47) magic
    If ( magic /= 20260612 ) Then
       Write(*,*) 'ERROR: bad magic in hit_file: ', magic
       Stop 'hit_file format mismatch (need v2, magic 20260612)'
    End If
    Read(47) n_planes_hit
    Read(47) ny_hit_f
    Read(47) nz_hit_f
    Read(47) dt_plane_hit
    Read(47) Lz_box_hit
    Read(47) Tu_hit, L11r_f, alphaL_f, A_vel_hit, alphat_f, &
             y_blend_hit, Uinf_f
    Allocate( y_hit_f(ny_hit_f) )
    Read(47) y_hit_f
    Close(47)

    ! data start offset (bytes, 1-based POS for stream access)
    hit_data_offset = 4_8*4_8 + 8_8*2_8 + 8_8*7_8 + 8_8*Int(ny_hit_f,8) + 1_8

    ! consistency checks
    ! true spanwise period of the code: dz*(nz_global-2)
    ! (z(1) and z(nz-1) are periodic images; z(nz) is a ghost)
    If ( Abs( Lz_box_hit - (z(2)-z(1))*Real(nz_global-2,8) ) > 1d-6*Lz_box_hit ) Then
       Write(*,*) 'ERROR: domain z-period=', (z(2)-z(1))*Real(nz_global-2,8), &
                  ' /= plane z-period=', Lz_box_hit
       Write(*,*) 'Planes are z-periodic over Lz_box: set boxsize Lz accordingly'
       Stop 'hit_file z-period mismatch'
    End If
    If ( y_hit_f(ny_hit_f) < yg_global(nyg_global-1) ) Then
       Write(*,*) 'WARNING: plane y range ', y_hit_f(ny_hit_f), &
                  ' < domain top ', yg_global(nyg_global-1), ' (will clamp)'
    End If

    If ( N_buffer_hit > n_planes_hit ) N_buffer_hit = n_planes_hit

    ! GPU-resident buffers on the staggered inlet grids
    Allocate( hit_buf_u(nyg, nzg, N_buffer_hit) )
    Allocate( hit_buf_v(ny , nzg, N_buffer_hit) )
    Allocate( hit_buf_w(nyg, nz , N_buffer_hit) )

    ! load first buffer (before the OpenACC data region: copied in by copyin)
    Call load_hit_buffer(1)

    If (myid==0) Then
       Write(*,*) '   planes: ', n_planes_hit, '  grid: ', ny_hit_f, 'x', nz_hit_f
       Write(*,*) '   dt_plane = ', dt_plane_hit, ' (', dt_plane_hit/Abs(CFL), ' steps)'
       Write(*,*) '   time covered = ', Real(n_planes_hit-1,8)*dt_plane_hit
       Write(*,*) '   buffer = ', N_buffer_hit, ' planes (', &
                  Real(N_buffer_hit,8)*Real(nyg*nzg+ny*nzg+nyg*nz,8)*8d0/1d9, ' GB)'
       Write(*,*) '   Tu = ', Tu_hit, '  y_blend = ', y_blend_hit
    End If

  End Subroutine init_hit_inflow

  !-------------------------------------------------------------------!
  ! Load planes [istart, istart+N_buffer-1] (1-based, clamped to end) !
  ! from disk and interpolate to the staggered inlet grids (CPU).     !
  !-------------------------------------------------------------------!
  Subroutine load_hit_buffer(istart)

    Integer(Int32), Intent(In) :: istart

    Integer(Int32) :: n_load, p, ip, j, k, j0, j1, k0, k1
    Integer(Int64) :: pos, rec_bytes
    Real   (Int64) :: yy, zz, fy, fz, wyj, wzk, dzf
    Real   (4), Allocatable, Dimension(:,:) :: raw_u, raw_v, raw_w
    ! precomputed y-weights for centers (yg) and faces (y)
    Integer(Int32), Allocatable, Dimension(:) :: jg0, jf0
    Real   (Int64), Allocatable, Dimension(:) :: wg, wf
    ! precomputed z-weights for centers (zg) and faces (z)
    Integer(Int32), Allocatable, Dimension(:) :: kc0, kc1, kf0v, kf1v
    Real   (Int64), Allocatable, Dimension(:) :: wc, wfz

    n_load = Min( N_buffer_hit, n_planes_hit - istart + 1 )
    ibuf_start_hit = istart
    rec_bytes = 3_8 * Int(ny_hit_f,8) * Int(nz_hit_f,8) * 4_8

    Allocate( raw_u(nz_hit_f, ny_hit_f) )   ! file order: iy slow, iz fast
    Allocate( raw_v(nz_hit_f, ny_hit_f) )
    Allocate( raw_w(nz_hit_f, ny_hit_f) )

    ! ---- y interpolation weights (clamped linear on y_hit_f) ----
    Allocate( jg0(nyg), wg(nyg), jf0(ny), wf(ny) )
    Do j = 1, nyg
       yy = Min( Max( yg(j), y_hit_f(1) ), y_hit_f(ny_hit_f) )
       j0 = 1
       Do j1 = 2, ny_hit_f
          If ( y_hit_f(j1) >= yy ) Then
             j0 = j1 - 1
             Exit
          End If
          If ( j1 == ny_hit_f ) j0 = ny_hit_f - 1
       End Do
       jg0(j) = j0
       wg (j) = ( yy - y_hit_f(j0) ) / ( y_hit_f(j0+1) - y_hit_f(j0) )
    End Do
    Do j = 1, ny
       yy = Min( Max( y(j), y_hit_f(1) ), y_hit_f(ny_hit_f) )
       j0 = 1
       Do j1 = 2, ny_hit_f
          If ( y_hit_f(j1) >= yy ) Then
             j0 = j1 - 1
             Exit
          End If
          If ( j1 == ny_hit_f ) j0 = ny_hit_f - 1
       End Do
       jf0(j) = j0
       wf (j) = ( yy - y_hit_f(j0) ) / ( y_hit_f(j0+1) - y_hit_f(j0) )
    End Do

    ! ---- z interpolation weights (periodic over Lz_box_hit) ----
    dzf = Lz_box_hit / Real(nz_hit_f,8)
    Allocate( kc0(nzg), kc1(nzg), wc(nzg), kf0v(nz), kf1v(nz), wfz(nz) )
    Do k = 1, nzg
       zz = Modulo( zg(k), Lz_box_hit )
       fz = zz / dzf
       k0 = Int(fz)
       wc(k)  = fz - Real(k0,8)
       kc0(k) = Mod(k0  , nz_hit_f) + 1
       kc1(k) = Mod(k0+1, nz_hit_f) + 1
    End Do
    Do k = 1, nz
       zz = Modulo( z(k), Lz_box_hit )
       fz = zz / dzf
       k0 = Int(fz)
       wfz(k)  = fz - Real(k0,8)
       kf0v(k) = Mod(k0  , nz_hit_f) + 1
       kf1v(k) = Mod(k0+1, nz_hit_f) + 1
    End Do

    ! ---- read + interpolate each plane ----
    Open(47, file=Trim(file_hit_planes), form='unformatted', &
         action='read', access='stream', status='old')
    Do p = 1, n_load
       ip  = istart + p - 1
       pos = hit_data_offset + Int(ip-1,8)*rec_bytes
       Read(47, POS=pos) raw_u, raw_v, raw_w

       ! u at (yg, zg)
       Do k = 1, nzg
          k0 = kc0(k); k1 = kc1(k); wzk = wc(k)
          Do j = 1, nyg
             j0 = jg0(j); wyj = wg(j)
             hit_buf_u(j,k,p) = &
               (1d0-wyj)*( (1d0-wzk)*Real(raw_u(k0,j0  ),8) + wzk*Real(raw_u(k1,j0  ),8) ) + &
                    wyj *( (1d0-wzk)*Real(raw_u(k0,j0+1),8) + wzk*Real(raw_u(k1,j0+1),8) )
          End Do
       End Do
       ! v at (y, zg)
       Do k = 1, nzg
          k0 = kc0(k); k1 = kc1(k); wzk = wc(k)
          Do j = 1, ny
             j0 = jf0(j); wyj = wf(j)
             hit_buf_v(j,k,p) = &
               (1d0-wyj)*( (1d0-wzk)*Real(raw_v(k0,j0  ),8) + wzk*Real(raw_v(k1,j0  ),8) ) + &
                    wyj *( (1d0-wzk)*Real(raw_v(k0,j0+1),8) + wzk*Real(raw_v(k1,j0+1),8) )
          End Do
       End Do
       ! w at (yg, z)
       Do k = 1, nz
          k0 = kf0v(k); k1 = kf1v(k); wzk = wfz(k)
          Do j = 1, nyg
             j0 = jg0(j); wyj = wg(j)
             hit_buf_w(j,k,p) = &
               (1d0-wyj)*( (1d0-wzk)*Real(raw_w(k0,j0  ),8) + wzk*Real(raw_w(k1,j0  ),8) ) + &
                    wyj *( (1d0-wzk)*Real(raw_w(k0,j0+1),8) + wzk*Real(raw_w(k1,j0+1),8) )
          End Do
       End Do
    End Do
    Close(47)

    If (myid==0) Write(*,'(a,i8,a,i8,a)') '   HIT buffer loaded: planes ', &
         istart, ' to ', istart+n_load-1, ' (interpolated to inlet grids)'

    Deallocate( raw_u, raw_v, raw_w, jg0, wg, jf0, wf )
    Deallocate( kc0, kc1, wc, kf0v, kf1v, wfz )

  End Subroutine load_hit_buffer

  !-------------------------------------------------------------------!
  ! Fill Ut/Vt/Wt_inlet from the HIT buffer at the current time t     !
  ! (linear interpolation between bracketing planes, GPU kernel).     !
  ! Reloads the buffer from disk when exhausted.                      !
  !-------------------------------------------------------------------!
  Subroutine update_hit_inlet_gpu

    Integer(Int32) :: ip0, ipw, p0, il0, j, k, n_loaded
    Real   (Int64) :: theta

    ! global plane interval (0-based). NO-WRAP CAMPAIGN (audit blueprint):
    ! the library must NEVER be recycled — running past its end is a fatal
    ! configuration error, not something to paper over with Mod().
    ip0   = Int( t / dt_plane_hit )
    theta = t/dt_plane_hit - Real(ip0,8)
    If ( ip0 >= n_planes_hit-1 ) Then
       Write(*,*) 'FATAL: HIT inflow library exhausted: t=', t, &
            ' ip0=', ip0, ' n_planes=', n_planes_hit
       Stop 'HIT library exhausted — no wrap allowed (extend library first)'
    End If
    ipw   = ip0
    p0    = ipw + 1                     ! 1-based file index of left plane

    ! reload if [p0, p0+1] not inside the buffer
    n_loaded = Min( N_buffer_hit, n_planes_hit - ibuf_start_hit + 1 )
    If ( p0 < ibuf_start_hit .or. p0+1 > ibuf_start_hit+n_loaded-1 ) Then
       Call load_hit_buffer(p0)
       !$acc update device(hit_buf_u, hit_buf_v, hit_buf_w)
    End If
    il0 = p0 - ibuf_start_hit + 1       ! local buffer index

    !$acc parallel loop collapse(2) default(present)
    Do k = 1, nzg
       Do j = 1, nyg
          Ut_inlet(j,k) = (1d0-theta)*hit_buf_u(j,k,il0) + theta*hit_buf_u(j,k,il0+1)
       End Do
    End Do
    !$acc parallel loop collapse(2) default(present)
    Do k = 1, nzg
       Do j = 1, ny
          Vt_inlet(j,k) = (1d0-theta)*hit_buf_v(j,k,il0) + theta*hit_buf_v(j,k,il0+1)
       End Do
    End Do
    !$acc parallel loop collapse(2) default(present)
    Do k = 1, nz
       Do j = 1, nyg
          Wt_inlet(j,k) = (1d0-theta)*hit_buf_w(j,k,il0) + theta*hit_buf_w(j,k,il0+1)
       End Do
    End Do

  End Subroutine update_hit_inlet_gpu

  !-------------------------------------------------------!
  !                                                       !
  !              Compute temporal component  for          !
  !                  boundary conditions                  !
  !                                                       !
  ! Ut_inlet = sum qu_inlet*exp( I*beta*n*z - I*omega*mt) !
  ! Vt_inlet = sum qv_inlet*exp( I*beta*n*z - I*omega*mt) !
  ! Wt_inlet = sum qw_inlet*exp( I*beta*n*z - I*omega*mt) !
  !                                                       !
  ! zmode_inlet = beta_inlet*n                            ! 
  ! tmode_inlet = omega_inlet*m                           !
  !                                                       !
  ! Input:  t                                             !
  ! Output: Ut_inlet, Vt_inlet, Wt_inlet                  !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine compute_temporal_inflow

    ! local variables
    Integer(Int32) :: j, k, n, m
    Real   (Int64) :: t_inlet
    Complex(Int64) :: Iu

    t_inlet = t         ! time    
    Iu      = (0d0,1d0) ! imaginary unit

    !----------------------------------------------------!
    ! Ut
    Ut_inlet = 0d0
    Do j = 1, nyg
       Do k = 1, nzg
          !
          Do n = 1, n_modes_inlet ! z-mode
             Do m = 1, m_modes_inlet ! t-mode
                Ut_inlet(j,k) = Ut_inlet(j,k) + & 
                Real( qu_inlet(j,n,m)*cdexp(Iu*zmode_inlet(n)*zg(k) - Iu*tmode_inlet(m)*t_inlet ) , 8)
             End Do
          End Do
          !
       End Do
    End Do

    !----------------------------------------------------!
    ! Vt
    Vt_inlet = 0d0
    Do j = 1, ny
       Do k = 1, nzg
          !
          Do n = 1, n_modes_inlet ! z-mode
             Do m = 1, m_modes_inlet ! t-mode
                Vt_inlet(j,k) = Vt_inlet(j,k) + & 
                Real( qv_inlet(j,n,m)*cdexp(Iu*zmode_inlet(n)*zg(k) - Iu*tmode_inlet(m)*t_inlet ) , 8)
             End Do
          End Do
          !
       End Do
    End Do

    !----------------------------------------------------!
    ! Wt
    Wt_inlet = 0d0
    Do j = 1, nyg
       Do k = 1, nz
          !
          Do n = 1, n_modes_inlet ! z-mode
             Do m = 1, m_modes_inlet ! t-mode
                Wt_inlet(j,k) = Wt_inlet(j,k) + & 
                Real( qw_inlet(j,n,m)*cdexp(Iu*zmode_inlet(n)*z(k) - Iu*tmode_inlet(m)*t_inlet ) , 8)
             End Do
          End Do
          !
       End Do
    End Do

  End Subroutine compute_temporal_inflow

  !-------------------------------------------------------------!
  !                                                             !
  !            Compute Neumann boundary conditions              !
  !    for pseudo-pressure when slip-wall model is active       !
  !                                                             !
  ! This has to be called every sub-step                        !
  !                                                             !
  ! Conditions top wall: (n->ny, ng->nyg)                       !
  !                                                             !
  !     V(n)   = V(  n)* -(p(ng)  -p(ng-1))/(yg(ng)  -yg(ng-1)) !
  !     V(n-1) = V(n-1)* -(p(ng-1)-p(ng-2))/(yg(ng-1)-yg(ng-2)) !
  !     V(n)   = alpha_y*(V(n)-V(n-1))/(y(n)-y(n-1))            !
  !     alpha_y -> Inf                                          !
  !                                                             !
  !     => p(ng)    = p_bn1*p(ng-2) + p_bcn*p(ng-1)             !
  !        p_bcn    = 1d0 + beta*Delta_r                        !
  !        p_bcn1   =     - beta*Delta_r                        !
  !        alpha_y  = - alpha_y -> same alpha_y each wall!      !
  !        Delta_r  = (yg(ng)-yg(ng-1))/(yg(ng-1)-yg(ng-2))     !
  !        beta     = alpha_y/( alpha_y - (y(n)-y(n-1)) )       !
  !                                                             !
  ! Equation for last interior points:                          !
  !                                                             !
  !   (b+a*p_bcn)*p(ng-1) + (c+a*p_bcn1)*p(ng-2) = rhs_p2(ng-1) !
  !                                                             !
  !                                                             !
  ! Assumed:                                                    !
  !      - alpha_y is not function of (x,z)                     !
  !      - V1* = alpha_y*(V2*-V1*)/(y(2)-y(1))                  !
  !        Otherwise the rhs for P must be modified             !
  !                                                             !
  !-------------------------------------------------------------!
  Subroutine compute_pseudo_pressure_bc_for_top_Neumann_bc

    ! local variables
    Real(Int64) :: Delta_r

    ! top wall
    Delta_r = ( yg(nyg) - yg(nyg-1) )/( yg(nyg-1) - yg(nyg-2) )

    Dyy(nyg-1,nyg-1) =  1d0 + Delta_r
    Dyy(nyg-1,nyg-2) = -Delta_r

  End Subroutine compute_pseudo_pressure_bc_for_top_Neumann_bc

End Module boundary_conditions
