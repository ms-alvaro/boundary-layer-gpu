!------------------------------------!
!   Module for boundary conditions   !
!------------------------------------!
Module boundary_conditions

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use Lund_rescaled_bc
  Use mpi
  Use ifport

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
    If ( top_boundary_flag.ne.0 ) Then
       ! U and W boundary condition at the top           
       Call apply_top_bc_y(U,W) 
       ! V boundary condition at the top    
       Call apply_Dirichlet_bc_y_top_BlowingSuction(V,top_boundary_flag)
    Else
       Call apply_top_bc_y_Falkner_Skan(U,V,W) 
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

    If ( step_beginning == 1 ) Then

       ! compute temporal component of the inflow
       Call compute_temporal_inflow
       
       ! variables at centers 
       Do j=1,nyg
          Do k=1,nzg
             U_(1,j,k) = U_inlet(j) + Amplitude_perturbations*(rand()-0.5d0)
          End Do
       End Do
       Do j=1,nyg
          Do k=1,nz
             W_(1,j,k) = W_inlet(j) + Amplitude_perturbations*(rand()-0.5d0)
          End Do
       End Do
       ! variables at faces
       Do j=1,ny
          Do k=1,nzg
             V_(1,j,k) = V_inlet(j) + Amplitude_perturbations*(rand()-0.5d0)
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
    Real   (Int64) :: Amp_per
    Integer(Int32) :: j, k

    If ( step_beginning == 1 ) Then

       ! compute temporal component of the inflow
       !Call compute_temporal_inflow

       ! perturbation 
       Amp_per = Amplitude_perturbations
       
       ! variables at centers 
       Do j=1,nyg
          Do k=1,nzg
             U_(1,j,k) = U_inlet(j) + Amp_per*(rand()-0.5d0)*(U_inlet(ny)-U_inlet(j))
          End Do
       End Do
       Do j=1,nyg
          Do k=1,nz
             W_(1,j,k) = W_inlet(j) + Amp_per*(rand()-0.5d0)*(U_inlet(ny)-U_inlet(j))
          End Do
       End Do
       ! variables at faces
       Do j=1,ny
          Do k=1,nzg
             V_(1,j,k) = V_inlet(j) + Amp_per*(rand()-0.5d0)*(U_inlet(ny)-U_inlet(j))
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
    Integer(Int32) :: j, k

    If ( step_beginning == 1 ) Then

       ! compute temporal component of the inflow
       !Call compute_temporal_inflow
       Call compute_rescaled_inflow(Uo_,Vo_,Wo_,Ut_inlet,Vt_inlet,Wt_inlet,delta_inlet,1)
       
       ! variables at centers 
       Do j=1,nyg
          Do k=1,nzg
             U_(1,j,k) = U_inlet(j) + Ut_inlet(j,k)
          End Do
       End Do
       Do j=1,nyg
          Do k=1,nz
             W_(1,j,k) = W_inlet(j) + Wt_inlet(j,k)
          End Do
       End Do
       ! variables at faces
       Do j=1,ny
          Do k=1,nzg
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
          F(i,ny,:) = Vbs_max*(x_bs-xg(i))/sigma_bs*dexp(.5d0 - ((x_bs-xg(i))/sigma_bs)**2d0 ) &
                      + phi_bs
       End Do     
    Elseif ( id==2 ) Then ! Abe 2017
       ! F defined at y faces
       Do i=1,nxg
          F(i,ny,:) = Vbs_max*(x_bs-xg(i))/sigma_bs*dexp(phi_bs - ((x_bs-xg(i))/sigma_bs)**2d0 )
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
  Subroutine apply_top_bc_y(U_,W_)

    Real(Int64), Dimension(:,:,:), Intent(InOut) :: U_, W_

    ! local variables
    Integer(Int32) :: i

    ! variables at centers
    Do i=1,nx
       U_(i,nyg,:) = U_top(i)
    End Do
    ! variables at faces
    Do i=1,nxg
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
          Open(5,file=file_blasius_own,form='unformatted',action='read',access='stream')
          Read(5) ny_profile
          If (ny_profile/=ny_global ) Then
             Write(*,*) 'file_blasius_own',file_blasius_own
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
       Open(4,file=file_blasius,form='formatted',action='read')
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
          
       If ( inflow_boundary_flag==1 .or. inflow_boundary_flag==3 .or. inflow_boundary_flag==4 ) Then

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
          Open(5,file=file_blasius_own,form='unformatted',action='read',access='stream')
          Read(5) ny_blasius
          If (ny_blasius/=ny_global ) Then 
             Write(*,*) 'file_blasius_own',file_blasius_own
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
