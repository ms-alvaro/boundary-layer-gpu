!----------------------------------------!
!    Module for subgrid-scale models     ! 
!----------------------------------------!
Module subgrid

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global 
  Use interpolation
  Use boundary_conditions
  Use mpi

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------------------!
  !                                                      !
  !             Select eddy viscosity model              !
  !                                                      ! 
  !------------------------------------------------------!
  Subroutine compute_eddy_viscosity(U_,V_,W_,avg_nu_t_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,    1), Intent(Out) :: avg_nu_t_
    Real(Int64), Dimension(nxg,nyg,nzm+2), Intent(Out) :: nu_t_

    If     ( LES_model==1 ) Then
       ! constant coef. Smagorinsky
       Call sgs_Smagorinsky(U_,V_,W_,avg_nu_t_,nu_t_)    
    Elseif ( LES_model==2 ) Then
       ! dynamic Smagorinsky (z-averaged)
       Call sgs_dynamic_Smagorinsky(U_,V_,W_,avg_nu_t_,nu_t_) 
    Elseif ( LES_model==3 ) Then
       ! dynamic Smagorinsky (no average)
       Call sgs_dynamic_Smagorinsky_noaverage(U_,V_,W_,avg_nu_t_,nu_t_) 
    Else
       ! No explicit model
       nu_t_     = 0d0
       avg_nu_t_ = 0d0
    End If

  End Subroutine compute_eddy_viscosity

  !--------------------------------------------------------!
  !      Filter field in homogeneous directions xz         !
  !                                                        !
  ! Simpson's rule (4th order):                            !
  !                                                        !
  !   Uf(0,0) = sum( fil(i,j)*U(i,j) )                     !
  !                                                        !
  ! Weights fil:                                           !
  !   1/36  1/9  1/36                                      !
  !   1/9   4/9  1/9                                       !
  !   1/36  1/9  1/36                                      !
  !                                                        !
  ! Note:                                                  !
  !   Must be called as filter_xz(U,Uf(2:n1-1,:,2:n3-1))   !
  !                                                        !
  !                                                        !
  ! Input:  U   (original flow field)                      !
  ! Output: Uf  (filtered flow field)                      !
  !                                                        ! 
  !--------------------------------------------------------!
  Subroutine filter_xz(U_,Uf_)

    Real(Int64), Intent(In)  :: U_  (:,:,:)
    Real(Int64), Intent(Out) :: Uf_ (:,:,:) 

    ! local variables
    Real   (Int64) :: fil(-1:1,-1:1)
    Integer(Int64) :: n(3), n1, n2, n3, i, j, k, ii, kk

    ! local size
    n  = Shape(U_)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ! filter size
    fil_size = 1.5874d0

    ! Set up filter coefficients
    fil(-1,-1) = 1d0/36d0
    fil( 1, 1) = 1d0/36d0
    fil(-1, 1) = 1d0/36d0
    fil( 1,-1) = 1d0/36d0

    fil( 0, 1) = 1d0/9d0
    fil( 1, 0) = 1d0/9d0
    fil( 0,-1) = 1d0/9d0
    fil(-1, 0) = 1d0/9d0

    fil( 0, 0) = 4d0/9d0

    ! Apply filter in homogeneous directions xz
    Uf_ = 0d0 
    ! loop for Uf(ii,:,kk)
    Do ii = 2, n1-1
       Do kk = 2, n3-1
          ! loop for filter
          Do i = -1, 1
             Do k = -1, 1
                ! the -1 shift is because filter_xz is called as filter_xz(U,Uf(2:end-1,:,2:end-1))
                Uf_(ii-1,1:n2,kk-1) = Uf_(ii-1,1:n2,kk-1) + U_(ii+i,1:n2,kk+k)*fil(i,k)
             End do
          End do
       End do
    End do
    
  End Subroutine filter_xz

  !--------------------------------------------------------!
  !             Filter field in xyz directions             !
  !                                                        !
  ! Simpson's rule (4th order):                            !
  !                                                        !
  !   Uf(0,0) = sum( fil(i,j)*U(i,j) )                     !
  !                                                        !
  ! Weights fil in x and z:                                !
  !   1/36  1/9  1/36                                      !
  !   1/9   4/9  1/9                                       !
  !   1/36  1/9  1/36                                      !
  !                                                        !
  ! Weights fil in y:                                      !
  !   1/6   2/3  1/6                                       !
  !                                                        !
  ! Note:                                                  !
  !   Must be called as filter_xzy(U,Uf(2:n1-1,:,2:n3-1))  !
  !                                                        !
  !                                                        !
  ! Input:  U   (original flow field)                      !
  ! Output: Uf  (filtered flow field)                      !
  !                                                        ! 
  !--------------------------------------------------------!
  Subroutine filter_xzy(U_,Uf_)

    Real(Int64), Intent(In)  :: U_  (:,:,:)
    Real(Int64), Intent(Out) :: Uf_ (:,:,:)

    ! local variables
    Real   (Int64) :: fil(-1:1,-1:1)
    Integer(Int64) :: n(3), n1, n2, n3, i, j, k, ii, kk

    ! local size
    n  = Shape(U_)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ! filter size
    fil_size = 2d0

    ! Set up filter coefficients 
    fil(-1,-1) = 1d0/36d0
    fil( 1, 1) = 1d0/36d0
    fil(-1, 1) = 1d0/36d0
    fil( 1,-1) = 1d0/36d0

    fil( 0, 1) = 1d0/9d0
    fil( 1, 0) = 1d0/9d0
    fil( 0,-1) = 1d0/9d0
    fil(-1, 0) = 1d0/9d0

    fil( 0, 0) = 4d0/9d0

    !Apply filter in homogeneous directions
    Uf_ = 0d0
    Do ii = 2, n1-1
       Do i = -1, 1
          Do kk = 2, n3-1
             Do k = -1,1
                Do j = 2, n2-1
                   Uf_(ii-1,j,kk-1) = Uf_(ii-1,j,kk-1) + (1d0/6d0 * U_(ii+i,j-1,kk+k) + 2d0/3d0 * U_(ii+i,j,kk+k) + 1d0/6d0 * U_(ii+i,j+1,kk+k)) * fil(i,k)
                End Do
                Uf_(ii-1,1,kk-1) = Uf_(ii-1,1,kk-1) + (2d0/3d0 * U_(ii+i,1,kk+k) + 1d0/3d0 * U_(ii+i,2,kk+k)) * fil(i,k)
                Uf_(ii-1,n2,kk-1) = Uf_(ii-1,n2,kk-1) + (1d0/3d0 * U_(ii+i,n2-1,kk+k) + 2d0/3d0 * U_(ii+i,n2,kk+k)) * fil(i,k)
             End do
          End do
       End do
    End do

  End Subroutine filter_xzy

  !------------------------------------------------------!
  !      Filter tensor in homogeneous directions xz      !
  !                                                      !
  ! Input:  T  (4th component is tensor position)        !
  ! Output: Tf                                           !
  !                                                      !
  !------------------------------------------------------!
  Subroutine filter_tensor_xz(T_,Tf_)

    Real(Int64), Intent(In)  :: T_  (:,:,:,:)
    Real(Int64), Intent(Out) :: Tf_ (:,:,:,:) 

    ! local variables
    Integer(Int64) :: n(4), n4, i

    ! local size
    n  = Shape(T_)
    n4 = n(4)

    ! filtering
    Do i = 1, n4
      Call filter_xz(T_(:,:,:,i),Tf_(:,:,:,i))
    End do 

  End Subroutine filter_tensor_xz

  !--------------------------------------------------------!
  !      Filter tensor in homogeneous directions xz and y  !
  !                                                        !
  ! Input:  T  (4th component is tensor position)          !
  ! Output: Tf                                             !
  !                                                        !
  !--------------------------------------------------------!
  Subroutine filter_tensor_xzy(T_,Tf_)

    Real(Int64), Intent(In)  :: T_  (:,:,:,:)
    Real(Int64), Intent(Out) :: Tf_ (:,:,:,:)

    ! local variables
    Integer(Int64) :: n(4), n4, i

    ! local size
    n  = Shape(T_)
    n4 = n(4)

    ! filtering
    Do i = 1, n4
      Call filter_xzy(T_(:,:,:,i),Tf_(:,:,:,i))
    End do

  End Subroutine filter_tensor_xzy

  !-----------------------------------------------------------!
  !                                                           !
  !      Compute dynamic Smagorinsky eddy-viscosity           !
  !                                                           !
  ! nu_t = -0.5 * <(Lij * Mij) / (Mij * Mij) * S>  (clipping) !
  ! Lij  = hat(ui*uj) - hat(ui)*hat(uj)                       !
  ! Mij  = fil_size^2 * |hat(S)| * hat(Sij) - hat(|S| * Sij)  !
  ! |S|  = sqrt(2 Sij Sij)                                    !
  ! fil_size = (2*dx2*dz*dy)^(1/3)/(dx*dz*dy)^(1/3)           !
  !                                                           !
  ! Tensors are organized in arrays as                        !
  !   ( 1 4 5 )                                               !
  !   ( 4 2 6 )                                               !
  !   ( 5 6 3 )                                               !
  ! where the number is the 4th component of the array        !
  !                                                           !
  ! Input:  U_,V_,W_ (velocities)                             !
  ! Output: nu_t, avg_nu_t (eddy viscosity and average in z)  !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine sgs_dynamic_Smagorinsky(U_,V_,W_,avg_nu_t_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,    1), Intent(Out) :: avg_nu_t_
    Real(Int64), Dimension(nxg,nyg,nzm+2), Intent(Out) :: nu_t_

    ! local variables
    Integer(Int32) :: i, j, k

    !------------------------------------------------------------------!
    ! Part 1: Compute Leonard term Lij = hat(ui*uj) - hat(ui)*hat(uj)
    ! interpolate velocity to cell centers (faces to centers)

    Call interpolate_x(U_,term  (2:nxg-1,:,:),1) 
    Call interpolate_y(V_,term_1(:,2:nyg-1,:),1) 
    Call interpolate_z(W_,term_2(:,:,2:nzg),1) 

    ! fill in missing values (periodicity)
    Call apply_Neumann_bc_x (term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! Lij (Sij as placeholder) = ui*uj at cell centers (why to nzg?)
    Sij(:,2:nyg-1,:,1) = term  (2:nxg,2:nyg-1,2:nzg) * term  (2:nxg,2:nyg-1,2:nzg)  ! u^2
    Sij(:,2:nyg-1,:,2) = term_1(2:nxg,2:nyg-1,2:nzg) * term_1(2:nxg,2:nyg-1,2:nzg)  ! v^2
    Sij(:,2:nyg-1,:,3) = term_2(2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! w^2
    Sij(:,2:nyg-1,:,4) = term  (2:nxg,2:nyg-1,2:nzg) * term_1(2:nxg,2:nyg-1,2:nzg)  ! uv
    Sij(:,2:nyg-1,:,5) = term  (2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! uw
    Sij(:,2:nyg-1,:,6) = term_1(2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! vw

    ten_buf(2:nxg,2:nyg,2:nzg,:) = Sij;
    Do i = 1,6
       Call apply_Neumann_bc_x (ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do
    ! Lij = hat(ui*uj)
    Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),Lij(2:nxg-1,2:nyg-1,2:nzg-1,:)) ! this stores filtered values in Lij(3:nxg-1,2:nyg-1,3:nzg-1)


    ! filter interpolated velocity
    Call filter_xzy(term  (1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,1)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,1)
    Call filter_xzy(term_1(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)
    Call filter_xzy(term_2(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,3)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)

    ! Lij = hat(ui*uj) - hat(ui)*hat(uj)
    ! Only valid on (2:nxg-1,2:nyg-1,2:nzg-1)
    Lij(:,:,:,1) = Lij(:,:,:,1) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,1) ! hat(u^2) - u^2
    Lij(:,:,:,2) = Lij(:,:,:,2) - Sij(2:nxg,2:nyg-1,2:nzg,2) * Sij(2:nxg,2:nyg-1,2:nzg,2) ! hat(v^2) - v^2
    Lij(:,:,:,3) = Lij(:,:,:,3) - Sij(2:nxg,2:nyg-1,2:nzg,3) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat(w^2) - w^2
    Lij(:,:,:,4) = Lij(:,:,:,4) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,2) ! hat( uv) - uv
    Lij(:,:,:,5) = Lij(:,:,:,5) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat( uw) - uw
    Lij(:,:,:,6) = Lij(:,:,:,6) - Sij(2:nxg,2:nyg-1,2:nzg,2) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat(u^2) - u^2

    !------------------------------------------------------------------!
    ! Part 2: Compute Mij = fil_size^2*|hat(S)|*hat(Sij) - hat(|S|*Sij) 

    ! filter velocity Same indices as U_
    Call filter_xzy(U_(1:nx,1:nyg,1:nzg),Uf(2:nx-1,1:nyg,2:nzg-1)) ! this stores filtered values in Uf(2:nx-1,2:nyg-1,2:nzg-1)
    Call filter_xzy(V_(1:nxg,1:ny,1:nzg),Vf(2:nxg-1,1:ny,2:nzg-1)) ! this stores filtered values in Vf(2:nxg-1,1:ny,2:nzg-1)
    Call filter_xzy(W_(1:nxg,1:nyg,1:nz),Wf(2:nxg-1,1:nyg,2:nz-1)) ! this stores filtered values in Wf(2:nxg-1,2:nyg-1,2:nz-1)

    ! fill in missing values (periodicity)
    Call apply_Neumann_bc_x (Uf,1)
    Call apply_periodic_bc_z(Uf,1)
    Call update_ghost_interior_planes(Uf,1)

    Call apply_Neumann_bc_x (Vf,2)
    Call apply_periodic_bc_z(Vf,2)
    Call update_ghost_interior_planes(Vf,2)

    Call apply_Neumann_bc_x (Wf,2)
    Call apply_periodic_bc_z(Wf,3)
    Call update_ghost_interior_planes(Wf,3)

    ! Compute hat(Sij) with filtered velocities
    Call compute_Sij(Uf,Vf,Wf,Sij,S)

    ! Mij = fil_size^2 * |hat(S)| * hat(Sij)
    Do i = 1,6
      Mij(:,:,:,i) = fil_size**2d0 * S * Sij(2:nxg-1,2:nyg-1,2:nzg-1,i)
    End Do

    ! Compute Sij with unfiltered velocities
    Call compute_Sij(U_,V_,W_,Sij,S)

    ! Compute hat(|S|Sij) and subtract from Mij
    Do i = 1,6
      term(2:nxg-1,2:nyg-1,2:nzg-1) = S * Sij(2:nxg-1,2:nyg-1,2:nzg-1,i)
      call apply_Neumann_bc_x (term,2)
      call apply_periodic_bc_z(term,4)
      call update_ghost_interior_planes(term,4)
      Call filter_xzy(term(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,i))
    End Do


    ! Mij = fil_size^2*|hat(S)|*hat(Sij) - hat(S*Sij)
    Mij = Mij - Sij(2:nxg-1,2:nyg-1,2:nzg-1,:)

    !------------------------------------------------------------------!
    ! Part 3: Compute eddy viscosity nu_t = Mij*Lij / Mij*Mij

    nu_t_ = 0d0
    Do i = 1,6
      If (i .le. 3) Then 
        nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) + Mij(:,:,:,i) * Lij(2:nxg-1,:,2:nzg-1,i)
      Else 
        nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) + 2d0 * Mij(:,:,:,i) * Lij(2:nxg-1,:,2:nzg-1,i)
      End if
    End Do

    ! compute <Mij*Lij>
    Do i = 2,nxg-1
      Do j = 2,nyg-1
         avg_nu_t_(i,j,1) = Sum( nu_t_(i,j,2:nzg-2) ) / Real((nzg-3), 8)
      End Do
    End Do

    nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = ( Mij(:,:,:,1)*Mij(:,:,:,1)  +      &
                                       Mij(:,:,:,2)*Mij(:,:,:,2)  +      & 
                                       Mij(:,:,:,3)*Mij(:,:,:,3)  +      & 
                                 2d0*( Mij(:,:,:,4)*Mij(:,:,:,4)  +      & 
                                       Mij(:,:,:,5)*Mij(:,:,:,5)  +      & 
                                       Mij(:,:,:,6)*Mij(:,:,:,6) ) )

    Do i = 2,nxg-1
      Do j = 2,nyg-1
         avg_nu_t_(i,j,1) = avg_nu_t_(i,j,1) / (Sum (nu_t_(i,j,2:nzg-2)) /Real((nzg-3),8) )
      End Do
    End Do

    Do i = 2,nxg-1
      Do j = 2,nyg-1
        nu_t_(i,j,2:nzg-1) = -0.5d0 * avg_nu_t_(i,j,1) * S(i,j,2:nzg-1) 
      End Do
    End Do

    Call update_ghost_interior_planes(nu_t_,4)


    ! clipping negative values
    nu_t_ = Max( nu_t_, 0d0 )

    ! boundary conditions ghost cell (must be done after clipping)
    If ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(nu_t_,2)
    Else
       Call apply_Neumann_bc_y  (nu_t_,2)
    End If

  End Subroutine sgs_dynamic_Smagorinsky

  !-----------------------------------------------------------!
  !                                                           !
  !      Compute dynamic Smagorinsky eddy-viscosity           !
  !                                                           !
  ! nu_t = -0.5 * (Lij * Mij) / (Mij * Mij) * S  (clipping)   !
  ! Lij  = hat(ui*uj) - hat(ui)*hat(uj)                       !
  ! Mij  = fil_size^2 * |hat(S)| * hat(Sij) - hat(|S| * Sij)  !
  ! |S|  = sqrt(2 Sij Sij)                                    !
  ! fil_size = (2*dx2*dz*dy)^(1/3)/(dx*dz*dy)^(1/3)           !
  !                                                           !
  ! Tensors are organized in arrays as                        !
  !   ( 1 4 5 )                                               !
  !   ( 4 2 6 )                                               !
  !   ( 5 6 3 )                                               !
  ! where the number is the 4th component of the array        !
  !                                                           !
  ! Input:  U_,V_,W_ (velocities)                             !
  ! Output: nu_t, avg_nu_t (eddy viscosity and average in xz) !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine sgs_dynamic_Smagorinsky_noaverage(U_,V_,W_,avg_nu_t_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(  1,nyg,  1), Intent(Out) :: avg_nu_t_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(Out) :: nu_t_

    ! local variables
    Integer(Int32) :: i, j, k

    !------------------------------------------------------------------!
    ! Part 1: Compute Leonard term Lij = hat(ui*uj) - hat(ui)*hat(uj)

    ! interpolate velocity to cell centers (faces to centers)
    Call interpolate_x(U_,term  (2:nxg-1,:,:),1) 
    Call interpolate_y(V_,term_1(:,2:nyg-1,:),1) 
    Call interpolate_z(W_,term_2(:,:,2:nzg  ),1) 

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! Lij (Sij as placeholder) = ui*uj at cell centers (why to nzg?)
    Sij(:,2:nyg-1,:,1) = term  (2:nxg,2:nyg-1,2:nzg) * term  (2:nxg,2:nyg-1,2:nzg)  ! u^2
    Sij(:,2:nyg-1,:,2) = term_1(2:nxg,2:nyg-1,2:nzg) * term_1(2:nxg,2:nyg-1,2:nzg)  ! v^2
    Sij(:,2:nyg-1,:,3) = term_2(2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! w^2
    Sij(:,2:nyg-1,:,4) = term  (2:nxg,2:nyg-1,2:nzg) * term_1(2:nxg,2:nyg-1,2:nzg)  ! uv
    Sij(:,2:nyg-1,:,5) = term  (2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! uw
    Sij(:,2:nyg-1,:,6) = term_1(2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! vw

    ten_buf(2:nxg,2:nyg,2:nzg,:) = Sij
    Do i = 1,6
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do

    ! Lij = hat(ui*uj)
    Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),Lij(2:nxg-1,2:nyg-1,2:nzg-1,:)) ! this stores filtered values in Lij(3:nxg-1,2:nyg-1,3:nzg-1)

!! fill in missing values ???
!    Lij(2,:,:,:) = Lij(nxg-1,:,:,:)
!    Lij(:,:,2,:) = Lij(:,:,nzg-1,:)

    ! filter interpolated velocity
    Call filter_xzy(term  (1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,1)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,1)
    Call filter_xzy(term_1(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)
    Call filter_xzy(term_2(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,3)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)

    ! Lij = hat(ui*uj) - hat(ui)*hat(uj)
    ! Only valid on (2:nxg-1,2:nyg-1,2:nzg-1)
    Lij(:,:,:,1) = Lij(:,:,:,1) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,1) ! hat(u^2) - u^2
    Lij(:,:,:,2) = Lij(:,:,:,2) - Sij(2:nxg,2:nyg-1,2:nzg,2) * Sij(2:nxg,2:nyg-1,2:nzg,2) ! hat(v^2) - v^2
    Lij(:,:,:,3) = Lij(:,:,:,3) - Sij(2:nxg,2:nyg-1,2:nzg,3) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat(w^2) - w^2
    Lij(:,:,:,4) = Lij(:,:,:,4) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,2) ! hat( uv) - uv
    Lij(:,:,:,5) = Lij(:,:,:,5) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat( uw) - uw
    Lij(:,:,:,6) = Lij(:,:,:,6) - Sij(2:nxg,2:nyg-1,2:nzg,2) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat(u^2) - u^2

    !------------------------------------------------------------------!
    ! Part 2: Compute Mij = fil_size^2*|hat(S)|*hat(Sij) - hat(|S|*Sij) 

    ! filter velocity Same indices as U_
    Call filter_xzy(U_(1:nx,1:nyg,1:nzg),Uf(2:nx-1,1:nyg,2:nzg-1)) ! this stores filtered values in Uf(2:nx-1,2:nyg-1,2:nzg-1)
    Call filter_xzy(V_(1:nxg,1:ny,1:nzg),Vf(2:nxg-1,1:ny,2:nzg-1)) ! this stores filtered values in Vf(2:nxg-1,1:ny,2:nzg-1)
    Call filter_xzy(W_(1:nxg,1:nyg,1:nz),Wf(2:nxg-1,1:nyg,2:nz-1)) ! this stores filtered values in Wf(2:nxg-1,2:nyg-1,2:nz-1)

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(Uf,1)
    Call apply_periodic_bc_z(Uf,1)
    Call update_ghost_interior_planes(Uf,1)

    Call apply_periodic_bc_x(Vf,2)
    Call apply_periodic_bc_z(Vf,2)
    Call update_ghost_interior_planes(Vf,2)

    Call apply_periodic_bc_x(Wf,2)
    Call apply_periodic_bc_z(Wf,3)
    Call update_ghost_interior_planes(Wf,3)

    ! Compute hat(Sij) with filtered velocities
    Call compute_Sij(Uf,Vf,Wf,Sij,S)

    ! Mij = fil_size^2 * |hat(S)| * hat(Sij)
    Do i = 1,6
      Mij(:,:,:,i) = fil_size**2d0 * S * Sij(2:nxg-1,2:nyg-1,2:nzg-1,i)
    End Do

    ! Compute Sij with unfiltered velocities
    Call compute_Sij(U_,V_,W_,Sij,S)

    ! Compute hat(|S|Sij) and subtract from Mij
    Do i = 1,6
      term(2:nxg-1,2:nyg-1,2:nzg-1) = S * Sij(2:nxg-1,2:nyg-1,2:nzg-1,i)
      Call apply_periodic_bc_x(term,2)
      Call apply_periodic_bc_z(term,4)
      Call update_ghost_interior_planes(term,4)
      Call filter_xzy(term(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,i))
    End Do

    ! Mij = fil_size^2*|hat(S)|*hat(Sij) - hat(S*Sij)
    Mij = Mij - Sij(2:nxg-1,2:nyg-1,2:nzg-1,:)

    !------------------------------------------------------------------!
    ! Part 3: Compute eddy viscosity nu_t = -0.5*(Lij*Mij)/(Mij*Mij)*S

    nu_t_ = 0d0
    Do i = 1,6
      If (i .le. 3) Then 
        nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) + Mij(:,:,:,i) * Lij(2:nxg-1,:,2:nzg-1,i)
      Else 
        nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) + 2d0 * Mij(:,:,:,i) * Lij(2:nxg-1,:,2:nzg-1,i)
      End if
    End Do
    nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) / & 
                                    ( Mij(:,:,:,1)*Mij(:,:,:,1)  +    &
                                      Mij(:,:,:,2)*Mij(:,:,:,2)  +    & 
                                      Mij(:,:,:,3)*Mij(:,:,:,3)  +    & 
                                2d0*( Mij(:,:,:,4)*Mij(:,:,:,4)  +    & 
                                      Mij(:,:,:,5)*Mij(:,:,:,5)  +    & 
                                      Mij(:,:,:,6)*Mij(:,:,:,6) ) )

    ! compute -0.5*Lij*Mij/Mij*Mij*|S|
    Do j = 2,nyg-1
       nu_t_(2:nxg-1,j,2:nzg-1) = -0.5d0 * nu_t_(2:nxg-1,j,2:nzg-1) * S(2:nxg-1,j,2:nzg-1)
    End Do

    ! apply boundary conditions in x and z
    Call apply_periodic_bc_x(nu_t_,4)
    Call apply_periodic_bc_z(nu_t_,4)
    Call update_ghost_interior_planes(nu_t_,4)

    ! clipping negative values
    nu_t_ = Max( nu_t_, 0d0 )

    ! boundary conditions ghost cell (must be done after clipping)
    If ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(nu_t_,2)
    Else
       Call apply_Neumann_bc_y  (nu_t_,2)
    End If

    ! not used here
    avg_nu_t_ = 0d0

  End Subroutine sgs_dynamic_Smagorinsky_noaverage

  !-----------------------------------------------------------!
  !                                                           !
  !  Compute constant coefficient Smagorinsky eddy-viscosity  !
  !       with van Driest damping function at the wall        !
  !                                                           !
  ! nu_t  = (Cs*Delta*f)^2 |S|                                !
  ! Cs    = 0.11 (usual range 0.1-0.2 )                       !
  ! Delta = (Cell V_lume)^1/3                                 !
  ! |S|   = sqrt(2 Sij Sij)                                   !
  !                                                           !
  ! van Driest damping function:                              !
  ! f = 1 - exp(-y+/25)                                       !
  ! + is wall-units computed from pressure gradient           !
  !                                                           !
  ! Tensors are organized in arrays as                        !
  !   ( 1 4 5 )                                               !
  !   ( 4 2 6 )                                               !
  !   ( 5 6 3 )                                               !
  ! where the number denotes the 4th component of the array   !
  !                                                           !
  ! Input:  U_,V_,W_ (velocities)                             !
  ! Output: nu_t, avg_nu_t (eddy viscosity and average in xz) !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine sgs_Smagorinsky(U_,V_,W_,avg_nu_t_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,    1), Intent(Out) :: avg_nu_t_
    Real(Int64), Dimension(nxg,nyg,nzm+2), Intent(Out) :: nu_t_

    ! local variables
    Real   (Int64) :: dum, maxerr, Cs, Delta, f, utau_t
    Integer(Int32) :: i, j, k

    ! Smagorinsky constant coefficient
    Cs = 0.11d0
    
    ! utau from pressure gradient
    utau_t = dPdx**0.5d0

    !------------------------------------------------------------------!
    ! Part 1: Compute rate-of-strain tensor Sij

    ! Compute rate-of-strain and S at cell centers
    Call compute_Sij(U_,V_,W_,Sij,S)

    !------------------------------------------------------------------!
    ! Part 2: Compute eddy viscosity nu_t = (Cs*Delta*f)^2|S| at cell centers
    Do j = 2, nyg-1
       If ( j<nyg/2 ) Then
          f = 1d0 - dexp( -( yg(j)-y(1)  )/25d0/(nu/utau_t) )
       Else
          f = 1d0 - dexp( -( y(ny)-yg(j) )/25d0/(nu/utau_t) )
       End if
       If ( Dirichlet_nu_t == 0 ) f = 1d0
       Delta                   = ( dx*dz*( y(j) - y(j-1) ) )**(1d0/3d0)
       nu_t_(2:nxg-1,j,2:nzg-1) = ( (Cs*Delta*f)**2d0 ) * S(2:nxg-1,j,2:nzg-1)
    End Do
 
    ! wall-parallel average
    Do j = 2, nyg-1
       avg_nu_t_(1,j,1) = Sum( nu_t_(2:nxg-1,j,2:nzg-1) )/Real( (nxg-2)*(nzg-2), 8)
    End Do

    ! clipping negative values
    avg_nu_t_ = Max( avg_nu_t_, 0d0 )

    ! boundary conditions ghost cell (must be done after clipping)
    If ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(avg_nu_t_,2)
    Else
       Call apply_Neumann_bc_y  (avg_nu_t_,2)
    End If

  End Subroutine sgs_Smagorinsky

  !-----------------------------------------------------------!
  !                                                           !
  !           Compute rate-of-strain at cell centers          !
  !                                                           !
  ! Input:  U, V, W, (flow velocities)                        !
  ! Output: Sij, S=sqrt(2Sij*Sij)                             !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine compute_Sij(U_,V_,W_,Sij_,S_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_    

    Real(Int64), Dimension(2:nxg,  2:nyg,  2:nzg,  6), Intent(Out) :: Sij_
    Real(Int64), Dimension(2:nxg-1,2:nyg-1,2:nzg-1),   Intent(Out) :: S_

    ! local variables
    Integer(Int32) :: i, j, k

    ! Compute diagonal Sij
    Do k = 2,nzg-1
      Do j = 2,nyg-1
         Do i = 2,nxg-1
          ! These are at cell center
          Sij_(i,j,k,1) = (U_(i,j,k) - U_(i-1,j,k))/(x(i)-x(i-1)) ! dU/dx
          Sij_(i,j,k,2) = (V_(i,j,k) - V_(i,j-1,k))/(y(j)-y(j-1)) ! dV/dy
          Sij_(i,j,k,3) = (W_(i,j,k) - W_(i,j,k-1))/(z(k)-z(k-1)) ! dW/dz
        End Do
      End Do
    End Do

    ! Compute off-diagonal Sij
    Do k = 2,nzg 
      Do j = 2,nyg
         Do i = 2,nxg
          ! These are at cell edges
          Sij_(i,j,k,4) = 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j-1,k))/(yg(j)-yg(j-1)) + (V_(i,j-1,k)-V_(i-1,j-1,k))/(xg(i)-xg(i-1))) ! 1/2*(dU/dy + dV/dx)
          Sij_(i,j,k,5) = 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j,k-1))/(zg(k)-zg(k-1)) + (W_(i,j,k-1)-W_(i-1,j,k-1))/(xg(i)-xg(i-1))) ! 1/2*(dU/dz + dW/dx)
          Sij_(i,j,k,6) = 0.5d0 * ((V_(i,j-1,k)-V_(i,j-1,k-1))/(zg(k)-zg(k-1)) + (W_(i,j,k-1)-W_(i,j-1,k-1))/(yg(j)-yg(j-1))) ! 1/2*(dV/dz + dW/dy)
        End Do
      End Do
    End Do

    ! Move values from cell edge to cell center
    Do k = 2,nzg-1
      Do j = 2,nyg-1
         Do i = 2,nxg-1       
          Sij_(i,j,k,4) = 0.25d0 * (Sij_(i,j,k,4) + Sij_(i+1,j,k,4) + Sij_(i,j+1,k,4) + Sij_(i+1,j+1,k,4)) 
          Sij_(i,j,k,5) = 0.25d0 * (Sij_(i,j,k,5) + Sij_(i+1,j,k,5) + Sij_(i,j,k+1,5) + Sij_(i+1,j,k+1,5))
          Sij_(i,j,k,6) = 0.25d0 * (Sij_(i,j,k,6) + Sij_(i,j+1,k,6) + Sij_(i,j,k+1,6) + Sij_(i,j+1,k+1,6))
        End Do
      End Do
    End Do

    ! Compute |S| = sqrt(2*Sij*Sij) (at cell centers)
    S_ = 0d0
    Do i = 1,6
      If (i .le. 3) Then
        S_ = S_ + 2d0 * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i) * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i)
      Else 
        S_ = S_ + 4d0 * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i) * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i)
      End If
    End Do
    S_ = S_**0.5d0

  End Subroutine compute_Sij

  !------------------------------------------------------------------------------!
  !                                                                              !
  !                    Compute rate-of-strain at V locations                     !
  !                     for the first two cell at the wall                       !
  !                                                                              !
  ! Assumed uniform mesh in y                                                    !
  !                                                                              !
  ! Input:  U, V, W, (flow velocities)                                           !
  ! Output: Sij(2:nxg,2:3,2:nzg) (bottom) and Sij(2:nxg,nyg-1:nyg,2:nzg) (top)   !
  !                                                                              !
  !------------------------------------------------------------------------------!
  Subroutine compute_Sij_at_V_location_wall(U_,V_,W_,Sij_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_    

    Real(Int64), Dimension(2:nxg,2:nyg,2:nzg,6), Intent(Out) :: Sij_

    ! local variables
    Integer(Int32) :: i, j, k, jv
    Integer(Int32) :: jindex(4), jj

    ! indices for diagonal part
    jindex(1) = 2
    jindex(2) = 3
    jindex(3) = nyg-1
    jindex(4) = nyg

    ! Compute diagonal Sij
    Do k = 2, nzg-1
       Do jj = 1, 4
          j  = jindex(jj)
          Do i = 2, nxg-1
             ! These are at cell center
             ! dU/dx
             Sij_(i,j,k,1) = 0.5d0*( (U_(i,j,k) - U_(i-1,j,k))/(x(i)-x(i-1)) + (U_(i,j-1,k) - U_(i-1,j-1,k))/(x(i)-x(i-1)) )
             ! dV/dy first order
             jv = j
             If ( j>nyg/2 ) jv = j-1
             Sij_(i,j,k,2) = (V_(i,jv,k) - V_(i,jv-1,k))/(y(jv)-y(jv-1)) 
             ! dW/dz
             Sij_(i,j,k,3) = 0.5d0*( (W_(i,j,k) - W_(i,j,k-1))/(z(k)-z(k-1)) + (W_(i,j-1,k) - W_(i,j-1,k-1))/(z(k)-z(k-1)) )
          End Do
       End Do
    End Do

    ! indices for off-diagonal part
    jindex(1) = 2
    jindex(2) = 3
    jindex(3) = nyg-1
    jindex(4) = nyg

    ! Compute off-diagonal Sij
    Do k = 2, nzg-1
       Do jj = 1, 4
          j  = jindex(jj)
          jv = j
          If ( j>nyg/2 ) jv = j-1
          Do i = 2, nxg-1
             ! These are at cell edges
             ! 1/2*(dU/dy + dV/dx)
                                     !(i,j)
             Sij_(i,j,k,4) = 0.5d0*( 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j-1,k))/(yg(j)-yg(j-1)) + (V_(i,jv-1,k)-V_(i-1,jv-1,k))/(xg(i)-xg(i-1))) + & 
                                     !(i+1,j)
                                     0.5d0 * ((U_(i-1+1,j,k)-U_(i-1+1,j-1,k))/(yg(j)-yg(j-1)) + (V_(i+1,jv-1,k)-V_(i-1+1,jv-1,k))/(xg(i+1)-xg(i-1+1))) )   
             ! 1/2*(dU/dz + dW/dx)
             Sij_(i,j,k,5) = 0.25d0*( &
                  !(i,k)
                  0.5d0 * ((U_(i-1,j,k)-U_(i-1,j,k-1))/(zg(k)-zg(k-1)) + (W_(i,j,k-1)-W_(i-1,j,k-1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k)
                  0.5d0 * ((U_(i-1+1,j,k)-U_(i-1+1,j,k-1))/(zg(k)-zg(k-1)) + (W_(i+1,j,k-1)-W_(i-1+1,j,k-1))/(xg(i+1)-xg(i-1+1))) + &  
                  !(i,k+1)
                  0.5d0 * ((U_(i-1,j,k+1)-U_(i-1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) + (W_(i,j,k-1+1)-W_(i-1,j,k-1+1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k+1)
                  0.5d0 * ((U_(i-1+1,j,k+1)-U_(i-1+1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) + (W_(i+1,j,k-1+1)-W_(i-1+1,j,k-1+1))/(xg(i+1)-xg(i-1+1))) )    

             ! 1/2*(dV/dz + dW/dy) 
                                     !(j,k)
             Sij_(i,j,k,6) = 0.5d0*( 0.5d0 * ((V_(i,jv-1,k)-V_(i,jv-1,k-1))/(zg(k)-zg(k-1)) + (W_(i,j,k-1)-W_(i,j-1,k-1))/(yg(j)-yg(j-1))) + & 
                                     !(j,k+1)
                                     0.5d0 * ((V_(i,jv-1,k+1)-V_(i,jv-1,k-1+1))/(zg(k+1)-zg(k-1+1)) + (W_(i,j,k-1+1)-W_(i,j-1,k-1+1))/(yg(j)-yg(j-1))) ) 
          End Do
       End Do
    End Do

    ! periodicity in x and z
    !Sij_(nxg,:,:,:) = Sij_(3,:,:,:)
    !Sij_(:,:,nzg,:) = Sij_(:,:,3,:)

    ten_buf(2:nxg,2:nyg,2:nzg,:) = Sij
    Do i = 1,6
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do
    Sij = ten_buf(2:nxg,2:nyg,2:nzg,:)

    
  End Subroutine compute_Sij_at_V_location_wall
  
End Module subgrid
