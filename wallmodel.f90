!------------------------------------!
!     Module for LES wall-models     !
!------------------------------------!
Module wallmodel

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use interpolation
  Use subgrid
  Use boundary_conditions
  Use Newton_solver
  Use functions_wallmodel

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------!
  !                                              !
  !              Select wall model               !
  !                                              !
  !----------------------------------------------!
  Subroutine compute_wall_model(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Integer(Int32) :: i, j
    Real   (Int64) :: T_alpha

    ! save alpha from time n
    alpha_xo = alpha_x
    alpha_yo = alpha_y
    alpha_zo = alpha_z

    ! select stress model
    If ( istress_model==1 ) Then
       ! law-of-the-wall model          
       Call compute_law_of_the_wall_model(U_)
    Else
       utau_model = utau_ref
    End If

    ! select boundary condition type
    If     ( iwall_model==1 ) Then
       ! momentum wall-model
       Call compute_alpha_momentum_model(U_,V_,W_)
    Elseif ( iwall_model==2 ) Then
       ! Boses's wall model
       Call compute_alpha_bose_model(U_,V_,W_)       
    Elseif ( iwall_model==3 ) Then
       ! Fitting best Robin BC at the wall
       Call compute_alpha_robin_fitting_model(U_,V_,W_)
    Elseif ( iwall_model==4 ) Then
       ! momentum equalibrium at the wall
       Call compute_alpha_momentum_wall_model(U_,V_,W_)
    Elseif ( iwall_model==5 ) Then
       ! set constant alpha_x = alpha_y = alpha_z
       Call compute_constant_alpha
    Elseif ( iwall_model==6 ) Then
       ! log-layer model, actual eddy viscosity
       Call compute_log_layer_wall_model(1)
    Elseif ( iwall_model==7 ) Then
       !  dynamic momentum equalibrium at the wall (for two filter)
       Call compute_alpha_momentum_wall_model_dynamic(U_,V_,W_)
    Elseif ( iwall_model==8 ) Then
       ! log-layer model, constant eddy viscosity
       Call compute_log_layer_wall_model(2)
    Elseif ( iwall_model==9 ) Then
       ! impose wall stress from turbulent Cf (with alpha_v)
       Call compute_alpha_v_from_Cf(U_,V_,W_)
    Elseif ( iwall_model==11 ) Then
       ! momentum equilibrium at the wall with slip_u/=0 and slip_v-w=0
       Call compute_alpha_u_momentum_wall_model(U_,V_,W_)
    Elseif ( iwall_model==12 ) Then
       ! impose wall stress from turbulent Cf (with alpha_u)
       Call compute_alpha_u_from_Cf(U_,V_,W_)
    Elseif ( iwall_model==13 ) Then
       ! impose wall stress from turbulent Cf (with v)
       Call compute_v_from_Cf(U_,V_,W_)
    Elseif ( iwall_model==14 ) Then
       ! 
       Call compute_v_from_Cf(U_,V_,W_)
    Elseif ( iwall_model==15 ) Then
       ! exact Neumann condition
       Call compute_exact_Neumann_from_Cf(U_,V_,W_)
    End If

    ! penetration vs. no penetration at the wall
    ! 1 -> penetration    (V/=0) 
    ! 0 -> no penetration (V =0)
    If ( penetration == 0 ) Then
       alpha_y = 0d0
    End If

    ! alpha limiter
    If ( .False. ) Then
       Do i=1,nx_global
          alpha_x(i,:,:) = Max( alpha_x(i,1,1), -0.1d0 )
          alpha_x(i,:,:) = Min( alpha_x(i,1,1),  0.1d0 )
       End Do       
       Do i=1,nxg_global
          alpha_y(i,:,:) = Max( alpha_y(i,1,1), -0.1d0 )
          alpha_y(i,:,:) = Min( alpha_y(i,1,1),  0.1d0 )
          alpha_z(i,:,:) = Max( alpha_z(i,1,1), -0.1d0 )
          alpha_z(i,:,:) = Min( alpha_z(i,1,1),  0.1d0 )
       End Do
    End If
    ! weighted alpha in time
    If ( .False. ) Then
       T_alpha = 1d0
       alpha_x = dt/T_alpha*alpha_x + (1d0-dt/T_alpha)*alpha_xo
       alpha_y = dt/T_alpha*alpha_y + (1d0-dt/T_alpha)*alpha_yo
       alpha_z = dt/T_alpha*alpha_z + (1d0-dt/T_alpha)*alpha_zo
    End If
    ! smooth in x
    If ( .False. ) Then
       Do j=1,10
          alpha_xo = alpha_x
          alpha_yo = alpha_y
          alpha_zo = alpha_z
          Do i=2,nx-1
             alpha_x(i,:,:) = 0.5d0*alpha_xo(i,:,:) + 0.25d0*alpha_xo(i-1,:,:) + 0.25d0*alpha_xo(i+1,:,:)
          End Do
          Do i=2,nxg-1
             alpha_y(i,:,:) = 0.5d0*alpha_yo(i,:,:) + 0.25d0*alpha_yo(i-1,:,:) + 0.25d0*alpha_yo(i+1,:,:)
          End Do
          Do i=2,nxg-1
             alpha_z(i,:,:) = 0.5d0*alpha_zo(i,:,:) + 0.25d0*alpha_zo(i-1,:,:) + 0.25d0*alpha_zo(i+1,:,:)
          End Do
       End Do
    End If
    
    ! compute boundary conditions for pseudo-pressure
    ! Not used for now
    !If ( iwall_model>0 ) Then
    !   Call compute_pseudo_pressure_bc_for_robin_bc
    !End If
    
  End Subroutine compute_wall_model

  !-----------------------------------------------------------!
  !                  Momentum wall-model                      !
  !     Compute slip lengths alpha_x, alpha_y and alpha_z     !
  !           for Ui = alpha_i dUi/dy at the wall             !
  !                                                           !
  ! Equations:                                                !
  !                                                           !
  !   alpha_x = <int_0^l(  hat(U)hat(V) + T_UV                !
  !                       - UV - tau_UV )dy>/<nu*a*dU/dy>     ! 
  !   alpha_y = alpha_x                                       !
  !   alpha_z = alpha_x                                       !
  !                                                           ! 
  ! Nomeclature:                                              ! 
  !                                                           ! 
  !   tau_UV = -2*nu_t*Suv                                    ! 
  !   T_UV   = -2*nu_t_hat*Suv_hat                            ! 
  !                                                           ! 
  !   <.> -> wall-parallel average                            ! 
  !   hat -> test filter                                      !
  !   l   -> integration length = int_len                     !
  !   a   -> 1-Delta_hat/Delta  = 1-fil_size                  !
  !                                                           !
  ! Note (not used):                                          !
  !                                                           !
  !   hat(U)hat(V) + hat(tau_UV) + Luv = hat(UV+tau_UV)       !
  !                                                           !
  ! Assumed instantaneously:                                  !
  !                                                           !
  !   2(    -UV       - tau_UV + nu*dU/dy     )_wall = dPx*Ly !
  !   2(-hat(U)hat(V) - T_UV   + nu*dhat(U)/dy)_wall = dPx*Ly !
  !                                                           !
  ! alpha_i is defined at cell edges (in xy plane)            !
  !                                                           !
  ! Input:  U_,V_,W_,avg_nu_t (flow velocities)               !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)          !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine compute_alpha_momentum_model(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64) :: dU, Err
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! Part 1: compute UV + tau_UV -> term_4 (at cell edges in xy plane)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate avg_nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(avg_nu_t,term_2(1:nxg,1:ny,1:1),2)

    ! interpolate U and V to cell edges
    Call interpolate_y(U_,term  ,2) 
    Call interpolate_x(V_,term_1,1)

    ! term_2 = UV + tau_UV at cell edges
    ! term   -> U
    ! term_1 -> V
    ! term_2 -> nu_t
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = term(1:nx,j,1:nzg)*term_1(1:nx,j,1:nzg) - 2d0*term_2(1,j,1)*term_3(1:nx,j,1:nzg)
    End Do

    !---------------------------------------------------------------------------!
    ! Part 2: compute hat(U)hat(V) + T_UV -> term_2 (at cell edges in xy plane)

    ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
    Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
    Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
    Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )

    ! apply periodicity in x and z
    Call apply_periodic_bc_x(Uff,1)
    Call apply_periodic_bc_z(Uff,2)
    Call apply_periodic_bc_x(Vff,2)
    Call apply_periodic_bc_z(Vff,2)
    Call apply_periodic_bc_x(Wff,2)
    Call apply_periodic_bc_z(Wff,1)

    ! compute eddy viscosity for filtered velocities
    Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t) ! nu_t should be change!

    ! interpolate avg_nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(avg_nu_t_hat,term_2(1:nxg,1:ny,1:1),2)

    ! compute Suv_hat at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (Uff(i,  j+1,k) - Uff(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (Vff(i+1,j  ,k) - Vff(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do

    ! interpolate Uff and Vff to cell edges
    Call interpolate_y(Uff,term  ,2) 
    Call interpolate_x(Vff,term_1,1)

    ! term_2 = hat(U)hat(V) + T_UV at cell edges
    ! term   -> Uff
    ! term_1 -> Vff
    ! term_2 -> nu_t_hat
    ! term_3 -> Suv_hat
    Do j = 1, ny
       term_2(1:nx,j,1:nzg) = term(1:nx,j,1:nzg)*term_1(1:nx,j,1:nzg) - 2d0*term_2(1,j,1)*term_3(1:nx,j,1:nzg)
    End Do
    
    !------------------------------------------------------------------------!
    ! Part 3: compute integral
    ! term_4 -> UV + tau_UV
    ! term_2 -> hat(U)hat(V) + T_UV
    ! alpha_x(:,1,:) -> integral for bottom wall (at cell edge in xy planes)
    ! alpha_x(:,2,:) -> integral for top    wall (at cell edge in xy planes)

    ! compute hat(U)hat(V) + T_UV - ( UV + tau_UV )
    term_2 = term_2 - term_4 
    
    ! integration from bottom wall with trapezoidal rule
    alpha_x(:,1,:) = 0d0
    Do j = 2, ny
       If ( (y(j)-y(1)) <= int_len ) Then
          alpha_x(:,1,2:nzg-1) = alpha_x(:,1,2:nzg-1) + 0.5d0*(term_2(1:nx,j,2:nzg-1) + term_2(1:nx,j-1,2:nzg-1))*(y(j)-y(j-1))
       End If
    End Do
    
    ! integration from top wall with trapezoidal rule (sign changed to later average alpha)
    alpha_x(:,2,:) = 0d0
    Do j = 1, ny-1
      If ( (y(ny)-y(j)) <= int_len ) Then
        alpha_x(:,2,2:nzg-1) = alpha_x(:,2,2:nzg-1) - 0.5d0*(term_2(1:nx,j,2:nzg-1) + term_2(1:nx,j+1,2:nzg-1))*(y(j+1)-y(j))
      End If
    End Do 

    ! averaging in xz plane
    alpha_x(:,:,:) = Sum( alpha_x(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 4: compute <dU/dy> at the wall
 
    ! compute term = dU/dy at each wall
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =   (U_(i,  2,k) - U_(i,    1,k)) / (yg(  2) - yg(    1))
          ! top wall (sign changed)
          term(i,2,k) = - (U_(i,nyg,k) - U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average dU/dy
    dU = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 5: compute alpha_i

    ! alpha_x = Integral_0^l( hat(U)hat(V)+T_UV - UV - tau_UV )/<nu*a*dU/dy>
    alpha_x(:,:,:) = alpha_x(:,:,:) / ( nu*( 1d0-fil_size )*dU )

    ! clipping 
    alpha_x = Max( alpha_x, 0d0 )

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)

  End Subroutine compute_alpha_momentum_model

  !--------------------------------------------------------------------!
  !                                                                    !
  !                        Boses' wall-model                           !
  !                                                                    !
  !          Compute slip lengths alpha_x, alpha_y and alpha_z         !
  !                  for ui = alpha_i dui/dy at the wall               !
  !                                                                    !
  ! Equation:                                                          !
  !                                                                    !
  !  alpha^2*( fil_size^2*dhat(ui)/dn*dhat(uj)/dn - dui/dn*duj/dn ) +  !
  !            Tij - hat(tau_ij) = hat(ui*uj) - ui*uj                  !
  !                                                                    !
  !  n->normal to the wall                                             !
  !                                                                    !
  ! Definitions:                                                       !
  !                                                                    !
  !  Lij    = hat(ui*uj) - ui*uj - Tij + hat(tau_ij)                   !
  !  Mij    = fil_size^2*hat(dui/dn)*hat(duj/dn) - dui/dn*duj/dn       !
  !  tau_ij = -2*nu_t*Sij                                              !
  !  Tij    = -2*nu_t_hat*hat(Sij)                                     !
  !                                                                    !
  ! alpha = alpha_x = alpha_y = alpha_z                                !
  ! alpha computed at V positions                                      !
  !                                                                    !
  ! Input:  U_,V_,W_,avg_nu_t (flow velocities)                        !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)                   !
  !                                                                    !
  !--------------------------------------------------------------------!
  Subroutine compute_alpha_bose_model(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Integer(Int32) :: i, j, k

    !--------------------------------------------------------------------------!
    ! Part 1: Compute Lij = hat(ui*uj) - ui*uj
    
    ! interpolate velocity to U and W to V location (at wall)
    call interpolate_x(U_,     term_1(2:nxg-1,:,:),1) 
    call interpolate_y(term_1, term  (:,1:ny,:),2)
    
    call interpolate_z(W_,     term_1(:,:,2:nz),1) 
    call interpolate_y(term_1, term_2(:,1:ny,:),2) 

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! Sij = ui*uj at v location (at wall)
    ! U at V location -> term  (2:nx-1,1:ny,:)
    ! W at V location -> term_2(:,1:ny,2:nzg-1)
    ! bottom wall 
    Sij(:,2:3,:,1) = term  (2:nxg,1:2,2:nzg) * term  (2:nxg,1:2,2:nzg)  ! u^2
    Sij(:,2:3,:,2) = V_    (2:nxg,1:2,2:nzg) * V_    (2:nxg,1:2,2:nzg)  ! v^2
    Sij(:,2:3,:,3) = term_2(2:nxg,1:2,2:nzg) * term_2(2:nxg,1:2,2:nzg)  ! w^2
    Sij(:,2:3,:,4) = term  (2:nxg,1:2,2:nzg) * V_    (2:nxg,1:2,2:nzg)  ! uv
    Sij(:,2:3,:,5) = term  (2:nxg,1:2,2:nzg) * term_2(2:nxg,1:2,2:nzg)  ! uw
    Sij(:,2:3,:,6) = V_    (2:nxg,1:2,2:nzg) * term_2(2:nxg,1:2,2:nzg)  ! vw
    ! top wall
    Sij(:,nyg-2:nyg-1,:,1) = term  (2:nxg,ny-1:ny,2:nzg) * term  (2:nxg,ny-1:ny,2:nzg) ! u^2
    Sij(:,nyg-2:nyg-1,:,2) = V_    (2:nxg,ny-1:ny,2:nzg) * V_    (2:nxg,ny-1:ny,2:nzg) ! v^2
    Sij(:,nyg-2:nyg-1,:,3) = term_2(2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! w^2
    Sij(:,nyg-2:nyg-1,:,4) = term  (2:nxg,ny-1:ny,2:nzg) * V_    (2:nxg,ny-1:ny,2:nzg) ! uv
    Sij(:,nyg-2:nyg-1,:,5) = term  (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! uw
    Sij(:,nyg-2:nyg-1,:,6) = V_    (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! vw

    ten_buf(2:nxg,2:nyg,2:nzg,1:6) = Sij;
    Do i = 1,6 
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do

    ! Compute Lij = hat(ui*uj) -> stored in Lij(3:nxg-1,2:nyg-1,3:nzg-1,:)
    Call filter_tensor_xzy(ten_buf(1:nxg,        2:5,1:nzg,:),Lij(2:nxg-1,        2:5,2:nzg-1,:)) 
    Call filter_tensor_xzy(ten_buf(1:nxg,nyg-4:nyg-1,1:nzg,:),Lij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)) 
    !Call filter_tensor_xzy(Sij(2:nxg,        2:5,2:nzg,:),Lij(3:nxg-1,        2:5,3:nzg-1,:)) 
    !Call filter_tensor_xzy(Sij(2:nxg,nyg-4:nyg-1,2:nzg,:),Lij(3:nxg-1,nyg-4:nyg-1,3:nzg-1,:)) 

    ! fill in missing values (periodicity) 
    !Lij(2,:,:,:) = Lij(nxg-1,:,:,:)
    !Lij(:,:,2,:) = Lij(:,:,nzg-1,:)

    ! Compute Lij = hat(ui*uj) - ui*uj at the wall 
    ! Valid for Lij(2:nxg-1,(/2,nyg-1/),2:nzg-1)
    Lij(2:nxg-1,    2,2:nzg-1,:) = Lij(2:nxg-1,    2,2:nzg-1,:) - Sij(2:nxg-1,    2,2:nzg-1,:)
    Lij(2:nxg-1,nyg-1,2:nzg-1,:) = Lij(2:nxg-1,nyg-1,2:nzg-1,:) - Sij(2:nxg-1,nyg-1,2:nzg-1,:)

    If ( Dirichlet_nu_t==0 ) Then
    !-----------------------------------------------------------------------------!
    ! Part 2: Compute -Tij + hat(tau_ij) at the wall and add it to Lij

       !----------compute hat(tau_ij)       
       ! compute Sij at V locations in the first two cells-> Sij(2:nxg,2:3,2:nzg) (bottom) and Sij(2:nxg,nyg-1:nyg,2:nzg) (top)
       Call compute_Sij_at_V_location_wall(U_,V_,W_,Sij)
       
       ! interpolate avg_nu_t from cell centers to cell faces (== cell edges because averaged in xz)
       Call interpolate_y(avg_nu_t,term_2(1:nxg,1:ny,1:1),2)
       
       ! tau_ij = -2*nu_t*Sij -> Sij
       Do i = 2, nxg
          Do k = 2, nzg
             Do j = 1, 6
                Sij(i,      2:3,k,j) = -2d0*term_2(1,    1:2,1)*Sij(i,      2:3,k,j) ! bottom 
                Sij(i,nyg-1:nyg,k,j) = -2d0*term_2(1,ny-1:ny,1)*Sij(i,nyg-1:nyg,k,j) ! top
             End Do
          End Do
       End Do
       

       ten_buf(2:nxg,2:nyg,2:nzg,1:6) = Sij;
       Do i = 1,6
          Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
          Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
          Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
       End Do
       ! hat(tau_ij) -> filter tau_ij 
       !     tau_ij  -> Sij(:,2:4,:) (bottom) and Sij(:,nyg-2:nyg  ,:)
       ! hat(tau_ij) -> Mij(:,2:4,:) (bottom) and Mij(:,nyg-3:nyg-1,:)
       Call filter_tensor_xzy(ten_buf(1:nxg,      2:4,1:nzg,:),Mij(2:nxg-1,        2:4,2:nzg-1,:)) 
       Call filter_tensor_xzy(ten_buf(1:nxg,nyg-2:nyg,1:nzg,:),Mij(2:nxg-1,nyg-3:nyg-1,2:nzg-1,:)) 
       !Call filter_tensor_xzy(Sij(2:nxg,      2:4,2:nzg,:),Mij(3:nxg-1,        2:4,3:nzg-1,:)) 
       !Call filter_tensor_xzy(Sij(2:nxg,nyg-2:nyg,2:nzg,:),Mij(3:nxg-1,nyg-3:nyg-1,3:nzg-1,:)) 
       
       ! periodicity
       !Mij(2,:,:,:) = Mij(nxg-1,:,:,:)
       !Mij(:,:,2,:) = Mij(:,:,nzg-1,:)    
       
       ! add it
       Lij(2:nxg-1,    2,2:nzg-1,:) = Lij(2:nxg-1,    2,2:nzg-1,:) + Mij(2:nxg-1,    2,2:nzg-1,:)
       Lij(2:nxg-1,nyg-1,2:nzg-1,:) = Lij(2:nxg-1,nyg-1,2:nzg-1,:) + Mij(2:nxg-1,nyg-1,2:nzg-1,:)
       
       !----------compute T_ij
       ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
       Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
       Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
       Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
       
       ! apply periodicity in x and z
       Call apply_periodic_bc_x(Uff,1)
       Call apply_periodic_bc_z(Uff,1)
       Call update_ghost_interior_planes(Uff,1)
       Call apply_periodic_bc_x(Vff,2)
       Call apply_periodic_bc_z(Vff,2)
       Call update_ghost_interior_planes(Vff,2)
       Call apply_periodic_bc_x(Wff,2)
       Call apply_periodic_bc_z(Wff,3)
       Call update_ghost_interior_planes(Wff,3)
       
       ! compute eddy viscosity for filtered velocities
       Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t) ! nu_t should be change
       
       ! interpolate avg_nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
       Call interpolate_y(avg_nu_t_hat,term_2(1:nxg,1:ny,1:1),2)    
       
       ! compute hat(Sij) at V locations in two cells-> Sij(2:nxg,2:3,2:nzg) (bottom) and Sij(2:nxg,nyg-1:nyg,2:nzg) (top)
       Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,Sij)
       
       ! Tij = -2*nu_t_hat*hat(Sij) -> Sij
       Do i = 2, nxg
          Do k = 2, nzg
             Do j = 1, 6
                Sij(i,  2,k,j) = -2d0*term_2(1, 1,1)*Sij(i,  2,k,j) ! bottom 
                Sij(i,nyg,k,j) = -2d0*term_2(1,ny,1)*Sij(i,nyg,k,j) ! top
             End Do
          End Do
       End Do
       
       ! add it
       Lij(2:nxg-1,    2,2:nzg-1,:) = Lij(2:nxg-1,    2,2:nzg-1,:) - Sij(2:nxg-1,  2,2:nzg-1,:)
       Lij(2:nxg-1,nyg-1,2:nzg-1,:) = Lij(2:nxg-1,nyg-1,2:nzg-1,:) - Sij(2:nxg-1,nyg,2:nzg-1,:)
    
    End If
   
    !-----------------------------------------------------------------------------!
    ! Part 3: Compute Mij = fil_size^2*hat(dui/dy)*hat(duj/dy) - dui/dy*duj/dy

    ! interpolate velocity to cell centers
    call interpolate_x(U_,term  (2:nxg-1,:,:),1) 
    call interpolate_z(W_,term_2(:,:,2:nz),1) 

    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! Compute dui/dn at V location (at wall)
    ! Sij(:,:,:,1) -> dU/dn
    ! Sij(:,:,:,2) -> dV/dn
    ! Sij(:,:,:,3) -> dW/dn 
    Do i = 2, nxg
      Do k = 2, nzg

        Sij(i,2,k,1)     = (term  (i,2,k)     - term  (i,1,k)  ) / (yg(2)-yg(1))            ! dU/dy at lower wall
        Sij(i,2,k,2)     = (V_    (i,2,k)     - V_    (i,1,k)  ) / (y (2)-y (1))            ! dV/dy at lower wall (1st order)
        Sij(i,2,k,3)     = (term_2(i,2,k)     - term_2(i,1,k)  ) / (yg(2)-yg(1))            ! dW/dy at lower wall

        Sij(i,nyg-1,k,1) = (term  (i,nyg-1,k) - term  (i,nyg,k)) / (yg(nyg)-yg(nyg-1))      ! dU/dy at upper wall
        Sij(i,nyg-1,k,2) = (V_    (i,ny -1,k) - V_    (i,ny ,k)) / (y ( ny)-y ( ny-1))      ! dV/dy at upper wall (1st order)
        Sij(i,nyg-1,k,3) = (term_2(i,nyg-1,k) - term_2(i,nyg,k)) / (yg(nyg)-yg(nyg-1))      ! dW/dy at upper wall

        Sij(i,3,k,1)     = (term  (i,3,k)     - term  (i,2,k)  ) / (yg(3)-yg(2))            ! dU/dy at lower wall+1
        Sij(i,3,k,2)     = (V_    (i,3,k)     - V_    (i,2,k)  ) / (y (3)-y (2))            ! dV/dy at lower wall+1 (1st order)
        Sij(i,3,k,3)     = (term_2(i,3,k)     - term_2(i,2,k)  ) / (yg(3)-yg(2))            ! dW/dy at lower wall+1

        Sij(i,nyg-2,k,1) = (term  (i,nyg-2,k) - term  (i,nyg-1,k)) / (yg(nyg-1)-yg(nyg-2))  ! dU/dy at upper wall-1
        Sij(i,nyg-2,k,2) = (V_    (i, ny-2,k) - V_    (i, ny-1,k)) / (y ( ny-1) -y ( ny-2)) ! dV/dy at upper wall-1 (1st order)
        Sij(i,nyg-2,k,3) = (term_2(i,nyg-2,k) - term_2(i,nyg-1,k)) / (yg(nyg-1)-yg(nyg-2))  ! dW/dy at upper wall-1

      End Do
    End Do

    ! hat(dui/dn) (NOTE: vs. dhat(ui)/dn)
    ! Sij(:,:,:,1) -> dU/dn
    ! Sij(:,:,:,2) -> dV/dn
    ! Sij(:,:,:,3) -> dW/dn 
    ! Sij(:,:,:,4) -> hat(dU/dn)
    ! Sij(:,:,:,5) -> hat(dV/dn)
    ! Sij(:,:,:,6) -> hat(dW/dn)
    ten_buf(2:nxg,2:nyg,2:nzg,1:6) = Sij;
    Do i = 1,6
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do

    Call filter_tensor_xzy(ten_buf(1:nxg,        2:5,1:nzg,1:3),Sij(2:nxg-1,        2:5,2:nzg-1,4:6))
    Call filter_tensor_xzy(ten_buf(1:nxg,nyg-4:nyg-1,1:nzg,1:3),Sij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,4:6))
    !Call filter_tensor_xzy(Sij(2:nxg,        2:5,2:nzg,1:3),Sij(3:nxg-1,        2:5,3:nzg-1,4:6))
    !Call filter_tensor_xzy(Sij(2:nxg,nyg-4:nyg-1,2:nzg,1:3),Sij(3:nxg-1,nyg-4:nyg-1,3:nzg-1,4:6))

    ! periodicity
    !Sij(2,:,:,:) = Sij(nxg-1,:,:,:)
    !Sij(:,:,2,:) = Sij(:,:,nzg-1,:)

    ! Mij = fil_size^2*hat(dui/dn)*hat(duj/dn) - dui/dn*duj/dn at wall
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! Diagonal elements
          Do j = 1, 3
             Mij(i,    2,k,j) = fil_size**2d0*Sij(i,    2,k,j+3)**2d0 - Sij(i,    2,k,j)**2d0 ! bottom
             Mij(i,nyg-1,k,j) = fil_size**2d0*Sij(i,nyg-1,k,j+3)**2d0 - Sij(i,nyg-1,k,j)**2d0 ! top
          End Do
          ! Off-diagonal elements
          Mij(i,    2,k,4) = fil_size**2d0*Sij(i,    2,k,4)*Sij(i,    2,k,5) - Sij(i,    2,k,1)*Sij(i,    2,k,2)
          Mij(i,nyg-1,k,4) = fil_size**2d0*Sij(i,nyg-1,k,4)*Sij(i,nyg-1,k,5) - Sij(i,nyg-1,k,1)*Sij(i,nyg-1,k,2)

          Mij(i,    2,k,5) = fil_size**2d0*Sij(i,    2,k,4)*Sij(i,    2,k,6) - Sij(i,    2,k,1)*Sij(i,    2,k,3)
          Mij(i,nyg-1,k,5) = fil_size**2d0*Sij(i,nyg-1,k,4)*Sij(i,nyg-1,k,6) - Sij(i,nyg-1,k,1)*Sij(i,nyg-1,k,3)

          Mij(i,    2,k,6) = fil_size**2d0*Sij(i,    2,k,5)*Sij(i,    2,k,6) - Sij(i,    2,k,2)*Sij(i,    2,k,3)
          Mij(i,nyg-1,k,6) = fil_size**2d0*Sij(i,nyg-1,k,5)*Sij(i,nyg-1,k,6) - Sij(i,nyg-1,k,2)*Sij(i,nyg-1,k,3)
       End Do
    End Do
    
    !-----------------------------------------------------------------------------!
    ! Part 4: Least squares solution term = Lij*Mij, term_1 = Mij*Mij

    term   = 0d0
    term_1 = 0d0
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! Diagonal elements
          Do j = 1, 3
             ! Lij*Mij
             term  (i, 1,k) = term  (i, 1,k) + Lij(i,    2,k,j)*Mij(i,    2,k,j) ! bottom
             term  (i,ny,k) = term  (i,ny,k) + Lij(i,nyg-1,k,j)*Mij(i,nyg-1,k,j) ! top
             ! Mij*Mij
             term_1(i, 1,k) = term_1(i, 1,k) + Mij(i,    2,k,j)*Mij(i,    2,k,j) ! bottom
             term_1(i,ny,k) = term_1(i,ny,k) + Mij(i,nyg-1,k,j)*Mij(i,nyg-1,k,j) ! top
          End Do
          ! Off-diagonal
          Do j = 4, 6
             ! Lij*Mij
             term  (i, 1,k) = term  (i, 1,k) + 2d0*Lij(i,    2,k,j)*Mij(i,    2,k,j) ! bottom
             term  (i,ny,k) = term  (i,ny,k) + 2d0*Lij(i,nyg-1,k,j)*Mij(i,nyg-1,k,j) ! top
             ! Mij*Mij
             term_1(i, 1,k) = term_1(i, 1,k) + 2d0*Mij(i,    2,k,j)*Mij(i,    2,k,j) ! bottom
             term_1(i,ny,k) = term_1(i,ny,k) + 2d0*Mij(i,nyg-1,k,j)*Mij(i,nyg-1,k,j) ! top
          End Do
       End Do
    End Do
    
    !-----------------------------------------------------------------------------!
    ! Part 5: Compute alpha
 
    alpha_x = Sum( term  (2:nxg-1,1,2:nzg-1) + term  (2:nxg-1,ny,2:nzg-1) )
    alpha_y = Sum( term_1(2:nxg-1,1,2:nzg-1) + term_1(2:nxg-1,ny,2:nzg-1) ) 
    
    If (myid == 0) Then
       Call MPI_Reduce(MPI_IN_PLACE,alpha_x(1,1,1),1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(MPI_IN_PLACE,alpha_y(1,1,1),1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    Else
       Call MPI_Reduce(alpha_x(1,1,1),0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(alpha_y(1,1,1),0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    End If

    Call Mpi_bcast (  alpha_x(1,1,1),1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  alpha_y(1,1,1),1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    ! Average over wall direction. Alpha^2 = <Lij*Mij> / <Mij*Mij>
    !alpha_x = Sum( term  (2:nxg-2,1,2:nzg-2) + term  (2:nxg-2,ny,2:nzg-2) ) / & 
    !          Sum( term_1(2:nxg-2,1,2:nzg-2) + term_1(2:nxg-2,ny,2:nzg-2) ) 
    alpha_x = alpha_x(1,1,1) / alpha_y(1,1,1)  
 
    ! Clipping if necessary
    alpha_x = Max( alpha_x,0d0 )
    alpha_x = alpha_x**0.5d0
    
    ! alpha_y
    alpha_y = alpha_x(2,1,2)
    
    ! alpha_z
    alpha_z = alpha_x(2,1,2) 
    
  End Subroutine compute_alpha_bose_model

  !--------------------------------------------------------------!
  !                                                              !
  !              Momentum balance at the wall                    !
  !                                                              !
  !     Compute slip lengths alpha_x, alpha_y and alpha_z        !
  !           for ui = alpha_i dui/dy at the wall                !
  !                                                              !
  ! alpha_x = alpha_y = alpha_z                                  !
  !                                                              !
  ! From x-momentum balance at the walls:                        !
  !       -2<UV+tau_UV> = dPdx*Ly - 2*nu*<dU/dy>                 !
  !    -> -2<UV> = dPdx*Ly - 2*<(nu+nu_t)*dU/dy>                 !
  !                                                              !
  ! Robin condition at the wall:                                 !
  !        U = alpha dU/dy, V = alpha dV/dy                      !
  !                                                              !
  ! Condition (everything is at the wall):                       !
  !      alpha^2 = (<nu*dU/dy>-dPdx*Ly/2-<tau_UV>)/<dU/dy*dV/dy> !
  !                                                              !
  ! Input:  U_,V_,W_,avg_nu_t (flow velocities)                  !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_alpha_momentum_wall_model(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64) :: dU, dUV, Err, mtau
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at cell edges in xy plane)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate avg_nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(avg_nu_t,term_2(1:nxg,1:ny,1:1),2)

    ! Compute term_4 = tau_UV at cell edges
    ! term_2 -> nu_t
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1,j,1)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    mtau = Sum( term_4(2:nx-1,1,2:nzg-1) - term_4(2:nx-1,ny,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall    

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k))/(yg(nyg) - yg(nyg-1))*(term_1(i,ny,k) - term_1(i,ny-1,k))/(y(ny) - y(ny-1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    dUV = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall    

    ! product
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average <dU/dy>
    dU = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 5: compute alpha_i

    alpha_x(:,:,:) = (nu*dU - dPdx*Ly/2d0 - mtau)/dUV

    ! clipping 
    alpha_x = Max( alpha_x(2,1,2), 0d0 )
    alpha_x = ( alpha_x(2,1,2) )**0.5d0

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)

  End Subroutine compute_alpha_momentum_wall_model

  !--------------------------------------------------------------!
  !                                                              !
  !       Alpha based on least-square fitting of Robin BC        !
  !                                                              !
  !      Compute slip lengths alpha_x, alpha_y and alpha_z       !
  !             for ui = alpha_i dui/dy at the wall              !
  !                                                              !
  ! F(C) = sum < (ui     - C*Delta    *dui/dy    )^2 +           !
  !              (ui_hat - C*Delta_hat*dui_hat/dy)^2 >           !
  ! i    = 1..3 -> U,V,W                                         !
  ! hat -> test filter                                           !
  !                                                              !
  ! Condition:                                                   !
  !    - dF(C)/dC = 0                                            !
  !                                                              !
  ! Equation:                                                    !
  !    - alpha=C*Delta= sum < ui*dui/dy + a*ui_hat*dui_hat/dy >/ !
  !                     sum < dui/dy^2  + a^2*dui_hat/dy^2    >  !
  !      a=Delta_hat/Delta=fil_size                              !
  !                                                              !
  ! alpha = alpha_x = alpha_y = alpha_z                          !
  !                                                              !
  ! Input:  U_,V_,W_ (flow velocities)                           !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine  compute_alpha_robin_fitting_model(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64) :: dU, dUff, Err
    Real   (Int64) :: num_U, num_V,num_W
    Real   (Int64) :: dem_U, dem_V,dem_W
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------------------------!
    ! Part 1: compute filtered velocities

    ! filter velocities
    Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
    Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
    Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )

    ! apply periodicity in x and z
    Call apply_periodic_bc_x(Uff,1)
    Call apply_periodic_bc_z(Uff,2)
    Call apply_periodic_bc_x(Vff,2)
    Call apply_periodic_bc_z(Vff,2)
    Call apply_periodic_bc_x(Wff,2)
    Call apply_periodic_bc_z(Wff,1)
    
    !---------------------------------------------------------!
    ! Part 2: compute numerator < ui*dui/dy + fil_size*ui_hat dui_hat/dy > at the wall
 
    ! interpolate U and Uff in y (centers to faces)
    Call interpolate_y(U_ ,term_1,2) 
    Call interpolate_y(Uff,term_2,2) 
   
    ! U
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          dU          = (   U_(i,2,k) -   U_(i,1,k)) / (yg(2) - yg(1))
          dUff        = (  Uff(i,2,k) -  Uff(i,1,k)) / (yg(2) - yg(1))
          term(i,1,k) = term_1(i,1,k)*dU + fil_size*term_2(i,1,k)*dUff
          ! top wall (sign changed to average)
          dU          = -(   U_(i,nyg,k) -   U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          dUff        = -(  Uff(i,nyg,k) -  Uff(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          term(i,2,k) = term_1(i,ny,k)*dU + fil_size*term_2(i,ny,k)*dUff
       End Do
    End Do
    num_U = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    ! No interpolattion for V (->first order)
    ! V
    term(:,1:2,:) = 0d0
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! bottom wall
          dU          = (   V_(i,2,k) -   V_(i,1,k)) / (y(2) - y(1))
          dUff        = (  Vff(i,2,k) -  Vff(i,1,k)) / (y(2) - y(1))
          term(i,1,k) = V(i,1,k)*dU + fil_size*Vff(i,1,k)*dUff
          ! top wall (sign changed to average)
          dU          = (   V_(i,ny,k) -   V_(i,ny-1,k)) / (y(ny) - y(ny-1))
          dUff        = (  Vff(i,ny,k) -  Vff(i,ny-1,k)) / (y(ny) - y(ny-1))
          term(i,2,k) = -V(i,ny,k)*dU - fil_size*Vff(i,ny,k)*dUff
       End Do
    End Do
    num_V = Sum( term(2:nxg-1,1:2,2:nzg-1) ) / Real( 2*(nxg-2)*(nzg-2), 8)

    ! interpolate W and Wff in y (centers to faces)
    Call interpolate_y(W_ ,term_1,2) 
    Call interpolate_y(Wff,term_2,2) 
   
    ! W
    term(:,1:2,:) = 0d0
    Do i = 2, nxg-1
       Do k = 2, nz-1
          ! bottom wall
          dU          = (   W_(i,2,k) -   W_(i,1,k)) / (yg(2) - yg(1))
          dUff        = (  Wff(i,2,k) -  Wff(i,1,k)) / (yg(2) - yg(1))
          term(i,1,k) = term_1(i,1,k)*dU + fil_size*term_2(i,1,k)*dUff
          ! top wall (sign changed to average)
          dU          = -(   W_(i,nyg,k) -   W_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          dUff        = -(  Wff(i,nyg,k) -  Wff(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          term(i,2,k) = term_1(i,ny,k)*dU + fil_size*term_2(i,ny,k)*dUff
       End Do
    End Do
    num_W = Sum( term(2:nxg-1,1:2,2:nz-1) ) / Real( 2*(nxg-2)*(nz-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute denomerator < dui/dy^2 + fil_size^2*ui_hat dui_hat/dy^2 > at the wall
    
    ! U
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          dU          = (   U_(i,2,k) -   U_(i,1,k)) / (yg(2) - yg(1))
          dUff        = (  Uff(i,2,k) -  Uff(i,1,k)) / (yg(2) - yg(1))
          term(i,1,k) = dU**2d0 + fil_size**2d0*dUff**2d0
          ! top wall (no need to change sign)
          dU          = (   U_(i,nyg,k) -   U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          dUff        = (  Uff(i,nyg,k) -  Uff(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          term(i,2,k) = dU**2d0 + fil_size**2d0*dUff**2d0
       End Do
    End Do
    dem_U = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    ! No interpolattion for V (->first order)
    ! V
    term(:,1:2,:) = 0d0
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! bottom wall
          dU          = (   V_(i,2,k) -   V_(i,1,k)) / (y(2) - y(1))
          dUff        = (  Vff(i,2,k) -  Vff(i,1,k)) / (y(2) - y(1))
          term(i,1,k) = dU**2d0 + fil_size**2d0*dUff**2d0
          ! top wall (no need to change sign)
          dU          = (   V_(i,ny,k) -   V_(i,ny-1,k)) / (y(ny) - y(ny-1))
          dUff        = (  Vff(i,ny,k) -  Vff(i,ny-1,k)) / (y(ny) - y(ny-1))
          term(i,2,k) = dU**2d0 + fil_size**2d0*dUff**2d0
       End Do
    End Do
    dem_V = Sum( term(2:nxg-1,1:2,2:nzg-1) ) / Real( 2*(nxg-2)*(nzg-2), 8)
   
    ! W
    term(:,1:2,:) = 0d0
    Do i = 2, nxg-1
       Do k = 2, nz-1
          ! bottom wall
          dU          = (   W_(i,2,k) -   W_(i,1,k)) / (yg(2) - yg(1))
          dUff        = (  Wff(i,2,k) -  Wff(i,1,k)) / (yg(2) - yg(1))
          term(i,1,k) = dU**2d0 + fil_size**2d0*dUff**2d0
          ! top wall (no need to change sign)
          dU          = (   W_(i,nyg,k) -   W_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          dUff        = (  Wff(i,nyg,k) -  Wff(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
          term(i,2,k) = dU**2d0 + fil_size**2d0*dUff**2d0
       End Do
    End Do
    dem_W = Sum( term(2:nxg-1,1:2,2:nz-1) ) / Real( 2*(nxg-2)*(nz-2), 8)

    !---------------------------------------------------------!
    ! Part 5: compute alpha_i

    ! alpha_x
    alpha_x(:,:,:) = (num_U+num_V+num_W)/(dem_U+dem_V+dem_W)

    ! clipping 
    alpha_x = Max( alpha_x, 0d0 )

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)

  End Subroutine compute_alpha_robin_fitting_model

  !---------------------------------------------------!
  !       Set constant alpha_x=alpha_y=alpha_z        !
  !---------------------------------------------------!
  Subroutine compute_constant_alpha

    Real(Int64) :: omega_alpha
    
    omega_alpha = 2d0*pi*freq_mult/(y(2)-y(1))*dPdx**0.5d0
    alpha_x = alpha_mean_x*( 1d0 + alpha_std*dsin(omega_alpha*t) ) 
    alpha_y = alpha_mean_y*( 1d0 + alpha_std*dsin(omega_alpha*t) ) 
    alpha_z = alpha_mean_z*( 1d0 + alpha_std*dsin(omega_alpha*t) ) 

  end Subroutine compute_constant_alpha

  !--------------------------------------------------------------------------!
  !                         Log-layer wall model                             !
  !                                                                          !
  ! Assumptions:                                                             !
  !                                                                          !
  !    Momentum balance at wall: nu_t*dU/dy = utau^2                         !
  !    Log-layer at at y(2)    : U = utau/kappa*ln(y(2)^+) + B*utau          !
  !                                                                          !
  ! Model:                                                                   !
  !                                                                          !
  ! U     = alpha*dU/dy                                                      !
  ! alpha = (nu_t+nu)/(kappa*utau)*ln(Delta*utau/nu*e^{B*kappa})             !
  !                                                                          !
  !                                                                          !
  ! utau  = sqrt(dPx)                                                        !
  ! kappa = 0.41                                                             !
  ! B     = 5.2                                                              !
  ! Delta = y(2)-y(1)                                                        !
  !                                                                          !
  !                                                                          !
  ! Input:  avg_nu_t (eddy viscosity)                                        !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)                         !
  !                                                                          !
  !--------------------------------------------------------------------------!
  Subroutine compute_log_layer_wall_model(iflag_log)

    Integer(Int32), Intent(In) :: iflag_log    

    ! local variables
    Real(Int64) :: Delta, kappa, B, e, utau_

    ! parameters
    kappa = 0.38d0
    B     = 5.20d0
    Delta = y(2)-y(1)    
    e     = dexp(1d0)
    utau_ = dPdx**0.5d0

    If ( iflag_log==1 ) Then
       ! Use actual eddy viscosity
    
       ! eddy viscosity at the wall (center to faces)
       Call interpolate_y(avg_nu_t,term_2(1:nxg,1:ny,1:1),2)
       
       ! compute alpha
       alpha_x = (term_2(1,1,1)+nu)/(kappa*utau_)* & 
                 dlog( Delta*utau_/nu*e**(B*kappa) )
       
    Else
       ! Use fixed eddy viscosity, nu_t = utau*kappa*y (y=Delta)
       
       ! compute alpha
       alpha_x = Delta*dlog( Delta*utau_/nu*e**(B*kappa) )

       ! change eddy viscosity boundary conditions to nu_t = kappa*utau*y
       avg_nu_t(:,  1,:) = -avg_nu_t(:,    2,:) + 2d0*kappa*utau_*Delta
       avg_nu_t(:,nyg,:) = -avg_nu_t(:,nyg-1,:) + 2d0*kappa*utau_*Delta

    end If

    alpha_y = alpha_x
    alpha_z = alpha_x

  End Subroutine compute_log_layer_wall_model

  !--------------------------------------------------------------!
  !                                                              !
  !            Dynamic Momentum balance at the wall              !
  !                                                              !
  !     Compute slip lengths alpha_x, alpha_y and alpha_z        !
  !           for ui = alpha_i dui/dy at the wall                !
  !                                                              !
  ! alpha_x = alpha_y = alpha_z                                  !
  !                                                              !
  ! From x-momentum balance at the walls:                        !
  !   impl: -2<UV+tau_UV> = dPdx*Ly - 2*nu*<dU/dy>               !
  !   test: -2<hat(U)hat(V)+T_UV> = dPdx*Ly - 2*nu*<dhat(U)/dy>  !
  !                                                              !
  ! Robin condition at the wall:                                 !
  !        U = alpha dU/dy, V = alpha dV/dy                      !
  !        hat(U) = alpha dhat(U)/dy, hat(V) = alpha dhat(V)/dy  !
  !                                                              !
  ! Condition (everything is at the wall):                       !
  !    alpha^2 = < tau_UV - T_UV + nu*(dhat(U)dy  - dUdy) > /    !
  !              (fil_size^2<dhat(U)dy*dhat(V)dy> - <dUdy*dVdy>) !
  !                                                              !
  ! Input:  U_,V_,W_,avg_nu_t (flow velocities)                  !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_alpha_momentum_wall_model_dynamic(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64) :: dU,     dUV,     mtau
    Real   (Int64) :: dU_hat, dUV_hat, mtau_hat
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! FOR LES FILTER

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at cell edges in xy plane)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate avg_nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(avg_nu_t,term_2(1:nxg,1:ny,1:1),2)

    ! Compute term_4 = tau_UV at cell edges
    ! term_2 -> nu_t
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1,j,1)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    mtau = Sum( term_4(2:nx-1,1,2:nzg-1) - term_4(2:nx-1,ny,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall    

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k))/(yg(nyg) - yg(nyg-1))*(term_1(i,ny,k) - term_1(i,ny-1,k))/(y(ny) - y(ny-1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    dUV = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall    

    ! product
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average <dU/dy>
    dU = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! FOR TEST FILTER

    ! filter velocities
    Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
    Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
    Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )

    ! apply periodicity in x and z
    Call apply_periodic_bc_x(Uff,1)
    Call apply_periodic_bc_z(Uff,2)
    Call apply_periodic_bc_x(Vff,2)
    Call apply_periodic_bc_z(Vff,2)
    Call apply_periodic_bc_x(Wff,2)
    Call apply_periodic_bc_z(Wff,1)

    !---------------------------------------------------------!
    ! Part 1: compute T_UV -> term_4 (at cell edges in xy plane)
   
    ! compute hat(Suv) at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (Uff(i,  j+1,k) - Uff(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (Vff(i+1,j  ,k) - Vff(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
   
    ! compute eddy viscosity for filtered velocities
    Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t) ! nu_t should be change! 

    ! interpolate avg_nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(avg_nu_t_hat,term_2(1:1,1:ny,1:1),2)

    ! Compute term_4 = T_UV at cell edges
    ! term_2 -> nu_t_hat
    ! term_3 -> Suv_hat
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1,j,1)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average T_UV
    mtau_hat = Sum( term_4(2:nx-1,1,2:nzg-1) - term_4(2:nx-1,ny,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dhat(U)/dy*dhat(V)/dy> at the wall    

    ! interpolate hat(V) to cell edges
    Call interpolate_x(Vff,term_1,1)

    ! product: dhat(U)/dy*dhat(V)/dy (first order)
    ! term_1 -> hat(V)
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (Uff(i,2,k) - Uff(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(Uff(i,nyg,k) - Uff(i,nyg-1,k))/(yg(nyg) - yg(nyg-1))*(term_1(i,ny,k) - term_1(i,ny-1,k))/(y(ny) - y(ny-1))
       End Do
    End Do
    
    ! average <dhat(U)/dy*dhat(V)/dy>
    dUV_hat = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dhat(U)/dy> at the wall    

    ! product
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (Uff(i,2,k) - Uff(i,1,k)) / (yg(2) - yg(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(Uff(i,nyg,k) - Uff(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average <dhat(U)/dy>
    dU_hat = Sum( term(2:nx-1,1:2,2:nzg-1) ) / Real( 2*(nx-2)*(nzg-2), 8)

    !---------------------------------------------------------!
    ! COMPUTE ALPHA_i

    alpha_x(:,:,:) = ( nu*(dU_hat-dU) + mtau - mtau_hat )/( fil_size**2d0*dUV_hat - dUV ) 

    ! clipping 
    alpha_x = Max( alpha_x(2,1,2), 0d0 )
    alpha_x = ( alpha_x(2,1,2) )**0.5d0

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)

  End Subroutine compute_alpha_momentum_wall_model_dynamic

  !--------------------------------------------------------------!
  !           Compute alpha from Cf turbulent correlation        !
  !                                                              !
  !  tau_wall = 0.5*Cf*Uinf^2                                    !
  !  tau_wall = (nu+nu_t)*dUmean/dy_wall - UVmean_wall           !
  !                                                              !
  !  1) Impose viscous stress lower than frac_vis:               !
  !     alpha_x = nu<U>/(frac_vis*tau_wall)                      !
  !                                                              !
  !  2) Impose correct stress at the wall:                       !
  !     alpha_y = (<nu*dU/dy>-tau_wall-<tau_UV>)/                !
  !                <alpha_x*dU/dy*dV/dy>                         !
  !                                                              !
  !  3) alpha_z = alpha_y                                        !
  !                                                              !
  ! frac_vis may be static or dynamic:                           !
  !     From box filtered DNS viscous stress                     !
  !     Assummed U_DNS_wall = 0                                  !
  !     frac = nu*U(Delta_f)/utau^2/Delta_f, Delta_f = 2*Delta   !
  !                                                              !
  ! F. M. White, Viscous Fluid Flow, 2005:                       !
  ! Cf = 0.020*Re_delta99^(-1/6)                                 !
  ! Cf = 0.027*Re_x^(-1/7)                                       !
  !                                                              !
  ! Flags:  iwall_model_nut, 0->no nu_t, 1->XXX, 2->nu dyn.frac  !
  ! Input:  U_, V_, W_, nu_t, frac_vis_wall_model                !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_alpha_u_from_Cf(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64), Dimension(nx) ::  Umean_wall,  Umean_wall_local 
    Real   (Int64), Dimension(nx) :: dUmean_wall, dUmean_wall_local
    Real   (Int64), Dimension(nx) ::   dUdV_wall,   dUdV_wall_local
    Real   (Int64), Dimension(nx) ::        Uref,        Uref_local 
    Real   (Int64), Dimension(nx) ::        mtau,        mtau_local
    Real   (Int64), Dimension(nx) :: tau_wall_ref, Cf
    Real   (Int64) :: frac_vis_dynamic, frac_vis, Err, Uinf
    Real   (Int64) :: alpha_y_1, alpha_y_2
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! compute Cf and tau_wall (maybe do it only once...)
    !
    Uinf         = U_(2,nyg,2)
    Cf           = 0.027d0*(Uinf*x/nu)**(-1d0/7d0) ! at x
    tau_wall_ref = 0.5d0*Cf*Uinf**2d0              ! at x

    !---------------------------------------------------------!   
    ! Reference viscous stress + model at the wall
    frac_vis = frac_vis_wall_model ! better if < 0.05

    !---------------------------------------------------------!
    ! Part 0: compute alpha_x
    Do i = 1,nx_global
       Umean_wall_local(i) = Sum( U_(i,1:2,2:nzg-1) )
    End Do
    Call MPI_Allreduce(Umean_wall_local, Umean_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)    
    Umean_wall = Umean_wall/Real( 2*(nzg_global-2), 8)

    If ( iwall_model_nut == 0 ) Then
       Do i = 1,nx_global       
          alpha_x(i,:,:) = nu*Umean_wall(i)/( frac_vis*tau_wall_ref(i) ) 
       End Do
    Elseif ( iwall_model_nut == 2 ) Then
       Do i = 1,nx_global
          Uref_local(i) = Sum( U_(i,3:4,2:nzg-1) )
       End Do
       Call MPI_Allreduce(Uref_local, Uref, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       Uref = Uref/Real( 2*(nzg_global-2), 8)
       Do i = 1,nx_global
          frac_vis_dynamic = nu*Uref(i)/( tau_wall_ref(i)*y(3) ) 
          alpha_x(i,:,:)   = nu*Umean_wall(i)/( frac_vis_dynamic*tau_wall_ref(i) ) 
       End Do
       frac_vis = frac_vis_dynamic
    End If
    ! avoid negative alpha_x 
    Do i = 1,nx_global
       If ( alpha_x(i,1,1)<0d0 ) Then
          alpha_x(i,:,:) = 1d-3
       End If
    End Do
    ! avoid zero alpha_x (note: there must be a better way of doing this)
    Do i = 1,nx_global
       If ( Abs(alpha_x(i,1,1))<1d-8 ) Then
          alpha_x(i,:,:) = 1d-3
       End If
    End Do

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at x)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from y cell centers to y cell faces 
    Call interpolate_y(nu_t,term_1(1:nxg,1:ny,1:nzg),2)
    ! interpolate nu_t from x cell centers to cell faces 
    Call interpolate_x(term_1(1:nxg,1:ny,1:nzg),term_2(1:nx,1:ny,1:nzg),2)

    ! Compute term_4 = tau_UV = -2*nu_t*Suv at cell edges
    ! term_2 -> nu_t 
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    Do i = 1,nx_global
       mtau_local(i) = Sum( term_4(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(mtau_local, mtau, nx,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( (nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall (at x)

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1,:) = 0d0
    Do i = 1, nx_global
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    Do i = 1,nx_global
       dUdV_wall_local(i) = Sum( term(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(dUdV_wall_local, dUdV_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUdV_wall = dUdV_wall/Real( nzg_global-2, 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall (at x)

    ! product
    term(:,1,:) = 0d0
    Do i = 1, nx_global
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
       End Do
    End Do
    
    ! average <dU/dy>
    Do i = 1, nx_global
       dUmean_wall_local(i) = Sum( term(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(dUmean_wall_local, dUmean_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUmean_wall = dUmean_wall/Real( nzg_global-2, 8)

    !---------------------------------------------------------!
    ! Part 4: compute alpha_y 
    Do i = 2, nxg_global-1
       alpha_y_1      = ( nu*dUmean_wall(i-1) - tau_wall_ref(i-1) - mtau(i-1) )/(alpha_x(i-1,1,2)*dUdV_wall(i-1))
       alpha_y_2      = ( nu*dUmean_wall(i  ) - tau_wall_ref(i  ) - mtau(i  ) )/(alpha_x(i  ,1,2)*dUdV_wall(i  ))
       alpha_y(i,:,:) = 0.5d0*( alpha_y_1 + alpha_y_2 )
    End Do
    alpha_y(  1,:,:) = alpha_y(    2,:,:)
    alpha_y(nxg,:,:) = alpha_y(nxg-1,:,:)

    ! avoid large alpha_y (note: there must be a better way of doing this)
    Do i = 1,nxg_global
       If ( Abs(alpha_y(i,1,1))>2d-3 ) Then
          alpha_y(i,:,:) = 2d-3*alpha_y(i,1,1)/Abs(alpha_y(i,1,1))
       End If
    End Do

    !---------------------------------------------------------!
    ! Part 5: compute alpha_z
    alpha_z = alpha_y

  end Subroutine compute_alpha_u_from_Cf

  !--------------------------------------------------------------!
  !                                                              !
  !      Momentum balance at the wall with only alpha_x          !
  !                                                              !
  !     Compute slip length alpha_x, alpha_y and alpha_z         !
  !              for dui/dy = alpha_i at the wall                !
  !                                                              !
  ! alpha_x /= 0                                                 ! 
  ! alpha_y  = 0                                                 !
  ! alpha_z  = 0                                                 !
  !                                                              !
  ! Equations:                                                   !
  ! -> alpha_x = utau^2/nu                                       !
  !                                                              !
  ! Input:  U_, V_, W_, nu_t,                                    !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_alpha_u_momentum_wall_model(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64), Dimension(nx) :: mtau, mtau_local
    Integer(Int32) :: i, j, k
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from y cell centers to y cell faces 
    Call interpolate_y(nu_t,term_1(1:nxg,1:ny,1:nzg),2)
    ! interpolate nu_t from x cell centers to cell faces 
    Call interpolate_x(term_1(1:nxg,1:ny,1:nzg),term_2(1:nx,1:ny,1:nzg),2)

    ! Compute term_4 = tau_UV = -2*nu_t*Suv at cell edges
    ! term_2 -> nu_t 
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    Do i = 1,nx_global
       mtau_local(i) = Sum( term_4(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(mtau_local, mtau, nx,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( (nzg_global-2), 8)
   
    ! compute alphas
    Do i = 1, nx_global
       If  ( Abs(utau_model(i))<1d-10 ) Then
          alpha_x(i,:,:) = 0d0
       Else
          alpha_x(i,:,:) = (utau_model(i)**2d0 + mtau(i))/nu
       End If
    End Do
    alpha_y = 0d0
    alpha_z = 0d0

    ! avoid zeros
    Do i = 1, nx_global
       alpha_x(i,:,:) = Max(alpha_x(i,1,1),0d0)
    End Do

  End Subroutine compute_alpha_u_momentum_wall_model

  !--------------------------------------------------------------!
  !           Compute alpha from Cf turbulent correlation        !
  !                                                              !
  !  tau_wall = 0.5*Cf*Uinf^2                                    !
  !  tau_wall = (nu+nu_t)*dUmean/dy_wall - UVmean_wall           !
  !                                                              !
  !  1) Impose viscous stress lower than frac_vis:               !
  !     alpha_x = nu<U>/(frac_vis*tau_wall)                      !
  !                                                              !
  !  2) Impose correct stress at the wall:                       !
  !     alpha_y = (<nu*dU/dy>-tau_wall-<tau_UV>)/                !
  !                <alpha_x*dU/dy*dV/dy>                         !
  !                                                              !
  !  3) alpha_z = alpha_y                                        !
  !                                                              !
  ! frac_vis may be static or dynamic:                           !
  !     From box filtered DNS viscous stress                     !
  !     Assummed U_DNS_wall = 0                                  !
  !     frac = nu*U(Delta_f)/utau^2/Delta_f, Delta_f = 2*Delta   !
  !                                                              !
  ! F. M. White, Viscous Fluid Flow, 2005:                       !
  ! Cf = 0.020*Re_delta99^(-1/6)                                 !
  ! Cf = 0.027*Re_x^(-1/7)                                       !
  !                                                              !
  ! Flags:  iwall_model_nut, 0->no nu_t, 1->XXX, 2->nu dyn.frac  !
  ! Input:  U_, V_, W_, nu_t, frac_vis_wall_model                !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_alpha_v_from_Cf(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64), Dimension(nx) ::  Umean_wall,  Umean_wall_local 
    Real   (Int64), Dimension(nx) :: dUmean_wall, dUmean_wall_local
    Real   (Int64), Dimension(nx) ::   dUdV_wall,   dUdV_wall_local
    Real   (Int64), Dimension(nx) ::        Uref,        Uref_local 
    Real   (Int64), Dimension(nx) ::        mtau,        mtau_local
    Real   (Int64), Dimension(nx) :: tau_wall_ref, Cf
    Real   (Int64) :: frac_vis_dynamic, frac_vis, Err, Uinf
    Real   (Int64) :: alpha_y_1, alpha_y_2
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! compute Cf and tau_wall (maybe do it only once...)
    !
    Uinf         = U_(2,nyg,2)
    Cf           = 0.027d0*(Uinf*x/nu)**(-1d0/7d0) ! at x
    tau_wall_ref = 0.5d0*Cf*Uinf**2d0              ! at x

    !---------------------------------------------------------!   
    ! Reference viscous stress + model at the wall
    frac_vis = frac_vis_wall_model ! better if < 0.05

    !---------------------------------------------------------!
    ! Part 0: compute alpha_x
    Do i = 1,nx_global
       Umean_wall_local(i) = Sum( U_(i,1:2,2:nzg-1) )
    End Do
    Call MPI_Allreduce(Umean_wall_local, Umean_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)    
    Umean_wall = Umean_wall/Real( 2*(nzg_global-2), 8)

    If ( iwall_model_nut == 0 ) Then
       Do i = 1,nx_global       
          alpha_x(i,:,:) = nu*Umean_wall(i)/( frac_vis*tau_wall_ref(i) ) 
       End Do
    Elseif ( iwall_model_nut == 2 ) Then
       Do i = 1,nx_global
          Uref_local(i) = Sum( U_(i,3:4,2:nzg-1) )
       End Do
       Call MPI_Allreduce(Uref_local, Uref, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       Uref = Uref/Real( 2*(nzg_global-2), 8)
       Do i = 1,nx_global
          frac_vis_dynamic = nu*Uref(i)/( tau_wall_ref(i)*y(3) ) 
          alpha_x(i,:,:)   = nu*Umean_wall(i)/( frac_vis_dynamic*tau_wall_ref(i) ) 
       End Do
       frac_vis = frac_vis_dynamic
    End If
    ! avoid negative alpha_x 
    Do i = 1,nx_global
       If ( alpha_x(i,1,1)<0d0 ) Then
          alpha_x(i,:,:) = 1d-3
       End If
    End Do
    ! avoid zero alpha_x (note: there must be a better way of doing this)
    Do i = 1,nx_global
       If ( Abs(alpha_x(i,1,1))<1d-8 ) Then
          alpha_x(i,:,:) = 1d-3
       End If
    End Do

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at x)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
!             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
!                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) )
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from y cell centers to y cell faces 
    Call interpolate_y(nu_t,term_1(1:nxg,1:ny,1:nzg),2)
    ! interpolate nu_t from x cell centers to cell faces 
    Call interpolate_x(term_1(1:nxg,1:ny,1:nzg),term_2(1:nx,1:ny,1:nzg),2)

    ! Compute term_4 = tau_UV = -2*nu_t*Suv at cell edges
    ! term_2 -> nu_t 
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    Do i = 1,nx_global
       mtau_local(i) = Sum( term_4(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(mtau_local, mtau, nx,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( (nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall (at x)

    ! interpolate V to cell edges
    Call interpolate_x(Vo,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1,:) = 0d0
    Do i = 1, nx_global
       Do k = 2, nzg-1
          ! bottom wall
!          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          term(i,1,k) =  0.5d0*(U_(i,2,k) + U_(i,1,k))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    Do i = 1,nx_global
       dUdV_wall_local(i) = Sum( term(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(dUdV_wall_local, dUdV_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUdV_wall = dUdV_wall/Real( nzg_global-2, 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall (at x)

    ! product
    term(:,1,:) = 0d0
    Do i = 1, nx_global
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
       End Do
    End Do
    
    ! average <dU/dy>
    Do i = 1, nx_global
       dUmean_wall_local(i) = Sum( term(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(dUmean_wall_local, dUmean_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUmean_wall = dUmean_wall/Real( nzg_global-2, 8)

    !---------------------------------------------------------!
    ! Part 4: compute alpha_y 
    Do i = 2, nxg_global-1
 !      alpha_y_1      = ( nu*dUmean_wall(i-1) - tau_wall_ref(i-1) - mtau(i-1) + Umean_wall(i)*beta_y )/(alpha_x(i-1,1,2)*dUdV_wall(i-1))
 !      alpha_y_2      = ( nu*dUmean_wall(i  ) - tau_wall_ref(i  ) - mtau(i  ) + Umean_wall(i)*beta_y )/(alpha_x(i  ,1,2)*dUdV_wall(i  ))
       alpha_y_1      = ( nu*dUmean_wall(i-1) - tau_wall_ref(i-1) - mtau(i-1) + Umean_wall(i-1)*beta_y )/(dUdV_wall(i-1))
       alpha_y_2      = ( nu*dUmean_wall(i  ) - tau_wall_ref(i  ) - mtau(i  ) + Umean_wall(i  )*beta_y )/(dUdV_wall(i  ))
       alpha_y(i,:,:) = 0.5d0*( alpha_y_1 + alpha_y_2 )
    End Do
    alpha_y(  1,:,:) = alpha_y(    2,:,:)
    alpha_y(nxg,:,:) = alpha_y(nxg-1,:,:)

    ! avoid large alpha_y (note: there must be a better way of doing this)
    !Do i = 1,nxg_global
    !   If ( Abs(alpha_y(i,1,1))>2d-3 ) Then
    !      alpha_y(i,:,:) = 2d-3*alpha_y(i,1,1)/Abs(alpha_y(i,1,1))
    !   End If
    !End Do

    !---------------------------------------------------------!
    ! Part 5: compute alpha_z
    alpha_z = alpha_y

  end Subroutine compute_alpha_v_from_Cf

  !--------------------------------------------------------------!
  !           Compute alpha from Cf turbulent correlation        !
  !                                                              !
  !  tau_wall = 0.5*Cf*Uinf^2                                    !
  !  tau_wall = (nu+nu_t)*dUmean/dy_wall - UVmean_wall           !
  !                                                              !
  !  1) Impose viscous stress lower than frac_vis:               !
  !     alpha_x = nu<U>/(frac_vis*tau_wall)                      !
  !                                                              !
  !  2) Impose correct stress at the wall:                       !
  !     V = a_sca*V^n - beta_sca (scaled velocity)               !
  !                                                              !
  !     a = (nu*dUmean_wall - tau_wall + beta*Umean_wall +       !
  !          <tau_model>)/<UV_wall>                              !
  !                                                              !
  !     alpha_y = < (V(1)+beta)/(V(2)-V(1))*dy >                 !
  !                                                              !
  !  3) alpha_z = alpha_y                                        !
  !                                                              !
  ! frac_vis may be static or dynamic:                           !
  !     From box filtered DNS viscous stress                     !
  !     Assummed U_DNS_wall = 0                                  !
  !     frac = nu*U(Delta_f)/utau^2/Delta_f, Delta_f = 2*Delta   !
  !                                                              !
  ! F. M. White, Viscous Fluid Flow, 2005:                       !
  ! Cf = 0.020*Re_delta99^(-1/6)                                 !
  ! Cf = 0.027*Re_x^(-1/7)                                       !
  !                                                              !
  ! Input:  U_, V_, W_, nu_t, frac_vis_wall_model                !
  ! Output: alpha_x, V_bottom, alpha_z                           !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_v_from_Cf(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64), Dimension(nx) ::   Umean_wall,  Umean_wall_local 
    Real   (Int64), Dimension(nx) ::  dUmean_wall, dUmean_wall_local
    Real   (Int64), Dimension(nx) ::      UV_wall,     UV_wall_local
    Real   (Int64), Dimension(nx) ::         Uref,        Uref_local 
    Real   (Int64), Dimension(nx) ::         mtau,        mtau_local
    Real   (Int64), Dimension(nx) :: tau_wall_ref, Cf
    Real   (Int64), Dimension(nxg)::        a_sca
    Real   (Int64) :: frac_vis_dynamic, frac_vis, Err, Uinf, T_mean
    Real   (Int64) :: a_sca_1, a_sca_2, beta_sca, beta_sca_local
    Integer(Int32) :: i, j, k, kk

    !---------------------------------------------------------!
    ! compute Cf and tau_wall (maybe do it only once...)
    !
    Uinf         = U_(2,nyg,2)
    Cf           = 0.027d0*(Uinf*x/nu)**(-1d0/7d0) ! at x
    tau_wall_ref = 0.5d0*Cf*Uinf**2d0              ! at x

    !---------------------------------------------------------!   
    ! Reference viscous stress + model at the wall
    frac_vis = frac_vis_wall_model ! better if < 0.05

    !---------------------------------------------------------!
    ! Part 0: compute alpha_x
    Do i = 1,nx_global
       Umean_wall_local(i) = Sum( U_(i,1:2,2:nzg-1) )
    End Do
    Call MPI_Allreduce(Umean_wall_local, Umean_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)    
    Umean_wall = Umean_wall/Real( 2*(nzg_global-2), 8)

    If ( iwall_model_nut == 0 ) Then
       Do i = 1,nx_global       
          alpha_x(i,:,:) = nu*Umean_wall(i)/( frac_vis*tau_wall_ref(i) ) 
       End Do
    Elseif ( iwall_model_nut == 2 ) Then
       Do i = 1,nx_global
          Uref_local(i) = Sum( U_(i,3:4,2:nzg-1) )
       End Do
       Call MPI_Allreduce(Uref_local, Uref, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       Uref = Uref/Real( 2*(nzg_global-2), 8)
       Do i = 1,nx_global
          frac_vis_dynamic = nu*Uref(i)/( tau_wall_ref(i)*y(3) ) 
          alpha_x(i,:,:)   = nu*Umean_wall(i)/( frac_vis_dynamic*tau_wall_ref(i) ) 
       End Do
       frac_vis = frac_vis_dynamic
    End If
    ! avoid negative alpha_x 
    Do i = 1,nx_global
       If ( alpha_x(i,1,1)<0d0 ) Then
          alpha_x(i,:,:) = 0d0
       End If
    End Do

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at x)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
!             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
!                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) )
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from y cell centers to y cell faces 
    Call interpolate_y(nu_t,term_1(1:nxg,1:ny,1:nzg),2)
    ! interpolate nu_t from x cell centers to cell faces 
    Call interpolate_x(term_1(1:nxg,1:ny,1:nzg),term_2(1:nx,1:ny,1:nzg),2)

    ! Compute term_4 = tau_UV = -2*nu_t*Suv at cell edges
    ! term_2 -> nu_t 
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    Do i = 1,nx_global
       mtau_local(i) = Sum( term_4(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(mtau_local, mtau, nx,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( (nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <U*V> at the wall (at x)

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: U*V
    ! term_1 -> V
    term(:,1,:) = 0d0
    Do i = 1, nx_global
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  0.5d0*(U_(i,2,k) + U_(i,1,k))*term_1(i,1,k)
       End Do
    End Do
    
    ! average <U*V>
    Do i = 1,nx_global
       UV_wall_local(i) = Sum( term(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(UV_wall_local, UV_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    UV_wall = UV_wall/Real( nzg_global-2, 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall (at x)

    ! product
    term(:,1,:) = 0d0
    Do i = 1, nx_global
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
       End Do
    End Do
    
    ! average <dU/dy>
    Do i = 1, nx_global
       dUmean_wall_local(i) = Sum( term(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(dUmean_wall_local, dUmean_wall, nx_global,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUmean_wall = dUmean_wall/Real( nzg_global-2, 8)

    ! accumulated variables
    If ( istep==1 ) Then
       dUmean_wall_T = dUmean_wall
       mtau_T        = mtau
       UV_wall_T     = UV_wall
    End If
    T_mean        = alpha_std
    dUmean_wall_T = dt/T_mean*dUmean_wall + (1d0-dt/T_mean)*dUmean_wall_T
    mtau_T        = dt/T_mean*mtau        + (1d0-dt/T_mean)*mtau_T
    UV_wall_T     = dt/T_mean*UV_wall     + (1d0-dt/T_mean)*UV_wall_T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dUmean_wall_T = dUmean_wall
    mtau_T        = mtau
    UV_wall_T     = UV_wall
    Do i = 1, nx
       If ( Abs(UV_wall_T(i))<1e-8 ) UV_wall_T(i) = 1d-8
    End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !---------------------------------------------------------!
    ! Part 4: compute a_sca (beta_sca assumed zero)
    Do i = 2, nxg_global-1
       a_sca_1  = ( nu*dUmean_wall_T(i-1) - tau_wall_ref(i-1) - mtau_T(i-1) )/(Umean_wall(i-1))
       a_sca_2  = ( nu*dUmean_wall_T(i  ) - tau_wall_ref(i  ) - mtau_T(i  ) )/(Umean_wall(i  ))
       a_sca(i) = 0.5d0*( a_sca_1 + a_sca_2 )
    End Do
    a_sca(  1) = a_sca(    2)
    a_sca(nxg) = a_sca(nxg-1)

    !---------------------------------------------------------!
    ! Part 5: V_bottom
    Do i = 1, nxg
       V_bottom(i,:) = a_sca(i) 
    End Do

    ! clipping
    Do i = 1, nxg
       Do k = 1,nzg
          If ( V_bottom(i,k)> 0.5d0 ) V_bottom(i,k) =  0.5d0
          If ( V_bottom(i,k)<-0.5d0 ) V_bottom(i,k) = -0.5d0
       End Do
    End Do
    
    kk = 1
    If ( myid==(nprocs-1) ) kk = 2 ! last periodic point in z excluded for mass conservation
    beta_sca_local = Sum( V_bottom(2:nxg-2,2:nzg-kk) )
    Call MPI_Allreduce(beta_sca_local,beta_sca,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    beta_sca = beta_sca / Real( (nxg-3)*(nzg_global-3) , 8 )
!!!!!!!!!!!!!!!!!!!!!
!    V_bottom = V_bottom - beta_sca    
!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    If (myid == 0) write(*,*) 'V_bottom', V_bottom(10:15,3)
    If (myid == 0 .and. Mod(istep,50)==0 ) write(*,*) 'beta_sca', beta_sca
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !---------------------------------------------------------!
    ! Part 6: estimated alpha_y (need to be averaged)
    Do i=1,nxg
       alpha_y(i,1,:) = V_bottom(i,:) !( V_(i,1,:) + beta_sca )/(V_(i,2,:)-V_(i,1,:))*(y(2)-y(1))
    End Do
    
    !---------------------------------------------------------!
    ! Part 5: compute alpha_z
    alpha_z = alpha_x/10d0

  end Subroutine compute_v_from_Cf

  !--------------------------------------------------------------------------!
  !                         Log-layer wall model                             !
  !                                                                          !
  !  Log-layer : <U(jref)> = utau/kappa*ln(y(jref)*utau/nu) + B*utau         !
  !              tau_wall  = utau^2                                          !
  !                                                                          !
  ! kappa = 0.38                                                             !
  ! B     = 5.2                                                              !
  ! Solved with iterative Newton method                                      !
  !                                                                          !
  ! Input:  U_                                                               !
  ! Output: utau_model                                                       !
  !                                                                          !
  !--------------------------------------------------------------------------!
  Subroutine compute_law_of_the_wall_model(U_)

    Real(Int64), Dimension(:,:,:), Intent(In) :: U_

    ! local variables
    Real   (Int64) :: utau_model0, utau_model1, Umean_model_local
    Integer(Int32) :: i, jref, iters, iters_max

    ! parameters
    jref        = 3
    kappa_model = 0.38d0
    B_model     = 5.20d0
    yg_model    = yg(jref)
    iters_max   = 0

    ! solve equation at each x-location
    Do i = 2, nx_global-1
       utau_model0       = utau_model(i)
       If ( utau_model(i) /= utau_model(i) ) utau_model(i) = 0d0
       Umean_model_local = Sum( U_(i,jref,2:nzg-1) )
       Call MPI_Allreduce(Umean_model_local,Umean_model,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       Umean_model = Umean_model/Real( nzg_global-2, 8) 
       Call Newton_iter_solver(f_law_of_wall, df_law_of_wall, utau_model0, utau_model1, iters, .false.)
       utau_model(i) = utau_model1
       If ( utau_model(i) /= utau_model(i) ) utau_model(i) = 0d0
       iters_max     = Max(iters_max,iters)
    End Do
    utau_model(1)         = utau_model(2)
    utau_model(nx_global) = utau_model(nx_global-1)

    If ( iters_max==30 .and. myid==0 ) Then
       Write(*,*) 'Warning: Newton solver not converged!'
    End If

  End Subroutine compute_law_of_the_wall_model

  !--------------------------------------------------------------!
  !        Compute alpha=dU/dy from Cf turbulent correlation     !
  !                                                              !
  !  tau_wall = 0.5*Cf*Uinf^2                                    !
  !                                                              !  
  !  Impose correct stress at the wall:                          !
  !   dU/dy = alpha_x = <0.5*Cf*Uinf^2/(nu + nut)>               !
  !                                                              !
  ! F. M. White, Viscous Fluid Flow, 2005:                       !
  ! Cf = 0.020*Re_delta99^(-1/6)                                 !
  ! Cf = 0.027*Re_x^(-1/7)                                       !
  !                                                              !
  ! Input:  U_, V_, W_, nu_t                                     !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_exact_Neumann_from_Cf(U_,V_,W_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    ! local variables
    Real   (Int64), Dimension(nx) :: mtau, mtau_local
    Real   (Int64), Dimension(nx) :: tau_wall_ref, Cf_ref
    Real   (Int64) :: Uinf
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! compute Cf and tau_wall (maybe do it only once...)
    Uinf         = 1d0 !U_(2,nyg,2)
    Cf_ref       = 0.027d0*(Uinf*x/nu)**(-1d0/7d0) ! at x
    tau_wall_ref = 0.5d0*Cf_ref*Uinf**2d0          ! at x

!!!!!!!!!!!!!!!!!!!!!!!!
!    if (myid==0) Then
!       write(*,*) '------------------Cf_ref ii'
!       Do i=1,nx
!          write(*,*) Cf_ref(i)
!       end Do
!    end if
!!!!!!!!!!!!!!!!!!!!!!!!

    !---------------------------------------------------------!
    ! interpolate nu_t from y cell centers to y cell faces
    Call interpolate_y(nu_t,term_1(1:nxg,1:ny,1:nzg),2)
    !Call interpolate_y(avg_nu_t,term_1(1:nxg,1:ny,1:nzg),2)
    ! interpolate nu_t from x cell centers to cell faces 
    Call interpolate_x(term_1(1:nxg,1:ny,1:nzg),term_2(1:nx,1:ny,1:nzg),2)
    
    !---------------------------------------------------------!
    ! Compute term_4 = 0.5*Cf*Uinf^2/(nu + nut) at cell centers
    !         term_4 = tau_wall/(nu+nut) 
    Do j = 1, 2
       Do k = 1, nzg
          !term_4(1:nx,j,k) = tau_wall_ref(1:nx)/(nu + term_2(1:nx,j,k))
          term_4(1:nx,j,k) = (nu + term_2(1:nx,j,k))/tau_wall_ref(1:nx)
       end Do
    End Do
    
    ! average tau_m in z
    Do i = 1,nx
       mtau_local(i) = Sum( term_4(i,1,2:nzg-1) )
    End Do
    Call MPI_Allreduce(mtau_local, mtau, nx, MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    !mtau = mtau/Real( nzg_global-2, 8)
    mtau = 1d0/(mtau/Real( nzg_global-2, 8))

    !---------------------------------------------------------!
    ! Part 1: compute alpha_x
    Do i=1,nx
       alpha_x(i,:,:) = mtau(i)
    End Do
    
    !---------------------------------------------------------!
    ! Part 2: compute alpha_y/z
    alpha_z = 0d0
    alpha_y = 0d0

  end Subroutine compute_exact_Neumann_from_Cf
  
  !-------------------------------------------------------------!
  !                                                             !
  !            Compute Neumann boundary conditions              !
  !    for pseudo-pressure when slip-wall model is active       !
  !                                                             !
  ! This has to be called every sub-step                        !
  !                                                             !
  ! Conditions bottom wall:                                     !
  !                                                             !
  !     V1 = V1*  -  (p2-p1)/(yg(2)-yg(1))                      !
  !     V2 = V2*  -  (p3-p2)/(yg(3)-yg(2))                      !
  !     V1 = alpha_y*(V2-V1)/(y (2)-y (1))                      !
  !                                                             !
  !     => p1 = p_b2*p2 + p_bc3*p3                              !
  !        p_bc2   = 1 + beta*Delta_r                           !
  !        p_bc3   =   - beta*Delta_r                           !
  !        Delta_r = ( yg(2)-yg(1) )/( yg(3)-yg(2) )            !
  !        beta    = alpha_y/(alpha_y + y(2)-y(1) )             !
  !                                                             !
  !     V* -> velocity without pressure                         !
  !                                                             !
  ! Equation for first interior points:                         !
  !                                                             !
  !    (a + c*p_bc3)*p3 + (b + c*p_bc2)*p2 = rhs_p2             !
  !                                                             !
  !                                                             !
  ! Conditions top wall: (n->ny, ng->nyg)                       !
  !                                                             !
  !     V(n)   = V(  n)* -(p(ng)  -p(ng-1))/(yg(ng)  -yg(ng-1)) !
  !     V(n-1) = V(n-1)* -(p(ng-1)-p(ng-2))/(yg(ng-1)-yg(ng-2)) !
  !     V(n)   = alpha_y*(V(n)-V(n-1))/(y(n)-y(n-1))            !
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
  ! NOT USED FOR NOW
  Subroutine compute_pseudo_pressure_bc_for_robin_bc

    ! local variables
    Real   (Int64) :: a, b, c
    Real   (Int64) :: beta, Delta_r, alphad
    Real   (Int64) :: p_bc2, p_bc3, p_bcn, p_bcn1
    Integer(Int64) :: j    

    ! bottom wall
    j        = 2 
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    Delta_r  = ( yg(2)-yg(1) )/( yg(3)-yg(2) )
    alphad   = alpha_y(2,2,2)
    beta     = alphad/(alphad + y(2)-y(1) )
    p_bc2    = 1d0 + beta*Delta_r
    p_bc3    =     - beta*Delta_r

    Dyy(2,2) = b + c*p_bc2
    Dyy(2,3) = a + c*p_bc3
    
    ! top wall
    j        = nyg-1
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    alphad   = -alpha_y(2,2,2)
    Delta_r  = ( yg(nyg) - yg(nyg-1) )/( yg(nyg-1) - yg(nyg-2) )
    beta     = alphad/( alphad - (y(ny)-y(ny-1)) )
    p_bcn    = 1d0 + beta*Delta_r
    p_bcn1   =     - beta*Delta_r
    
    Dyy(nyg-1,nyg-1) = b + a*p_bcn
    Dyy(nyg-1,nyg-2) = c + a*p_bcn1
    
  End Subroutine compute_pseudo_pressure_bc_for_robin_bc
      
End Module wallmodel
