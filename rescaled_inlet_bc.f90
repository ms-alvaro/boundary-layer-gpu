!--------------------------------------------------------!
!      Lund's rescaling inlet boundary condition         !
!--------------------------------------------------------!
Module Lund_rescaled_bc

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global, Only : nx, nxg, ny, nyg, nz, nzg, nzg_global,     &
                     i_rescale, yg, y, nu, dt, istep, nmonitor, &
                     Umean_resc_T,  Vmean_resc_T, T_resc,       &
                     Umean_inlet_T, Vmean_inlet_T
  Use mpi

  ! prevent implicit typing
  Implicit None
  
  ! set private variables
  Private
  Public  :: compute_rescaled_inflow

Contains

  !-------------------------------------------------------------------!
  !            Compute inflow veloties from Lund's rescaling          !
  !                                                                   !
  !  U     = Umean + Ufluc                                            !
  !  gamma = utau(x_inlet)/utau(x_resc)                               !
  !                                                                   !
  !  Means:                                                           !
  !     Uinner        = utau(x)*f_1(y+)                               !
  !     Uinf - Uouter = utau(x)*f_2(y/delta)                          !
  !     Vinner        = Uinf*f_3(y+)                                  !
  !     Vouter        = Uinf*f_4(y/delta)                             !
  !                                                                   !
  !     Umean_inner_inlet  = gamma*Umean_resc(y+)                     !
  !     Umean_outer_inlet  = gamma*Umean_resc(y/delta)+(1-gamma)*Uinf !
  !     Vmean_inner_inlet  = Vmean_resc(y+)                           !
  !     Vmean_outer_inlet  = Vmean_resc(y/delta)                      !
  !                                                                   !
  !  Fluctuations:                                                    !
  !     Ufluc_inner = utau(x)*g(y+)                                   !
  !     Ufluc_outer = utau(x)*h(y/delta)                              !
  !                                                                   !
  !     Ufluc_inner_inlet = gamma*Ufluc_resc(y+)                      !
  !     Ufluc_outer_inlet = gamma*Ufluc_resc(y/delta)                 !
  !                                                                   !
  !     (same for V and W)                                            !
  !                                                                   !
  !  Consistent utau at inlet: (Ludwing-Tillmann)                     !
  !     utau_inlet = utau_resc*(theta_resc/theta_inlet)^{1/(2(n-1))}  !
  !     n = 5                                                         !
  !                                                                   !
  !  Means are weighted in time:                                      !
  !     Umean^(n+1) = dt/T_resc*<U^(n+1)>_z + (1-dt/T_resc)*Umean^n   !
  !                                                                   !
  !  Input:  U, V, W, delta_inlet                                     !
  !  Output: U_inlet, V_inlet, W_inlet                                !  
  !                                                                   !
  !-------------------------------------------------------------------!
  Subroutine compute_rescaled_inflow(U,V,W,U_inlet,V_inlet,W_inlet,delta_inlet,iflag_fluc)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In)  :: U
    Real(Int64), Dimension(nxg,ny ,nzg), Intent(In)  :: V
    Real(Int64), Dimension(nxg,nyg,nz ), Intent(In)  :: W
    Real(Int64), Intent(In) :: delta_inlet

    Real(Int64), Dimension( nyg,nzg), Intent(Out) :: U_inlet
    Real(Int64), Dimension( ny ,nzg), Intent(Out) :: V_inlet
    Real(Int64), Dimension( nyg,nz ), Intent(Out) :: W_inlet    

    Integer(Int32), Intent(In) :: iflag_fluc
    
    ! local variables
    Integer(Int32)                  :: j, jref
    Real   (Int64)                  :: utau_inlet, utau_resc, delta_resc, w1, w0
    Real   (Int64)                  :: alpha, b, gamma, UVmean_wall, UVmean_wall_local, dUmean_wall 
    Real   (Int64)                  :: Uinf, Umean99, theta_resc, theta_inlet, delta_inlet2
        ! rescaled and inlet
    Real(Int64), Dimension(nyg,  1) :: Umean_resc, Umean_inlet, temp
    Real(Int64), Dimension(ny ,  1) :: Vmean_resc, Vmean_inlet
    Real(Int64), Dimension(nyg,  1) :: Umean_resc_local, Umean_inlet_local
    Real(Int64), Dimension(ny ,  1) :: Vmean_resc_local, Vmean_inlet_local
    Real(Int64), Dimension(nyg,nzg) :: U_resc, Ufluc_resc, Ufluc_inlet
    Real(Int64), Dimension(ny ,nzg) :: V_resc, Vfluc_resc, Vfluc_inlet
    Real(Int64), Dimension(nyg,nz ) :: W_resc, Wfluc_resc, Wfluc_inlet
        ! inner and outer
    Real(Int64), Dimension(nyg,  1) :: Umean_inner, Umean_outer
    Real(Int64), Dimension(ny ,  1) :: Vmean_inner, Vmean_outer
    Real(Int64), Dimension(nyg,nzg) :: Ufluc_inner, Ufluc_outer
    Real(Int64), Dimension(ny ,nzg) :: Vfluc_inner, Vfluc_outer
    Real(Int64), Dimension(nyg,nz ) :: Wfluc_inner, Wfluc_outer
        ! y meshes
    Real(Int64), Dimension(nyg)     :: ygp_inlet, ygp_resc, etag_inlet, etag_resc, Weig
    Real(Int64), Dimension(ny )     :: yp_inlet , yp_resc , eta_inlet , eta_resc , Wei
    
    !---------------------------------------------------------------------!
    ! store velocity planes at x rescaled location
    U_resc = U(i_rescale,:,:)
    V_resc = V(i_rescale,:,:)
    W_resc = W(i_rescale,:,:)

    ! store velocity planes at x inlet location
    U_inlet = U(1,:,:)
    V_inlet = V(1,:,:)
    W_inlet = W(1,:,:)

    !---------------------------------------------------------------------!
    ! compute means
    Do j=1,nyg
       Umean_resc_local (j,1) = Sum(U_resc (j,2:nzg-1))
       Umean_inlet_local(j,1) = Sum(U_inlet(j,2:nzg-1))
    End Do
    Do j=1,ny
       Vmean_resc_local (j,1) = Sum(V_resc (j,2:nzg-1))
       Vmean_inlet_local(j,1) = Sum(V_inlet(j,2:nzg-1))
    End Do
    
    Call MPI_Allreduce(Umean_resc_local, Umean_resc, nyg,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_Allreduce(Vmean_resc_local, Vmean_resc, ny ,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_Allreduce(Umean_inlet_local,Umean_inlet,nyg,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_Allreduce(Vmean_inlet_local,Vmean_inlet,ny ,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    
    Umean_resc  = Umean_resc /Real(nzg_global-2,8)
    Vmean_resc  = Vmean_resc /Real(nzg_global-2,8)
    Umean_inlet = Umean_inlet/Real(nzg_global-2,8)
    Vmean_inlet = Vmean_inlet/Real(nzg_global-2,8)
    
    ! garranz: this may be the problem when we change the BC
    !Uinf = Umean_resc(nyg,1)
    
    ! Let us change pick a lower point (this does not work, either)
    Uinf = Umean_resc(nyg-3,1)
    
    ! Let us fix it to 1
    !Uinf = 1d0

    ! time average
    ! NOTE: should be changed to read from file
    !       in order to properly restart the simulation
    Umean_resc_T  = dt/T_resc*Umean_resc  + (1d0-dt/T_resc)*Umean_resc_T
    Vmean_resc_T  = dt/T_resc*Vmean_resc  + (1d0-dt/T_resc)*Vmean_resc_T

    Umean_inlet_T = dt/T_resc*Umean_inlet + (1d0-dt/T_resc)*Umean_inlet_T
    Vmean_inlet_T = dt/T_resc*Vmean_inlet + (1d0-dt/T_resc)*Vmean_inlet_T

    !---------------------------------------------------------------------!    
    ! fluctuating part
    Do j = 1,nyg
       Ufluc_resc (j,:) = U_resc (j,:) - Umean_resc_T (j,1)
       Ufluc_inlet(j,:) = U_inlet(j,:) - Umean_inlet_T(j,1)
       Wfluc_resc (j,:) = W_resc (j,:)
       Wfluc_inlet(j,:) = W_inlet(j,:)
    End Do    
    Do j = 1,ny
       Vfluc_resc (j,:) = V_resc (j,:) - Vmean_resc_T(j,1)
       Vfluc_inlet(j,:) = V_inlet(j,:) - Vmean_inlet (j,1)
    End Do

    !---------------------------------------------------------------------!
    ! compute delta_99 at i_rescale
    Umean99 = 0.99d0*Uinf
    jref    = 0
    Do j = 1,nyg
       If ( Umean_resc_T(j,1)>=Umean99 ) then
          jref = j
          Exit
       End If
    End Do
    If ( jref < 2 ) Then 
       jref = nyg/2
!       If ( myid==0 ) Write(*,*) 'WARNING: no j for delta_resc'
    End If
    w1         = ( Umean99 - Umean_resc_T(jref-1,1) )/( Umean_resc_T(jref,1) - Umean_resc_T(jref-1,1) )
    w0         = 1d0 - w1
    delta_resc = w1*yg(jref) + w0*yg(jref-1)
    If ( delta_resc<0 ) Then
       delta_resc = Abs(delta_resc)
!       If ( myid==0 ) Write(*,*) 'WARNING: delta_resc<0'
    End If

    ! compute delta_99 at inlet
    jref = 0
    Do j = 1,nyg
       If ( Umean_inlet_T(j,1)>=Umean99 ) then
          jref = j
          Exit
       End If
    End Do
    If ( jref < 2 ) Then 
       jref = nyg/2
!       If ( myid==0 ) Write(*,*) 'WARNING: no j for delta_inlet2'
    End If
    w1           = ( Umean99 - Umean_inlet_T(jref-1,1) )/( Umean_inlet_T(jref,1) - Umean_inlet_T(jref-1,1) )
    w0           = 1d0 - w1
    delta_inlet2 = w1*yg(jref) + w0*yg(jref-1)
    If ( delta_inlet2<0 ) Then
       delta_inlet2 = Abs(delta_inlet2)
!       If ( myid==0 ) Write(*,*) 'WARNING: delta_inlet2<0'
    End If

    !---------------------------------------------------------------------!
    ! compute momentum thickness at i_rescale
    ! trapezoidal rule
    ! Q: should I use the ghost cells?
    theta_resc = 0d0
    temp       = Umean_resc_T/Uinf*(1d0 - Umean_resc_T/Uinf)
    Do j = 2,nyg
       theta_resc = theta_resc + 0.5d0*( temp(j,1) + temp(j-1,1) )*( yg(j)-yg(j-1) )
    End Do
    If ( theta_resc<0 ) Then
       theta_resc = Abs(theta_resc)
!       If ( myid==0 ) Write(*,*) 'WARNING: theta_resc<0'
    End If

    ! compute momentum thickness at inlet
    theta_inlet = 0d0
    temp        = Umean_inlet_T/Uinf*(1d0 - Umean_inlet_T/Uinf)
    Do j = 2, nyg
       theta_inlet = theta_inlet + 0.5d0*( temp(j,1) + temp(j-1,1) )*( yg(j)-yg(j-1) )
    End Do
    If ( theta_inlet<0 ) Then
       theta_inlet = Abs(theta_inlet)
!       If ( myid==0 ) Write(*,*) 'WARNING: theta_inlet<0'
    End If

    !---------------------------------------------------------------------!
    ! compute utau at i_rescale
    UVmean_wall_local = Sum( 0.5d0*(U_resc(1,2:nzg-1)+U_resc(2,2:nzg-1))*V_resc(1,2:nzg-1) )  
    Call MPI_Allreduce(UVmean_wall_local,UVmean_wall,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
   
    UVmean_wall = UVmean_wall/Real(nzg_global-2,8)
    dUmean_wall = ( Umean_resc_T(2,1)-Umean_resc_T(1,1) )/( yg(2)-yg(1) )
    utau_resc   = Abs( nu*dUmean_wall - UVmean_wall )**0.5d0

    !---------------------------------------------------------------------!
    ! estimate correct utau at inlet (from Ludwig-Tillmann correlation)
    utau_inlet = utau_resc*(theta_resc/theta_inlet)**(1d0/8d0)
    gamma      = utau_inlet/utau_resc    
    
    !---------------------------------------------------------------------!
    ! y and yg wall-normal grid in plus units
    ygp_inlet  = utau_inlet*yg/nu
    yp_inlet   = utau_inlet*y /nu
    ygp_resc   = utau_resc *yg/nu
    yp_resc    = utau_resc *y /nu

    ! y and yg wall-normal grid in delta_99 units
    etag_inlet = yg/delta_inlet
    eta_inlet  = y /delta_inlet
    etag_resc  = yg/delta_resc
    eta_resc   = y /delta_resc
    
    !---------------------------------------------------------------------!
    ! compute inner/outer weights
    alpha = 4d0
    b     = 0.2d0
    Weig  = 0.5d0*(1d0+dtanh(alpha*(etag_inlet-b)/((1d0-2d0*b)*etag_inlet+b))/dtanh(alpha))
    Wei   = 0.5d0*(1d0+dtanh(alpha*(eta_inlet -b)/((1d0-2d0*b)*eta_inlet +b))/dtanh(alpha))

    !---------------------------------------------------------------------!        
    ! compute inner and outer components 

    ! U mean
    Call interp_vel(Umean_resc_T, Umean_inner,  ygp_resc,  ygp_inlet)
    Call interp_vel(Umean_resc_T, Umean_outer, etag_resc, etag_inlet)
    Umean_inner = gamma*Umean_inner
    Umean_outer = gamma*Umean_outer + (1d0-gamma)*Uinf
    
    ! V mean
    Call interp_vel(Vmean_resc_T, Vmean_inner,  yp_resc,  yp_inlet)
    Call interp_vel(Vmean_resc_T, Vmean_outer, eta_resc, eta_inlet)
    
    ! U fluctuations
    Call interp_vel(Ufluc_resc, Ufluc_inner,  ygp_resc,  ygp_inlet)
    Call interp_vel(Ufluc_resc, Ufluc_outer, etag_resc, etag_inlet)
    Ufluc_inner = gamma*Ufluc_inner
    Ufluc_outer = gamma*Ufluc_outer
    
    ! V fluctuations
    Call interp_vel(Vfluc_resc, Vfluc_inner,  yp_resc,  yp_inlet)
    Call interp_vel(Vfluc_resc, Vfluc_outer, eta_resc, eta_inlet)
    Vfluc_inner = gamma*Vfluc_inner
    Vfluc_outer = gamma*Vfluc_outer
    
    ! W fluctuations
    Call interp_vel(Wfluc_resc, Wfluc_inner,  ygp_resc,  ygp_inlet)
    Call interp_vel(Wfluc_resc, Wfluc_outer, etag_resc, etag_inlet)
    Wfluc_inner = gamma*Wfluc_inner
    Wfluc_outer = gamma*Wfluc_outer    

    !---------------------------------------------------------------------!        
    ! compute the inlet flow field 
    If ( iflag_fluc==1 ) Then
       Umean_inner = 0d0
       Umean_outer = 0d0
       Vmean_inner = 0d0
       Vmean_outer = 0d0
    End If
    Do j = 1,nyg
       U_inlet(j,:) = (Umean_inner(j,1)+Ufluc_inner(j,:))*(1d0-Weig(j)) + (Umean_outer(j,1)+Ufluc_outer(j,:))*Weig(j)
       W_inlet(j,:) = (                 Wfluc_inner(j,:))*(1d0-Weig(j)) + (                 Wfluc_outer(j,:))*Weig(j)
    End Do
    Do j = 1,ny
       V_inlet(j,:) = (Vmean_inner(j,1)+Vfluc_inner(j,:))*(1d0-Wei (j)) + (Vmean_outer(j,1)+Vfluc_outer(j,:))*Wei (j)
    End Do
    
  End Subroutine compute_rescaled_inflow

  !---------------------------------------------------------!
  !                Interpolate velocity in y                !
  !---------------------------------------------------------!
  Subroutine interp_vel(U_,U_interp,y_,y_interp)

    Real(Int64), Dimension(:,:), Intent(In)  :: U_
    Real(Int64), Dimension(:)  , Intent(In)  :: y_
    Real(Int64), Dimension(:)  , Intent(In)  :: y_interp
    Real(Int64), Dimension(:,:), Intent(Out) :: U_interp
    
    ! local variables
    Integer(Int32) :: n(2), ni(2), n1, ni1, i, j, k
    Real   (Int64) :: coef

    ! get sizes
    n   = Shape(U_)
    n1  = n(1)
    ni  = Shape(U_interp)
    ni1 = ni(1)

    If ( n1/=ni1 ) Stop 'Error! ni1/=ni'

    ! default value is top bc
    Do i = 2, ni1
       U_interp(i,:) = U_(n1,:)
    End Do
    
    ! interpolate
    Do i = 2, ni1 ! i=1 is ghost cell
       j = 0
       Do k = 1, n1
          If ( y_(k) > y_interp(i) ) Then
             j = k
             Exit
          End If
       End Do   
       If ( j==1 ) Stop 'Error! j==1'
       If ( j >0 ) Then 
          coef          = ( y_(j)-y_interp(i) )/( y_(j)-y_(j-1) ) 
          U_interp(i,:) = coef*U_(j-1,:) + (1d0-coef)*U_(j,:)
       End If
    End Do
       
  End Subroutine interp_vel

End Module Lund_rescaled_bc
