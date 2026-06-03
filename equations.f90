!--------------------------------------------------------!
! Module to compute right-hand side of Navier-Stokes eq. !
!--------------------------------------------------------!
Module equations

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : x, xm, xg, y, ym, yg, z, zm, zg, term_1, &
                              term_2, term, nx, nxg, ny, nyg, nz, nzg, &
                              nu, dPdx, dPdz, yg_m, nu_t, in1, in2,    &
                              weight_y_0, weight_y_1, Omega_z
  Use interpolation

  ! prevent implicit typing
  Implicit None

Contains

  !--------------------------------------------------------------------!
  !                      Compute RHS for du/dt                         !
  !                                                                    !
  ! du/dt = -du^2/dx - duv/dy - duw/dz + div(nu grad(u)) + 2*Omega_z*v !
  !--------------------------------------------------------------------!
  Subroutine compute_rhs_u(U_,V_,W_,rhs_u)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nx-1,2:nyg-1,2:nzg-1), Intent(Out) :: rhs_u

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2

    !-------------compute convective terms----------------!

    ! 1)---------compute -du^2/dx--------------

    ! interpolate u in x (faces to centers)
    Call interpolate_x(U_,term_1,in1)

    ! -du^2/dx (squaring fused into derivative)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             rhs_u(i,j,k) = -(term_1(i,j,k)**2-term_1(i-1,j,k)**2)/( xg(i+1) - xg(i) )
          End Do
       End Do
    End Do

    ! 2)-----------compute -duv/dy--------------

    ! interpolate u in y (centers to faces)
    Call interpolate_y(U_,term_1,in2)

    ! interpolate v in x (centers to faces)
    Call interpolate_x(V_,term_2,in2)

    ! -duv/dy (multiply fused into derivative, accumulate into rhs_u)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             rhs_u(i,j,k) = rhs_u(i,j,k) &
                - ( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j-1,k)*term_2(i,j-1,k) )/( y(j) - y(j-1) )
          End Do
       End Do
    End Do

    ! 3)--------------compute -duw/dz--------------

    ! interpolate u in z (centers to faces)
    Call interpolate_z(U_,term_1,in2)

    ! interpolate w in x (centers to faces)
    Call interpolate_x(W_,term_2,in2)

    ! -duw/dz (multiply fused, accumulate into rhs_u)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             rhs_u(i,j,k) = rhs_u(i,j,k) &
                - ( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j,k-1)*term_2(i,j,k-1) )/( z(k) - z(k-1) )
          End Do
       End Do
    End Do

    !--------------compute viscous terms------------!

    ! 4)--------compute div( nu grad(u) )------------

    ! first derivation in y: du/dy goes to yg_m(1:ny)
    !$acc parallel loop collapse(2) default(present)
    Do i = 1, nx
       Do k = 1, nzg
          term_1(i,1:nyg-1,k) = ( U_(i,2:nyg,k) - U_(i,1:nyg-1,k) )/( yg(2:nyg) - yg(1:nyg-1) )
       End Do
    End Do
    ! interpolate du/dy from yg_m(1:ny) to faces y(1:ny)
    !Call interpolate_y_2nd( yg_m, term_1, y, term_2 )
    !$acc kernels default(present)
    term_2 = term_1 ! -> uncomment for no interpolation
    !$acc end kernels

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dx_1,dx_2,dx_3,dy_3,dz_1,dz_2,dz_3,nu_x1,nu_x2,nu_y1,nu_y2,nu_z1,nu_z2)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1

             ! grid spacings
             dz_1 = zg(  k) - zg(k-1)
             dz_2 = zg(k+1) - zg(k  )
             dz_3 = z(k) - z(k-1)
             dy_3 = y(j) - y(j-1)
             dx_1 = x(  i) - x(i-1)
             dx_2 = x(i+1) - x(i  )
             dx_3 = xg(i+1) - xg(i)

             ! total viscosity at x locations
             nu_x1  = nu + nu_t(i  ,j,k)
             nu_x2  = nu + nu_t(i+1,j,k)

             ! total viscosity at y locations
             nu_y1 = nu + 0.5d0*( weight_y_0(j-1)*nu_t(i,  j-1,k) + weight_y_1(j-1)*nu_t(i,  j  ,k) + &
                                  weight_y_0(j-1)*nu_t(i+1,j-1,k) + weight_y_1(j-1)*nu_t(i+1,j  ,k) )
             nu_y2 = nu + 0.5d0*( weight_y_0(j  )*nu_t(i,  j,  k) + weight_y_1(j  )*nu_t(i,  j+1,k) + &
                                  weight_y_0(j  )*nu_t(i+1,j,  k) + weight_y_1(j  )*nu_t(i+1,j+1,k) )

             ! total viscosity at z locations
             nu_z1 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k-1)+nu_t(i+1,j,k)+nu_t(i+1,j,k-1))
             nu_z2 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i+1,j,k)+nu_t(i+1,j,k+1))

             ! viscous term (accumulated directly into rhs_u with pressure gradient)
             rhs_u(i,j,k) = rhs_u(i,j,k) + dPdx + &
                           1d0/dx_3*(nu_x2*1d0/dx_2*(U_(i+1,j,k)-U_(i,j,k)) - nu_x1*1d0/dx_1*(U_(i,j,k)-U_(i-1,j,k)) ) + &
                           1d0/dy_3*(nu_y2*term_2(i,j,k)                    - nu_y1*term_2(i,j-1,k) )                  + &
                           1d0/dz_3*(nu_z2*1d0/dz_2*(U_(i,j,k+1)-U_(i,j,k)) - nu_z1*1d0/dz_1*(U_(i,j,k)-U_(i,j,k-1)) )

       End Do
    End Do
    End Do

    !-------------------Rotating force--------------------!
    If ( Abs(Omega_z)>1d-6 ) Then

       ! interpolate v in y (faces to centers)
       call interpolate_y(V_,term,1)
       ! interpolate v in x (center to faces)
       call interpolate_x(term,term_1,2)
       !$acc kernels default(present)
       rhs_u = rhs_u + 2d0*Omega_z*term_1(2:nx-1,1:ny-1,2:nzg-1)
       !$acc end kernels

    End If

    !-----------------last plane is dummy----------------!
    !$acc parallel loop collapse(2) default(present)
    Do k=2,nzg-1
       Do j=2,nyg-1
          rhs_u(nx-1,j,k) = 0d0
       End Do
    End Do

  End Subroutine compute_rhs_u

  !--------------------------------------------------------------------!
  !                       Compute RHS for dv/dt                        !
  !                                                                    !
  ! dv/dt = -duv/dx - dv^2/dy - dvw/dz + div(nu grad(v)) - 2*Omega_z*u !
  !--------------------------------------------------------------------!
  Subroutine compute_rhs_v(U_,V_,W_,rhs_v)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nxg-1,2:ny-1,2:nzg-1), Intent(Out) :: rhs_v

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2

    !-------------compute convective terms----------------!

    ! 1)---------compute -dv^2/dy--------------

    ! interpolate v in y (faces to centers)
    Call interpolate_y(V_,term_1,in1)

    ! -dv^2/dy (squaring fused into derivative)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             rhs_v(i,j,k) = -(term_1(i,j,k)**2-term_1(i,j-1,k)**2)/( yg(j+1) - yg(j) )
          End Do
       End Do
    End Do

    ! 2)-----------compute -duv/dx--------------

    ! interpolate u in y (centers to faces)
    Call interpolate_y(U_,term_1,in2)

    ! interpolate v in x (centers to faces)
    Call interpolate_x(V_,term_2,in2)

    ! -duv/dx (multiply fused, accumulate into rhs_v)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             rhs_v(i,j,k) = rhs_v(i,j,k) &
                - ( term_1(i,j,k)*term_2(i,j,k) - term_1(i-1,j,k)*term_2(i-1,j,k) )/( x(i) - x(i-1) )
          End Do
       End Do
    End Do

    ! 3)--------------compute -dvw/dz--------------

    ! interpolate v in z (centers to faces)
    Call interpolate_z(V_,term_1,in2)

    ! interpolate w in y (centers to faces)
    Call interpolate_y(W_,term_2,in2)

    ! -dvw/dz (multiply fused, accumulate into rhs_v)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             rhs_v(i,j,k) = rhs_v(i,j,k) &
                - ( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j,k-1)*term_2(i,j,k-1) )/( z(k) - z(k-1) )
          End Do
       End Do
    End Do

    !--------------compute viscous terms------------!

    ! 4)-----------compute div( nu grad(v) )----------

    ! interpolate eddy viscosity to faces

    ! second order remain, no need to interpolate
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dx_1,dx_2,dx_3,dy_1,dy_2,dy_3,dz_1,dz_2,dz_3,nu_x1,nu_x2,nu_y1,nu_y2,nu_z1,nu_z2)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1

             ! grid spacings
             dz_1 = zg(  k) - zg(k-1)
             dz_2 = zg(k+1) - zg(k  )
             dz_3 = z(k) - z(k-1)
             dy_1 = y(  j) - y(j-1)
             dy_2 = y(j+1) - y(j  )
             dy_3 = yg(j+1) - yg(j)
             dx_1 = xg(  i) - xg(i-1)
             dx_2 = xg(i+1) - xg(i  )
             dx_3 = x(i) - x(i-1)

             ! eddy viscosity at x locations
             nu_x1 = nu + 0.5d0*( weight_y_0(j)*nu_t(i-1,j,k) + weight_y_1(j)*nu_t(i-1,j+1,k) + &
                                  weight_y_0(j)*nu_t(i  ,j,k) + weight_y_1(j)*nu_t(i  ,j+1,k) )
             nu_x2 = nu + 0.5d0*( weight_y_0(j)*nu_t(i,  j,k) + weight_y_1(j)*nu_t(i,  j+1,k) + &
                                  weight_y_0(j)*nu_t(i+1,j,k) + weight_y_1(j)*nu_t(i+1,j+1,k) )

             ! eddy viscosity at the centers
             nu_y1 = nu + nu_t(i,j  ,k)
             nu_y2 = nu + nu_t(i,j+1,k)

             ! eddy viscosity at z locations
             nu_z1 = nu + 0.5d0*( weight_y_0(j)*nu_t(i,j,k-1) + weight_y_1(j)*nu_t(i,j+1,k-1) + &
                                  weight_y_0(j)*nu_t(i,j,k  ) + weight_y_1(j)*nu_t(i,j+1,k  ) )
             nu_z2 = nu + 0.5d0*( weight_y_0(j)*nu_t(i,j,  k) + weight_y_1(j)*nu_t(i,j+1,  k) + &
                                  weight_y_0(j)*nu_t(i,j,k+1) + weight_y_1(j)*nu_t(i,j+1,k+1) )

             ! viscous term (accumulated directly into rhs_v)
             rhs_v(i,j,k) = rhs_v(i,j,k) + &
                           1d0/dx_3*(nu_x2*1d0/dx_2*(V_(i+1,j,k) - V_(i,j,k)) - nu_x1*1d0/dx_1*(V_(i,j,k)-V_(i-1,j,k)) ) + &
                           1d0/dy_3*(nu_y2*1d0/dy_2*(V_(i,j+1,k) - V_(i,j,k)) - nu_y1*1d0/dy_1*(V_(i,j,k)-V_(i,j-1,k)) ) + &
                           1d0/dz_3*(nu_z2*1d0/dz_2*(V_(i,j,k+1) - V_(i,j,k)) - nu_z1*1d0/dz_1*(V_(i,j,k)-V_(i,j,k-1)) )

       End Do
    End Do
    End Do

    !-------------------Rotating force--------------------!
    If ( Abs(Omega_z)>1d-6 ) Then

       ! interpolate u in x (faces to centers)
       call interpolate_x(U_,term,1)
       ! interpolate u in y (center to faces)
       call interpolate_y(term,term_1,2)
       !$acc kernels default(present)
       rhs_v = rhs_v - 2d0*Omega_z*term_1(1:nx-1,2:ny-1,2:nzg-1)
       !$acc end kernels

    End If

    !-----------------last plane is dummy----------------!
    !$acc parallel loop collapse(2) default(present)
    Do k=2,nzg-1
       Do j=2,ny-1
          rhs_v(nxg-1,j,k) = 0d0
       End Do
    End Do

  End Subroutine compute_rhs_v

  !-------------------------------------------------------!
  !                Compute RHS for dw/dt                  !
  !                                                       !
  !  dw/dt = -duw/dx - dvw/dy - dw^2/dz + div(nu grad(w)) !
  !-------------------------------------------------------!
  Subroutine compute_rhs_w(U_,V_,W_,rhs_w)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nxg-1,2:nyg-1,2:nz-1), Intent(Out) :: rhs_w

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr, w0, w1, yd
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2

    !-------------compute convective terms---------------!

    ! 1)---------compute -dw^2/dz--------------

    ! interpolate w in z (faces to centers)
    Call interpolate_z(W_,term_1,in1)

    ! -dw^2/dz (squaring fused into derivative)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             rhs_w(i,j,k) = -(term_1(i,j,k)**2-term_1(i,j,k-1)**2)/( zg(k+1) - zg(k) )
          End Do
       End Do
    End Do

    ! 2)-----------compute -duw/dx--------------

    ! interpolate u in z (centers to faces)
    Call interpolate_z(U_,term_1,in2)

    ! interpolate w in x (centers to faces)
    Call interpolate_x(W_,term_2,in2)

    ! -duw/dx (multiply fused, accumulate into rhs_w)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             rhs_w(i,j,k) = rhs_w(i,j,k) &
                - ( term_1(i,j,k)*term_2(i,j,k) - term_1(i-1,j,k)*term_2(i-1,j,k) )/( x(i) - x(i-1) )
          End Do
       End Do
    End Do

    ! 3)--------------compute -dvw/dy--------------

    ! interpolate v in z (centers to faces)
    Call interpolate_z(V_,term_1,in2)

    ! interpolate w in y (centers to faces)
    Call interpolate_y(W_,term_2,in2)

    ! -dvw/dy (multiply fused, accumulate into rhs_w)
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             rhs_w(i,j,k) = rhs_w(i,j,k) &
                - ( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j-1,k)*term_2(i,j-1,k) )/( y(j) - y(j-1) )
          End Do
       End Do
    End Do

    !--------------compute viscous terms------------!

    ! 4)----------compute div( nu grad(w) )-----------
    ! first derivation in y: dw/dy -> term_2 (need to preserve before term_1 is reused)
    !$acc parallel loop collapse(3) default(present)
    Do k = 1, nz
       Do j = 1, nyg-1
          Do i = 1, nxg
             term_2(i,j,k) = ( W_(i,j+1,k) - W_(i,j,k) )/( yg(j+1) - yg(j) )
          End Do
       End Do
    End Do

    ! interpolate eddy viscosity to faces
    Call interpolate_y(nu_t(2:nxg-1,1:nyg,2:nzg-1),term_1(2:nxg-1,1:ny,2:nzg-1),in2)

    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dx_1,dx_2,dx_3,dy_3,dz_1,dz_2,dz_3,nu_x1,nu_x2,nu_y1,nu_y2,nu_z1,nu_z2)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1

             ! grid spacings
             dz_1 = z(  k) - z(k-1)
             dz_2 = z(k+1) - z(k  )
             dz_3 = zg(k+1) - zg(k)
             dy_3 = y(j) - y(j-1)
             dx_1 = xg(  i) - xg(i-1)
             dx_2 = xg(i+1) - xg(i  )
             dx_3 = x(i) - x(i-1)

             ! eddy viscosity at x locations
             nu_x1 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i-1,j,k)+nu_t(i-1,j,k+1))
             nu_x2 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i+1,j,k)+nu_t(i+1,j,k+1))

             ! eddy viscosity at y locations
             nu_y1 = nu + 0.5d0*( weight_y_0(j-1)*nu_t(i,j-1,  k) + weight_y_1(j-1)*nu_t(i,j  ,  k) + &
                                  weight_y_0(j-1)*nu_t(i,j-1,k+1) + weight_y_1(j-1)*nu_t(i,j  ,k+1) )
             nu_y2 = nu + 0.5d0*( weight_y_0(j  )*nu_t(i,j  ,  k) + weight_y_1(j  )*nu_t(i,j+1,  k) + &
                                  weight_y_0(j  )*nu_t(i,j  ,k+1) + weight_y_1(j  )*nu_t(i,j+1,k+1) )

             ! eddy viscosity at z locations
             nu_z1 = nu + nu_t(i,j,k)
             nu_z2 = nu + nu_t(i,j,k+1)

             ! viscous term (accumulated directly into rhs_w with pressure gradient)
             rhs_w(i,j,k) = rhs_w(i,j,k) + dPdz + &
                           1d0/dx_3*(nu_x2*1d0/dx_2*(W_(i+1,j,k) - W_(i,j,k)) - nu_x1*1d0/dx_1*(W_(i,j,k) - W_(i-1,j,k)) ) + &
                           1d0/dy_3*(nu_y2*term_2(i,j,k)                      - nu_y1*term_2(i,j-1,k) )                    + &
                           1d0/dz_3*(nu_z2*1d0/dz_2*(W_(i,j,k+1) - W_(i,j,k)) - nu_z1*1d0/dz_1*(W_(i,j,k) - W_(i,j,k-1)) )

       End Do
    End Do
    End Do

    !-----------------last plane is dummy----------------!
    !$acc parallel loop collapse(2) default(present)
    Do k=2,nz-1
       Do j=2,nyg-1
          rhs_w(nxg-1,j,k) = 0d0
       End Do
    End Do

  End Subroutine compute_rhs_w

End Module equations
