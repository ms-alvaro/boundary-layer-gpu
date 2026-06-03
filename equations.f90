!--------------------------------------------------------!
! Module to compute right-hand side of Navier-Stokes eq. !
!--------------------------------------------------------!
Module equations

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : x, xm, xg, y, ym, yg, z, zm, zg, term_1, &
                              term_2, term, nx, nxg, ny, nyg, nz, nzg, &
                              nu, dPdx, dPdz, yg_m, nu_t, in1, in2,    &
                              weight_y_0, weight_y_1, Omega_z, LES_model
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

    ! All convective terms for u, fully inlined (no intermediate arrays)
    ! 1) -du^2/dx: interp_x(U,centers) then d/dx
    !    u_c(i) = 0.5*(U(i)+U(i+1)), du^2/dx = (u_c(i)^2 - u_c(i-1)^2)/dx
    ! 2) -duv/dy: interp_y(U,faces)*interp_x(V,faces) then d/dy
    !    u_f(j) = w0*U(j)+w1*U(j+1), v_f(i) = 0.5*(V(i)+V(i+1))
    ! 3) -duw/dz: interp_z(U,faces)*interp_x(W,faces) then d/dz
    !    u_f(k) = 0.5*(U(k)+U(k+1)), w_f(i) = 0.5*(W(i)+W(i+1))
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(dx_1)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             ! 1) -du^2/dx (interp_x faces->centers = simple average)
             rhs_u(i,j,k) = -( (0.5d0*(U_(i,j,k)+U_(i+1,j,k)))**2 &
                              - (0.5d0*(U_(i-1,j,k)+U_(i,j,k)))**2 ) / (xg(i+1)-xg(i))
             ! 2) -duv/dy (interp_y centers->faces = weighted, interp_x centers->faces)
             rhs_u(i,j,k) = rhs_u(i,j,k) &
                - ( (weight_y_0(j)*U_(i,j,k)+weight_y_1(j)*U_(i,j+1,k)) &
                    * 0.5d0*(V_(i,j,k)+V_(i+1,j,k)) &
                  - (weight_y_0(j-1)*U_(i,j-1,k)+weight_y_1(j-1)*U_(i,j,k)) &
                    * 0.5d0*(V_(i,j-1,k)+V_(i+1,j-1,k)) ) / (y(j)-y(j-1))
             ! 3) -duw/dz (interp_z centers->faces = simple avg, interp_x centers->faces)
             rhs_u(i,j,k) = rhs_u(i,j,k) &
                - ( 0.5d0*(U_(i,j,k)+U_(i,j,k+1)) * 0.5d0*(W_(i,j,k)+W_(i+1,j,k)) &
                  - 0.5d0*(U_(i,j,k-1)+U_(i,j,k)) * 0.5d0*(W_(i,j,k-1)+W_(i+1,j,k-1)) ) / (z(k)-z(k-1))
          End Do
       End Do
    End Do

    !--------------compute viscous terms (fused into convective kernel for DNS)-----!

    If ( LES_model == 0 ) Then
       ! DNS: nu_t = 0, viscous = nu * Laplacian(U) — inline into single kernel
       ! Rewrite the convective kernel to include viscous + pressure gradient
       ! Note: we OVERWRITE rhs_u here (convective already computed above)
       ! Actually we need to ADD viscous to the existing rhs_u from convective
       !$acc parallel loop collapse(3) default(present) &
       !$acc private(dx_1,dx_2,dx_3,dy_3,dz_1,dz_2,dz_3)
       Do k=2,nzg-1
          Do j=2,nyg-1
             Do i=2,nx-1
                dx_1 = x(i) - x(i-1)
                dx_2 = x(i+1) - x(i)
                dx_3 = xg(i+1) - xg(i)
                dy_3 = y(j) - y(j-1)
                dz_1 = zg(k) - zg(k-1)
                dz_2 = zg(k+1) - zg(k)
                dz_3 = z(k) - z(k-1)
                rhs_u(i,j,k) = rhs_u(i,j,k) + dPdx + nu * ( &
                   (1d0/dx_2*(U_(i+1,j,k)-U_(i,j,k)) - 1d0/dx_1*(U_(i,j,k)-U_(i-1,j,k)))/dx_3 + &
                   ((U_(i,j+1,k)-U_(i,j,k))/(yg(j+1)-yg(j)) - (U_(i,j,k)-U_(i,j-1,k))/(yg(j)-yg(j-1)))/dy_3 + &
                   (1d0/dz_2*(U_(i,j,k+1)-U_(i,j,k)) - 1d0/dz_1*(U_(i,j,k)-U_(i,j,k-1)))/dz_3 )
             End Do
          End Do
       End Do
    Else
       ! LES: variable nu_t, need interpolated du/dy
       !$acc parallel loop collapse(3) default(present)
       Do k = 1, nzg
          Do j = 1, nyg-1
             Do i = 1, nx
                term_2(i,j,k) = ( U_(i,j+1,k) - U_(i,j,k) )/( yg(j+1) - yg(j) )
             End Do
          End Do
       End Do

       !$acc parallel loop collapse(3) default(present) &
       !$acc private(dx_1,dx_2,dx_3,dy_3,dz_1,dz_2,dz_3,nu_x1,nu_x2,nu_y1,nu_y2,nu_z1,nu_z2)
       Do k=2,nzg-1
          Do j=2,nyg-1
             Do i=2,nx-1
                dz_1 = zg(k) - zg(k-1)
                dz_2 = zg(k+1) - zg(k)
                dz_3 = z(k) - z(k-1)
                dy_3 = y(j) - y(j-1)
                dx_1 = x(i) - x(i-1)
                dx_2 = x(i+1) - x(i)
                dx_3 = xg(i+1) - xg(i)
                nu_x1 = nu + nu_t(i,j,k)
                nu_x2 = nu + nu_t(i+1,j,k)
                nu_y1 = nu + 0.5d0*(weight_y_0(j-1)*nu_t(i,j-1,k) + weight_y_1(j-1)*nu_t(i,j,k) + &
                                     weight_y_0(j-1)*nu_t(i+1,j-1,k) + weight_y_1(j-1)*nu_t(i+1,j,k))
                nu_y2 = nu + 0.5d0*(weight_y_0(j)*nu_t(i,j,k) + weight_y_1(j)*nu_t(i,j+1,k) + &
                                     weight_y_0(j)*nu_t(i+1,j,k) + weight_y_1(j)*nu_t(i+1,j+1,k))
                nu_z1 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k-1)+nu_t(i+1,j,k)+nu_t(i+1,j,k-1))
                nu_z2 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i+1,j,k)+nu_t(i+1,j,k+1))
                rhs_u(i,j,k) = rhs_u(i,j,k) + dPdx + &
                   (nu_x2/dx_2*(U_(i+1,j,k)-U_(i,j,k)) - nu_x1/dx_1*(U_(i,j,k)-U_(i-1,j,k)))/dx_3 + &
                   (nu_y2*term_2(i,j,k) - nu_y1*term_2(i,j-1,k))/dy_3 + &
                   (nu_z2/dz_2*(U_(i,j,k+1)-U_(i,j,k)) - nu_z1/dz_1*(U_(i,j,k)-U_(i,j,k-1)))/dz_3
             End Do
          End Do
       End Do
    End If

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

    ! All convective terms for v, fully inlined
    ! 1) -dv^2/dy: interp_y(V,centers=0.5avg) then d/dy
    ! 2) -duv/dx: interp_y(U,faces=weighted)*interp_x(V,faces=0.5avg) then d/dx
    ! 3) -dvw/dz: interp_z(V,faces=0.5avg)*interp_y(W,faces=weighted) then d/dz
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             ! 1) -dv^2/dy (interp_y faces->centers = 0.5 average)
             rhs_v(i,j,k) = -( (0.5d0*(V_(i,j,k)+V_(i,j+1,k)))**2 &
                              - (0.5d0*(V_(i,j-1,k)+V_(i,j,k)))**2 ) / (yg(j+1)-yg(j))
             ! 2) -duv/dx (interp_y centers->faces=weighted, interp_x centers->faces=0.5avg)
             rhs_v(i,j,k) = rhs_v(i,j,k) &
                - ( (weight_y_0(j)*U_(i,j,k)+weight_y_1(j)*U_(i,j+1,k)) &
                    * 0.5d0*(V_(i,j,k)+V_(i+1,j,k)) &
                  - (weight_y_0(j)*U_(i-1,j,k)+weight_y_1(j)*U_(i-1,j+1,k)) &
                    * 0.5d0*(V_(i-1,j,k)+V_(i,j,k)) ) / (x(i)-x(i-1))
             ! 3) -dvw/dz (interp_z centers->faces=0.5avg, interp_y centers->faces=weighted)
             rhs_v(i,j,k) = rhs_v(i,j,k) &
                - ( 0.5d0*(V_(i,j,k)+V_(i,j,k+1)) &
                    * (weight_y_0(j)*W_(i,j,k)+weight_y_1(j)*W_(i,j+1,k)) &
                  - 0.5d0*(V_(i,j,k-1)+V_(i,j,k)) &
                    * (weight_y_0(j)*W_(i,j,k-1)+weight_y_1(j)*W_(i,j+1,k-1)) ) / (z(k)-z(k-1))
          End Do
       End Do
    End Do

    !--------------compute viscous terms------------!

    If ( LES_model == 0 ) Then
       ! DNS: nu * Laplacian(V) inlined
       !$acc parallel loop collapse(3) default(present) &
       !$acc private(dx_1,dx_2,dx_3,dy_1,dy_2,dy_3,dz_1,dz_2,dz_3)
       Do k=2,nzg-1
          Do j=2,ny-1
             Do i=2,nxg-1
                dx_1 = xg(i) - xg(i-1)
                dx_2 = xg(i+1) - xg(i)
                dx_3 = x(i) - x(i-1)
                dy_1 = y(j) - y(j-1)
                dy_2 = y(j+1) - y(j)
                dy_3 = yg(j+1) - yg(j)
                dz_1 = zg(k) - zg(k-1)
                dz_2 = zg(k+1) - zg(k)
                dz_3 = z(k) - z(k-1)
                rhs_v(i,j,k) = rhs_v(i,j,k) + nu * ( &
                   (1d0/dx_2*(V_(i+1,j,k)-V_(i,j,k)) - 1d0/dx_1*(V_(i,j,k)-V_(i-1,j,k)))/dx_3 + &
                   (1d0/dy_2*(V_(i,j+1,k)-V_(i,j,k)) - 1d0/dy_1*(V_(i,j,k)-V_(i,j-1,k)))/dy_3 + &
                   (1d0/dz_2*(V_(i,j,k+1)-V_(i,j,k)) - 1d0/dz_1*(V_(i,j,k)-V_(i,j,k-1)))/dz_3 )
             End Do
          End Do
       End Do
    Else
       ! LES: variable nu_t
       !$acc parallel loop collapse(3) default(present) &
       !$acc private(dx_1,dx_2,dx_3,dy_1,dy_2,dy_3,dz_1,dz_2,dz_3,nu_x1,nu_x2,nu_y1,nu_y2,nu_z1,nu_z2)
       Do k=2,nzg-1
          Do j=2,ny-1
             Do i=2,nxg-1
                dz_1 = zg(k) - zg(k-1)
                dz_2 = zg(k+1) - zg(k)
                dz_3 = z(k) - z(k-1)
                dy_1 = y(j) - y(j-1)
                dy_2 = y(j+1) - y(j)
                dy_3 = yg(j+1) - yg(j)
                dx_1 = xg(i) - xg(i-1)
                dx_2 = xg(i+1) - xg(i)
                dx_3 = x(i) - x(i-1)
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
    End If

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

    ! All convective terms for w, fully inlined
    ! 1) -dw^2/dz: interp_z(W,centers=0.5avg) then d/dz
    ! 2) -duw/dx: interp_z(U,faces=0.5avg)*interp_x(W,faces=0.5avg) then d/dx
    ! 3) -dvw/dy: interp_z(V,faces=0.5avg)*interp_y(W,faces=weighted) then d/dy
    !$acc parallel loop collapse(3) default(present)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             ! 1) -dw^2/dz (interp_z faces->centers = 0.5 average)
             rhs_w(i,j,k) = -( (0.5d0*(W_(i,j,k)+W_(i,j,k+1)))**2 &
                              - (0.5d0*(W_(i,j,k-1)+W_(i,j,k)))**2 ) / (zg(k+1)-zg(k))
             ! 2) -duw/dx (interp_z centers->faces=0.5avg, interp_x centers->faces=0.5avg)
             rhs_w(i,j,k) = rhs_w(i,j,k) &
                - ( 0.5d0*(U_(i,j,k)+U_(i,j,k+1)) * 0.5d0*(W_(i,j,k)+W_(i+1,j,k)) &
                  - 0.5d0*(U_(i-1,j,k)+U_(i-1,j,k+1)) * 0.5d0*(W_(i-1,j,k)+W_(i,j,k)) ) / (x(i)-x(i-1))
             ! 3) -dvw/dy (interp_z centers->faces=0.5avg, interp_y centers->faces=weighted)
             rhs_w(i,j,k) = rhs_w(i,j,k) &
                - ( 0.5d0*(V_(i,j,k)+V_(i,j,k+1)) &
                    * (weight_y_0(j)*W_(i,j,k)+weight_y_1(j)*W_(i,j+1,k)) &
                  - 0.5d0*(V_(i,j-1,k)+V_(i,j-1,k+1)) &
                    * (weight_y_0(j-1)*W_(i,j-1,k)+weight_y_1(j-1)*W_(i,j,k)) ) / (y(j)-y(j-1))
          End Do
       End Do
    End Do

    !--------------compute viscous terms------------!

    If ( LES_model == 0 ) Then
       ! DNS: nu * Laplacian(W) + dPdz inlined
       !$acc parallel loop collapse(3) default(present) &
       !$acc private(dx_1,dx_2,dx_3,dy_3,dz_1,dz_2,dz_3)
       Do k=2,nz-1
          Do j=2,nyg-1
             Do i=2,nxg-1
                dx_1 = xg(i) - xg(i-1)
                dx_2 = xg(i+1) - xg(i)
                dx_3 = x(i) - x(i-1)
                dy_3 = y(j) - y(j-1)
                dz_1 = z(k) - z(k-1)
                dz_2 = z(k+1) - z(k)
                dz_3 = zg(k+1) - zg(k)
                rhs_w(i,j,k) = rhs_w(i,j,k) + dPdz + nu * ( &
                   (1d0/dx_2*(W_(i+1,j,k)-W_(i,j,k)) - 1d0/dx_1*(W_(i,j,k)-W_(i-1,j,k)))/dx_3 + &
                   ((W_(i,j+1,k)-W_(i,j,k))/(yg(j+1)-yg(j)) - (W_(i,j,k)-W_(i,j-1,k))/(yg(j)-yg(j-1)))/dy_3 + &
                   (1d0/dz_2*(W_(i,j,k+1)-W_(i,j,k)) - 1d0/dz_1*(W_(i,j,k)-W_(i,j,k-1)))/dz_3 )
             End Do
          End Do
       End Do
    Else
       ! LES: variable nu_t
       !$acc parallel loop collapse(3) default(present)
       Do k = 1, nz
          Do j = 1, nyg-1
             Do i = 1, nxg
                term_2(i,j,k) = ( W_(i,j+1,k) - W_(i,j,k) )/( yg(j+1) - yg(j) )
             End Do
          End Do
       End Do
       Call interpolate_y(nu_t(2:nxg-1,1:nyg,2:nzg-1),term_1(2:nxg-1,1:ny,2:nzg-1),in2)
       !$acc parallel loop collapse(3) default(present) &
       !$acc private(dx_1,dx_2,dx_3,dy_3,dz_1,dz_2,dz_3,nu_x1,nu_x2,nu_y1,nu_y2,nu_z1,nu_z2)
       Do k=2,nz-1
          Do j=2,nyg-1
             Do i=2,nxg-1
                dz_1 = z(k) - z(k-1)
                dz_2 = z(k+1) - z(k)
                dz_3 = zg(k+1) - zg(k)
                dy_3 = y(j) - y(j-1)
                dx_1 = xg(i) - xg(i-1)
                dx_2 = xg(i+1) - xg(i)
                dx_3 = x(i) - x(i-1)
                nu_x1 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i-1,j,k)+nu_t(i-1,j,k+1))
                nu_x2 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i+1,j,k)+nu_t(i+1,j,k+1))
                nu_y1 = nu + 0.5d0*(weight_y_0(j-1)*nu_t(i,j-1,k)+weight_y_1(j-1)*nu_t(i,j,k) + &
                                     weight_y_0(j-1)*nu_t(i,j-1,k+1)+weight_y_1(j-1)*nu_t(i,j,k+1))
                nu_y2 = nu + 0.5d0*(weight_y_0(j)*nu_t(i,j,k)+weight_y_1(j)*nu_t(i,j+1,k) + &
                                     weight_y_0(j)*nu_t(i,j,k+1)+weight_y_1(j)*nu_t(i,j+1,k+1))
                nu_z1 = nu + nu_t(i,j,k)
                nu_z2 = nu + nu_t(i,j,k+1)
                rhs_w(i,j,k) = rhs_w(i,j,k) + dPdz + &
                   (nu_x2/dx_2*(W_(i+1,j,k)-W_(i,j,k)) - nu_x1/dx_1*(W_(i,j,k)-W_(i-1,j,k)))/dx_3 + &
                   (nu_y2*term_2(i,j,k) - nu_y1*term_2(i,j-1,k))/dy_3 + &
                   (nu_z2/dz_2*(W_(i,j,k+1)-W_(i,j,k)) - nu_z1/dz_1*(W_(i,j,k)-W_(i,j,k-1)))/dz_3
             End Do
          End Do
       End Do
    End If

    !-----------------last plane is dummy----------------!
    !$acc parallel loop collapse(2) default(present)
    Do k=2,nz-1
       Do j=2,nyg-1
          rhs_w(nxg-1,j,k) = 0d0
       End Do
    End Do

  End Subroutine compute_rhs_w

End Module equations
