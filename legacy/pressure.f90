!------------------------------------------------!
!      Module for computing actual pressure      !
!------------------------------------------------!
Module pressure

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use interpolation
  Use equations
  Use boundary_conditions

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------------!
  !              Compute actual pressure               !
  !----------------------------------------------------!
  Subroutine compute_pressure

    Call compute_pressure_BC

    Call compute_pressure_rhs

    Call solve_poisson_equation_for_pressure

  End Subroutine compute_pressure

  !-------------------------------------------------------!
  !     Compute right-hand size for pressure equation     !
  !                                                       !
  !   Equation:                                           !
  !    dvel/dt      = - grad(P) + rhs_vel                 !
  !    Laplacian(p) = rhs_p = div(rhs_vel)                !
  !                                                       !
  ! Input:  U, V, W                                       !
  ! Output: rhs_p                                         !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine compute_pressure_rhs

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: maxerr

    ! compute rhs terms
    Call compute_rhs_u( U, V, W, rhs_uf(2:nx-1, 2:nyg-1,2:nzg-1) )
    Call compute_rhs_v( U, V, W, rhs_vf(2:nxg-1,2:ny-1, 2:nzg-1) )
    Call compute_rhs_w( U, V, W, rhs_wf(2:nxg-1,2:nyg-1,2:nz-1 ) )

    ! compute rhs at the boundaries
    Call apply_periodic_bc_x   ( rhs_uf, 1 )
    Call compute_rhs_v_boundary( rhs_vf    )
    Call apply_periodic_bc_z   ( rhs_wf, 3 )

    ! compute divergence of rhs terms
    !$acc kernels default(present)
    rhs_p = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) default(present)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             rhs_p(i,j,k)  = ( rhs_uf(i,j,k) - rhs_uf(i-1,j,k) )/( x(i)-x(i-1) ) + & ! d rhs_u/dx
                             ( rhs_vf(i,j,k) - rhs_vf(i,j-1,k) )/( y(j)-y(j-1) ) + & ! d rhs_v/dy
                             ( rhs_wf(i,j,k) - rhs_wf(i,j,k-1) )/( z(k)-z(k-1) )     ! d rhs_w/dz
          End Do
       End Do
    End Do

  End Subroutine compute_pressure_rhs

  !----------------------------------------------------!
  !              Solve pressure equation               !
  !----------------------------------------------------!
  Subroutine solve_poisson_equation_for_pressure

    Integer(Int32) :: i, j, k, i_global, k_global, info
    Real   (Int64) :: maxerr

    External :: zgesv

    ! Transfer data from GPU to CPU for FFT operations
    !$acc update self(rhs_p, V, bc_1, bc_2)

    ! 2D Fourier transform interior points rhs_p
    ! NOTE: use plane as both in/out (in-place) to avoid nvfortran pointer alias bug
    Do j = 2, nyg-1
       plane = dcmplx( rhs_p(2:nxp+1,j,2:nzp+1) ) ! nxp+1 = nxg-2, nzp+1 = nzg-2
       Call fftw_mpi_execute_dft(plan_d,plane,plane)
       rhs_p_hat(:,j,:) = plane(1:mx+1, 1:mz+1)
    End Do

    ! 2D Fourier transform boundary conditions
    plane = dcmplx( bc_1(2:nxp+1,2:nzp+1) )
    Call fftw_mpi_execute_dft(plan_d,plane,plane)
    bc_1_hat = plane(1:mx+1, 1:mz+1)

    plane = dcmplx( bc_2(2:nxp+1,2:nzp+1) )
    Call fftw_mpi_execute_dft(plan_d,plane,plane)
    bc_2_hat = plane(1:mx+1, 1:mz+1)

    ! solve for each mode
    Do k = 0, mz
       Do i = 0, mx
          ! mapping to x-mode
          i_global = imode_map_fft( i, k )
          ! mapping to z-mode
          k_global = kmode_map_fft( i, k )
          ! form matrix
          Do j = 2, nyg-1 ! diagonal
             D(j) = Dyy(j,j) + kxx(i_global) + kzz(k_global)
          End Do
          Do j = 2, nyg-2 ! lower diagonal
             DL(j) = Dyy(j+1,j)
          End Do
          Do j = 2, nyg-2 ! upper diagonal
             DU(j) = Dyy(j,j+1)
          End Do
          ! remove singularity 00 mode (set a reference pressure)
          if ( i_global==0 .And. k_global==0 ) D(2) = 3d0/2d0*D(2)
          ! rhs with boundary conditions for pressure
          rhs_aux        = rhs_p_hat(i,:,k)
          rhs_aux(    2) = rhs_aux(    2) + coef_bc_1*(yg(  2)-yg(    1))*bc_1_hat(i,k)
          rhs_aux(nyg-1) = rhs_aux(nyg-1) - coef_bc_2*(yg(nyg)-yg(nyg-1))*bc_2_hat(i,k)
          ! solve M*u = rhs (solution stored in rhs_p_hat)
          Call Zgtsv( nr, nrhs, DL, D, DU, rhs_aux, nr, info)
          rhs_p_hat(i,:,k) = rhs_aux
       End Do
    End Do

    ! 2D inverse Fourier transform
    ! NOTE: use plane as both in/out (in-place) to avoid nvfortran pointer alias bug
    Do j = 2, nyg-1
       plane = (0d0,0d0)
       plane(1:mx+1, 1:mz+1) = rhs_p_hat(:,j,:)
       Call fftw_mpi_execute_dft(plan_i,plane,plane)
       P(2:nxg-2,j,2:nzg-2) = plane/Real(nxp_global*nzp_global,8)
    End Do

    ! Transfer P back to GPU
    !$acc update device(P)

    ! extend values to ghost cells in y
    !$acc kernels default(present)
    P(:,  1,:) = P(:,    2,:) - ( yg(  2)-yg(    1) )*bc_1
    P(:,nyg,:) = P(:,nyg-1,:) + ( yg(nyg)-yg(nyg-1) )*bc_2
    !$acc end kernels

    ! apply periodicity in x and z
    Call apply_periodic_xz_pressure

    ! update ghost interior planes
    Call update_ghost_interior_planes_pressure

  End Subroutine solve_poisson_equation_for_pressure

  !--------------------------------------------------!
  ! Periodic boundary conditions for pseudo-pressure !
  !--------------------------------------------------!
  Subroutine apply_periodic_xz_pressure

    ! apply periodicity in x (All processors, no MPI needed)
    !$acc kernels default(present)
    P ( nxg-1, :, : ) = P ( 2, :, : )
    !$acc end kernels

    ! apply periodicity in z (Only first and last processor, MPI needed)
    !$acc update self(P)
    If     ( myid==0 ) Then
       buffer_p = P ( 2:nxg-1, :, 2 )
       ! send data to nprocs-1
       Call Mpi_send(buffer_p, (nxg-2)*(nyg-2), MPI_real8, nprocs-1, 0, &
             MPI_COMM_WORLD,ierr)
    Elseif ( myid==nprocs-1 ) Then
       ! receive data from 0
       Call Mpi_recv(buffer_p, (nxg-2)*(nyg-2), MPI_real8, 0, 0, &
            MPI_COMM_WORLD,istat,ierr)
       P ( 2:nxg-1, :, nzp+1+1 ) = buffer_p
    End If
    !$acc update device(P)

  End Subroutine apply_periodic_xz_pressure

  !------------------------------------------------!
  !    Update ghost interior planes for pressure   !
  !------------------------------------------------!
  Subroutine update_ghost_interior_planes_pressure

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom

    !----------------------update P-----------------------!
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
    !$acc update self(P(:,:,2))
    buffer_ps = P(:,:,2)  ! send buffer
    Call Mpi_sendrecv(buffer_ps, (nxg-2)*(nyg-2), Mpi_real8, sendto, tagto,        &
         buffer_pr, (nxg-2)*(nyg-2), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
         istat, ierr)
    If ( myid/=nprocs-1 ) P(:,:,nzg) = buffer_pr ! received buffer
    !$acc update device(P(:,:,nzg))

  End Subroutine update_ghost_interior_planes_pressure

  !----------------------------------------------------------!
  !          Compute rhs for V at the y-boundaries           !
  !----------------------------------------------------------!
  Subroutine compute_rhs_v_boundary(rhs_v)

    Real(Int64), Dimension(1:nxg,1:ny,1:nzg), Intent(InOut) :: rhs_v

    !$acc kernels default(present)
    ! bottom wall
    rhs_v(2:nxg-1, 1,2:nzg-1) = bc_1

    ! top wall
    rhs_v(2:nxg-1,ny,2:nzg-1) = bc_2
    !$acc end kernels

  End Subroutine compute_rhs_v_boundary

  !----------------------------------------------------------!
  !     Compute boundary conditions for actual pressure      !
  !----------------------------------------------------------!
  Subroutine compute_pressure_BC

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: nui, dnu, dV, ddV, dzeta, g, g2, maxerr

    dzeta  = 1d0 ! arbitrary
    maxerr = 0d0

    !-------------------------------------------------!
    ! Part 1: bottom wall

    ! metric factor for first derivative dy/dzeta
    g  = ( -3d0*y(1) + 4d0*y(2) - y(3) )/(2d0*dzeta)

    ! metric factor for second derivative d^2y/dzeta^2
    g2 = ( 2d0*y(1) - 5d0*y(2) + 4d0*y(3) - 1d0*y(4) )/dzeta**2d0

    !$acc parallel loop collapse(2) default(present) private(nui,dnu,dV,ddV)
    Do k = 2, nzg-1
       Do i = 2, nxg-1

          ! interpolated nu at faces
          nui = nu + 0.5d0*( avg_nu_t(1,1,1) + avg_nu_t(1,2,1) )

          ! derivative of nu at the wall
          dnu = ( avg_nu_t(1,2,1) - avg_nu_t(1,1,1) ) / ( yg(2)-yg(1) )

          ! derivative of V at the wall
          dV = ( -3d0*V(i,1,k) + 4d0*V(i,2,k) - V(i,3,k) )/(2d0*dzeta)*1d0/g

          ! second derivative of V at the wall
          ddV = ( 2d0*V(i,1,k) - 5d0*V(i,2,k) + 4d0*V(i,3,k) - 1d0*V(i,4,k) )/dzeta**2d0
          ddV = ( ddV - dV*g2 )/g**2d0

          ! boundary conditions
          bc_1(i,k) = dnu*dV + nui*ddV

       End Do
    End Do

    !-------------------------------------------------!
    ! Part 2: top wall

    ! metric factor for first derivative dy/dzeta
    g  = ( 3d0*y(ny) - 4d0*y(ny-1) + y(ny-2) )/(2d0*dzeta)

    ! metric factor for second derivative d^2y/dzeta^2
    g2 = ( 2d0*y(ny) - 5d0*y(ny-1) + 4d0*y(ny-2) - 1d0*y(ny-3) )/dzeta**2d0

    !$acc parallel loop collapse(2) default(present) private(nui,dnu,dV,ddV)
    Do k = 2, nzg-1
       Do i = 2, nxg-1

          ! interpolated nu at faces
          nui = nu + 0.5d0*( avg_nu_t(1,nyg,1) + avg_nu_t(1,nyg-1,1) )

          ! derivative of nu at the wall
          dnu = ( avg_nu_t(1,nyg,1) - avg_nu_t(1,nyg-1,1) ) / ( yg(nyg)-yg(nyg-1) )

          ! derivative of V at the wall
          dV = ( 3d0*V(i,ny,k) - 4d0*V(i,ny-1,k) + V(i,ny-2,k) )/(2d0*dzeta)*1d0/g

          ! second derivative of V at the wall
          ddV = ( 2d0*V(i,ny,k) - 5d0*V(i,ny-1,k) + 4d0*V(i,ny-2,k) - 1d0*V(i,ny-3,k) )/dzeta**2d0
          ddV = ( ddV - dV*g2 )/g**2d0

          ! boundary conditions
          bc_2(i,k) = dnu*dV + nui*ddV

       End Do
    End Do

  End Subroutine compute_pressure_BC

End Module pressure
