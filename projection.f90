!------------------------------------------------!
!      Module for fractional step method         !
!------------------------------------------------!
Module projection

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use cufft_solver

  ! prevent implicit typing
  Implicit None

Contains

  !--------------------------------------------!
  !    Compute incompressible velocity with    !
  !         fractional step method             !
  !--------------------------------------------!
  Subroutine compute_projection_step

    Call compute_pseudo_pressure_rhs

    Call solve_poisson_equation

    Call project_velocity

  End Subroutine compute_projection_step

  !------------------------------------------------------!
  ! Compute right-hand size for pseudo-pressure equation !
  !------------------------------------------------------!
  Subroutine compute_pseudo_pressure_rhs

    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) default(present)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             rhs_p(i,j,k)  = ( U(i,j,k) - U(i-1,j,k) )/( x(i)-x(i-1) ) + &
                             ( V(i,j,k) - V(i,j-1,k) )/( y(j)-y(j-1) ) + &
                             ( W(i,j,k) - W(i,j,k-1) )/( z(k)-z(k-1) )
          End Do
       End Do
    End Do

  End Subroutine compute_pseudo_pressure_rhs

  !----------------------------------------------------!
  !             Solve pseudo-pressure equation         !
  !----------------------------------------------------!
  Subroutine solve_poisson_equation

    Integer(Int32) :: ii, i, j, k, k_global, i_global, info
    Real   (Int64) :: dum, dumref, maxerr
    Integer :: istat

    External :: zgesv

    If ( use_cufft ) Then
       Call solve_poisson_gpu
    Else
       Call solve_poisson_cpu
    End If

    ! periodic boundary conditions in z
    Call apply_periodic_z_pressure

    ! update ghost interior planes
    Call update_ghost_interior_planes_pressure

    ! save pressure for statistics
    If ( rk_step == 1 ) Then
       !$acc kernels default(present)
       P( 2:nxg-1, 2:nyg-1, 2:nzg-1 ) = rhs_p/(dt*rk2_coef(1,1))
       P(  1,:,:) = P(    2,:,:)
       P(nxg,:,:) = P(nxg-1,:,:)
       P(:,  1,:) = P(:,    2,:)
       P(:,nyg,:) = P(:,nyg-1,:)
       P(:,:,  1) = P(:,:,    2)
       P(:,:,nzg) = P(:,:,nzg-1)
       !$acc end kernels
    End If

  End Subroutine solve_poisson_equation

  !----------------------------------------------------!
  !     GPU pressure solver using cuFFT               !
  !     All operations on GPU, no host transfers       !
  !----------------------------------------------------!
  Subroutine solve_poisson_gpu

    Integer(Int32) :: ii, i, j, k, i_global, k_global, info, jj
    Integer :: istat

    External :: zgesv

    ! Forward cosine+Fourier transform per y-slice (on GPU)
    Do j = 2, nyg-1
       ! Pack rhs_p into plane_short
       !$acc kernels default(present)
       plane_short = dcmplx( rhs_p ( 2:nxp+1, j, 2:nzp+1 ) )
       plane_gpu   = (0d0,0d0)
       !$acc end kernels

       ! Cosine extension
       !$acc parallel loop default(present)
       Do ii = 1, Int(nxp)
          plane_gpu(2*ii,   1:nzp) = plane_short(ii, 1:nzp)
          plane_gpu(4*Int(nxp)-2*ii+2, 1:nzp) = plane_short(ii, 1:nzp)
       End Do

       ! cuFFT forward (in-place on GPU)
       !$acc host_data use_device(plane_gpu)
       istat = cufftExecZ2Z(cufft_plan_fwd, plane_gpu, plane_gpu, CUFFT_FORWARD)
       !$acc end host_data

       ! Extract modes (0-indexed modes -> 1-indexed plane_gpu)
       !$acc kernels default(present)
       rhs_p_hat(0:mx, j, 0:mz) = plane_gpu(1:mx+1, 1:mz+1)
       !$acc end kernels
    End Do

    ! Tridiagonal solve per mode — on CPU (small data transfer)
    !$acc update self(rhs_p_hat)
    Do k = 0, Int(mz)
       Do i = 0, Int(mx)
          i_global = imode_map( i )
          k_global = kmode_map( k )
          Do jj = 2, nyg-1
             D(jj) = Dyy(jj,jj) + kxx(i_global) + kzz(k_global)
          End Do
          Do jj = 2, nyg-2
             DL(jj) = Dyy(jj+1,jj)
          End Do
          Do jj = 2, nyg-2
             DU(jj) = Dyy(jj,jj+1)
          End Do
          If ( i_global==0 .And. k_global==0 ) D(2) = 3d0/2d0*D(2)
          rhs_aux = rhs_p_hat(i, :, k)
          Call Zgtsv( nr, nrhs, DL, D, DU, rhs_aux, nr, info)
          rhs_p_hat(i, :, k) = rhs_aux
       End Do
    End Do
    !$acc update device(rhs_p_hat)

    ! Inverse cosine+Fourier transform per y-slice (on GPU)
    Do j = 2, nyg-1
       ! Build full spectrum for inverse cosine transform
       !$acc kernels default(present)
       plane_gpu = (0d0,0d0)
       plane_gpu(1:nxp, 1:nzp) = rhs_p_hat(0:nxp-1, j, 0:mz)
       !$acc end kernels

       !$acc parallel loop default(present)
       Do ii = Int(nxp)-1, 0, -1
          plane_gpu(Int(nxp)+2-1-Int(nxp)+1+ii, 1:nzp) = -rhs_p_hat(Int(nxp)-1-ii, j, 0:mz)
       End Do

       !$acc parallel loop default(present)
       Do ii = 1, Int(nxp)-1
          plane_gpu(2*Int(nxp)+2-1+ii-1, 1:nzp) = -rhs_p_hat(ii, j, 0:mz)
       End Do

       !$acc parallel loop default(present)
       Do ii = Int(nxp)-1, 1, -1
          plane_gpu(4*Int(nxp)-Int(nxp)+2-1+Int(nxp)-1-ii, 1:nzp) = rhs_p_hat(ii, j, 0:mz)
       End Do

       ! cuFFT inverse (in-place on GPU)
       !$acc host_data use_device(plane_gpu)
       istat = cufftExecZ2Z(cufft_plan_inv, plane_gpu, plane_gpu, CUFFT_INVERSE)
       !$acc end host_data

       ! Extract real part from even positions
       !$acc parallel loop default(present)
       Do ii = 1, Int(nxp)
          rhs_p( ii+1, j, 2:nzp+1 ) = Real(plane_gpu(2*ii,:), 8)/Real(nxpe_global*nzp_global, 8)
       End Do
    End Do

  End Subroutine solve_poisson_gpu

  !----------------------------------------------------!
  !     CPU pressure solver using FFTW-MPI            !
  !     Fallback for multi-rank runs                  !
  !----------------------------------------------------!
  Subroutine solve_poisson_cpu

    Integer(Int32) :: ii, i, j, k, k_global, i_global, info

    External :: zgesv

    ! Transfer rhs_p from GPU to CPU for FFT
    !$acc update self(rhs_p)

    ! Forward cosine+Fourier transform
    Do j = 2, nyg-1
       plane_short = dcmplx( rhs_p ( 2:nxp+1, j, 2:nzp+1 ) )
       plane       = (0d0,0d0)
       ii = 1
       Do i=2,2*nxp,2
          plane(i,:) = plane_short(ii,:)
          ii = ii + 1
       End Do
       ii = 1
       Do i=4*nxp,2*nxp+2,-2
          plane(i,:) = plane_short(ii,:)
          ii = ii + 1
       End Do
       Call fftw_mpi_execute_dft(plan_d, plane, plane)
       rhs_p_hat(0:mx, j, 0:mz) = plane(1:mx+1, 1:mz+1)
    End Do

    ! Tridiagonal solve
    Do k = 0, mz
       Do i = 0, mx
          i_global = imode_map( i )
          k_global = kmode_map( k )
          Do j = 2, nyg-1
             D(j) = Dyy(j,j) + kxx(i_global) + kzz(k_global)
          End Do
          Do j = 2, nyg-2
             DL(j) = Dyy(j+1,j)
          End Do
          Do j = 2, nyg-2
             DU(j) = Dyy(j,j+1)
          End Do
          If ( i_global==0 .And. k_global==0 ) D(2) = 3d0/2d0*D(2)
          rhs_aux = rhs_p_hat (i, :, k)
          Call Zgtsv( nr, nrhs, DL, D, DU, rhs_aux, nr, info)
          rhs_p_hat (i, :, k) = rhs_aux
       End Do
    End Do

    ! Inverse cosine+Fourier transform
    Do j = 2, nyg-1
       plane               = (0d0,0d0)
       plane(1:nxp, 1:nzp) = rhs_p_hat(0:nxp-1,j,0:mz)
       ii = nxp+2-1
       Do i=nxp-1,1-1,-1
          plane(ii+1, 1:nzp) = -rhs_p_hat(i,j,0:mz)
          ii = ii + 1
       End Do
       ii = 2-1
       Do i=2*nxp+2-1,3*nxp-1
          plane(i+1, 1:nzp) = -rhs_p_hat(ii,j,0:mz)
          ii = ii + 1
       End Do
       ii = 4*nxp-nxp+2-1
       Do i=nxp-1,2-1,-1
          plane(ii+1, 1:nzp) = rhs_p_hat(i,j,0:mz)
          ii = ii + 1
       End Do
       Call fftw_mpi_execute_dft(plan_i, plane, plane)
       ii = 2
       Do i=2,2*nxp,2
          rhs_p( ii, j, 2:nzp+1 ) = plane(i,:)/Real( nxpe_global*nzp_global, 8)
          ii = ii + 1
       End Do
    End Do

    ! Transfer rhs_p back to GPU
    !$acc update device(rhs_p)

  End Subroutine solve_poisson_cpu

  !--------------------------------------------------!
  ! Periodic boundary conditions for pseudo-pressure !
  !--------------------------------------------------!
  Subroutine apply_periodic_z_pressure

    !$acc kernels default(present)
    rhs_p( nxg-1, :, : ) = rhs_p( nxg-2, :, : )
    !$acc end kernels

    If ( nprocs==1 ) Then
       !$acc kernels default(present)
       rhs_p( 2:nxg-1, :, nzp+1+1 ) = rhs_p ( 2:nxg-1, :, 2 )
       !$acc end kernels
       Return
    End If

    !$acc update self(rhs_p)
    If     ( myid==0 ) Then
       buffer_p = rhs_p ( 2:nxg-1, :, 2 )
       Call Mpi_send(buffer_p, (nxg-2)*(nyg-2), MPI_real8, nprocs-1, 0, &
             MPI_COMM_WORLD,ierr)
    Elseif ( myid==nprocs-1 ) Then
       Call Mpi_recv(buffer_p, (nxg-2)*(nyg-2), MPI_real8, 0, 0, &
            MPI_COMM_WORLD,istat,ierr)
       rhs_p ( 2:nxg-1, :, nzp+1+1 ) = buffer_p
    End If
    !$acc update device(rhs_p)

  End Subroutine apply_periodic_z_pressure

  !------------------------------------------------!
  !    Update ghost interior planes for pressure   !
  !------------------------------------------------!
  Subroutine update_ghost_interior_planes_pressure

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom

    If ( nprocs==1 ) Return

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
    !$acc update self(rhs_p(:,:,2))
    buffer_ps = rhs_p(:,:,2)
    Call Mpi_sendrecv(buffer_ps, (nxg-2)*(nyg-2), Mpi_real8, sendto, tagto,        &
         buffer_pr, (nxg-2)*(nyg-2), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
         istat, ierr)
    If ( myid/=nprocs-1 ) rhs_p(:,:,nzg) = buffer_pr
    !$acc update device(rhs_p(:,:,nzg))

  End Subroutine update_ghost_interior_planes_pressure

  !------------------------------------------------!
  !             Project velocity field             !
  !------------------------------------------------!
  Subroutine project_velocity

    Integer(Int32) :: i, j, k

    !$acc parallel loop default(present)
    Do i = 2, nx-1
       U(i,2:nyg-1,2:nzg-1) = U(i,2:nyg-1,2:nzg-1) - &
            ( rhs_p(i+1,2:nyg-1,2:nzg-1) - rhs_p(i,2:nyg-1,2:nzg-1) )/( xg(i+1) - xg(i) )
    End Do

    !$acc parallel loop default(present)
    Do j = 2, ny-1
       V(2:nxg-1,j,2:nzg-1) = V(2:nxg-1,j,2:nzg-1) - &
            ( rhs_p(2:nxg-1,j+1,2:nzg-1) - rhs_p(2:nxg-1,j,2:nzg-1) )/( yg(j+1) - yg(j) )
    End Do

    !$acc parallel loop default(present)
    Do k = 2, nz-1
       W(2:nxg-1,2:nyg-1,k) = W(2:nxg-1,2:nyg-1,k) - &
            ( rhs_p(2:nxg-1,2:nyg-1,k+1) - rhs_p(2:nxg-1,2:nyg-1,k) )/( zg(k+1) - zg(k) )
    End Do

  End Subroutine project_velocity

  !----------------------------------------------------!
  !   Compute maximum divergence for interior points   !
  !----------------------------------------------------!
  Subroutine check_divergence(max_divergence)

    Real(Int64), Intent(Out) :: max_divergence

    Real   (Int64) :: max_divergence_local, div
    Integer(Int32) :: i, j, k

    div = 0d0
    !$acc parallel loop collapse(3) default(present) reduction(max:div)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-2
             div  = Max( div, Abs( ( U(i,j,k) - U(i-1,j,k) )/( x(i)-x(i-1) ) + &
                                   ( V(i,j,k) - V(i,j-1,k) )/( y(j)-y(j-1) ) + &
                                   ( W(i,j,k) - W(i,j,k-1) )/( z(k)-z(k-1) )   ) )
          End Do
       End Do
    End Do

    max_divergence_local = div
    Call MPI_Reduce(max_divergence_local,max_divergence,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

  End Subroutine check_divergence

End module projection
