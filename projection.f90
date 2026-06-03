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

    Integer(Int32) :: ii, i, j, k, b, jj, nm
    Integer(Int32) :: nxp_i, nx2_i, nzp_i, mx_i, mz_i
    Integer :: istat
    Real(Int64) :: fft_norm

    nm    = Int(nyg) - 2
    nxp_i = Int(nxp)
    nx2_i = 2 * nxp_i
    nzp_i = Int(nzp)
    mx_i  = Int(mx) + 1    ! = nxp for nprocs=1
    mz_i  = Int(mz) + 1    ! = nzp for nprocs=1
    fft_norm = 1d0 / Real(2*nxp_i*nzp_i, 8)

    ! ======== FORWARD: symmetric pack → 2D FFT → twiddle extract ========

    ! Pack DCT-II symmetric extension: g[0..N-1]=f, g[2N-1..N]=f[0..N-1]
    !$acc parallel loop collapse(3) default(present)
    Do b = 1, nm
       Do k = 1, nzp_i
          Do ii = 1, nxp_i
             plane_gpu(ii,          k, b) = dcmplx(rhs_p(ii+1, b+1, k+1))
             plane_gpu(nx2_i+1-ii,  k, b) = dcmplx(rhs_p(ii+1, b+1, k+1))
          End Do
       End Do
    End Do

    ! 2D batched forward FFT (2*nxp × nzp, nm batches)
    !$acc host_data use_device(plane_gpu)
    istat = cufftExecZ2Z(cufft_plan_2d, plane_gpu, plane_gpu, CUFFT_FORWARD)
    !$acc end host_data

    ! Extract modes with DCT-II twiddle
    !$acc parallel loop collapse(3) default(present)
    Do b = 1, nm
       Do k = 1, mz_i
          Do i = 1, mx_i
             rhs_hat_gpu(i, b, k) = plane_gpu(i, k, b) * dct_twiddle(i-1)
          End Do
       End Do
    End Do

    ! ======== THOMAS TRIDIAGONAL SOLVE ========
    !$acc parallel loop collapse(2) default(present)
    Do k = 1, mz_i
       Do i = 1, mx_i
          Do jj = 2, nm
             rhs_hat_gpu(i, jj, k) = rhs_hat_gpu(i, jj, k) &
                - thomas_dl_fact(i, k, jj-1) * rhs_hat_gpu(i, jj-1, k)
          End Do
          rhs_hat_gpu(i, nm, k) = rhs_hat_gpu(i, nm, k) * thomas_d_pivot(i, k, nm)
          Do jj = nm-1, 1, -1
             rhs_hat_gpu(i, jj, k) = (rhs_hat_gpu(i, jj, k) &
                - thomas_du(jj) * rhs_hat_gpu(i, jj+1, k)) * thomas_d_pivot(i, k, jj)
          End Do
       End Do
    End Do

    ! ======== INVERSE: z-IFFT first (makes data real), then x-IDCT ========

    ! Step 1: z-IFFT on rhs_hat_gpu (nzp-point, batched over nxp*nm)
    ! After this, rhs_hat_gpu(i,b,k) is in z-physical space
    !$acc host_data use_device(rhs_hat_gpu)
    istat = cufftExecZ2Z(cufft_plan_zi, rhs_hat_gpu, rhs_hat_gpu, CUFFT_INVERSE)
    !$acc end host_data

    ! Step 2: un-twiddle + Hermitian extend into plane_gpu for x-IFFT
    ! Data is now real in z-physical, so Hermitian in x is valid
    !$acc parallel loop collapse(3) default(present)
    Do b = 1, nm
       Do k = 1, nzp_i
          Do ii = 1, nx2_i
             plane_gpu(ii, k, b) = (0d0, 0d0)
          End Do
       End Do
    End Do
    ! First nxp modes with inverse twiddle
    !$acc parallel loop collapse(3) default(present)
    Do b = 1, nm
       Do k = 1, nzp_i
          Do ii = 1, nxp_i
             plane_gpu(ii, k, b) = rhs_hat_gpu(ii, b, k) * conjg(dct_twiddle(ii-1))
          End Do
       End Do
    End Do
    ! Hermitian extension: plane(2N+2-ii) = conj(plane(ii)) for ii=2..nxp
    !$acc parallel loop collapse(3) default(present)
    Do b = 1, nm
       Do k = 1, nzp_i
          Do ii = 2, nxp_i
             plane_gpu(nx2_i+2-ii, k, b) = conjg(plane_gpu(ii, k, b))
          End Do
       End Do
    End Do

    ! Step 3: x-IFFT (1D, batched over nzp*nm)
    !$acc host_data use_device(plane_gpu)
    istat = cufftExecZ2Z(cufft_plan_xi, plane_gpu, plane_gpu, CUFFT_INVERSE)
    !$acc end host_data

    ! Extract physical-space result
    !$acc parallel loop collapse(3) default(present)
    Do b = 1, nm
       Do k = 1, nzp_i
          Do ii = 1, nxp_i
             rhs_p(ii+1, b+1, k+1) = Real(plane_gpu(ii, k, b), 8) * fft_norm
          End Do
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

    !$acc parallel loop collapse(3) default(present)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             U(i,j,k) = U(i,j,k) - (rhs_p(i+1,j,k) - rhs_p(i,j,k))/(xg(i+1) - xg(i))
          End Do
       End Do
    End Do

    !$acc parallel loop collapse(3) default(present)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             V(i,j,k) = V(i,j,k) - (rhs_p(i,j+1,k) - rhs_p(i,j,k))/(yg(j+1) - yg(j))
          End Do
       End Do
    End Do

    !$acc parallel loop collapse(3) default(present)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             W(i,j,k) = W(i,j,k) - (rhs_p(i,j,k+1) - rhs_p(i,j,k))/(zg(k+1) - zg(k))
          End Do
       End Do
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
