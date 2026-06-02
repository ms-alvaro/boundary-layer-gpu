!------------------------------------------------!
!  Module for cuFFT-based 2D transforms          !
!  Replaces FFTW-MPI for single-rank runs        !
!  Batched: all y-planes in one call             !
!------------------------------------------------!
Module cufft_solver

  Use iso_fortran_env, Only : Int32, Int64
  Use, Intrinsic :: iso_c_binding
  Use cufft

  Implicit None

  ! Single cuFFT plan (used for both forward and inverse)
  Integer :: cufft_plan = 0

  ! Number of batches (= nyg-2 interior y-planes)
  Integer :: cufft_nbatch = 0

  ! Flag: use cuFFT (true when nprocs==1)
  Logical :: use_cufft = .False.

  Private
  Public :: use_cufft, cufft_plan, cufft_nbatch
  Public :: cufft_init_plans, cufft_destroy_plans
  Public :: cufftExecZ2Z, CUFFT_FORWARD, CUFFT_INVERSE, CUFFT_SUCCESS

Contains

  !--------------------------------------------!
  !  Initialize cuFFT batched Z2Z plan         !
  !  plane_gpu(nx_dim, nz_dim, nbatch)         !
  !--------------------------------------------!
  Subroutine cufft_init_plans(nz_dim, nx_dim, nbatch)

    Integer, Intent(In) :: nz_dim, nx_dim, nbatch
    Integer :: istat
    Integer :: rank_fft, n_arr(2), inembed(2), onembed(2)
    Integer :: istride, idist, ostride, odist

    cufft_nbatch = nbatch
    rank_fft = 2

    ! cuFFT uses C-order: n_arr(1)=slow dim, n_arr(2)=fast dim
    ! Fortran column-major plane_gpu(nx_dim, nz_dim, nbatch):
    !   element (i,k,b) at offset (i-1) + (k-1)*nx_dim + (b-1)*nx_dim*nz_dim
    n_arr(1) = nz_dim
    n_arr(2) = nx_dim

    istride = 1
    idist   = nx_dim * nz_dim
    inembed(1) = nz_dim
    inembed(2) = nx_dim

    ostride = istride
    odist   = idist
    onembed = inembed

    istat = cufftPlanMany(cufft_plan, rank_fft, n_arr, &
         inembed, istride, idist, onembed, ostride, odist, &
         CUFFT_Z2Z, nbatch)
    If (istat /= CUFFT_SUCCESS) Then
       Write(*,*) 'Error: cufftPlanMany failed, status=', istat
       Write(*,*) '  nz=', nz_dim, ' nx=', nx_dim, ' batch=', nbatch
       Stop
    End If

    use_cufft = .True.

  End Subroutine cufft_init_plans

  Subroutine cufft_destroy_plans
    Integer :: istat
    If (use_cufft) Then
       istat = cufftDestroy(cufft_plan)
       use_cufft = .False.
    End If
  End Subroutine cufft_destroy_plans

End Module cufft_solver
