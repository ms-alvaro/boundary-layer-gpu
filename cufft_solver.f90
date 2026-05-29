!------------------------------------------------!
!  Module for cuFFT-based 2D transforms          !
!  Replaces FFTW-MPI for single-rank runs        !
!------------------------------------------------!
Module cufft_solver

  Use iso_fortran_env, Only : Int32, Int64
  Use, Intrinsic :: iso_c_binding
  Use cufft

  Implicit None

  ! cuFFT plan handles
  Integer :: cufft_plan_fwd = 0
  Integer :: cufft_plan_inv = 0

  ! Flag: use cuFFT (true when nprocs==1)
  Logical :: use_cufft = .False.

  Private
  Public :: use_cufft, cufft_plan_fwd, cufft_plan_inv
  Public :: cufft_init_plans, cufft_destroy_plans
  Public :: cufftExecZ2Z, CUFFT_FORWARD, CUFFT_INVERSE, CUFFT_SUCCESS

Contains

  !--------------------------------------------!
  !     Initialize cuFFT 2D Z2Z plans         !
  !--------------------------------------------!
  Subroutine cufft_init_plans(nz_dim, nx_dim)

    Integer, Intent(In) :: nz_dim, nx_dim
    Integer :: istat

    ! Create 2D Z2Z (complex-to-complex double precision) plans
    ! cuFFT uses C-order dimensions like FFTW: (nz, nx)
    istat = cufftPlan2d(cufft_plan_fwd, nz_dim, nx_dim, CUFFT_Z2Z)
    If (istat /= CUFFT_SUCCESS) Then
       Write(*,*) 'Error: cufftPlan2d forward failed, status=', istat
       Stop
    End If

    istat = cufftPlan2d(cufft_plan_inv, nz_dim, nx_dim, CUFFT_Z2Z)
    If (istat /= CUFFT_SUCCESS) Then
       Write(*,*) 'Error: cufftPlan2d inverse failed, status=', istat
       Stop
    End If

    use_cufft = .True.

  End Subroutine cufft_init_plans

  !--------------------------------------------!
  !        Destroy cuFFT plans                !
  !--------------------------------------------!
  Subroutine cufft_destroy_plans

    Integer :: istat

    If (use_cufft) Then
       istat = cufftDestroy(cufft_plan_fwd)
       istat = cufftDestroy(cufft_plan_inv)
       use_cufft = .False.
    End If

  End Subroutine cufft_destroy_plans

End Module cufft_solver
