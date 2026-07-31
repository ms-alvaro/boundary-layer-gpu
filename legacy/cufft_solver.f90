!------------------------------------------------!
!  Module for cuFFT-based transforms              !
!  Forward: 2D batched (2*nxp × nzp)             !
!  Inverse: split z-IFFT + x-IFFT                !
!------------------------------------------------!
Module cufft_solver

  Use iso_fortran_env, Only : Int32, Int64
  Use, Intrinsic :: iso_c_binding
  Use cufft

  Implicit None

  Integer :: cufft_plan_2d = 0   ! 2D C2C for forward (2*nxp × nzp × nm batches)
  Integer :: cufft_plan_zi = 0   ! 1D C2C for inverse z (on rhs_hat_gpu)
  Integer :: cufft_plan_xi = 0   ! 1D C2C for inverse x (on plane_gpu)
  Logical :: use_cufft = .False.

  Private
  Public :: use_cufft, cufft_plan_2d, cufft_plan_zi, cufft_plan_xi
  Public :: cufft_init_plans, cufft_destroy_plans
  Public :: cufftExecZ2Z, CUFFT_FORWARD, CUFFT_INVERSE, CUFFT_SUCCESS

Contains

  Subroutine cufft_init_plans(nxp_in, nzp_in, nm_in)
    Integer, Intent(In) :: nxp_in, nzp_in, nm_in
    Integer :: istat, n2(2), ne(2), n1(1), ne1(1)
    Integer :: nx2

    nx2 = 2 * nxp_in

    ! Plan 1: 2D C2C forward — plane_gpu(2*nxp, nzp, nm)
    n2(1) = nzp_in;  n2(2) = nx2
    ne(1) = nzp_in;  ne(2) = nx2
    istat = cufftPlanMany(cufft_plan_2d, 2, n2, &
         ne, 1, nx2*nzp_in, ne, 1, nx2*nzp_in, CUFFT_Z2Z, nm_in)
    If (istat /= CUFFT_SUCCESS) Then
       Write(*,*) 'cufftPlanMany 2D failed:', istat; Stop
    End If

    ! Plan 2: 1D C2C inverse z — rhs_hat_gpu(nxp, nm, nzp)
    ! stride = nxp*nm (between consecutive z-values)
    ! dist = 1, batch = nxp*nm
    n1(1) = nzp_in; ne1(1) = nzp_in
    istat = cufftPlanMany(cufft_plan_zi, 1, n1, &
         ne1, nxp_in*nm_in, 1, ne1, nxp_in*nm_in, 1, CUFFT_Z2Z, nxp_in*nm_in)
    If (istat /= CUFFT_SUCCESS) Then
       Write(*,*) 'cufftPlanMany zi failed:', istat; Stop
    End If

    ! Plan 3: 1D C2C inverse x — plane_gpu(2*nxp, nzp, nm)
    ! stride = 1, dist = 2*nxp, batch = nzp*nm
    n1(1) = nx2; ne1(1) = nx2
    istat = cufftPlanMany(cufft_plan_xi, 1, n1, &
         ne1, 1, nx2, ne1, 1, nx2, CUFFT_Z2Z, nzp_in*nm_in)
    If (istat /= CUFFT_SUCCESS) Then
       Write(*,*) 'cufftPlanMany xi failed:', istat; Stop
    End If

    use_cufft = .True.
  End Subroutine cufft_init_plans

  Subroutine cufft_destroy_plans
    Integer :: istat
    If (use_cufft) Then
       istat = cufftDestroy(cufft_plan_2d)
       istat = cufftDestroy(cufft_plan_zi)
       istat = cufftDestroy(cufft_plan_xi)
       use_cufft = .False.
    End If
  End Subroutine cufft_destroy_plans

End Module cufft_solver
