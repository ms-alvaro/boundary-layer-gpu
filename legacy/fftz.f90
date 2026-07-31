!------------------------------------------------!
! Module to compute fast Fourier Transform in z  !
!------------------------------------------------!
Module fftz

  ! Modules
  Use iso_fortran_env, Only : Int32, Int64

  ! prevent implicit typing
  Implicit None

  ! module private variables
  Integer(Int32) :: N
  Integer(Int64) :: plan_fftz, plan_ifftz
  Real   (Int64), Allocatable, Dimension(:) :: line
  Complex(Int64), Allocatable, Dimension(:) :: line_hat
  
  ! set private variables
  Private 
  Public  :: initialize_fftz, finalize_fftz, &
             compute_fft_z, compute_ifft_z

  ! FFTW
  Include 'fftw3.f'

  ! declarations
Contains  

  !--------------------------------------------!
  !           Initialize fftz module           !
  !--------------------------------------------!
  Subroutine initialize_fftz(N_)

    Integer(Int32), Intent(In) :: N_

    N = N_

    ! allocate variables
    Allocate(line(N),line_hat(N/2+1))
    
    ! create fft plans
    Call dfftw_plan_dft_r2c_1d(plan_fftz ,N,line,line_hat,FFTW_ESTIMATE)
    Call dfftw_plan_dft_c2r_1d(plan_ifftz,N,line_hat,line,FFTW_ESTIMATE)

  End Subroutine initialize_fftz

  !--------------------------------------------!
  !           Initialize fftz module           !
  !--------------------------------------------!
  Subroutine finalize_fftz
   
    ! destroy fft plan
    Call dfftw_destroy_plan(plan_fftz)
    Call dfftw_destroy_plan(plan_ifftz)
    
  End Subroutine finalize_fftz
 
  !---------------------------------------------!
  !     Compute fft in z for 2D xz planes       !
  !                                             !
  ! Input:  U     (physical space)              !
  ! Output: U_hat (Fourier  space)              !
  !                                             !
  ! Sizes: U(1:nx,1:nzg), U_hat(1:nx,1:nzg/2+1) !
  !                                             !
  ! Wavenumbers: kz=2*pi/Lz*[0 1 2 ... nzg/2+1] !
  ! Wavelengths: lambda_z = 2pi/kz              !
  !---------------------------------------------!
  Subroutine compute_fft_z(U,U_hat)
    
    Real   (Int64), Dimension(:,:), Intent(In)  :: U
    Complex(Int64), Dimension(:,:), Intent(Out) :: U_hat

    ! local variables
    Integer(Int32), Dimension(2) :: nn
    Integer(Int32) :: n1, n2, i

    ! get size
    nn = Shape(U)
    n1 = nn(1)
    n2 = nn(2)

    ! sanity check
    If ( n2/=N ) Stop "Error, different size!"

    ! take fft
    Do i=1,n1
       line = U(i,:)
       Call dfftw_execute_dft_r2c(plan_fftz, line, line_hat)
       U_hat(i,:) = line_hat
    End Do
       
  End Subroutine compute_fft_z

  !---------------------------------------------!
  !  Compute inverse fft in z for 2D xz planes  !
  !                                             !
  ! Intput: U_hat (Fourier  space)              !
  ! Output: U     (physical space)              !
  !                                             !
  ! Sizes: U(1:nx,1:nzg), U_hat(1:nx,1:nzg/2+1) !
  !---------------------------------------------!
  Subroutine compute_ifft_z(U,U_hat)

    Complex(Int64), Dimension(:,:), Intent(In)  :: U_hat    
    Real   (Int64), Dimension(:,:), Intent(Out) :: U

    ! local variables
    Integer(Int32), Dimension(2) :: nn
    Integer(Int32) :: n1, n2, i

    ! get size
    nn = Shape(U_hat)
    n1 = nn(1)
    n2 = nn(2)

    ! sanity check
    If ( n2/=(N/2+1) ) Stop "Error, different size!"

    ! take ifft
    Do i=1,n1
       line_hat = U_hat(i,:)
       Call dfftw_execute_dft_c2r(plan_ifftz, line_hat, line)
       U(i,:) = line/Real(N,8)
    End Do
       
  End Subroutine compute_ifft_z

End Module fftz
