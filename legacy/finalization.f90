!-----------------------------------------!
!      Module to finalize FFT and etc     !
!-----------------------------------------!
Module finalization
  
  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use cufft_solver
  
  ! prevent implicit typing
  Implicit None
  
Contains

  !-------------------------------------------!
  !        finalize FFTW plans and MPI        !
  !-------------------------------------------!
  Subroutine finalize
  
    ! finalize cuFFT
    Call cufft_destroy_plans
    ! finalize FFTW
    Call dfftw_destroy_plan(plan_d)
    Call dfftw_destroy_plan(plan_i)
    
    ! finalized MPI
    Call Mpi_barrier (MPI_COMM_WORLD,ierr)
    Call Mpi_finalize(ierr)
    
    If ( myid==0 ) Then
       Write(*,*) 'Done!'
    End If
    
  End Subroutine finalize
  
End Module finalization
