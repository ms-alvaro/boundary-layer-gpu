!--------------------------------------------!
! Module for computing some basic statistics !
!--------------------------------------------!
Module statistics

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use interpolation
  Use input_output
  Use equations
  Use boundary_conditions

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------------!
  !                                                    !
  !      Compute some basic statistics on the fly      !
  !       Statsitics computed at x and y positions     !
  !                                                    !
  !----------------------------------------------------!
  Subroutine compute_statistics    

    Integer(Int32) :: ii, jj, j, jref
    Real   (Int64) :: Uinf, theta_inlet, Umean99, w0, w1, T_mean
    Real   (Int64) :: tau_sgs_wall(nx_global)
    Real   (Int64), Dimension(ny) :: temp_1d

    ! if pressure not computed 
    pressure_computed = .False.

    ! statistics computed at grid y -> U and W interpolated    
    If ( Mod(istep,nstats)==0 .Or. istep==1 ) Then

       ! Compute actual pressure (should be called first, uses term_1,...)
       !Call compute_pressure
       ! now computed in projection.f90
       pressure_computed = .True.

       ! interpolate W in x and y -> term
       Call interpolate_x(   W,term_1,1)
       Call interpolate_y(term_1,term,2)
       
       ! interpolate U in y -> term_1
       Call interpolate_y(U,term_1,2)

       ! interpolate V in x -> term_2
       Call interpolate_x(V,term_2,1)

       ! interpolate P in y -> term_3
      ! Call interpolate_y(P,term,2)

       ! Transfer interpolated fields from GPU to CPU for statistics
       !$acc update self(term,term_1,term_2)

       ! compute local statistics
       Do ii=1,nx
          Do jj=1,ny
             
             Umean  (ii,jj)   = Sum(term_1 (ii, jj, 2:nzg-1) )
             Vmean  (ii,jj)   = Sum(term_2 (ii, jj, 2:nzg-1) )
             Wmean  (ii,jj)   = Sum(term   (ii, jj, 1:nz-2 ) )

             U2mean (ii,jj)   = Sum(term_1 (ii, jj, 2:nzg-1)**2d0 )
             V2mean (ii,jj)   = Sum(term_2 (ii, jj, 2:nzg-1)**2d0 )
             W2mean (ii,jj)   = Sum(term   (ii, jj, 1:nz-2 )**2d0 )             
             
             UVmean (ii,jj)   = Sum( term_1(ii,jj,2:nzg-1)*term_2(ii,jj,2:nzg-1) )
             
             Pmean  (ii,jj)   = Sum( P (ii, jj, 2:nzg-1)      )
             P2mean (ii,jj)   = Sum( P (ii, jj, 2:nzg-1)**2d0 )

             nu_t_mean(ii,jj) = Sum( nu_t(ii, jj, 2:nzg-1) )
             
          End Do
       End Do

       ! reduce statatistics between processors      
       IF ( myid==0 ) Then

          Call MPI_Reduce(MPI_IN_PLACE,Umean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,Vmean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,Wmean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(MPI_IN_PLACE,U2mean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,V2mean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,W2mean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(MPI_IN_PLACE,UVmean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

          Call MPI_Reduce(MPI_IN_PLACE, Pmean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,P2mean,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

          Call MPI_Reduce(MPI_IN_PLACE,nu_t_mean,nxg*nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
                    
       Else

          Call MPI_Reduce(Umean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(Vmean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(Wmean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(U2mean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(V2mean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(W2mean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(UVmean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

          Call MPI_Reduce( Pmean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(P2mean,0,nx*ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

          Call MPI_Reduce(nu_t_mean,0,nxg*nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
       End If

       !-------------------------------------------------------!
       !     These statistics are only good for processor 0    !
       !-------------------------------------------------------!
       Umean  = Umean/Real( nzg_global-2, 8)
       Vmean  = Vmean/Real( nzg_global-2, 8)
       Wmean  = Wmean/Real(  nz_global-2, 8)

       U2mean = U2mean/Real( nzg_global-2, 8)
       V2mean = V2mean/Real( nzg_global-2, 8)
       W2mean = W2mean/Real(  nz_global-2, 8)

       UVmean = UVmean/Real( nzg_global-2 ,8)

       Pmean  =  Pmean/Real( nzg_global-2 ,8)
       P2mean = P2mean/Real( nzg_global-2 ,8)

       nu_t_mean = nu_t_mean/Real( nzg_global-2 ,8)

       ! Mean derivative at the wall
       Do ii=1,nx
          Uaux_1_local(ii) = Sum( U(ii,1,2:nzg-1) )
          Uaux_2_local(ii) = Sum( U(ii,2,2:nzg-1) )
       End Do
       Uaux_1 = 0d0
       Uaux_2 = 0d0
       Call MPI_Reduce(Uaux_1_local,Uaux_1,nx,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(Uaux_2_local,Uaux_2,nx,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Uaux_1    = Uaux_1/Real( nzg_global-2, 8 )
       Uaux_2    = Uaux_2/Real( nzg_global-2, 8 )
       dUdy_wall = ( Uaux_2 - Uaux_1 )/( yg(2) - yg(1) )

       ! Mean Reynolds stress at the wall
       UV_wall = UVmean(:,1) 

       ! SGS stress at the wall
       Do ii=1,nx_global
          tau_sgs_wall(ii) = 0.25d0*( nu_t_mean(ii,1) + nu_t_mean(ii,2) + nu_t_mean(ii+1,1) + nu_t_mean(ii+1,2) )*dUdy_wall(ii)
       End Do

       ! Skin friction coefficient Cf = tau_w/(1/2*rho*U_inf^2)
       Uinf        = 1d0 ! U(1,nyg,1)
       Cf          = 2d0*( -UV_wall + nu*dUdy_wall + tau_sgs_wall )/Uinf**2d0
       utau_wall   = ( 0.5d0*Uinf**2d0*Abs(Cf) )**0.5d0
       If (istep == 1) utau_wall_T = utau_wall
       !!!!!!!!!!!!!!!!!!!!!!!!
       T_mean      = alpha_std 
       !!!!!!!!!!!!!!!!!!!!!!!!
       utau_wall_T = dt/T_mean*utau_wall + (1d0-dt/T_mean)*utau_wall_T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       If (myid==0) Then
!          Write(*,*) '--------------------------------Cf'
!          Do ii=1,nx
!             Write(*,*) Cf(ii), 0.027d0*(x(ii)/nu)**(-1d0/7d0)
!          end Do
!       End If
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! momentum boundary layer thickness
       theta_inlet = 0d0
       temp_1d     = Umean(1,:)/Uinf*(1d0 - Umean(1,:)/Uinf)
       Do j = 2,ny
          theta_inlet = theta_inlet + 0.5d0*( temp_1d(j) + temp_1d(j-1) )*( y(j)-y(j-1) )
       End Do

       ! 99% boundary layer thickness
       Umean99 = 0.99d0*Uinf
       jref    = 0
       Do j = 1,ny
          If ( Umean(1,j)>=Umean99 ) then
             jref = j
             Exit
          End If
       End Do
       If ( jref<2 ) jref = 2 
       w1                = ( Umean99 - Umean(1,jref-1) )/( Umean(1,jref) - Umean(1,jref-1) )
       w0                = 1d0 - w1
       delta99_inlet_ins = w1*y(jref) + w0*y(jref-1)

       ! Reynolds numbers at inlet
       Rex_inlet     =       x_global(1)*Uinf/nu
       Retheta_inlet =       theta_inlet*Uinf/nu
       Redelta_inlet = delta99_inlet_ins*Uinf/nu

       ! write statistics
!!!!!!!!!!!!!!!!!!!
       Call output_statistics       
!!!!!!!!!!!!!!!!!!!

       ! Sanity check
       ! NaN check (disabled: nvfortran 24.3 ICE workaround)
       !If ( Any( ieee_is_nan(U) ) ) Stop 'Error U NaNs!'
       !If ( Any( ieee_is_nan(V) ) ) Stop 'Error V NaNs!'
       !If ( Any( ieee_is_nan(W) ) ) Stop 'Error W NaNs!'
       
    End If
   
  End Subroutine compute_statistics

  !----------------------------------------------------!
  !                                                    !
  !       Compute z-modes and save them for latter     !
  !                                                    !
  !----------------------------------------------------!
  ! deleted, should be re-adapted
  Subroutine compute_statistics_z_modes

    Use fftz

    Real   (Int64) :: nu_original
    Integer(Int32) :: j, iproc, nzge
    
  End Subroutine compute_statistics_z_modes

End Module statistics
