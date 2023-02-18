!--------------------------------------------!
! Module to monitor status of the simulation !
!--------------------------------------------!
Module monitor

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use projection, Only : check_divergence

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------!
  ! Output some key values during simulation !
  !------------------------------------------!
  Subroutine output_monitor

    Real(Int64) ::  maxU,  maxV,  maxW, maxNut, local_sum
    Real(Int64) :: meanU, meanV, meanW, meanNut
    Real(Int64) :: max_divergence

    If ( Mod(istep,nmonitor)==0 .Or. istep==1 ) Then

       ! compute mean values
       local_sum = Sum ( U(2:nx-1,2:nyg-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,meanU,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       local_sum = Sum ( V(2:nxg-1,2:ny-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,meanV,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       local_sum = Sum ( W(2:nxg-1,2:nyg-1,2:nz-1) )
       Call MPI_Reduce (local_sum,meanW,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       local_sum = Sum ( nu_t(2:nxg-1,2:nyg-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,meanNut,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       ! compute maximum values
       local_sum = Maxval ( U(2:nx-1,2:nyg-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,maxU,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       local_sum = Maxval ( V(2:nxg-1,2:ny-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,maxV,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       local_sum = Maxval ( W(2:nxg-1,2:nyg-1,2:nz-1) )
       Call MPI_Reduce (local_sum,maxW,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       local_sum = Maxval ( nu_t(2:nxg-1,2:nyg-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,maxNut,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       Call check_divergence(max_divergence)

       ! end measure time per step
       time2 = MPI_WTIME()

       ! processor 0 shows the results
       If ( myid==0 ) Then

          Write(*,*) 'step number :', istep
          Write(*,*) 'time        :', t
          Write(*,*) 'time step   :', dt
          
          Write(*,*) ' '
          Write(*,*) 'Mean Cf:    :', Sum(Cf)/Real(nx,8)

          Write(*,*) ' '          
          Write(*,*) 'Maximum U   :', maxU
          Write(*,*) 'Maximum V   :', maxV
          Write(*,*) 'Maximum W   :', maxW
          
          Write(*,*) 'Mean U      :', meanU/Real( nxm_global*nym_global*nzm_global, 8 )
          Write(*,*) 'Mean V      :', meanV/Real( nxm_global*nym_global*nzm_global, 8 )
          Write(*,*) 'Mean W      :', meanW/Real( nxm_global*nym_global*nzm_global, 8 )

          If ( LES_model>0 ) Then
             Write(*,*) ' '
             Write(*,*) "Maximum nu_t: ", maxNut
             Write(*,*) "Average nu_t: ", meanNut/Real( nxm_global*nym_global*nzm_global, 8 )
          End If

          Write(*,*) ' '
          Write(*,*) 'Mean pressure gradient in x :', dPdx
          Write(*,*) 'Mean pressure gradient in y :', dPdy
          Write(*,*) 'Spanwise rotation Omega_z   :', Omega_z

          Write(*,*) ' '
          If     ( iwall_model==1 ) Then
             Write(*,*) 'Wall model : Integral momentum balance'
          Elseif ( iwall_model==2 ) Then
             Write(*,*) 'Wall model : Boses'
          Elseif ( iwall_model==3 ) Then
             Write(*,*) 'Wall model : Robin BC fitting'
          Elseif ( iwall_model==4 ) Then
             Write(*,*) 'Wall model : Momentum balance at the wall'
          Elseif ( iwall_model==5 ) Then
             Write(*,*) 'Wall model : constant alpha'
             Write(*,*) '     Std               : ',alpha_std
             Write(*,*) '     freq. Multiplier  : ',freq_mult
          Elseif ( iwall_model==6 ) Then
             Write(*,*) 'Wall model : log-layer model (actual nu_t)'
          Elseif ( iwall_model==7 ) Then
             Write(*,*) 'Wall model : Dynamic Momentum balance at the wall'
          Elseif ( iwall_model==8 ) Then
             Write(*,*) 'Wall model : log-layer model (fixed nu_t)'
          Elseif ( iwall_model==9 ) Then
             Write(*,*) 'Wall model : imposed Cf + dynamic viscous fraction'
          Elseif ( iwall_model==11 ) Then
             Write(*,*) 'Wall model : Momentum balance at the wall with alpha_u only'
          Elseif ( iwall_model==13 ) Then
             Write(*,*) 'Wall model : imposed Cf + dynamic viscous fraction (with v)'
          Elseif ( iwall_model==15 ) Then
             Write(*,*) 'Wall model : imposed Cf with Neumann condition'
          Else
             Write(*,*) 'No wall model'
          End If

          Write(*,*) ' '
          If ( iwall_model > 0 ) Then
             If     ( istress_model==0 ) Then
                Write(*,*) 'No stress model : using Cf_ref'
                Write(*,*) '        utau wall      :',Sum( utau_wall  )/Real(nx_global,8)
                Write(*,*) '        utau wall avg  :',Sum( utau_wall_T)/Real(nx_global,8)
                Write(*,*) '        utau reference :',Sum( utau_ref   )/Real(nx_global,8)
             Elseif ( istress_model==1 ) Then
                Write(*,*) 'Law-of-the-wall model:'
                Write(*,*) '        utau wall      :',Sum( utau_wall  )/Real(nx_global,8)
                Write(*,*) '        utau model     :',Sum( utau_model )/Real(nx_global,8)
                Write(*,*) '        utau reference :',Sum( utau_ref   )/Real(nx_global,8)
             End If
          End If

          If ( iwall_model>0 ) Then
             Write(*,*) ' '
             Write(*,*) 'Wall boundary condition:'
             Write(*,*) '     beta_y            : ', beta_y
             Write(*,*) '     alpha_x     Mean  : ', Sum(    alpha_x(:,1,2) )/Real(  nx_global , 8)
             Write(*,*) '                 Max   : ', Maxval( alpha_x(:,1,2) )
             Write(*,*) '                 Min   : ', Minval( alpha_x(:,1,2) )
             Write(*,*) '     alpha_y     Mean  : ', Sum(    alpha_y(:,1,2) )/Real(  nxg_global , 8)
             Write(*,*) '                 Max   : ', Maxval( alpha_y(:,1,2) )
             Write(*,*) '                 Min   : ', Minval( alpha_y(:,1,2) )
             Write(*,*) '     alpha_z     Mean  : ', Sum(    alpha_z(:,1,2) )/Real(  nxg_global , 8)
             Write(*,*) '                 Max   : ', Maxval( alpha_z(:,1,2) )
             Write(*,*) '                 Min   : ', Minval( alpha_z(:,1,2) )
          End If
          
          Write(*,*) ' '
          Write(*,*) 'Inflow Reynolds numbers:'
          Write(*,*) '                Re_x        : ', Rex_inlet
          Write(*,*) '                Re_theta    : ', Retheta_inlet
          Write(*,*) '                Re_delta    : ', Redelta_inlet
          Write(*,*) '                delta99_i   : ', delta99_inlet_ins

          Write(*,*) ' '
          write(*,*) 'Maximum divergence          :', max_divergence
          write(*,*) 'Elapsed time (s)            :', time2-time1
          
          Write(*,*) '------------------------------------------------------'

       End If

       ! start measure time per step
       time1 = MPI_WTIME()

       Call Mpi_barrier(MPI_COMM_WORLD,ierr)

    End If

  End Subroutine output_monitor

  !------------------------------------------!
  !   Output summary of initial parameters   !
  !------------------------------------------!
  Subroutine summary

    If ( myid==0 ) Then

       Write(*,*) '------------------------------------------------------------'
       Write(*,*) '              Summary of initial parameters                 '
       Write(*,*) ' '

       Write(*,*) ' '
       Write(*,*) 'Number of processors:',nprocs

       Write(*,*) ' '
       If     ( itime_step==1 ) Then
          Write(*,*) 'Numerical integration: Explicit Euler'
       Elseif ( itime_step==2 ) Then
          Write(*,*) 'Numerical integration: Explicit RK2'
       Else
          Write(*,*) 'Numerical integration: Explicit RK3'
       End If

       Write(*,*) ' '
       Write(*,*) 'Input  file : ',Trim(filein)
       Write(*,*) 'Output file : ',Trim(fileout)

       Write(*,*) ' '
       Write(*,*) 'Inflow parameters: '
       If    ( inflow_boundary_flag == 1 ) Then
          Write(*,*) '     Blasius profile + temporal perturbations'
          Write(*,*) '     Blasius  file  : ', Trim(file_inflow)
          Write(*,*) '     temporal file  : ', Trim(file_temporal_inlet)
          Write(*,*) '     beta           : ',  beta_inlet
          Write(*,*) '     omega          : ', omega_inlet
       Elseif ( inflow_boundary_flag == 2 ) Then
          Write(*,*) '     Blasius profile from file + temporal perturbations'
          Write(*,*) '     Blasius  file2 : ', Trim(file_inflow)
          Write(*,*) '     temporal file  : ', Trim(file_temporal_inlet)
          Write(*,*) '     beta           : ',  beta_inlet
          Write(*,*) '     omega          : ', omega_inlet
       Elseif ( inflow_boundary_flag == 3 ) Then
          Write(*,*) '     Lunds rescaling'
          Write(*,*) '     rescaling plane index: ',i_rescale
          Write(*,*) '     delta 99 at inlet:     ',delta_inlet
          Write(*,*) '     time for averaging:    ',T_resc
          If ( i_rescale < 2 .or. i_rescale>nx_global-1 ) Stop 'Error, i_rescale out of range!'             
       Elseif ( inflow_boundary_flag == 4 ) Then
          Write(*,*) '     Blasius profile + temporal perturbations + random'
          Write(*,*) '     temporal file  : ', Trim(file_temporal_inlet)
          Write(*,*) '     plus random perturbations of amplitude: ', Amplitude_perturbations
       Elseif ( inflow_boundary_flag == 5 ) Then 
          Write(*,*) '     Lund fixed mean and rescaled fluctuations'
          Write(*,*) '     turbulent mean inflow file : ',Trim(file_inflow)
          Write(*,*) '     rescaling plane index:     : ',i_rescale
          !Write(*,*) '     delta 99 at inlet:         : ',delta_inlet
          !Write(*,*) '     time for averaging:        : ',T_resc
       End If
       
       Write(*,*) ' '
       Write(*,*) 'Top bc parameters: '
       If ( top_boundary_flag == 1 ) Then
          Write(*,*) '    Blasius velocity'
       Else
          Write(*,*) '    Falkner-Skan velocity with beta: ', beta_hartree
          Write(*,*) '    assuming Uinf==1 and x_origin==-1'
       End If
       
       Write(*,*) ' '
       Write(*,*) 'nu        :', nu
       Write(*,*) 'Rex inlet :', x_global(1)*1d0/nu
       Write(*,*) 'CFL       :', CFL
       If ( CFL<0 ) Write(*,*) 'fixed dt  :', dt_period
       Write(*,*) 'time      :', t

       Write(*,*) ' '
       Write(*,*) 'dPdx    :', dPdx
       Write(*,*) 'dPdz    :', dPdz
       Write(*,*) 'Omega_z :', Omega_z
       
       Write(*,*) ' '
       Write(*,*) 'nsteps   :', nsteps
       Write(*,*) 'nsave    :', nsave
       Write(*,*) 'nstats   :', nstats
       Write(*,*) 'nmonitor :', nmonitor

       Write(*,*) ' '
       Write(*,*) 'statistics z modes: '
       Write(*,*) '    nstats zmodes :', nstats_zmodes
       Write(*,*) '    i_stat, x     :', i_stat, x (i_stat)
       Write(*,*) '    j_stat, y     :', j_stat, yg(j_stat)
       Write(*,*) '    Delta_i_stat  :', Delta_i_stat
       Write(*,*) '    modes saved   :', nzu_first_modes

       Write(*,*) ' '       
       Write(*,*) 'Lx,Ly,Lz:'
       Write(*,*) Lx,Ly,Lz
       
       Write(*,*) ' '
       Write(*,*) 'nx,nxg,nxm :',nx_global,nxg_global,nxm_global
       Write(*,*) 'ny,nyg,nym :',ny_global,nyg_global,nym_global
       Write(*,*) 'nz,nzg,nzm :',nz_global,nzg_global,nzm_global
       Write(*,*) 'nz,nzg,nzm (local) :',nz,nzg,nzm
       
       Write(*,*) ' '
       Write(*,*) 'xg(1),xg(end) :',xg_global(1), xg_global(nxg_global)
       Write(*,*) 'x (1),x (end) :', x_global(1), x_global ( nx_global)
       
       Write(*,*) 'yg(1),yg(end) :',yg_global(1), yg_global(nyg_global)
       Write(*,*) 'y (1),y (end) :', y_global(1), y_global ( ny_global)
       
       Write(*,*) 'zg(1),xz(end) :',zg_global(1), zg_global(nzg_global)
       Write(*,*) 'z (1),z (end) :', z_global(1), z_global ( nz_global)

       Write(*,*) ' '
       If     ( LES_model==1 ) Then
          Write(*,*) 'LES model : constant coefficient Smagorinsky'
       Elseif ( LES_model==2 ) Then
          Write(*,*) 'LES model : dynamic Smagorinsky, z-averaged'
       Elseif ( LES_model==3 ) Then
          Write(*,*) 'LES model : dynamic Smagorinsky, no-average'
       Else
          Write(*,*) 'No LES model'
       End If

       Write(*,*) ' '
       If     ( iwall_model==1 ) Then
          Write(*,*) 'Wall model : Integral momentum balance'
       Elseif ( iwall_model==2 ) Then
          Write(*,*) 'Wall model : Boses'
       Elseif ( iwall_model==3 ) Then
          Write(*,*) 'Wall model : Robin BC fitting'
       Elseif ( iwall_model==4 ) Then
          Write(*,*) 'Wall model : Momentum balance at the wall'
       Elseif ( iwall_model==5 ) Then
          Write(*,*) 'Wall model : constant alpha'
          Write(*,*) '     Mean-x            : ',alpha_mean_x
          Write(*,*) '     Mean-y            : ',alpha_mean_y
          Write(*,*) '     Mean-z            : ',alpha_mean_z
          Write(*,*) '     Std               : ',alpha_std
          Write(*,*) '     freq. Multiplier  : ',freq_mult
       Elseif ( iwall_model==6 ) Then
          Write(*,*) 'Wall model : log-layer model (actual nu_t)'
       Elseif ( iwall_model==7 ) Then
          Write(*,*) 'Wall model : Dynamic Momentum balance at the wall'
       Elseif ( iwall_model==8 ) Then
          Write(*,*) 'Wall model : log-layer model (fixed nu_t)'
       Elseif ( iwall_model==9 ) Then
          Write(*,*) 'Wall model : imposed Cf + dynamic viscous fraction'
       Elseif ( iwall_model==11 ) Then
          Write(*,*) 'Wall model : Momentum balance at the wall with alpha_u only'
       Elseif ( iwall_model==13 ) Then
          Write(*,*) 'Wall model : imposed Cf + dynamic viscous fraction (with v)'
       Elseif ( iwall_model==15 ) Then
          Write(*,*) 'Wall model : imposed Cf with Neumann condition'
       Else
          Write(*,*) 'No wall model'
       End If

       Write(*,*) ' '
       If     ( istress_model==0 ) Then
          Write(*,*) 'No stress model : using utau_ref'
       Elseif ( istress_model==1 ) Then
          Write(*,*) 'Stress model : law-of-the-wall'
       Else
          Stop 'Error: unknown wall stress model'
       End If

       Write(*,*) ' '
       Write(*,*) 'Separation parameters (Vtop)'
       Write(*,*) '        Vbs_max  :',Vbs_max
       Write(*,*) '        sigma_bs :',sigma_bs
       Write(*,*) '        phi_bs   :',phi_bs
       Write(*,*) '        x_bs     :',x_bs
       
       Write(*,*) ' '
       Write(*,*) '------------------------------------------------------------'
    
    End If
    
  End Subroutine summary

End Module monitor
