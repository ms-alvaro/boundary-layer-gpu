!--------------------------------------!
!          Module for I/O              !
!--------------------------------------!
Module input_output

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use ifport
  Use pressure
  Use boundary_conditions, Only : apply_periodic_bc_z

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------!
  !  Read input parameters from a txt file !
  !----------------------------------------!
  Subroutine read_input_parameters

    Character(200) :: dummy_line
    Real(Int64)    :: Rossby_plus, utau_

    ! initialize default variables
    nstep_init = 0
    Rossby_plus = 0d0
    Amplitude_perturbations = 0d0
    beta_hartree = 0d0

    ! processor 0 reads the data
    If ( myid==0 ) Then
       call initparams()
      
       ! WARNING: nstats_zmodes not read

       !Read(*,*) dummy_line
       !Read(*,*) Rossby_plus

       !Read(*,*) dummy_line
       !Read(*,*) Amplitude_perturbations

       !Read(*,*) dummy_line
       !Read(*,*) beta_hartree

       ! reference utau
       utau_   = dPdx**0.5d0

       ! nominal rotation
       Omega_z = Rossby_plus*utau_/2d0
       
    End If

    ! broadcast data to all processors
    Call Mpi_bcast ( nx_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ny_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nz_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (     CFL,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nu,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    dPdx,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    dPdz,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Omega_z,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (       LES_model,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     iwall_model,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( iwall_model_nut,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( frac_vis_wall_model,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (        alpha_mean_x,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (        alpha_mean_y,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (        alpha_mean_z,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (           alpha_std,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (           freq_mult,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (         delta_inlet,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (              T_resc,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (          nstep_init,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (              nsteps,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (               nsave,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (              nstats,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nstats_zmodes,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (            nmonitor,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (         random_init,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (inflow_boundary_flag,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   top_boundary_flag,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (          itime_step,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      Dirichlet_nu_t,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (           i_rescale,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       istress_model,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (Amplitude_perturbations,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (           beta_hartree,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    !  Vbs_max, x_bs, sigma_bs, phi_bs, dummy
    Call Mpi_bcast (  Vbs_max,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     x_bs,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( sigma_bs,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   phi_bs,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    freq_mult,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    
    Call Mpi_bcast (    Lx_rand,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    Ly_rand,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    Lz_rand,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( alpha_rand,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    ! small check
    If ( i_rescale>nx_global ) Stop 'Error! i_rescale>nx_global'
    
  End Subroutine read_input_parameters

  !------------------------------------------------!
  !    Generates an initial condition for channel  !
  !                                                !
  ! Output: U,V,W                                  !
  !                                                !
  !------------------------------------------------!
  Subroutine init_flow
  
    Integer(Int32) :: ll, ii, jj, kk, i, j, k
    Real   (Int64) :: alpha

    If (myid==0) Then
       Write(*,*) 'Generating random initial condition'
    End If

    ! start time
    t = 0d0

    ! xmesh
    ! x_global(1) = x0 = 1 -> distance to the leading edge
    Do i=1,nx_global
       x_global(i) = Real(i-1,8)
    End Do
    x_global = x_global - x_global(1)
    x_global = 1d0 + x_global/Maxval(x_global)*Lx_rand 
    
    ! zmesh
    Do k=1,nz_global
       z_global(k) = Real(k-1,8)
    End Do
    z_global = z_global/z_global(nz_global-1)*Lz_rand

    ! ymesh
    ! Note: delta99/x0 = 4.91/sqrt(Rex0) = 0.0155 for Rex0 = 1e5
    Do i=1,ny_global
       y_global(i) = Real(i-1,8)
    End Do
    y_global = y_global - y_global(1)
    y_global = y_global/Maxval(y_global)
    !If ( iwall_model == 0 ) Then
       alpha = alpha_rand
       Do i=1,ny_global
          y_global(i) = dsinh(alpha*y_global(i))/dsinh(alpha)
       End Do
    !End If
    y_global = y_global - y_global(1)
    y_global = y_global/Maxval(y_global)*Ly_rand

    ! U
    U = 1d0
!    Do jj=1,ny_global/10
!       U(:,jj,:) = ( y_global(jj)/y_global(ny_global/10) )**2d0
!    end Do
    Do ii=1,nx_global
       Do jj=1,nyg_global
          Do kk=1,nzg
             U(ii,jj,kk) =  U(ii,jj,kk) + 0.1*(rand()-0.5)*( real(nyg_global-jj)/real(nyg_global) )**8d0
          End Do
       End Do
    End Do

    ! V
    V = 0d0
    Do ii=1,nxg_global
       Do jj=1,ny_global
          Do kk=1,nzg
             V(ii,jj,kk) = V(ii,jj,kk) + 0.1*(rand()-0.5)*( real(nyg_global-jj)/real(nyg_global) )**8d0
          End Do
       End Do
    End Do

    ! W
    W = 0d0
    Do ii=1,nxg_global
       Do jj=1,ny_global
          Do kk=1,nz
             W(ii,jj,kk) = W(ii,jj,kk) + 0.1*(rand()-0.5)*( real(nyg_global-jj)/real(nyg_global) )**8d0
          End Do
       End Do
    End Do

    If ( myid==0 ) Then
       Write(*,*) 'Random initial condition:'
       Write(*,*) '   Lx     : ',Lx_rand
       Write(*,*) '   Ly     : ',Ly_rand
       Write(*,*) '   Lz     : ',Lz_rand
       Write(*,*) '   alpha  : ',alpha_rand

       Write(*,*) '   Max U : ',MaxVal(U)
       Write(*,*) '   Max V : ',MaxVal(V)
       Write(*,*) '   Max W : ',MaxVal(W)
       
       Write(*,*) '   Mean U : ',sum(U)/Real(nx_global*nyg_global*nzg_global,8)
       Write(*,*) '   Mean V : ',sum(V)/Real(nxg_global*ny_global*nzg_global,8)
       Write(*,*) '   Mean W : ',sum(W)/Real(nxg_global*nyg_global*nz_global,8)
    End If
 
  End Subroutine init_flow

  !--------------------------------------------!
  !    Read binary snapshot: mesh, U,V and W   !
  !                                            !
  ! Input:  filein                             !
  ! Output: U,V,W,x,y,z                        !
  !                                            !
  !--------------------------------------------!
  Subroutine read_input_data

    Integer(Int32) ::  nx_global_f,  ny_global_f,  nz_global_f, iproc, nze, nzge
    Integer(Int32) :: nxm_global_f, nym_global_f, nzm_global_f, nn(3), ndum, nyg_resc
    Integer(Int64) :: pos_header, nsize_U, nsize_V
    Real   (Int64) :: nu_dummy
    Character(200) :: filein_resc
    logical        :: iostat_res
    
    ! processor 0 Reads the all the data
    If ( myid==0 ) Then

       Write(*,*) 'reading ',Trim(Adjustl(filein)),'...'
       Open(1,file=filein,access='stream',form='unformatted',action='Read')       

       ! metadata
       Read(1) t, nu_dummy
       
       ! mesh
       Read(1) nx_global_f
       If ( nx_global_f/=nx_global ) Stop 'nx_f/=nx'
       Read(1) x_global
       
       Read(1) ny_global_f
       If ( ny_global_f/=ny_global ) Stop 'ny_f/=ny'
       Read(1) y_global
       
       Read(1) nz_global_f
       If ( nz_global_f/=nz_global ) Stop 'nz_f/=nz'
       Read(1) z_global
       
       Read(1) nxm_global_f
       If ( nxm_global_f/=nxm_global ) Stop 'nxm_f/=nxm'
       Read(1) xm_global
       
       Read(1) nym_global_f
       If ( nym_global_f/=nym_global ) Stop 'nym_f/=nym'
       Read(1) ym_global
       
       Read(1) nzm_global_f
       If ( nzm_global_f/=nzm_global ) Stop 'nzm_f/=nzm'
       Read(1) zm_global

       ! get header position and size
       Inquire(1,pos=pos_header)
       pos_header = pos_header - 1
       nsize_U    = nx_global*nyg_global*nzg_global*8
       nsize_V    = nxg_global*ny_global*nzg_global*8
                     
    End If

    ! U
    If ( myid==0 ) Then
       ! read dummy
       Read(1) nn 
       ! read data for processor 0
       nzge = kg2_global(myid) - kg1_global(myid) + 1
       Read(1) U(:,:,1:nzge)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             ndum = fseek(1,-2*nx_global*nyg_global*8,seek_cur) ! ghost cell
             Read(1) Uo(:,:,1:nzge)
             Call Mpi_send(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
          Else ! especial case: U has different size for last processor
             ndum = fseek(1,-2*nx_global*nyg_global*8,seek_cur) ! ghost cell
             Read(1) Uoo(:,:,1:nzge)
             Call Mpi_send(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(U,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! V
    If ( myid==0 ) Then
       ! go to correct position. I dont know, if I dont do this it gets lost
       ndum = fseek(1,pos_header+3*4+nsize_U,seek_set)
       ! read dummy
       Read(1) nn
       ! read data for processor 0
       nzge = kg2_global(myid) - kg1_global(myid) + 1
       Read(1) V(:,:,1:nzge)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             ndum = fseek(1,-2*nxg_global*ny_global*8,seek_cur) ! ghost cell
             Read(1) Vo(:,:,1:nzge) 
             Call Mpi_send(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
          Else ! especial case: V has different size for last processor
             ndum = fseek(1,-2*nxg_global*ny_global*8,seek_cur) ! ghost cell
             Read(1) Voo(:,:,1:nzge) 
             Call Mpi_send(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(V,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! W
    If ( myid==0 ) Then
       ! go to correct position. I dont know, if I dont do this it gets lost
       ndum = fseek(1,pos_header+3*4+nsize_U+3*4+nsize_V,seek_set)
       ! read dummy
       Read(1) nn
       ! read data for processor 0
       nzge = k2_global(myid) - k1_global(myid) + 1
       Read(1) W(:,:,1:nzge)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             ndum = fseek(1,-2*nxg_global*nyg_global*8,seek_cur) ! ghost cell
             Read(1) Wo(:,:,1:nzge)
             Call Mpi_send(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
          Else ! especial case: W has different size for last processor
             ndum = fseek(1,-2*nxg_global*nyg_global*8,seek_cur) ! ghost cell
             Read(1) Woo(:,:,1:nzge)
             Call Mpi_send(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(W,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! close file
    If (myid==0) Then
       Close(1)
    End If

    ! read Umean and Vmean for recycling (if exists)
    If ( myid==0 .And. inflow_boundary_flag==3 ) Then
       filein_resc = Trim(Adjustl(filein))//'.mean.rescaling'
       Inquire(file = filein_resc, exist=iostat_res)
       If ( iostat_res ) Then
          Write(*,*) 'reading ',Trim(Adjustl(filein_resc)),'...'
          Open(2,file=filein_resc,access='stream',form='unformatted',action='Read',status='old',convert='big_endian')
          Read(2) nyg_resc
          If ( nyg_resc/=nyg_global ) Then 
             Write(*,*) nyg_resc, nyg_global
             Stop 'Error! nyg_resc/=nyg_global'
          End If
          Read(2) Umean_resc_To
          Read(2) Vmean_resc_To
          Close(2)
       Else
          Umean_resc_To = 0d0
          Vmean_resc_To = 0d0
       End If
    End If
    Call Mpi_bcast ( Umean_resc_To,nyg_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Vmean_resc_To, ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    ! send data to all other processors
    ! mesh
    Call Mpi_bcast ( x_global,nx_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( y_global,ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( z_global,nz_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( xm_global,nxm_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ym_global,nym_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( zm_global,nzm_global,MPI_real8,0,MPI_COMM_WORLD,ierr )   

    Call Mpi_bcast ( t,1,MPI_real8,0,MPI_COMM_WORLD,ierr ) 

    ! set solution for zero step
    Uo = U
    Vo = V
    Wo = W

    ! Sanity check
    If ( Any( Isnan(U) ) ) Stop 'Error U NaNs!'
    write (*,*) 'no nans here'
    If ( Any( Isnan(V) ) ) Stop 'Error V NaNs!'
    If ( Any( Isnan(W) ) ) Stop 'Error W NaNs!'
    
  End Subroutine read_input_data

  !--------------------------------------------!
  !    write binary snapshot: mesh, U,V and W  !
  !                                            !
  ! Input: U,V,W,x,y,z,xm,ym,zm                !
  ! Output: fileout                            !
  !                                            !
  !--------------------------------------------!
  Subroutine output_data

    Character(200)   :: fname, fileout_resc
    Character(8)     :: ext
    Integer  (Int32) :: iproc, nze, nzge
    
    If ( Mod(istep,nsave)==0 ) then

       ! P
       !If (pressure_computed==.False.) Then
       !   Call compute_pressure
       !End If

       ! processor 0 writes the data
       If ( myid==0 ) Then
          
          Write(ext,'(I0.8)') istep + nstep_init
          
          fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))
          Write(*,*) 'writting ',Trim(Adjustl(fname))
          Open(1,file=fname,access='stream',form='unformatted',action='write')
          
          ! metadata
          Write(1) t, nu

          ! mesh
          Write(1) Shape(x_global), x_global
          Write(1) Shape(y_global), y_global
          Write(1) Shape(z_global), z_global
          
          Write(1) Shape(xm_global), xm_global
          Write(1) Shape(ym_global), ym_global
          Write(1) Shape(zm_global), zm_global          
         
       End If

       ! U
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(U,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
       Else
          ! write U size
          Write(1) nx_global,nyg_global,nzg_global
          ! processor 0 writes its data
          Write(1) U(:,:,1:nzg-1) 
          ! processor 0 receives and writes rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Uo(:,:,2:nzge-1)
             Else
                Call Mpi_recv(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Uoo(:,:,2:nzge)
             End If
          End Do
       Endif

       ! V
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(V,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
       Else
          ! write V size
          Write(1) nxg_global,ny_global,nzg_global
          ! processor 0 writes its data
          Write(1) V(:,:,1:nzg-1)
          ! processor 0 receives and write rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Vo(:,:,2:nzge-1)
             Else
                Call Mpi_recv(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Voo(:,:,2:nzge)
             End If
          End Do
       Endif

       ! W
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(W,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
       Else
          ! write W size
          Write(1) nxg_global,nyg_global,nz_global
          ! processor 0 writes its data
          Write(1) W(:,:,1:nz-1)
          ! processor 0 receives and writes rest data
          Do iproc = 1, nprocs-1
             nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Wo(:,:,2:nze-1)
             Else
                Call Mpi_recv(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Woo(:,:,2:nze)
             End If
          End Do
       Endif

       ! nu_t
       If ( myid/=0 ) Then
          ! data from processor n>0
          Call Mpi_send(nu_t,nxg*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
       Else
          ! write nu_t size
          Write(1) nxg_global,nyg_global,nzg_global
          ! processor 0 writes its data
          Write(1) nu_t(:,:,1:nzg-1)
          ! processor 0 receives and write rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(nu_to,nxg*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) nu_to(:,:,2:nzge-1)
             Else
                Call Mpi_recv(nu_too,nxg*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) nu_too(:,:,2:nzge)
             End If
          End Do
       Endif

       ! P (to be checked)
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(P,nxg*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
       Else
          ! write P size
          Write(1) nxg_global,nyg_global,nzg_global
          ! processor 0 writes its data
          Write(1) P
          ! processor 0 receives and write rest data
          Do iproc = 1, nprocs-1
             nze = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Po,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Po
             Else
                Call Mpi_recv(Po,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                Write(1) Po
             End If
          End Do
       Endif
                    
       ! close file
       If (myid==0) Then
          Close(1)
       End If

       ! save means for Lund's rescaling
       If ( myid==0 .And. inflow_boundary_flag==3 ) Then
          fileout_resc = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.mean.rescaling'
          Open(3,file=fileout_resc,access='stream',form='unformatted',action='Write')
          Write(3) nyg_global
          Write(3) Umean_resc_T
          Write(3) Vmean_resc_T
          Close(3)
       End If
       
    End If
       
  End Subroutine output_data

  !----------------------------------------------!
  !   Write some basic statistics in a txt file  !
  !----------------------------------------------!
  Subroutine output_statistics

    Character(200) :: fname,my_format
    Character(8)   :: ext,ext_ny,ext_nx
    Integer(Int32) :: ii

    If ( myid==0 ) Then

       ! create file name
       Write(ext,'(I0.8)') istep + nstep_init

       Write(ext_nx,'(I8)') nx
       Write(ext_ny,'(I8)') ny
              
       fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.stats.txt'
       Write(*,*) 'writting ',Trim(Adjustl(fname))
       Open(3,file=fname,form='formatted',action='write') 
       !
       Write(3,'(A,2F15.8,4I)') '%',t, nu, nx_global, ny_global, nz_global, istep

       my_format = '('//Trim(Adjustl(ext_nx))//'F15.8)'
       Write(3,my_format) Cf

       my_format = '('//Trim(Adjustl(ext_ny))//'F15.8)'
       Write(3,my_format) y

       !
       Do ii=1,nx
          Write(3,my_format) Umean(ii,:)
       End Do
       Do ii=1,nx
          Write(3,my_format) Vmean(ii,:)
       End Do
       Do ii=1,nx
          Write(3,my_format) Wmean(ii,:)
       End Do
       !
       Do ii=1,nx
          Write(3,my_format) U2mean(ii,:)
       End Do
       Do ii=1,nx
          Write(3,my_format) V2mean(ii,:)
       End Do
       Do ii=1,nx
          Write(3,my_format) W2mean(ii,:)
       End Do
       !
       Do ii=1,nx
          Write(3,my_format) Pmean(ii,:)
       End Do
       Do ii=1,nx
          Write(3,my_format) P2mean(ii,:)
       End Do
       !
       Close(3)

    End If

  End Subroutine output_statistics

  !----------------------------------------------!
  !   Write zmodes statistics on a binary file   !
  !----------------------------------------------!
  Subroutine output_statistics_zmodes

    Character(200) :: fname
    Character(8)   :: ext

    If ( myid==0 ) Then

       ! create file name
       Write(ext,'(I8)') istep + nstep_init
       fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.zmodes.dat'
       Write(*,*) 'writting ',Trim(Adjustl(fname))
       
       ! write data
       Open(6,file=fname,action='write',form='unformatted',access='stream')
       Write(6) t, beta_inlet, omega_inlet
       Write(6) nxu_reduced, nyu_reduced, nzu_reduced, nzu_first_modes
       Write(6) xu_reduced, yu_reduced, zu_reduced
       Write(6) Real (U_reduced_hat_z,8)
       Write(6) Dimag(U_reduced_hat_z)
       Write(6) Real (rhs_u_reduced_hat_z,8)
       Write(6) Dimag(rhs_u_reduced_hat_z)
       Close(6)

    End If

  End Subroutine output_statistics_zmodes

End Module input_output
