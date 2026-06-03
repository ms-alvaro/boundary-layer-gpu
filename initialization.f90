!------------------------------------------------!
! Module with initialization of global variables !
!------------------------------------------------!
Module initialization

  ! Modules
  Use, Intrinsic :: iso_c_binding
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use input_output
  Use fftz
  Use cufft_solver
  Use boundary_conditions, Only : apply_periodic_bc_z
 
  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  !----------------------------------------!
  !         Initialize everything          !
  !----------------------------------------!
  Subroutine initialize
    
    Integer(Int32) :: i, j, k, kk, nzpe, pos, ipos, nze, nzme
    Real   (Int64) :: dy1, dy2, det, a, b, c, r, alpha_o_local, alpha_o_global
    Integer(Int32), Dimension(:,:), Allocatable :: A_kmodes, A_kmodes_local

    !---------------------first initialize MPI-------------------!
    call Mpi_init(ierr)
    call Mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call Mpi_comm_rank(MPI_COMM_WORLD,   myid, ierr)

    If (myid==0) Then 
       Write(*,*) '----------------------------------------------------------------------'       
       Write(*,*) ' '
       Write(*,*) '    My boundary layer ^^, parallel version 0.7.2 with separation      '
       Write(*,*) ' '
       Write(*,*) '----------------------------------------------------------------------'
    End If

    !--------------read parameters from standard input-----------!
    If ( myid==0 ) Write(*,*) 'reading parameters...' 
    Call read_input_parameters

    ! time
    ! now read from input
    
    !-------------------grid definitions-------------------------!
    Allocate (  k1_global(0:nprocs-1),  k2_global(0:nprocs-1) )
    Allocate ( kg1_global(0:nprocs-1), kg2_global(0:nprocs-1) )

    ! restrictions for FFTW-MPI mapping (only enforced for multi-rank)
    If ( nprocs > 1 ) Then
       If ( Mod( nx_global, 2 )/=0 )      Stop 'Error: nx must be even for MPI'
       If ( Mod( nz_global, 2 )/=0 )      Stop 'Error: nz must be even for MPI'
       If ( Mod( nz_global-2, nprocs)/=0 ) Stop 'nz-2 should be divisible by nprocs'
    End If

    ! number of interior z-planes per processor based on fftw decomposition
    nslices_z = Nint( Real((nz_global-2))/Real(nprocs) ) 
    
    ! restriction for MPI boundaries
    !If ( nslices_z<2 ) Stop 'Error: nslices_z must be at least 2' 
    If ( nslices_z<1 ) Stop 'Error: nslices_z must be at least 1' 

    ! domain decomposition. Must be consistent with fftw
    Do i = 0, nprocs-1
       ! range index for faces in each processor
       k1_global(i)  = i*nslices_z  + 1
       k2_global(i)  = k1_global(i) + nslices_z + 1
       ! range index for centers in each processor
       kg1_global(i) = i*nslices_z   + 1
       kg2_global(i) = kg1_global(i) + nslices_z + 1
    End Do    

    ! remaining planes in last processor
    k2_global (nprocs-1) = nz_global 
    kg2_global(nprocs-1) = nz_global + 1

    ! face points
    nx = nx_global
    ny = ny_global
    nz = k2_global(myid) - k1_global(myid) + 1 

    ! middle points
    nxm_global = nx_global - 1
    nym_global = ny_global - 1
    nzm_global = nz_global - 1

    nxm = nx - 1
    nym = ny - 1
    nzm = kg2_global(myid) - kg1_global(myid) + 1 - 2  
    
    ! middle points + ghost cells
    nxg_global = nxm_global + 2
    nyg_global = nym_global + 2
    nzg_global = nzm_global + 2

    nxg = nxm + 2
    nyg = nym + 2
    nzg = kg2_global(myid) - kg1_global(myid) + 1 

    ! size for last proccesor nz and nzm -> nze and nzme
    nze  = nz
    nzme = nzm
    Call Mpi_bcast (  nze,1,MPI_integer,nprocs-1,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nzme,1,MPI_integer,nprocs-1,MPI_COMM_WORLD,ierr )
   
    ! Allocate main arrays
    If ( myid==0 ) Write(*,*) 'allocating...' 
    Allocate ( x_global (  nx_global),  y_global (  ny_global),  z_global (  nz_global)  )
    Allocate ( xm_global( nxm_global),  ym_global( nym_global),  zm_global( nzm_global)  )
    Allocate ( xg_global(nxm_global+2), yg_global(nym_global+2), zg_global(nzm_global+2) )

    Allocate (  x (  nx),  y (  ny),  z (  nz) )
    Allocate (  xm( nxm),  ym( nym),  zm( nzm) )
    Allocate ( xg(nxm+2), yg(nym+2), zg(nzm+2) )

    Allocate ( yg_m (nyg-1) )
    Allocate ( yg_mm(nyg-2) )
    
    ! global interior + boundary + ghost points
    Allocate (U (    nx, nym+2, nzm+2) )
    Allocate (V ( nxm+2,    ny, nzm+2) )
    Allocate (W ( nxm+2, nym+2,    nz) )
    Allocate (P ( nxm+2, nym+2, nzm+2) )

    Allocate (Uo  (    nx,  nym+2, nzm+2) )
    Allocate (Vo  ( nxm+2,     ny, nzm+2) )
    Allocate (Wo  ( nxm+2,  nym+2,    nz) )
    Allocate (Po  ( nxm+2,  nym+2, nzm+2) )

    Allocate (Uoo (    nx,  nym+2, nzme+2) ) ! z-planes modified for I/O
    Allocate (Voo ( nxm+2,     ny, nzme+2) )
    Allocate (Woo ( nxm+2,  nym+2,    nze) )

    ! Auxiliary arrays
    Allocate ( term   ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_1 ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_2 ( nxg, nyg, nzm+2 ) ) 
    If (iwall_model>0) Then
       Allocate ( term_3 ( nxg, nyg, nzm+2 ) ) 
       Allocate ( term_4 ( nxg, nyg, nzm+2 ) ) 
    End If
    Allocate (V_bottom( nxm+2, nzm+2) )

    ! RHS: interior points only
    Allocate ( rhs_uo ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) ) 
    Allocate ( rhs_vo ( 2:nxg-1, 2:ny-1,  2:nzg-1 ) )
    Allocate ( rhs_wo ( 2:nxg-1, 2:nyg-1, 2:nz-1  ) )
    Allocate ( rhs_p  ( 2:nxg-1, 2:nyg-1, 2:nzg   ) ) ! ONE EXTRA PLANE IN Z FOR GHOST CELL
    If (.False.) Then ! used for computing pressure from eq. 
       Allocate ( rhs_uf  ( 1:nx  ,  1:nyg  , 1:nzg   ) ) 
       Allocate ( rhs_vf  ( 1:nxg  , 1:ny  ,  1:nzg   ) )
       Allocate ( rhs_wf  ( 1:nxg  , 1:nyg  , 1:nz    ) )
    End If

    ! Lund's rescaling
    Allocate(Umean_resc_To(nyg,1),Vmean_resc_To(ny,1))
    Umean_resc_To = 0d0
    Vmean_resc_To = 0d0
    
    ! read data 
    If ( myid==0 ) Write(*,*) 'initializing initial condition...' 
    If ( random_init==1 ) Then
       Call init_flow
    Else
       Call read_input_data
    End If

    ! definie global grids from x_global, y_global and z_global (face to centers)
    ! local faces
    x = x_global
    y = y_global
    z = z_global( k1_global(myid):k2_global(myid) )

    ! global interior centers
    Do i = 1, nxm_global
       xm_global(i) = 0.5d0*( x_global(i) + x_global(i+1) )
    End Do
    Do j=1,nym_global
       ym_global(j) = 0.5d0*( y_global(j) + y_global(j+1) )
    End Do
    Do k=1,nzm_global
       zm_global(k) = 0.5d0*( z_global(k) + z_global(k+1) )
    End Do

    ! local interior centers
    xm = xm_global
    ym = ym_global
    zm = zm_global( kg1_global(myid):kg2_global(myid)-2 )

    ! global 
    xg_global(2:nxm_global+1) = xm_global
    xg_global(1)              = xm_global(1)          - 2d0*(xm_global(1)-x_global(1))
    xg_global(nxm_global+2)   = xm_global(nxm_global) + 2d0*(x_global(nx_global)-xm_global(nxm_global))
    
    yg_global(2:nym_global+1) = ym_global
    yg_global(1)              = ym_global(1)          - 2d0*(ym_global(1)-y_global(1))
    yg_global(nym_global+2)   = ym_global(nym_global) + 2d0*(y_global(ny_global)-ym_global(nym_global))
    
    zg_global(2:nzm_global+1) = zm_global
    zg_global(1)              = zm_global(1)          - 2d0*(zm_global(1)-z_global(1))
    zg_global(nzm_global+2)   = zm_global(nzm_global) + 2d0*(z_global(nz_global)-zm_global(nzm_global))    

    xg = xg_global
    yg = yg_global
    zg = zg_global( kg1_global(myid):kg2_global(myid) )

    ! middle points for yg (.not. equal to y in general)
    yg_m = 0.5d0*( yg(2:nyg) + yg(1:nyg-1) )

    ! middle points for yg_m (.not. equal to ym in general)
    yg_mm = 0.5d0*( yg_m(2:nyg-1) + yg_m(1:nyg-2) )

    ! local minimum grid size for CFL
    dxmin = Minval ( xg_global(2:nxg_global) - xg_global(1:nxg_global-1) )
    dymin = Minval ( yg_global(2:nyg_global) - yg_global(1:nyg_global-1) )
    dzmin = Minval ( zg_global(2:nzg_global) - zg_global(1:nzg_global-1) )

    ! total domain size
    Lx = x_global(nx_global) - x_global(1)
    Ly = y_global(ny_global) - y_global(1)
    Lz = z_global(nz_global) - z_global(1)

    !--------------------------Boundary conditions--------------------------!
    ! local velocity, initial z-planes
    Allocate ( buffer_ui(nx,nyg,2:3), buffer_vi(nxg,ny,2:3), buffer_wi(nxg,nyg), buffer_ci(nxg,nyg,2:3) )
    ! local velocity, ending  z-planes
    Allocate ( buffer_ue(nx,nyg),     buffer_ve(nxg,ny),     buffer_we(nxg,nyg), buffer_ce(nxg,nyg) )
    ! local pressure z-plane
    Allocate ( buffer_p(2:nxg-1,2:nyg-1) ) 

    !------------------------Interior communications------------------------!
    Allocate ( buffer_us(nx ,nyg), buffer_ur(nx ,nyg) )
    Allocate ( buffer_vs(nxg, ny), buffer_vr(nxg, ny) )
    Allocate ( buffer_ws(nxg,nyg), buffer_wr(nxg,nyg) )
    Allocate ( buffer_ps(2:nxg-1,2:nyg-1), buffer_pr(2:nxg-1,2:nyg-1) ) 

    ! force periodicity in z in case initial condition is not
    Call apply_periodic_bc_z(U,1)
    Call apply_periodic_bc_z(V,2)
    Call apply_periodic_bc_z(W,3)

    !---------------------------Fourier transform---------------------------!
    ! initialize MPI FFTW
    If ( myid==0 ) Write(*,*) 'initializing FFT...' 
    Call fftw_mpi_init()

    ! Fourier constant grid spacing
    dx = dxmin
    dz = dzmin

    ! length for periodic domain
    Lxp = Lx - dx 
    Lzp = Lz - dz

    ! global points for periodic domain in physical space
    nxp_global = nxm_global - 1
    nzp_global = nzm_global - 1

    ! global indices for fourier modes starting from 0
    mx_global = nxp_global - 1
    mz_global = nzp_global - 1

    ! Get local sizes:
    ! local data size in x direction
    nxp = nxp_global
    mx  =  mx_global
    ! extended values. Factor of 4 for cosine transform
    nxpe        = 4*nxp
    nxpe_global = 4*nxp_global
    ! local data size in z direction (note dimension reversal)
    alloc_local = fftw_mpi_local_size_2d(nzp_global, nxpe_global, MPI_COMM_WORLD, nzp, local_k_offset)
    mz  = nzp - 1

    ! sanity check and restrictions in fftw
    If ( (nzp/=nzm .And. myid/=nprocs-1) .Or. (nzp/=nzm-1 .And. myid==nprocs-1) ) Then 
       Write(*,*) nzp,nzm
       Stop 'Error: something wrong in FFTW size'
    End If

    ! allocate variables
    cplane_fft = fftw_alloc_complex(alloc_local)
    Call c_f_pointer(cplane_fft,plane,[nxpe,nzp])
    plane_hat(0:,0:) => plane
    Allocate ( rhs_p_hat ( 0:mx, 2:nyg-1, 0:mz ) )
    Allocate ( rhs_aux   ( 2:nyg-1 ) )
    Allocate ( plane_short(nxp,nzp) )
   
    ! create MPI plan for forward DFT (note dimension reversal and transposed_out/in and x4 for cosine transform)
    ! uses imode_map_fft and kmode_map_fft
!    plan_d = fftw_mpi_plan_dft_2d( nzp_global, nxpe_global, plane, plane_hat,          & 
!             MPI_COMM_WORLD,  FFTW_FORWARD, ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_OUT) ) 
!    plan_i = fftw_mpi_plan_dft_2d( nzp_global, nxpe_global, plane_hat, plane,          & 
!             MPI_COMM_WORLD, FFTW_BACKWARD, ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_IN)  ) 

    ! create MPI plan for forward DFT (note dimension reversal and x4 for cosine transform)
    ! uses imode_map, kmode_map
    ! NOTE: use plane for both in/out (in-place plans) to avoid nvfortran
    !       pointer alias bug with plane_hat(0:,0:) => plane
    plan_d = fftw_mpi_plan_dft_2d( nzp_global, nxpe_global, plane, plane, &
             MPI_COMM_WORLD,  FFTW_FORWARD, FFTW_ESTIMATE )
    plan_i = fftw_mpi_plan_dft_2d( nzp_global, nxpe_global, plane, plane, &
             MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )

    ! Initialize cuFFT plans for single-rank GPU solver
    If ( nprocs == 1 ) Then
       ! Batched cuFFT: all nyg-2 interior y-planes at once
       Call cufft_init_plans(Int(nzp_global), Int(nxpe_global), Int(nyg-2))
       Allocate( plane_gpu(nxpe, nzp, nyg-2) )
       Allocate( rhs_hat_gpu(mx+1, nyg-2, mz+1) )
       If (myid==0) Write(*,*) 'cuFFT batched plans created:', nyg-2, 'y-planes'
    Else
       Allocate( plane_gpu(1, 1, 1) )  ! dummy
       Allocate( rhs_hat_gpu(1, 1, 1) )  ! dummy
    End If

    ! global Fourier coeficients with modified wave-number for the second derivative
    Allocate ( kxx(0:mx_global), kzz(0:mz_global) ) 

    ! modified wave-number for cosine transform
    kxx = 0d0
    Do i = 0, mx_global
       kxx(i) = 2d0*( dcos(pi*Real(i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0  
    End do

    ! modified wave-number for Fourier transform
    kzz = 0d0
    Do k = 0, Ceiling( Real(nzp_global)/2d0 )
       kzz(k) = 2d0*( dcos(2d0*pi*Real(k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0  
    End do
    Do k = Ceiling( Real(nzp_global)/2d0 )+1, mz_global
       kzz(k) = 2d0*( dcos(2d0*pi*Real(-nzp_global+k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
    End do
    
    ! MPI mapping for z-modes: from local to global without transposed_out/in
    Allocate ( imode_map(0:mx) ) 
    Allocate ( kmode_map(0:mz) ) 
    imode_map = 0
    kmode_map = 0
    Do i = 0, mx
       imode_map(i) = i
    End Do
    Do k = 0, mz
       kmode_map(k) = k + myid*nslices_z
    End Do

    ! FFTW+MPI mapping for x and z-modes when using FFTW with transposed_out/in
    ! from local to global
    ! this needs (mz_global+1)*(mx_global+1)/nprocs to be an integer 
    If ( Mod((mz_global+1)*(mx_global+1),nprocs)/=0 ) Stop 'Error: (mz_global+1)*(mx_global+1)/nprocs should be an integer'
    Allocate ( imode_map_fft(0:mx_global,0:mz) ) 
    Allocate ( kmode_map_fft(0:mx_global,0:mz) ) 
    Do i = 0, mx_global
       Do k = 0, mz            
          pos = i + (mx_global+1)*k + (mz_global+1)*(mx_global+1)/nprocs*myid
          imode_map_fft(i,k) = Floor( Real(pos/(mz_global+1)) )
          kmode_map_fft(i,k) = Mod  ( pos, mz_global+1 )
          ! sanity check
       end Do
    End Do

    ! Sanity check for FFTW mapping
    Allocate(A_kmodes      (0:mx_global,0:mz_global))
    Allocate(A_kmodes_local(0:mx_global,0:mz_global))
    A_kmodes       = 0
    A_kmodes_local = 0
    Do i = 0, mx_global
       Do k = 0, mz            
          A_kmodes_local( imode_map_fft(i,k), kmode_map_fft(i,k) ) =  A_kmodes_local(imode_map_fft(i,k), kmode_map_fft(i,k) ) + 1
       end Do
    End Do
    Call MPI_AllReduce(A_kmodes_local,A_kmodes,(mx_global+1)*(mz_global+1),MPI_integer,MPI_sum,MPI_COMM_WORLD,ierr)
    If ( Any(A_kmodes>1) .Or. Any(A_kmodes==0) ) Stop 'Error: wrong combination of nx, nz and processors'
    Deallocate(A_kmodes)
    Deallocate(A_kmodes_local)
     
    !------------------------Tridiagonal linear solver-------------------------!
    If ( myid==0 ) Write(*,*) 'initializing pressure solver...'
    Allocate ( pivot(nyg) )
    Allocate ( Dyy(2:nyg-1,2:nyg-1), M(2:nyg-1,2:nyg-1) )
    Allocate ( D(2:nyg-1), DL(2:nyg-2), DU(2:nyg-2) )
     
    ! second derivative matrix for pressure (full data in y assumed)
    Dyy = 0d0
    Do j=3,nyg-2

       a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
       b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
       c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 

       Dyy(j,j+1) = a
       Dyy(j,j-1) = c 
       Dyy(j,j  ) = b

    End Do

    ! Boundary conditions for pressure (full data in y assumed)
    j = 2
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    ! Neumann in V: p(1)==p(2) 
    Dyy(2,2)   = b + c 
    Dyy(2,3)   = a
    coef_bc_1  = c

    j = nyg-1
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    ! Neumann in V: p(nyg)==p(nyg-1) 
    Dyy(nyg-1,nyg-1) = a + b
    Dyy(nyg-1,nyg-2) = c
    coef_bc_2        = a

    ! Precompute Thomas LU factorization for all (i,k) modes
    If ( nprocs == 1 ) Then
       Block
          Integer(Int32) :: ii, kk, jj, nm_loc, ig, kg
          Complex(Int64) :: dd
          nm_loc = Int(nyg) - 2
          ! Transposed layout: (mx+1, mz+1, nm) for coalesced GPU access
          Allocate( thomas_dl_fact(mx+1, mz+1, nm_loc) )
          Allocate( thomas_d_pivot(mx+1, mz+1, nm_loc) )
          Allocate( thomas_du(nm_loc) )
          ! DU is the same for all modes
          Do jj = 1, nm_loc-1
             thomas_du(jj) = dcmplx(Dyy(jj+1, jj+2))
          End Do
          thomas_du(nm_loc) = (0d0, 0d0)
          ! Factorize each (i,k) mode (stored transposed)
          Do kk = 1, Int(mz)+1
             Do ii = 1, Int(mx)+1
                ig = imode_map(ii-1)
                kg = kmode_map(kk-1)
                Do jj = 1, nm_loc
                   thomas_d_pivot(ii, kk, jj) = dcmplx(Dyy(jj+1, jj+1) + kxx(ig) + kzz(kg))
                End Do
                Do jj = 1, nm_loc-1
                   thomas_dl_fact(ii, kk, jj) = dcmplx(Dyy(jj+2, jj+1))
                End Do
                thomas_dl_fact(ii, kk, nm_loc) = (0d0, 0d0)
                If ( ig==0 .And. kg==0 ) thomas_d_pivot(ii, kk, 1) = 3d0/2d0 * thomas_d_pivot(ii, kk, 1)
                ! Forward elimination (factorize)
                Do jj = 2, nm_loc
                   dd = thomas_dl_fact(ii, kk, jj-1) / thomas_d_pivot(ii, kk, jj-1)
                   thomas_dl_fact(ii, kk, jj-1) = dd
                   thomas_d_pivot(ii, kk, jj) = thomas_d_pivot(ii, kk, jj) - dd * thomas_du(jj-1)
                End Do
                ! Invert D_pivot for multiply instead of divide at runtime
                Do jj = 1, nm_loc
                   thomas_d_pivot(ii, kk, jj) = (1d0, 0d0) / thomas_d_pivot(ii, kk, jj)
                End Do
             End Do
          End Do
          If (myid==0) Write(*,*) 'Thomas LU precomputed:', nm_loc, 'x', Int(mx)+1, 'x', Int(mz)+1
       End Block
    End If

    Allocate ( bc_1(2:nxg-1,2:nzg-1), bc_2(2:nxg-1,2:nzg-1) )
    Allocate ( bc_1_hat(0:mx,0:mz),   bc_2_hat(0:mx,0:mz)   )

    ! some parameters for linear solver
    nr   = nym
    nrhs = 1

    !--------------------interpolation weights--------------------!
    in1 = 1
    in2 = 2 
    If ( in2==1 ) Write(*,*) 'Conservative interpolations'
    Allocate ( weight_y_0(ny), weight_y_1(ny) )
    weight_y_0 = ( yg(2:nyg) - y(1:ny) ) / ( yg(2:nyg) - yg(1:nyg-1)  )
    weight_y_1 = 1d0 - weight_y_0

    !------------------------statistics---------------------------!
    Allocate ( Cf(nx), dUdy_wall(nx), UV_wall(nx) )
    Allocate ( Uaux_1(nx), Uaux_2(nx), Uaux_1_local(nx), Uaux_2_local(nx) )
    Allocate (  Umean(nx,ny),  Vmean(nx,ny),  Wmean(nx,ny),  Pmean(nx,ny)                )
    Allocate ( U2mean(nx,ny), V2mean(nx,ny), W2mean(nx,ny), UVmean(nx,ny), P2mean(nx,ny) )    
    Allocate ( nu_t_mean(nxg,nyg) )
    Umean     = 0d0
    Vmean     = 0d0
    Wmean     = 0d0
    Pmean     = 0d0
    U2mean    = 0d0
    V2mean    = 0d0
    W2mean    = 0d0
    UVmean    = 0d0
    P2mean    = 0d0
    nu_t_mean = 0d0

    !------------------------Runge-Kutta 2-------------------------!
    Allocate( rk2_coef(2,2), rk2_t(0:2) ) 

    rk2_t(0)      =  0d0
    rk2_t(1)      =  1d0/2d0
    rk2_t(2)      =  1d0

    rk2_coef      =  0d0
    rk2_coef(1,1) =  1d0/2d0
    rk2_coef(2,1) =  0d0
    rk2_coef(2,2) =  1d0

    !------------------------Runge-Kutta 3-------------------------!
    If ( myid==0 ) Write(*,*) 'initializing time integration...' 
    Allocate( rk_coef(3,3), rk_t(0:3) ) 

    rk_t(0)      =  0d0
    rk_t(1)      =  8d0/15d0
    rk_t(2)      =  2d0/3d0
    rk_t(3)      =  1d0

    rk_coef      =  0d0
    rk_coef(1,1) =  8d0/15d0
    rk_coef(2,1) =  1d0/4d0
    rk_coef(2,2) =  5d0/12d0
    rk_coef(3,1) =  1d0/4d0
    rk_coef(3,2) =  0d0
    rk_coef(3,3) =  3d0/4d0

    Allocate ( Fu1 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu2 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    If (itime_step==3) Then
       Allocate ( Fu3 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    End If

    Allocate ( Fv1 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv2 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    If (itime_step==3) Then
       Allocate ( Fv3 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    End If

    Allocate ( Fw1 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw2 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    If (itime_step==3) Then
       Allocate ( Fw3 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    End If

    !----------------------sgs model------------------------------!

    ! Interior points only
    If (LES_model>0) Then
       Allocate( Lij    (2:nxg  ,2:nyg-1,2:nzm+2 ,6) ) ! These need one more points in each periodic direction (for filtering)
       Allocate( Mij    (2:nxg-1,2:nyg-1,2:nzm+1, 6) )
       Allocate( Sij    (2:nxg  ,2:nyg  ,2:nzm+2, 6) ) ! These need one more points in each direction (for filtering)
       Allocate( S      (2:nxg-1,2:nyg-1,2:nzm+1)    )

       ! Tensor buffer
       Allocate( ten_buf(1:nxg  ,1:nyg  ,1:nzm+2  ,6) )

       ! Filtered velocities
       Allocate( Uf (1:nx ,1:nyg,1:nzm+2) )
       Allocate( Vf (1:nxg,1:ny ,1:nzm+2) )
       Allocate( Wf (1:nxg,1:nyg,1:nz ) )
       
       Allocate( Uff (1:nx ,1:nyg,1:nzm+2) )
       Allocate( Vff (1:nxg,1:ny ,1:nzm+2) )
       Allocate( Wff (1:nxg,1:nyg,1:nz ) )
    End If
       
    ! Eddy viscosity
    Allocate( nu_t        (1:nxg , 1:nyg, 1:nzg ) )
    Allocate( avg_nu_t    (1:nxg , 1:nyg, 1     ) )
    Allocate( avg_nu_t_hat(1     , 1:nyg, 1     ) )
    Allocate( nu_to       ( nxm+2, 1:nyg, nzm+2 ) )
    Allocate( nu_too      ( nxm+2, 1:nyg, nzme+2) )

    ! Never use Dirichlet for nu_t
    Dirichlet_nu_t = 0
 
    !--------------------Integral Wall-model-------------------!

    ! ui = alpha_i dui/dy
    Allocate( alpha_x(1:nx ,1:2,1:nzg) )
    Allocate( alpha_y(1:nxg,1:2,1:nzg) )
    Allocate( alpha_z(1:nxg,1:2,1:nz ) )    

    Allocate( alpha_xo(1:nx ,1:2,1:nzg) )
    Allocate( alpha_yo(1:nxg,1:2,1:nzg) )
    Allocate( alpha_zo(1:nxg,1:2,1:nz ) ) 

    Allocate( dUmean_wall_T(nx_global), mtau_T(nx_global), UV_wall_T(nx_global) ) 

    ! estimate initial alphas for time average
    Do i=1,nx_global
       alpha_o_local = Sum( 0.5d0*( U(i,2,:) + U(i,1,:) )/( (U(i,2,:) - U(i,1,:))/(yg(2)-yg(1)) ) )
       Call MPI_Allreduce(alpha_o_local, alpha_o_global,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       alpha_xo(i,:,:) = alpha_o_global/Real( nzg_global , 8)
    End Do
    Do i=1,nxg_global
       alpha_o_local = Sum( V(i,1,:) /( (V(i,2,:) - V(i,1,:))/(y(2)-y(1)) ) )
       Call MPI_Allreduce(alpha_o_local, alpha_o_global,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       alpha_yo(i,:,:) = alpha_o_global/Real( nzg_global , 8)

       alpha_o_local = Sum( 0.5d0*( W(i,2,:) + W(i,1,:) )/( (W(i,2,:) - W(i,1,:))/(yg(2)-yg(1)) ) )
       Call MPI_Allreduce(alpha_o_local, alpha_o_global,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       alpha_zo(i,:,:) = alpha_o_global/Real( nz_global , 8)
    End Do

    ! always activate penetration (V/=0 at the wall)
    penetration = 1 

    !------Blasius/turbulent solution for boundary conditions for bc-----!    
    If ( inflow_boundary_flag == 3 .or. inflow_boundary_flag == 5 ) Then
       If ( myid==0 ) Write(*,*) 'initializing turbulent boundary conditions...' 
       Call compute_turbulent_solution_for_bc
    Else
       If ( myid==0 ) Write(*,*) 'initializing blasius boundary conditions...'
       Call compute_blasius_solution_for_bc
    End If

    If ( myid==0 ) Write(*,*) 'initializating some statistics...'
    !--------------Means for Lund's rescaling inflow-----------!    
    Allocate(Umean_resc_T       (nyg,1), Vmean_resc_T       (ny,1))
    Allocate(Umean_resc_T_local (nyg,1), Vmean_resc_T_local (ny,1))
    Allocate(Umean_inlet_T      (nyg,1), Vmean_inlet_T      (ny,1))
    Allocate(Umean_inlet_T_local(nyg,1), Vmean_inlet_T_local(ny,1))

    If ( Sum(Umean_resc_To)>1e-4 ) Then
       ! used those read from file
       Umean_resc_T = Umean_resc_To
       Vmean_resc_T = Vmean_resc_To
    Else
       ! initialize for first time
       Do j=1,nyg
          Umean_resc_T_local (j,1) = Sum(U(i_rescale,j,2:nzg-1))
          Umean_inlet_T_local(j,1) = Sum(U(        1,j,2:nzg-1))
       End Do
       Do j=1,ny
          Vmean_resc_T_local (j,1) = Sum(V(i_rescale,j,2:nzg-1))
          Vmean_inlet_T_local(j,1) = Sum(V(        1,j,2:nzg-1))
       End Do
       
       Umean_resc_T(:,1) = 0d0
       Vmean_resc_T(:,1) = 0d0
       Call MPI_Allreduce(Umean_resc_T_local, Umean_resc_T, nyg,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       Call MPI_Allreduce(Vmean_resc_T_local, Vmean_resc_T, ny ,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
       Umean_resc_T = Umean_resc_T/Real(nzg_global-2,8)
       Vmean_resc_T = Vmean_resc_T/Real(nzg_global-2,8)
    End If

    Umean_inlet_T = 0d0
    Vmean_inlet_T = 0d0
    Call MPI_Allreduce(Umean_inlet_T_local, Umean_inlet_T, nyg,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_Allreduce(Vmean_inlet_T_local, Vmean_inlet_T, ny ,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Umean_inlet_T = Umean_inlet_T/Real(nzg_global-2,8)
    Vmean_inlet_T = Vmean_inlet_T/Real(nzg_global-2,8)

    !---------------statistics for z modes---------------------!
    ! delta_inlet/x0 reference = 0.0115
    j_stat = 0
    Do j = 1, nyg_global
       !If ( yg_global(j)>3d0*0.0115 ) Then
       If ( yg_global(j)>0.1d0 ) Then
          j_stat = j
          Exit
       End If
    End Do
    If ( j_stat==0 ) Then 
       j_stat = ny
       If (myid==0) write(*,*) 'Warning! j_stat not found'
    ENd If

    ! last station 3.0*x0
    i_stat = 0
    Do i = 1, nx_global
       !If ( x_global(i)>3.0d0 ) Then
       If ( x_global(i)>4.0d0 ) Then
          i_stat = i
          Exit
       End If
    End Do
    If ( i_stat==0 ) Then 
       i_stat = nx
       If (myid==0) write(*,*) 'Warning! i_stat not found'
    End If
    ! save 100 point in x
    !Delta_i_stat = Nint( Real(i_stat)/Real(100) )
    !i_stat       = Delta_i_stat*100
    !
    !Delta_i_stat = Nint( Real(i_stat)/Real(1000) )
    !i_stat       = Delta_i_stat*1000
    !
    Delta_i_stat = Nint( Real(i_stat)/Real(50) )
    i_stat       = Delta_i_stat*50
    !
    i_stat       = Min(i_stat,nx_global)
    If ( Delta_i_stat==0 ) Stop 'Error!: Delta_i_stat=0'

    nxu_reduced = 0
    Do i = 1, i_stat, Delta_i_stat
       nxu_reduced = nxu_reduced + 1
    End Do
    nyu_reduced = 0
    Do j = 1, j_stat
       nyu_reduced = nyu_reduced + 1
    End Do
    nzu_reduced     = (nzg_global-2) - 2 + 1
    nzu_modes       = nzu_reduced/2+1
    nzu_first_modes = Min(6,nzg_global-4)
    If ( Mod(nzu_reduced,2)/=0 .And. nprocs > 1 ) Stop 'Error! nzu_reduced must be even'
    If ( Mod(nzu_reduced,2)/=0 ) nzu_reduced = nzu_reduced - 1  ! trim for odd case
    Call initialize_fftz(nzu_reduced)

    Allocate( xu_reduced(nxu_reduced) )
    Allocate( yu_reduced(nyu_reduced) )
    Allocate( zu_reduced(nzu_reduced) )

    xu_reduced = x_global(1:i_stat:Delta_i_stat)
    yu_reduced = yg_global(1:j_stat)
    zu_reduced = zg_global(2:nzg_global-2)

    Allocate( U_reduced      (nxu_reduced, nyu_reduced, nzu_reduced) )
    Allocate( U_reduced_hat_z(nxu_reduced, nyu_reduced, 0:nzu_first_modes-1) )

    Allocate( rhs_u_reduced      (nxu_reduced, nyu_reduced, nzu_reduced) )
    Allocate( rhs_u_reduced_hat_z(nxu_reduced, nyu_reduced, 0:nzu_first_modes-1) )

    Allocate( U_plane      (nxu_reduced, nzu_reduced) )
    Allocate( U_plane_hat_z(nxu_reduced, 0:nzu_modes-1 ) )

    !----------------------wall stress model-------------------!
    Allocate(utau_model (nx_global))
    Allocate(utau_wall  (nx_global))
    Allocate(utau_wall_T(nx_global))
    Allocate(utau_ref   (nx_global))
    ! utau reference for turbulent flows
    utau_ref    = ( 0.5d0*0.027d0*(x_global/nu)**(-1d0/7d0) )**0.5d0 
    utau_model  = utau_ref ! initial guess
    utau_wall   = 0d0      ! actual utau at the wall 
    utau_wall_T = 0d0      ! averaged actual utau at the wall 
    beta_y      = 0d0
    
    !-------------------------Done-----------------------------!
    Call Mpi_barrier(MPI_COMM_WORLD,ierr)

    ! Measure time
    time1 = MPI_WTIME()

    If ( myid==0 ) Write(*,*) 'initialization done...' 
    
  End Subroutine initialize
  
End Module initialization
