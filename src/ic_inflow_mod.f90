!------------------------------------------------------------------------------
! ic_inflow_mod — Blasius profiles, initial condition, inflow/top BC tables
!
! Part of the CUDA Fortran port of the boundary-layer DNS solver
! (see docs/CUDA_PORT_DESIGN.md). Single GPU, nprocs = 1, no MPI.
!
! Reference provenance (verbatim numerics):
!   * input_output.f90, subroutine init_flow (lines ~176-277):
!       Blasius initial condition for random_init == 1 (IC part only; the
!       face-grid generation of init_flow lives in grid_mod:grid_generate).
!       The xg_init/yg_init center+ghost grids that init_flow derives locally
!       are recomputed here with the exact reference expressions (bitwise
!       identical to grid_mod's xg/yg, which use the same formulas).
!   * boundary_conditions.f90, subroutine compute_blasius_solution_for_bc
!     (lines ~1452-1797):
!       PART 1 — inlet profiles U_inlet/V_inlet/W_inlet and top profiles
!                U_top/V_top/W_top from the self-similar Blasius solution.
!       PART 2 — temporal inflow-mode tables read from file_temporal_inlet
!                (qu/qv/qw_inlet complex modes, zmode/tmode_inlet wavenumbers)
!                interpolated to the present y-grids, or the no-modes path
!                when the input holds the '**TO DO**' placeholder.
!
! Blasius normalization (U_inf = 1, x0 = distance to leading edge):
!       U/U_inf = f'(eta)
!       V/U_inf = 1/sqrt(2*Rex) * ( eta*f'(eta) - f(eta) )
!       W/U_inf = 0
!       eta     = y * sqrt( U_inf/(2*nu*x) ),   Rex = U_inf*x/nu
! The file 'file_inflow' (formatted) holds: n_source, then eta(:), f(:),
! f'(:) each as one list-directed record.
!
! Temporal inflow modes (used by bc_kernels on the device every substep):
!       Ut_inlet(j,k) = sum_{n,m} Re[ qu_inlet(j,n,m)
!                        * exp( i*zmode_inlet(n)*zg(k) - i*tmode_inlet(m)*t ) ]
! (and likewise Vt with y-faces, Wt with z-faces). zmode_inlet holds
! beta_inlet*n, tmode_inlet holds omega_inlet*m, both premultiplied here.
!
! Ownership contract:
!   * This module OWNS the inflow/top tables: host arrays with the exact
!     reference names (global.f90) plus device copies suffixed `_d`, filled
!     by inflow_tables_to_device().
!   * Ut/Vt/Wt_inlet are the per-substep temporal-perturbation planes; they
!     are recomputed ON THE DEVICE by bc_kernels when n_modes_inlet > 0 and
!     must remain zero when n_modes_inlet == 0. The reference relies on
!     `!$acc create` leaving them zero in that case (uninitialized device
!     memory that happens to be zero); the port zero-initializes explicitly.
!   * Restart reading (random_init == 0, input_output.f90:read_input_data)
!     is NOT implemented here: it is snapshot-format I/O and belongs to
!     io_mod (assumed public routine `read_restart()`, see report). The
!     driver selects generate_initial_condition() vs. read_restart().
!
! Call order (driver): grid_generate -> fields_allocate ->
!   generate_initial_condition (or io_mod read_restart) ->
!   compute_blasius_solution_for_bc -> inflow_tables_to_device.
!------------------------------------------------------------------------------
module ic_inflow_mod

  use precision_mod, only: dp, pi
  use param_mod,     only: nu, CFL, file_inflow, file_temporal_inlet, &
                           inflow_boundary_flag, &
                           Lx_rand, Ly_rand, Lz_rand, alpha_rand
  use grid_mod,      only: nx, ny, nz, nxm, nym, nxg, nyg, nzg, x, y, xg, yg
  use field_mod,     only: U, V, W
  use cudafor

  implicit none
  private

  !----------------------------- public API ----------------------------------
  public :: generate_initial_condition    ! Blasius IC (random_init == 1)
  public :: compute_blasius_solution_for_bc ! inlet/top profiles + mode tables
  public :: inflow_tables_to_device       ! host tables -> device copies

  !------------------- inlet / top boundary profiles (host) ------------------
  ! (reference names, global.f90; allocated by compute_blasius_solution_for_bc)
  real(dp), allocatable, public :: U_inlet(:)  ! (nyg) u at inlet, y-centers
  real(dp), allocatable, public :: V_inlet(:)  ! (ny ) v at inlet, y-faces
  real(dp), allocatable, public :: W_inlet(:)  ! (nyg) w at inlet, y-centers
  real(dp), allocatable, public :: U_top(:)    ! (nx ) u at top, x-faces
  real(dp), allocatable, public :: V_top(:)    ! (nxg) v at top, x-centers
  real(dp), allocatable, public :: W_top(:)    ! (nxg) w at top, x-centers

  !------------- temporal inflow-perturbation planes (host mirrors) ----------
  ! Zero when n_modes_inlet == 0; otherwise recomputed on the device by
  ! bc_kernels every substep (host copies stay zero after setup).
  real(dp), allocatable, public :: Ut_inlet(:,:) ! (nyg,nzg)
  real(dp), allocatable, public :: Vt_inlet(:,:) ! (ny ,nzg)
  real(dp), allocatable, public :: Wt_inlet(:,:) ! (nyg,nz )

  !--------------------- temporal inflow-mode tables (host) ------------------
  integer,  public :: ny_inlet      = 0      ! wall-normal points in mode file
  integer,  public :: n_modes_inlet = 0      ! number of spanwise (z) modes
  integer,  public :: m_modes_inlet = 0      ! number of temporal modes
  real(dp), public :: beta_inlet    = 0.0_dp ! fundamental spanwise wavenumber
  real(dp), public :: omega_inlet   = 1.0_dp ! fundamental angular frequency
  real(dp), public :: dt_period     = 1.0_dp ! dt as a fraction of the period

  real(dp),    allocatable, public :: ymesh_inlet(:) ! (ny_inlet) mode-file y
  real(dp),    allocatable, public :: zmode_inlet(:) ! (n_modes) beta_inlet*n
  real(dp),    allocatable, public :: tmode_inlet(:) ! (m_modes) omega_inlet*m
  complex(dp), allocatable, public :: qu_inlet(:,:,:) ! (nyg,n_modes,m_modes)
  complex(dp), allocatable, public :: qv_inlet(:,:,:) ! (ny ,n_modes,m_modes)
  complex(dp), allocatable, public :: qw_inlet(:,:,:) ! (nyg,n_modes,m_modes)

  !--------------------------- device copies ---------------------------------
  ! Same names + `_d` (design-doc convention); filled by inflow_tables_to_device.
  real(dp), device, allocatable, public :: U_inlet_d(:), V_inlet_d(:), W_inlet_d(:)
  real(dp), device, allocatable, public :: U_top_d(:),   V_top_d(:),   W_top_d(:)
  real(dp), device, allocatable, public :: Ut_inlet_d(:,:), Vt_inlet_d(:,:), Wt_inlet_d(:,:)
  real(dp), device, allocatable, public :: zmode_inlet_d(:), tmode_inlet_d(:)
  complex(dp), device, allocatable, public :: qu_inlet_d(:,:,:)
  complex(dp), device, allocatable, public :: qv_inlet_d(:,:,:)
  complex(dp), device, allocatable, public :: qw_inlet_d(:,:,:)

contains

  !----------------------------------------------------------------------------
  ! generate_initial_condition — Blasius initial condition (random_init == 1).
  !
  ! Fills the field_mod host mirrors U, V, W with the self-similar Blasius
  ! solution evaluated on the staggered grid and sets t = 0:
  !   U(x-faces, y-centers) = U_inf * f'(eta),
  !   V(x-centers, y-faces) = U_inf/sqrt(2*Rex) * (eta*f' - f),
  !   W = 0 (2D laminar flow),
  ! with linear interpolation of (f, f') in eta from the file table.
  !
  ! Verbatim transcription of input_output.f90:init_flow (IC part). The face
  ! grids are generated beforehand by grid_mod:grid_generate with the same
  ! reference numerics; the xg_init/yg_init center+ghost grids are recomputed
  ! here exactly as init_flow derives them (bitwise identical to xg/yg).
  ! NOTE: this IC path contains NO random_number calls in the reference; the
  ! generated field is z-uniform, so the reference's post-IC periodic-z fixup
  ! (initialization.f90:244-246) is a no-op and is not repeated here.
  ! The caller then pushes U,V,W to the device (fields_to_device).
  !----------------------------------------------------------------------------
  subroutine generate_initial_condition(t)

    real(dp), intent(out) :: t   !< simulation start time (set to 0)

    integer :: ii, jj, k
    integer :: n_source, j0, j1
    real(dp), allocatable :: eta_source(:), f_source(:), df_source(:)
    real(dp) :: eta_local, w0, w1, U_inf, Rex_ref
    real(dp), allocatable :: xg_init(:), yg_init(:)

    write(*,*) 'Generating initial condition'

    ! start time
    t = 0.0_dp

    !---- center grids (yg_init, xg_init) needed for the IC ----
    ! Recomputed locally exactly as init_flow does (input_output.f90:179-191);
    ! same expressions and operation order as grid_mod's xg/yg.
    ! yg_init: center points + ghost for U,W (y-centers)
    allocate( yg_init(nyg) )
    do jj = 1, nym
       yg_init(jj+1) = 0.5_dp*( y(jj) + y(jj+1) )
    end do
    yg_init(1)   = yg_init(2)     - 2.0_dp*( yg_init(2) - y(1) )
    yg_init(nyg) = yg_init(nym+1) + 2.0_dp*( y(ny) - yg_init(nym+1) )

    ! xg_init: center points + ghost for V (x-centers)
    allocate( xg_init(nxg) )
    do ii = 1, nxm
       xg_init(ii+1) = 0.5_dp*( x(ii) + x(ii+1) )
    end do
    xg_init(1)   = xg_init(2)     - 2.0_dp*( xg_init(2) - x(1) )
    xg_init(nxg) = xg_init(nxm+1) + 2.0_dp*( x(nx) - xg_init(nxm+1) )

    !---- read Blasius self-similar solution ----
    U_inf = 1.0_dp
    call read_blasius_source(n_source, eta_source, f_source, df_source)
    write(*,*) '   Blasius IC: read ', n_source, ' points from ', trim(file_inflow)

    !---- initialize U with Blasius: U(x,y) = U_inf * f'(eta) ----
    ! U is at x-faces (x), y-centers (yg_init)
    U = U_inf
    do ii = 1, nx
       do jj = 1, nyg
          eta_local = yg_init(jj) * sqrt( U_inf/(2.0_dp*nu*x(ii)) )
          j0 = 0
          do k = 2, n_source
             if ( eta_source(k) > eta_local ) then
                j0 = k - 1
                j1 = k
                w1 = (eta_local - eta_source(j0))/(eta_source(j1) - eta_source(j0))
                w0 = 1.0_dp - w1
                exit
             end if
          end do
          if ( j0 > 0 ) then
             U(ii,jj,:) = U_inf * (w0*df_source(j0) + w1*df_source(j1))
          else
             U(ii,jj,:) = U_inf
          end if
       end do
    end do

    !---- initialize V with Blasius: V(x,y) = U_inf/sqrt(2*Rex)*(eta*f'-f) ----
    ! V is at x-centers (xg_init), y-faces (y)
    V = 0.0_dp
    do ii = 1, nxg
       Rex_ref = U_inf*xg_init(ii)/nu
       do jj = 1, ny
          eta_local = y(jj) * sqrt( U_inf/(2.0_dp*nu*xg_init(ii)) )
          j0 = 0
          do k = 2, n_source
             if ( eta_source(k) > eta_local ) then
                j0 = k - 1
                j1 = k
                w1 = (eta_local - eta_source(j0))/(eta_source(j1) - eta_source(j0))
                w0 = 1.0_dp - w1
                exit
             end if
          end do
          if ( j0 > 0 ) then
             V(ii,jj,:) = U_inf/sqrt(2.0_dp*Rex_ref) * &
                ( w0*(eta_source(j0)*df_source(j0)-f_source(j0)) + &
                  w1*(eta_source(j1)*df_source(j1)-f_source(j1)) )
          else
             V(ii,jj,:) = U_inf/sqrt(2.0_dp*Rex_ref) * &
                (eta_source(n_source)*df_source(n_source)-f_source(n_source))
          end if
       end do
    end do

    ! W = 0 (2D laminar flow)
    W = 0.0_dp

    deallocate(eta_source, f_source, df_source, yg_init, xg_init)

    write(*,*) 'Blasius initial condition:'
    write(*,*) '   Lx     : ', Lx_rand
    write(*,*) '   Ly     : ', Ly_rand
    write(*,*) '   Lz     : ', Lz_rand
    write(*,*) '   alpha  : ', alpha_rand

    write(*,*) '   Max U : ', maxval(U)
    write(*,*) '   Max V : ', maxval(V)
    write(*,*) '   Max W : ', maxval(W)

    write(*,*) '   Mean U : ', sum(U)/real(nx *nyg*nzg, dp)
    write(*,*) '   Mean V : ', sum(V)/real(nxg*ny *nzg, dp)
    write(*,*) '   Mean W : ', sum(W)/real(nxg*nyg*nz , dp)

  end subroutine generate_initial_condition

  !----------------------------------------------------------------------------
  ! compute_blasius_solution_for_bc — inlet/top profiles + temporal mode tables.
  !
  ! Verbatim transcription of boundary_conditions.f90:
  ! compute_blasius_solution_for_bc (single rank; MPI broadcasts dropped).
  !
  ! PART 1 — Blasius profiles for the boundary conditions:
  !   U_inlet(nyg), V_inlet(ny), W_inlet(nyg): self-similar solution at the
  !     inlet station x(1)/xg(1) (linear interpolation in eta of the file
  !     table; U -> 1 and V -> its edge value above the tabulated range).
  !   U_top(nx) = 1; V_top(nxg) = Blasius edge value of V at each x-station
  !     (eta -> infinity limit, last table point); W_top(nxg) = 0.
  !   The reference's U_outlet/V_outlet/W_outlet are NOT ported: they are
  !   never used by the GPU BC path (convective outflow uses Uo), dead code.
  !
  ! PART 2 — temporal inflow modes:
  !   If file_temporal_inlet is empty or holds the literal '**TO DO**'
  !   placeholder (tested exactly as the reference: first four characters
  !   '**TO'; param_mod's parser preserves this), there are no perturbations:
  !   Ut/Vt/Wt_inlet = 0, n_modes_inlet = m_modes_inlet = 0, beta_inlet = 0,
  !   omega_inlet = 1, dt_period = 1, and zero-size mode tables are allocated
  !   so the device mirroring is well defined.
  !   Otherwise the stream-unformatted mode file is read:
  !     ny_inlet, n_modes_inlet, m_modes_inlet (int4), beta_inlet,
  !     omega_inlet (real8), then ymesh_inlet, zmode_inlet, tmode_inlet and
  !     the real/imaginary parts of the (qu,qv,qw) mode shapes on ymesh_inlet.
  !   Wavenumbers are premultiplied (zmode = beta*n, tmode = omega*m), the
  !   complex modes are linearly interpolated in y to the present yg (u,w
  !   at centers) and y (v at faces) grids — zero above the file's y-range —
  !   and dt_period = 2*pi/(|CFL|*omega_inlet)/1000 (fixed-dt convention).
  !
  ! Requires grid_generate() (x, y, xg, yg) and read_input_parameters().
  !----------------------------------------------------------------------------
  subroutine compute_blasius_solution_for_bc()

    ! solution from file (different size than the mesh)
    real(dp), allocatable :: eta_source(:), f_source(:), df_source(:)

    ! temporal-mode work arrays (reference qu_inlet_r/_i/_o etc., local here)
    real(dp),    allocatable :: qu_inlet_r(:,:,:), qv_inlet_r(:,:,:), qw_inlet_r(:,:,:)
    real(dp),    allocatable :: qu_inlet_i(:,:,:), qv_inlet_i(:,:,:), qw_inlet_i(:,:,:)
    complex(dp), allocatable :: qu_inlet_o(:,:,:), qv_inlet_o(:,:,:), qw_inlet_o(:,:,:)

    real(dp) :: eta_local, Rex0, Rex_ref, w0, w1, U_inf
    integer  :: j, jj, j0, j1, i_ref, n_source
    integer  :: funit, ios

    !----------------------------------------------------------------------!
    ! PART 1: compute Blasius

    ! set U_inf
    U_inf = 1.0_dp

    ! allocate boundary velocities (nprocs = 1: nyg_global = nyg, ...)
    allocate( U_inlet(nyg), V_inlet(ny), W_inlet(nyg) )
    allocate( U_top(nx), V_top(nxg), W_top(nxg) )

    ! read self-similar Blasius solution
    call read_blasius_source(n_source, eta_source, f_source, df_source)

    ! inflow Rex (informational in the reference; kept)
    Rex0 = 1.0_dp*x(1)/nu

    if ( inflow_boundary_flag==1 .or. inflow_boundary_flag==3 .or. &
         inflow_boundary_flag==4 ) then

       ! Generate own Blasius
       write(*,*) 'Generating own Blasius for inlet'

       ! compute solution at inlet
       ! U
       U_inlet = 1.0_dp ! old version: U_inlet = 0d0
       i_ref   = 1
       Rex_ref = 1.0_dp*x(i_ref)/nu
       do j = 1, nyg
          eta_local = yg(j)*(U_inf/(2.0_dp*nu*x(i_ref)))**0.5_dp
          j0 = 0
          do jj = 2, n_source
             if ( eta_source(jj) > eta_local ) then
                j0 = jj - 1
                j1 = jj
                w1 = ( eta_local - eta_source(j0) )/( eta_source(j1) - eta_source(j0) )
                w0 = 1.0_dp - w1
                exit
             end if
          end do
          if ( j0 > 0 ) then
             U_inlet(j) = w0*df_source(j0) + w1*df_source(j1)
          else
             U_inlet(j) = 1.0_dp
          end if
       end do
       ! V
       V_inlet = 0.0_dp
       i_ref   = 1
       Rex_ref = 1.0_dp*xg(i_ref)/nu
       do j = 1, ny
          eta_local = y(j)*(U_inf/(2.0_dp*nu*xg(i_ref)))**0.5_dp
          j0 = 0
          do jj = 2, n_source
             if ( eta_source(jj) > eta_local ) then
                j0 = jj - 1
                j1 = jj
                w1 = ( eta_local - eta_source(j0) )/( eta_source(j1) - eta_source(j0) )
                w0 = 1.0_dp - w1
                exit
             end if
          end do
          if ( j0 > 0 ) then
             V_inlet(j) = w0*1.0_dp/(2.0_dp*Rex_ref)**0.5_dp*(eta_source(j0)*df_source(j0)-f_source(j0)) + &
                  w1*1.0_dp/(2.0_dp*Rex_ref)**0.5_dp*(eta_source(j1)*df_source(j1)-f_source(j1))
          else
             V_inlet(j) = maxval(V_inlet)
          end if
       end do
       ! W
       W_inlet = 0.0_dp

    else if ( inflow_boundary_flag==2 ) then

       ! reference: read Blasius profile from a binary file. Not ported
       ! (check_supported rejects inflow_flag /= 1 before reaching here).
       stop 'inflow_flag = 2 (Blasius profile from file) is not ported'

    else
       stop 'Error! inflow_boundary_flag unknown'
    end if

    ! compute solution at top
    ! U
    U_top = 1.0_dp
    ! V
    do i_ref = 1, nxg
       Rex_ref      = 1.0_dp*xg(i_ref)/nu
       j            = ny
       eta_local    = y(j)*(U_inf/(2.0_dp*nu*xg(i_ref)))**0.5_dp
       j0           = n_source
       V_top(i_ref) = 1.0_dp/(2.0_dp*Rex_ref)**0.5_dp*(eta_source(j0)*df_source(j0)-f_source(j0))
    end do
    ! W
    W_top = 0.0_dp

    deallocate(eta_source, f_source, df_source)

    !----------------------------------------------------------------------!
    ! PART 2: prepare temporal component for boundary conditions

    ! allocate inlet temporal component
    allocate( Ut_inlet(nyg,nzg) )
    allocate( Vt_inlet(ny ,nzg) )
    allocate( Wt_inlet(nyg,nz ) )

    ! Skip temporal modes if no file is provided (pure Blasius, no
    ! perturbations). The '**TO' test is the reference's placeholder check
    ! (boundary_conditions.f90:1653) and relies on param_mod's parser not
    ! left-adjusting the value after '=' — do not "fix" either side.
    if ( len_trim(file_temporal_inlet)==0 .or. &
         file_temporal_inlet(1:4)=='**TO' ) then
       write(*,*) 'No temporal inflow file -> pure Blasius (no perturbations)'
       Ut_inlet = 0.0_dp
       Vt_inlet = 0.0_dp
       Wt_inlet = 0.0_dp
       n_modes_inlet = 0
       m_modes_inlet = 0
       beta_inlet    = 0.0_dp
       omega_inlet   = 1.0_dp
       dt_period     = 1.0_dp
       ! zero-size tables so inflow_tables_to_device stays well defined
       allocate( zmode_inlet(0), tmode_inlet(0) )
       allocate( qu_inlet(nyg,0,0), qv_inlet(ny,0,0), qw_inlet(nyg,0,0) )
       return
    end if

    ! read sizes and wavenumbers (stream binary: 3 x int4, 2 x real8)
    open(newunit=funit, file=file_temporal_inlet, form='unformatted', &
         action='read', access='stream', iostat=ios)
    if ( ios /= 0 ) then
       write(*,*) 'ERROR: cannot open timeinflow_file: ', trim(file_temporal_inlet)
       stop 'missing temporal inflow file'
    end if
    read(funit)      ny_inlet
    read(funit) n_modes_inlet
    read(funit) m_modes_inlet
    read(funit)    beta_inlet
    read(funit)   omega_inlet

    ! time multiple of the period (fixed-dt convention: dt = |CFL|)
    dt_period = 2.0_dp*pi/(abs(CFL)*omega_inlet)/real(1000,dp)

    ! allocate y-mesh, wavenumbers and coefficients
    allocate( ymesh_inlet(     ny_inlet) ) ! wall-normal points
    allocate( zmode_inlet(n_modes_inlet) ) ! number of spanwise modes
    allocate( tmode_inlet(m_modes_inlet) ) ! number of temporal modes

    allocate( qu_inlet(nyg,n_modes_inlet,m_modes_inlet) ) ! u modes y-interpolated
    allocate( qv_inlet(ny ,n_modes_inlet,m_modes_inlet) ) ! v modes y-interpolated
    allocate( qw_inlet(nyg,n_modes_inlet,m_modes_inlet) ) ! w modes y-interpolated

    ! temporal arrays
    ! real part
    allocate( qu_inlet_r(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! u modes from source
    allocate( qv_inlet_r(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! v modes from source
    allocate( qw_inlet_r(ny_inlet,n_modes_inlet,m_modes_inlet) ) ! w modes from source
    ! imaginary part
    allocate( qu_inlet_i(ny_inlet,n_modes_inlet,m_modes_inlet) )
    allocate( qv_inlet_i(ny_inlet,n_modes_inlet,m_modes_inlet) )
    allocate( qw_inlet_i(ny_inlet,n_modes_inlet,m_modes_inlet) )
    ! complex
    allocate( qu_inlet_o(ny_inlet,n_modes_inlet,m_modes_inlet) )
    allocate( qv_inlet_o(ny_inlet,n_modes_inlet,m_modes_inlet) )
    allocate( qw_inlet_o(ny_inlet,n_modes_inlet,m_modes_inlet) )

    ! read
    read(funit) ymesh_inlet
    read(funit) zmode_inlet
    read(funit) tmode_inlet
    read(funit) qu_inlet_r
    read(funit) qu_inlet_i
    read(funit) qv_inlet_r
    read(funit) qv_inlet_i
    read(funit) qw_inlet_r
    read(funit) qw_inlet_i
    close(funit)

    ! build wavenumber
    zmode_inlet =  beta_inlet*zmode_inlet
    tmode_inlet = omega_inlet*tmode_inlet

    ! complex modes
    qu_inlet_o = cmplx( qu_inlet_r, qu_inlet_i, kind=dp )
    qv_inlet_o = cmplx( qv_inlet_r, qv_inlet_i, kind=dp )
    qw_inlet_o = cmplx( qw_inlet_r, qw_inlet_i, kind=dp )

    ! interpolate qu and qw to yg (u,w at y-centers)
    qu_inlet = (0.0_dp, 0.0_dp)
    qw_inlet = (0.0_dp, 0.0_dp)
    do j = 1, nyg
       j1 = 0
       do jj = 2, ny_inlet
          if ( ymesh_inlet(jj) >= yg(j) ) then
             j1 = jj
             exit
          end if
          if ( jj == ny_inlet ) j1 = ny_inlet
       end do
       if ( j1 == ny_inlet ) then
          qu_inlet(j,:,:) = 0.0_dp
          qw_inlet(j,:,:) = 0.0_dp
       else
          j0 = j1 - 1
          qu_inlet(j,:,:) = qu_inlet_o(j0,:,:) + &
                            (yg(j)-ymesh_inlet(j0))*(qu_inlet_o(j1,:,:)-qu_inlet_o(j0,:,:))/(ymesh_inlet(j1)-ymesh_inlet(j0))
          qw_inlet(j,:,:) = qw_inlet_o(j0,:,:) + &
                            (yg(j)-ymesh_inlet(j0))*(qw_inlet_o(j1,:,:)-qw_inlet_o(j0,:,:))/(ymesh_inlet(j1)-ymesh_inlet(j0))
       end if
    end do

    ! interpolate qv to present y mesh (v at y-faces)
    qv_inlet = (0.0_dp, 0.0_dp)
    do j = 1, ny
       j1 = 0
       do jj = 2, ny_inlet
          if ( ymesh_inlet(jj) >= y(j) ) then
             j1 = jj
             exit
          end if
          if ( jj == ny_inlet ) j1 = ny_inlet
       end do
       if ( j1 == ny_inlet ) then
          qv_inlet(j,:,:) = 0.0_dp
       else
          j0 = j1 - 1
          qv_inlet(j,:,:) = qv_inlet_o(j0,:,:) + &
                            (y(j)-ymesh_inlet(j0))*(qv_inlet_o(j1,:,:)-qv_inlet_o(j0,:,:))/(ymesh_inlet(j1)-ymesh_inlet(j0))
       end if
    end do

    ! deallocate
    deallocate(qu_inlet_r,qv_inlet_r,qw_inlet_r)
    deallocate(qu_inlet_i,qv_inlet_i,qw_inlet_i)
    deallocate(qu_inlet_o,qv_inlet_o,qw_inlet_o)

    ! perturbation planes start at zero; bc_kernels rebuilds them on the
    ! device every substep from the mode tables and the current time
    Ut_inlet = 0.0_dp
    Vt_inlet = 0.0_dp
    Wt_inlet = 0.0_dp

  end subroutine compute_blasius_solution_for_bc

  !----------------------------------------------------------------------------
  ! inflow_tables_to_device — allocate the device copies (first call) and copy
  ! all inflow/top tables host -> device. Call once after
  ! compute_blasius_solution_for_bc(); the tables are time-invariant (the
  ! time-dependent Ut/Vt/Wt_inlet_d planes are OWNED here but recomputed in
  ! place on the device by bc_kernels when n_modes_inlet > 0 — when there are
  ! no modes they keep the zeros copied here). Assignments between host and
  ! device allocatables are synchronous copies (CUDA Fortran semantics).
  !----------------------------------------------------------------------------
  subroutine inflow_tables_to_device()

    if (.not. allocated(U_inlet_d)) then
       allocate( U_inlet_d(nyg), V_inlet_d(ny), W_inlet_d(nyg) )
       allocate( U_top_d(nx), V_top_d(nxg), W_top_d(nxg) )
       allocate( Ut_inlet_d(nyg,nzg), Vt_inlet_d(ny,nzg), Wt_inlet_d(nyg,nz) )
       ! The CUDA runtime rejects zero-byte device allocations ("0 bytes
       ! requested"), so the mode tables exist on the device only when there
       ! are modes; bc_kernels touches them only when n_modes_inlet > 0.
       if ( n_modes_inlet > 0 ) then
          allocate( zmode_inlet_d(n_modes_inlet), tmode_inlet_d(m_modes_inlet) )
          allocate( qu_inlet_d(nyg,n_modes_inlet,m_modes_inlet) )
          allocate( qv_inlet_d(ny ,n_modes_inlet,m_modes_inlet) )
          allocate( qw_inlet_d(nyg,n_modes_inlet,m_modes_inlet) )
       end if
    end if

    U_inlet_d = U_inlet
    V_inlet_d = V_inlet
    W_inlet_d = W_inlet
    U_top_d   = U_top
    V_top_d   = V_top
    W_top_d   = W_top

    ! zero (no modes) or zeroed scratch (recomputed on device each substep)
    Ut_inlet_d = Ut_inlet
    Vt_inlet_d = Vt_inlet
    Wt_inlet_d = Wt_inlet

    if ( n_modes_inlet > 0 ) then
       zmode_inlet_d = zmode_inlet
       tmode_inlet_d = tmode_inlet
       qu_inlet_d    = qu_inlet
       qv_inlet_d    = qv_inlet
       qw_inlet_d    = qw_inlet
    end if

  end subroutine inflow_tables_to_device

  !----------------------------------------------------------------------------
  ! read_blasius_source — read the self-similar Blasius table (file_inflow).
  !
  ! Formatted file: n_source, then eta(:), f(:), f'(:) as list-directed
  ! records. Shared by the IC and the BC-profile builder, which in the
  ! reference read the same file with identical inline code
  ! (input_output.f90:196-202, boundary_conditions.f90:1475-1483).
  !----------------------------------------------------------------------------
  subroutine read_blasius_source(n_source, eta_source, f_source, df_source)

    integer,               intent(out) :: n_source
    real(dp), allocatable, intent(out) :: eta_source(:), f_source(:), df_source(:)

    integer :: funit, ios

    open(newunit=funit, file=file_inflow, form='formatted', action='read', &
         iostat=ios)
    if ( ios /= 0 ) then
       write(*,*) 'ERROR: cannot open inflow_file: ', trim(file_inflow)
       stop 'missing Blasius profile file'
    end if
    read(funit,*) n_source
    allocate( eta_source(n_source), f_source(n_source), df_source(n_source) )
    read(funit,*) eta_source
    read(funit,*)   f_source
    read(funit,*)  df_source
    close(funit)

  end subroutine read_blasius_source

end module ic_inflow_mod
