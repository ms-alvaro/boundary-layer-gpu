!------------------------------------------------------------------------------
! param_mod — input-file parsing and runtime parameters
!
! Part of the CUDA Fortran port of the boundary-layer DNS solver
! (see docs/CUDA_PORT_DESIGN.md). Single GPU, nprocs = 1, no MPI.
!
! Provenance (reference files, transcribed faithfully):
!   - params.f90            subroutine initparams: keyword list, parse order,
!                           '-i <file>' command-line handling, output-directory
!                           creation, error messages.
!   - miscel.f90            getvar / get_dbl / get_dbl_arr / get_int /
!                           get_int_arr / get_str: the '.turbb' key = value
!                           parser. Transcribed exactly, including the fact
!                           that the value string keeps any characters right
!                           after '=' unshifted (no adjustl). This is load-
!                           bearing: downstream code recognizes the literal
!                           placeholder '**TO DO**' for timeinflow_file by
!                           testing file_temporal_inlet(1:4) == '**TO', which
!                           only works because the harness input files write
!                           `timeinflow_file =**TO DO**` with no space and
!                           getvar does not left-adjust the value.
!   - input_output.f90      subroutine read_input_parameters: defaults set
!                           before parsing (nstep_init, nstep_init_input,
!                           Rossby_plus, Amplitude_perturbations) and the
!                           Omega_z computation. MPI broadcasts dropped.
!   - initialization.f90    derived grid-point counts (nxm/nxg families) for
!                           the single-rank case, and the startup banner.
!   - monitor.f90           subroutine summary: parameter echo (the grid-
!                           coordinate lines that need x/y/z arrays are left
!                           to the module that owns the grids).
!
! Time-step convention (compute_dt, time_integration.f90): CFL > 0 means a
! true CFL condition; CFL < 0 means a fixed time step dt = |CFL| (the CFL
! computation is skipped entirely). param_mod only stores CFL; the convention
! is applied in timestep_mod.
!
! Input format ('.turbb'): free-form lines `keyword = value [! comment]`.
! Unknown lines are ignored; a keyword that is absent leaves the target
! variable at its default. Strings must not be quoted.
!------------------------------------------------------------------------------
module param_mod

  use precision_mod, only: dp
  use mpi_mod,       only: nprocs, is_root

  implicit none
  private

  !----------------------------- public API ---------------------------------
  public :: read_input_parameters   ! parse '-i <file>', fill all parameters
  public :: check_supported         ! stop if an unported feature is requested
  public :: print_banner            ! startup banner (initialization.f90)
  public :: param_summary           ! parameter echo   (monitor.f90 summary)

  ! parameters file name (from '-i <file>')
  public :: paramsfilename

  ! grid-point counts (single rank: local = global)
  public :: nx_global, ny_global, nz_global
  public :: nxm_global, nym_global, nzm_global
  public :: nxg_global, nyg_global, nzg_global
  ! NOTE: the derived local counts nx..nzg are intentionally NOT public —
  ! grid_mod owns them (set by grid_generate); everything downstream
  ! use-associates the grid_mod versions. The private copies below are
  ! kept only for param-internal consistency checks.

  ! physics / numerics
  public :: nu, CFL, dPdx, dPdz, Omega_z, itime_step

  ! run control
  public :: nsteps, nsave, nstats, nmonitor
  public :: nstep_init, nstep_init_input, random_init

  ! files
  public :: filein, fileout, file_inflow, file_temporal_inlet
  public :: file_hit_planes, N_buffer_hit, file_ygrid
  public :: boxout_every, boxout_start, n_boxout, dir_boxout
  public :: boxout_i0, boxout_i1, boxout_is, boxout_jmax

  ! boundary conditions
  public :: inflow_boundary_flag, top_boundary_flag
  public :: Amplitude_perturbations
  public :: Vbs_max, x_bs, sigma_bs, phi_bs

  ! Lund rescaling inflow (parsed, not ported — see check_supported)
  public :: i_rescale, delta_inlet, T_resc

  ! LES / wall-model flags (parsed, not ported — see check_supported)
  public :: LES_model, iwall_model, iwall_model_nut, frac_vis_wall_model
  public :: istress_model, Dirichlet_nu_t

  ! slip-length wall-model parameters ('alphas' keyword)
  public :: alpha_mean_x, alpha_mean_y, alpha_mean_z, alpha_std, freq_mult

  ! random-IC box size and y-stretching ('boxsize' keyword)
  public :: Lx_rand, Ly_rand, Lz_rand, alpha_rand

  !--------------------------- module constants -----------------------------
  integer, parameter :: IO_TMP     = 10   !< scratch unit for the params file
  integer, parameter :: ndim       = 3    !< number of dimensions
  integer, parameter :: baselength = 256  !< length of the params-file name

  !--------------------------- module variables -----------------------------
  ! All defaults are the effective reference defaults: the reference relies
  ! on zero-initialized module storage for keywords that are absent from the
  ! input file (plus the explicit defaults set in read_input_parameters of
  ! input_output.f90, applied below before parsing).

  character(baselength) :: paramsfilename = ' ' !< input file from '-i <file>'

  ! face-point counts (from 'nxyz'); single rank: local = global
  integer :: nx_global = 0, ny_global = 0, nz_global = 0
  integer :: nx  = 0, ny  = 0, nz  = 0
  ! center-point counts: n*m = n* - 1
  integer :: nxm_global = 0, nym_global = 0, nzm_global = 0
  integer :: nxm = 0, nym = 0, nzm = 0
  ! center points + ghost cells: n*g = n*m + 2
  integer :: nxg_global = 0, nyg_global = 0, nzg_global = 0
  integer :: nxg = 0, nyg = 0, nzg = 0

  real(dp) :: nu   = 0.0_dp  !< kinematic viscosity
  real(dp) :: CFL  = 0.0_dp  !< CFL number; CFL < 0 means fixed dt = |CFL|
  real(dp) :: dPdx = 0.0_dp  !< imposed streamwise pressure gradient
  real(dp) :: dPdz = 0.0_dp  !< imposed spanwise   pressure gradient
  real(dp) :: Omega_z = 0.0_dp !< spanwise rotation (Rossby_plus is hard-zero
                               !! in the reference, so Omega_z ends up 0)

  integer :: itime_step = 0  !< 'RKscheme': 1 Euler (unsupported), 2 RK2, 3 RK3

  integer :: nsteps   = 0    !< total number of time steps
  integer :: nsave    = 0    !< snapshot output period (steps)
  integer :: nstats   = 0    !< statistics output period (steps)
  integer :: nmonitor = 0    !< monitor output period (steps)

  integer :: nstep_init       = 0    !< step-number offset for output names
  integer :: nstep_init_input = -45  !< 'init_step'; -45 sentinel = not given
                                     !! (restart then takes the step from the
                                     !! input flow-field file, io_mod)
  integer :: random_init = 0         !< 'init_rand': 1 = random/Blasius IC,
                                     !! 0 = read initial field from 'filein'

  character(200) :: filein  = ' '            !< input flow field (restart)
  character(200) :: fileout = ' '            !< output flow-field prefix
  character(200) :: file_inflow = ' '        !< Blasius / mean-profile file
  ! HIT plane inflow (inflow_flag=6): time-resolved y-z planes from a
  ! precursor HIT simulation, preprocessed by preprocess_planes.py (v2).
  ! (tbl-gpu working tree / legacy global.f90+params.f90 additions.)
  character(200) :: file_hit_planes = ' ' !< binary planes file ('hit_file')
  integer        :: N_buffer_hit = 1000   !< planes held on the GPU
  character(200) :: file_ygrid = ' '      !< optional wall-normal grid file

  ! Subvolume ("box") output for causal-analysis data campaigns: up to 8
  ! boxes of u,v,w saved as float32 every boxout_every steps (absolute-step
  ! gated so the cadence is continuous across chain restarts).
  integer        :: boxout_every = 0      !< steps between box dumps (0 = off)
  integer        :: boxout_start = 0      !< absolute step to begin output
  integer        :: n_boxout = 0          !< number of boxes (<= 8)
  character(200) :: dir_boxout = './boxout' !< output directory (must exist)
  integer        :: boxout_i0(8) = 0, boxout_i1(8) = 0 !< x-index window
  integer        :: boxout_is(8) = 0      !< x-index stride
  integer        :: boxout_jmax(8) = 0    !< y-index cap (1..jmax)

  character(200) :: file_temporal_inlet = ' '!< temporal inflow modes file;
                                             !! '**TO DO**' placeholder means
                                             !! "none" (pure Blasius inflow)

  integer :: inflow_boundary_flag = 0 !< 1 Blasius(+temporal modes) — only
                                      !! supported value; 2 file, 3/5 Lund,
                                      !! 4 Blasius+random
  integer :: top_boundary_flag    = 0 !< 0 impose velocity — only supported
                                      !! value; 1 Coleman2018, 2 Abe2017,
                                      !! 3 Falkner-Skan

  real(dp) :: Amplitude_perturbations = 0.0_dp !< 'Amplitude' (inflow_flag=4)

  ! blowing/suction top-BC shape parameters (top_flag = 1 or 2)
  real(dp) :: Vbs_max  = 0.0_dp !< 'Vmax'  maximum vertical velocity
  real(dp) :: x_bs     = 0.0_dp !< 'x0'    blowing-to-suction position
  real(dp) :: sigma_bs = 0.0_dp !< 'sigma' gaussian width
  real(dp) :: phi_bs   = 0.0_dp !< 'phi'   net-flux parameter

  ! Lund rescaling (inflow_flag = 3 or 5)
  integer  :: i_rescale   = 0      !< 'Lund_ix' recycling-plane grid index
  real(dp) :: delta_inlet = 0.0_dp !< 'Lund_deltai' delta99 at the inlet
  real(dp) :: T_resc      = 0.0_dp !< 'Lund_T' rescaling averaging time

  ! LES / wall models
  integer  :: LES_model       = 0      !< 'LES' 0 none, 1 Smag, 2 DSM
  integer  :: iwall_model     = 0      !< 'WM'  0 none, 1..15 wall models
  integer  :: iwall_model_nut = 0      !< 'WMnutflag' (WM = 9 only)
  real(dp) :: frac_vis_wall_model = 0.0_dp !< 'WMnut' viscous fraction
  integer  :: istress_model   = 0      !< 'TauwModel' 0 utau_ref, 1 log law
  integer  :: Dirichlet_nu_t  = 0      !< 'nutBC' 0 Neumann, 1 Dirichlet

  ! slip-length wall-model parameters ('alphas' = 5 values, optional)
  real(dp) :: alpha_mean_x = 0.0_dp
  real(dp) :: alpha_mean_y = 0.0_dp
  real(dp) :: alpha_mean_z = 0.0_dp
  real(dp) :: alpha_std    = 0.0_dp
  real(dp) :: freq_mult    = 0.0_dp

  ! random-IC domain: [Lx Ly Lz alpha] from 'boxsize'
  real(dp) :: Lx_rand    = 0.0_dp !< domain length in x
  real(dp) :: Ly_rand    = 0.0_dp !< domain height in y
  real(dp) :: Lz_rand    = 0.0_dp !< domain width  in z
  real(dp) :: alpha_rand = 0.0_dp !< y-stretching factor (1 = uniform)

contains

  !--------------------------------------------------------------------------
  !> Startup banner (initialization.f90, subroutine initialize).
  subroutine print_banner()
    write(*,*) '----------------------------------------------------------------------'
    write(*,*) ' '
    write(*,*) '    My boundary layer ^^, CUDA Fortran single-GPU port                '
    write(*,*) ' '
    write(*,*) '----------------------------------------------------------------------'
  end subroutine print_banner

  !--------------------------------------------------------------------------
  !> Parse the command line and the '.turbb' input file, then derive the
  !! grid-point counts. Faithful transcription of params.f90:initparams
  !! (keyword names, parse order, error messages, directory creation) plus
  !! the pre-parse defaults and Omega_z from input_output.f90:
  !! read_input_parameters. MPI broadcasts are dropped (single rank).
  subroutine read_input_parameters()

    logical :: f, f1
    integer :: i, Nout, nslices_z
    character(len=32) :: arg
    integer  :: nxyz(ndim)         ! [nx ny nz]
    real(dp) :: boxsize(ndim+1)    ! [Lx Ly Lz alpha]
    real(dp) :: alphas(ndim+3)     ! alpha_mean_[xyz], std, freq_mult (5 read)
    real(dp) :: Rossby_plus, utau_

    write(*,*) 'reading parameters...'

    ! defaults set before parsing (input_output.f90:read_input_parameters)
    nstep_init              = 0
    nstep_init_input        = -45
    Rossby_plus             = 0.0_dp
    Amplitude_perturbations = 0.0_dp

    nxyz    = 0
    boxsize = 0.0_dp
    alphas  = 0.0_dp

    !---------------- read name of parameters file ('-i <file>') ------------
    f = .false.
    i = 0
    do
      call get_command_argument(i, arg)
      if (len_trim(arg) == 0) exit

      select case (arg)
      case ('-i')
        call get_command_argument(i+1, paramsfilename)
        f = .true.   ! input file found
      end select
      i = i + 1
    end do
    if (f .eqv. .false.) then
      stop 'No input file provided. Use *.exe -i <input_file>'
    end if

    inquire(file=paramsfilename, exist=f)
    if (f .eqv. .false.) then
      write(*,*) 'input file: ', adjustl(trim(paramsfilename)), ' does not exist'
      stop
    end if

    !---------------- keyword parsing (order as in initparams) --------------
    call get_dbl('nu'           , nu              , f)
    call get_dbl('CFL'          , CFL             , f)

    call get_dbl('Vmax'         , Vbs_max         , f)
    call get_dbl('x0'           , x_bs            , f)
    call get_dbl('sigma'        , sigma_bs        , f)
    call get_dbl('phi'          , phi_bs          , f)

    call get_int('nsteps'       , nsteps          , f)
    call get_int('nsave'        , nsave           , f)
    call get_int('nstats'       , nstats          , f)
    call get_int('nmonitor'     , nmonitor        , f)

    call get_int_arr('nxyz'     , nxyz , ndim     , f)
    call get_dbl_arr('boxsize'  , boxsize, ndim+1 , f)

    call get_int('inflow_flag'  , inflow_boundary_flag, f)
    call get_int('top_flag'     , top_boundary_flag, f)

    call get_str('inflow_file'  , file_inflow     , f)
    if (f .eqv. .false.) then
      stop ' ERROR: you must specify a BL profile in inflow_file '
    end if

    call get_int('Lund_ix'      , i_rescale       , f)
    call get_dbl('Lund_deltai'  , delta_inlet     , f)
    call get_dbl('Lund_T'       , T_resc          , f)

    call get_dbl('dPdx'         , dPdx            , f)
    call get_dbl('dPdz'         , dPdz            , f)

    call get_int('LES'          , LES_model       , f)
    call get_int('WM'           , iwall_model     , f)
    call get_int('TauwModel'    , istress_model   , f)
    call get_int('nutBC'        , Dirichlet_nu_t  , f)

    ! wall model.  NOTE: the reference passes the real(8) variable
    ! frac_vis_wall_model to get_int (an implicit-interface type mismatch);
    ! with an explicit interface that is illegal, so 'WMnut' is read with
    ! get_dbl here. WM > 0 is rejected by check_supported anyway.
    call get_int('WMnutflag'    , iwall_model_nut , f)
    call get_dbl('WMnut'        , frac_vis_wall_model, f)

    call get_dbl('Amplitude'    , Amplitude_perturbations, f)

    call get_int('init_step'    , nstep_init_input, f)
    call get_int('init_rand'    , random_init     , f)
    call get_str('filein'       , filein          , f1)
    if ( (random_init == 1) .and. (f1 .eqv. .true.) ) then
      stop ' ERROR: init_rand = 1, but filein exists'
    end if
    ! Assign initial step number if given (-45 = sentinel "not given";
    ! a restart then takes the step count from the input flow field, io_mod)
    if (nstep_init_input /= -45) then
      nstep_init = nstep_init_input
    end if

    call get_str('fileout'      , fileout         , f)

    ! Create the parent directory of the output files if it does not exist.
    ! Faithful to initparams: scan backwards for the last '/', mkdir -p the
    ! prefix. The reference loop walks past index 1 (undefined behavior) if
    ! fileout has no '/'; here that case stops with a clear message instead.
    Nout = len(trim(adjustl(fileout)))  ! length of fileout string

    i = Nout
    f = .false.
    do while (f .eqv. .false.)
      if (i < 1) then
        stop ' ERROR: fileout must contain a directory component, e.g. ./data/BL '
      end if
      if (fileout(i:i) == '/') then
        f = .true.
      end if
      i = i - 1
    end do

    write(*,*) fileout, Nout, fileout(1:i)

    inquire(file=trim(adjustl(fileout(1:i))), exist=f1)
    if (f1 .eqv. .false.) then
      call execute_command_line('mkdir -p '//trim(adjustl(fileout(1:i))))
    end if

    call get_str('timeinflow_file', file_temporal_inlet, f)

    ! HIT plane inflow (inflow_flag=6)
    call get_str('hit_file'     , file_hit_planes , f)
    if (.not. f) file_hit_planes = ''
    call get_int('N_buffer_hit' , N_buffer_hit    , f)
    if (.not. f) N_buffer_hit = 1000
    ! optional wall-normal grid from file (blended-sinh)
    call get_str('ygrid_file'   , file_ygrid      , f)
    if (.not. f) file_ygrid = ''

    ! subvolume (box) output for causal-analysis campaigns
    call get_int('boxout_every' , boxout_every    , f)
    if (.not. f) boxout_every = 0
    call get_int('boxout_start' , boxout_start    , f)
    if (.not. f) boxout_start = 0
    call get_int('boxout_n'     , n_boxout        , f)
    if (.not. f) n_boxout = 0
    call get_str('boxout_dir'   , dir_boxout      , f)
    if (.not. f) dir_boxout = './boxout'
    if (n_boxout > 8) stop ' ERROR: boxout_n > 8'
    if (n_boxout > 0) then
       call get_int_arr('boxout_i0'  , boxout_i0  , n_boxout, f)
       if (.not. f) stop ' ERROR: boxout_n>0 needs boxout_i0'
       call get_int_arr('boxout_i1'  , boxout_i1  , n_boxout, f)
       if (.not. f) stop ' ERROR: boxout_n>0 needs boxout_i1'
       call get_int_arr('boxout_is'  , boxout_is  , n_boxout, f)
       if (.not. f) stop ' ERROR: boxout_n>0 needs boxout_is'
       call get_int_arr('boxout_jmax', boxout_jmax, n_boxout, f)
       if (.not. f) stop ' ERROR: boxout_n>0 needs boxout_jmax'
    end if

    call get_int('RKscheme'     , itime_step      , f)

    call get_dbl_arr('alphas'   , alphas, ndim+2  , f)
    if (f .eqv. .true.) then
      alpha_mean_x = alphas(1)
      alpha_mean_y = alphas(2)
      alpha_mean_z = alphas(3)
      alpha_std    = alphas(4)
      freq_mult    = alphas(5)
    end if

    nx_global = nxyz(1)
    ny_global = nxyz(2)
    nz_global = nxyz(3)

    Lx_rand    = boxsize(1)
    Ly_rand    = boxsize(2)
    Lz_rand    = boxsize(3)
    alpha_rand = boxsize(4)

    !------------- nominal rotation (input_output.f90, verbatim) ------------
    utau_   = dPdx**0.5d0
    Omega_z = Rossby_plus*utau_/2d0

    !------------- derived point counts (initialization.f90, nprocs=1) ------
    ! number of interior z-planes (single rank owns all of them)
    nslices_z = nz_global - 2
    if (nslices_z < 1) stop 'Error: nslices_z must be at least 1'

    ! face points (local = global on a single rank)
    nx = nx_global
    ny = ny_global
    nz = nz_global

    ! middle (center) points
    nxm_global = nx_global - 1
    nym_global = ny_global - 1
    nzm_global = nz_global - 1

    nxm = nx - 1
    nym = ny - 1
    nzm = nz - 1

    ! middle points + ghost cells
    nxg_global = nxm_global + 2
    nyg_global = nym_global + 2
    nzg_global = nzm_global + 2

    nxg = nxm + 2
    nyg = nym + 2
    nzg = nzm + 2

  end subroutine read_input_parameters

  !--------------------------------------------------------------------------
  !> Stop with a clear message if the input requests a feature that the
  !! Phase 2 CUDA port does not implement (docs/CUDA_PORT_DESIGN.md):
  !! LES, wall models, wall-stress models, inflow flags other than 1
  !! (Blasius + temporal modes), top flags other than 0 (imposed velocity),
  !! and Euler time stepping. Restart from 'filein' (init_rand = 0) IS
  !! supported (io_mod), but requires a file name.
  subroutine check_supported()

    if (LES_model < 0 .or. LES_model > 3) then
      write(*,*) 'ERROR: LES = ', LES_model, ' requested; supported: 0 (DNS),', &
                 ' 1 (constant Smagorinsky), 2 (dynamic, z-averaged),', &
                 ' 3 (dynamic, no average).'
      stop 'unsupported feature: LES model'
    end if

    if (iwall_model < 0 .or. iwall_model == 10 .or. iwall_model > 15) then
      write(*,*) 'ERROR: WM = ', iwall_model, ' is not a defined wall model', &
                 ' (supported: 0-9, 11-15; the reference dispatcher has no 10).'
      stop 'unsupported feature: wall model'
    end if

    if (istress_model /= 0 .and. istress_model /= 1) then
      write(*,*) 'ERROR: TauwModel = ', istress_model, &
                 ' requested; supported: 0 (utau_ref) and 1 (log law).'
      stop 'unsupported feature: wall-stress model'
    end if

    if (inflow_boundary_flag /= 1 .and. inflow_boundary_flag /= 3 .and. &
        inflow_boundary_flag /= 5 .and. inflow_boundary_flag /= 6) then
      write(*,*) 'ERROR: inflow_flag = ', inflow_boundary_flag, &
                 ' requested, but only inflow_flag = 1 (Blasius + temporal modes),', &
                 ' 3/5 (Lund rescaling) and 6 (Blasius + HIT planes) are ported.'
      stop 'unsupported feature: inflow boundary condition'
    end if
    if (inflow_boundary_flag == 3 .or. inflow_boundary_flag == 5) then
      if (i_rescale < 2 .or. i_rescale > nx_global-1) then
        write(*,*) 'ERROR: inflow_flag = 3/5 needs 2 <= Lund_ix <= nx-1, got ', i_rescale
        stop 'invalid Lund_ix'
      end if
      if (T_resc <= 0.0_dp) then
        write(*,*) 'ERROR: inflow_flag = 3/5 needs Lund_T > 0, got ', T_resc
        stop 'invalid Lund_T'
      end if
      if (delta_inlet <= 0.0_dp) then
        write(*,*) 'ERROR: inflow_flag = 3/5 needs Lund_deltai > 0, got ', delta_inlet
        stop 'invalid Lund_deltai'
      end if
      if (inflow_boundary_flag == 5 .and. len_trim(file_inflow) == 0) then
        write(*,*) 'ERROR: inflow_flag = 5 requires inflow_file (turbulent mean profile).'
        stop 'missing inflow_file'
      end if
    end if
    if (inflow_boundary_flag == 6 .and. len_trim(file_hit_planes) == 0) then
      write(*,*) 'ERROR: inflow_flag = 6 requires hit_file in the input.'
      stop 'missing hit_file'
    end if

    if (top_boundary_flag /= 0 .and. top_boundary_flag /= 1 .and. &
        top_boundary_flag /= 2 .and. top_boundary_flag /= 4) then
      write(*,*) 'ERROR: top_flag = ', top_boundary_flag, &
                 ' requested, but only top_flag = 0 (imposed velocity),', &
                 ' 1/2 (blowing-suction lid) and 4 (zero-shear tangential) are ported.'
      stop 'unsupported feature: top boundary condition'
    end if
    if ((top_boundary_flag == 1 .or. top_boundary_flag == 2) .and. &
        sigma_bs == 0.0_dp) then
      write(*,*) 'ERROR: top_flag = 1/2 needs sigma /= 0 (Gaussian width).'
      stop 'invalid sigma for blowing/suction top'
    end if

    if (itime_step == 1) then
      write(*,*) 'ERROR: RKscheme = 1 (Euler) is not supported by the CUDA port; use RKscheme = 2 or 3.'
      stop 'unsupported feature: Euler time stepping'
    end if
    if (itime_step /= 2 .and. itime_step /= 3) then
      write(*,*) 'ERROR: RKscheme = ', itime_step, ' is invalid; use RKscheme = 2 or 3.'
      stop 'invalid RKscheme'
    end if

    if (random_init /= 1 .and. len_trim(filein) == 0) then
      write(*,*) 'ERROR: init_rand = 0 (restart) but no filein was given.'
      stop 'missing filein for restart'
    end if

    ! ---- multi-rank (P5.1) restrictions ----
    if (nprocs > 1) then
      ! legacy z-slab constraints (initialization.f90:56-59)
      if (mod(nx_global, 2) /= 0)          stop 'Error: nx must be even for MPI'
      if (mod(nz_global, 2) /= 0)          stop 'Error: nz must be even for MPI'
      if (mod(nz_global-2, nprocs) /= 0)   stop 'nz-2 should be divisible by nprocs'
    end if

  end subroutine check_supported

  !--------------------------------------------------------------------------
  !> Echo of the run parameters, ported loosely from monitor.f90:summary.
  !! Only the lines that depend exclusively on input parameters are printed
  !! here; the grid-coordinate lines (Lx/Ly/Lz, xg(1)... etc.) belong to the
  !! module that owns the grids and are printed by the io/monitor module.
  subroutine param_summary()

    write(*,*) '------------------------------------------------------------'
    write(*,*) '              Summary of initial parameters                 '
    write(*,*) ' '

    write(*,*) ' '
    write(*,*) 'Single GPU, no MPI (CUDA Fortran port)'

    write(*,*) ' '
    if     ( itime_step == 1 ) then
      write(*,*) 'Numerical integration: Explicit Euler'
    elseif ( itime_step == 2 ) then
      write(*,*) 'Numerical integration: Explicit RK2'
    else
      write(*,*) 'Numerical integration: Explicit RK3'
    end if

    write(*,*) ' '
    write(*,*) 'Input  file : ', trim(filein)
    write(*,*) 'Output file : ', trim(fileout)

    write(*,*) ' '
    write(*,*) 'Inflow parameters: '
    if     ( inflow_boundary_flag == 1 ) then
      write(*,*) '     Blasius profile + temporal perturbations'
      write(*,*) '     Blasius  file  : ', trim(file_inflow)
      write(*,*) '     temporal file  : ', trim(file_temporal_inlet)
    elseif ( inflow_boundary_flag == 2 ) then
      write(*,*) '     Blasius profile from file + temporal perturbations'
    elseif ( inflow_boundary_flag == 3 ) then
      write(*,*) '     Lunds rescaling'
    elseif ( inflow_boundary_flag == 4 ) then
      write(*,*) '     Blasius profile + temporal perturbations + random'
    elseif ( inflow_boundary_flag == 5 ) then
      write(*,*) '     Lund fixed mean and rescaled fluctuations'
    end if

    write(*,*) ' '
    write(*,*) 'Top bc parameters: '
    if     ( top_boundary_flag == 0 ) then
      write(*,*) '    Impose velocity'
    elseif ( top_boundary_flag == 1 ) then
      write(*,*) '    Coleman 2018'
    elseif ( top_boundary_flag == 2 ) then
      write(*,*) '    Abe 2017'
    elseif ( top_boundary_flag == 3 ) then
      write(*,*) '    Falkner-Skan velocity'
    end if

    write(*,*) ' '
    write(*,*) 'nu        :', nu
    write(*,*) 'CFL       :', CFL
    if ( CFL < 0.0_dp ) write(*,*) 'fixed dt  :', -CFL

    write(*,*) ' '
    write(*,*) 'dPdx    :', dPdx
    write(*,*) 'dPdz    :', dPdz
    write(*,*) 'Omega_z :', Omega_z

    write(*,*) ' '
    write(*,*) 'nsteps   :', nsteps
    write(*,*) 'nsave    :', nsave
    write(*,*) 'nstats   :', nstats
    write(*,*) 'nmonitor :', nmonitor

    write(*,*) ' '
    write(*,*) 'nx,nxg,nxm :', nx_global, nxg_global, nxm_global
    write(*,*) 'ny,nyg,nym :', ny_global, nyg_global, nym_global
    write(*,*) 'nz,nzg,nzm :', nz_global, nzg_global, nzm_global

    write(*,*) ' '
    if     ( LES_model == 1 ) then
      write(*,*) 'LES model : constant coefficient Smagorinsky'
    elseif ( LES_model == 2 ) then
      write(*,*) 'LES model : dynamic Smagorinsky, z-averaged'
    elseif ( LES_model == 3 ) then
      write(*,*) 'LES model : dynamic Smagorinsky, no-average'
    else
      write(*,*) 'No LES model'
    end if

    write(*,*) ' '
    if ( iwall_model == 0 ) then
      write(*,*) 'No wall model'
    else
      write(*,*) 'Wall model :', iwall_model
    end if

    write(*,*) ' '
    if     ( istress_model == 0 ) then
      write(*,*) 'No stress model : using utau_ref'
    elseif ( istress_model == 1 ) then
      write(*,*) 'Stress model : law-of-the-wall'
    end if

    write(*,*) ' '
    write(*,*) 'Separation parameters (Vtop)'
    write(*,*) '        Vbs_max  :', Vbs_max
    write(*,*) '        sigma_bs :', sigma_bs
    write(*,*) '        phi_bs   :', phi_bs
    write(*,*) '        x_bs     :', x_bs

    write(*,*) ' '
    write(*,*) '------------------------------------------------------------'

  end subroutine param_summary

  !==========================================================================
  ! '.turbb' key = value parser — faithful transcription of miscel.f90.
  !==========================================================================

  !> Scan the parameters file for `envvarname = value` and return the value
  !! string, with any trailing '!' comment stripped. Returns blank if the
  !! keyword is absent. Faithful to miscel.f90:getvar — the value keeps the
  !! characters immediately after '=' unshifted (trailing trim only, NO
  !! adjustl); downstream placeholder checks rely on this.
  subroutine getvar(envvarname, string)

    character(len=*),   intent(in)  :: envvarname
    character(len=256), intent(out) :: string

    character(len=256) :: line_str
    integer            :: eq_pos, ios
    logical            :: exists

    string = ' '

    ! check parameters file for the variable
    inquire(file=trim(adjustl(paramsfilename)), exist=exists)

    if (exists) then
      open(IO_TMP, file=trim(adjustl(paramsfilename)), status='old', iostat=ios)
      if (ios == 0) then
        do
          read(IO_TMP, '(a)', iostat=ios) line_str
          if (ios /= 0) exit
          eq_pos = scan(line_str, '=')

          if (trim(adjustl(line_str(1:eq_pos-1))) == envvarname) then
            string = trim(line_str(eq_pos+1:))
            eq_pos = scan(string, '!')
            if (eq_pos /= 0) then
              string = trim(string(1:eq_pos-1))
            end if
            exit
          end if
        end do
      end if
      close(IO_TMP)
    end if

  end subroutine getvar

  !> Read one double for keyword envvarname; found=.false. leaves var as is.
  subroutine get_dbl(envvarname, var, found)
    character(len=*), intent(in)    :: envvarname
    real(dp),         intent(inout) :: var
    logical,          intent(out)   :: found
    character(len=256) :: string

    call getvar(envvarname, string)
    if (string /= ' ') then
      read(string, *) var
      found = .true.
    else
      found = .false.
    end if
  end subroutine get_dbl

  !> Read n doubles for keyword envvarname; found=.false. leaves var as is.
  subroutine get_dbl_arr(envvarname, var, n, found)
    character(len=*), intent(in)    :: envvarname
    integer,          intent(in)    :: n
    real(dp),         intent(inout) :: var(n)
    logical,          intent(out)   :: found
    character(len=256) :: string

    call getvar(envvarname, string)
    if (string /= ' ') then
      read(string, *) var(1:n)
      found = .true.
    else
      found = .false.
    end if
  end subroutine get_dbl_arr

  !> Read one integer for keyword envvarname; found=.false. leaves var as is.
  subroutine get_int(envvarname, var, found)
    character(len=*), intent(in)    :: envvarname
    integer,          intent(inout) :: var
    logical,          intent(out)   :: found
    character(len=256) :: string

    call getvar(envvarname, string)
    if (string /= ' ') then
      read(string, *) var
      found = .true.
    else
      found = .false.
    end if
  end subroutine get_int

  !> Read n integers for keyword envvarname; found=.false. leaves var as is.
  subroutine get_int_arr(envvarname, var, n, found)
    character(len=*), intent(in)    :: envvarname
    integer,          intent(in)    :: n
    integer,          intent(inout) :: var(n)
    logical,          intent(out)   :: found
    character(len=256) :: string

    call getvar(envvarname, string)
    if (string /= ' ') then
      read(string, *) var
      found = .true.
    else
      found = .false.
    end if
  end subroutine get_int_arr

  !> Read a string for keyword envvarname; found=.false. leaves var as is.
  !! The value is transferred with an '(a)' read, so leading blanks (if any)
  !! after '=' are preserved, exactly as in the reference.
  subroutine get_str(envvarname, var, found)
    character(len=*), intent(in)    :: envvarname
    character(len=*), intent(inout) :: var
    logical,          intent(out)   :: found
    character(len=256) :: string

    call getvar(envvarname, string)
    if (string /= ' ') then
      read(string, '(a)') var
      found = .true.
    else
      found = .false.
    end if
  end subroutine get_str

end module param_mod
