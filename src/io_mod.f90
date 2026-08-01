!------------------------------------------------------------------------------
! io_mod — snapshots, statistics, restart and monitor output (host-only)
!
! Part of the CUDA Fortran port of the boundary-layer DNS solver
! (see docs/CUDA_PORT_DESIGN.md). Single GPU, nprocs = 1, no MPI.
!
! Provenance (reference files, transcribed faithfully):
!   * input_output.f90   output_data        -> output_snapshot + write_restart
!                        read_input_data    -> read_restart
!   * statistics.f90     compute_statistics -> output_stats (computation)
!   * input_output.f90   output_statistics  -> output_stats (.stats.txt writer)
!   * monitor.f90        output_monitor     -> output_monitor
!   * interpolation.f90  interpolate_x/interpolate_y (the two variants that
!                        compute_statistics calls) -> private host helpers
!
! Design decisions for the port (details flagged in the integration report):
!   * The reference downloads the interpolated fields from the GPU and then
!     reduces on the host. Here the interpolations themselves are HOST loops
!     acting on the U,V,W,P host mirrors, which are refreshed from the device
!     via field_mod's fields_from_device() at gated steps only (same numbers:
!     the interpolation formulas are transcribed verbatim and evaluated in
!     double precision either way).
!   * MPI_WTIME is replaced by system_clock (int64 count / count_rate).
!   * All MPI_Reduce calls collapse to their single-rank local values.
!   * The snapshot writer is byte-identical to the reference layout,
!     INCLUDING the nu_t block (identically zero in this DNS-only port; the
!     reference writes nu_t = 0 for LES_model = 0) and the full P block.
!   * The reference's '.restart' file is a SYMLINK to the latest snapshot,
!     recreated after every snapshot write; write_restart() reproduces that
!     (including the reference's inquire/rm/ln -s sequence).
!
! Assumed external interfaces (modules not owned by this file; the design doc
! assigns these symbols there but they were not on disk when io_mod was
! written — integrator must verify):
!   * timestep_mod (timestep_mod.cuf): public host scalars
!         real(dp) :: t            ! current simulation time
!         real(dp) :: dt           ! current time step
!         real(dp) :: dt_min_cfl   ! min CFL-limited dt (printed when CFL < 0)
!     (reference: global.f90 declares t, dt, dt_min_cfl; time_integration.f90
!     advances t and sets dt/dt_min_cfl, so the port's timestep_mod owns them.)
!   * poisson_mod (poisson_mod.cuf): public subroutine
!         check_divergence(max_divergence)   ! real(dp), intent(out)
!     (reference: projection.f90, same name and signature.)
!
! Driver contract (main.f90):
!   * call io_init() once after grid_generate() + fields_allocate() and before
!     the time loop (starts the elapsed-time clock, like the reference's
!     time1 = MPI_WTIME() at the end of initialize).
!   * inside the loop call, in this order (reference main.f90 order):
!       output_stats(istep); output_monitor(istep); output_snapshot(istep)
!     so that Cf and the inflow Reynolds numbers printed by the monitor are
!     the ones computed at the same step (both trigger at istep == 1).
!   * for a restart run (init_rand = 0): call read_restart() after
!     grid_generate() + fields_allocate() and BEFORE grid_to_device() /
!     fields_to_device(); read_restart only fills host arrays and t.
!------------------------------------------------------------------------------
module io_mod

  use iso_fortran_env, only: int32, int64
  use ieee_arithmetic, only: ieee_is_nan
  use precision_mod,   only: dp
  use param_mod,       only: nu, CFL, dPdx, Omega_z, alpha_std,               &
                             nsave, nstats, nmonitor,                          &
                             nstep_init, nstep_init_input,                     &
                             filein, fileout,                                  &
                             nx_global, ny_global, nz_global,                  &
                             nxm_global, nym_global, nzm_global,               &
                             nyg_global, nzg_global,                           &
                             boxout_every, boxout_start, n_boxout,             &
                             dir_boxout, boxout_i0, boxout_i1, boxout_is,      &
                             boxout_jmax, inflow_boundary_flag, LES_model
  use grid_mod,        only: nx, ny, nz, nxg, nyg, nzg,                        &
                             x, y, z, xm, ym, zm, yg,                          &
                             weight_y_0, weight_y_1,                           &
                             z_global, zm_global, k1, k2, kg1, nzm
  use field_mod,       only: U, V, W, P, fields_from_device
  use lund_inflow_mod, only: Umean_resc_T, Vmean_resc_T, &
                             Umean_resc_To, Vmean_resc_To
  use les_mod,         only: nu_t, les_nut_to_host
  use timestep_mod,    only: t, dt, dt_min_cfl
  use poisson_mod,     only: check_divergence
  use mpi_mod,         only: nprocs, myid, is_root, allreduce_max, &
                             allreduce_sum, allreduce_sum_arr, &
                             gather_slabs_host, gather_planes_r4

  implicit none
  private

  !----------------------------- public API ---------------------------------
  public :: io_init          ! subroutine io_init()
                             !   allocate statistics work arrays, start the
                             !   elapsed-time clock; call once after setup
  public :: output_snapshot  ! subroutine output_snapshot(istep)
                             !   integer, intent(in) :: istep
                             !   gated on mod(istep,nsave)==0: write the binary
                             !   snapshot fileout.XXXXXXXX and refresh the
                             !   '.restart' symlink
  public :: output_stats     ! subroutine output_stats(istep)
                             !   integer, intent(in) :: istep
                             !   gated on mod(istep,nstats)==0 .or. istep==1:
                             !   compute z-averaged statistics + Cf(x) and
                             !   write fileout.XXXXXXXX.stats.txt
  public :: output_monitor   ! subroutine output_monitor(istep)
                             !   integer, intent(in) :: istep
                             !   gated on mod(istep,nmonitor)==0 .or. istep==1:
                             !   print the monitor block (verbatim formats)
  public :: output_boxes     ! subroutine output_boxes(istep)
                             !   subvolume float32 dumps for causal-analysis
                             !   campaigns (self-gated on the ABSOLUTE step)
  public :: write_restart    ! subroutine write_restart()
                             !   recreate the fileout.restart symlink pointing
                             !   at the last snapshot written (no-op before the
                             !   first snapshot); called by output_snapshot
  public :: read_restart     ! subroutine read_restart()
                             !   read the binary snapshot 'filein' (grids,
                             !   U,V,W, t, and the init step) into host arrays

  !--------------------------- module constants -----------------------------
  ! dPdy exists in the reference global.f90 but is never read from the input
  ! nor assigned in the DNS path (always 0); the monitor still prints it.
  real(dp), parameter :: dPdy = 0.0_dp

  !--------------------------- module state ---------------------------------
  logical :: io_initialized = .false.

  ! elapsed-time clock (reference: time1/time2 via MPI_WTIME; here
  ! system_clock). time1 is (re)set at io_init and after every monitor print.
  real(dp) :: time1 = 0.0_dp

  ! host-mirror refresh bookkeeping: fields are downloaded from the device at
  ! most once per step, shared by stats / monitor / snapshot at gated steps.
  integer :: last_refresh_step = -huge(1)

  ! last snapshot written (target of the '.restart' symlink)
  character(200) :: last_snapshot_fname = ' '
  logical        :: have_snapshot = .false.

  !------------------- statistics state (statistics.f90) --------------------
  ! z-averaged first/second moments at (x-face, y-face) points, and wall
  ! quantities; shapes from initialization.f90:487-491 (single rank).
  real(dp), allocatable :: Umean(:,:),  Vmean(:,:),  Wmean(:,:),  Pmean(:,:)
  real(dp), allocatable :: U2mean(:,:), V2mean(:,:), W2mean(:,:), UVmean(:,:)
  real(dp), allocatable :: P2mean(:,:), nu_t_mean(:,:)
  real(dp), allocatable :: Cf(:), dUdy_wall(:), UV_wall(:)
  real(dp), allocatable :: Uaux_1(:), Uaux_2(:)
  real(dp), allocatable :: utau_wall(:), utau_wall_T(:)

  ! interpolation work arrays (reference term/term_1/term_2, host-only here;
  ! shapes from initialization.f90:141-143, nzm+2 = nzg at nprocs = 1)
  real(dp), allocatable :: term(:,:,:), term_1(:,:,:), term_2(:,:,:)

  ! inflow Reynolds numbers computed by output_stats, printed by the monitor
  real(dp) :: Rex_inlet         = 0.0_dp
  real(dp) :: Retheta_inlet     = 0.0_dp
  real(dp) :: Redelta_inlet     = 0.0_dp
  real(dp) :: delta99_inlet_ins = 0.0_dp

contains

  !---------------------------------------------------------------------------
  ! io_init — allocate the statistics work arrays and start the elapsed-time
  ! clock. Requires grid_generate() (sizes nx..nzg set). Idempotent.
  !
  ! Reference: initialization.f90:486-501 (statistics allocation + zeroing),
  ! 740-748 (utau arrays), 755 (time1 = MPI_WTIME()).
  !---------------------------------------------------------------------------
  subroutine io_init()

    if (io_initialized) return

    allocate ( Cf(nx), dUdy_wall(nx), UV_wall(nx) )
    allocate ( Uaux_1(nx), Uaux_2(nx) )
    allocate (  Umean(nx,ny),  Vmean(nx,ny),  Wmean(nx,ny),  Pmean(nx,ny) )
    allocate ( U2mean(nx,ny), V2mean(nx,ny), W2mean(nx,ny), UVmean(nx,ny) )
    allocate ( P2mean(nx,ny) )
    allocate ( nu_t_mean(nxg,nyg) )
    allocate ( utau_wall(nx), utau_wall_T(nx) )
    allocate ( term  (nxg,nyg,nzg) )
    allocate ( term_1(nxg,nyg,nzg) )
    allocate ( term_2(nxg,nyg,nzg) )

    Umean     = 0.0_dp
    Vmean     = 0.0_dp
    Wmean     = 0.0_dp
    Pmean     = 0.0_dp
    U2mean    = 0.0_dp
    V2mean    = 0.0_dp
    W2mean    = 0.0_dp
    UVmean    = 0.0_dp
    P2mean    = 0.0_dp
    nu_t_mean = 0.0_dp   ! DNS-only port: nu_t == 0 identically (LES_model = 0)

    Cf        = 0.0_dp
    dUdy_wall = 0.0_dp
    UV_wall   = 0.0_dp
    Uaux_1    = 0.0_dp
    Uaux_2    = 0.0_dp

    utau_wall   = 0.0_dp   ! actual utau at the wall (initialization.f90:747)
    utau_wall_T = 0.0_dp   ! averaged actual utau at the wall

    term   = 0.0_dp
    term_1 = 0.0_dp
    term_2 = 0.0_dp

    ! start measure time (reference: time1 = MPI_WTIME() at end of initialize)
    time1 = wtime()

    io_initialized = .true.

  end subroutine io_init

  !---------------------------------------------------------------------------
  ! wtime — wall-clock seconds from system_clock (replaces MPI_WTIME; the
  ! port has no MPI). int64 count avoids wrap; only differences are used.
  !---------------------------------------------------------------------------
  function wtime() result(secs)
    real(dp) :: secs
    integer(int64) :: cnt, rate
    call system_clock(cnt, rate)
    secs = real(cnt, dp)/real(rate, dp)
  end function wtime

  !---------------------------------------------------------------------------
  ! refresh_host_fields — download U,V,W,P from the device into the host
  ! mirrors, at most once per step (stats/monitor/snapshot share the copy,
  ! like the single gated '!$acc update self(U,V,W,P)' in the reference
  ! main.f90:99-101).
  !---------------------------------------------------------------------------
  subroutine refresh_host_fields(istep)
    integer, intent(in) :: istep
    if (istep /= last_refresh_step) then
      call fields_from_device()
      last_refresh_step = istep
    end if
  end subroutine refresh_host_fields

  !===========================================================================
  ! (a) Binary snapshot writer + restart symlink
  !===========================================================================

  !---------------------------------------------------------------------------
  ! output_snapshot — write the binary snapshot, byte-identical to the
  ! reference output_data (input_output.f90:496-680) at nprocs = 1.
  !
  ! Gate: mod(istep,nsave) == 0 (as reference). File name:
  ! fileout.XXXXXXXX with XXXXXXXX = I0.8 of istep + nstep_init.
  !
  ! Record layout (stream access, unformatted, native little endian):
  !   t, nu                      2 x real64
  !   -73, istep+nstep_init      2 x int32   (magic marker + step number)
  !   nx,  x(nx)                 int32, nx  x real64
  !   ny,  y(ny)                 int32, ny  x real64
  !   nz,  z(nz)                 int32, nz  x real64
  !   nxm, xm(nxm)               int32, nxm x real64
  !   nym, ym(nym)               int32, nym x real64
  !   nzm, zm(nzm)               int32, nzm x real64
  !   nx,nyg,nzg,   U(:,:,1:nzg-1)   3 x int32, nx*nyg*(nzg-1)  x real64
  !   nxg,ny,nzg,   V(:,:,1:nzg-1)   3 x int32, nxg*ny*(nzg-1)  x real64
  !   nxg,nyg,nz,   W(:,:,1:nz-1)    3 x int32, nxg*nyg*(nz-1)  x real64
  !   nxg,nyg,nzg,  nu_t(:,:,1:nzg-1) = 0      (DNS: eddy viscosity is zero;
  !                                             block kept for byte identity)
  !   nxg,nyg,nzg,  P(:,:,:)         3 x int32, nxg*nyg*nzg     x real64
  !                                  (FULL P, all nzg planes — reference
  !                                   writes 'P', not 'P(:,:,1:nzg-1)')
  ! then the '.restart' symlink is refreshed (write_restart).
  !---------------------------------------------------------------------------
  subroutine output_snapshot(istep)

    integer, intent(in) :: istep

    character(200) :: fname
    character(8)   :: ext
    integer        :: unit_s
    real(dp), allocatable :: nu_t_zero(:,:,:)

    if (.not. io_initialized) call io_init()

    if ( mod(istep,nsave) == 0 ) then

      call refresh_host_fields(istep)

      if (nprocs > 1) then
         call snapshot_write_global(istep)
         return
      end if

      write(ext,'(I0.8)') istep + nstep_init

      fname = trim(adjustl(fileout))//'.'//trim(adjustl(ext))
      write(*,*) 'writting ', trim(adjustl(fname))
      open(newunit=unit_s, file=fname, access='stream', form='unformatted', &
           action='write')

      ! metadata
      write(unit_s) t, nu
      write(unit_s) -73, istep + nstep_init  ! magic number for compatibility
                                             ! with older versions

      ! mesh (shape() of a rank-1 array = one default int32)
      write(unit_s) shape(x), x
      write(unit_s) shape(y), y
      write(unit_s) shape(z), z

      write(unit_s) shape(xm), xm
      write(unit_s) shape(ym), ym
      write(unit_s) shape(zm), zm

      ! U (single rank: reference rank-0 branch, planes 1:nzg-1)
      write(unit_s) nx, nyg, nzg
      write(unit_s) U(:,:,1:nzg-1)

      ! V
      write(unit_s) nxg, ny, nzg
      write(unit_s) V(:,:,1:nzg-1)

      ! W
      write(unit_s) nxg, nyg, nz
      write(unit_s) W(:,:,1:nz-1)

      ! nu_t — identically zero in this DNS-only port (the reference writes
      ! nu_t computed by compute_eddy_viscosity, which is 0 for LES_model=0);
      ! the block is written anyway so the file is byte-identical.
      write(unit_s) nxg, nyg, nzg
      if (LES_model > 0) then
        call les_nut_to_host()
        write(unit_s) nu_t(:,:,1:nzg-1)
      else
        allocate ( nu_t_zero(nxg,nyg,nzg-1) )
        nu_t_zero = 0.0_dp
        write(unit_s) nu_t_zero
        deallocate ( nu_t_zero )
      end if

      ! P — the reference writes the FULL array (all nzg planes)
      write(unit_s) nxg, nyg, nzg
      write(unit_s) P

      close(unit_s)

      ! remember the file and refresh the '.restart' symlink
      last_snapshot_fname = fname
      have_snapshot       = .true.
      call write_restart()

      ! save means for Lund's rescaling (reference input_output.f90:733-741,
      ! root only, inflow_flag = 3; NATIVE endianness — the reference's
      ! writer was native too, its big-endian read matched Intel-era files)
      if ( myid == 0 .and. inflow_boundary_flag == 3 ) then
        block
          integer :: unit_l
          open(newunit=unit_l, file=trim(adjustl(fname))//'.mean.rescaling', &
               access='stream', form='unformatted', action='write')
          write(unit_l) nyg_global
          write(unit_l) Umean_resc_T
          write(unit_l) Vmean_resc_T
          close(unit_l)
        end block
      end if

    end if

  end subroutine output_snapshot

  !---------------------------------------------------------------------------
  ! snapshot_write_global — P5.1 multi-rank snapshot: assemble the global
  ! fields on rank 0 (host slab gather; overlapping planes are identical
  ! because ghosts are exchanged every substep) and write the reference
  ! format at GLOBAL sizes. Byte-identical to a single-rank snapshot of the
  ! same state.
  !---------------------------------------------------------------------------
  subroutine snapshot_write_global(istep)

    integer, intent(in) :: istep

    character(200) :: fname
    character(8)   :: ext
    integer        :: unit_s
    real(dp), allocatable :: Ug(:,:,:), Vg(:,:,:), Wg(:,:,:), Pg(:,:,:)
    real(dp), allocatable :: nu_t_zero(:,:,:), nutg(:,:,:)

    if (is_root()) then
       allocate ( Ug(nx,  nyg, nzg_global) )
       allocate ( Vg(nxg, ny,  nzg_global) )
       allocate ( Wg(nxg, nyg, nz_global ) )
       allocate ( Pg(nxg, nyg, nzg_global) )
       if (LES_model > 0) allocate ( nutg(nxg, nyg, nzg_global) )
    else
       allocate ( Ug(1,1,1), Vg(1,1,1), Wg(1,1,1), Pg(1,1,1) )
       if (LES_model > 0) allocate ( nutg(1,1,1) )
    end if

    call gather_slabs_host(U, nx,  nyg, nzg, kg1, Ug, nzg_global)
    call gather_slabs_host(V, nxg, ny,  nzg, kg1, Vg, nzg_global)
    call gather_slabs_host(W, nxg, nyg, nz,  k1,  Wg, nz_global )
    call gather_slabs_host(P, nxg, nyg, nzg, kg1, Pg, nzg_global)
    if (LES_model > 0) then
       call les_nut_to_host()
       call gather_slabs_host(nu_t, nxg, nyg, nzg, kg1, nutg, nzg_global)
    end if

    if (is_root()) then
       write(ext,'(I0.8)') istep + nstep_init
       fname = trim(adjustl(fileout))//'.'//trim(adjustl(ext))
       write(*,*) 'writting ', trim(adjustl(fname))
       open(newunit=unit_s, file=fname, access='stream', form='unformatted', &
            action='write')

       write(unit_s) t, nu
       write(unit_s) -73, istep + nstep_init

       write(unit_s) nx, x
       write(unit_s) ny, y
       write(unit_s) nz_global, z_global
       write(unit_s) nx-1, xm
       write(unit_s) ny-1, ym
       write(unit_s) nz_global-1, zm_global

       write(unit_s) nx, nyg, nzg_global
       write(unit_s) Ug(:,:,1:nzg_global-1)
       write(unit_s) nxg, ny, nzg_global
       write(unit_s) Vg(:,:,1:nzg_global-1)
       write(unit_s) nxg, nyg, nz_global
       write(unit_s) Wg(:,:,1:nz_global-1)

       write(unit_s) nxg, nyg, nzg_global
       if (LES_model > 0) then
         write(unit_s) nutg(:,:,1:nzg_global-1)
       else
         allocate ( nu_t_zero(nxg,nyg,nzg_global-1) )
         nu_t_zero = 0.0_dp
         write(unit_s) nu_t_zero
         deallocate ( nu_t_zero )
       end if

       write(unit_s) nxg, nyg, nzg_global
       write(unit_s) Pg

       close(unit_s)

       last_snapshot_fname = fname
       have_snapshot       = .true.
       call write_restart()

       ! (audit fix) Lund companion also in the multi-rank path
       ! (reference output_data writes it for inflow_flag = 3 on root)
       if ( inflow_boundary_flag == 3 ) then
         block
           integer :: unit_l
           open(newunit=unit_l, file=trim(adjustl(fname))//'.mean.rescaling', &
                access='stream', form='unformatted', action='write')
           write(unit_l) nyg_global
           write(unit_l) Umean_resc_T
           write(unit_l) Vmean_resc_T
           close(unit_l)
         end block
       end if
    end if

    deallocate (Ug, Vg, Wg, Pg)
    if (LES_model > 0) deallocate (nutg)

  end subroutine snapshot_write_global

  !---------------------------------------------------------------------------
  ! write_restart — recreate the '<fileout>.restart' symlink pointing at the
  ! last snapshot written. Faithful port of the symlink block at the end of
  ! the reference output_data (input_output.f90:652-664): inquire whether the
  ! link resolves, 'rm' it if so, then 'ln -s <snapshot> <fileout>.restart'
  ! ('call system' replaced by execute_command_line). No-op until the first
  ! snapshot exists.
  !
  ! The link lives in the same directory as the snapshot, so the target must
  ! be the snapshot BASENAME (a cwd-relative path would dangle from that
  ! directory). 'ln -sf' overwrites so the link always tracks the newest
  ! file. (Bug fix imported from the tbl-gpu working tree.)
  !---------------------------------------------------------------------------
  subroutine write_restart()

    character(200) :: fname_symlnk
    character(400) :: string_link

    if (.not. have_snapshot) return

    fname_symlnk = trim(adjustl(fileout))//'.'//'restart'

    string_link = 'ln -sf '// &
         trim(adjustl(last_snapshot_fname( &
              index(last_snapshot_fname, '/', back=.true.)+1: )))// &
         ' '//trim(fname_symlnk)
    call execute_command_line( string_link )

  end subroutine write_restart

  !---------------------------------------------------------------------------
  ! output_boxes — subvolume ("box") output for causal-analysis campaigns.
  ! Verbatim port of legacy/input_output.f90 output_boxes (tbl-gpu working
  ! tree). Every boxout_every steps (gated on the ABSOLUTE step so the
  ! cadence survives chain restarts), writes u,v,w as float32 on up to 8
  ! index-window boxes to per-launch stream files:
  !     <dir_boxout>/box<b>_<first-abs-step>.bin
  ! File: int32 header (magic 20260709, version, id, i0, i1, is, jmax,
  ! npx, npy, npz, every, ncomp) + f8 x(sel), y(1:jmax), z(1:NKB), then per
  ! record: f8 t, int32 absstep, r4 u,v,w (npx, jmax, NKB). Components stay
  ! on their native staggered grids over the same index window.
  !---------------------------------------------------------------------------
  subroutine output_boxes(istep)

    integer, intent(in) :: istep

    ! z planes saved: one full period when the grid has it, else all
    ! interior planes (the header's npz field makes the format
    ! self-describing, so readers adapt automatically)
    integer :: NKB
    integer :: b, i, j, k, ii, absstep, npx, jm
    integer :: g0, g1, nploc, pl
    character(300) :: fname
    character(8)   :: ext
    integer,       save :: box_unit(8) = -1
    logical,       save :: box_open(8) = .false.
    logical,       save :: checked = .false.
    real(4), allocatable, save :: loc4(:,:,:,:)   ! (npx, jm, nploc, 3)
    real(4), allocatable, save :: glb4(:,:,:,:)   ! (npx, jm, NKB, 3) root

    if ( boxout_every <= 0 .or. n_boxout <= 0 ) return
    NKB = min( 256, nz_global - 2 )
    absstep = istep + nstep_init
    if ( absstep < boxout_start ) return
    if ( mod(absstep, boxout_every) /= 0 ) return

    if (.not. io_initialized) call io_init()

    ! this rank's window of the global box z-planes 1..NKB (faces):
    ! global face g = k1 + p - 1; overlapping planes are ghost-synced, so
    ! assembly-by-overwrite on the root is exact.
    g0 = max( 1, k1 )
    g1 = min( NKB, k2 )
    nploc = max( 0, g1 - g0 + 1 )

    if ( .not. checked ) then
       ii = 0; j = 0
       do b = 1, n_boxout
          if ( boxout_i1(b) > nx .or. boxout_jmax(b) > ny .or. NKB > nz_global ) &
               stop 'ERROR: boxout window exceeds grid'
          if ( boxout_is(b) < 1 .or. boxout_i0(b) < 1 ) stop 'ERROR: bad boxout spec'
          ii = max( ii, (boxout_i1(b)-boxout_i0(b))/boxout_is(b) + 1 )
          j  = max( j, boxout_jmax(b) )
       end do
       allocate( loc4(ii, j, max(nploc,1), 3) )   ! sized once to the largest box
       if (is_root()) allocate( glb4(ii, j, NKB, 3) )
       checked = .true.
    end if

    ! fresh fields on the host (shared with stats/monitor/snapshot)
    call refresh_host_fields(istep)

    do b = 1, n_boxout
       npx = (boxout_i1(b) - boxout_i0(b))/boxout_is(b) + 1
       jm  = boxout_jmax(b)

       if ( is_root() .and. .not. box_open(b) ) then
          write(ext,'(I0.8)') absstep
          fname = trim(dir_boxout)//'/box'//char(48+b)//'_'//ext//'.bin'
          open(newunit=box_unit(b), file=trim(fname), access='stream', &
               form='unformatted', action='write', status='replace')
          write(box_unit(b)) 20260709, 1, b, boxout_i0(b), boxout_i1(b), &
                             boxout_is(b), jm, npx, jm, NKB, boxout_every, 3
          write(box_unit(b)) x(boxout_i0(b):boxout_i1(b):boxout_is(b))
          write(box_unit(b)) y(1:jm)
          write(box_unit(b)) z_global(1:NKB)
          box_open(b) = .true.
          write(*,'(a,i2,a,a)') '   boxout: opened box ', b, ' -> ', trim(fname)
       end if

       ! local slab contribution (local plane pl = global g0..g1)
       do k = 1, nploc
          pl = g0 - k1 + k       ! local z index of global plane g0+k-1
          do j = 1, jm
             ii = 0
             do i = boxout_i0(b), boxout_i1(b), boxout_is(b)
                ii = ii + 1
                loc4(ii,j,k,1) = real( U(i,j,pl), 4 )
                loc4(ii,j,k,2) = real( V(i,j,pl), 4 )
                loc4(ii,j,k,3) = real( W(i,j,pl), 4 )
             end do
          end do
       end do

       ! assemble the global box on the root (identity copy at nprocs = 1)
       do k = 1, 3
          if (is_root()) then
             call gather_planes_r4(loc4(1:npx,1:jm,1:max(nploc,1),k), npx, jm, &
                                   nploc, g0, glb4(1:npx,1:jm,:,k), NKB)
          else
             call gather_planes_r4(loc4(1:npx,1:jm,1:max(nploc,1),k), npx, jm, &
                                   nploc, g0, loc4(1:npx,1:jm,1:1,k), 1)
          end if
       end do

       if (is_root()) then
          write(box_unit(b)) t, absstep
          write(box_unit(b)) glb4(1:npx,1:jm,1:NKB,1)
          write(box_unit(b)) glb4(1:npx,1:jm,1:NKB,2)
          write(box_unit(b)) glb4(1:npx,1:jm,1:NKB,3)
          flush(box_unit(b))
       end if
    end do

  end subroutine output_boxes

  !---------------------------------------------------------------------------
  ! read_restart — read a binary snapshot (grids, U, V, W, time, init step)
  ! from 'filein' into the HOST arrays. Port of read_input_data
  ! (input_output.f90:288-487) at nprocs = 1, minus MPI, minus the Lund
  ! '.mean.rescaling' companion file (read back for inflow_flag = 3).
  !
  ! Reads: t; the -73 marker + step number (taken as nstep_init when the
  ! input file gave no init_step, i.e. nstep_init_input == -45); the six
  ! grids (each checked against the input-file sizes with the reference's
  ! error messages) overwriting grid_mod's x,y,z,xm,ym,zm; then U,V,W.
  !
  ! FIELD PLANE COUNTS — deliberate deviation from the reference reader:
  ! the reference read_input_data at nprocs = 1 reads nzg z-planes per field
  ! and seeks with nzg-plane offsets, but its own writer stores only nzg-1
  ! planes (U,V) and nz-1 planes (W) — the single-rank restart path in the
  ! reference is internally inconsistent (it was only ever exercised through
  ! multi-rank files). This port reads exactly what output_snapshot writes:
  ! U(:,:,1:nzg-1), V(:,:,1:nzg-1), W(:,:,1:nz-1). The remaining last plane
  ! of each field is a periodic z-ghost/image plane and is rebuilt by the
  ! periodic-z boundary condition before the first RHS evaluation.
  !
  ! Post-conditions / driver duties: host U,V,W and grids are filled and t is
  ! set; the driver must afterwards call grid_to_device() and
  ! fields_to_device() (the reference's 'Uo = U' zero-step copy is subsumed
  ! by the per-step Uo save in timestep_mod).
  !---------------------------------------------------------------------------
  subroutine read_restart()

    integer(int32) :: nx_global_f, ny_global_f, nz_global_f
    integer(int32) :: nxm_global_f, nym_global_f, nzm_global_f, nn(3)
    real(dp)       :: nu_dummy
    character(200) :: err_msg
    integer        :: unit_r

    write(*,*) 'reading ', trim(adjustl(filein)), '...'
    open(newunit=unit_r, file=filein, access='stream', form='unformatted', &
         action='read')

    ! metadata
    read(unit_r) t, nu_dummy

    ! mesh
    read(unit_r) nx_global_f

    if ( nx_global_f == -73 ) then
      read(unit_r) ny_global_f  ! this is the time step, but we only use it
                                ! if we haven't provided an input init_step
      if ( nstep_init_input == -45 ) then
        nstep_init = ny_global_f
        write(*,*) 'Reading init step from file: ', nstep_init
      end if
      read(unit_r) nx_global_f  ! read the actual nx_global in file
    end if

    write(err_msg,'(A,I5,A,I5,A)') 'nx_f(', nx_global_f, ')/=nx(', nx_global, ')'
    if ( nx_global_f /= nx_global ) stop trim(err_msg)
    read(unit_r) x

    read(unit_r) ny_global_f
    write(err_msg,'(A,I5,A,I5,A)') 'ny_f(', ny_global_f, ')/=ny(', ny_global, ')'
    if ( ny_global_f /= ny_global ) stop trim(err_msg)
    read(unit_r) y

    read(unit_r) nz_global_f
    write(err_msg,'(A,I5,A,I5,A)') 'nz_f(', nz_global_f, ')/=nz(', nz_global, ')'
    if ( nz_global_f /= nz_global ) stop trim(err_msg)
    ! z record is GLOBAL; every rank reads it and slices its slab
    read(unit_r) z_global
    z = z_global(k1:k2)

    read(unit_r) nxm_global_f
    write(err_msg,'(A,I5,A,I5,A)') 'nxm_f(', nxm_global_f, ')/=nxm(', nxm_global, ')'
    if ( nxm_global_f /= nxm_global ) stop trim(err_msg)
    read(unit_r) xm

    read(unit_r) nym_global_f
    write(err_msg,'(A,I5,A,I5,A)') 'nym_f(', nym_global_f, ')/=nym(', nym_global, ')'
    if ( nym_global_f /= nym_global ) stop trim(err_msg)
    read(unit_r) ym

    read(unit_r) nzm_global_f
    write(err_msg,'(A,I5,A,I5,A)') 'nzm_f(', nzm_global_f, ')/=nzm(', nzm_global, ')'
    if ( nzm_global_f /= nzm_global ) stop trim(err_msg)
    read(unit_r) zm_global
    zm = zm_global(kg1:kg1+nzm-1)

    ! Field blocks. The file stores GLOBAL fields (nzg_global-1 z-planes
    ! for the z-centered U,V; nz_global-1 for the z-faced W). Every rank
    ! seeks directly to its own z-window (POS= stream addressing, Int64) —
    ! fully parallel reads, no MPI. Planes past the stored range (the top
    ! ghost planes) stay zero and are rebuilt by the first BC application,
    ! exactly like the single-rank read of 1:nzg-1.
    block
      integer(int64) :: pos_u, pos_v, pos_w, pu, pv, pw
      integer :: nu_pl, nv_pl, nw_pl
      inquire(unit_r, pos=pos_u)          ! start of the U size record
      pu = pos_u + 12_int64
      pos_v = pu + int(nx, int64)*int(nyg, int64)*int(nzg_global-1, int64)*8_int64
      pv = pos_v + 12_int64
      pos_w = pv + int(nxg, int64)*int(ny, int64)*int(nzg_global-1, int64)*8_int64
      pw = pos_w + 12_int64

      nu_pl = min( nzg, nzg_global-1 - kg1 + 1 )
      nv_pl = nu_pl
      nw_pl = min( nz, nz_global-1 - k1 + 1 )

      read(unit_r, pos=pu + int(kg1-1, int64)*int(nx, int64)*int(nyg, int64)*8_int64) &
           U(:,:,1:nu_pl)
      read(unit_r, pos=pv + int(kg1-1, int64)*int(nxg, int64)*int(ny, int64)*8_int64) &
           V(:,:,1:nv_pl)
      read(unit_r, pos=pw + int(k1-1, int64)*int(nxg, int64)*int(nyg, int64)*8_int64) &
           W(:,:,1:nw_pl)
    end block

    ! close file (nu_t and P blocks are not read, as in the reference)
    close(unit_r)

    ! sanity check
    if ( any( ieee_is_nan(U) ) ) stop 'Error U NaNs!'
    if ( any( ieee_is_nan(V) ) ) stop 'Error V NaNs!'
    if ( any( ieee_is_nan(W) ) ) stop 'Error W NaNs!'

    ! read Umean/Vmean for Lund recycling if the companion snapshot exists
    ! (reference input_output.f90:505-527, inflow_flag = 3). Missing file
    ! -> zeros, and lund_init_means falls back to the IC spanwise means.
    ! Every rank reads (no bcast needed; y is undecomposed). NATIVE
    ! endianness: matches this port's writer above.
    if ( inflow_boundary_flag == 3 ) then
      block
        integer :: unit_l, ios_l, nyg_resc
        logical :: have_resc
        inquire(file=trim(adjustl(filein))//'.mean.rescaling', exist=have_resc)
        if ( have_resc ) then
          write(*,*) 'reading ', trim(adjustl(filein))//'.mean.rescaling', '...'
          open(newunit=unit_l, file=trim(adjustl(filein))//'.mean.rescaling', &
               access='stream', form='unformatted', action='read', &
               status='old', iostat=ios_l)
          read(unit_l) nyg_resc
          if ( nyg_resc /= nyg_global ) then
            write(*,*) nyg_resc, nyg_global
            stop 'Error! nyg_resc/=nyg_global'
          end if
          read(unit_l) Umean_resc_To
          read(unit_l) Vmean_resc_To
          close(unit_l)
        else
          Umean_resc_To = 0.0_dp
          Vmean_resc_To = 0.0_dp
        end if
      end block
    end if

  end subroutine read_restart

  !===========================================================================
  ! (b) Statistics: z-averaged means, Cf(x), .stats.txt writer
  !===========================================================================

  !---------------------------------------------------------------------------
  ! output_stats — compute the on-the-fly statistics and write the
  ! .stats.txt file. Port of statistics.f90:compute_statistics (which calls
  ! input_output.f90:output_statistics at the end) at nprocs = 1.
  !
  ! Gate: mod(istep,nstats) == 0 .or. istep == 1 (as reference).
  !
  ! Interpolations (reference calls GPU interpolation.f90 routines and
  ! downloads the results; here the same formulas run as host loops on the
  ! refreshed host mirrors — identical arithmetic):
  !   term_1 <- interpolate_x(W, di=1)   x-faces average, W to x-centers
  !   term   <- interpolate_y(term_1, 2) weighted centers->y-faces (W at y-faces)
  !   term_1 <- interpolate_y(U, 2)      weighted centers->y-faces (U at y-faces)
  !   term_2 <- interpolate_x(V, 1)      x average (V at x-faces)
  ! Then z-averaged moments at (x-face, y-face) points, the wall derivative
  ! dU/dy from the first two U planes, the wall Reynolds stress, the (zero)
  ! SGS wall stress, Cf, and the inlet momentum/99% thicknesses and Reynolds
  ! numbers. All loop bodies verbatim from statistics.f90:61-201.
  !---------------------------------------------------------------------------
  subroutine output_stats(istep)

    integer, intent(in) :: istep

    integer  :: ii, jj, j, jref
    real(dp) :: Uinf, theta_inlet, Umean99, w0, w1, T_mean
    real(dp) :: tau_sgs_wall(nx)
    real(dp) :: temp_1d(ny)

    if (.not. io_initialized) call io_init()

    ! statistics computed at grid y -> U and W interpolated
    if ( mod(istep,nstats) == 0 .or. istep == 1 ) then

      call refresh_host_fields(istep)
      if (LES_model > 0) call les_nut_to_host()

      ! interpolate W in x and y -> term
      call interp_x_avg(W, term_1)
      call interp_y_weighted(term_1, term)

      ! interpolate U in y -> term_1
      call interp_y_weighted(U, term_1)

      ! interpolate V in x -> term_2
      call interp_x_avg(V, term_2)

      ! compute local statistics (statistics.f90:61-80; single rank, so the
      ! MPI_Reduce stage is the identity)
      do ii = 1, nx
        do jj = 1, ny

          Umean (ii,jj) = sum( term_1(ii, jj, 2:nzg-1) )
          Vmean (ii,jj) = sum( term_2(ii, jj, 2:nzg-1) )
          Wmean (ii,jj) = sum( term  (ii, jj, 1:nz-2 ) )

          U2mean(ii,jj) = sum( term_1(ii, jj, 2:nzg-1)**2.0_dp )
          V2mean(ii,jj) = sum( term_2(ii, jj, 2:nzg-1)**2.0_dp )
          W2mean(ii,jj) = sum( term  (ii, jj, 1:nz-2 )**2.0_dp )

          UVmean(ii,jj) = sum( term_1(ii,jj,2:nzg-1)*term_2(ii,jj,2:nzg-1) )

          Pmean (ii,jj) = sum( P(ii, jj, 2:nzg-1)        )
          P2mean(ii,jj) = sum( P(ii, jj, 2:nzg-1)**2.0_dp )

          if (LES_model > 0) then
            nu_t_mean(ii,jj) = sum( nu_t(ii, jj, 2:nzg-1) )
          end if

        end do
      end do

      ! P5.3: cross-rank reduction of the partial z-sums. The local sum
      ! ranges above partition the global z exactly (overlapping-slab
      ! bookkeeping: centers 2:nzg-1 and faces 1:nz-2 tile the globe), so
      ! summing the per-rank partials and dividing by the *_global counts
      ! below reproduces the single-rank statistics exactly.
      if (nprocs > 1) then
         call allreduce_sum_arr(Umean,  nx*ny)
         call allreduce_sum_arr(Vmean,  nx*ny)
         call allreduce_sum_arr(Wmean,  nx*ny)
         call allreduce_sum_arr(U2mean, nx*ny)
         call allreduce_sum_arr(V2mean, nx*ny)
         call allreduce_sum_arr(W2mean, nx*ny)
         call allreduce_sum_arr(UVmean, nx*ny)
         call allreduce_sum_arr(Pmean,  nx*ny)
         call allreduce_sum_arr(P2mean, nx*ny)
         if (LES_model > 0) call allreduce_sum_arr(nu_t_mean, nxg*nyg)
      end if

      ! z-averages (statistics.f90:122-135)
      Umean  = Umean /real( nzg_global-2, dp )
      Vmean  = Vmean /real( nzg_global-2, dp )
      Wmean  = Wmean /real(  nz_global-2, dp )

      U2mean = U2mean/real( nzg_global-2, dp )
      V2mean = V2mean/real( nzg_global-2, dp )
      W2mean = W2mean/real(  nz_global-2, dp )

      UVmean = UVmean/real( nzg_global-2, dp )

      Pmean  = Pmean /real( nzg_global-2, dp )
      P2mean = P2mean/real( nzg_global-2, dp )

      if (LES_model > 0) nu_t_mean = nu_t_mean/real( nzg_global-2, dp )

      ! mean derivative at the wall (statistics.f90:138-148)
      do ii = 1, nx
        Uaux_1(ii) = sum( U(ii,1,2:nzg-1) )
        Uaux_2(ii) = sum( U(ii,2,2:nzg-1) )
      end do
      if (nprocs > 1) then
         call allreduce_sum_arr(Uaux_1, nx)
         call allreduce_sum_arr(Uaux_2, nx)
      end if
      Uaux_1    = Uaux_1/real( nzg_global-2, dp )
      Uaux_2    = Uaux_2/real( nzg_global-2, dp )
      dUdy_wall = ( Uaux_2 - Uaux_1 )/( yg(2) - yg(1) )

      ! mean Reynolds stress at the wall
      UV_wall = UVmean(:,1)

      ! SGS stress at the wall (statistics.f90:154-156; nu_t_mean == 0 in
      ! this DNS-only port, expression kept verbatim)
      do ii = 1, nx
        tau_sgs_wall(ii) = 0.25_dp*( nu_t_mean(ii,1) + nu_t_mean(ii,2) + &
                           nu_t_mean(ii+1,1) + nu_t_mean(ii+1,2) )*dUdy_wall(ii)
      end do

      ! skin friction coefficient Cf = tau_w/(1/2*rho*U_inf^2)
      Uinf        = 1.0_dp ! U(1,nyg,1)
      Cf          = 2.0_dp*( -UV_wall + nu*dUdy_wall + tau_sgs_wall )/Uinf**2.0_dp
      utau_wall   = ( 0.5_dp*Uinf**2.0_dp*abs(Cf) )**0.5_dp
      if (istep == 1) utau_wall_T = utau_wall
      T_mean      = alpha_std
      utau_wall_T = dt/T_mean*utau_wall + (1.0_dp - dt/T_mean)*utau_wall_T

      ! momentum boundary layer thickness (statistics.f90:178-182)
      theta_inlet = 0.0_dp
      temp_1d     = Umean(1,:)/Uinf*(1.0_dp - Umean(1,:)/Uinf)
      do j = 2, ny
        theta_inlet = theta_inlet + 0.5_dp*( temp_1d(j) + temp_1d(j-1) )*( y(j)-y(j-1) )
      end do

      ! 99% boundary layer thickness (statistics.f90:185-196)
      Umean99 = 0.99_dp*Uinf
      jref    = 0
      do j = 1, ny
        if ( Umean(1,j) >= Umean99 ) then
          jref = j
          exit
        end if
      end do
      if ( jref < 2 ) jref = 2
      w1                = ( Umean99 - Umean(1,jref-1) )/( Umean(1,jref) - Umean(1,jref-1) )
      w0                = 1.0_dp - w1
      delta99_inlet_ins = w1*y(jref) + w0*y(jref-1)

      ! Reynolds numbers at inlet (statistics.f90:199-201; x == x_global)
      Rex_inlet     = x(1)*Uinf/nu
      Retheta_inlet = theta_inlet*Uinf/nu
      Redelta_inlet = delta99_inlet_ins*Uinf/nu

      ! write statistics (all ranks hold identical reduced values;
      ! only the root writes the file)
      if (is_root()) call write_stats_file(istep)

    end if

  end subroutine output_stats

  !---------------------------------------------------------------------------
  ! write_stats_file — the .stats.txt writer, verbatim port of
  ! input_output.f90:output_statistics (formats and record order unchanged;
  ! the harness read_stats parser depends on them):
  !   line 1 : '%', t, nu (2F15.8), nx_global, ny_global, nz_global, istep (4I8)
  !   line 2 : Cf(1:nx)      (nxF15.8)
  !   line 3 : x(1:nx)       (nxF15.8)
  !   line 4 : y(1:ny)       (nyF15.8)
  !   then, one x-station per line in (nyF15.8):
  !   Umean, Vmean, Wmean, U2mean, V2mean, W2mean, Pmean, P2mean (nx lines each)
  !---------------------------------------------------------------------------
  subroutine write_stats_file(istep)

    integer, intent(in) :: istep

    character(200) :: fname, my_format
    character(8)   :: ext, ext_ny, ext_nx
    integer        :: ii, unit_t

    ! create file name
    write(ext,'(I0.8)') istep + nstep_init

    write(ext_nx,'(I8)') nx
    write(ext_ny,'(I8)') ny

    fname = trim(adjustl(fileout))//'.'//trim(adjustl(ext))//'.stats.txt'
    write(*,*) 'writting ', trim(adjustl(fname))
    open(newunit=unit_t, file=fname, form='formatted', action='write')
    !
    write(unit_t,'(A,2F15.8,4I8)') '%', t, nu, nx_global, ny_global, nz_global, istep

    my_format = '('//trim(adjustl(ext_nx))//'F15.8)'
    write(unit_t,my_format) Cf

    my_format = '('//trim(adjustl(ext_nx))//'F15.8)'
    write(unit_t,my_format) x

    my_format = '('//trim(adjustl(ext_ny))//'F15.8)'
    write(unit_t,my_format) y

    !
    do ii = 1, nx
      write(unit_t,my_format) Umean(ii,:)
    end do
    do ii = 1, nx
      write(unit_t,my_format) Vmean(ii,:)
    end do
    do ii = 1, nx
      write(unit_t,my_format) Wmean(ii,:)
    end do
    !
    do ii = 1, nx
      write(unit_t,my_format) U2mean(ii,:)
    end do
    do ii = 1, nx
      write(unit_t,my_format) V2mean(ii,:)
    end do
    do ii = 1, nx
      write(unit_t,my_format) W2mean(ii,:)
    end do
    !
    do ii = 1, nx
      write(unit_t,my_format) Pmean(ii,:)
    end do
    do ii = 1, nx
      write(unit_t,my_format) P2mean(ii,:)
    end do
    !
    close(unit_t)

  end subroutine write_stats_file

  !---------------------------------------------------------------------------
  ! interp_x_avg — linear interpolation in x, faces -> centers (uniform mesh
  ! assumed). Host transcription of interpolation.f90:interpolate_x (which
  ! applies the same formula for either di; statistics calls it with di = 1).
  !   ui = 0 everywhere, then
  !   ui(1:n1-1,1:n2,1:n3) = 0.5*( u(1:n1-1,:,:) + u(2:n1,:,:) )
  ! with (n1,n2,n3) = shape(u); ui may be larger than u (term arrays).
  !---------------------------------------------------------------------------
  subroutine interp_x_avg(u, ui)

    real(dp), intent(in)  :: u (:,:,:)
    real(dp), intent(out) :: ui(:,:,:)

    integer :: n1, n2, n3

    n1 = size(u,1)
    n2 = size(u,2)
    n3 = size(u,3)

    ui = 0.0_dp
    ui(1:n1-1, 1:n2, 1:n3) = 0.5_dp*( u(1:n1-1,:,:) + u(2:n1,:,:) )

  end subroutine interp_x_avg

  !---------------------------------------------------------------------------
  ! interp_y_weighted — linear interpolation in y, centers -> faces with the
  ! stretched-grid weights weight_y_0/weight_y_1(ny) from grid_mod. Host
  ! transcription of interpolation.f90:interpolate_y, di = 2 branch:
  !   ui = 0 everywhere, then for each (i1,i3)
  !   ui(i1,1:n2-1,i3) = weight_y_0*u(i1,1:n2-1,i3) + weight_y_1*u(i1,2:n2,i3)
  ! with (n1,n2,n3) = shape(u); n2-1 = ny for all callers here.
  !---------------------------------------------------------------------------
  subroutine interp_y_weighted(u, ui)

    real(dp), intent(in)  :: u (:,:,:)
    real(dp), intent(out) :: ui(:,:,:)

    integer :: n1, n2, n3, i1, i3

    n1 = size(u,1)
    n2 = size(u,2)
    n3 = size(u,3)

    ui = 0.0_dp
    do i1 = 1, n1
      do i3 = 1, n3
        ui(i1, 1:n2-1, i3) = weight_y_0*u(i1,1:n2-1,i3) + weight_y_1*u(i1,2:n2,i3)
      end do
    end do

  end subroutine interp_y_weighted

  !===========================================================================
  ! (c) Monitor
  !===========================================================================

  !---------------------------------------------------------------------------
  ! output_monitor — print the per-interval monitor block. Port of
  ! monitor.f90:output_monitor at nprocs = 1 with LES_model = 0 and
  ! iwall_model = 0 (the only configuration the port accepts — the LES and
  ! wall-model print branches are unreachable and therefore not transcribed;
  ! the blank lines the reference emits on this path are preserved exactly).
  !
  ! Gate: mod(istep,nmonitor) == 0 .or. istep == 1 (as reference).
  !
  ! All format strings are VERBATIM from monitor.f90 (the validation harness
  ! regex-parses 'Mean Cf:     :', 'Maximum U    :',
  ! 'Maximum divergence          :', 'Elapsed time (s)            :').
  ! Interior means/maxima are host reductions over the refreshed mirrors
  ! (identical to the reference, which also reduces on the host after the
  ! gated device download). The divergence comes from poisson_mod's
  ! check_divergence, computed on the device. MPI_WTIME -> system_clock.
  !
  ! Note: 'Mean Cf' and the inflow Reynolds numbers are the values last
  ! computed by output_stats — call output_stats before output_monitor
  ! (reference main.f90 order); both trigger at istep == 1.
  !---------------------------------------------------------------------------
  subroutine output_monitor(istep)

    integer, intent(in) :: istep

    real(dp) :: maxU, maxV, maxW
    real(dp) :: meanU, meanV, meanW
    real(dp) :: max_divergence
    real(dp) :: time2

    if (.not. io_initialized) call io_init()

    if ( mod(istep,nmonitor) == 0 .or. istep == 1 ) then

      call refresh_host_fields(istep)

      ! compute mean values. The overlapping-slab bookkeeping makes the
      ! local interiors an exact partition of the global interior, so the
      ! cross-rank sums with the *_global denominators are exact.
      meanU = sum( U(2:nx-1, 2:nyg-1, 2:nzg-1) )
      meanV = sum( V(2:nxg-1, 2:ny-1, 2:nzg-1) )
      meanW = sum( W(2:nxg-1, 2:nyg-1, 2:nz-1) )
      call allreduce_sum(meanU)
      call allreduce_sum(meanV)
      call allreduce_sum(meanW)

      ! compute maximum values (global across ranks)
      maxU = maxval( U(2:nx-1, 2:nyg-1, 2:nzg-1) )
      maxV = maxval( V(2:nxg-1, 2:ny-1, 2:nzg-1) )
      maxW = maxval( W(2:nxg-1, 2:nyg-1, 2:nz-1) )
      call allreduce_max(maxU)
      call allreduce_max(maxV)
      call allreduce_max(maxW)

      call check_divergence(max_divergence)   ! already globally reduced

      ! end measure time per step
      time2 = wtime()

      if (.not. is_root()) then
         time1 = wtime()
         return
      end if

      write(*,*) 'step number  :', istep
      write(*,*) 'time         :', t
      write(*,*) 'time step    :', dt
      if (CFL < 0.0_dp) write(*,*) 'min time step:', dt_min_cfl

      write(*,*) ' '
      write(*,*) 'Mean Cf:     :', sum(Cf)/real(nx,dp)

      write(*,*) ' '
      write(*,*) 'Maximum U    :', maxU
      write(*,*) 'Maximum V    :', maxV
      write(*,*) 'Maximum W    :', maxW

      write(*,*) 'Mean U       :', meanU/real( nxm_global*nym_global*nzm_global, dp )
      write(*,*) 'Mean V       :', meanV/real( nxm_global*nym_global*nzm_global, dp )
      write(*,*) 'Mean W       :', meanW/real( nxm_global*nym_global*nzm_global, dp )

      ! (LES_model > 0 nu_t block unreachable in this port)

      write(*,*) ' '
      write(*,*) 'Mean pressure gradient in x :', dPdx
      write(*,*) 'Mean pressure gradient in y :', dPdy
      write(*,*) 'Spanwise rotation Omega_z   :', Omega_z

      write(*,*) ' '
      write(*,*) 'No wall model'   ! iwall_model == 0 always (check_supported)

      write(*,*) ' '
      ! (iwall_model > 0 utau and alpha blocks unreachable in this port)

      write(*,*) ' '
      write(*,*) 'Inflow Reynolds numbers:'
      write(*,*) '                Re_x        : ', Rex_inlet
      write(*,*) '                Re_theta    : ', Retheta_inlet
      write(*,*) '                Re_delta    : ', Redelta_inlet
      write(*,*) '                delta99_i   : ', delta99_inlet_ins

      write(*,*) ' '
      write(*,*) 'Maximum divergence          :', max_divergence
      write(*,*) 'Elapsed time (s)            :', time2-time1

      write(*,*) '------------------------------------------------------'

      ! start measure time per step
      time1 = wtime()

    end if

  end subroutine output_monitor

end module io_mod
