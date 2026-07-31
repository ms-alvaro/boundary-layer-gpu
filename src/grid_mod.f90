!------------------------------------------------------------------------------
! grid_mod — grid generation and metrics (host + device), single-GPU port
!
! Purpose:
!   Generates the staggered finite-difference grid and all metric arrays used
!   by the DNS kernels, and mirrors them to the GPU.
!
! Reference provenance (verbatim numerics):
!   * face-grid generation (x, y with sinh stretching and alpha parameter, z):
!       input_output.f90, subroutine init_flow      (lines ~148-174)
!   * center grids, ghost-extended grids, yg_m/yg_mm, min spacings, domain
!     lengths, dx/dz Fourier spacings, interpolation weights:
!       initialization.f90, subroutine initialize   (lines ~84-104, 174-227,
!                                                    253-259, 478-484)
!
! Single-GPU simplifications (contract: docs/CUDA_PORT_DESIGN.md):
!   * nprocs == 1, no MPI: the reference local/global split collapses
!     (x == x_global, nz == nz_global, ...). Only the LOCAL names are kept,
!     since those are the names the reference kernels use.
!   * Restart runs read the grid from file (io_mod); grid_generate() implements
!     the generated-grid path (random_init == 1) of the reference.
!
! Grid layout (staggered, 2nd-order FD):
!   x(nx),  y(ny),  z(nz)    face coordinates
!   xm(nxm), ym(nym), zm(nzm) cell-center coordinates (midpoints of faces)
!   xg(nxg), yg(nyg), zg(nzg) centers extended with one ghost point per side
!   yg_m(nyg-1)               midpoints of yg   (NOT equal to y in general)
!   yg_mm(nyg-2)              midpoints of yg_m (NOT equal to ym in general)
!   weight_y_0/1(ny)          linear interpolation weights, centers -> y-faces
!
! Device copies carry the same names with an `_d` suffix and are filled by
! grid_to_device(); kernels use-associate them from this module.
!------------------------------------------------------------------------------
module grid_mod

  use precision_mod, only: dp
  use mpi_mod,       only: myid, nprocs
  use param_mod,     only: nx_global, ny_global, nz_global, file_ygrid, &
                           Lx_rand, Ly_rand, Lz_rand, alpha_rand
  use cudafor

  implicit none
  private

  !----------------------------- grid sizes ----------------------------------
  ! face points / center points / centers + ghost cells
  ! (nprocs = 1: local sizes equal global sizes; initialization.f90:84-104)
  integer, public :: nx,  ny,  nz    ! face points
  integer, public :: nxm, nym, nzm   ! center points        (n*m = n* - 1)
  integer, public :: nxg, nyg, nzg   ! centers + 2 ghosts   (n*g = n*m + 2)

  !--------------------------- host grid arrays ------------------------------
  ! z-slab decomposition (P5.1; legacy initialization.f90:52-110 bookkeeping):
  ! rank r owns global face planes k1..k2 and center(+ghost) planes kg1..kg2;
  ! neighboring slabs overlap by two index planes by construction. nprocs = 1
  ! collapses to k1 = 1, k2 = nz_global (identity).
  integer, public :: nslices_z = 0
  integer, public :: k1 = 1, k2 = 0, kg1 = 1, kg2 = 0
  real(dp), allocatable, public :: z_global(:), zm_global(:)  ! for IO/Poisson

  real(dp), allocatable, public :: x(:),  y(:),  z(:)    ! faces
  real(dp), allocatable, public :: xm(:), ym(:), zm(:)   ! centers
  real(dp), allocatable, public :: xg(:), yg(:), zg(:)   ! centers + ghosts

  ! midpoints of yg and of yg_m (initialization.f90:214-217)
  real(dp), allocatable, public :: yg_m(:), yg_mm(:)

  ! interpolation weights centers -> y-faces (initialization.f90:483-484)
  real(dp), allocatable, public :: weight_y_0(:), weight_y_1(:)

  !--------------------------- scalar metrics --------------------------------
  ! total domain size (initialization.f90:225-227)
  real(dp), public :: Lx, Ly, Lz
  ! minimum grid spacings, used by the CFL condition (initialization.f90:220-222)
  real(dp), public :: dxmin, dymin, dzmin
  ! constant Fourier grid spacings for the pressure solver: dx = dxmin,
  ! dz = dzmin (initialization.f90:254-255)
  real(dp), public :: dx, dz

  !--------------------------- device copies ---------------------------------
  ! Same names + _d suffix (design-doc convention); filled by grid_to_device().
  real(dp), device, allocatable, public :: x_d(:),  y_d(:),  z_d(:)
  real(dp), device, allocatable, public :: xm_d(:), ym_d(:), zm_d(:)
  real(dp), device, allocatable, public :: xg_d(:), yg_d(:), zg_d(:)
  real(dp), device, allocatable, public :: yg_m_d(:), yg_mm_d(:)
  real(dp), device, allocatable, public :: weight_y_0_d(:), weight_y_1_d(:)

  public :: grid_generate, grid_to_device

contains

  !----------------------------------------------------------------------------
  ! grid_generate — host-side grid generation and metric computation.
  !
  ! Sets the size scalars, allocates all host grid arrays and computes them
  ! with numerics transcribed verbatim from the reference:
  !   1. face grids x, y (sinh stretching with alpha = alpha_rand), z
  !      [input_output.f90, init_flow]
  !   2. center grids xm, ym, zm; ghost-extended grids xg, yg, zg;
  !      yg_m, yg_mm; dxmin/dymin/dzmin; Lx/Ly/Lz; dx/dz;
  !      weight_y_0/weight_y_1 [initialization.f90, initialize]
  !
  ! Requires nx_global/ny_global/nz_global and Lx_rand/Ly_rand/Lz_rand/
  ! alpha_rand to have been set by param_mod input parsing.
  !----------------------------------------------------------------------------
  subroutine grid_generate()

    integer  :: i, j, k
    real(dp) :: alpha

    !------------------------- grid sizes (nprocs = 1) ------------------------
    ! face points (initialization.f90:84-86; nprocs=1 -> nz = nz_global)
    nx = nx_global
    ny = ny_global

    ! ---- z-slab decomposition (identity at nprocs = 1) ----
    nslices_z = (nz_global - 2)/nprocs
    if (nslices_z < 1) stop 'Error: nslices_z must be at least 1'
    k1  = myid*nslices_z + 1
    k2  = k1 + nslices_z + 1
    kg1 = myid*nslices_z + 1
    kg2 = kg1 + nslices_z + 1
    if (myid == nprocs-1) then
       k2  = nz_global
       kg2 = nz_global + 1
    end if
    nz = k2 - k1 + 1

    ! center points (initialization.f90:93-99). In z the counts come from
    ! the center-slab bookkeeping: nzm = kg2-kg1-1, nzg = kg2-kg1+1 — for
    ! interior ranks nzg = nz (NOT nz+1); the last rank (and nprocs = 1)
    ! carries one extra center-ghost plane (nzg = nz+1).
    nxm = nx - 1
    nym = ny - 1
    nzm = kg2 - kg1 + 1 - 2

    ! centers + ghost cells (initialization.f90:102-108)
    nxg = nxm + 2
    nyg = nym + 2
    nzg = kg2 - kg1 + 1

    !----------------------------- allocation ---------------------------------
    ! (initialization.f90:118-124, 482)
    allocate ( x (nx),    y (ny),    z (nz)    )
    allocate ( xm(nxm),   ym(nym),   zm(nzm)   )
    allocate ( xg(nxm+2), yg(nym+2), zg(nzm+2) )
    allocate ( yg_m (nyg-1) )
    allocate ( yg_mm(nyg-2) )
    allocate ( weight_y_0(ny), weight_y_1(ny) )

    !------------------- face grids (input_output.f90:148-174) ----------------
    ! xmesh: x(1) = x0 = 1 -> distance to the leading edge
    do i = 1, nx
       x(i) = real(i-1, dp)
    end do
    x = x - x(1)
    x = 1.0_dp + x/maxval(x)*Lx_rand

    ! zmesh (uniform; last face duplicates the periodic image). Built
    ! GLOBALLY, then sliced to this rank's slab (legacy generates global
    ! grids in init_flow and slices in initialization.f90).
    allocate ( z_global(nz_global), zm_global(nz_global-1) )
    do k = 1, nz_global
       z_global(k) = real(k-1, dp)
    end do
    z_global = z_global/z_global(nz_global-1)*Lz_rand
    do k = 1, nz_global-1
       zm_global(k) = 0.5_dp*( z_global(k) + z_global(k+1) )
    end do
    z = z_global(k1:k2)

    ! ymesh: either read from a file ('ygrid_file', e.g. a blended-sinh grid
    ! for HIT freestream resolution, inflow_flag=6 campaigns) or generated as
    ! a sinh-stretched coordinate rescaled to Ly (the reference default).
    ! (Port of the tbl-gpu init_flow branch, legacy/input_output.f90.)
    if ( len_trim(file_ygrid) > 0 ) then
       call read_ygrid_from_file()
    else
       ! Note: delta99/x0 = 4.91/sqrt(Rex0) = 0.0155 for Rex0 = 1e5
       do i = 1, ny
          y(i) = real(i-1, dp)
       end do
       y = y - y(1)
       y = y/maxval(y)
       alpha = alpha_rand
       do i = 1, ny
          y(i) = sinh(alpha*y(i))/sinh(alpha)
       end do
       y = y - y(1)
       y = y/maxval(y)*Ly_rand
    end if

    !------------------ center grids (initialization.f90:181-194) -------------
    do i = 1, nxm
       xm(i) = 0.5_dp*( x(i) + x(i+1) )
    end do
    do j = 1, nym
       ym(j) = 0.5_dp*( y(j) + y(j+1) )
    end do
    ! local z-centers are a SLICE of the global centers: local center m
    ! corresponds to global center kg1+m-1 (matches the field-plane mapping
    ! local plane p <-> global center plane kg1+p-1).
    do k = 1, nzm
       zm(k) = zm_global(kg1 + k - 1)
    end do

    !--------------- ghost-extended grids (initialization.f90:197-211) --------
    ! interior = centers; ghosts mirror the boundary faces at twice the
    ! face-to-first-center distance (keeps the wall face at the midpoint).
    xg(2:nxm+1) = xm
    xg(1)       = xm(1)   - 2.0_dp*( xm(1) - x(1)  )
    xg(nxm+2)   = xm(nxm) + 2.0_dp*( x(nx) - xm(nxm) )

    yg(2:nym+1) = ym
    yg(1)       = ym(1)   - 2.0_dp*( ym(1) - y(1)  )
    yg(nym+2)   = ym(nym) + 2.0_dp*( y(ny) - ym(nym) )

    ! zg: local center+ghost planes are the global center grid evaluated at
    ! global indices kg1-1 .. kg1+nzg-2 shifted by the ghost convention. At
    ! the global ends the reference ghost extrapolation applies; interior
    ! slab edges take the neighboring GLOBAL center directly (identical for
    ! the uniform z of this solver, and consistent across ranks).
    do k = 1, nzg
       block
         integer :: gc
         gc = kg1 + k - 2      ! global center index for local zg slot k
         if (gc >= 1 .and. gc <= nz_global-1) then
            zg(k) = zm_global(gc)
         else if (gc < 1) then
            zg(k) = zm_global(1) - 2.0_dp*( zm_global(1) - z_global(1) )
         else
            zg(k) = zm_global(nz_global-1) &
                  + 2.0_dp*( z_global(nz_global) - zm_global(nz_global-1) )
         end if
       end block
    end do

    !------------- yg midpoint hierarchies (initialization.f90:214-217) -------
    ! middle points for yg (.not. equal to y in general)
    yg_m = 0.5_dp*( yg(2:nyg) + yg(1:nyg-1) )

    ! middle points for yg_m (.not. equal to ym in general)
    yg_mm = 0.5_dp*( yg_m(2:nyg-1) + yg_m(1:nyg-2) )

    !------------- minimum grid sizes for CFL (initialization.f90:220-222) ----
    dxmin = minval ( xg(2:nxg) - xg(1:nxg-1) )
    dymin = minval ( yg(2:nyg) - yg(1:nyg-1) )
    dzmin = minval ( zg(2:nzg) - zg(1:nzg-1) )
    ! (uniform z: local and global minima coincide on every rank)

    !------------------ total domain size (initialization.f90:225-227) --------
    Lx = x(nx) - x(1)
    Ly = y(ny) - y(1)
    Lz = z_global(nz_global) - z_global(1)

    !------- Fourier constant grid spacings (initialization.f90:254-255) ------
    dx = dxmin
    dz = dzmin

    !------------- interpolation weights (initialization.f90:483-484) ---------
    ! linear weights for center -> y-face interpolation:
    !   q_face(j) = weight_y_0(j)*q_center(j) + weight_y_1(j)*q_center(j+1)
    weight_y_0 = ( yg(2:nyg) - y(1:ny) ) / ( yg(2:nyg) - yg(1:nyg-1) )
    weight_y_1 = 1.0_dp - weight_y_0

  end subroutine grid_generate

  !----------------------------------------------------------------------------
  ! grid_to_device — allocate the device grid arrays and copy the host grids
  ! to the GPU. Call once after grid_generate() (grids are time-invariant).
  ! Assignments of host arrays to device allocatables perform synchronous
  ! host-to-device copies (CUDA Fortran semantics).
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! read_ygrid_from_file — custom wall-normal face grid from an ASCII file:
  ! one value per line, blank lines and lines starting with '#' skipped
  ! (port of legacy/input_output.f90 read_ygrid_from_file, verbatim behavior).
  !----------------------------------------------------------------------------
  subroutine read_ygrid_from_file()

    integer            :: j, ios, unit_y
    character(len=256) :: line

    write(*,*) '   Reading y-grid from: ', trim(file_ygrid)

    open(newunit=unit_y, file=trim(file_ygrid), form='formatted', &
         action='read', status='old')
    j = 0
    do while ( j < ny )
       read(unit_y, '(a)', iostat=ios) line
       if ( ios /= 0 ) exit
       line = adjustl(line)
       if ( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
       j = j + 1
       read(line, *) y(j)
    end do
    close(unit_y)

    if ( j /= ny ) then
       write(*,*) 'ERROR: ygrid_file has ', j, ' values, need ny_global=', ny
       stop 'ygrid_file size mismatch'
    end if

    write(*,*) '   y-grid: [', y(1), ',', y(ny), ']'
    write(*,*) '   dy_wall = ', y(2)-y(1), ' dy_top = ', y(ny)-y(ny-1)

  end subroutine read_ygrid_from_file

  subroutine grid_to_device()

    if (.not. allocated(x_d)) then
       allocate ( x_d (nx),  y_d (ny),  z_d (nz)  )
       allocate ( xm_d(nxm), ym_d(nym), zm_d(nzm) )
       allocate ( xg_d(nxg), yg_d(nyg), zg_d(nzg) )
       allocate ( yg_m_d (nyg-1) )
       allocate ( yg_mm_d(nyg-2) )
       allocate ( weight_y_0_d(ny), weight_y_1_d(ny) )
    end if

    x_d  = x
    y_d  = y
    z_d  = z
    xm_d = xm
    ym_d = ym
    zm_d = zm
    xg_d = xg
    yg_d = yg
    zg_d = zg
    yg_m_d  = yg_m
    yg_mm_d = yg_mm
    weight_y_0_d = weight_y_0
    weight_y_1_d = weight_y_1

  end subroutine grid_to_device

end module grid_mod
