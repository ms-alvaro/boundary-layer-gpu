!------------------------------------------------------------------------------
! field_mod — state arrays of the DNS timestep (device) + host mirrors
!
! Purpose:
!   Owns every array the RK2/RK3 timestep touches, as module-level
!   `device, allocatable` variables with the reference names suffixed `_d`
!   (contract: docs/CUDA_PORT_DESIGN.md). Host mirrors exist only for
!   U, V, W, P — needed for I/O (snapshots, restart) and setup (IC, inflow).
!
! Reference provenance (shapes transcribed from allocation code):
!   * global.f90            — array declarations
!   * initialization.f90    — allocation shapes:
!       U/V/W/P and Uo/Vo/Wo                 lines ~126-134
!       rhs_uo/rhs_vo/rhs_wo/rhs_p           lines ~151-154
!       RK stage arrays Fu1..Fw3             lines ~532-554
!
! Single-GPU simplifications (nprocs = 1, no MPI):
!   * nym+2 = nyg, nzm+2 = nzg, so the reference shapes read
!       U(nx,nyg,nzg), V(nxg,ny,nzg), W(nxg,nyg,nz), P(nxg,nyg,nzg).
!   * Uoo/Voo/Woo (rank-local I/O staging), MPI buffers, and Po (only used
!     by the dead true-pressure solve) are NOT ported.
!   * Inflow/top-BC tables live in ic_inflow_mod; Poisson work arrays
!     (plane, rhs_hat, Thomas factors, twiddles) live in poisson_mod.
!
! Array roles in the timestep (time_integration.f90):
!   U_d/V_d/W_d/P_d      current velocity components and pressure
!   Uo_d/Vo_d/Wo_d       velocity at the beginning of the step (RK base state,
!                        also used by the convective-outflow BC)
!   rhs_uo/vo/wo_d       Navier-Stokes RHS of the Euler path in the reference;
!                        kept for RHS staging (equations.f90 output shape)
!   rhs_p_d              divergence of the intermediate velocity (Poisson RHS);
!                        ONE EXTRA PLANE IN Z FOR GHOST CELL (reference note)
!   Fu1..Fu3_d (v,w)     RK substep RHS storage; the *3 arrays are only
!                        needed by RK3 (itime_step == 3) and are allocated
!                        (1,1,1) otherwise, exactly as the reference does.
!
! Index conventions (staggered grid): RHS and stage arrays are allocated on
! INTERIOR points only, with the reference's explicit lower bounds of 2
! (e.g. rhs_uo_d(2:nx-1, 2:nyg-1, 2:nzg-1)); kernels must index accordingly.
!------------------------------------------------------------------------------
module field_mod

  use precision_mod, only: dp
  use param_mod,     only: itime_step
  use grid_mod,      only: nx, ny, nz, nxg, nyg, nzg
  use cudafor

  implicit none
  private

  !--------------------- host mirrors (I/O and setup only) --------------------
  real(dp), allocatable, public :: U(:,:,:)   ! (nx , nyg, nzg) u at x-faces
  real(dp), allocatable, public :: V(:,:,:)   ! (nxg, ny , nzg) v at y-faces
  real(dp), allocatable, public :: W(:,:,:)   ! (nxg, nyg, nz ) w at z-faces
  real(dp), allocatable, public :: P(:,:,:)   ! (nxg, nyg, nzg) p at centers

  !------------------------- device state arrays ------------------------------
  ! velocities and pressure
  real(dp), device, allocatable, public :: U_d(:,:,:), V_d(:,:,:), W_d(:,:,:)
  real(dp), device, allocatable, public :: P_d(:,:,:)

  ! beginning-of-step velocities (RK base state / convective outflow BC)
  real(dp), device, allocatable, public :: Uo_d(:,:,:), Vo_d(:,:,:), Wo_d(:,:,:)

  ! Navier-Stokes RHS staging (interior points only)
  real(dp), device, allocatable, public :: rhs_uo_d(:,:,:)
  real(dp), device, allocatable, public :: rhs_vo_d(:,:,:)
  real(dp), device, allocatable, public :: rhs_wo_d(:,:,:)

  ! Poisson RHS: div(u*) at centers, one extra ghost plane in z
  real(dp), device, allocatable, public :: rhs_p_d(:,:,:)

  ! Runge-Kutta stage arrays (interior points only)
  real(dp), device, allocatable, public :: Fu1_d(:,:,:), Fu2_d(:,:,:), Fu3_d(:,:,:)
  real(dp), device, allocatable, public :: Fv1_d(:,:,:), Fv2_d(:,:,:), Fv3_d(:,:,:)
  real(dp), device, allocatable, public :: Fw1_d(:,:,:), Fw2_d(:,:,:), Fw3_d(:,:,:)

  public :: fields_allocate, fields_to_device, fields_from_device

contains

  !----------------------------------------------------------------------------
  ! fields_allocate — allocate host mirrors and all device state arrays.
  !
  ! Shapes transcribed from initialization.f90 (nprocs = 1: nym+2 = nyg,
  ! nzm+2 = nzg, nze = nz). Requires grid_generate() to have set the size
  ! scalars in grid_mod, and itime_step to have been read by param_mod.
  ! All arrays are zero-initialized for deterministic contents (the IC /
  ! restart reader then fills U, V, W; P is rebuilt by the first projection).
  !----------------------------------------------------------------------------
  subroutine fields_allocate()

    !----------------- host mirrors (initialization.f90:126-129) --------------
    allocate ( U(nx , nyg, nzg) )
    allocate ( V(nxg, ny , nzg) )
    allocate ( W(nxg, nyg, nz ) )
    allocate ( P(nxg, nyg, nzg) )
    U = 0.0_dp
    V = 0.0_dp
    W = 0.0_dp
    P = 0.0_dp

    !----------------- device state: velocities and pressure ------------------
    allocate ( U_d(nx , nyg, nzg) )
    allocate ( V_d(nxg, ny , nzg) )
    allocate ( W_d(nxg, nyg, nz ) )
    allocate ( P_d(nxg, nyg, nzg) )
    U_d = 0.0_dp
    V_d = 0.0_dp
    W_d = 0.0_dp
    P_d = 0.0_dp

    !----------------- beginning-of-step velocities (131-133) -----------------
    allocate ( Uo_d(nx , nyg, nzg) )
    allocate ( Vo_d(nxg, ny , nzg) )
    allocate ( Wo_d(nxg, nyg, nz ) )
    Uo_d = 0.0_dp
    Vo_d = 0.0_dp
    Wo_d = 0.0_dp

    !----------------- RHS: interior points only (151-154) --------------------
    allocate ( rhs_uo_d( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    allocate ( rhs_vo_d( 2:nxg-1, 2:ny-1,  2:nzg-1 ) )
    allocate ( rhs_wo_d( 2:nxg-1, 2:nyg-1, 2:nz-1  ) )
    allocate ( rhs_p_d ( 2:nxg-1, 2:nyg-1, 2:nzg   ) ) ! ONE EXTRA PLANE IN Z FOR GHOST CELL
    rhs_uo_d = 0.0_dp
    rhs_vo_d = 0.0_dp
    rhs_wo_d = 0.0_dp
    rhs_p_d  = 0.0_dp

    !----------------- Runge-Kutta stage arrays (532-554) ---------------------
    allocate ( Fu1_d( 2:nx-1, 2:nyg-1, 2:nzg-1 ) )
    allocate ( Fu2_d( 2:nx-1, 2:nyg-1, 2:nzg-1 ) )
    if (itime_step == 3) then
       allocate ( Fu3_d( 2:nx-1, 2:nyg-1, 2:nzg-1 ) )
    else
       allocate ( Fu3_d( 1:1, 1:1, 1:1 ) )
    end if

    allocate ( Fv1_d( 2:nxg-1, 2:ny-1, 2:nzg-1 ) )
    allocate ( Fv2_d( 2:nxg-1, 2:ny-1, 2:nzg-1 ) )
    if (itime_step == 3) then
       allocate ( Fv3_d( 2:nxg-1, 2:ny-1, 2:nzg-1 ) )
    else
       allocate ( Fv3_d( 1:1, 1:1, 1:1 ) )
    end if

    allocate ( Fw1_d( 2:nxg-1, 2:nyg-1, 2:nz-1 ) )
    allocate ( Fw2_d( 2:nxg-1, 2:nyg-1, 2:nz-1 ) )
    if (itime_step == 3) then
       allocate ( Fw3_d( 2:nxg-1, 2:nyg-1, 2:nz-1 ) )
    else
       allocate ( Fw3_d( 1:1, 1:1, 1:1 ) )
    end if

    Fu1_d = 0.0_dp
    Fu2_d = 0.0_dp
    Fu3_d = 0.0_dp
    Fv1_d = 0.0_dp
    Fv2_d = 0.0_dp
    Fv3_d = 0.0_dp
    Fw1_d = 0.0_dp
    Fw2_d = 0.0_dp
    Fw3_d = 0.0_dp

  end subroutine fields_allocate

  !----------------------------------------------------------------------------
  ! fields_to_device — copy the host mirrors U, V, W, P to the device.
  ! Used after setup (IC / restart read) and whenever the host modifies the
  ! state. Synchronous host-to-device copies (CUDA Fortran assignment).
  !----------------------------------------------------------------------------
  subroutine fields_to_device()

    U_d = U
    V_d = V
    W_d = W
    P_d = P

  end subroutine fields_to_device

  !----------------------------------------------------------------------------
  ! fields_from_device — copy the device state U, V, W, P back to the host
  ! mirrors. Used before I/O (snapshots, statistics, restart writes).
  ! Synchronous device-to-host copies (CUDA Fortran assignment).
  !----------------------------------------------------------------------------
  subroutine fields_from_device()

    U = U_d
    V = V_d
    W = W_d
    P = P_d

  end subroutine fields_from_device

end module field_mod
