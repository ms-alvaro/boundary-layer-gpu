!=======================================================================
! wm_bridge.f90 — glue between the CUDA solver and the VENDORED
! reference wall-model stack (legacy_* files: modules global,
! mpi_legacy, interpolation, boundary_conditions, subgrid,
! Newton_solver/functions_wallmodel, wallmodel — copied verbatim from
! the pre-port CPU solver).
!
! The wall models are research-grade host algorithms (full-field test
! filters, Newton solves, wall-plane averages). Rather than re-deriving
! ~2000 lines as kernels, the reference code RUNS AS-IS on the host:
! once per RK substep (only when iwall_model > 0 or istress_model > 0)
! the bridge copies the device fields into the vendored `global`
! arrays, calls the reference compute_wall_model, and uploads the
! resulting alpha_x/y/z slip-length planes (+ V_bottom for WM=13, and
! nu_t, which several models CLOBBER with the filtered-field eddy
! viscosity — a reference quirk faithfully preserved) back to the
! device for the wall-BC kernels.
!
! Cost: a few full-field D2H/H2D copies + host compute per substep —
! slow relative to the 10 ms/step DNS path, but bit-faithful. Optimize
! per-model later if a model becomes production-critical.
!=======================================================================
module wm_bridge

  use precision_mod, only: dp
  use param_mod,     only: iwall_model_p => iwall_model,               &
                           istress_model_p => istress_model,           &
                           LES_model_p => LES_model,                   &
                           Dirichlet_nu_t_p => Dirichlet_nu_t,         &
                           iwall_model_nut_p => iwall_model_nut,       &
                           frac_vis_wall_model_p => frac_vis_wall_model, &
                           nu_p => nu, dPdx_p => dPdx, dPdz_p => dPdz, &
                           nx_global_p => nx_global,                   &
                           ny_global_p => ny_global,                   &
                           nz_global_p => nz_global,                   &
                           nxg_global_p => nxg_global,                 &
                           nyg_global_p => nyg_global,                 &
                           nzg_global_p => nzg_global,                 &
                           nxm_global_p => nxm_global,                 &
                           nym_global_p => nym_global,                 &
                           nzm_global_p => nzm_global,                 &
                           amx_p => alpha_mean_x, amy_p => alpha_mean_y, &
                           amz_p => alpha_mean_z, astd_p => alpha_std,  &
                           fmult_p => freq_mult
  use grid_mod,      only: nx_g => nx, ny_g => ny, nz_g => nz,         &
                           nxg_g => nxg, nyg_g => nyg, nzg_g => nzg,   &
                           nzm_g => nzm,                               &
                           x_g => x, y_g => y, z_g => z,               &
                           xg_g => xg, yg_g => yg, zg_g => zg,         &
                           w0_g => weight_y_0, w1_g => weight_y_1,     &
                           dx_g => dx, dz_g => dz,                     &
                           z_glob => z_global, zm_glob => zm_global
  use mpi_mod,       only: myid_p => myid, nprocs_p => nprocs
  use field_mod,     only: U_d, V_d, W_d, Uo_d, Vo_d, Wo_d
  use les_mod,       only: nu_t_d, les_avg => avg_nu_t
  use global   ! the vendored reference declaration module
  use mpi_legacy, only: myid, nprocs, nslices_z, k1_global, k2_global, &
                        kg1_global, kg2_global
  use wallmodel,  only: compute_wall_model
  use cudafor

  implicit none
  private

  public :: wm_active, wm_allocate, wm_update

  ! device copies of the wall-model outputs for the BC kernels
  real(dp), device, allocatable, public :: alpha_x_d(:,:,:), &
                                           alpha_y_d(:,:,:), &
                                           alpha_z_d(:,:,:), &
                                           V_bottom_d(:,:)

contains

  logical function wm_active()
    wm_active = ( iwall_model_p > 0 .or. istress_model_p > 0 )
  end function wm_active

  !--------------------------------------------------------------------
  ! wm_allocate — size the vendored `global` arrays exactly as the
  ! reference initialization.f90 does, and copy the static state
  ! (grids, weights, parameters). Call once after grid_generate.
  !--------------------------------------------------------------------
  subroutine wm_allocate()

    integer :: r

    if (.not. wm_active()) return

    ! ---- mpi_legacy state ----
    myid   = myid_p
    nprocs = nprocs_p
    nslices_z = (nz_global_p-2)/nprocs_p
    allocate( k1_global(0:nprocs_p-1),  k2_global(0:nprocs_p-1) )
    allocate( kg1_global(0:nprocs_p-1), kg2_global(0:nprocs_p-1) )
    do r = 0, nprocs_p-1
       ! reference initialization.f90 decomposition bookkeeping
       k1_global(r)  = r*nslices_z + 1
       k2_global(r)  = k1_global(r) + (nslices_z+2) - 1
       kg1_global(r) = r*nslices_z + 1
       kg2_global(r) = kg1_global(r) + (nslices_z+3) - 1
    end do

    ! ---- scalar/parameter state ----
    nx  = nx_g;  ny  = ny_g;  nz  = nz_g
    nxg = nxg_g; nyg = nyg_g; nzg = nzg_g
    nxm = nx_g-1; nym = ny_g-1; nzm = nzm_g
    nx_global  = nx_global_p;  ny_global  = ny_global_p
    nz_global  = nz_global_p
    nxg_global = nxg_global_p; nyg_global = nyg_global_p
    nzg_global = nzg_global_p
    nxm_global = nxm_global_p; nym_global = nym_global_p
    nzm_global = nzm_global_p
    nu   = nu_p
    dPdx = dPdx_p; dPdz = dPdz_p
    LES_model      = LES_model_p
    iwall_model    = iwall_model_p
    istress_model  = istress_model_p
    iwall_model_nut = iwall_model_nut_p
    frac_vis_wall_model = frac_vis_wall_model_p
    Dirichlet_nu_t = 0          ! reference hard-override
    penetration    = 1          ! reference initialization.f90:545
    alpha_mean_x = amx_p; alpha_mean_y = amy_p; alpha_mean_z = amz_p
    alpha_std    = astd_p; freq_mult = fmult_p
    int_len        = 0.0_dp     ! never assigned in the reference input path
    beta_y         = 0.0_dp

    ! ---- grids ----
    allocate( x(nx), y(ny), z(nz), xg(nxg), yg(nyg), zg(nzg) )
    x = x_g; y = y_g; z = z_g; xg = xg_g; yg = yg_g; zg = zg_g
    allocate( x_global(nx_global), y_global(ny_global), z_global(nz_global) )
    x_global = x_g; y_global = y_g; z_global = z_glob
    allocate( weight_y_0(ny), weight_y_1(ny) )
    weight_y_0 = w0_g; weight_y_1 = w1_g
    dx = dx_g; dz = dz_g

    ! ---- fields + scratch (reference initialization sizes) ----
    allocate( U(nx,nyg,nzg), V(nxg,ny,nzg), W(nxg,nyg,nz) )
    allocate( Uo(nx,nyg,nzg), Vo(nxg,ny,nzg), Wo(nxg,nyg,nz) )
    U = 0.0_dp; V = 0.0_dp; W = 0.0_dp
    Uo = 0.0_dp; Vo = 0.0_dp; Wo = 0.0_dp
    allocate( term  (nxg,nyg,nzm+2), term_1(nxg,nyg,nzm+2) )
    allocate( term_2(nxg,nyg,nzm+2), term_3(nxg,nyg,nzm+2) )
    allocate( term_4(nxg,nyg,nzm+2) )
    term = 0.0_dp; term_1 = 0.0_dp; term_2 = 0.0_dp
    term_3 = 0.0_dp; term_4 = 0.0_dp
    allocate( nu_t(nxg,nyg,nzg), avg_nu_t(nxg,nyg,1), avg_nu_t_hat(nxg,nyg,1) )
    nu_t = 0.0_dp; avg_nu_t = 0.0_dp; avg_nu_t_hat = 0.0_dp

    ! LES scratch used by the vendored subgrid module (reference sizes)
    allocate( Lij(2:nxg,2:nyg-1,2:nzm+2,6), Mij(2:nxg-1,2:nyg-1,2:nzm+1,6) )
    allocate( Sij(2:nxg,2:nyg,2:nzm+2,6),   S(2:nxg-1,2:nyg-1,2:nzm+1) )
    allocate( ten_buf(1:nxg,1:nyg,1:nzm+2,6) )
    allocate( Uf(1:nx,1:nyg,1:nzm+2),  Vf(1:nxg,1:ny,1:nzm+2),  Wf(1:nxg,1:nyg,1:nzm+2) )
    allocate( Uff(1:nx,1:nyg,1:nzm+2), Vff(1:nxg,1:ny,1:nzm+2), Wff(1:nxg,1:nyg,1:nzm+2) )
    Lij = 0.0_dp; Mij = 0.0_dp; Sij = 0.0_dp; S = 0.0_dp; ten_buf = 0.0_dp
    Uf = 0.0_dp; Vf = 0.0_dp; Wf = 0.0_dp
    Uff = 0.0_dp; Vff = 0.0_dp; Wff = 0.0_dp

    ! wall-model state
    allocate( alpha_x(1:nx,1:2,1:nzg), alpha_y(1:nxg,1:2,1:nzg), &
              alpha_z(1:nxg,1:2,1:nz) )
    allocate( alpha_xo(1:nx,1:2,1:nzg), alpha_yo(1:nxg,1:2,1:nzg), &
              alpha_zo(1:nxg,1:2,1:nz) )
    alpha_x = 0.0_dp; alpha_y = 0.0_dp; alpha_z = 0.0_dp
    alpha_xo = 0.0_dp; alpha_yo = 0.0_dp; alpha_zo = 0.0_dp
    allocate( V_bottom(nxm+2, nzm+2) ); V_bottom = 0.0_dp
    allocate( utau_model(nx_global), utau_wall(nx_global), &
              utau_wall_T(nx_global), utau_ref(nx_global) )
    utau_ref    = ( 0.5_dp*0.027_dp*(x_global/nu)**(-1.0_dp/7.0_dp) )**0.5_dp
    utau_model  = utau_ref
    utau_wall   = 0.0_dp
    utau_wall_T = 0.0_dp

    ! communication buffers (reference initialization.f90:228-238)
    allocate( buffer_ui(nx,nyg,2:3), buffer_vi(nxg,ny,2:3), &
              buffer_wi(nxg,nyg),    buffer_ci(nxg,nyg,2:3) )
    allocate( buffer_ue(nx,nyg), buffer_ve(nxg,ny), &
              buffer_we(nxg,nyg), buffer_ce(nxg,nyg) )
    allocate( buffer_us(nx,nyg),  buffer_ur(nx,nyg) )
    allocate( buffer_vs(nxg,ny),  buffer_vr(nxg,ny) )
    allocate( buffer_ws(nxg,nyg), buffer_wr(nxg,nyg) )

    ! device outputs for the wall-BC kernels
    allocate( alpha_x_d(1:nx,1:2,1:nzg), alpha_y_d(1:nxg,1:2,1:nzg), &
              alpha_z_d(1:nxg,1:2,1:nz), V_bottom_d(nxm+2, nzm+2) )
    alpha_x_d = 0.0_dp; alpha_y_d = 0.0_dp; alpha_z_d = 0.0_dp
    V_bottom_d = 0.0_dp

  end subroutine wm_allocate

  !--------------------------------------------------------------------
  ! wm_update — run the reference wall model on the host for the
  ! current substep state. first_substep: the step's base state is
  ! still in U_d/V_d/W_d (the fused Uo save has not happened yet).
  !--------------------------------------------------------------------
  subroutine wm_update(dt_in, t_in, first_substep)

    real(dp), intent(in) :: dt_in, t_in
    logical,  intent(in) :: first_substep

    if (.not. wm_active()) return

    dt = dt_in
    t  = t_in

    ! device -> vendored host state
    U = U_d; V = V_d; W = W_d
    if (first_substep) then
       Uo = U; Vo = V; Wo = W
    else
       Uo = Uo_d; Vo = Vo_d; Wo = Wo_d
    end if
    nu_t = nu_t_d_or_zero()
    if (allocated(les_avg)) avg_nu_t(:,:,1) = les_avg

    call compute_wall_model(U, V, W)

    ! vendored host outputs -> device
    alpha_x_d  = alpha_x
    alpha_y_d  = alpha_y
    alpha_z_d  = alpha_z
    V_bottom_d = V_bottom
    ! several models clobber nu_t with the filtered-field eddy
    ! viscosity (reference quirk) — push the host copy back
    if (LES_model_p > 0) nu_t_d = nu_t

  end subroutine wm_update

  function nu_t_d_or_zero() result(a)
    real(dp) :: a(nxg,nyg,nzg)
    if (LES_model_p > 0) then
       a = nu_t_d
    else
       a = 0.0_dp
    end if
  end function nu_t_d_or_zero

end module wm_bridge
