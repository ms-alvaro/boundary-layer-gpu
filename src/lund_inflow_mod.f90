!=======================================================================
! lund_inflow_mod.f90 — Lund–Wu–Squires recycling inflow for the CUDA
! Fortran port (inflow_flag = 3 full rescaling, = 5 fluctuations only).
!
! Reference provenance
! --------------------
! Port of Module Lund_rescaled_bc (rescaled_inlet_bc.f90, identical in
! the pre-port master and the gpu-openacc production branch) plus the
! EMA-state initialization block of initialization.f90:556-591 and the
! per-step dispatch of boundary_conditions.f90:41-50, 371-456.
! Loop bodies are transcribed VERBATIM (same expressions, same operation
! order; d0 literals -> _dp) per docs/CUDA_PORT_DESIGN.md.
!
! Execution model
! ---------------
! The recycling planes are tiny (nyg x nzg), so the whole computation
! runs on the HOST, exactly as the reference wrote it: once per step
! (reference guard: step_beginning == 1) the driver calls
! lund_update_inlet(dt), which
!   1. copies the i = i_rescale and i = 1 y-z planes of U_d/V_d/W_d to
!      host buffers (strided D2H, ~100 KB total),
!   2. runs compute_rescaled_inflow verbatim,
!   3. uploads the result into Ut_inlet_d/Vt_inlet_d/Wt_inlet_d, which
!      bc_kernels' bc_inflow_kernel already applies every substep as
!      U(1,j,k) = U_inlet(j) + Ut_inlet(j,k).
! The call site is the top of the RK step, where U_d still holds the
! step's base state — the same field the reference passes as Uo/Vo/Wo.
!
! Flag semantics at the inlet (must match bc_inflow_kernel's addition):
!   * flag 3 (apply_inflow_bc_x_rescaling): the reference OVERWRITES
!     U(1,j,:) = Ut_inlet(j,:) — so ic_inflow_mod zeroes U_inlet for
!     flag 3 and the full rescaled plane is carried by Ut_inlet.
!   * flag 5 (apply_inflow_bc_x_Turbulent_rescaled_flu): the reference
!     adds fluctuation planes to the file-read turbulent mean profile
!     U_inlet(j), with the spanwise index of W SWAPPED:
!       W(1,j,k) = W_inlet(j) + Wt_inlet(j, mod(k+16,nz)+1)
!     (U and V compute the swapped index but use k — reference quirk,
!     preserved). The swap is applied host-side before upload.
!
! Multi-rank: the spanwise means are local sums over the rank's owned
! z-centers reduced with allreduce_sum_arr, divided by the GLOBAL
! center count (nzg_global-2) — the reference's exact MPI pattern. The
! flag-5 spanwise swap uses the LOCAL k and nz, as the reference does.
!
! Restart continuity: Umean_resc_To/Vmean_resc_To are filled by io_mod
! from the '.mean.rescaling' companion snapshot (inflow_flag = 3, as in
! the production branch); when absent (or flag 5) the EMA state starts
! from the IC's spanwise means (initialization.f90:567-583). The
! reference's companion branch reads Umean_inlet_T from uninitialized
! locals; the port always initializes it from the IC means instead.
!=======================================================================
module lund_inflow_mod

  use precision_mod, only: dp
  use grid_mod,      only: ny, nyg, nz, nzg, y, yg
  use param_mod,     only: i_rescale, delta_inlet, T_resc, nu, nzg_global, &
                           inflow_boundary_flag
  use mpi_mod,       only: allreduce_sum_arr
  use field_mod,     only: U_d, V_d, W_d
  use ic_inflow_mod, only: Ut_inlet_d, Vt_inlet_d, Wt_inlet_d

  implicit none
  private

  public :: lund_allocate       ! allocate state + plane buffers
  public :: lund_init_means     ! EMA state from companion or IC means
  public :: lund_update_inlet   ! once per step: planes -> Ut/Vt/Wt_inlet_d

  ! EMA state (reference global.f90: Umean_resc_T etc., shape (n,1))
  real(dp), allocatable, public :: Umean_resc_T(:,:),  Vmean_resc_T(:,:)
  real(dp), allocatable, public :: Umean_inlet_T(:,:), Vmean_inlet_T(:,:)
  ! companion-snapshot carriers (read by io_mod before lund_init_means)
  real(dp), allocatable, public :: Umean_resc_To(:,:), Vmean_resc_To(:,:)

  ! host plane buffers (i = i_rescale and i = 1 y-z planes)
  real(dp), allocatable :: U_resc(:,:), V_resc(:,:), W_resc(:,:)
  real(dp), allocatable :: U_inl(:,:),  V_inl(:,:),  W_inl(:,:)
  real(dp), allocatable :: Wt_swap(:,:)

contains

  !--------------------------------------------------------------------
  ! lund_allocate — reference initialization.f90:159-163, 558-561.
  !--------------------------------------------------------------------
  subroutine lund_allocate()

    allocate( Umean_resc_T (nyg,1), Vmean_resc_T (ny,1) )
    allocate( Umean_inlet_T(nyg,1), Vmean_inlet_T(ny,1) )
    allocate( Umean_resc_To(nyg,1), Vmean_resc_To(ny,1) )
    Umean_resc_To = 0.0_dp
    Vmean_resc_To = 0.0_dp

    allocate( U_resc(nyg,nzg), V_resc(ny,nzg), W_resc(nyg,nz) )
    allocate( U_inl (nyg,nzg), V_inl (ny,nzg), W_inl (nyg,nz) )
    allocate( Wt_swap(nyg,nz) )

  end subroutine lund_allocate

  !--------------------------------------------------------------------
  ! lund_init_means — reference initialization.f90:563-591. Called once
  ! after the initial condition (or restart) is resident in U_d/V_d.
  !--------------------------------------------------------------------
  subroutine lund_init_means()

    integer :: j

    ! pull the two planes of the freshly initialized field
    U_resc = U_d(i_rescale,:,:)
    V_resc = V_d(i_rescale,:,:)
    U_inl  = U_d(1,:,:)
    V_inl  = V_d(1,:,:)

    if ( sum(Umean_resc_To) > 1e-4_dp ) then
       ! use the companion-snapshot means (reference lines 563-566)
       Umean_resc_T = Umean_resc_To
       Vmean_resc_T = Vmean_resc_To
    else
       ! initialize for the first time (reference lines 567-583)
       do j = 1, nyg
          Umean_resc_T(j,1) = sum(U_resc(j,2:nzg-1))
       end do
       do j = 1, ny
          Vmean_resc_T(j,1) = sum(V_resc(j,2:nzg-1))
       end do
       call allreduce_sum_arr(Umean_resc_T(:,1), nyg)
       call allreduce_sum_arr(Vmean_resc_T(:,1), ny)
       Umean_resc_T = Umean_resc_T/real(nzg_global-2, dp)
       Vmean_resc_T = Vmean_resc_T/real(nzg_global-2, dp)
    end if

    ! inlet means always start from the IC (reference lines 586-591)
    do j = 1, nyg
       Umean_inlet_T(j,1) = sum(U_inl(j,2:nzg-1))
    end do
    do j = 1, ny
       Vmean_inlet_T(j,1) = sum(V_inl(j,2:nzg-1))
    end do
    call allreduce_sum_arr(Umean_inlet_T(:,1), nyg)
    call allreduce_sum_arr(Vmean_inlet_T(:,1), ny)
    Umean_inlet_T = Umean_inlet_T/real(nzg_global-2, dp)
    Vmean_inlet_T = Vmean_inlet_T/real(nzg_global-2, dp)

  end subroutine lund_init_means

  !--------------------------------------------------------------------
  ! lund_update_inlet — once per step, at the top of the RK step (U_d
  ! still holds the step's base state = the reference's Uo/Vo/Wo).
  ! Fills Ut/Vt/Wt_inlet_d for bc_inflow_kernel.
  !--------------------------------------------------------------------
  subroutine lund_update_inlet(dt)

    real(dp), intent(in) :: dt

    real(dp) :: Ut(nyg,nzg), Vt(ny,nzg), Wt(nyg,nz)
    integer  :: j, k, ks

    ! device planes -> host (reference passes full Uo/Vo/Wo; only these
    ! two planes are ever read)
    U_resc = U_d(i_rescale,:,:)
    V_resc = V_d(i_rescale,:,:)
    W_resc = W_d(i_rescale,:,:)
    U_inl  = U_d(1,:,:)
    V_inl  = V_d(1,:,:)
    W_inl  = W_d(1,:,:)

    if ( inflow_boundary_flag == 3 ) then
       ! full rescaling (reference apply_inflow_bc_x_rescaling)
       call compute_rescaled_inflow(Ut, Vt, Wt, dt, 0)
       Wt_swap = Wt
    else
       ! fluctuations only (apply_inflow_bc_x_Turbulent_rescaled_flu):
       ! means zeroed inside (iflag_fluc = 1); spanwise swap for W only
       call compute_rescaled_inflow(Ut, Vt, Wt, dt, 1)
       do k = 1, nz
          ks = mod( k + 16, nz ) + 1
          Wt_swap(:,k) = Wt(:,ks)
       end do
    end if

    ! upload; bc_inflow_kernel adds U_inlet(j) (zero for flag 3) every
    ! substep, so the plane persists across the step like the reference
    Ut_inlet_d = Ut
    Vt_inlet_d = Vt
    Wt_inlet_d = Wt_swap

  end subroutine lund_update_inlet

  !--------------------------------------------------------------------
  ! compute_rescaled_inflow — VERBATIM port of rescaled_inlet_bc.f90's
  ! compute_rescaled_inflow, operating on the module's extracted planes
  ! (U_resc/U_inl replace U(i_rescale,:,:)/U(1,:,:) etc.). Outputs the
  ! inlet planes into Ut/Vt/Wt (the reference's U_inlet/V_inlet/W_inlet
  ! outputs, which the caller stores as Ut_inlet/Vt_inlet/Wt_inlet).
  !--------------------------------------------------------------------
  subroutine compute_rescaled_inflow(Ut, Vt, Wt, dt, iflag_fluc)

    real(dp), intent(out) :: Ut(nyg,nzg), Vt(ny,nzg), Wt(nyg,nz)
    real(dp), intent(in)  :: dt
    integer,  intent(in)  :: iflag_fluc

    ! local variables (reference declarations, shapes preserved)
    integer  :: j, jref, k
    real(dp) :: utau_inlet, utau_resc, delta_resc, w1, w0
    real(dp) :: alpha, b, gamma, UVmean_wall, dUmean_wall, UVwall1(1)
    real(dp) :: Uinf, Umean99, theta_resc, theta_inlet, delta_inlet2
        ! rescaled and inlet means
    real(dp) :: Umean_resc(nyg,1), Umean_inlet(nyg,1), temp(nyg,1)
    real(dp) :: Vmean_resc(ny ,1), Vmean_inlet(ny ,1)
        ! fluctuations
    real(dp) :: Ufluc_resc(nyg,nzg), Ufluc_inlet(nyg,nzg)
    real(dp) :: Vfluc_resc(ny ,nzg), Vfluc_inlet(ny ,nzg)
    real(dp) :: Wfluc_resc(nyg,nz ), Wfluc_inlet(nyg,nz )
        ! inner and outer
    real(dp) :: Umean_inner(nyg,1), Umean_outer(nyg,1)
    real(dp) :: Vmean_inner(ny ,1), Vmean_outer(ny ,1)
    real(dp) :: Ufluc_inner(nyg,nzg), Ufluc_outer(nyg,nzg)
    real(dp) :: Vfluc_inner(ny ,nzg), Vfluc_outer(ny ,nzg)
    real(dp) :: Wfluc_inner(nyg,nz ), Wfluc_outer(nyg,nz )
        ! y meshes
    real(dp) :: ygp_inlet(nyg), ygp_resc(nyg), etag_inlet(nyg), etag_resc(nyg), Weig(nyg)
    real(dp) :: yp_inlet (ny ), yp_resc (ny ), eta_inlet (ny ), eta_resc (ny ), Wei (ny )

    !---------------------------------------------------------------------!
    ! compute means (reference lines 108-126; local z sum + allreduce)
    do j = 1, nyg
       Umean_resc (j,1) = sum(U_resc(j,2:nzg-1))
       Umean_inlet(j,1) = sum(U_inl (j,2:nzg-1))
    end do
    do j = 1, ny
       Vmean_resc (j,1) = sum(V_resc(j,2:nzg-1))
       Vmean_inlet(j,1) = sum(V_inl (j,2:nzg-1))
    end do

    call allreduce_sum_arr(Umean_resc (:,1), nyg)
    call allreduce_sum_arr(Vmean_resc (:,1), ny)
    call allreduce_sum_arr(Umean_inlet(:,1), nyg)
    call allreduce_sum_arr(Vmean_inlet(:,1), ny)

    Umean_resc  = Umean_resc /real(nzg_global-2, dp)
    Vmean_resc  = Vmean_resc /real(nzg_global-2, dp)
    Umean_inlet = Umean_inlet/real(nzg_global-2, dp)
    Vmean_inlet = Vmean_inlet/real(nzg_global-2, dp)

    Uinf = Umean_resc(nyg,1)

    ! time average (reference lines 130-137)
    Umean_resc_T  = dt/T_resc*Umean_resc  + (1.0_dp-dt/T_resc)*Umean_resc_T
    Vmean_resc_T  = dt/T_resc*Vmean_resc  + (1.0_dp-dt/T_resc)*Vmean_resc_T

    Umean_inlet_T = dt/T_resc*Umean_inlet + (1.0_dp-dt/T_resc)*Umean_inlet_T
    Vmean_inlet_T = dt/T_resc*Vmean_inlet + (1.0_dp-dt/T_resc)*Vmean_inlet_T

    !---------------------------------------------------------------------!
    ! fluctuating part (reference lines 140-150; note Vfluc_inlet uses the
    ! INSTANTANEOUS Vmean_inlet — reference quirk, preserved)
    do j = 1, nyg
       Ufluc_resc (j,:) = U_resc(j,:) - Umean_resc_T (j,1)
       Ufluc_inlet(j,:) = U_inl (j,:) - Umean_inlet_T(j,1)
       Wfluc_resc (j,:) = W_resc(j,:)
       Wfluc_inlet(j,:) = W_inl (j,:)
    end do
    do j = 1, ny
       Vfluc_resc (j,:) = V_resc(j,:) - Vmean_resc_T(j,1)
       Vfluc_inlet(j,:) = V_inl (j,:) - Vmean_inlet (j,1)
    end do

    !---------------------------------------------------------------------!
    ! compute delta_99 at i_rescale (reference lines 153-172)
    Umean99 = 0.99_dp*Uinf
    jref    = 0
    do j = 1, nyg
       if ( Umean_resc_T(j,1) >= Umean99 ) then
          jref = j
          exit
       end if
    end do
    if ( jref < 2 ) jref = nyg/2
    w1         = ( Umean99 - Umean_resc_T(jref-1,1) )/( Umean_resc_T(jref,1) - Umean_resc_T(jref-1,1) )
    w0         = 1.0_dp - w1
    delta_resc = w1*yg(jref) + w0*yg(jref-1)
    if ( delta_resc < 0.0_dp ) delta_resc = abs(delta_resc)

    ! compute delta_99 at inlet (reference lines 174-192)
    jref = 0
    do j = 1, nyg
       if ( Umean_inlet_T(j,1) >= Umean99 ) then
          jref = j
          exit
       end if
    end do
    if ( jref < 2 ) jref = nyg/2
    w1           = ( Umean99 - Umean_inlet_T(jref-1,1) )/( Umean_inlet_T(jref,1) - Umean_inlet_T(jref-1,1) )
    w0           = 1.0_dp - w1
    delta_inlet2 = w1*yg(jref) + w0*yg(jref-1)
    if ( delta_inlet2 < 0.0_dp ) delta_inlet2 = abs(delta_inlet2)

    !---------------------------------------------------------------------!
    ! momentum thickness at i_rescale and inlet (reference lines 195-217)
    theta_resc = 0.0_dp
    temp       = Umean_resc_T/Uinf*(1.0_dp - Umean_resc_T/Uinf)
    do j = 2, nyg
       theta_resc = theta_resc + 0.5_dp*( temp(j,1) + temp(j-1,1) )*( yg(j)-yg(j-1) )
    end do
    if ( theta_resc < 0.0_dp ) theta_resc = abs(theta_resc)

    theta_inlet = 0.0_dp
    temp        = Umean_inlet_T/Uinf*(1.0_dp - Umean_inlet_T/Uinf)
    do j = 2, nyg
       theta_inlet = theta_inlet + 0.5_dp*( temp(j,1) + temp(j-1,1) )*( yg(j)-yg(j-1) )
    end do
    if ( theta_inlet < 0.0_dp ) theta_inlet = abs(theta_inlet)

    !---------------------------------------------------------------------!
    ! utau at i_rescale (reference lines 220-226)
    UVwall1(1) = sum( 0.5_dp*(U_resc(1,2:nzg-1)+U_resc(2,2:nzg-1))*V_resc(1,2:nzg-1) )
    call allreduce_sum_arr(UVwall1, 1)
    UVmean_wall = UVwall1(1)/real(nzg_global-2, dp)
    dUmean_wall = ( Umean_resc_T(2,1)-Umean_resc_T(1,1) )/( yg(2)-yg(1) )
    utau_resc   = abs( nu*dUmean_wall - UVmean_wall )**0.5_dp

    !---------------------------------------------------------------------!
    ! consistent utau at inlet, Ludwig-Tillmann (reference lines 229-231)
    utau_inlet = utau_resc*(theta_resc/theta_inlet)**(1.0_dp/8.0_dp)
    gamma      = utau_inlet/utau_resc

    !---------------------------------------------------------------------!
    ! wall-normal grids in plus and delta99 units (reference lines 234-244)
    ygp_inlet  = utau_inlet*yg/nu
    yp_inlet   = utau_inlet*y /nu
    ygp_resc   = utau_resc *yg/nu
    yp_resc    = utau_resc *y /nu

    etag_inlet = yg/delta_inlet
    eta_inlet  = y /delta_inlet
    etag_resc  = yg/delta_resc
    eta_resc   = y /delta_resc

    !---------------------------------------------------------------------!
    ! inner/outer weights (reference lines 247-251)
    alpha = 4.0_dp
    b     = 0.2_dp
    Weig  = 0.5_dp*(1.0_dp+tanh(alpha*(etag_inlet-b)/((1.0_dp-2.0_dp*b)*etag_inlet+b))/tanh(alpha))
    Wei   = 0.5_dp*(1.0_dp+tanh(alpha*(eta_inlet -b)/((1.0_dp-2.0_dp*b)*eta_inlet +b))/tanh(alpha))

    !---------------------------------------------------------------------!
    ! inner and outer components (reference lines 253-282)

    ! U mean
    call interp_vel(Umean_resc_T, Umean_inner,  ygp_resc,  ygp_inlet)
    call interp_vel(Umean_resc_T, Umean_outer, etag_resc, etag_inlet)
    Umean_inner = gamma*Umean_inner
    Umean_outer = gamma*Umean_outer + (1.0_dp-gamma)*Uinf

    ! V mean
    call interp_vel(Vmean_resc_T, Vmean_inner,  yp_resc,  yp_inlet)
    call interp_vel(Vmean_resc_T, Vmean_outer, eta_resc, eta_inlet)

    ! U fluctuations
    call interp_vel(Ufluc_resc, Ufluc_inner,  ygp_resc,  ygp_inlet)
    call interp_vel(Ufluc_resc, Ufluc_outer, etag_resc, etag_inlet)
    Ufluc_inner = gamma*Ufluc_inner
    Ufluc_outer = gamma*Ufluc_outer

    ! V fluctuations
    call interp_vel(Vfluc_resc, Vfluc_inner,  yp_resc,  yp_inlet)
    call interp_vel(Vfluc_resc, Vfluc_outer, eta_resc, eta_inlet)
    Vfluc_inner = gamma*Vfluc_inner
    Vfluc_outer = gamma*Vfluc_outer

    ! W fluctuations
    call interp_vel(Wfluc_resc, Wfluc_inner,  ygp_resc,  ygp_inlet)
    call interp_vel(Wfluc_resc, Wfluc_outer, etag_resc, etag_inlet)
    Wfluc_inner = gamma*Wfluc_inner
    Wfluc_outer = gamma*Wfluc_outer

    !---------------------------------------------------------------------!
    ! inlet flow field (reference lines 284-298)
    if ( iflag_fluc == 1 ) then
       Umean_inner = 0.0_dp
       Umean_outer = 0.0_dp
       Vmean_inner = 0.0_dp
       Vmean_outer = 0.0_dp
    end if
    do j = 1, nyg
       Ut(j,:) = (Umean_inner(j,1)+Ufluc_inner(j,:))*(1.0_dp-Weig(j)) + (Umean_outer(j,1)+Ufluc_outer(j,:))*Weig(j)
       Wt(j,:) = (                 Wfluc_inner(j,:))*(1.0_dp-Weig(j)) + (                 Wfluc_outer(j,:))*Weig(j)
    end do
    do j = 1, ny
       Vt(j,:) = (Vmean_inner(j,1)+Vfluc_inner(j,:))*(1.0_dp-Wei (j)) + (Vmean_outer(j,1)+Vfluc_outer(j,:))*Wei (j)
    end do

  end subroutine compute_rescaled_inflow

  !--------------------------------------------------------------------
  ! interp_vel — VERBATIM port of rescaled_inlet_bc.f90:interp_vel.
  !--------------------------------------------------------------------
  subroutine interp_vel(U_, U_interp, y_, y_interp)

    real(dp), intent(in)  :: U_(:,:)
    real(dp), intent(in)  :: y_(:)
    real(dp), intent(in)  :: y_interp(:)
    real(dp), intent(out) :: U_interp(:,:)

    integer  :: n(2), ni(2), n1, ni1, i, j, k
    real(dp) :: coef

    n   = shape(U_)
    n1  = n(1)
    ni  = shape(U_interp)
    ni1 = ni(1)

    if ( n1 /= ni1 ) stop 'Error! ni1/=ni'

    ! default value is top bc
    do i = 2, ni1
       U_interp(i,:) = U_(n1,:)
    end do

    ! interpolate
    do i = 2, ni1 ! i=1 is ghost cell
       j = 0
       do k = 1, n1
          if ( y_(k) > y_interp(i) ) then
             j = k
             exit
          end if
       end do
       if ( j == 1 ) stop 'Error! j==1'
       if ( j > 0 ) then
          coef          = ( y_(j)-y_interp(i) )/( y_(j)-y_(j-1) )
          U_interp(i,:) = coef*U_(j-1,:) + (1.0_dp-coef)*U_(j,:)
       end if
    end do

  end subroutine interp_vel

end module lund_inflow_mod
