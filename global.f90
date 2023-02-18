!-----------------------------------------!
! Module with all shared global variables !
!-----------------------------------------!
Module global

  ! General Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use, Intrinsic :: iso_c_binding

  ! prevent implicit typing
  Implicit None

  ! FFTW
  Include 'fftw3-mpi.f03'  

  !------------------Declarations-----------------!

  ! step number
  Integer(Int32) :: istep, rk_step, itime_step
  Real   (Int64) :: time1, time2

  ! constants
  Real(Int64) :: pi = 4d0*datan(1d0)  

  ! files
  Character(200) :: filein, fileout
  Integer(Int32) :: nsave, nmonitor

  ! set random initial condition
  Integer(Int32) :: random_init

  ! domain size
  Real(Int64) :: Lx, Lz, Ly, Lxp, Lzp
  Real(Int64) :: Lx_rand, Ly_rand, Lz_rand, alpha_rand

  ! steps
  Integer(Int32) :: nsteps, nstep_init
  Real   (Int64) :: dt, t, dt_period

  ! viscosity
  Real(Int64) :: nu

  ! Reynolds numbers and delta99
  Real(Int64) :: Rex_inlet, Retheta_inlet, Redelta_inlet, delta99_inlet_ins, Retau_inlet

  ! global face points
  Integer(Int32) :: nx_global, ny_global, nz_global

  ! local face points
  Integer(Int32) :: nx, ny, nz

  ! global center points
  Integer(Int32) :: nxm_global, nym_global, nzm_global

  ! local center points
  Integer(Int32) :: nxm, nym, nzm

  ! global center points + ghost cells
  Integer(Int32) :: nxg_global, nyg_global, nzg_global

  ! local center points + ghost cells
  Integer(Int32) :: nxg, nyg, nzg

  ! global grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x_global, y_global, z_global

  ! local grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x, y, z

  ! global grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm_global, ym_global, zm_global

  ! local grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm, ym, zm

  ! global grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg_global, yg_global, zg_global

  ! local grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg, yg, zg

  ! middle points for yg->yg_m and yg_m->yg_mm
  Real(Int64), Allocatable, Dimension(:) :: yg_m, yg_mm

  ! local velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: U,V,W,P
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uo,Vo,Wo,Po
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uoo,Voo,Woo

  ! local auxiliary 
  Real(Int64), Allocatable, Dimension(:,:,:) :: term_1, term_2, term_3, term_4, term

  ! local rhs for velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:)   :: px_bottom, px_top
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_p
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uo, rhs_vo, rhs_wo
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uf, rhs_vf, rhs_wf

  ! local rhs for pressure in Fourier
  Complex(Int64), Dimension(:,:,:), Allocatable :: rhs_p_hat
  Complex(Int64), Dimension(:),     Allocatable :: rhs_aux

  ! local auxiliary arrays for MPI_sendrev boundary conditions
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ui, buffer_vi, buffer_ci
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ue, buffer_ve, buffer_we, buffer_wi, buffer_ce
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_p

  ! inlet, outlet and top boundary conditions
  Real   (Int64), Allocatable, Dimension(:)   :: U_inlet,  U_outlet, U_top
  Real   (Int64), Allocatable, Dimension(:)   :: V_inlet,  V_outlet, V_top
  Real   (Int64), Allocatable, Dimension(:)   :: W_inlet,  W_outlet, W_top
  Real   (Int64), Allocatable, Dimension(:,:) :: Ut_inlet, Vt_inlet, Wt_inlet, V_bottom

  ! Lund's rescaling boundary condition
  Real   (Int64), Allocatable, Dimension(:,:) :: Umean_resc_To,       Vmean_resc_To 
  Real   (Int64), Allocatable, Dimension(:,:) :: Umean_resc_T,        Vmean_resc_T 
  Real   (Int64), Allocatable, Dimension(:,:) :: Umean_resc_T_local,  Vmean_resc_T_local
  Real   (Int64), Allocatable, Dimension(:,:) :: Umean_inlet_T,       Vmean_inlet_T 
  Real   (Int64), Allocatable, Dimension(:,:) :: Umean_inlet_T_local, Vmean_inlet_T_local
  Integer(Int32) :: i_rescale, step_beginning
  Real   (Int64) :: delta_inlet, T_resc

  ! local auxiliary arrays for MPI_sendrev interior planes
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_us, buffer_ur
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_vs, buffer_vr
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_ws, buffer_wr
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_ps, buffer_pr
  
  ! local auxiliary planes for FFTW
  Type(C_PTR) :: cplane_fft
  Complex(C_DOUBLE_COMPLEX), Pointer,     Dimension(:,:) :: plane, plane_hat
  Complex(C_DOUBLE_COMPLEX), Allocatable, Dimension(:,:) :: plane_short

  ! Fourier points and wave numbers 
  Integer(C_INTPTR_T) :: nxp_global, nxpe_global, nzp_global, local_k_offset
  Integer(C_INTPTR_T) :: nxp, nxpe, nzp
  Integer(C_INTPTR_T) :: mx_global, mz_global
  Integer(C_INTPTR_T) :: mx, mz
  Real   (Int64)      :: dx, dz
  Real   (Int64), Dimension(:), Allocatable :: kxx, kzz

  ! Mappings for fft modes
  Integer(Int64), Dimension(:),   Allocatable :: imode_map, kmode_map
  Integer(Int64), Dimension(:,:), Allocatable :: imode_map_fft, kmode_map_fft
  
  ! FFTW plans
  Integer(C_INTPTR_T) :: alloc_local
  Type   (C_PTR)      :: plan_d, plan_i

  ! finite differences (second derivative)
  Real(Int64) :: ddx1, ddx2, ddx3
  Real(Int64) :: ddy1, ddy2, ddy3
  Real(Int64) :: ddz1, ddz2, ddz3

  ! linear solver
  Integer (Int32) :: nr, nrhs
  Integer (Int32), Dimension(:),   Allocatable :: pivot  
  Complex (Int64), Dimension(:),   Allocatable :: D, DL, DU
  Complex (Int64), Dimension(:,:), Allocatable :: M, Dyy

  ! pressure gradients
  Real(Int64) :: dPdx, dPdy, dPdz
    
  ! CFL parameters
  Real(Int64) :: CFL, dxmin, dymin, dzmin

  ! interpolation weights 
  Integer(Int32) :: in1, in2
  Real(Int64), Dimension(:), Allocatable :: weight_y_0, weight_y_1

  ! actual pressure boundary conditions
  Real   (Int64) :: coef_bc_1, coef_bc_2, Uc
  Real   (Int64), Dimension(:,:), Allocatable :: bc_1,     bc_2
  Complex(Int64), Dimension(:,:), Allocatable :: bc_1_hat, bc_2_hat
  Logical(Int32) :: pressure_computed

  ! statistics
  Integer(Int32) :: nstats
  Real   (Int64), Dimension(:,:), Allocatable ::  Umean,  Vmean,  Wmean, Pmean
  Real   (Int64), Dimension(:,:), Allocatable :: U2mean, V2mean, W2mean, UVmean, P2mean, nu_t_mean
  Real   (Int64), Dimension(:),   Allocatable :: Cf, dUdy_wall, UV_wall     ! auxiliary arrays 
  Real   (Int64), Dimension(:),   Allocatable :: Uaux_1_local, Uaux_2_local ! auxiliary arrays 
  Real   (Int64), Dimension(:),   Allocatable :: Uaux_1,       Uaux_2       ! auxiliary arrays 

  ! Runge-Kutta 3 coefficients and buffers
  Real(Int64), Dimension(:),     Allocatable :: rk_t, rk2_t
  Real(Int64), Dimension(:,:),   Allocatable :: rk_coef, rk2_coef
  Real(Int64), Dimension(:,:,:), Allocatable :: Fu1, Fu2, Fu3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fv1, Fv2, Fv3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fw1, Fw2, Fw3

  ! sgs model
  Integer(Int32) :: LES_model, Dirichlet_nu_t
  Real   (Int64), Allocatable, Dimension(:,:,:,:) :: Lij, Mij, Sij, ten_buf
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: Uf, Vf, Wf, Uff, Vff, Wff, S
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: avg_nu_t, avg_nu_t_hat
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: nu_t, nu_to, nu_too
  Real   (Int64) :: fil_size

  ! channel spanwise rotation
  Real(Int64) :: Omega_z

  ! slip length wall-model
  Integer(Int32) :: penetration, iwall_model, iwall_model_nut, istress_model
  Real   (Int64) :: alpha_mean_x, alpha_mean_y, alpha_mean_z, freq_mult, alpha_std, int_len
  Real   (Int64) :: mean_alpha_x, mean_alpha_y, mean_alpha_z, beta_y
  Real   (Int64), Allocatable, Dimension(:,:,:) :: alpha_x , alpha_y , alpha_z
  Real   (Int64), Allocatable, Dimension(:,:,:) :: alpha_xo, alpha_yo, alpha_zo
  Real   (Int64), Allocatable, Dimension(:)     :: dUmean_wall_T, mtau_T, UV_wall_T
  Real   (Int64) :: frac_vis_wall_model, frac_vis

  ! stress model
  Real   (Int64) :: Umean_model, kappa_model, B_model, yg_model
  Real   (Int64), Allocatable, Dimension(:) :: utau_model, utau_ref, utau_wall, utau_wall_T

  ! boundary conditions flags and parameters
  Integer(Int32) :: inflow_boundary_flag, top_boundary_flag
  Real   (Int64) :: Amplitude_perturbations
  Real   (Int64) :: beta_hartree ! for Falkner-Skan 

  ! Blasius boundary layer for boundary conditions
  ! input file
  Character(200) :: file_inflow, file_temporal_inlet
  Integer(Int32) :: ny_inlet, n_modes_inlet, m_modes_inlet
  Real   (Int64) :: beta_inlet, omega_inlet
  Real   (Int64), Allocatable, Dimension(:)     :: ymesh_inlet, zmode_inlet, tmode_inlet
  Complex(Int64), Allocatable, Dimension(:,:,:) :: qu_inlet,   qv_inlet,   qw_inlet
  Complex(Int64), Allocatable, Dimension(:,:,:) :: qu_inlet_o, qv_inlet_o, qw_inlet_o ! temporal
  Real   (Int64), Allocatable, Dimension(:,:,:) :: qu_inlet_r, qv_inlet_r, qw_inlet_r ! temporal
  Real   (Int64), Allocatable, Dimension(:,:,:) :: qu_inlet_i, qv_inlet_i, qw_inlet_i ! temporal

  ! statistics for z modes
  Integer(Int32) :: i_stat, j_stat, Delta_i_stat, nstats_zmodes
  Integer(Int32) :: nxu_reduced, nyu_reduced, nzu_reduced, nzu_modes, nzu_first_modes
  Real   (Int64), Allocatable, Dimension(:)     :: xu_reduced, yu_reduced, zu_reduced
  Real   (Int64), Allocatable, Dimension(:,:)   :: U_plane     
  Complex(Int64), Allocatable, Dimension(:,:)   :: U_plane_hat_z
  Real   (Int64), Allocatable, Dimension(:,:,:) :: U_reduced, rhs_u_reduced
  Complex(Int64), Allocatable, Dimension(:,:,:) :: U_reduced_hat_z, rhs_u_reduced_hat_z

  ! top boundary condition for V
  Real(Int64) :: Vbs_max, x_bs, sigma_bs, phi_bs
  
End Module global
