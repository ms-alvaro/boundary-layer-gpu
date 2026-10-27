# Physics and models

This document describes what the solver in `src/` computes and which input
keys select each option. Every statement comes from the code on `main`. The
file names in parentheses point to the module that implements each item.
Input keys are written as they appear in the `.turbb` file and are parsed by
`src/param_mod.f90`.

Contents

1. [Governing equations and numerics](#1-governing-equations-and-numerics)
2. [Boundary conditions](#2-boundary-conditions)
3. [Large-eddy simulation](#3-large-eddy-simulation)
4. [Wall models](#4-wall-models)
5. [Outputs and statistics](#5-outputs-and-statistics)
6. [Features that exist only on branches](#6-features-that-exist-only-on-branches)

## Input file conventions

An input file is a list of lines `key = value`. Text after `!` is a comment.
Unknown keys are ignored. A key that is absent leaves its variable at the
default, which is zero or blank unless a table below says otherwise. Strings
are not quoted.

The parser keeps every character that follows `=` and does not remove leading
blanks (`param_mod.f90`, subroutine `getvar`). Numbers are unaffected. String
values must therefore be written without a blank after `=`, as in
`fileout =./data/BL` and `timeinflow_file =**TO DO**`. The shipped cases in
`cases/` follow this rule. The placeholder test for `timeinflow_file`
compares the first four characters of the value with `**TO`
(`ic_inflow_mod.f90`), so `timeinflow_file = **TO DO**` with a blank is read
as a file name and the run stops.

The parser enforces three requirements. `inflow_file` must be present.
`fileout` must contain a directory component, and the directory is created if
it does not exist. `filein` must be present when `init_rand = 0`. All other
keys default silently.

## 1. Governing equations and numerics

### Equations and units

The solver integrates the incompressible Navier-Stokes equations for a
constant-density fluid with kinematic viscosity $\nu$ (`rhs_kernels.cuf`),

$$
\frac{\partial u_i}{\partial t} + \frac{\partial u_i u_j}{\partial x_j}
= -\frac{\partial p}{\partial x_i}
+ \frac{\partial}{\partial x_j}\left[(\nu+\nu_t)\left(\frac{\partial u_i}{\partial x_j}\right)\right]
+ f_i,
\qquad \frac{\partial u_j}{\partial x_j} = 0 .
$$

Here $u_i = (u, v, w)$ are the velocity components along the streamwise
($x$), wall-normal ($y$) and spanwise ($z$) directions, $p$ is the pressure
divided by the density, $\nu_t$ is the eddy viscosity (zero in DNS, see
section 3) and $f_i = (\mathrm{dPdx}, 0, \mathrm{dPdz})$ is a uniform body
force set by the keys `dPdx` and `dPdz`. The convective term is discretized in
divergence form, and the viscous term in flux form with the viscosity
evaluated on the cell faces. The rotation term of the reference solver is not
ported; `Omega_z` is always zero.

The code is dimensionless. The free-stream velocity is $U_\infty = 1$
(`ic_inflow_mod.f90`). The unit of length is the distance from the plate
leading edge to the inlet, so the first face coordinate is $x(1) = 1$
(`grid_mod.f90`). The only physical input is the viscosity key `nu`, which
is $\nu/(U_\infty x_0)$ with $x_0$ the inlet distance. The inlet Reynolds
number is therefore

$$
\mathrm{Re}_{x,\mathrm{inlet}} = \frac{U_\infty x(1)}{\nu} = \frac{1}{\texttt{nu}} ,
$$

and this is the value printed as `Re_x` by the monitor (`io_mod.f90`). The
Blasius similarity variable used for the initial condition and the inflow
profile is $\eta = y\sqrt{U_\infty/(2\nu x)}$, with $u/U_\infty = f'(\eta)$
and $v/U_\infty = (\eta f' - f)/\sqrt{2\,\mathrm{Re}_x}$.

### Grid

The key `nxyz = nx ny nz` gives the number of face points in each direction.
Cell centers number one fewer, and the arrays carry one ghost point on each
side of the centers. The grid is staggered: $u$ lives on $x$ faces, $v$ on
$y$ faces, $w$ on $z$ faces and $p$ on cell centers. Finite differences are
second order on this arrangement. Interpolations from centers to $y$ faces
use the linear weights of the stretched grid; all other interpolations are
plain averages.

The key `boxsize = Lx Ly Lz alpha` sets the domain and the wall-normal
stretching (`grid_mod.f90`). The $x$ faces are uniform from $1$ to $1 + L_x$.
The $z$ faces are uniform with period $L_z$: the face $z(nz-1)$ is the
periodic image of $z(1)$, and $z(nz)$ is a ghost, so the spanwise period
contains $nz - 2$ cells. The $y$ faces follow

$$
y_j = L_y\,\frac{\sinh(\alpha\, s_j)}{\sinh \alpha},\qquad s_j = \frac{j-1}{ny-1},
$$

so cells grow away from the wall. The ratio of the last to the first spacing
approaches $\cosh\alpha$ for fine grids; $\alpha = 3$, the value of the
shipped cases, gives a ratio near 10. The formula is not uniform for
$\alpha = 1$ and is undefined for $\alpha = 0$. The key `ygrid_file` replaces
this formula by a list of $ny$ face coordinates read from a text file, one
value per line, with blank lines and lines starting with `#` ignored.

When the run starts from a snapshot (`init_rand = 0`) the face and center
coordinates stored in the file replace the generated ones (`io_mod.f90`,
subroutine `read_restart`). The ghost-extended grids, the minimum spacings
used by the Poisson solver and the interpolation weights are not recomputed
and remain those of `boxsize`, so `boxsize` must describe the same domain as
the snapshot.

The pressure solver uses the uniform spacings in $x$ and $z$, so the Poisson
transform lengths are $nx - 2$ and $nz - 2$. The half-length cosine transform
requires $nx - 2$ to be even. For multi-GPU runs $nx$ and $nz$ must be even
and $nz - 2$ must be divisible by the number of ranks.

### Time integration

The key `RKscheme` selects the explicit scheme (`timestep_mod.cuf`). Writing
$F_s$ for the right-hand side evaluated at substep $s$, each substep forms
$u_s = u^n + \Delta t \sum_q c_{sq} F_q$ on the interior points and then
applies the boundary conditions, the projection and the boundary conditions
again.

| `RKscheme` | scheme | coefficients $c_{sq}$ | substep times $t_s - t^n$ |
|---|---|---|---|
| 1 | forward Euler | $c_{11} = 1$ | $\Delta t$ |
| 2 | two-stage Runge-Kutta (midpoint) | $c_{11} = 1/2$; $c_{21} = 0$, $c_{22} = 1$ | $\Delta t/2$, $\Delta t$ |
| 3 | three-stage Runge-Kutta | $c_{11} = 8/15$; $c_{21} = 1/4$, $c_{22} = 5/12$; $c_{31} = 1/4$, $c_{32} = 0$, $c_{33} = 3/4$ | $8\Delta t/15$, $2\Delta t/3$, $\Delta t$ |

The Euler step is executed as the last substep of the three-stage table, so
its boundary conditions and projection use the substep-3 coefficients.

The key `CFL` sets the time step. A negative value means a fixed step
$\Delta t = \lvert\texttt{CFL}\rvert$, and the stability estimate is skipped.
All shipped cases use this form. A positive value computes

$$
\Delta t = \texttt{CFL}\,\min\!\left(\Delta t_{\mathrm{conv}},\;
\frac{\min(\Delta x_{\min}^2, \Delta y_{\min}^2, \Delta z_{\min}^2)}{\nu}\right),
$$

where $\Delta t_{\mathrm{conv}}$ is the smallest over the three directions of
the maximum of $\Delta x_i/\lvert u_i\rvert$ over the field. The code takes
the maximum of this ratio, as the reference solver did, and the ratio is
infinite wherever a velocity component is exactly zero. In that case the
viscous estimate sets the step. In multi-GPU runs the minimum over ranks is
taken.

### Fractional step and Poisson solver

After the explicit advance and the first boundary-condition pass, the
intermediate velocity $u^*$ is projected (`poisson_mod.cuf`),

$$
\nabla^2 \phi = \nabla\cdot u^*,\qquad u = u^* - \nabla\phi ,
$$

so that the divergence of $u$ on the interior cells is zero to round-off. The
pressure stored for output is $P = \phi/(\Delta t/2)$, taken at substep 1
of every step regardless of the scheme. The monitor prints the maximum
interior divergence as a check.

The Poisson equation is diagonalized by a type-II discrete cosine transform
in $x$ and a discrete Fourier transform in $z$, both with the modified
wavenumbers of the second-order Laplacian on the uniform spacings, and the
remaining tridiagonal system in $y$ is solved for each mode pair by the
Thomas algorithm with precomputed LU factors. The cosine transform
corresponds to a homogeneous Neumann condition for $\phi$ at the inlet and
outlet planes. The tridiagonal matrix carries the Neumann closures
$\phi(1) = \phi(2)$ at the wall and $\phi(nyg) = \phi(nyg-1)$ at the top. The
singular mean mode is pinned by scaling the first pivot by $3/2$, as in the
reference solver. The transforms run in cuFFT with a half-length real
cosine transform (Makhoul packing), and the whole projection stays on the
GPU. Multi-GPU runs use a $z$-slab decomposition with an all-to-all
transpose of the cosine coefficients; `docs/MULTI_GPU_DESIGN.md` describes
that pipeline.

The only periodic direction is $z$. The streamwise direction has an inflow
and an outflow boundary, and the wall-normal direction has a wall and a top
boundary. Section 2 describes each.

### Order of operations in a substep

Each substep runs, in this order: the eddy viscosity from the current
velocity when `LES > 0`, the wall model when `WM > 0` or `TauwModel > 0`,
the three right-hand sides, the advance, the boundary conditions, the
projection, and the boundary conditions again. The Lund recycling planes and
the random inflow planes are refreshed once per step before substep 1.

## 2. Boundary conditions

The boundary-condition kernels are in `src/bc_kernels.cuf` and the inflow
and top tables in `src/ic_inflow_mod.f90`, `src/lund_inflow_mod.f90` and
`src/hit_inflow_mod.cuf`. The kernels run in the order inflow, outflow,
spanwise periodicity, top, wall, global mass correction.

### 2.1 Inflow (plane $i = 1$)

Every inflow option writes the inlet plane as

$$
u(1,j,k) = U_{\mathrm{inlet}}(j) + U_t(j,k),
$$

and likewise for $v$ and $w$, where $U_{\mathrm{inlet}}$ is a steady profile
and $U_t$ a time-dependent plane. The key `inflow_flag` selects how the two
parts are built.

| `inflow_flag` | steady profile | time-dependent part | keys | files |
|---|---|---|---|---|
| 1 | Blasius similarity solution evaluated at $x(1)$ | sum of spanwise-temporal Fourier modes, or none | `inflow_file`, `timeinflow_file` | Blasius table (formatted); mode file (stream binary) or the placeholder `**TO DO**` |
| 2 | $u$ and $v$ profiles read from a binary file | same as 1 | `inflow_file`, `timeinflow_file` | profile file (stream binary: `int32 ny`, `U(nyg)`, `V(ny)` as `real64`, endianness detected); the Blasius table is read from `blasius_solution.dat` in the working directory |
| 3 | zero | the full plane at $i = $ `Lund_ix`, rescaled to inlet units every step | `Lund_ix`, `Lund_deltai`, `Lund_T` | `<snapshot>.mean.rescaling` companion on restart (optional) |
| 4 | Blasius as in 1 | uniform random noise redrawn every step | `inflow_file`, `Amplitude` | Blasius table |
| 5 | turbulent mean profile read from a binary file | fluctuations at $i = $ `Lund_ix`, rescaled every step | `inflow_file`, `Lund_ix`, `Lund_deltai`, `Lund_T` | profile file (stream binary, big-endian: `int32 ny`, `U(nyg)`, `V(ny)`) |
| 6 | Blasius as in 1 | planes from a precursor homogeneous-turbulence library, interpolated in time | `inflow_file`, `hit_file`, `N_buffer_hit`, `ygrid_file` (optional) | Blasius table; plane library (format v2) |

Any other value stops the run. The paragraphs below give the details of
each option.

Option 1 reads the Blasius table from `inflow_file`. The file is formatted
and holds the number of points $n$, then $\eta(1{:}n)$, $f(1{:}n)$ and
$f'(1{:}n)$, each as one list-directed record. The profiles are linear
interpolations of $f'$ and $\eta f' - f$ in $\eta$; above the tabulated
range $u$ is $1$ and $v$ keeps its edge value. If `timeinflow_file` holds a
file name, the time-dependent plane is

$$
U_t(j,k) = \sum_{n=1}^{N}\sum_{m=1}^{M}
\mathrm{Re}\left[\hat q_u(j,n,m)\,
\exp\!\big(i\,n\beta\, z_k - i\, m\omega\, t\big)\right],
$$

with $\beta$ the fundamental spanwise wavenumber, $\omega$ the fundamental
angular frequency, and $\hat q_u$ complex mode shapes interpolated to the
wall-normal grid. The sum is evaluated on the GPU at every substep time. The
mode file is a stream binary with `int32 ny_inlet`, `int32 N`, `int32 M`,
`real64 beta`, `real64 omega`, then `ymesh(ny_inlet)`, `zmode(N)`,
`tmode(M)`, and the real and imaginary parts of the $u$, $v$ and $w$ shapes,
each of size `(ny_inlet, N, M)`. `cases/transition/generate_temporal_modes_local.py`
writes this format. The value `**TO DO**` means no modes, and the inflow is
then the steady Blasius profile.

Option 2 differs from option 1 only in the steady profile. The binary profile
is tried in native byte order and then in big-endian order, and its `ny`
must equal the run's `ny`. The time-dependent modes of option 1 are applied
if `timeinflow_file` names a file.

Options 3 and 5 are the recycling method of Lund, Wu and Squires. Once per
step the $y$-$z$ planes at the recycling station `Lund_ix` and at the inlet
are copied to the host, rescaled, and uploaded as $U_t, V_t, W_t$. The
spanwise means at both stations are averaged in time with the exponential
filter

$$
\langle u\rangle_T^{n+1} = \frac{\Delta t}{T}\,\langle u\rangle^{n+1}
+ \left(1 - \frac{\Delta t}{T}\right)\langle u\rangle_T^{n},
$$

with $T$ the key `Lund_T`. The friction velocity at the recycling station is
$u_{\tau,r} = \sqrt{\lvert\nu\,\partial_y\langle u\rangle_T - \langle u v\rangle\rvert_w}$,
the one at the inlet follows the Ludwig-Tillmann relation
$u_{\tau,i} = u_{\tau,r}(\theta_r/\theta_i)^{1/8}$ with $\theta$ the momentum
thickness of the filtered means, and $\gamma = u_{\tau,i}/u_{\tau,r}$. The
inner part of every profile is interpolated in $y^+ = y u_\tau/\nu$ and the
outer part in $y/\delta$, where $\delta$ at the recycling station is the
99 percent thickness of $\langle u\rangle_T$ and $\delta$ at the inlet is the
key `Lund_deltai`. The two parts are blended with the weight

$$
W(\eta) = \frac{1}{2}\left[1 + \frac{\tanh\big(4(\eta-0.2)/(0.6\,\eta+0.2)\big)}{\tanh 4}\right],
\qquad \eta = y/\delta_i ,
$$

and the means are scaled as $\gamma\langle u\rangle$ (inner) and
$\gamma\langle u\rangle + (1-\gamma)U_\infty$ (outer). Option 3 imposes the
full rescaled plane and sets $U_{\mathrm{inlet}} = 0$. Option 5 zeroes the
rescaled means and adds the rescaled fluctuations to the mean profile read
from `inflow_file`; its $w$ plane is shifted in $z$ by 16 cells,
$W_t(j,k) \leftarrow W_t(j,\,\mathrm{mod}(k+16, nz)+1)$, as in the reference
solver. For both options the top profile is the turbulent displacement
velocity $V_{\mathrm{top}}(x) = 0.0160\,\mathrm{Re}_x^{-1/7}$. The run checks
$2 \le$ `Lund_ix` $\le nx-1$, `Lund_T` $> 0$ and `Lund_deltai` $> 0$. For
option 3 every snapshot is accompanied by `<snapshot>.mean.rescaling` holding
the filtered means, and a restart reads it back when it exists; otherwise the
means restart from the spanwise averages of the initial field.

Option 4 draws, once per step on the host, $U_t(j,k) = A\,(1 - U_{\mathrm{inlet}}(j))\,(r - 1/2)$
with $r$ uniform in $[0, 1)$ and $A$ the key `Amplitude`, and likewise for
$V_t$ and $W_t$. The factor $1 - U_{\mathrm{inlet}}$ confines the noise to the
boundary layer. The random sequence is not reproducible across compilers.

Option 6 reads a library of $y$-$z$ velocity planes produced by a
preprocessing pipeline that is not part of this repository. The file
(`hit_file`) is a stream binary with the header `int32 magic = 20260612`,
`int32 n_planes`, `int32 ny_hit`, `int32 nz_hit`, `real64 dt_plane`,
`real64 Lz_box`, seven `real64` descriptors (`Tu, L11r, alphaL, A_vel, alphat, y_blend, Uinf`),
and `real64 y_hit(ny_hit)`; then every plane stores $u$, $v$, $w$ as
`real32` arrays with $z$ fastest. The planes are interpolated to the inlet
grids on the host (linear and clamped in $y$, linear and periodic in $z$) and
held on the GPU in a buffer of `N_buffer_hit` planes (default 1000), which is
reloaded from disk when the simulation time leaves it. At every substep the
inlet plane is the linear interpolation in time between the two planes that
bracket $t$. The library period `Lz_box` must equal the domain period
$\Delta z\,(nz-2)$ to a relative tolerance of $10^{-6}$. When $t$ passes the
last plane the run stops; the library is never recycled. The key
`ygrid_file` was added for this option so that the wall-normal grid can
resolve the free-stream turbulence.

### 2.2 Wall ($y = 0$)

Without a wall model (`WM = 0`) the wall is a no-slip plate. The ghost rows
of the center-located components are set antisymmetric,
$u(i,1,k) = -u(i,2,k)$ and $w(i,1,k) = -w(i,2,k)$, so the interpolated wall
value vanishes, and the wall-normal velocity on the wall face is
$v(i,1,k) = 0$. There is no input key for wall blowing or suction; the only
nonzero wall-normal velocity at the wall comes from wall model 13.

With a wall model the slip lengths $\alpha_x, \alpha_y, \alpha_z$ returned by
the model (section 4) replace the no-slip condition. The table lists the
wall conditions by `WM` value.

| `WM` | $u$ ghost row | $v$ on the wall face | $w$ ghost row |
|---|---|---|---|
| 0 | $u(1) = -u(2)$ | $0$ | $w(1) = -w(2)$ |
| 11, 15 | $u(1) = u(2) - \alpha_x\,(y_{g,2} - y_{g,1})$ | $0$ | $w(1) = -w(2)$ |
| 13 | Robin with $\alpha_x$ | $V_{\mathrm{bottom}}(i,k)$ from the model | Robin with $\alpha_z$ |
| all other values | Robin with $\alpha_x$ | pseudo-Robin with $\alpha_y$ | Robin with $\alpha_z$ |

The Robin condition $u = \alpha\,\partial u/\partial y$ is imposed through the
ghost row,

$$
u(i,1,k) = \frac{2\alpha/\Delta_g - 1}{2\alpha/\Delta_g + 1}\,u(i,2,k),
\qquad \Delta_g = y_{g,2} - y_{g,1},
$$

which reduces to the antisymmetric no-slip ghost for $\alpha = 0$. For
`WM = 11` and `15` the quantity $\alpha_x$ is the wall velocity gradient
itself, and the ghost row imposes it as a Neumann condition. The pseudo-Robin
condition for $v$ uses the field at the beginning of the step,
$v(i,1,k) = \alpha_y\,(v^n(i,2,k) - v^n(i,1,k))/(y_2 - y_1)$, and then
subtracts the plane average of the result so that no net mass enters through
the wall.

### 2.3 Top ($y = L_y$)

The key `top_flag` selects the condition on the top plane, which is the last
row of each array ($j = nyg$ for $u$ and $w$, the face $j = ny$ for $v$).

| `top_flag` | $u$ | $v$ | $w$ | extra keys |
|---|---|---|---|---|
| 0 | $u = U_{\mathrm{top}}(i)$ | $v = V_{\mathrm{top}}(i)$ | $w = 0$ | none |
| 1 | zero spanwise vorticity: $u(i,nyg) = u(i,nyg-1) + \frac{\Delta y}{\Delta x}\,[v(i+1,ny) - v(i,ny)]$ | $v = V_{bs}(i)$ (Coleman et al. 2018) | $\partial w/\partial y = 0$ | `Vmax`, `x0`, `sigma`, `phi` |
| 2 | as 1 | $v = V_{bs}(i)$ (Abe 2017) | $\partial w/\partial y = 0$ | `Vmax`, `x0`, `sigma`, `phi` |
| 3 | $u = U_{FS}(i)$ (Falkner-Skan) | $v = 0$ | $\partial w/\partial y = 0$ | `beta_hartree` |
| 4 | $\partial u/\partial y = 0$ | $v = V_{\mathrm{top}}(i)$ | $\partial w/\partial y = 0$ | none |

Any other value stops the run. The profile $U_{\mathrm{top}} = 1$ for every
inflow option. The displacement velocity $V_{\mathrm{top}}(i)$ is the Blasius
edge value $V = U_\infty(\eta f' - f)_{\eta_{\max}}/\sqrt{2\,\mathrm{Re}_x}$
at each station for `inflow_flag` 1, 2, 4 and 6, and
$0.0160\,\mathrm{Re}_x^{-1/7}$ for `inflow_flag` 3 and 5. It is also the
value used by `top_flag = 4`, so options 0 and 4 share the same wall-normal
flux through the top and differ only in the tangential condition. Option 4
lets free-stream turbulence pass the lid and is the one used with the
precursor-turbulence inflow.

Options 1 and 2 impose a steady blowing-and-suction profile on the lid
(`ic_inflow_mod.f90`, subroutine `compute_blowsuction_top`). With
$\xi = (x_0 - x_{g,i})/\sigma$, $V_{\max}$ the key `Vmax` and $\varphi$ the
key `phi`,

$$
V_{bs}(i) = \sqrt{2}\,V_{\max}\,\xi\,e^{\,1/2 - \xi^2} + \varphi
\quad(\texttt{top\_flag} = 1),
\qquad
V_{bs}(i) = \sqrt{2}\,V_{\max}\,\xi\,e^{\,\varphi - \xi^2}
\quad(\texttt{top\_flag} = 2).
$$

The first form is equation (2.1) of Coleman, Rumsey and Spalart (2018) and
the second equation (2.1) of Abe (2017). The key `sigma` must be nonzero.
For $i \le$ `Lund_ix` the $u$ and $w$ top rows are reset to
$U_{\mathrm{top}}$ and $W_{\mathrm{top}} = 0$, so the inlet region keeps a
Dirichlet top; when `Lund_ix` is absent or zero this reset is skipped. The
reference solver also wrote one row past the $v$ array here; the port keeps
the blowing-and-suction profile over the whole lid, which is what that write
amounted to after the wall condition overwrote it.

Option 3 imposes the Falkner-Skan free stream
$U_{FS}(i) = C\,(x_{g,i} + 1)^{\beta/(2-\beta)}$ with $C$ chosen so that
$U_{FS} = 1$ at the inlet and $\beta$ the Hartree parameter from the key
`beta_hartree` (default 0, which gives a flat plate). The profile is
evaluated at the center coordinates $x_{g,i}$ and stored on the $u$ faces,
as in the reference solver. The convective outflow below uses
$U_{FS}(nx)$ as its velocity.

### 2.4 Outflow (planes $i = nx-1$ and $i = nx$)

The outflow is a first-order convective condition. With $u^n$ the velocity at
the beginning of the step, $U_c = U_{\mathrm{top}}(nx)$ the convection
velocity, and $\Delta x_o$ the spacing between the last two interior
stations,

$$
u(nx-1,j,k) = u^n(nx-1,j,k) - c_s\,\Delta t\,U_c\,
\frac{u^n(nx-1,j,k) - u^n(nx-2,j,k)}{\Delta x_o},
\qquad u(nx,j,k) = u(nx-1,j,k),
$$

and the same for $v$ and $w$. The factor $c_s$ is the substep time fraction
of the three-stage table, $8/15$, $2/3$, $1$, for every value of `RKscheme`,
which reproduces the reference solver. The last plane is a copy of the
outflow plane and its right-hand side is never used.

After the wall condition, the net volume flux through the four non-periodic
boundaries is removed by a uniform shift of the outflow plane
(`bc_kernels.cuf`, mass-conservation block),

$$
\Delta U = \frac{Q_{\mathrm{in}} - Q_{\mathrm{out}} + Q_{\mathrm{wall}} - Q_{\mathrm{top}}}{L_y\,(nzg - 3)},
\qquad u(nx-1,j,k) \leftarrow u(nx-1,j,k) + \Delta U ,
$$

where each $Q$ is the sum of the normal velocity times the face length over
the plane. The sums are fixed-order device reductions, so the result is
bitwise reproducible.

### 2.5 Spanwise ($z$)

The spanwise direction is periodic with period $L_z$ over $nz - 2$ cells.
The center-located arrays ($u$, $v$, $p$) carry three image planes,
$F(1) = F(nzg-2)$, $F(nzg-1) = F(2)$, $F(nzg) = F(3)$, and the face-located
$w$ carries two, $w(1) = w(nz-1)$, $w(nz) = w(2)$. On several GPUs the
domain is split into $z$ slabs and the image planes become halo exchanges
(`mpi_mod.f90`).

## 3. Large-eddy simulation

The key `LES` selects the subgrid-scale model (`les_mod.cuf`). The eddy
viscosity $\nu_t$ is defined at cell centers and enters the momentum
equations through the flux-form viscous term, with $\nu + \nu_t$
interpolated to each face as in the reference `equations.f90`. The eddy
viscosity is recomputed on the GPU from the current velocity before the
right-hand side of every substep. All LES kernels run on the device; the
host receives $\nu_t$ only at output steps and, for the wall models, the
wall-parallel average.

| `LES` | model | coefficient | averaging | ghost treatment in $x$ |
|---|---|---|---|---|
| 0 | none (DNS) | $\nu_t = 0$ | | |
| 1 | constant-coefficient Smagorinsky | $C_s = 0.11$ | none | |
| 2 | dynamic Smagorinsky | Germano identity, least squares | $\langle L_{ij}M_{ij}\rangle_z / \langle M_{ij}M_{ij}\rangle_z$ over the local $z$ slab | Neumann |
| 3 | dynamic Smagorinsky, no averaging | Germano identity, least squares | pointwise ratio | periodic |

Any other value stops the run.

For `LES = 1`,

$$
\nu_t = (C_s \Delta)^2\,\lvert S\rvert,\qquad
\Delta = (\Delta x\,\Delta z\,\Delta y_j)^{1/3},\qquad
\lvert S\rvert = \sqrt{2 S_{ij} S_{ij}} ,
$$

where $\Delta y_j = y_j - y_{j-1}$ is the local cell height and $S_{ij}$ is
the rate-of-strain tensor at the cell center. The diagonal components are
differences across the cell; the off-diagonal components are formed on cell
edges and averaged to the center. The reference carries a van Driest damping
factor $f = 1 - \exp(-y^+/25)$ with $u_\tau = \sqrt{\mathrm{dPdx}}$, but the
reference hard-sets `Dirichlet_nu_t = 0` after the input is read and the
factor is then $f = 1$; the port keeps this override. The key `nutBC` is
parsed but has no effect. The wall-parallel average $\langle\nu_t\rangle(y)$
is computed for the wall models and clipped at zero.

For `LES = 2` and `3` the test filter is the three-point Simpson filter
applied in $x$, $z$ and $y$. In $x$ and $z$ its two-dimensional weights are
$1/36$, $1/9$ and $4/9$ (the product of two one-dimensional
$(1/6, 2/3, 1/6)$ stencils). In $y$ the stencil is $(1/6, 2/3, 1/6)$ on
interior rows and the one-sided blends $(2/3, 1/3)$ and $(1/3, 2/3)$ on the
first and last row of the filtered window. The ratio of test to grid filter
width is $\hat\Delta/\Delta = 2$. With $\hat{(\cdot)}$ the test filter, the
Leonard term and the model tensor are

$$
L_{ij} = \widehat{u_i u_j} - \hat u_i \hat u_j,\qquad
M_{ij} = \left(\frac{\hat\Delta}{\Delta}\right)^{2}\lvert\hat S\rvert\,\hat S_{ij}
- \widehat{\lvert S\rvert S_{ij}} ,
$$

with the velocities interpolated to cell centers before the products are
formed, and

$$
\nu_t = -\frac{1}{2}\,\frac{L_{ij}M_{ij}}{M_{kl}M_{kl}}\,\lvert S\rvert .
$$

For `LES = 2` the numerator and denominator are averaged over the spanwise
cells of the slab owned by each rank before the ratio is formed, so on
several GPUs the coefficient depends on the rank count. For `LES = 3` the
ratio is taken pointwise. In both cases $\nu_t$ is clipped at zero over the
whole array after the ratio, and the wall-normal ghost rows are Neumann
copies. For `LES = 2` the filtered fields take Neumann ghosts in $x$; for
`LES = 3` they take periodic ghosts in $x$, as the reference code did, even
though the domain is not periodic in $x$. There is no wall damping in the
dynamic models. For `LES = 2` the raw, unclipped ratio
$\langle L_{ij}M_{ij}\rangle/\langle M_{ij}M_{ij}\rangle$ is what the wall
models receive as the average eddy viscosity; for `LES = 3` they receive
zero.

The LES models were compared with the CPU reference from a shared restart:
`LES = 1` over 200 steps with velocity differences near $4\times10^{-15}$,
`LES = 3` over 200 steps near $2\times10^{-15}$, and `LES = 2` over 10 steps
near $6\times10^{-15}$ in the interior (commits `318fd7d` and `5b25f5a`).
The bypass-transition case shown in the README was run with `LES = 2`.

## 4. Wall models

### Keys and execution path

| key | values | meaning |
|---|---|---|
| `WM` | 0, 1 to 9, 11 to 15 | wall-model option; 10 is not defined and any other value stops the run |
| `TauwModel` | 0, 1 | wall-stress estimate used by `WM = 11`: 0 takes $u_{\tau,\mathrm{ref}}(x) = \sqrt{\tfrac{1}{2}\,0.027\,(x/\nu)^{-1/7}}$, 1 solves the logarithmic law by Newton iteration |
| `WMnutflag` | 0, 2 | for `WM` 9, 12, 13, 14: 0 uses the fixed viscous fraction `WMnut`, 2 estimates it from the velocity at the third and fourth cells |
| `WMnut` | real | viscous fraction $f_{\mathrm{vis}}$ of the wall stress carried by the resolved gradient (`WM` 9, 12, 13, 14) |
| `alphas` | five reals | `WM = 5`: $\bar\alpha_x$, $\bar\alpha_y$, $\bar\alpha_z$, relative amplitude $a$, frequency multiplier $n_f$ |

The wall-model code is the reference `wallmodel.f90` compiled unchanged on
the host (`src/legacy_wallmodel.f90`, with `legacy_subgrid.f90`,
`legacy_interpolation.f90`, `legacy_bc.f90` and `legacy_newton.f90`). The
bridge `src/wm_bridge.f90` runs once per substep when `WM > 0` or
`TauwModel > 0`: it copies $u$, $v$, $w$, the beginning-of-step fields and
$\nu_t$ from the GPU to the host arrays of the reference, calls
`compute_wall_model`, and copies the slip-length planes $\alpha_x$,
$\alpha_y$, $\alpha_z$, the wall velocity $V_{\mathrm{bottom}}$ and $\nu_t$
back to the GPU. Several models overwrite $\nu_t$ with the eddy viscosity of
the test-filtered field, and that overwrite is carried back too. The wall
kernels of section 2.2 then apply the planes. This path moves several full
fields across the PCIe bus at every substep and is slow compared with the
DNS loop.

The models were written for a channel. Models 1, 2, 3, 4 and 7 average their
wall quantities over both the bottom and the top rows of the domain, and the
log-law models 6 and 8 use $u_\tau = \sqrt{\mathrm{dPdx}}$. In this solver
only the bottom boundary receives the resulting condition, and `dPdx` is
zero in a boundary layer. The top-row contributions and the zero friction
velocity are inherited from the channel configuration of the reference.

In the models below $\langle\cdot\rangle$ is a wall-parallel average,
$\tau_{uv} = -2\langle\nu_t\rangle S_{uv}$ is the modeled shear stress,
$\hat{(\cdot)}$ is the test filter of section 3 with ratio
$r = \hat\Delta/\Delta = 2$, $\Delta_1 = y_2 - y_1$ is the first cell height,
and $C_f^{\mathrm{ref}} = 0.027\,\mathrm{Re}_x^{-1/7}$ with
$\tau_w^{\mathrm{ref}} = \tfrac{1}{2} C_f^{\mathrm{ref}} U_\infty^2$ is the
flat-plate correlation of White. All gradients are first-order differences
between the wall ghost row and the first cell.

### Options

| `WM` | model | what it computes | matching location | returned condition |
|---|---|---|---|---|
| 1 | momentum integral | $\alpha_x = \langle\int_0^{l}(\hat u\hat v + T_{uv} - uv - \tau_{uv})\,dy\rangle / [\nu(1-r)\langle\partial_y u\rangle_w]$; the integration length $l$ is a variable the input path never sets, so the integral and $\alpha_x$ are zero | wall and the band $y \le l$ | Robin, $\alpha_y = \alpha_z = \alpha_x$ |
| 2 | dynamic slip (Bose-type) | $\alpha^2 = \langle L_{ij}M_{ij}\rangle/\langle M_{ij}M_{ij}\rangle$ with $L_{ij} = \widehat{u_iu_j} - u_iu_j - T_{ij} + \hat\tau_{ij}$ and $M_{ij} = r^2\,\widehat{\partial_n u_i}\,\widehat{\partial_n u_j} - \partial_n u_i\,\partial_n u_j$, clipped at zero, square root | first two cells at the wall (and top) | Robin, single $\alpha$ |
| 3 | least-squares Robin fit | $\alpha = \sum_i\langle u_i\,\partial_y u_i + r\,\hat u_i\,\partial_y\hat u_i\rangle / \sum_i\langle(\partial_y u_i)^2 + r^2(\partial_y\hat u_i)^2\rangle$, clipped at zero | wall face | Robin, single $\alpha$ |
| 4 | momentum balance at the wall | $\alpha^2 = (\nu\langle\partial_y u\rangle - \mathrm{dPdx}\,L_y/2 - \langle\tau_{uv}\rangle)/\langle\partial_y u\,\partial_y v\rangle$, clipped, square root | wall face | Robin, single $\alpha$ |
| 5 | prescribed slip | $\alpha_i = \bar\alpha_i\,[1 + a\sin(\omega t)]$, $\omega = 2\pi n_f\sqrt{\mathrm{dPdx}}/\Delta_1$, from `alphas` | none | Robin |
| 6 | log layer with resolved $\nu_t$ | $\alpha = \dfrac{\langle\nu_t\rangle_w + \nu}{\kappa u_\tau}\ln\!\left(\dfrac{\Delta_1 u_\tau}{\nu}e^{\kappa B}\right)$, $\kappa = 0.38$, $B = 5.2$, $u_\tau = \sqrt{\mathrm{dPdx}}$ | first cell | Robin, single $\alpha$ |
| 7 | dynamic momentum balance | $\alpha^2 = [\nu(\langle\partial_y\hat u\rangle - \langle\partial_y u\rangle) + \langle\tau_{uv}\rangle - \langle T_{uv}\rangle]/[r^2\langle\partial_y\hat u\,\partial_y\hat v\rangle - \langle\partial_y u\,\partial_y v\rangle]$, clipped, square root | wall face | Robin, single $\alpha$ |
| 8 | log layer with fixed $\nu_t$ | $\alpha = \Delta_1\ln\!\left(\dfrac{\Delta_1 u_\tau}{\nu}e^{\kappa B}\right)$ and the ghost rows of $\langle\nu_t\rangle$ set to $\nu_t = \kappa u_\tau\Delta_1$ at the wall | first cell | Robin, single $\alpha$ |
| 9 | wall stress from $C_f$ through $v$ | $\alpha_x(x) = \nu\langle u\rangle_w/(f_{\mathrm{vis}}\tau_w^{\mathrm{ref}})$ floored at $10^{-3}$; $\alpha_y(x) = [\nu\langle\partial_y u\rangle - \tau_w^{\mathrm{ref}} - \langle\tau_{uv}\rangle + \beta_y\langle u\rangle_w]/\langle u\,\partial_y v\rangle$ with $\beta_y$ the wall-plane mean removed by the $v$ condition | wall face, per $x$ station | Robin $u$, $w$; pseudo-Robin $v$; $\alpha_z = \alpha_y$ |
| 11 | Neumann from the momentum balance | $\alpha_x(x) = (u_{\tau,\mathrm{model}}^2 + \langle\tau_{uv}\rangle)/\nu$ with $u_{\tau,\mathrm{model}}$ from `TauwModel`; floored at zero | wall face, per $x$ station | Neumann $\partial_y u = \alpha_x$; $v = 0$; $w$ no-slip |
| 12 | wall stress from $C_f$ through $u$ | $\alpha_x$ as in 9; $\alpha_y(x) = [\nu\langle\partial_y u\rangle - \tau_w^{\mathrm{ref}} - \langle\tau_{uv}\rangle]/(\alpha_x\langle\partial_y u\,\partial_y v\rangle)$ capped at $\lvert\alpha_y\rvert \le 2\times10^{-3}$ | wall face, per $x$ station | Robin $u$, $w$; pseudo-Robin $v$; $\alpha_z = \alpha_y$ |
| 13 | wall transpiration from $C_f$ | $\alpha_x$ as in 9 (floored at zero); $V_{\mathrm{bottom}}(x) = [\nu\langle\partial_y u\rangle - \tau_w^{\mathrm{ref}} - \langle\tau_{uv}\rangle]/\langle u\rangle_w$ clipped to $\lvert V\rvert \le 0.5$; $\alpha_z = \alpha_x/10$ | wall face, per $x$ station | Robin $u$, $w$; Dirichlet $v = V_{\mathrm{bottom}}$ |
| 14 | same routine as 13 | identical quantities to 13 | as 13 | the wall kernels do not single out 14, so $v$ takes the pseudo-Robin condition with $\alpha_y = V_{\mathrm{bottom}}$ instead of the Dirichlet value |
| 15 | exact Neumann from $C_f$ | $\alpha_x(x) = \tau_w^{\mathrm{ref}}(x)\big/\langle\nu + \nu_t\rangle_w$ with the spanwise average of $\nu + \nu_t$ in the first cell, $U_\infty = 1$ | wall face, per $x$ station | Neumann $\partial_y u = \alpha_x$; $v = 0$; $w$ no-slip |

For models 9, 12 and 13 the free-stream velocity in $C_f^{\mathrm{ref}}$ is
taken as the top-row value $u(2, nyg, 2)$. With `WMnutflag = 2` the viscous
fraction is $f_{\mathrm{vis}} = \nu\langle u\rangle_{3,4}/(\tau_w^{\mathrm{ref}} y_3)$,
where $\langle u\rangle_{3,4}$ is the average of the third and fourth cells.
Models 13 and 14 keep running averages of their wall quantities with the time
constant `alpha_std`, but the code then overwrites the averages with the
instantaneous values, so the time constant has no effect.

The option `TauwModel = 1` solves, at every $x$ station and for the velocity
at the third cell $y_{g,3}$,

$$
\frac{u_\tau}{\kappa}\ln\!\left(\frac{y_{g,3}\,u_\tau}{\nu}\right) + B\,u_\tau = \langle u\rangle(y_{g,3}),
\qquad \kappa = 0.38,\; B = 5.2,
$$

by Newton iteration with at most 30 iterations and tolerance $10^{-10}$,
starting from the previous value (`legacy_newton.f90`). A warning is printed
when the iteration count reaches 30. The result is used only by `WM = 11`.

### What is validated

The validation record is in the commit messages `e1e765a` and `efd413b` on
`main`. Validation means agreement with the pre-port CPU solver started from
the same restart file, after two memory errors in the reference were
repaired in both codes (periodic-$z$ buffer identifiers in models 1, 3 and 7,
and the allocation of the filtered average eddy viscosity).

| `WM` | laminar probe, 40 steps, 2 ranks | turbulent probe, 20 steps, 2 ranks |
|---|---|---|
| 1, 2 | machine precision ($\sim 10^{-15}$) | machine precision ($1.3\times10^{-15}$) |
| 3, 5, 7, 11, 15, and `TauwModel = 1` | machine precision | not reported |
| 4 | division by zero in both codes | $2.9\times10^{-10}$ (its ratio amplifies round-off; the rank-resolved terms agree to $10^{-12}$) |
| 6, 8, 9, 12, 13, 14 | division by zero in both codes | division by zero in both codes |

Models 6 and 8 divide by $u_\tau = \sqrt{\mathrm{dPdx}}$, which is zero in
every boundary-layer input, and models 9, 12, 13 and 14 divide by
wall-plane averages that vanish on the probes used. These six options are
ported and run identically to the reference as far as the reference runs,
but no case in this repository exercises them to a physical result. No wall
model is physically validated against boundary-layer data here; the record
above is code-to-code agreement only. The automated test `tests/validate.py`
covers the laminar DNS case and does not run a wall model.

## 5. Outputs and statistics

All output is handled by `src/io_mod.f90`. File names are built from the key
`fileout` and the absolute step number, which is the loop step plus the
offset `nstep_init`. The offset is the key `init_step` when given, otherwise
the step stored in the restart file, otherwise zero.

### Monitor

Every `nmonitor` steps, and at step 1, the root rank prints the step, the
time, the time step, the mean skin-friction coefficient over the plate, the
maximum and mean of each velocity component over the interior, the imposed
pressure gradients, the inlet Reynolds numbers $\mathrm{Re}_x$,
$\mathrm{Re}_\theta$, $\mathrm{Re}_\delta$ and $\delta_{99}$ computed by the
statistics, the maximum interior divergence and the elapsed wall-clock time
since the previous monitor. The line `No wall model` is printed regardless
of `WM`; the monitor was not extended when the wall models were added. With
`BL_PROFILE=1` in the environment the driver also prints every 100 steps the
average time of the first substep split into copy, right-hand side,
advance, boundary conditions and projection; the fences this needs slow the
run.

### Statistics

Every `nstats` steps, and at step 1, the solver writes
`<fileout>.<step>.stats.txt`. The statistics are spanwise averages at the
points ($x$ face, $y$ face): $u$ is interpolated to $y$ faces with the
stretched-grid weights, $v$ to $x$ faces, $w$ to $x$ faces and then $y$
faces, and $p$ is averaged as stored. The wall velocity gradient is the
difference between the first cell and the ghost row, and the skin-friction
coefficient is

$$
C_f(x) = \frac{2}{U_\infty^2}\left(-\langle u v\rangle_w + \nu\,\langle\partial_y u\rangle_w + \tau_{\mathrm{sgs},w}\right),
\qquad U_\infty = 1,
$$

where $\tau_{\mathrm{sgs},w}$ is the mean eddy viscosity near the wall times
the wall gradient and is zero in DNS. The inlet momentum thickness and
$\delta_{99}$ come from the mean profile at $i = 1$. The file is formatted:
the first line holds `%`, $t$, $\nu$ (`2F15.8`) and $nx$, $ny$, $nz$, the
step (`4I8`); then $C_f(1{:}nx)$, $x(1{:}nx)$ and $y(1{:}ny)$, one record
each; then, one $x$ station per line with $ny$ values, the blocks
$\langle u\rangle$, $\langle v\rangle$, $\langle w\rangle$,
$\langle u^2\rangle$, $\langle v^2\rangle$, $\langle w^2\rangle$,
$\langle p\rangle$, $\langle p^2\rangle$. The mean eddy viscosity is
accumulated but not written.

### Snapshots and restart

Every `nsave` steps the solver writes the binary snapshot `<fileout>.<step>`
with the step as eight digits, and refreshes the symbolic link
`<fileout>.restart` to point at it. The file is a stream of native-endian
records: $t$, $\nu$ (`real64`); the marker `-73` and the step (`int32`);
the six grids, each as its length (`int32`) followed by the coordinates
(`real64`); then $u$, $v$, $w$, $\nu_t$ and $p$, each as three `int32`
dimensions followed by the `real64` array. The $u$, $v$ and $\nu_t$ blocks
store the first $nzg - 1$ spanwise planes, $w$ the first $nz - 1$, and $p$
all $nzg$ planes. The $\nu_t$ block is zero in DNS. Multi-GPU runs write one
global file of the same layout.

A run starts from a snapshot when `init_rand = 0` and `filein` names the
file. The reader checks the six grid lengths against `nxyz`, overwrites the
generated grids with the stored ones, reads $t$, $u$, $v$ and $w$, and
rebuilds the spanwise image planes. The $\nu_t$ and $p$ blocks are not read.
Snapshots are interchangeable between rank counts. With `init_rand = 1` the
initial field is the Blasius solution on the whole domain at $t = 0$, with
$w = 0$, and `filein` must be absent.

### Box output

The keys `boxout_every`, `boxout_start`, `boxout_n`, `boxout_dir`,
`boxout_i0`, `boxout_i1`, `boxout_is` and `boxout_jmax` write up to eight
index windows of $u$, $v$, $w$ in single precision every `boxout_every`
steps once the absolute step reaches `boxout_start`. Box $b$ spans the
$x$ indices `boxout_i0(b)` to `boxout_i1(b)` with stride `boxout_is(b)`, the
$y$ indices 1 to `boxout_jmax(b)` and the first $\min(256, nz-2)$ spanwise
faces. Each box writes one stream file `<boxout_dir>/box<b>_<first-step>.bin`
per launch with a twelve-integer header (magic `20260709`, version 1, box
id, `i0`, `i1`, `is`, `jmax`, the three array sizes, `every`, 3), the
coordinates of the window as `real64`, and then per record $t$ (`real64`),
the step (`int32`) and the three `real32` arrays. The directory must exist.
The gate is on the absolute step, so the cadence continues across restarts.

### Strided three-dimensional snapshots

The keys `snap3d_every`, `snap3d_stride` (default 2) and `snap3d_dir`
(default `./snapshots`) write the three velocity components over the whole
domain in single precision with the given stride in $x$ and $z$ and full
resolution in $y$, every `snap3d_every` absolute steps. Each file
`<snap3d_dir>/snap_<step>.bin` starts with the magic `20260803`, the three
array shapes, the stride and the step (`int32`), then $t$ (`real64`) and
the three `real32` arrays on their native staggered grids. A one-time file
`snap_grid.bin` holds the strided coordinates. On `main` the fields are
gathered to the root rank in double precision and strided there; the branch
in section 6 strides them on the GPU.

### Environment variables and build options

`BL_PROFILE=1` enables the substep timers described under Monitor.
`BL_LES_DBG` set to any value prints one Leonard and model tensor sample per
substep for the dynamic models. `BL_POISSON_GATHER=1` selects the gathered
Poisson solve on several ranks instead of the distributed transpose. A build
with `make DEBUG=1` adds bounds checks and a launch-error check after every
kernel.

## 6. Features that exist only on branches

Three remote branches carry code that `main` does not have. The commands
`git log main..origin/<branch>` and `git diff main origin/<branch>` show the
exact changes.

### snap3d-devpack

One commit on top of `main` (`6249b6d`). It adds a flat plate whose leading
edge lies inside the domain, and moves the striding of the `snap3d`
snapshots to the GPU. The new keys are optional; without them the code paths
are those of `main`.

| key | meaning |
|---|---|
| `x_le` | streamwise position of the leading edge of a zero-thickness plate. For $x < x_{le}$ the bottom is a free-slip symmetry plane ($u(1) = u(2)$, $w(1) = w(2)$, $v = 0$); for $x \ge x_{le}$ it is the no-slip plate. The inlet mean is uniform plus the potential-flow correction of the plate displacement, and the Blasius initial condition and the top displacement velocity are measured from $x_{le}$. Only with `WM = 0`. |
| `le_open` | 1 opens the bottom upstream of $x_{le}$ for the precursor-turbulence inflow (`inflow_flag = 6`): the fluctuations crossing $y = 0$ are the unmirrored library planes delayed by the Taylor time $(x - x_{\mathrm{inlet}})/U_\infty$, and the net volume flux through the open patch is removed every step. |
| `le_ramp` | $n > 0$ blends the bottom condition linearly from free slip to no slip over the $n$ cells before $x_{le}$ and takes the open-bottom weight from 1 to 0 over the same cells; 0 gives a sharp edge. |

The new module `src/snap3d_pack.cuf` strides the fields on the device and
gathers only the stored samples, so the bytes written equal those of the
host path on `main`. The commit message records its use for a free-stream
turbulence bypass-transition DNS on a $5672\times239\times194$ grid on two
A100 GPUs and a zero-perturbation control that stays laminar to
$\mathrm{Re}_x = 10^6$ with $C_f$ within 0.42 percent of Blasius.

### top-potential-bc

This branch contains `snap3d-devpack` plus three commits (`2e1a809`,
`1c58e66`, `71c2a16`) that add a potential-flow matching condition on the top
plane, implemented in the new module `src/topbc_mod.cuf` and activated by
`top_potential` with `top_flag = 4`. The fluctuation about the steady
potential flow is matched to an irrotational field that decays above the
lid: every top-plane mode with wavenumber $k$ satisfies
$\partial_y \hat q = -\lvert k\rvert\,\hat q$ for $q = u, v, w$. Discretely,
the lid value is the value one cell below times the per-cell decay factor
of a harmonic mode of the solver's Laplacian, using the modified wavenumbers
of the pressure solve; $x$ is extended evenly about the inlet and outlet
planes and $z$ is periodic. The lid row of $v$ is computed at the start of
each substep from the projected field and held through the projection.

| key | values | meaning |
|---|---|---|
| `top_potential` | 0 (default), 1, 2 | 0 off; 1 match every mode $(k_x, k_z) \ne (0,0)$; 2 match only the spanwise-uniform modes |
| `top_potential_long` | 0 (default), 1, 2 | treatment of the spanwise-uniform modes with $k_x L_y < 1$: 0 keeps the `top_flag = 4` lid for them, 1 matches them, 2 matches them with a first-order relaxation |
| `top_potential_flux` | 0, 1 (default) | the mean mode of the lid: 1 lets the total lid flux float ($\partial_y v = 0$ for the mean), 0 pins it to $V_{\mathrm{top}}$ |
| `top_potential_tau` | real (default 0.5) | relaxation time for `top_potential_long = 2`; the filter state is stored in `<snapshot>.topbc` |
| `top_potential_vmean` | file | a row $\delta V(1{:}nxg)$ added to $V_{\mathrm{top}}$, the mean excess transpiration of the actual layer from a previous run (`tools/make_top_vmean.py`) |
| `top_potential_ref` | file | a reference top row $U_{\mathrm{ref}}(1{:}nx)$ from a steady laminar solution of the same grid (`tools/make_top_reference.py`); otherwise the analytic potential flow of the plate displacement is used |

The sub-keys are rejected when `top_potential` is absent, and
`top_potential` requires `top_flag = 4`. The commit messages record that the
option-off path is byte-identical to `6249b6d`, that `top_potential_long = 1`
and `2` were measured on the production configuration and rejected for box-wide
free-stream flicker, and that `top_potential_long = 0` with
`top_potential_vmean` is the production setting. Multi-GPU runs all-gather
the top rows.

### gpu-openacc

This branch holds the pre-port solver: the CPU MPI code and its OpenACC
port in a flat directory of Fortran files (`boundary_conditions.f90`,
`equations.f90`, `projection.f90`, `wallmodel.f90`, `subgrid.f90`,
`rescaled_inlet_bc.f90` and the rest), built with gfortran and FFTW or with
nvfortran. `main` dropped this tree when it became CUDA-only
(`9fc23f0`). Beyond the history that `main` shares, the branch has two
commits: the fatal stop when the precursor-turbulence library is exhausted
(`84f3669`, which replaced a silent recycling of the planes) and three check
figures from a production run (`db6d2c7`). The CUDA port on `main` includes
the no-wrap stop. The branch `cuda-port` has no commits beyond `main`.
