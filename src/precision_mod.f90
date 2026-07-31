!------------------------------------------------------------------------------
! precision_mod — floating-point kind, constants, launch-extent helper
!
! Part of the CUDA Fortran port of the boundary-layer DNS solver
! (see docs/CUDA_PORT_DESIGN.md).
!
! Provenance:
!   - dp mirrors the reference code's Real(Int64) working precision
!     (global.f90 uses iso_fortran_env Int64 kind values for 8-byte reals).
!   - pi mirrors global.f90: `pi = 4d0*datan(1d0)`. The literal below rounds
!     to the identical IEEE-754 double (0x400921FB54442D18), so numerics are
!     bit-identical to the reference.
!   - ceil_div is the small helper suggested by the design document for
!     computing 1-D kernel launch extents (number of thread blocks needed to
!     cover n items with block size b). Host-only, no CUDA dependencies.
!------------------------------------------------------------------------------
module precision_mod

  use iso_fortran_env, only: real64

  implicit none
  private

  public :: dp, pi, ceil_div

  !> Working precision of the whole solver: IEEE double (8-byte real).
  integer, parameter :: dp = real64

  !> pi to full double precision; bit-identical to 4d0*datan(1d0)
  !! used by the reference (global.f90).
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

contains

  !> Ceiling division: smallest integer nb with nb*b >= n.
  !!
  !! Used to compute launch extents for 1-D/2-D/3-D kernel grids, e.g.
  !!   grid = dim3( ceil_div(ihi-ilo+1, 64), ceil_div(jhi-jlo+1, 4), khi-klo+1 )
  !! per the binding thread-mapping convention (tBlock = dim3(64,4,1)).
  pure integer function ceil_div(n, b)
    integer, intent(in) :: n  !< number of items to cover (n >= 0)
    integer, intent(in) :: b  !< block size (b >= 1)
    ceil_div = (n + b - 1) / b
  end function ceil_div

end module precision_mod
