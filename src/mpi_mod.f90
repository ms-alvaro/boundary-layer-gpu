!------------------------------------------------------------------------------
! mpi_mod — multi-GPU (z-slab) communication layer, phase P5.1
! (docs/MULTI_GPU_DESIGN.md). One MPI rank per GPU, CUDA-aware OpenMPI:
! all exchange buffers are DEVICE arrays passed straight to MPI (the NVHPC
! OpenMPI build is CUDA-aware; z-planes and z-plane PAIRS of the solver's
! arrays are contiguous, so no packing or staging is needed at all — the
! legacy host bounce buffers existed only because the CPU code was not
! CUDA-aware).
!
! Semantics are verbatim ports of the legacy exchanges:
!   * halo_exchange_z   <- boundary_conditions.f90 update_ghost_interior_planes
!       send plane (last-1) up   -> neighbor's plane 1        (centers & faces)
!       send plane 2        down -> neighbor's plane (last)
!     where "last" = nzg for center-staggered fields (U, V and pressure) and
!     nz for the face-staggered W.
!   * periodic_exchange_z <- boundary_conditions.f90 apply_periodic_bc_z
!       centers: rank 0 sends planes {2,3} -> last writes {nzg-1, nzg};
!                last sends plane nzg-2    -> rank 0 writes plane 1.
!       faces:   rank 0 sends plane 2      -> last writes plane nz;
!                last sends plane nz-1     -> rank 0 writes plane 1.
!   * periodic_pressure_z: rank 0 sends planes {2,3} of rhs_p -> last rank
!       writes its {nzg-1, nzg} (the two planes the gathered solve does not
!       fill on the last slab; single-rank equivalent: k_periodic_z copies).
!
! Rank placement: cudaSetDevice(local rank). Launch with
!   CUDA_VISIBLE_DEVICES=<gpus> mpirun -np <P> ./boundary_layer_cuda -i case
! Single-rank runs work unchanged (every routine below is a no-op or a local
! copy is used instead by the caller); the executable also runs WITHOUT
! mpirun (OpenMPI singleton mode).
!------------------------------------------------------------------------------
module mpi_mod

  use mpi
  use cudafor
  use precision_mod, only: dp

  implicit none
  private

  public :: mpi_initialize, mpi_finalize_run
  public :: myid, nprocs, is_root
  public :: halo_exchange_z, periodic_exchange_z
  public :: allreduce_sum, allreduce_max
  public :: gather_dev, scatter_dev
  public :: alltoall_dev
  public :: pressure_exchange_z
  public :: gather_slabs_host
  public :: allreduce_sum_arr
  public :: gather_planes_r4
  public :: barrier

  integer, public :: STAG_CENTER = 1   ! z-center staggering (U, V, rhs_p)
  integer, public :: STAG_FACE   = 2   ! z-face staggering (W)

  integer :: myid = 0, nprocs = 1
  integer :: ierr

contains

  logical function is_root()
    is_root = (myid == 0)
  end function is_root

  subroutine mpi_initialize()
    integer :: istat, ndev
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    ! one GPU per rank, by local ordinal (single node; CUDA_VISIBLE_DEVICES
    ! selects the physical cards)
    istat = cudaGetDeviceCount(ndev)
    if (nprocs > ndev) then
       if (myid == 0) write(*,*) 'ERROR: ', nprocs, ' ranks but only ', &
                                 ndev, ' visible GPUs'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    istat = cudaSetDevice(mod(myid, ndev))
  end subroutine mpi_initialize

  subroutine mpi_finalize_run()
    call MPI_Finalize(ierr)
  end subroutine mpi_finalize_run

  subroutine barrier()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine barrier

  !---------------------------------------------------------------------------
  ! halo_exchange_z — interior-neighbor ghost update for a device field
  ! F_d(n1, n2, n3loc): plane n3loc-1 goes up into the upper neighbor's
  ! plane 1; plane 2 goes down into the lower neighbor's plane n3loc.
  ! (For both staggering types the LOCAL last index plays the same role;
  ! the caller passes n3loc = nzg for centers, nz for faces.)
  !---------------------------------------------------------------------------
  subroutine halo_exchange_z(F_d, n1, n2, n3loc, tag)
    integer, intent(in) :: n1, n2, n3loc, tag
    real(dp), device, intent(inout) :: F_d(n1, n2, n3loc)
    integer :: istat_mpi(MPI_STATUS_SIZE)
    integer :: up, down, istat

    if (nprocs == 1) return
    istat_mpi = 0
    ! CUDA-aware MPI is NOT stream-aware: kernels still writing F_d must
    ! finish before its device pointer is handed to MPI (the legacy code's
    ! '!$acc update self' implied this synchronization).
    istat = cudaDeviceSynchronize()

    up   = merge(MPI_PROC_NULL, myid+1, myid == nprocs-1)
    down = merge(MPI_PROC_NULL, myid-1, myid == 0)

    ! upward: my plane n3loc-1 -> upper's plane 1
    call MPI_Sendrecv(F_d(1,1,n3loc-1), n1*n2, MPI_REAL8, up,   tag,   &
                      F_d(1,1,1),       n1*n2, MPI_REAL8, down, tag,   &
                      MPI_COMM_WORLD, istat_mpi, ierr)
    ! downward: my plane 2 -> lower's plane n3loc
    call MPI_Sendrecv(F_d(1,1,2),       n1*n2, MPI_REAL8, down, tag+100, &
                      F_d(1,1,n3loc),   n1*n2, MPI_REAL8, up,   tag+100, &
                      MPI_COMM_WORLD, istat_mpi, ierr)
  end subroutine halo_exchange_z

  !---------------------------------------------------------------------------
  ! periodic_exchange_z — periodic wrap between rank 0 and the last rank.
  ! Centers (stag = STAG_CENTER, n3loc = nzg): rank 0 sends its planes {2,3}
  ! (contiguous pair) into the last rank's {n3loc-1, n3loc}; the last rank
  ! sends its plane n3loc-2 into rank 0's plane 1.
  ! Faces (stag = STAG_FACE, n3loc = nz): single planes, 2 -> last's n3loc
  ! and n3loc-1 -> rank 0's plane 1.
  !---------------------------------------------------------------------------
  subroutine periodic_exchange_z(F_d, n1, n2, n3loc, stag, tag)
    integer, intent(in) :: n1, n2, n3loc, stag, tag
    real(dp), device, intent(inout) :: F_d(n1, n2, n3loc)
    integer :: istat_mpi(MPI_STATUS_SIZE)
    integer :: last, istat

    if (nprocs == 1) return
    istat_mpi = 0
    istat = cudaDeviceSynchronize()   ! see halo_exchange_z
    last = nprocs - 1

    if (stag == STAG_CENTER) then
       if (myid == 0) then
          call MPI_Sendrecv(F_d(1,1,2),        2*n1*n2, MPI_REAL8, last, tag,   &
                            F_d(1,1,1),          n1*n2, MPI_REAL8, last, tag+1, &
                            MPI_COMM_WORLD, istat_mpi, ierr)
       else if (myid == last) then
          call MPI_Sendrecv(F_d(1,1,n3loc-2),    n1*n2, MPI_REAL8, 0, tag+1, &
                            F_d(1,1,n3loc-1),  2*n1*n2, MPI_REAL8, 0, tag,   &
                            MPI_COMM_WORLD, istat_mpi, ierr)
       end if
    else
       if (myid == 0) then
          call MPI_Sendrecv(F_d(1,1,2),          n1*n2, MPI_REAL8, last, tag,   &
                            F_d(1,1,1),          n1*n2, MPI_REAL8, last, tag+1, &
                            MPI_COMM_WORLD, istat_mpi, ierr)
       else if (myid == last) then
          call MPI_Sendrecv(F_d(1,1,n3loc-1),    n1*n2, MPI_REAL8, 0, tag+1, &
                            F_d(1,1,n3loc),      n1*n2, MPI_REAL8, 0, tag,   &
                            MPI_COMM_WORLD, istat_mpi, ierr)
       end if
    end if
  end subroutine periodic_exchange_z

  !---------------------------------------------------------------------------
  ! Host-scalar reductions (legacy MPI_Allreduce equivalents).
  !---------------------------------------------------------------------------
  subroutine allreduce_sum(v)
    real(dp), intent(inout) :: v
    real(dp) :: r
    if (nprocs == 1) return
    call MPI_Allreduce(v, r, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    v = r
  end subroutine allreduce_sum

  subroutine allreduce_max(v)
    real(dp), intent(inout) :: v
    real(dp) :: r
    if (nprocs == 1) return
    call MPI_Allreduce(v, r, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    v = r
  end subroutine allreduce_max

  !---------------------------------------------------------------------------
  ! gather_dev / scatter_dev — fixed-count root gather/scatter of DEVICE
  ! buffers (CUDA-aware). Used by the P5.1 Poisson gather: every rank owns
  ! an equal-sized compact chunk (nxp*ns*nm reals).
  !---------------------------------------------------------------------------
  subroutine gather_dev(send_d, count, recv_d)
    integer, intent(in) :: count
    real(dp), device, intent(in)    :: send_d(count)
    real(dp), device, intent(inout) :: recv_d(count*nprocs)
    integer :: istat
    istat = cudaDeviceSynchronize()   ! see halo_exchange_z
    call MPI_Gather(send_d, count, MPI_REAL8, recv_d, count, MPI_REAL8, &
                    0, MPI_COMM_WORLD, ierr)
  end subroutine gather_dev

  subroutine scatter_dev(send_d, count, recv_d)
    integer, intent(in) :: count
    real(dp), device, intent(in)    :: send_d(count*nprocs)
    real(dp), device, intent(inout) :: recv_d(count)
    integer :: istat
    istat = cudaDeviceSynchronize()   ! see halo_exchange_z
    call MPI_Scatter(send_d, count, MPI_REAL8, recv_d, count, MPI_REAL8, &
                     0, MPI_COMM_WORLD, ierr)
  end subroutine scatter_dev

  !---------------------------------------------------------------------------
  ! alltoall_dev — device-buffer MPI_Alltoall with fixed per-rank block size
  ! (the P5.2 Poisson transposes). Blocks are packed by destination rank in
  ! the send buffer and land ordered by source rank in the receive buffer.
  !---------------------------------------------------------------------------
  subroutine alltoall_dev(send_d, count, recv_d)
    integer, intent(in) :: count
    real(dp), device, intent(in)    :: send_d(count*nprocs)
    real(dp), device, intent(inout) :: recv_d(count*nprocs)
    integer :: istat
    istat = cudaDeviceSynchronize()   ! see halo_exchange_z
    call MPI_Alltoall(send_d, count, MPI_REAL8, recv_d, count, MPI_REAL8, &
                      MPI_COMM_WORLD, ierr)
  end subroutine alltoall_dev

  !---------------------------------------------------------------------------
  ! pressure_exchange_z — ghost planes of the pseudo-pressure after the
  ! gathered solve fills each rank's OWNED planes (buffer planes 1..ns,
  ! i.e. actual z-indices 2..ns+1 of rhs_p(2:, 2:, 2:nzg)). The buffer
  ! passed here is rhs_p viewed as nplanes = nzg-1 planes (actual 2..nzg):
  !   * interior: each rank sends buffer plane 1 DOWN; receives from UP
  !     into buffer plane nplanes  (legacy update_ghost_interior_planes_
  !     pressure: neighbor's first owned plane -> my top ghost).
  !   * periodic: rank 0 sends buffer planes {1,2} (contiguous) -> last
  !     rank's {nplanes-1, nplanes} (legacy apply_periodic_z_pressure).
  !---------------------------------------------------------------------------
  subroutine pressure_exchange_z(F_d, n1, n2, nplanes, tag)
    integer, intent(in) :: n1, n2, nplanes, tag
    real(dp), device, intent(inout) :: F_d(n1, n2, nplanes)
    integer :: istat_mpi(MPI_STATUS_SIZE)
    integer :: up, down, last, istat

    if (nprocs == 1) return
    istat_mpi = 0
    istat = cudaDeviceSynchronize()   ! see halo_exchange_z
    last = nprocs - 1
    up   = merge(MPI_PROC_NULL, myid+1, myid == last)
    down = merge(MPI_PROC_NULL, myid-1, myid == 0)

    call MPI_Sendrecv(F_d(1,1,1),       n1*n2, MPI_REAL8, down, tag,   &
                      F_d(1,1,nplanes), n1*n2, MPI_REAL8, up,   tag,   &
                      MPI_COMM_WORLD, istat_mpi, ierr)

    if (myid == 0) then
       call MPI_Send(F_d(1,1,1), 2*n1*n2, MPI_REAL8, last, tag+1, &
                     MPI_COMM_WORLD, ierr)
    else if (myid == last) then
       call MPI_Recv(F_d(1,1,nplanes-1), 2*n1*n2, MPI_REAL8, 0, tag+1, &
                     MPI_COMM_WORLD, istat_mpi, ierr)
    end if
  end subroutine pressure_exchange_z

  !---------------------------------------------------------------------------
  ! gather_slabs_host — assemble a global field on rank 0 from per-rank HOST
  ! slabs. Rank r's slab covers global planes start_r .. start_r+n3loc_r-1;
  ! slabs overlap by construction and overlapping planes hold identical
  ! values (ghosts are exchanged), so later writes simply overwrite.
  !---------------------------------------------------------------------------
  subroutine gather_slabs_host(loc, n1, n2, n3loc, start, glob, n3glob)
    integer,  intent(in)    :: n1, n2, n3loc, start, n3glob
    real(dp), intent(in)    :: loc(n1, n2, n3loc)
    real(dp), intent(inout) :: glob(n1, n2, n3glob)
    integer :: istat_mpi(MPI_STATUS_SIZE)
    integer :: r, meta(2)
    real(dp), allocatable :: buf(:,:,:)

    if (nprocs == 1) then
       glob(:,:,start:start+n3loc-1) = loc
       return
    end if
    istat_mpi = 0

    if (myid == 0) then
       glob(:,:,start:start+n3loc-1) = loc
       do r = 1, nprocs-1
          call MPI_Recv(meta, 2, MPI_INTEGER, r, 900, MPI_COMM_WORLD, &
                        istat_mpi, ierr)
          allocate (buf(n1, n2, meta(2)))
          call MPI_Recv(buf, n1*n2*meta(2), MPI_REAL8, r, 901, &
                        MPI_COMM_WORLD, istat_mpi, ierr)
          glob(:,:,meta(1):meta(1)+meta(2)-1) = buf
          deallocate (buf)
       end do
    else
       meta = [start, n3loc]
       call MPI_Send(meta, 2, MPI_INTEGER, 0, 900, MPI_COMM_WORLD, ierr)
       call MPI_Send(loc, n1*n2*n3loc, MPI_REAL8, 0, 901, MPI_COMM_WORLD, ierr)
    end if
  end subroutine gather_slabs_host

  !---------------------------------------------------------------------------
  ! allreduce_sum_arr — in-place elementwise sum-reduce of a host array
  ! (the P5.3 statistics: per-rank partial z-sums -> global sums).
  !---------------------------------------------------------------------------
  subroutine allreduce_sum_arr(a, n)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: a(n)
    if (nprocs == 1) return
    call MPI_Allreduce(MPI_IN_PLACE, a, n, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
  end subroutine allreduce_sum_arr

  !---------------------------------------------------------------------------
  ! gather_planes_r4 — assemble float32 z-plane windows on rank 0 (the P5.3
  ! box output). Rank r's contribution covers global planes gstart_r ..
  ! gstart_r+nploc_r-1 of glob's third dimension; overlapping planes hold
  ! identical values (ghost-synced fields) and are simply overwritten.
  !---------------------------------------------------------------------------
  subroutine gather_planes_r4(loc4, n1, n2, nploc, gstart, glob4, n3glob)
    integer, intent(in)    :: n1, n2, nploc, gstart, n3glob
    real(4), intent(in)    :: loc4(n1, n2, nploc)
    real(4), intent(inout) :: glob4(n1, n2, n3glob)
    integer :: istat_mpi(MPI_STATUS_SIZE)
    integer :: r, meta(2)
    real(4), allocatable :: buf(:,:,:)

    if (nprocs == 1) then
       if (nploc > 0) glob4(:,:,gstart:gstart+nploc-1) = loc4
       return
    end if
    istat_mpi = 0

    if (myid == 0) then
       if (nploc > 0) glob4(:,:,gstart:gstart+nploc-1) = loc4
       do r = 1, nprocs-1
          call MPI_Recv(meta, 2, MPI_INTEGER, r, 910, MPI_COMM_WORLD, &
                        istat_mpi, ierr)
          if (meta(2) > 0) then
             allocate (buf(n1, n2, meta(2)))
             call MPI_Recv(buf, n1*n2*meta(2), MPI_REAL4, r, 911, &
                           MPI_COMM_WORLD, istat_mpi, ierr)
             glob4(:,:,meta(1):meta(1)+meta(2)-1) = buf
             deallocate (buf)
          end if
       end do
    else
       meta = [gstart, nploc]
       call MPI_Send(meta, 2, MPI_INTEGER, 0, 910, MPI_COMM_WORLD, ierr)
       if (nploc > 0) &
          call MPI_Send(loc4, n1*n2*nploc, MPI_REAL4, 0, 911, &
                        MPI_COMM_WORLD, ierr)
    end if
  end subroutine gather_planes_r4

end module mpi_mod
