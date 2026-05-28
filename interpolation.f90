!------------------------------------------!
!        Module for interpolation          !
!------------------------------------------!
Module interpolation

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : weight_y_0, weight_y_1

  ! prevent implicit typing
  Implicit None

Contains

  !-------------------------------------------------------!
  !            Linear interpolation of u in x             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  ! Parameter: di-> 1=faces to center  (normal   average) !
  !                 2=centers to faces (weighted average) !
  !-------------------------------------------------------!
  ! uniform mesh assumed
  Subroutine interpolate_x(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer(Int32) :: n(3), n1, n2, n3

    n  = shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    !$acc kernels default(present)
    ui = 0d0
    ui( 1:n1-1, 1:n2, 1:n3) = 0.5d0*( u(1:n1-1,:,:) + u(2:n1,:,:) )
    !$acc end kernels

  End Subroutine interpolate_x

  !-------------------------------------------------------!
  !            Linear interpolation of u in y             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  ! Parameter: di-> 1=faces to center  (normal   average) !
  !                 2=centers to faces (weighted average) !
  !-------------------------------------------------------!
  Subroutine interpolate_y(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer (Int32) :: n(3), n1, n2, n3, i1, i3

    n  = shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    !$acc kernels default(present)
    ui = 0d0
    !$acc end kernels
    If ( di==1 ) Then
       ! faces to centers (normal average)
       !$acc kernels default(present)
       ui(1:n1, 1:n2-1, 1:n3) = 0.5d0*( u(:,1:n2-1,:) + u(:,2:n2,:) )
       !$acc end kernels
    Elseif ( di==2 ) Then
       ! centers to faces (weighted average)
       !$acc parallel loop collapse(2) default(present)
       Do i1 = 1, n1
          Do i3 = 1, n3
             ui(i1, 1:n2-1, i3) = weight_y_0*u(i1,1:n2-1,i3) + weight_y_1*u(i1,2:n2,i3)
          End Do
       End Do
    Else
       Stop 'Error: invalid interpolation'
    End If

  End Subroutine interpolate_y

  !-------------------------------------------------------!
  !            Linear interpolation of u in z             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  ! Parameter: di-> 1=faces to center  (normal   average) !
  !                 2=centers to faces (weighted average) !
  !-------------------------------------------------------!
  ! uniform mesh assumed
  Subroutine interpolate_z(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer(Int32) :: n(3), n1, n2, n3

    n  = shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    !$acc kernels default(present)
    ui = 0d0
    ui(1:n1, 1:n2, 1:n3-1) = 0.5d0*( u(:,:,1:n3-1) + u(:,:,2:n3) )
    !$acc end kernels

  End Subroutine interpolate_z


  !-------------------------------------------------------!
  !        General second order interpolation in y        !
  !                                                       !
  ! Input : y,u,yi                                        !
  ! Output: ui                                            !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine interpolate_y_2nd(y,u,yi,ui)

    Real    (Int64), Intent(In)  :: u (:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Real    (Int64), Intent(In)  :: y(:), yi(:)

    Integer(Int32) :: n(3), n1, n2, n3, j, nn(1)
    Real   (Int64) :: w0, w1, w2, yref

    n  = shape(u)
    n1 = n(1)
    n3 = n(3)

    nn = shape(yi)
    n2 = nn(1)

    !$acc kernels default(present)
    ui = 0d0
    !$acc end kernels
    ! middle points
    !$acc parallel loop default(present) private(yref, w0, w1, w2)
    Do j = 2, n2-1
       yref = yi(j)
       w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
       w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j)   - y(j-1) )*( y(j)   - y(j+1) ) )
       w2   = ( yref - y(j-1) )*( yref - y(j)   )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
       ui(1:n1, j, 1:n3) = w0*u(1:n1,j-1,1:n3) + w1*u(1:n1,j,1:n3) +  w2*u(1:n1,j+1,1:n3)
    End Do

    ! first point
    j    = 2
    yref = yi(1)
    w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
    w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j  ) - y(j-1) )*( y(j)   - y(j+1) ) )
    w2   = ( yref - y(j-1) )*( yref - y(j  ) )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
    !$acc kernels default(present)
    ui(1:n1, 1, 1:n3) = w0*u(1:n1,j-1,1:n3) + w1*u(1:n1,j,1:n3) +  w2*u(1:n1,j+1,1:n3)
    !$acc end kernels

    ! last point
    j    = n2-1
    yref = yi(n2)
    w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
    w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j  ) - y(j-1) )*( y(j)   - y(j+1) ) )
    w2   = ( yref - y(j-1) )*( yref - y(j  ) )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
    !$acc kernels default(present)
    ui(1:n1, n2, 1:n3) = w0*u(1:n1,j-1,1:n3) + w1*u(1:n1,j,1:n3) +  w2*u(1:n1,j+1,1:n3)
    !$acc end kernels

  End Subroutine interpolate_y_2nd

  !-------------------------------------------------------!
  !        General second order interpolation in x        !
  !                                                       !
  ! Input : y,u,yi                                        !
  ! Output: ui                                            !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine interpolate_x_2nd(y,u,yi,ui)

    Real    (Int64), Intent(In)  :: u (:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Real    (Int64), Intent(In)  :: y(:), yi(:)

    Integer(Int32) :: n(3), n1, n2, n3, j, nn(1)
    Real   (Int64) :: w0, w1, w2, yref

    n  = shape(u)
    n2 = n(2)
    n3 = n(3)

    nn = shape(yi)
    n1 = nn(1)

    !$acc kernels default(present)
    ui = 0d0
    !$acc end kernels
    ! middle points
    !$acc parallel loop default(present) private(yref, w0, w1, w2)
    Do j = 2, n1-1
       yref = yi(j)
       w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
       w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j)   - y(j-1) )*( y(j)   - y(j+1) ) )
       w2   = ( yref - y(j-1) )*( yref - y(j)   )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
       ui(j, 1:n2, 1:n3) = w0*u(j-1,1:n2,1:n3) + w1*u(j,1:n2,1:n3) +  w2*u(j+1,1:n2,1:n3)
    End Do

    ! first point
    j    = 2
    yref = yi(1)
    w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
    w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j  ) - y(j-1) )*( y(j)   - y(j+1) ) )
    w2   = ( yref - y(j-1) )*( yref - y(j  ) )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
    !$acc kernels default(present)
    ui(1, 1:n2, 1:n3) = w0*u(j-1,1:n2,1:n3) + w1*u(j,1:n2,1:n3) +  w2*u(j+1,1:n2,1:n3)
    !$acc end kernels

    ! last point
    j    = n1-1
    yref = yi(n1)
    w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
    w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j  ) - y(j-1) )*( y(j)   - y(j+1) ) )
    w2   = ( yref - y(j-1) )*( yref - y(j  ) )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
    !$acc kernels default(present)
    ui(n1, 1:n2, 1:n3) = w0*u(j-1,1:n2,1:n3) + w1*u(j,1:n2,1:n3) +  w2*u(j+1,1:n2,1:n3)
    !$acc end kernels

  End Subroutine interpolate_x_2nd

  !-------------------------------------------------------!
  !        General second order interpolation in z        !
  !                                                       !
  ! Input : y,u,yi                                        !
  ! Output: ui                                            !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine interpolate_z_2nd(y,u,yi,ui)

    Real    (Int64), Intent(In)  :: u (:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Real    (Int64), Intent(In)  :: y(:), yi(:)

    Integer(Int32) :: n(3), n1, n2, n3, j, nn(1)
    Real   (Int64) :: w0, w1, w2, yref

    n  = shape(u)
    n1 = n(1)
    n2 = n(2)

    nn = shape(yi)
    n3 = nn(1)

    !$acc kernels default(present)
    ui = 0d0
    !$acc end kernels
    ! middle points
    !$acc parallel loop default(present) private(yref, w0, w1, w2)
    Do j = 2, n3-1
       yref = yi(j)
       w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
       w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j)   - y(j-1) )*( y(j)   - y(j+1) ) )
       w2   = ( yref - y(j-1) )*( yref - y(j)   )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
       ui(1:n1, 1:n2, j) = w0*u(1:n1,1:n2,j-1) + w1*u(1:n1,1:n2,j) + w2*u(1:n1,1:n2,j+1)
    End Do

    ! first point
    j    = 2
    yref = yi(1)
    w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
    w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j  ) - y(j-1) )*( y(j)   - y(j+1) ) )
    w2   = ( yref - y(j-1) )*( yref - y(j  ) )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
    !$acc kernels default(present)
    ui(1:n1, 1:n2, 1) = w0*u(1:n1,1:n2,j-1) + w1*u(1:n1,1:n2,j) + w2*u(1:n1,1:n2,j+1)
    !$acc end kernels

    ! last point
    j    = n3-1
    yref = yi(n3)
    w0   = ( yref - y(j)   )*( yref - y(j+1) )/( ( y(j-1) - y(j)   )*( y(j-1) - y(j+1) ) )
    w1   = ( yref - y(j-1) )*( yref - y(j+1) )/( ( y(j  ) - y(j-1) )*( y(j)   - y(j+1) ) )
    w2   = ( yref - y(j-1) )*( yref - y(j  ) )/( ( y(j+1) - y(j-1) )*( y(j+1) - y(j)   ) )
    !$acc kernels default(present)
    ui( 1:n1, 1:n2, n3) = w0*u(1:n1,1:n2,j-1) + w1*u(1:n1,1:n2,j) + w2*u(1:n1,1:n2,j+1)
    !$acc end kernels

  End Subroutine interpolate_z_2nd

End Module interpolation
