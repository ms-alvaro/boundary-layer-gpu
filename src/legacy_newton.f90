!-----------------------------------------------!
!                                               !
!-----------------------------------------------!
Module Newton_solver

  ! module parameters:
  Implicit none
  Integer,      Parameter :: maxiter = 30
  Real(kind=8), Parameter :: tol     = 1.d-10

Contains

  Subroutine Newton_iter_solver(f, fp, x0, x, iters, debug)
    
    ! Estimate the zero of f(x) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!) 
    !   the number of iterations iters
     
    Implicit none
    
    Real(kind=8), Intent(in)  :: x0
    Real(kind=8), External    :: f, fp
    Logical,      Intent(in)  :: debug
    Real(kind=8), Intent(out) :: x
    Integer,      Intent(out) :: iters

    ! Declare any local variables:
    Real(kind=8) :: deltax, fx, fxprime
    Integer      :: k

    ! initial guess
    x = x0

    If (debug) Then      
       Print 11, x
11     Format('Initial guess: x = ', e22.15)
    Endif

    ! Newton iteration to find a zero of f(x) 
    Do k=1,maxiter

        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)

        If ( abs(fx) < tol ) then
           Exit  ! jump out of do loop
        Endif
        
        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

        If (debug) Then
           Print 12, k,x
12         Format('After', i3, ' iterations, x = ', e22.15)
        Endif
     Enddo
     
     If (k > maxiter .and. .False.) Then
        ! might not have converged        
        fx = f(x)
        If ( Abs(fx) > tol ) then
           Print *, 'Warning: utau_model not converged'
        Endif
     Endif
     
     ! number of iterations taken:
     iters = k-1
     
   End subroutine Newton_iter_solver
   
 End module Newton_solver
 
 !-----------------------------------------------!
 !                                               !
 !-----------------------------------------------!
 Module functions_wallmodel
   
   ! Use 
   Use global, Only : Umean_model, kappa_model, B_model, yg_model, nu
   
   Implicit None
   
 Contains
   
   !---------------------------------------------------------!
   ! F = utau/kappa*ln(y(jref)*utau/nu) + B*utau - Umean = 0 !
   !---------------------------------------------------------!
   Real(kind=8) function f_law_of_wall(utau_model)
     Implicit none
     Real(kind=8), Intent(in) :: utau_model
     
     f_law_of_wall = utau_model/kappa_model*dlog( yg_model*utau_model/nu ) + B_model*utau_model - Umean_model
     
   End function f_law_of_wall
   
   !------------------------------------------------------!
   ! dF = 1/kappa*ln(y*utau/nu) + 1/kappa + B             !
   !------------------------------------------------------!  
   Real(kind=8) function df_law_of_wall(utau_model)
     Implicit none
     Real(kind=8), Intent(in) :: utau_model
     
     df_law_of_wall = 1d0/kappa_model*dlog(yg_model*utau_model/nu) + 1d0/kappa_model + B_model
     
   End function df_law_of_wall
   
End module functions_wallmodel
 
