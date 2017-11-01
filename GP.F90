module GP
#include "definition.h"

  use grid_data
  use sim_data
  use linalg

  


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Quadrature Rules for Integrated Kernel!!!!!!

  function intg_kernel(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f
    f = quad_exact(x, y)
    return
  end function intg_kernel

  function quad_cross(x, t) result(f)
    implicit none
    real, intent(IN) :: x, t
    real :: f
    f = 0.
    f = 0.5*sim_sigdel*SQRT(2.*PI)*int_egrand(x,t)
    return
  end function quad_cross

  function intg_predvec(x) result(T)
    implicit none
    real, intent(IN) :: x
    real, dimension(2) :: T

    T(1) = quad_cross(x, -0.5)
    T(2) = quad_cross(x,  0.5)
    return
  end function intg_predvec

  

 
  function SE(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r
    r = abs(x-y)
    f = EXP( -0.5*(r/sim_sigdel)**2 )
    return
  end function SE

  
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! SE Kernel !!!!!!!!!!!!!!!!!!!

  function int_egrand(x, t) result(f)
    implicit none
    real, intent(IN) :: x, t
    real :: f, sigdel
    sigdel = sim_sigdel*SQRT(2.)

    f = ERF( (x + .5 - t)/sigdel) - ERF( (x - .5 - t)/sigdel )
    !f = ERF( (x + .5*gr_dx - t)/sigdel) - ERF( (x - .5*gr_dx - t)/sigdel )

  end function int_egrand

  function quad_exact(x1,x2) result(Integ)
    
    !exact quadrature, only good for SE kernel
    real, intent(IN) :: x1, x2

    real :: Integ, yxp, yxn, yxm, sigdel
    sigdel = sim_sigdel*SQRT(2.)

    yxp = (x1 - x2 + 1.)/sigdel
    yxn = (x1      -x2)/sigdel
    yxm = (x1 - x2 -1.)/sigdel
    
    
    Integ = 0.5*SQRT(PI)*(sigdel)**2 *( yxp*ERF(yxp) + yxm*ERF(yxm) &
         - 2.*( yxn*ERF(yxn) + 1./SQRT(PI) *EXP(-yxn**2) ) &
         + 1./SQRT(PI) * ( EXP(-yxp**2) + exp(-yxm**2) ) )
    return
  end function quad_exact

  

  function cross_cor(x) result(T)
    !returns the cross-correlation between the left and right states and the cell centered at x
    !see eq 24
    implicit none
    real, intent(IN) :: x
    real, dimension(2) :: T
    real :: sigdel
    sigdel = sim_sigdel*SQRT(2.)
    !sim_sigdel = sim_sigma/gr_dx
    T(1) = int_egrand(x, -.5)
    T(2) = int_egrand(x, 0.5)
    T = T*.5*sigdel*SQRT(PI)

  end function cross_cor

  function K(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f

    f = 0.
    f = SE(x,y)
    return
  end function K

  
end module GP
