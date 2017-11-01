subroutine gp_eigens()
  !this needs to calculate the eigensystem for each of the GPWENO stencils
  !for use in calculating the smoothness indicators

#include "definition.h"

  use sim_data
  use grid_data

  implicit none

  real, allocatable, dimension(:,:) :: C
  integer :: N, LDA, LWORK, INFO, i, j

  real, allocatable, dimension(:) :: W, WORK

  !recons on cell avgs
  N = gr_radius + 1

  LDA = N
  LWORK = 66*N
  allocate(C(N,N))
  allocate(W(N))
  allocate(WORK(LWORK))

  do i = 1, N
     do j = 1, N
        C(i,j) = SE(REAL(i),REAL(j))
     end do
  end do


  call DSYEV('V', 'L', N, C, LDA, W, WORK, LWORK, INFO)
  if (INFO < 0) then
     print *, "error in dsyev, info != 1"
     stop
  end if
  !allocate eigensystem vars for GP
  allocate(sim_gp_evals(N  )); sim_gp_evals = W
  allocate(sim_gp_evecs(N,N)); sim_gp_evecs = C
  
  deallocate(C)
  deallocate(W)
  deallocate(WORK)

  
  return

contains

function RQ(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r, alpha
    alpha = 1.e-8
    r = abs(x-y)
    f = (1. + r**2/(2.*sim_RQ_alpha*sim_sigdel**2))**(-sim_RQ_alpha)
    return
  end function RQ
  
function quad_exact(x1,x2) result(Integ)
  use sim_data, only: sim_sigdel
  implicit none
  !exact quadrature, only good for SE kernel
  real, intent(IN) :: x1, x2

  real :: Integ, yxp, yxn, yxm, sigdel
  sigdel = sim_sigdel*SQRT(2.)
  !sigdel = 4.*(2.*REAL(gr_radius)+1)

  yxp = (x1 - x2 + 1.)/sigdel
  yxn = (x1      -x2)/sigdel
  yxm = (x1 - x2 -1.)/sigdel


  Integ = 0.5*SQRT(PI)*(sigdel)**2 *( yxp*ERF(yxp) + yxm*ERF(yxm) &
       - 2.*( yxn*ERF(yxn) + 1./SQRT(PI) *EXP(-yxn**2) ) &
       + 1./SQRT(PI) * ( EXP(-yxp**2) + exp(-yxm**2) ) )
  return
end function quad_exact

function SE(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r
    r = abs(x-y)
    f = EXP( -0.5*(r/sim_sigdel)**2 )
    return
  end function SE

  function mat_3h(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r, rt3li
    r = abs(x-y)
    rt3li = SQRT(3.)/sim_sigdel
    f = (1 + rt3li*r)*EXP(-rt3li*r)
    return
  end function mat_3h

  function mat_5h(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r, rt5li
    r = abs(x-y)
    rt5li = SQRT(5.)/sim_sigdel
    f = (1 + rt5li*r + 5./3. * (r/sim_sigdel)**2)*EXP(-rt5li*r)
    return
  end function mat_5h

  function mat_7h(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r, rt7li
    r = abs(x-y)
    rt7li = SQRT(7.)/sim_sigdel

    f = ( 1 + rt7li*r + 14./5.*(r/sim_sigdel)**2 + 7.*SQRT(7.)/15.*(r/sim_sigdel)**3 )*EXP(-rt7li*r)

    return
  end function mat_7h

end subroutine gp_eigens
