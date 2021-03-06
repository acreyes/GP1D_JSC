subroutine soln_WENO(dt, V)

#include "definition.h"

  use grid_data
  use sim_data
  use eigensystem

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  integer :: i, s, var, varL, varR
  real :: RHS, delV, dot
  logical :: conservative
  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig0, leig0
  real, dimension(NSYS_VAR, 3) :: wlp, wlm, beta, wbarp, wbarm
  real, dimension(3          ) :: gamp, gamm, Plp, Plm
  real, dimension(3       , 3) :: coeffm, coeffp
  real, dimension(5          ) :: Stencil
  real, dimension(NSYS_VAR   ) :: vL, vR, vim2, vim1, vi, vip1, vip2, wsumm, wsump, C0, C1, C2, vecL, vecR


  !initialize w/ zeros
  vim2 = 0.; vim1 = 0.
  vi   = 0.
  vip1 = 0.; vip2 = 0.
  vL   = 0.; vR   = 0.
  
  wsumm = 0.; wsump = 0.
  wlp   = 0.; wlm   = 0.
  wbarp = 0.; wbarm = 0.
  beta  = 0.

  conservative = .false.
  
  gamp(1) = .1
  gamp(2) = .6
  gamp(3) = .3

  gamm(1) = .3
  gamm(2) = .6
  gamm(3) = .1

  coeffp(1, 1:3) = (/ 2., -7., 11. /)
  coeffp(2, 1:3) = (/-1.,  5., 2.  /)
  coeffp(3, 1:3) = (/ 2.,  5.,-1.  /)

  coeffm(1, 1:3) = (/-1.,  5., 2. /)
  coeffm(2, 1:3) = (/ 2.,  5.,-1. /)
  coeffm(3, 1:3) = (/11., -7., 2. /)

  coeffp = coeffp/6.
  coeffm = coeffm/6.
  

  do i = gr_ibeg-1, gr_iend+1
     !copy cell-centered values to left and right states
     gr_vL(DENS_VAR:NUMB_VAR,i) = V(DENS_VAR:NUMB_VAR,i)
     gr_vR(DENS_VAR:NUMB_VAR,i) = V(DENS_VAR:NUMB_VAR,i)

     !get eigen information
     call eigenvalues(V(DENS_VAR:GAME_VAR,i), lambda)
     call mhd_eigenvectors(V(DENS_VAR:GAME_VAR,i), conservative, reig0, leig0)

     if (sim_charLimiting) then
        !do WENO reconstruction on char variables
        do var = 1, NUMB_WAVE
           vim2(var) = dot_product(V(DENS_VAR:PRES_VAR, i-2), leig0(DENS_VAR:PRES_VAR, var))
           vim1(var) = dot_product(V(DENS_VAR:PRES_VAR, i-1), leig0(DENS_VAR:PRES_VAR, var))
           vi(  var) = dot_product(V(DENS_VAR:PRES_VAR, i  ), leig0(DENS_VAR:PRES_VAR, var))
           vip1(var) = dot_product(V(DENS_VAR:PRES_VAR, i+1), leig0(DENS_VAR:PRES_VAR, var))
           vip2(var) = dot_product(V(DENS_VAR:PRES_VAR, i+2), leig0(DENS_VAR:PRES_VAR, var))
        end do
     else
        !do WENO reconstruction on primitive variables
        vim2(DENS_VAR:PRES_VAR) = V(DENS_VAR:PRES_VAR, i-2)
        vim1(DENS_VAR:PRES_VAR) = V(DENS_VAR:PRES_VAR, i-1)
        vi(  DENS_VAR:PRES_VAR) = V(DENS_VAR:PRES_VAR, i  )
        vip1(DENS_VAR:PRES_VAR) = V(DENS_VAR:PRES_VAR, i+1)
        vip2(DENS_VAR:PRES_VAR) = V(DENS_VAR:PRES_VAR, i+2)
     end if

     !Calculate the smoothness indicators, beta
     beta(DENS_VAR:PRES_VAR,1) = 13./12.*(vim2(DENS_VAR:PRES_VAR) - 2.*vim1(DENS_VAR:PRES_VAR) &
          + vi(DENS_VAR:PRES_VAR))**2 &
          + .25*(vim2(DENS_VAR:PRES_VAR) - 4.*vim1(DENS_VAR:PRES_VAR) + 3.*vi(DENS_VAR:PRES_VAR))**2
     
     beta(DENS_VAR:PRES_VAR,2) = 13./12.*(vim1(DENS_VAR:PRES_VAR) - 2.*vi(DENS_VAR:PRES_VAR) &
          + vip1(DENS_VAR:PRES_VAR))**2  &
          + .25*(vim1(DENS_VAR:PRES_VAR) - vip1(DENS_VAR:PRES_VAR))**2
     
     beta(DENS_VAR:PRES_VAR,3) = 13./12.*(vi(DENS_VAR:PRES_VAR) - 2.*vip1(DENS_VAR:PRES_VAR) &
          + vip2(DENS_VAR:PRES_VAR))**2 &
          + .25*(3.*vi(DENS_VAR:PRES_VAR) - 4.*vip1(DENS_VAR:PRES_VAR) + vip2(DENS_VAR:PRES_VAR))**2

     !calculate weight factors
     do s = 1,3
        wbarp(DENS_VAR:PRES_VAR, s) = gamp(s)/((sim_WENeps+beta(DENS_VAR:PRES_VAR, s))**sim_mval)
        wbarm(DENS_VAR:PRES_VAR, s) = gamm(s)/((sim_WENeps+beta(DENS_VAR:PRES_VAR, s))**sim_mval)
     end do

     wsumm(DENS_VAR:PRES_VAR) = wbarm(DENS_VAR:PRES_VAR,1) + wbarm(DENS_VAR:PRES_VAR,2) &
          + wbarm(DENS_VAR:PRES_VAR,3)
     wsump(DENS_VAR:PRES_VAR) = wbarp(DENS_VAR:PRES_VAR, 1) + wbarp(DENS_VAR:PRES_VAR, 2) &
          + wbarp(DENS_VAR:PRES_VAR, 3)
     do s = 1,3
        wlp(DENS_VAR:PRES_VAR, s) = wbarp(DENS_VAR:PRES_VAR,s)/wsump
        wlm(DENS_VAR:PRES_VAR, s) = wbarm(DENS_VAR:PRES_VAR,s)/wsumm
     end do

      !calculate initial left and right states
     do var = DENS_VAR, PRES_VAR
        !make stencil for weno
        stencil = (/ vim2(var), vim1(var), vi(var), vip1(var), vip2(var) /)
        !calculate weno polynomials at interfaces
        Plp(1) = dot_product(stencil(1:3), coeffp(1, 1:3))
        Plp(2) = dot_product(stencil(2:4), coeffp(2, 1:3))
        Plp(3) = dot_product(stencil(3:5), coeffp(3, 1:3))

        Plm(1) = dot_product(stencil(1:3), coeffm(1, 1:3))
        Plm(2) = dot_product(stencil(2:4), coeffm(2, 1:3))
        Plm(3) = dot_product(stencil(3:5), coeffm(3, 1:3))

        vecL(var) = dot_product(wlm(var, 1:3), Plm(1:3))
        vecR(var) = dot_product(wlp(var, 1:3), Plp(1:3))
        
        
     end do
     do var = DENS_VAR, PRES_VAR
        if (sim_charLimiting) then
           !project char vars back onto primitive vars
           vL(var) = dot_product(reig0(var,1:NUMB_WAVE),vecL(1:NUMB_WAVE))
           vR(var) = dot_product(reig0(var,1:NUMB_WAVE),vecR(1:NUMB_WAVE))
        else
           vL(var) = vecL(var)
           vR(var) = vecR(var)
        end if
     end do

     gr_vL(DENS_VAR:PRES_VAR, i) = vL(DENS_VAR:PRES_VAR)
     gr_vR(DENS_VAR:PRES_VAR, i) = vR(DENS_VAR:PRES_VAR)


  end do
end subroutine soln_WENO
