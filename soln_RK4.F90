subroutine soln_RK4(dt, j, Vj, Flux)
  !performs the jth step of the RK4 algorithm
  !Total Flux = (k1 + 2k2 + 2k3 + k4)/6
#include "definition.h"
  
  use grid_data
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: j
  real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: Vj
  real, dimension(NSYS_VAR,gr_imax), intent(INOUT) :: Flux

  real :: A, F, dtx
  integer :: i
  real, dimension(NUMB_VAR,gr_imax) :: Uj

  if (j == 1 .OR. j == 2) then
     A = 0.5
  else 
     A = 1.
  end if

  if (j == 1 .OR. j == 4) then
     F = 1./6.
  else
     F = 1./3.
  end if

  dtx = A*dt/gr_dx
  
  call soln_reconstruct(dt, Vj)
  call soln_getFlux()

  if (j .NE. 4) then
     !update cons variables to jth step only if not 4th step
     do i = gr_ibeg, gr_iend
        Uj(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
             dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))
     end do


     !! get updated primitive vars from the updated conservative vars
     do i = gr_ibeg, gr_iend
        ! Eos is automatically callled inside cons2prim
        call cons2prim(Uj(DENS_VAR:ENER_VAR,i),Vj(DENS_VAR:GAME_VAR,i))
     end do
     
  end if

  !finally add to total flux
  Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) + &
       F*gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)
  
end subroutine soln_RK4
