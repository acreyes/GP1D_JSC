subroutine soln_ReconEvolveAvg(dt)

#include "definition.h"  

  use grid_data
  use sim_data
  use bc

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax) :: Vj
  real, dimension(NSYS_VAR,gr_imax) :: Flux
  integer :: j

  
  !RK4 in time
  Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = 0.
  !initialize Vj as V0
  Vj(DENS_VAR:NUMB_VAR,gr_i0:gr_imax) = gr_V(DENS_VAR:NUMB_VAR,gr_i0:gr_imax)
  !do rk4 steps
  do j = 1, 4
     !each step adds the flux at the kj'th state with the appropriate weight
     !to the total flux
     call soln_RK4(dt, j, Vj, Flux)
     call bc_apply(Vj)
  end do
  gr_flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax) = Flux(DENS_VAR:ENER_VAR,gr_i0:gr_imax)
     
  return
end subroutine soln_ReconEvolveAvg
