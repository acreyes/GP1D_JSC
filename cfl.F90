subroutine cfl(dt)

#include "definition.h"
  
  use grid_data
  use sim_data, only: sim_cfl, sim_Bx

  implicit none
  real, intent(OUT) :: dt
  integer :: i
  real :: maxSpeed, lambda, a, Ca, aCa2, Cf

  maxSpeed = -1e30
  !! update conservative vars
  do i = gr_ibeg, gr_iend
     a = sqrt(gr_V(GAMC_VAR,i)*gr_V(PRES_VAR,i)/gr_V(DENS_VAR,i))
     
     !alfven speed
     Ca = SQRT(sim_Bx*sim_Bx/gr_V(DENS_VAR,i))
     aCa2 = a*a + (sim_Bx*sim_Bx + gr_V(MAGY_VAR,i)*gr_V(MAGY_VAR,i) + gr_V(MAGZ_VAR,i)*gr_V(MAGZ_VAR,i))/gr_V(DENS_VAR,i)
     !fast/slow waves
     Cf = 1./sqrt(2.) * sqrt( aCa2 + sqrt( aCa2*aCa2 - 4.*a*a*Ca*Ca ) )

     lambda=(abs(gr_V(VELX_VAR,i)) + Cf)
     maxSpeed=max(maxSpeed,lambda)
  end do

  ! cfl
  dt = sim_cfl*gr_dx/maxSpeed

  return

end subroutine cfl
