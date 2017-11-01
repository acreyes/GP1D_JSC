subroutine soln_reconstruct(dt, V)

#include "definition.h"  

  use grid_data
  use sim_data

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  
  if (sim_order == 5) then
     call soln_WENO(dt, V)
  elseif (sim_order == 9) then
     call soln_GP2(dt, V)
  end if

  return
end subroutine soln_reconstruct
