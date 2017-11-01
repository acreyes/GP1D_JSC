module eos

#include "definition.h"
  
  use grid_data
  use sim_data, only : sim_smallPres

contains

  subroutine eos_all()
    implicit none

    integer :: i
    real :: eint, ekin,pres,velx
    

    ! dens-eint mode:
    ! This assumes that conservative vars are known
    ! (but not prim vars yet) to evaluate pres
    ! == inputs : dens & eint
    ! == outputs: pres
    ! This routine fills out pressure also the GC regions
    !I dont' think this gets used anywhere anymore...
    
    do i =  gr_i0,gr_imax
       ! updated velx, ekin, eint
       velx = gr_U(MOMX_VAR,i)/gr_U(DENS_VAR,i)
       ekin = 0.5*velx*velx*gr_U(DENS_VAR,i)
       eint = max(gr_U(ENER_VAR,i) - ekin,sim_smallPres)

       ! now let's get a new pressure
       pres = (gr_V(GAME_VAR,i)-1.)*eint
       gr_V(PRES_VAR,i) = max(sim_smallPres,pres)
    end do
  end subroutine eos_all


  
  subroutine eos_cell(dens,eint,game,pres)
    implicit none
    real, intent(IN) :: dens,eint,game
    real, intent(OUT):: pres
    
    ! ideal gas law
    ! eint = pres/dens/(game-1)
    pres = max((game-1.)*dens*eint,sim_smallPres)

    
  end subroutine eos_cell

end module eos
