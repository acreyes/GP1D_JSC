subroutine soln_getFlux()

#include "definition.h"  

  use grid_data
  use sim_data

  implicit none
  integer :: i

  if (sim_riemann == 'hll') then
     do i = gr_ibeg, gr_iend+1
        call hll(gr_vR(DENS_VAR:GAME_VAR,i-1),&
                 gr_vL(DENS_VAR:GAME_VAR,i  ),&
                 gr_flux(DENS_VAR:ENER_VAR,i))
        !print *, i, gr_flux(MOMX_VAR, i)
!!$        print *, gr_VR(1:7,i)
!!$        print *, gr_V(1:7,i)
!!$        print *, gr_vL(1:7,i)
!!$        print*,
     enddo

  elseif (sim_riemann == 'roe') then
     do i = gr_ibeg, gr_iend+1
        call roe(gr_vR(DENS_VAR:GAME_VAR,i-1),&
                 gr_vL(DENS_VAR:GAME_VAR,i  ),&
                 gr_flux(DENS_VAR:ENER_VAR,i))
     enddo
  elseif (sim_riemann == 'hllc') then
     do i = gr_ibeg, gr_iend+1
        call hllc(gr_vR(DENS_VAR:GAME_VAR,i-1),&
                 gr_vL(DENS_VAR:GAME_VAR,i  ),&
                 gr_flux(DENS_VAR:ENER_VAR,i))
        !print *, gr_flux(:,i), gr_xCoord(i)
        !print *, gr_vL(1:NSYS_VAR, i), gr_xCoord(i)
     enddo
  endif



  return
end subroutine soln_getFlux
