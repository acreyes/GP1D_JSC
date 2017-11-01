subroutine soln_update(dt)

#include "definition.h"
  
  use grid_data
  use sim_data
  use primconsflux, only : cons2prim

  implicit none
  real, intent(IN) :: dt
  integer :: i
  real :: dtx

  dtx = dt/gr_dx

  !! update conservative vars
  do i = gr_ibeg, gr_iend
     gr_U(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
          dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))
!!$     print *, (gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))!, gr_xCoord(i)
!!$     print *, gr_vR(1:7,i-1)
!!$     print *, gr_vL(1:7,i)
!!$     print *, gr_flux(:,i)
!!$     print*,
!!$     !print *, i, gr_U(:,i)
  end do


  !! get updated primitive vars from the updated conservative vars
  do i = gr_ibeg, gr_iend
     ! Eos is automatically callled inside cons2prim
!!$     if (isnan(gr_U(DENS_VAR,i))) then
!!$        print *, 'fdiff',(gr_flux(DENS_VAR:ENER_VAR,i+1) - gr_flux(DENS_VAR:ENER_VAR,i))!, gr_xCoord(i)
!!$       
!!$        print *, 'FL',gr_flux(:,i)
!!$        print *, 'vR',gr_vR(1:7,i-1)
!!$        print *, 'vi',gr_V(1:7,i)
!!$        print *, 'vL',gr_vL(1:7,i)
!!$        print *, 'FR',gr_flux(:,i+1)
!!$        print *, 'vR',gr_vR(1:7,i)
!!$        print *, 'vi+1',gr_V(1:7,i+1)
!!$        print *, 'vL',gr_vL(1:7,i+1)
!!$     end if
     call cons2prim(gr_U(DENS_VAR:ENER_VAR,i),gr_V(DENS_VAR:GAME_VAR,i))
     
     !print *, i, gr_v(1:7,i)
  end do
  

  return
end subroutine soln_update
