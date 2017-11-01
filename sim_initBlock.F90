subroutine sim_initBlock()

#include "definition.h"
  
  use sim_data
  use grid_data, only : gr_V,gr_U,gr_i0,gr_imax,gr_xCoord,gr_dx,gr_ngc, gr_xbeg
  use primconsflux, only : prim2cons
  
  implicit none

  integer :: i
  real :: ekin, eint, x, del, gamm, dens, v2, b2
  
  ! generate x-coordinate
  do i = gr_i0,gr_imax
     gr_xCoord(i) = gr_xbeg + (real(i-gr_ngc)-0.5)*gr_dx
  end do

  do i = gr_i0,gr_imax
     if (sim_icType == 'BrioWu') then
        if (gr_xCoord(i) < sim_shockLoc) then
           !left state
           gr_V(DENS_VAR,i) = 1.
           gr_V(VELX_VAR,i) = 0.
           gr_V(VELY_VAR,i) = 0.
           gr_V(VELZ_VAR,i) = 0.

           gr_V(MAGY_VAR,i) = 1.
           gr_V(MAGZ_VAR,i) = 0.
           gr_V(PRES_VAR,i) = 1.
        else
           !right state
           gr_V(DENS_VAR,i) = 0.125
           gr_V(VELX_VAR,i) = 0.
           gr_V(VELY_VAR,i) = 0.
           gr_V(VELZ_VAR,i) = 0.

           gr_V(MAGY_VAR,i) = -1.
           gr_V(MAGZ_VAR,i) = 0.
           gr_V(PRES_VAR,i) = 0.1
        end if
        elseif (sim_icType == 'sod') then
           if (gr_xCoord(i) < sim_shockLoc) then
              !left state
              gr_V(DENS_VAR,i) = 1.
              gr_V(VELX_VAR,i) = 0.
              gr_V(VELY_VAR,i) = 0.
              gr_V(VELZ_VAR,i) = 0.
              
              gr_V(MAGY_VAR,i) = 0.
              gr_V(MAGZ_VAR,i) = 0.
              gr_V(PRES_VAR,i) = 1.
           else
              !right state
              gr_V(DENS_VAR,i) = 0.125
              gr_V(VELX_VAR,i) = 0.
              gr_V(VELY_VAR,i) = 0.
              gr_V(VELZ_VAR,i) = 0.
              
              gr_V(MAGY_VAR,i) = 0.
              gr_V(MAGZ_VAR,i) = 0.
              gr_V(PRES_VAR,i) = 0.1
           end if
        elseif (sim_icType == 'shu') then
           !do IC for shu-osher problem
           !transform the domain [0,1] onto the one give for the problem:[-4.5,4.5]
           x = gr_xCoord(i) - 4.5
           gr_V(VELY_VAR:MAGZ_VAR,i) = 0.
           if (x < -4.) then
              !left state
              gr_V(DENS_VAR,i) = 3.857143
              gr_V(VELX_VAR,i) = 2.629369
              gr_V(PRES_VAR,i) = 10.33333
           else
              gr_V(DENS_VAR,i) = 1 + .2*SIN(5.*x)
              gr_V(VELX_VAR,i) = 0.
              gr_V(PRES_VAR,i) = 1.
           end if
        elseif (sim_icType == 'shock') then
           !uses sim_pres to get pressure
            if (gr_xCoord(i) < sim_shockLoc) then
              gr_V(DENS_VAR,i) = sim_densL
              gr_V(VELX_VAR,i) = sim_velxL
              gr_V(VELY_VAR,i) = sim_velyL
              gr_V(VELZ_VAR,i) = sim_velzL
              gr_V(MAGY_VAR,i) = sim_magyL
              gr_V(MAGZ_VAR,i) = sim_magzL
              gr_V(PRES_VAR,i) = sim_presL
           else
              gr_V(DENS_VAR,i) = sim_densR
              gr_V(VELX_VAR,i) = sim_velxR
              gr_V(VELY_VAR,i) = sim_velyR
              gr_V(VELZ_VAR,i) = sim_velzR
              gr_V(MAGY_VAR,i) = sim_magyR
              gr_V(MAGZ_VAR,i) = sim_magzR
              gr_V(PRES_VAR,i) = sim_presR
           end if
        elseif (sim_icType == 'gauss') then
           x = gr_xCoord(i) - 0.5
           del = 0.1
           gr_V(DENS_VAR,i) = 1. + EXP(-x**2/(del)**2)
           gr_V(VELX_VAR,i) = 1.
           gr_V(PRES_VAR,i) = 1./sim_gamma

        elseif (sim_icType == 'RJ') then
           !uses energy to determine pressure
           !this is to make it easier to run tests from Ryu and Jones paper
           if (gr_xCoord(i) < sim_shockLoc) then
              gr_V(DENS_VAR,i) = sim_densL
              gr_V(VELX_VAR,i) = sim_velxL
              gr_V(VELY_VAR,i) = sim_velyL
              gr_V(VELZ_VAR,i) = sim_velzL
              gr_V(MAGY_VAR,i) = sim_magyL
              gr_V(MAGZ_VAR,i) = sim_magzL
              gr_U(ENER_VAR,i) = sim_enerL
           else
              gr_V(DENS_VAR,i) = sim_densR
              gr_V(VELX_VAR,i) = sim_velxR
              gr_V(VELY_VAR,i) = sim_velyR
              gr_V(VELZ_VAR,i) = sim_velzR
              gr_V(MAGY_VAR,i) = sim_magyR
              gr_V(MAGZ_VAR,i) = sim_magzR
              gr_U(ENER_VAR,i) = sim_enerR
           end if
           v2 = dot_product(gr_V(VELX_VAR:VELZ_VAR,i),gr_V(VELX_VAR:VELZ_VAR,i))
           b2 = sim_bx*sim_bx + dot_product(gr_V(MAGY_VAR:MAGZ_VAR,i),gr_V(MAGY_VAR:MAGZ_VAR,i))
           gr_V(PRES_VAR,i) = max((sim_gamma - 1.)*(gr_U(ENER_VAR,i) - 0.5*(gr_V(DENS_VAR,i)*v2 + b2)),sim_smallPres)
           print *, gr_V(PRES_VAR,i), (sim_gamma - 1.)*(gr_U(ENER_VAR,i) - 0.5*(gr_V(DENS_VAR,i)*v2 + b2)), gr_U(ENER_VAR,i)
        
     end if
     gr_V(GAMC_VAR,i) = sim_gamma
     gr_V(GAME_VAR,i) = sim_gamma
     gr_V(EINT_VAR,i) = gr_V(PRES_VAR,i)/(gr_V(GAME_VAR,i)-1.)/gr_V(DENS_VAR,i)
     
  end do

  ! also initialize conservative vars
  do i = gr_i0,gr_imax
     call prim2cons(gr_V(:,i), gr_U(DENS_VAR:ENER_VAR,i))
  end do

  
end subroutine sim_initBlock
