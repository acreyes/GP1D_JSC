subroutine hll(vL,vR,Flux)

#include "definition.h"  

  use grid_data
  use sim_data, only : sim_Bx
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR
  real :: sL,sR,aL2,aR2, Cax2, Ca2, CfL, CfR, Bx, By, Bz
  
  call prim2flux(vL,FL)
  call prim2flux(vR,FR)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)


!!$  print *, uL
!!$  print *, uR
!!$  print *,
  
  ! left and right sound speed a
  aL2 = vL(GAMC_VAR)*vL(PRES_VAR)/vL(DENS_VAR)
  aR2 = vR(GAMC_VAR)*vR(PRES_VAR)/vR(DENS_VAR)

  Bx = sim_Bx
  By = vL(MAGY_VAR)
  Bz = vL(MAGZ_VAR)

  Ca2  = Bx*Bx/vL(DENS_VAR)
  Cax2 = (Bx*Bx + By*By + Bz*Bz)/vL(DENS_VAR)

  CfL = sqrt(0.5*( aL2 + Cax2 + sqrt( (aL2+Cax2)*(aL2+Cax2) - 4.*aL2*Ca2 ) ))

  By = vR(MAGY_VAR)
  Bz = vR(MAGZ_VAR)

  Ca2  = Bx*Bx/vR(DENS_VAR)
  Cax2 = (Bx*Bx + By*By + Bz*Bz)/vR(DENS_VAR)

  CfR = sqrt(0.5*( aR2 + Cax2 + sqrt( (aR2+Cax2)*(aR2+Cax2) - 4.*aR2*Ca2 ) ))

  sL = min(vL(VELX_VAR) - CfL, vR(VELX_VAR) - cfR)
  sR = max(vL(VELX_VAR) + CfL, vR(VELX_VAR) + cfR)

  ! numerical flux
  if (sL >= 0.) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR)
  elseif ( (sL < 0.) .and. (sR >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = (    sR*FL(DENS_VAR:ENER_VAR) &
                                   -sL*FR(DENS_VAR:ENER_VAR) &
                               +sR*sL*(uR(DENS_VAR:ENER_VAR) &
                                      -uL(DENS_VAR:ENER_VAR)))/(sR-sL)
  else
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR)
  endif
  !print *, sL, sR, uR(2)-uL(2)
!!$  print *, FL
!!$  print *, Flux
!!$  print *, FR
!!$  print*,

  return
end subroutine hll
