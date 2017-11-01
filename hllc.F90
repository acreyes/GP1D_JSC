subroutine hllc(vL,vR,Flux)

#include "definition.h"  

  use grid_data
  use sim_data, only : sim_bx
  use primconsflux, only : prim2flux,prim2cons
  use eigensystem, only : eigenvalues


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR,Uhll,UstarL,UstarR
  real, dimension(NUMB_WAVE) :: lambdaL, lambdaR
  real, dimension(3) :: Bhll, bL, bR
  real :: sL,sR,aL,aR,qstar,pStar,pTotL, pTotR
  real :: numerL, denomL, numerR, denomR
  real :: dStarL, dStarR, bxstar

  real :: a2, bbx, bby, bbz, bb2, CfL, CfR, sqrtdi

  call prim2flux(vL,FL)
  call prim2flux(vR,FR)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)

  bL(1) = sim_bx
  bL(2:3) = vL(MAGY_VAR:MAGZ_VAR)
  bR(1) = sim_bx
  bR(2:3) = vR(MAGY_VAR:MAGZ_VAR)

  !total pressure w/ mag contribution
  pTotL = vL(PRES_VAR) + 0.5*dot_product(bL(1:3),bL(1:3))
  pTotR = vR(PRES_VAR) + 0.5*dot_product(bR(1:3),bR(1:3))
  !left state first
  sqrtdi = 1./sqrt(vL(DENS_VAR))
  bbx = sim_bx*sqrtdi
  bby = vL(MAGY_VAR)*sqrtdi
  bbz = vL(MAGZ_VAR)*sqrtdi

  bb2 = bbx*bbx + bby*bby + bbz*bbz
  a2 = vL(GAMC_VAR)*vL(PRES_VAR)/vL(DENS_VAR)

  CfL = sqrt( 0.5*( a2 + bb2 + sqrt( (a2+bb2)**2 - 4.*a2*bbx*bbx)))

  !right state
  sqrtdi = 1./sqrt(vR(DENS_VAR))
  bbx = sim_bx*sqrtdi
  bby = vR(MAGY_VAR)*sqrtdi
  bbz = vR(MAGZ_VAR)*sqrtdi

  bb2 = bbx*bbx + bby*bby + bbz*bbz
  a2 = vR(GAMC_VAR)*vR(PRES_VAR)/vR(DENS_VAR)

  CfR = sqrt( 0.5*( a2 + bb2 + sqrt( (a2+bb2)**2 - 4.*a2*bbx*bbx)))

  !get fastest left and right going waves
  sL = min(vL(VELX_VAR) - CfL, vR(VELX_VAR) - CfR)
  sR = max(vL(VELX_VAR) + CfR, vR(VELX_VAR) + CfR)
  
  
  ! Get HLL states for later use
  if (sL > 0.) then
     Uhll(DENS_VAR:ENER_VAR) = uL(DENS_VAR:ENER_VAR)
  elseif ((sL <= 0.) .and. (sR >= 0.)) then
     Uhll(DENS_VAR:ENER_VAR) = &
          ( sR*uR(DENS_VAR:ENER_VAR) &
           -sL*uL(DENS_VAR:ENER_VAR) &
             - FR(DENS_VAR:ENER_VAR) &
             + FL(DENS_VAR:ENER_VAR)&
             )/(sR - sL)
  else
     Uhll(DENS_VAR:ENER_VAR) = uR(DENS_VAR:ENER_VAR)
  endif

  
  Bhll(1) = sim_bx
  Bhll(2:3) = Uhll(MAGY_VAR:MAGZ_VAR)
  

  


  
  
  ! Get qstar
  qstar = vR(DENS_VAR)*vR(VELX_VAR)*(sR-vR(VELX_VAR)) &
       -vL(DENS_VAR)*vL(VELX_VAR)*(sL-vL(VELX_VAR)) &
       +pTotL - pTotR
  qstar = qstar - bL(1)*bL(1) + bR(1)*bR(1)
  qstar = qstar/( vR(DENS_VAR)*(sR-vR(VELX_VAR)) &
       -vL(DENS_VAR)*(sL-vL(VELX_VAR)))

  ! Convenient parameters
  numerL = sL-vL(VELX_VAR)
  denomL = sL-qstar
  numerR = sR-vR(VELX_VAR)
  denomR = sR-qstar

  

  ! Get pStar
  pStar = vL(DENS_VAR)*numerL*(qstar-vL(VELX_VAR)) + pTotL - bL(1)*bL(1) + Bhll(1)*Bhll(1)


  ! density
  dStarL = uL(DENS_VAR)*numerL/denomL
  dStarR = uR(DENS_VAR)*numerR/denomR

  ! left and right star regions
  UstarL(DENS_VAR) = dStarL
  UstarL(MOMX_VAR) = dStarL*qstar
  UstarL(MOMY_VAR) = uL(MOMY_VAR)*numerL/denomL
  UstarL(MOMZ_VAR) = uL(MOMZ_VAR)*numerL/denomL
  UstarL(ENER_VAR) = uL(ENER_VAR)*numerL/denomL+&
       (pStar*qstar - pTotL*vL(VELX_VAR))/denomL

  UstarR(DENS_VAR) = dStarR
  UstarR(MOMX_VAR) = dStarR*qstar
  UstarR(MOMY_VAR) = uR(MOMY_VAR)*numerR/denomR
  UstarR(MOMZ_VAR) = uR(MOMZ_VAR)*numerR/denomR
  UstarR(ENER_VAR) = uR(ENER_VAR)*numerR/denomR+&
       (pStar*qstar - pTotR*vR(VELX_VAR))/denomR

  !add MHD terms
  UstarL(MOMY_VAR) = UstarL(MOMY_VAR) - ( Bhll(1)*Bhll(2) - bL(1)*bL(2) )/denomL
  UstarL(MOMZ_VAR) = UstarL(MOMZ_VAR) - ( Bhll(1)*Bhll(3) - bL(1)*bL(3) )/denomL
  UstarL(ENER_VAR) = UstarL(ENER_VAR) &
       -((Bhll(1)*(dot_product(Bhll(:),Uhll(MOMX_VAR:MOMZ_VAR)))/Uhll(DENS_VAR)) &
       -(  bL(1)*(dot_product(  bL(:),  vL(VELX_VAR:VELZ_VAR)))               ))/denomL

  UstarR(MOMY_VAR) = UstarR(MOMY_VAR) - ( Bhll(1)*Bhll(2) - bR(1)*bR(2) )/denomR
  UstarR(MOMZ_VAR) = UstarR(MOMZ_VAR) - ( Bhll(1)*Bhll(3) - bR(1)*bR(3) )/denomR
  UstarR(ENER_VAR) = UstarR(ENER_VAR) &
       -((Bhll(1)*dot_product(Bhll(:),Uhll(MOMX_VAR:MOMZ_VAR))/Uhll(DENS_VAR)) &
       -(  bR(1)*dot_product(  bR(:),  vR(VELX_VAR:VELZ_VAR))               ))/denomR


  UstarL(MAGY_VAR:MAGZ_VAR) = Uhll(MAGY_VAR:MAGZ_VAR)
  UstarR(MAGY_VAR:MAGZ_VAR) = Uhll(MAGY_VAR:MAGZ_VAR)




  ! numerical flux
  if (sL >= 0.) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR)

  elseif ( (sL < 0.) .and. (qstar >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR) &
          + sL*(UstarL(DENS_VAR:ENER_VAR) - uL(DENS_VAR:ENER_VAR))
  elseif ( (qstar < 0.) .and. (sR >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR) &
          + sR*(UstarR(DENS_VAR:ENER_VAR) - uR(DENS_VAR:ENER_VAR))
  else
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR)
  endif

!!$  print *, sL, qstar, sR
!!$  print *, FL
!!$  print *, Flux
!!$  print *, FR
!!$  print *,



  return
end subroutine hllc
