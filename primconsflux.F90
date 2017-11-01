module primconsflux

#include "definition.h"
  
  use grid_data
  use sim_data, only : sim_gamma, sim_smallPres, sim_Bx
  use eos, only : eos_cell
  
contains

  subroutine prim2cons(V,U)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: U

    real :: ekin, eint, emag

    U(DENS_VAR) = V(DENS_VAR)
    U(MOMX_VAR) = V(DENS_VAR)*V(VELX_VAR)
    U(MOMY_VAR) = V(DENS_VAR)*V(VELY_VAR)
    U(MOMZ_VAR) = V(DENS_VAR)*V(VELZ_VAR)

    U(MAGY_VAR) = V(MAGY_VAR)
    U(MAGZ_VAR) = V(MAGZ_VAR)

    ekin = 0.5*V(DENS_VAR)*dot_product(V(VELX_VAR:VELZ_VAR), V(VELX_VAR:VELZ_VAR))
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    emag = 0.5*(dot_product(V(MAGY_VAR:MAGZ_VAR),V(MAGY_VAR:MAGZ_VAR)) + sim_Bx*sim_Bx)
    U(ENER_VAR) = ekin + eint + emag
         
  end subroutine prim2cons


  subroutine cons2prim(U,V)
    implicit none
    real, dimension(NSYS_VAR), intent(IN)  :: U
    real, dimension(NUMB_VAR), intent(OUT) :: V
    real :: eint, ekin, emag, pres, B2

    B2 = dot_product(U(MAGY_VAR:MAGZ_VAR),U(MAGY_VAR:MAGZ_VAR)) + sim_Bx*sim_Bx
    
    V(DENS_VAR) = U(DENS_VAR)
    V(VELX_VAR) = U(MOMX_VAR)/U(DENS_VAR)
    V(VELY_VAR) = U(MOMY_VAR)/U(DENS_VAR)
    V(VELZ_VAR) = U(MOMZ_VAR)/U(DENS_VAR)

    V(MAGY_VAR) = U(MAGY_VAR)
    V(MAGZ_VAR) = U(MAGZ_VAR)
    
    ekin = 0.5*V(DENS_VAR)*dot_product(V(VELX_VAR:VELZ_VAR), V(VELX_VAR:VELZ_VAR))
    emag = 0.5*B2
    eint = max(U(ENER_VAR) - ekin - emag, sim_smallPres) !eint=rho*e
    eint = eint/U(DENS_VAR)
    ! get pressure by calling eos
    call eos_cell(U(DENS_VAR),eint,sim_gamma,pres)
    V(PRES_VAR) = pres

    V(EINT_VAR) = eint*U(DENS_VAR)
    V(GAMC_VAR) = sim_gamma
    V(GAME_VAR) = sim_gamma
    
  end subroutine cons2prim

  subroutine prim2flux(V,Flux)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: Flux

    real :: ekin,eint,ener, Ptot, B2, v2

    B2 = dot_product(V(MAGY_VAR:MAGZ_VAR),V(MAGY_VAR:MAGZ_VAR)) + sim_Bx*sim_Bx
    v2 = dot_product(V(VELX_VAR:VELZ_VAR),V(VELX_VAR:VELZ_VAR))

    !total pressure is mag + hydro
    Ptot = V(PRES_VAR) + 0.5*B2
!    print *, Ptot
    
    Flux(DENS_VAR) = V(DENS_VAR)*V(VELX_VAR)
    
    Flux(MOMX_VAR) = Flux(DENS_VAR)*V(VELX_VAR) + Ptot - sim_Bx*sim_Bx
    Flux(MOMY_VAR) = Flux(DENS_VAR)*V(VELY_VAR) - sim_Bx*V(MAGY_VAR)
    Flux(MOMZ_VAR) = Flux(DENS_VAR)*V(VELZ_VAR) - sim_Bx*V(MAGZ_VAR)

    Flux(MAGY_VAR) = V(MAGY_VAR)*V(VELX_VAR) - sim_Bx*V(VELY_VAR)
    Flux(MAGZ_VAR) = V(MAGZ_VAR)*V(VELX_VAR) - sim_Bx*V(VELZ_VAR)
    
    ekin = 0.5*V(DENS_VAR)*v2
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    ener = ekin + eint + 0.5*B2
    Flux(ENER_VAR) = V(VELX_VAR)*(ener + Ptot) - &
         sim_Bx*(sim_Bx*V(VELX_VAR) + V(MAGY_VAR)*V(VELY_VAR) + V(MAGZ_VAR)*V(VELZ_VAR))
    
  end subroutine prim2flux

  subroutine cons2flux
    implicit none
  end subroutine cons2flux
  
end module primconsflux
