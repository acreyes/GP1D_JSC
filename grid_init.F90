subroutine grid_init()

#include "definition.h"

  use grid_data
  use sim_data
  use read_initFile

  implicit none

  call read_initFileInt ('slug.init','gr_nx',   gr_nx)
  call read_initFileInt ('slug.init','gr_radius',   gr_radius)
  call read_initFileInt ('slug.init','gr_ngc',  gr_ngc)
  call read_initFileReal('slug.init','gr_xbeg', gr_xbeg)
  call read_initFileReal('slug.init','gr_xend', gr_xend)

  ! allocate cell coordinates
  allocate(gr_xCoord(gr_nx+2*gr_ngc)); gr_xCoord = 0.0

  ! the first and the last interior cell index
  gr_ibeg = gr_ngc + 1
  gr_iend = gr_ngc + gr_nx

  ! the min and max of the entire cell index
  gr_i0   = 1
  gr_imax = 2*gr_ngc + gr_nx

  ! grid delta
  gr_dx = (gr_xend - gr_xbeg)/gr_nx

  ! allocate grid variables
  allocate(gr_U(NSYS_VAR,gr_imax)); gr_U = 0.
  allocate(gr_V(NUMB_VAR,gr_imax)); gr_V = 0.
  allocate(gr_W(NSYS_VAR,gr_imax)); gr_W = 0.

  ! allocate grid Riemann states
  allocate(gr_vL(NUMB_VAR,gr_imax)); gr_vL = 0.
  allocate(gr_vR(NUMB_VAR,gr_imax)); gr_vR = 0.

  ! allocate grid fluxes
  allocate(gr_flux(NSYS_VAR,gr_imax)); gr_flux = 0.

  ! allocate grid eigensystem
  allocate(gr_eigval(NUMB_WAVE,gr_imax)); gr_eigval = 0.
  allocate(gr_leigvc(NSYS_VAR,NUMB_WAVE,gr_imax)); gr_leigvc = 0.
  allocate(gr_reigvc(NSYS_VAR,NUMB_WAVE,gr_imax)); gr_reigvc = 0.


  ! allocate GP variables
  allocate(gr_GPv(2*gr_radius+1   )); gr_GPv = 0.
  allocate(gr_GPZ(2, 2*gr_radius+1)); gr_GPZ = 0.
  allocate(gr_GPv2(2*gr_radius    )); gr_GPv2 = 0.
  allocate(gr_GPZ2(2, 2*gr_radius )); gr_GPZ2 = 0.

!!$  allocate(gr_GPvk(3, 3   ) ); gr_GPvk = 0.
!!$  allocate(gr_GPZk(2, 3, 3) ); gr_GPZk = 0.
!!$  allocate(gr_GPvk2(3, 3   ) ); gr_GPvk = 0.
!!$  allocate(gr_GPZk2(2, 3, 3) ); gr_GPZk = 0.

  allocate(gr_GPvk(gr_radius+1, gr_radius+1   ) ); gr_GPvk = 0.
  allocate(gr_GPZk(2, gr_radius+1, gr_radius+1) ); gr_GPZk = 0.
  allocate(gr_GPvk2(gr_radius+1, gr_radius+1   ) ); gr_GPvk = 0.
  allocate(gr_GPZk2(2, gr_radius+1, gr_radius+1) ); gr_GPZk = 0.

  allocate(gr_gp_Zvecs( gr_radius+1,gr_radius+1,gr_radius+1)); gr_gp_Zvecs = 0.
  allocate(gr_gp_Zvecs2(gr_radius+1,gr_radius+1,gr_radius+1)); gr_gp_Zvecs2 = 0.
  
  allocate(gr_U2_16(2)); gr_U2_16 = 0.
  allocate(gr_U2(2)); gr_U2 = 0.

  allocate(gr_GP_w2(2,gr_radius+1)); gr_GP_w2 = 0.
  allocate(gr_GP_w( 2,gr_radius+1)); gr_GP_w  = 0.

  allocate( gr_GP_Kki(gr_radius+1,gr_radius+1) ); gr_GP_Kki = 0.
  allocate( gr_GP_Kki2(gr_radius+1,gr_radius+1) ); gr_GP_Kki2 = 0.

  allocate( gr_GP_char_Kki(gr_radius+1,gr_radius+1) ); gr_GP_char_Kki = 0.
  allocate( gr_GP_char_Kki2(gr_radius+1,gr_radius+1) ); gr_GP_char_Kki2 = 0.
  return
end subroutine grid_init
