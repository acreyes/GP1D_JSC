module bc

#include "definition.h"

  use grid_data
  use sim_data, only : sim_bcType
  implicit none

contains

  subroutine bc_apply(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    if (sim_bcType == 'outflow') then
       call bc_outflow(V)
    elseif (sim_bcType == 'user') then
       call bc_user(V)
    elseif (sim_bcType == 'reflect') then
       call bc_reflect(V)
    elseif (sim_bcType == 'periodic') then
       call bc_periodic(V)
    endif

  end subroutine bc_apply

  subroutine bc_periodic(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    real, dimension(NUMB_VAR) :: Vbeg, Vend
    integer :: i

    do i = 1,gr_ngc
       Vbeg(1:NUMB_VAR) = V(1:NUMB_VAR, gr_ngc + i)
       Vend(1:NUMB_VAR) = V(1:NUMB_VAR, gr_iend - i+1)
       
       V(1:NUMB_VAR, gr_iend+i) = Vbeg(1:NUMB_VAR)
       V(1:NUMB_VAR, gr_ibeg-i) = Vend(1:NUMB_VAR)
    end do
  end subroutine bc_periodic

  subroutine bc_user(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    !BC for shu-osher problem.
    !don't do anything!
  end subroutine bc_user
  
  subroutine bc_outflow(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    integer :: i
 
    do i = 1, gr_ngc
       ! on the left GC
       V(1:NUMB_VAR,i) = V(1:NUMB_VAR,i+1) 

       ! on the right GC
       V(1:NUMB_VAR,gr_imax+1-i) = V(1:NUMB_VAR,gr_imax-i)
    end do

    return
  end subroutine bc_outflow

  subroutine bc_reflect(V)
    implicit none
    real, dimension(NUMB_VAR,gr_imax), intent(INOUT) :: V
    integer :: i,k0,k1

    do i = 1, gr_ngc
       k0 = 2*gr_ngc+1
       k1 = gr_iend-gr_ngc

       ! on the left GC
       V(       :,i) = V(       :,k0-i)
       V(VELX_VAR,i) =-V(VELX_VAR,k0-i)

       ! on the right GC
       V(       :,k1+k0-i) = V(         :,k1+i)
       V(VELX_VAR,k1+k0-i) =-V(  VELX_VAR,k1+i)
    end do

    return
  end subroutine bc_reflect


end module bc
