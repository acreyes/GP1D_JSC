subroutine soln_GP2(dt, V)
#include "definition.h"

  use grid_data
  use sim_data
  use eigensystem
  use WENO

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax), intent(IN) :: V

  integer :: i, var, s, N, k, R
  logical :: conservative

  real :: sum_wbar
  real, dimension(NSYS_VAR) :: vecL, vecR, vL, vR, C0, C1, C2

  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig0, leig0
  real, dimension(2*gr_radius+1) :: G, un
  real, dimension(gr_radius+1) :: vLk, vRk, beta_k
  real, dimension(2, gr_radius+1) :: weights, wbar, beta



  conservative = .false.
  R = gr_radius
  N = 2*R+1


  do i = gr_ibeg-1,gr_iend+1
     
     !copy cell-centered values to left and right states
     gr_vL(DENS_VAR:NUMB_VAR,i) = V(DENS_VAR:NUMB_VAR,i)
     gr_vR(DENS_VAR:NUMB_VAR,i) = V(DENS_VAR:NUMB_VAR,i)

     !get eigen information
     call eigenvalues(V(DENS_VAR:GAME_VAR,i), lambda)
     call mhd_eigenvectors(V(DENS_VAR:GAME_VAR,i), conservative, reig0, leig0)
     !this is the global stencil
     do var = DENS_VAR, PRES_VAR
        do s = 1, N
           if (sim_charLimiting) then
              G(s) = dot_product(V(DENS_VAR:PRES_VAR, i+s-gr_radius-1), leig0(DENS_VAR:PRES_VAR, var))
           else
              G(s) = V(var, i+(s-gr_radius-1))
           end if
        end do

        !now we begin the GP reconstruction

        !now lets compute the GP prediction over the smaller candidate stencils
        do k = 1, R+1
           vLk(k) = dot_product(gr_GPzk2(1, :, k), G(k:k+R))
           vRk(k) = dot_product(gr_GPzk2(2, :, k), G(k:k+R))
        end do

        
        call gp_betas(G, R, beta_k)

        do k = 1, R+1
           wbar(1, k) = gr_GP_w2(1, k)/(1.e-36+beta_k(k))**sim_mval
           wbar(2, k) = gr_GP_w2(2, k)/(1.e-36+beta_k(k))**sim_mval
        end do
        
        sum_wbar = SUM(wbar(1,:))
        weights(1,:) = wbar(1,:)/sum_wbar
        sum_wbar = SUM(wbar(2,:))
        weights(2,:) = wbar(2,:)/sum_wbar

        !WENO-M weights
        do k = 1, R+1
           do s = 1, 2
              wbar(s, k) = weights(s,k) * ( gr_GP_w2(s,k) + gr_GP_w2(s,k)**2 - &
                   3.*gr_GP_w2(s,k)*weights(s,k) + weights(s,k)**2)/( gr_GP_w2(s,k)**2 + &
                   weights(s,k)*(1.-2.*gr_GP_w2(s,k)) )
           end do
        end do

        sum_wbar = SUM(wbar(1,:))
        weights(1,:) = wbar(1,:)/sum_wbar
        sum_wbar = SUM(wbar(2,:))
        weights(2,:) = wbar(2,:)/sum_wbar
        

        vecL(var) = dot_product(weights(1,:), vLk(:))
        vecR(var) = dot_product(weights(2,:), vRk(:))
       
     end do !var
      
     do var = DENS_VAR, PRES_VAR
        if (sim_charLimiting) then
           !project char vars back onto prim vars
           vL(var) = dot_product(reig0(var,1:NUMB_WAVE),vecL(1:NUMB_WAVE))
           vR(var) = dot_product(reig0(var,1:NUMB_WAVE),vecR(1:NUMB_WAVE))
        else
           vL(var) = vecL(var)
           vR(var) = vecR(var)
        end if
     end do !var

     gr_vL(DENS_VAR:PRES_VAR, i) = vL(DENS_VAR:PRES_VAR)
     gr_vR(DENS_VAR:PRES_VAR, i) = vR(DENS_VAR:PRES_VAR)
  end do !i



end subroutine soln_GP2
