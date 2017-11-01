module WENO

#include "definition.h"

contains

  subroutine gp_betas(V, R, beta)
    use sim_data , only: sim_gp_evals, sim_gp_evecs
    use grid_data, only: gr_dx, gr_gp_Zvecs2
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer :: i, N, s
    real    :: sum_b
    real, dimension(R+1) :: gamms, f
    real, dimension(R+1) :: beta_vec

    N = R+1
    !loop over stencils
    !primitive recons. stencil begins w/ 0.
    do s = 1 ,N
       !build stencil for reconstruction by primitive vars
       do i = 1, N
          f(i) = dot_product(gr_gp_Zvecs2(:,i,s),V(s:s+R))
       end do
       !loop over eigen values
       do i = 1, N
          !print *, sim_gp_evecs(:,i)
          gamms(i) = dot_product(f(:), sim_gp_evecs(:,i))
          !gamms(i) = dot_product(V(s:s+R), sim_gp_evecs(:,i))!/gr_dx**(1)!/sim_gp_evals(i)
          !gamms(i) = dot_product(f(:), sim_gp_evecs(:,i))
       end do
       !gamms = gamms/sqrt(dot_product(gamms,gamms))
       beta(s) = 0.
       do i = 1, N
          !beta_vec(i) = (abs(gamms(i))**2)/(sim_gp_evals(i))
          beta(s) = beta(s) + abs(gamms(i)**2)/(sim_gp_evals(i))
       end do
!!$       sum_b = sqrt(dot_product( beta_vec, beta_vec ))
!!$       beta(s) = sum(beta_vec(:))
       !beta(s) =  maxval(beta_vec(:))
    end do
    

!!$    call betas(V,R,weno_beta)
!!$    print *, sim_gp_evals
!!$    print *, sim_gp_evecs
!!$    print *, gamms
!!$    print *, beta
!!$    print *, weno_beta
!!$    print *, V
    !print *, 


    return
  end subroutine gp_betas

  subroutine betas(V, R, beta)
    !subroutine to calculate the smoothness-indicators for a WENO scheme on a 2R+1 point stencil
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer :: i

    i = R+1

    select case(R)
    case(1)
       beta(1) = (V(i  )-V(i-1))**2
       beta(2) = (V(i+1)-V(i  ))**2
    case(2)
       beta(1) = 13./12.*(V(i-2) - 2.*V(i-1) + V(i  ) )**2 + 0.25*(   V(i-2) - 4.*V(i-1) + 3.*V(i  ) )**2
       beta(2) = 13./12.*(V(i-1) - 2.*V(i  ) + V(i+1) )**2 + 0.25*(   V(i-1)             -    V(i+1) )**2
       beta(3) = 13./12.*(V(i  ) - 2.*V(i+1) + V(i+2) )**2 + 0.25*(3.*V(i  ) - 4.*V(i+1) +    V(i+2) )**2
    case(3)
       beta(1) = V(i-3)*(  547.*V(i-3) -  3882.*V(i-2) + 4642.*V(i-1) - 1854.*V(i  )) + &
                 V(i-2)*( 7043.*V(i-2) - 17246.*V(i-1) + 7042.*V(i  )               ) + &
                 V(i-1)*(11003.*V(i-1) -  9402.*V(i  )                              ) + 2107.*V(i  )**2
       beta(2) = V(i-2)*(  267.*V(i-2) -  1642.*V(i-1) + 1602.*V(i  ) -  494.*V(i+1)) + &
                 V(i-1)*( 2843.*V(i-1) -  5966.*V(i  ) + 1922.*V(i+1)               ) + &
                 V(i  )*( 3443.*V(i  ) -  2522.*V(i+1)                              ) +  547.*V(i+1)**2
       beta(3) = V(i-1)*(  547.*V(i-1) -  2522.*V(i  ) + 1922.*V(i+1) -  494.*V(i+2)) + &
                 V(i  )*( 3443.*V(i  ) -  5966.*V(i+1) + 1602.*V(i+2)               ) + &
                 V(i+1)*( 2843.*V(i+1) -  1642.*V(i+2)                              ) +  267.*V(i+2)**2
       beta(4) = V(i  )*( 2107.*V(i  ) -  9402.*V(i+1) + 7042.*V(i+2) - 1854.*V(i+3)) + &
                 V(i+1)*(11003.*V(i+1) - 17246.*V(i+2) + 4642.*V(i+3)               ) + &
                 V(i+2)*( 7043.*V(i+2) -  3882.*V(i+3)                              ) +  547.*V(i+3)**2
    case DEFAULT
       beta = 1.
    end select

    return
  end subroutine betas

end module WENO
