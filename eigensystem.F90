module eigensystem

#include "definition.h"
  
  use grid_data
  use sim_data, only : sim_Bx

contains

  subroutine eigenvalues(V,lambda)
    implicit none

    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NUMB_WAVE),intent(OUT) :: lambda

    real :: a, u, Cf, Cs, Ca, aCa2

    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VELX_VAR)

    !alfven speed
    Ca = SQRT(sim_Bx*sim_Bx/V(DENS_VAR))
    aCa2 = a*a + (sim_Bx*sim_Bx + V(MAGY_VAR)*V(MAGY_VAR) + V(MAGZ_VAR)*V(MAGZ_VAR))/V(DENS_VAR)
    !fast/slow waves
    Cf = 1./sqrt(2.) * sqrt( aCa2 + sqrt( aCa2*aCa2 - 4.*a*a*Ca*Ca ) )
    Cs = 1./sqrt(2.) * sqrt( aCa2 - sqrt( aCa2*aCa2 - 4.*a*a*Ca*Ca ) )
    
    lambda(SHOCKLEFT) = u - Cf
    lambda(ALFVNLEFT) = u - Ca
    lambda(SLOWWLEFT) = u - Cs
    lambda(CTENTROPY) = u
    lambda(SLOWWRGHT) = u + Cs
    lambda(ALFVNRGHT) = u + Ca
    lambda(SHOCKRGHT) = u + Cf
    
    return
  end subroutine eigenvalues

  subroutine eigen_params(V, Cf, Ca, Cs, a, A_f, A_s, beta)
    !There are a variety of sources that describe the eigensystem for MHD in 1D
    !I ended up using that described in Roe & Balsara, (Roe & Balsara, SIAM, 1996)
    !It looks like this is what is used in Flash, and it is also the simplest
    !description I have seen.
    
    implicit none

    real, dimension(NUMB_VAR), intent(IN) :: V
    
    real, intent(OUT) :: Cf, Ca, Cs, a, A_f, A_s
    real, dimension(2), intent(OUT) :: beta

    real :: a2, Cs2, Cf2, A_s2, A_f2, u, u2
    real :: abx, bbx, bby, bbz, bb2, bbT, sqrtdi, eps

    Ca = 0.
    Cs = 0.
    Cf = 0.
    A_f = 0.
    A_s = 0.
    beta = 0.
    eps = 1.e-14

    u = V(VELX_VAR)
    a2 = V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR)
    a = sqrt(a2)

    sqrtdi = 1./sqrt(V(DENS_VAR))
    bbx = sim_bx*sqrtdi
    abx = abs(bbx)
    bby = V(MAGY_VAR)*sqrtdi
    bbz = V(MAGZ_VAR)*sqrtdi
    bb2 = bbx*bbx + bby*bby + bbz*bbz

    Ca = bbx
    bbT = sqrt(bby*bby + bbz*bbz)
    u2 = dot_product(V(VELX_VAR:VELZ_VAR),V(VELX_VAR:VELZ_VAR))

    Cf2 = 0.5*(a2 + bb2 + sqrt((a2-bb2)*(a2-bb2) + 4.*a2*bbT*bbT ) )
    Cs2 = a2*Ca*Ca/Cf2

    A_f2 = (a2-Cs2)/(Cf2-Cs2)
    A_s2 = 1. - A_f2

    !renormalization coefficients
    !this is more directly from ROE & Balsara
!!$    if (bbT < eps) then
!!$       !leads to cases II, III, IV, V
!!$       if (abx < eps) then
!!$          !case II
!!$          Cf2 = a2
!!$          Cs2 = 0.
!!$          A_f2 = 1.
!!$          A_s2 = 0.
!!$       elseif (abx < a - eps) then
!!$          !case III
!!$          Cf2 = a2
!!$          Cs2 = bbx*bbx
!!$          A_f2 = 1.
!!$          A_s2 = 0.
!!$       elseif (abx - a > eps) then
!!$          !case IV
!!$          Cf2 = bbx*bbx
!!$          Cs2 = a2
!!$          A_f2 = 0.
!!$          A_s2 = 1.
!!$       elseif (abx - a < eps) then
!!$          !case V
!!$          A_f2 = 0.5
!!$          A_s2 = 0.5
!!$       end if
!!$
!!$    elseif (abx < eps) then
!!$       !case I
!!$       Cf2 = a2 + bbT*bbT
!!$       Cs2 = a2*bbx*bbx/Cf2
!!$       A_f2 = a2/Cf2
!!$       A_s2 = bbT*bbT/Cf2
!!$
!!$    end if
       

    !this is from FLASH, and I don't really understand it
    if (bb2 > 0.) then
       if (abs(Cf2-Cs2) > 1.e-16*a2) then
          A_f2 = min(1., max(0.,(a2-Cs2)/(Cf2-Cs2) ) )
          A_s2 = 1. - A_f2
          !print *, 'A',A_f2, A_s2, (a2-Cs2)/(Cf2-Cs2)
       else
          A_f2 = 0.5
          A_s2 = 0.5
       end if

    else
       A_f2 = 1.
       A_s2 = 0.
       Cf2 = a2
       Cs2 = 0.
       Ca = 0.

    end if

    !indeterminante betas
    if (bbT > 0.) then
       beta = (/bby, bbz/)/bbT
    else
       beta = sqrt(0.5)
    end if

    Cf  = sqrt(Cf2)
    Cs  = sqrt(Cs2)
    A_f = sqrt(A_f2)
    A_s = sqrt(A_s2)
    a   = sqrt(a2)
    

    

    return
  end subroutine eigen_params

  subroutine mhd_eigenvectors(V, conservative, reig, leig)
    !There are a variety of sources that describe the eigensystem for MHD in 1D
    !I ended up using that described in Roe & Balsara, (Roe & Balsara, SIAM, 1996)
    !It looks like this is what is used in Flash, and it is also the simplest
    !description I have seen.
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical, intent(IN) :: conservative !so far I'm only doing primitive eigenvectors
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig, leig

    !we're gonna get these guys with eigen_params
    real               :: Cf, Ca, Cs, a, A_f, A_s
    real, dimension(2) :: beta

    !these will be useful
    real :: sqrtd, sqrtdi, signB, k, dinv, a2inv, a2, u2
    real :: a_f1, a_f2, a_f3, a_s1, a_s2, a_s3, dot

    integer :: varL, varR, ii, jj

    real, dimension(NSYS_VAR,NSYS_VAR)  :: JacQ, JacQi
    real, dimension(NSYS_VAR,NUMB_WAVE) :: LeigTemp
    real, dimension(NUMB_WAVE,NSYS_VAR) :: ReigTemp

    !lets initialize w/ 0 and only fill in non-zero ones
    leig = 0.
    reig = 0.
    
    call eigen_params(V, Cf, Ca, Cs, a, A_f, A_s, beta)

    dinv = 1./V(DENS_VAR)
    k = 1. - V(GAME_VAR)                                      ! 1 - gamma
    u2 = dot_product(V(VELX_VAR:VELZ_VAR),V(VELX_VAR:VELZ_VAR))
    sqrtd  = sqrt(V(DENS_VAR))
    sqrtdi = 1./sqrtd
    signB  = sign(1., sim_bx)                                     !Ca = sim_bx/sqrt(dens)
    a2     = V(GAMC_VAR)*V(PRES_VAR)*dinv
    a2inv  = 1./a2
    

    a_f1   = A_f*Cf*signB
    a_f2   = A_f*a*sqrtdi
    a_f3   = A_f*a*sqrtd

    a_s1   = A_s*Cs*signB
    a_s2   = A_s*a*sqrtdi
    a_s3   = A_s*a*sqrtd

!!$    print *, A_f, A_s
!!$    print *, a_f1, a_f2, a_f3
!!$    print *, a_s1, a_s2, a_s3

    !primitive eigenvectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!lets start with the left eigenvectors!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
    
    !Fast waves
    !leig(DENS_VAR,SHOCKLEFT) = 0. 
    leig(VELX_VAR,SHOCKLEFT) = -A_f*Cf 
    leig(VELY_VAR,SHOCKLEFT) = a_s1*beta(1)
    leig(VELZ_VAR,SHOCKLEFT) = a_s1*beta(2)
    leig(MAGY_VAR,SHOCKLEFT) = a_s2*beta(1)
    leig(MAGZ_VAR,SHOCKLEFT) = a_s2*beta(2)
    leig(PRES_VAR,SHOCKLEFT) = A_f*dinv
    ! SHOCKRGHT same up to some - signs
    leig(VELX_VAR:VELZ_VAR, SHOCKRGHT) = -leig(VELX_VAR:VELZ_VAR, SHOCKLEFT)
    leig(MAGY_VAR:PRES_VAR, SHOCKRGHT) =  leig(MAGY_VAR:PRES_VAR, SHOCKLEFT)


    !Alfven waves
    !leig(DENS_VAR,ALFVNLEFT) = 0.
    !leig(VELX_VAR,ALFVNLEFT) = 0.
    leig(VELY_VAR,ALFVNLEFT) = -beta(2)
    leig(VELZ_VAR,ALFVNLEFT) =  beta(1)
    leig(MAGY_VAR,ALFVNLEFT) = -beta(2)*sqrtdi
    leig(MAGZ_VAR,ALFVNLEFT) =  beta(1)*sqrtdi
    !leig(PRES_VAR,ALFVNLEFT) = 0.

    leig(VELY_VAR:VELZ_VAR,ALFVNRGHT) = -leig(VELY_VAR:VELZ_VAR, ALFVNLEFT)
    leig(MAGY_VAR:MAGZ_VAR,ALFVNRGHT) =  leig(MAGY_VAR:MAGZ_VAR, ALFVNLEFT)

    !SLOW waves
    !leig(DENS_VAR,SLOWWLEFT) = 0.
    leig(VELX_VAR,SLOWWLEFT) = -A_s*Cs
    leig(VELY_VAR,SLOWWLEFT) = -a_f1*beta(1)
    leig(VELZ_VAR,SLOWWLEFT) = -a_f1*beta(2)
    leig(MAGY_VAR,SLOWWLEFT) = -a_f2*beta(1)
    leig(MAGZ_VAR,SLOWWLEFT) = -a_f2*beta(2)
    leig(PRES_VAR,SLOWWLEFT) = A_s*dinv

    leig(VELX_VAR:VELZ_VAR,SLOWWRGHT) = -leig(VELX_VAR:VELZ_VAR,SLOWWLEFT)
    leig(MAGY_VAR:PRES_VAR,SLOWWRGHT) =  leig(MAGY_VAR:PRES_VAR,SLOWWLEFT)


    ! Entropy Wave
    leig(DENS_VAR,CTENTROPY) =  1.
    !leig(VELX_VAR,CTENTROPY) = 0.
    !leig(VELY_VAR,CTENTROPY) = 0
    !leig(VELZ_VAR,CTENTROPY) = 0
    !leig(MAGY_VAR,CTENTROPY) = 0
    !leig(MAGZ_VAR,CTENTROPY) = 0
    leig(PRES_VAR,CTENTROPY) = -1.*a2inv

    !scale fast and slow waves by 1/2a2
    leig(:,SHOCKLEFT) = 0.5*a2inv*leig(:,SHOCKLEFT)
    leig(:,SHOCKRGHT) = 0.5*a2inv*leig(:,SHOCKRGHT)

    leig(:,SLOWWRGHT) = 0.5*a2inv*leig(:,SLOWWRGHT)
    leig(:,SLOWWLEFT) = 0.5*a2inv*leig(:,SLOWWLEFT)

    !scale Alfven waves by 1/2
    leig(:,ALFVNLEFT) = 0.5*leig(:,ALFVNLEFT)
    leig(:,ALFVNRGHT) = 0.5*leig(:,ALFVNRGHT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!right eigenvectors!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !fast waves
    reig(DENS_VAR,SHOCKLEFT) =  A_f*V(DENS_VAR)
    reig(VELX_VAR,SHOCKLEFT) = -A_f*Cf
    reig(VELY_VAR,SHOCKLEFT) =  a_s1*beta(1)
    reig(VELZ_VAR,SHOCKLEFT) =  a_s1*beta(2)
    reig(MAGY_VAR,SHOCKLEFT) =  a_s3*beta(1)
    reig(MAGZ_VAR,SHOCKLEFT) =  a_s3*beta(2)
    reig(PRES_VAR,SHOCKLEFT) =  A_f*V(DENS_VAR)*a2

    reig(DENS_VAR         ,SHOCKRGHT) =  reig(DENS_VAR         ,SHOCKLEFT)
    reig(VELX_VAR:VELZ_VAR,SHOCKRGHT) = -reig(VELX_VAR:VELZ_VAR,SHOCKLEFT)
    reig(MAGY_VAR:PRES_VAR,SHOCKRGHT) =  reig(MAGY_VAR:PRES_VAR,SHOCKLEFT)

    !slow waves 
    reig(DENS_VAR,SLOWWLEFT) =  A_s*V(DENS_VAR)
    reig(VELX_VAR,SLOWWLEFT) = -A_s*Cs
    reig(VELY_VAR,SLOWWLEFT) = -a_f1*beta(1)
    reig(VELZ_VAR,SLOWWLEFT) = -a_f1*beta(2)
    reig(MAGY_VAR,SLOWWLEFT) = -a_f3*beta(1)
    reig(MAGZ_VAR,SLOWWLEFT) = -a_f3*beta(2)
    reig(PRES_VAR,SLOWWLEFT) = A_s*V(DENS_VAR)*a2

    reig(DENS_VAR         ,SLOWWRGHT) =  reig(DENS_VAR         ,SLOWWLEFT)
    reig(VELX_VAR:VELZ_VAR,SLOWWRGHT) = -reig(VELX_VAR:VELZ_VAR,SLOWWLEFT)
    reig(MAGY_VAR:PRES_VAR,SLOWWRGHT) =  reig(MAGY_VAR:PRES_VAR,SLOWWLEFT)

    !alfven waves
    !reig(DENS_VAR,ALFVNLEFT) = 0.
    !reig(VELX_VAR,ALFVNLEFT) = 0.
    reig(VELY_VAR,ALFVNLEFT) = -beta(2)
    reig(VELZ_VAR,ALFVNLEFT) =  beta(1)
    reig(MAGY_VAR,ALFVNLEFT) = -beta(2)*sqrtd
    reig(MAGZ_VAR,ALFVNLEFT) =  beta(1)*sqrtd
    !reig(PRES_VAR,ALFVNLEFT) = 0.

    reig(VELY_VAR:VELZ_VAR,ALFVNRGHT) = -reig(VELY_VAR:VELZ_VAR,ALFVNLEFT)
    reig(MAGY_VAR:MAGZ_VAR,ALFVNRGHT) =  reig(MAGY_VAR:MAGZ_VAR,ALFVNLEFT)

    reig(DENS_VAR,CTENTROPY) = 1.
    !reig(VELX_VAR,CTENTROPY) = 0.
    !reig(VELY_VAR,CTENTROPY) = 0.
    !reig(VELZ_VAR,CTENTROPY) = 0.
    !reig(MAGY_VAR,CTENTROPY) = 0.
    !reig(MAGZ_VAR,CTENTROPY) = 0.
    !reig(PRES_VAR,CTENTROPY) = 0.




!normalization check
!!$    do varL = 1, NUMB_WAVE
!!$       dot = dot_product(leig(:,varL),reig(:,varL))
!!$       if (abs(dot-1.) > 1.e-4) then
!!$          print *, 'not unity', dot, varL
!!$       end if
!!$    end do
!!$
!!$    do varL = 1, NUMB_WAVE
!!$       if (varL < NUMB_WAVE) then
!!$          dot = dot_product(leig(:,varL),reig(:,varL+1))
!!$          if (abs(dot) > 1.e-4) then
!!$             print *, 'not 0', dot, varL, varL+1
!!$          end if
!!$       else
!!$          dot = dot_product(leig(:,varL),reig(:,varL-1))
!!$          if (abs(dot) > 1.e-4 ) then
!!$             print *, 'not 0',dot, varL, varL-1
!!$          end if
!!$       end if
!!$    end do

    if (conservative) then
       !convert from eigenvectors in primitive vars to cons. vars using jacobian and inverse

       JacQ  = 0.
       JacQi = 0.

       !make Q
       JacQ (DENS_VAR,DENS_VAR)          = 1.
       
       JacQ (VELX_VAR,DENS_VAR)          = V(VELX_VAR)
       JacQ (VELX_VAR,VELX_VAR)          = V(DENS_VAR)
       
       JacQ (VELY_VAR,DENS_VAR)          = V(VELY_VAR)
       JacQ (VELY_VAR,VELY_VAR)          = V(DENS_VAR)
       
       JacQ (VELZ_VAR,DENS_VAR)          = V(VELZ_VAR)
       JacQ (VELZ_VAR,VELZ_VAR)          = V(DENS_VAR)
       
       JacQ (MAGY_VAR,MAGY_VAR)          = 1.
       JacQ (MAGZ_VAR,MAGZ_VAR)          = 1.
       
       JacQ (PRES_VAR,DENS_VAR:PRES_VAR) =(/0.5*u2, V(DENS_VAR)*V(VELX_VAR),&
                                                    V(DENS_VAR)*V(VELY_VAR), &
                                                    V(DENS_VAR)*V(VELZ_VAR),&
                                                    V(MAGY_VAR), &
                                                    V(MAGZ_VAR), &
                                                    -1./K /)

       !make Q inverse
       JacQi(DENS_VAR,DENS_VAR)          = 1.
       JacQi(DENS_VAR,VELX_VAR)          = -V(VELX_VAR)*dinv
       JacQi(DENS_VAR,VELY_VAR)          = -V(VELY_VAR)*dinv
       JacQi(DENS_VAR,VELZ_VAR)          = -V(VELZ_VAR)*dinv
       
       JacQi(VELX_VAR,VELX_VAR)          = dinv
       JacQi(VELY_VAR,VELY_VAR)          = dinv
       JacQi(VELZ_VAR,VELZ_VAR)          = dinv
       
       JacQi(MAGY_VAR,MAGY_VAR)          = 1.
       JacQi(MAGZ_VAR,MAGZ_VAR)          = 1.
       
       JacQi(DENS_VAR:PRES_VAR,PRES_VAR) =(/-0.5*u2,V(VELX_VAR),&
                                                    V(VELY_VAR), &
                                                    V(VELZ_VAR),&
                                                    V(MAGY_VAR), &
                                                    V(MAGZ_VAR), &
                                                    -1. /)*k

       !right eigenvectors
       do jj = 1, NUMB_WAVE
          do ii = 1, NSYS_VAR
             ReigTemp(ii,jj) = dot_product(JacQ(ii,:), reig(:,jj))
          end do
       end do
       reig = ReigTemp

       !left eigenvectors
       do jj = 1, NUMB_WAVE
          do ii = 1, NSYS_VAR
             LeigTemp(ii,jj) = dot_product(leig(:,jj),JacQi(ii,:))
          end do
       end do
       leig = LeigTemp



    end if


    return
  end subroutine mhd_eigenvectors


  
  


 

  
end module eigensystem
