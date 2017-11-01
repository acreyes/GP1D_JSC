subroutine sim_GPinit2()
  !want to initialize all the GP parameters that do not depend on the data, and only on the grid geometries.
  !I am using sig/del = 2
  !This includes calculating the following:
  !****** the covariance matrix C 
  !****** the cross-correlation matrix T 
  !****** The corresponding vectors and matrices to be used in actual calculations (see eqs 30-32) v & Z
  !******  ****** C.v  = u  (gr_GPv)
  !******  ****** C.Z* = T* (gr_GPZ)

  !for the case of general WENO like candidate stencils for the global stencil radius R
  !the global stencil will be of size 2R+1
  !then each of the smaller stencils will be of size R+1 --> 3 for the case of WENO5
  !and there will be R+1 of these smaller stencils (l) so l = 1, R+1
  !each stencil S_l will extend from i-R+(l-1) to i+(l-1) --> R=1 for WENO5 and l
#include "definition.h"

  use grid_data
  use linalg
  use GP

  implicit none
  real, dimension(2*gr_radius+1,2*gr_radius+1) :: C, L, W, Vt, Sigmai
  real, dimension(2,2*gr_radius+1) :: T
  real, dimension(2*gr_radius+1,2) :: B
  real, dimension(2*gr_radius+1)   :: stencil, u, S


  real, dimension(gr_radius+1, gr_radius+1) :: Ck, Lk, Wk, Vtk, Sigmaik
  real, dimension(2,gr_radius+1)  :: Tk
  real, dimension(gr_radius+1, 2) :: Bk
  real, dimension(gr_radius+1)    :: stencilk, uk, Sk

  real, dimension(2*gr_radius+2, gr_radius+1) :: Zmat
  real, dimension(2*gr_radius+2)    :: Zvec
  real, dimension(gr_radius+1)    :: vec, ul, Pk
  real, dimension(2*gr_radius+1)    :: un

  real, dimension(gr_radius+1,gr_radius+1)  :: K_char, L_char
  real, dimension(gr_radius+1)    :: unl, char_stencil
  
  integer :: LR,i,j,N,M, ROW, COL, R
  real :: small, sigma, sigdel
  
  !initialize
  R = gr_radius
  N = 2*R+1
  C = 0.
  T = 0.
  u = 1.

  do i = 1, N
     stencil(i) =  REAL(i - R - 1)     
  end do

  !first thing is to calculate the covariance matrix according to eq. 15
  !since C is symmetric only bother with one side of the diaganol
  do i = 1,N
     do j = 1,N
        C(i,j) = intg_kernel(stencil(i), stencil(j))     
     end do
    T(1:2, i) = intg_predvec(stencil(i))
  end do
  !now we need to solve the linear eqns for v & Z (see eqs 30-32)
  
  
  call chol(C, N, L)
  call solve_Axb(C, gr_GPv, u, L, N)
  call solve_CZT(C, gr_GPZ, T, L, N)

  !this should only be needed for calculating the conditional probability i.e. I haven't found a need for it yet
  gr_U2_16(1) = K(0.5,0.5)   - dot_product(gr_GPZ(1,:), T(1,:))
  gr_U2_16(2) = K(-0.5,-0.5) - dot_product(gr_GPZ(2,:), T(2,:))

  !now lets solve for the gp vars over the local stencils (size = R+1)
  N = R+1
  Ck = 0.
  Tk = 0.
  uk = 1.

  do m = 1,N
     !make the m-th stencil 
     !these are just the R+1 eno stencils
     !this should extend from -R+m-1 to m-1
     do i = 1, N
        stencil(i) = REAL( i-1 - R + m-1 )
     end do
     !now lets calculate the covariance matrix

     do i = 1,N
        do j = 1,N
           Ck(i,j) = intg_kernel(stencil(i), stencil(j))
        end do
        Tk(1:2, i) = intg_predvec(stencil(i))
     end do

     !now we do cholesky
     call chol(Ck, N, Lk)
     call solve_Axb(Ck, gr_GPvk(:,m), uk, Lk, N)
     call solve_CZT(Ck, gr_GPZk(:,:,m), Tk, Lk, N)
     !z-vectors for eigen system
     do i = 1, N
        do j = 1, N
           Pk(j) = quad_cross(stencil(j), stencil(i))
        end do
        call solve_Axb(Ck, gr_gp_Zvecs(:,i,m), Pk(:), Lk, N)
     end do

  end do

  !now we have the Z vectors for all of our stencils
  !lets now compute the linear weights
  !we do this by solving Ax=b, using least squares for the overdetermined system
  !----where x are the linear weights
  !----A is going to be the matrix of Z_k vectors
  !----b is the Z vector for the larger stencil

  !lets first make Z matrix
  !Zmat = 0.
  ROW = 2*R+1
  COL = R+1
  ul = 1.
  un = 1.
  
  do LR = 1, 2
     Zmat = 0.
     do m = 1, COL
        Zmat(m:m+R,m) = gr_GPZk(LR, :, m)
        Zmat(2*R+2,m) = dot_product(gr_GPZk(LR,:,m), ul)
     end do
     
     Zvec(1:ROW) = gr_GPZ(LR,:)
     Zvec(2*R+2)   = dot_product(gr_GPZ(LR,:),un)

     call LSTSQ(ROW+1, COL, Zmat, gr_GP_w(LR,:), Zvec)
  end do

  !last thing to do is calculate K-inverse for ths smaller stencils
  !we need this in order to compute the marginal-likelyhood later
  N = R+1
  sigma  = sim_sigma
  sigdel = sim_sigdel
  
  do i = 1,N
     do j = 1,N
        Ck(i,j) =  intg_kernel(stencil(i), stencil(j))
     end do
  end do
  do i = 1, N
     vec = 0.
     vec(i) = 1.
     call solve_Axb(Ck, gr_GP_Kki(:,i), vec, Lk, N)

  end do

  
end subroutine sim_GPinit2

