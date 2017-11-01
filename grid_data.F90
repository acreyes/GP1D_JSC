module grid_data
  implicit none
  real, allocatable, dimension(:), save :: gr_xCoord
  real, save :: gr_xbeg,gr_xend,gr_dx
  integer, save :: gr_i0,gr_ibeg,gr_iend,gr_imax,gr_ngc,gr_nx,gr_radius

  real, allocatable, dimension(:,:) :: gr_U ! conservative vars
  real, allocatable, dimension(:,:) :: gr_V ! primitive vars
  real, allocatable, dimension(:,:) :: gr_W ! characteristic vars

  real, allocatable, dimension(:,:) :: gr_vL   ! left Riemann states
  real, allocatable, dimension(:,:) :: gr_vR   ! right Riemann states
  real, allocatable, dimension(:,:) :: gr_flux ! fluxes

  real, allocatable, dimension(:,:)   :: gr_eigval ! eigenvalues
  real, allocatable, dimension(:,:,:) :: gr_reigvc ! right eigenvectors
  real, allocatable, dimension(:,:,:) :: gr_leigvc ! left  eigenvectors

  !GP vars
  !quad precision GP vars
  real(KIND=16), allocatable, dimension(:  ) :: gr_GPv
  real(KIND=16), allocatable, dimension(:,:) :: gr_GPZ
  
  real(KIND=16), allocatable, dimension(:,:  ) :: gr_GPvk
  real(KIND=16), allocatable, dimension(:,:,:) :: gr_GPZk, gr_gp_Zvecs

  real(KIND=8), allocatable, dimension(:  ) :: gr_GPv2
  real(KIND=8), allocatable, dimension(:,:) :: gr_GPZ2

  real(KIND=8), allocatable, dimension(:,:  ) :: gr_GPvk2
  real(KIND=8), allocatable, dimension(:,:,:) :: gr_GPZk2, gr_gp_Zvecs2

  real(KIND=16), allocatable, dimension(:) :: gr_U2_16
  real(KIND=8 ), allocatable, dimension(:)  :: gr_U2

  real(KIND=8 ), allocatable, dimension(:, :) :: gr_GP_w2
  real(KIND=16), allocatable, dimension(:, :) :: gr_GP_w

  real(KIND=16), allocatable, dimension(:, :) :: gr_GP_Kki
  real(KIND=8 ), allocatable, dimension(:, :) :: gr_GP_Kki2

  real(KIND=16), allocatable, dimension(:, :) :: gr_GP_char_Kki
  real(KIND=8 ), allocatable, dimension(:, :) :: gr_GP_char_Kki2

  !double precision GP vars

!!$  real(KIND=8), allocatable, dimension(:  ) :: gr_GPv
!!$  real(KIND=8), allocatable, dimension(:,:) :: gr_GPZ
!!$  
!!$  real(KIND=8), allocatable, dimension(:,:  ) :: gr_GPvk
!!$  real(KIND=8), allocatable, dimension(:,:,:) :: gr_GPZk
!!$
!!$  real(KIND=8), allocatable, dimension(:  ) :: gr_GPv2
!!$  real(KIND=8), allocatable, dimension(:,:) :: gr_GPZ2
!!$
!!$  real(KIND=8), allocatable, dimension(:,:  ) :: gr_GPvk2
!!$  real(KIND=8), allocatable, dimension(:,:,:) :: gr_GPZk2
!!$
!!$  real(KIND=8), allocatable, dimension(:) :: gr_U2_16
!!$  real(KIND=8 ), allocatable, dimension(:)  :: gr_U2
!!$
!!$  real(KIND=8), allocatable, dimension(:, :) :: gr_GP_w2
!!$  real(KIND=8 ), allocatable, dimension(:, :) :: gr_GP_w

end module grid_data
