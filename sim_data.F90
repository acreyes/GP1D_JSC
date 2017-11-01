module sim_data

#include "definition.h"

  implicit none

  !! numerics
  real, save :: sim_cfl, sim_tmax, sim_outputIntervalTime, sim_WENeps, sim_dt, sim_mval, sim_Bx
  integer, save :: sim_order, sim_nStep, sim_Torder!, sim_mval
  character(len=MAX_STRING_LENGTH), save :: sim_name, sim_limiter,sim_riemann, sim_WENO, sim_mono
  logical, save :: sim_charLimiting, sim_RK, sim_fixDt, sim_nlim, sim_shockDetect, sim_upwindTVD

  !! ICs
  real, save :: sim_gamma
  real, save :: sim_densL,sim_velxL,sim_velyL,sim_velzL,sim_presL,sim_enerL,sim_magyL,sim_magzL
  real, save :: sim_densR,sim_velxR,sim_velyR,sim_velzR,sim_presR,sim_enerR,sim_magyR,sim_magzR
  real, save :: sim_shockLoc
  real, save :: sim_smallPres
  character(len=MAX_STRING_LENGTH), save :: sim_icType

  !! BCs
  character(len=MAX_STRING_LENGTH), save :: sim_bcType

  !! IO
  integer, save :: sim_ioNfreq
  real,    save :: sim_ioTfreq

  !! GP
  character(len=MAX_STRING_LENGTH), save :: sim_quad, sim_gp_kernel
  real, save :: sim_sigdel, sim_sigma, sim_matern_nu, sim_RQ_alpha
  real, allocatable, dimension(:) :: sim_gp_evals
  real, allocatable, dimension(:,:) :: sim_gp_evecs

end module sim_data
