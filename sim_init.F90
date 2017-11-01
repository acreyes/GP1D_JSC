subroutine sim_init()

#include "definition.h"

  use sim_data
  use read_initFile

  
  implicit none

  call read_initFileInt ('slug.init','sim_order',  sim_order)
  call read_initFileInt ('slug.init','sim_Torder',  sim_Torder)
  call read_initFileInt ('slug.init','sim_nstep',  sim_nStep)
  call read_initFileReal('slug.init','sim_dt',    sim_dt)
  call read_initFileReal('slug.init','sim_cfl',    sim_cfl)
  call read_initFileReal('slug.init','sim_tmax',   sim_tmax)
  call read_initFileReal('slug.init','sim_WENeps',   sim_WENeps)
  call read_initFileReal('slug.init','sim_outputIntervalTime',sim_outputIntervalTime)
  call read_initFileChar('slug.init','sim_riemann',sim_riemann)
  call read_initFileChar('slug.init','sim_limiter',sim_limiter)
  call read_initFileChar('slug.init','sim_name',sim_name)
  call read_initFileChar('slug.init','sim_WENO',sim_WENO)
  call read_initFileBool('slug.init','sim_charLimiting',sim_charLimiting)
  call read_initFileBool('slug.init','sim_RK',sim_RK)
  call read_initFileBool('slug.init','sim_fixDt',sim_fixDt)
  call read_initFileBool('slug.init','sim_nlim',sim_nlim)
  call read_initFileBool('slug.init','sim_shockDetect',sim_shockDetect)
  call read_initFileBool('slug.init','sim_upwindTVD',sim_upwindTVD)
  call read_initFileChar('slug.init','sim_mono',sim_mono)

  call read_initFileChar('slug.init','sim_icType',   sim_icType)
  call read_initFileReal('slug.init','sim_Bx'    ,   sim_Bx)
  call read_initFileReal('slug.init','sim_densL',    sim_densL)
  call read_initFileReal('slug.init','sim_velxL',    sim_velxL)
  call read_initFileReal('slug.init','sim_velyL',    sim_velyL)
  call read_initFileReal('slug.init','sim_velzL',    sim_velzL)
  call read_initFileReal('slug.init','sim_presL',    sim_presL)
  call read_initFileReal('slug.init','sim_enerL',    sim_enerL)
  call read_initFileReal('slug.init','sim_magyL',    sim_magyL)
  call read_initFileReal('slug.init','sim_magzL',    sim_magzL)
  call read_initFileReal('slug.init','sim_densR',    sim_densR)
  call read_initFileReal('slug.init','sim_velxR',    sim_velxR)
  call read_initFileReal('slug.init','sim_velyR',    sim_velyR)
  call read_initFileReal('slug.init','sim_velzR',    sim_velzR)
  call read_initFileReal('slug.init','sim_enerR',    sim_enerR)
  call read_initFileReal('slug.init','sim_presR',    sim_presR)
  call read_initFileReal('slug.init','sim_magyR',    sim_magyR)
  call read_initFileReal('slug.init','sim_magzR',    sim_magzR)
  call read_initFileReal('slug.init','sim_gamma',    sim_gamma)
  call read_initFileReal('slug.init','sim_shockLoc', sim_shockLoc)
  call read_initFileReal('slug.init','sim_smallPres', sim_smallPres)

  call read_initFileChar('slug.init','sim_bcType',sim_bcType)

  call read_initFileReal('slug.init','sim_ioTfreq',  sim_ioTfreq)
  call read_initFileInt ('slug.init','sim_ioNfreq',  sim_ioNfreq)
  !call read_initFileInt ('slug.init','sim_mval',  sim_mval)
  call read_initFileReal ('slug.init','sim_mval',  sim_mval)

  call read_initFileChar('slug.init','sim_quad'     ,  sim_quad)
  call read_initFileChar('slug.init','sim_gp_kernel',  sim_gp_kernel)
  call read_initFileReal('slug.init','sim_sigma'    ,  sim_sigma)
  call read_initFileReal('slug.init','sim_sigdel'   ,  sim_sigdel)
  call read_initFileReal('slug.init','sim_matern_nu',  sim_matern_nu)
  call read_initFileReal('slug.init','sim_RQ_alpha',  sim_RQ_alpha)

  call sim_initBlock()

end subroutine sim_init
