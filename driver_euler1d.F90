program driver_euler1d

#include "definition.h"
  
  use sim_data
  use grid_data
  use io
  use bc
  use eos, only : eos_all

  implicit none

  real :: t,dt, small_dt, dt_temp
  integer :: nStep,ioCounter,ioTimeFreqCounter,zerodt
  real :: ioCheckTime

  zerodt = 1
  small_dt = 1.e-10
  t = 0.
  nStep = 0
  ioCounter = 0
  ioTimeFreqCounter = 0
  
  
  ! grid_init should be called first before sim_init
  call grid_init()
  call sim_init()
  
  

  write(*,*)''
  write(*,*)'================================================='
  write(*,*)'                GP1D MHD Code                    '
  write(*,*)'      Written by Adam Reyes & Prof. Dongwook Lee              '
  write(*,*)'================================================='
  write(*,*)''

  
  ! write the initial condition
  write(*,*)''
  write(*,*)'       Initial condition was written!            '
  write(*,*)'================================================='
  write(*,*)'   Steps      Time              dt               '
  write(*,*)'================================================='
  write(*,*)''
  call io_writeOutput(t, nStep,ioCounter)

  if (sim_order == 9) then
     !we're doing GP and need to initialize
     if (sim_sigdel == 0.) then
        sim_sigdel = sim_sigma/gr_dx
     else
        sim_sigma = sim_sigdel*gr_dx
     end if
     
     call sim_GPinit2
     !sim_GPinit2 calculates GP parameters to quad precision
     !truncate to double precision for use in the rest of simulation
     gr_GPV2 = gr_GPV
     gr_GPZ2 = gr_GPZ
     gr_GPvk2 = gr_GPvk
     gr_GPZk2 = gr_GPZk
     gr_U2    = gr_U2_16

     gr_GP_w2 = gr_GP_w
     gr_GP_Kki2 = gr_GP_Kki

     gr_GP_w2(1,:) = gr_GP_w2(1,:)/sum(gr_GP_w2(1,:))
     gr_GP_w2(2,:) = gr_GP_w2(2,:)/sum(gr_GP_w2(2,:))
     gr_gp_Zvecs2(:,:,:) = gr_gp_Zvecs(:,:,:)
     call gp_eigens

     
  end if

  do while ( (t < sim_tmax))
     
     if (sim_fixDt) then
        dt = sim_dt
     else
        call cfl(dt)
     end if
     dt_temp = 2.**nstep*small_dt
     if (dt_temp < dt) then
        dt = dt_temp
     end if
     !check to see if there is a reason to stop
     if (  sim_nlim .and. (nStep .ge. sim_nstep)) then
        exit
     elseif ( abs(t - sim_tmax) .le. dt ) then
        dt = abs(t - sim_tmax)
        !exit
     end if

     call soln_ReconEvolveAvg(dt)
     call soln_update(dt)


     ! call BC on primitive vars
     call bc_apply(gr_V)
     
     ! write outputs every ioNfreq cycle or ioTfreq cycle
     ioCheckTime = sim_ioTfreq*real(ioTimeFreqCounter+1)
     if (t-dt < ioCheckTime .and. t>ioCheckTime) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
        ioCounter = ioCounter + 1
        ioTimeFreqCounter = ioTimeFreqCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif

     if (sim_ioNfreq > 0) then
     if (mod(nStep, sim_ioNfreq) == 0) then
        write(*,*)''
        write(*,*)' Output no.',ioCounter+1, 'has been written      '
        write(*,*)'================================================='
        write(*,*)'   Steps      Time              dt               '
        write(*,*)'================================================='
        write(*,*)''
        ioCounter = ioCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif
     endif
     
     ! update your time and step count
     t = t + dt
     nStep = nStep + 1

     write(*,900)nstep,t,dt
     if (dt .le. 0.) then
        !catch for negative time step
        zerodt = 0
        exit
     end if
  enddo


  !! Let's write the final result before exiting
  write(*,*)''
  write(*,*)' Final output no.',ioCounter+1, 'has been written'
  write(*,*)'================================================='
  write(*,*)'        The final tmax has reached, bye!         '
  write(*,*)'================================================='
  write(*,*)''
  call io_writeOutput(t, nStep,ioCounter+1)

  !! finalize and deallocate memories
  call grid_finalize()

900 format(1x,i5,e16.8,1x,e16.8)
  
  return
end program driver_euler1d
