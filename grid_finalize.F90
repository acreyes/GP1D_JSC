subroutine grid_finalize()

  use grid_data!, only : gr_xCoord, gr_U, gr_V, gr_W, gr_eigval, gr_leigvc, gr_reigvc

  implicit none

  if (allocated(gr_xCoord) .eqv. .true.) deallocate(gr_xCoord)
  if (allocated(gr_U) .eqv. .true.) deallocate(gr_U)
  if (allocated(gr_V) .eqv. .true.) deallocate(gr_V)
  if (allocated(gr_W) .eqv. .true.) deallocate(gr_W)

  if (allocated(gr_vR) .eqv. .true.) deallocate(gr_vR)
  if (allocated(gr_vL) .eqv. .true.) deallocate(gr_vL)
  if (allocated(gr_flux) .eqv. .true.) deallocate(gr_flux)

  if (allocated(gr_eigval) .eqv. .true.) deallocate(gr_eigval)
  if (allocated(gr_leigvc) .eqv. .true.) deallocate(gr_leigvc)
  if (allocated(gr_reigvc) .eqv. .true.) deallocate(gr_reigvc)

  
  return
end subroutine grid_finalize
