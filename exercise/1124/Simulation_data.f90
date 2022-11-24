module Simulation_data
    implicit none
    integer, parameter :: imax   = 500   ! number of points in the x direction
    integer, parameter :: ibuf   = 1     ! number of ghost zones for BC.
    integer, parameter :: istart = 1     ! starting point
    integer, parameter :: iend   = imax  ! end point

    real, parameter  :: c        = 1.0   ! velocity
    real, parameter  :: xmin     = 0.0   ! left position
    real, parameter  :: xmax     = 1.0   ! right position
    real, parameter  :: tend     = 0.7   ! final time

    real, parameter  :: cfl      = 0.4   ! cfl number
    real, save       :: dx

    real, dimension(istart-ibuf:iend+ibuf), save :: u, uold, x

    integer,parameter  :: io_interval = 10

end module Simulation_data
