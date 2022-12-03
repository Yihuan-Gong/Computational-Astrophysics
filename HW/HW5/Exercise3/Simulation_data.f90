module Simulation_data
    implicit none
    integer, parameter :: imax   = 128   ! number of points in the x and y direction
    integer, parameter :: ibuf   = 1     ! number of ghost zones for BC.
    integer, parameter :: istart = 1     ! starting point
    integer, parameter :: iend   = imax  ! end point

    real, parameter  :: cx       = 1.0   ! velocity x comonent
    real, parameter  :: cy       = 1.0   ! velocity y comonent
    real, parameter  :: xmin     = 0.0   ! left position
    real, parameter  :: xmax     = 1.0   ! right position
    real, parameter  :: ymin     = 0.0   ! bottom position
    real, parameter  :: ymax     = 1.0   ! top position
    real, parameter  :: tend     = 2.5   ! final time

    real, parameter  :: cfl      = 0.4   ! cfl number
    real, save       :: dx, dy

    real, dimension(istart-ibuf:iend+ibuf), save :: x, y
    real, dimension(istart-ibuf:iend+ibuf, istart-ibuf:iend+ibuf), save :: u, uold

    integer,parameter  :: io_interval = 10

end module Simulation_data
