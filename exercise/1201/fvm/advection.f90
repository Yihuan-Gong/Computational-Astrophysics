!===================================================================
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.11.20
!
! The simpliest first order finite-volume method for calculating 
!
! the advection equation: u_t = - c u_x
!
!===================================================================
program advection

    implicit none

    call initial()

    call evolution()

    call finalize()

end program advection
