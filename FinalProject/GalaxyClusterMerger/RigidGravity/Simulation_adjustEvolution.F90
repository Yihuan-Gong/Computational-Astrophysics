!!****f* source/Simulation/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
  
  use Simulation_data

  use Grid_interface, ONLY: Grid_getCellCoords, Grid_getDeltas, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_fillGuardCells, &
       Grid_getBlkIndexLimits

  use Eos_interface, ONLY: Eos_wrapped

  use IO_interface, ONLY: IO_setScalar

!! Hack
  use Driver_data, ONLY : dr_dt, dr_dtOld

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "UHD.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real, dimension(MDIM)             :: del

  real, dimension(:), allocatable :: x, y, z

  real :: bub_rad1, bub_rad2, sphere_mass, lsphere_mass, cell_volume, Einj
  integer :: i, j, k, blockID, lb, ierr, myPE
  logical, save :: gcell = .true.
  logical, save :: first_call = .true.
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: csize(3), jet_dir, vel_x, vel_y, vel_z

  logical :: blow_bubbles, fire_jets

  real, parameter :: onethird = 1./3.
  real, parameter :: onesixth = 1./6.

  real :: rr1, rr2, dtOld, dtNew, phi, theta
  real :: xc, yc, zc, atot1, atot2
  real :: wterm, woldterm, rad2, dx_l, dy_l, dz_l
  real :: dx_r, dy_r, dz_r, Mj, Pj, Ej, cell_area

  real, external :: interpolate

  dtOld = dr_dtOld
  dtNew = dr_dt

  call MPI_Comm_rank(MPI_COMM_WORLD, myPE, ierr)

  ! First update the position of the cluster(s)

  if (.not. sim_testSingleCluster) then

     ! The subcluster travels on a "Keplerian" orbit until it leaves
     ! the other boundary, after which it goes on a constant-velocity
     ! trajectory

     xc = sim_xCtr(2) - sim_xCtr(1)
     yc = sim_yCtr(2) - sim_yCtr(1)
     zc = sim_zCtr(2) - sim_zCtr(1)

     rr1 = sqrt(xc*xc + yc*yc + zc*zc)

     if (rr1 <= sim_rMax1) then
        atot1 = -interpolate(grav1,numPoints1,r1,rr1)
     else
        atot1 = -sim_Newton*sim_totMass1/(rr1*rr1)
     endif
     if (rr1 <= sim_rMax2) then
        atot2 = -interpolate(grav2,numPoints2,r2,rr1)
     else
        atot2 = -sim_Newton*sim_totMass2/(rr1*rr1)
     endif

     sim_axCtr(2) = atot1*(xc/rr1)
     sim_ayCtr(2) = atot1*(yc/rr1)
     sim_azCtr(2) = atot1*(zc/rr1)

     sim_axCtr(1) = atot2*(-xc/rr1)
     sim_ayCtr(1) = atot2*(-yc/rr1)
     sim_azCtr(1) = atot2*(-zc/rr1)

     if (myPE == 0) then
        open(9, file='sub_trajectory.dat', position='append')
        write(9, '(10(E14.7, T1))') stime, sim_xCtr(2), sim_yCtr(2), sim_zCtr(2), sim_vxCtr(2), &
             sim_vyCtr(2), sim_vzCtr(2), sim_axCtr(2), sim_ayCtr(2), sim_azCtr(2)
        close(9)
     endif

     if (first_call .and. ((.not. sim_restart) .or. sim_setupSubclusterWithRestart)) then
        wterm = 0.5*dtNew
        woldterm = 0.
        first_call = .false.
     else
        wterm = 0.5*dtNew + onethird*dtOld + onesixth*dtNew**2/dtOld
        woldterm = onesixth*(dtOld**2 - dtNew**2)/dtOld
     endif

     sim_vxCtr(2) = sim_vxCtr(2) + wterm*sim_axCtr(2) + woldterm*sim_oaxCtr(2)
     sim_vyCtr(2) = sim_vyCtr(2) + wterm*sim_ayCtr(2) + woldterm*sim_oayCtr(2)
     sim_vzCtr(2) = sim_vzCtr(2) + wterm*sim_azCtr(2) + woldterm*sim_oazCtr(2)

     sim_oxCtr(2) = sim_xCtr(2)
     sim_oyCtr(2) = sim_yCtr(2)
     sim_ozCtr(2) = sim_zCtr(2)

     sim_xCtr(2) = sim_xCtr(2) + dtNew*sim_vxCtr(2)
     sim_yCtr(2) = sim_yCtr(2) + dtNew*sim_vyCtr(2)
     sim_zCtr(2) = sim_zCtr(2) + dtNew*sim_vzCtr(2)

     call IO_setScalar("cluster 2 x", sim_xCtr(2))
     call IO_setScalar("cluster 2 y", sim_yCtr(2))
     call IO_setScalar("cluster 2 z", sim_zCtr(2))
     call IO_setScalar("cluster 2 vx", sim_vxCtr(2))
     call IO_setScalar("cluster 2 vy", sim_vyCtr(2))
     call IO_setScalar("cluster 2 vz", sim_vzCtr(2))
     call IO_setScalar("cluster 2 ax", sim_axCtr(2))
     call IO_setScalar("cluster 2 ay", sim_ayCtr(2))
     call IO_setScalar("cluster 2 az", sim_azCtr(2))

     sim_oaxCtr(2) = sim_axCtr(2)
     sim_oayCtr(2) = sim_ayCtr(2)
     sim_oazCtr(2) = sim_azCtr(2)

     if (.not. sim_mainClusterFixed) then

        if (myPE == 0) then
           open(9, file='main_trajectory.dat', position='append')
           write(9, '(10(E14.7, T1))') stime, sim_xCtr(1), sim_yCtr(1), sim_zCtr(1), sim_vxCtr(1), &
                sim_vyCtr(1), sim_vzCtr(1), sim_axCtr(1), sim_ayCtr(1), sim_azCtr(1)
           close(9)
        endif
        
        sim_vxCtr(1) = sim_vxCtr(1) + wterm*sim_axCtr(1) + woldterm*sim_oaxCtr(1)
        sim_vyCtr(1) = sim_vyCtr(1) + wterm*sim_ayCtr(1) + woldterm*sim_oayCtr(1)
        sim_vzCtr(1) = sim_vzCtr(1) + wterm*sim_azCtr(1) + woldterm*sim_oazCtr(1)

        sim_oxCtr(1) = sim_xCtr(1)
        sim_oyCtr(1) = sim_yCtr(1)
        sim_ozCtr(1) = sim_zCtr(1)

        sim_xCtr(1) = sim_xCtr(1) + dtNew*sim_vxCtr(1)
        sim_yCtr(1) = sim_yCtr(1) + dtNew*sim_vyCtr(1)
        sim_zCtr(1) = sim_zCtr(1) + dtNew*sim_vzCtr(1)

        call IO_setScalar("cluster 1 x", sim_xCtr(1))
        call IO_setScalar("cluster 1 y", sim_yCtr(1))
        call IO_setScalar("cluster 1 z", sim_zCtr(1))
        call IO_setScalar("cluster 1 vx", sim_vxCtr(1))
        call IO_setScalar("cluster 1 vy", sim_vyCtr(1))
        call IO_setScalar("cluster 1 vz", sim_vzCtr(1))
        call IO_setScalar("cluster 1 ax", sim_axCtr(1))
        call IO_setScalar("cluster 1 ay", sim_ayCtr(1))
        call IO_setScalar("cluster 1 az", sim_azCtr(1))

        sim_oaxCtr(1) = sim_axCtr(1)
        sim_oayCtr(1) = sim_ayCtr(1)
        sim_oazCtr(1) = sim_azCtr(1)

     endif

  endif

  fire_jets = sim_useJets
  fire_jets = sim_useJets .and. stime > sim_jetStart .and. stime < sim_jetStart + sim_jetDuration

  if (fire_jets) then

     theta = sim_jetAngle
     phi = 2*PI*(stime-sim_jetStart)/sim_jetPeriod
     
     do lb = 1, blkcnt

        blockID = blklst(lb)
        
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        
        csize(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
        csize(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
        csize(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

        allocate(x(csize(1)), y(csize(2)), z(csize(3)))
        
        call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, csize(1))
        call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, csize(2))
        call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, csize(3))
     
        call Grid_getDeltas(blockID, del)
        call Grid_getBlkPtr(blockID, solnData)

        cell_volume = del(1)*del(2)*del(3)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 select case(sim_jetAxis) 
                 case(DIR_X)
                    dx_l = x(i)-(sim_xCtr(1)-sim_jetHeight)
		    dx_r = x(i)-(sim_xCtr(1)+sim_jetHeight)
                    if (abs(dx_l) < 0.5*del(DIR_X)) then
                       jet_dir = -1
                    else if (abs(dx_r) < 0.5*del(DIR_X)) then
                       jet_dir = 1
                    else
                       jet_dir = 0
                    endif
                    vel_x = VELY_VAR
                    vel_y = VELZ_VAR
                    vel_z = VELX_VAR
                    rad2 = (y(j)-sim_yCtr(1))**2+(z(k)-sim_zCtr(1))**2
                    cell_area = del(2)*del(3) 
                 case(DIR_Y)
                    dy_l = y(j)-(sim_yCtr(1)-sim_jetHeight)
		    dy_r = y(j)-(sim_yCtr(1)+sim_jetHeight)
                    if (abs(dy_l) < 0.5*del(DIR_Y)) then
                       jet_dir = -1
                    else if (abs(dy_r) < 0.5*del(DIR_Y)) then
                       jet_dir = 1
                    else
                       jet_dir = 0
                    endif
                    vel_x = VELX_VAR
                    vel_y = VELZ_VAR
                    vel_z = VELY_VAR
                    rad2 = (x(i)-sim_xCtr(1))**2+(z(k)-sim_zCtr(1))**2
                    cell_area = del(1)*del(3)
                 case(DIR_Z)
                    dz_l = z(k)-(sim_zCtr(1)-sim_jetHeight)
		    dz_r = z(k)-(sim_zCtr(1)+sim_jetHeight)
                    if (abs(dz_l) < 0.5*del(DIR_Z)) then
                       jet_dir = -1
                    else if (abs(dz_r) < 0.5*del(DIR_Z)) then
                       jet_dir = 1
                    else
                       jet_dir = 0
                    endif
		    vel_x = VELX_VAR
		    vel_y = VELY_VAR
                    vel_z = VELZ_VAR
                    rad2 = (x(i)-sim_xCtr(1))**2+(y(j)-sim_yCtr(1))**2
                    cell_area = del(1)*del(2)
                 end select

                 if (jet_dir /= 0 .and. rad2 < sim_jetRadius**2) then

                    Mj = sim_jetDensity*sim_jetVelocity*dt*cell_area/cell_volume
                    Pj = jet_dir*Mj*sim_jetVelocity
                    Ej = Mj*sim_jetEnergy

                    solnData(vel_x,i,j,k) = solnData(vel_x,i,j,k)*solnData(DENS_VAR,i,j,k) + &
                         Pj*sin(theta)*cos(phi)
                    solnData(vel_y,i,j,k) = solnData(vel_y,i,j,k)*solnData(DENS_VAR,i,j,k) + &
                         Pj*sin(theta)*sin(phi)
                    solnData(vel_z,i,j,k) = solnData(vel_z,i,j,k)*solnData(DENS_VAR,i,j,k) + &
                         Pj*cos(theta)
                    solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + &
                         Ej
                    solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)+Mj

                    solnData(vel_x,i,j,k) = solnData(vel_x,i,j,k) / solnData(DENS_VAR,i,j,k)
                    solnData(vel_y,i,j,k) = solnData(vel_y,i,j,k) / solnData(DENS_VAR,i,j,k)
                    solnData(vel_z,i,j,k) = solnData(vel_z,i,j,k) / solnData(DENS_VAR,i,j,k)
                    solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) / solnData(DENS_VAR,i,j,k)
                    solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + & 
                                     0.5*(solnData(VELX_VAR,i,j,k)**2. + &
                                          solnData(VELY_VAR,i,j,k)**2. + &
                                          solnData(VELZ_VAR,i,j,k)**2.) 

                 endif

              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID, solnData)
        
        deallocate(x)
        deallocate(y)
        deallocate(z)
        
     enddo

  endif

  blow_bubbles = sim_insertBubbles
  blow_bubbles = blow_bubbles .and. sim_bubbleMethod == 2
  blow_bubbles = stime > sim_bubbleStart .and. stime < sim_bubbleStart + sim_bubbleDuration .and. blow_bubbles

  if (.not. blow_bubbles) return

  lsphere_mass = 0.0

  do lb = 1, blkcnt

     blockID = blklst(lb)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     csize(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     csize(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     csize(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

     allocate(x(csize(1)), y(csize(2)), z(csize(3)))

     call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, csize(1))
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, csize(2))
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, csize(3))

     call Grid_getDeltas(blockID, del)

     call Grid_getBlkPtr(blockID, solnData)

     cell_volume = del(1)*del(2)*del(3)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              bub_rad1 = sqrt((x(i)-sim_bubbleX(1))**2 + &
                   (y(j)-sim_bubbleY(1))**2 + &
                   (z(k)-sim_bubbleZ(1))**2)
              bub_rad2 = sqrt((x(i)-sim_bubbleX(2))**2 + &
                   (y(j)-sim_bubbleY(2))**2 + &
                   (z(k)-sim_bubbleZ(2))**2)

              if (bub_rad1 <= sim_bubbleRadius .or. bub_rad2 <= sim_bubbleRadius) then	
                 lsphere_mass = lsphere_mass + solnData(DENS_VAR,i,j,k)*cell_volume
              endif

           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID, solnData)

     deallocate(x)
     deallocate(y)
     deallocate(z)

  enddo

  call mpi_allreduce(lsphere_mass, sphere_mass, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  do lb = 1, blkcnt

     blockID = blklst(lb)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     csize(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     csize(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     csize(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

     allocate(x(csize(1)), y(csize(2)), z(csize(3)))

     call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, csize(1))
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, csize(2))
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, csize(3))

     call Grid_getBlkPtr(blockID, solnData)
  
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              bub_rad1 = sqrt((x(i)-sim_bubbleX(1))**2 + &
                   (y(j)-sim_bubbleY(1))**2 + &
                   (z(k)-sim_bubbleZ(1))**2)
              bub_rad2 = sqrt((x(i)-sim_bubbleX(2))**2 + &
                   (y(j)-sim_bubbleY(2))**2 + &
                   (z(k)-sim_bubbleZ(2))**2)

              if (bub_rad1 <= sim_bubbleRadius .or. bub_rad2 <= sim_bubbleRadius) then
                               
                 Einj = sim_bubblePower*dt/sphere_mass

                 solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + Einj

#ifdef ENER_VAR
                 solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + & 
                      0.5*(solnData(VELX_VAR,i,j,k)**2. + &
                      solnData(VELY_VAR,i,j,k)**2. + &
                      solnData(VELZ_VAR,i,j,k)**2.)
#endif

              endif
                 
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID, solnData)

     deallocate(x)
     deallocate(y)
     deallocate(z)

     call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)

  enddo

  call Grid_fillGuardCells(CENTER, ALLDIR, eosMode=MODE_DENS_EI, doEos=.true.)

end subroutine Simulation_adjustEvolution
