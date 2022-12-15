!!****if* source/Simulation/SimulationMain/magnetoHD/ClusterSloshingMagneticNew/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2),intent(in):: pos, 
!!                      integer, intent(in) :: sweepDir, 
!!                      integer, intent(in) :: blockID, 
!!                      integer, intent(in) :: numCells, 
!!                      real(numCells),intent(out) :: grav, 
!!                      integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  test against sweep_x,
!!                                 sweep_y, and sweep_z
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, potentialIndex)

  use Simulation_data

  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getCellCoords

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer, intent(IN),optional        :: potentialIndex
  real            :: blockSize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  
  integer         :: ii, iimin, iimax, lb
  real            :: gpot(numCells), delxinv
  real, parameter :: onesixth = 1.e0/6.e0
  integer         :: potVar
  real, dimension(2) :: xCtr, yCtr, zCtr, axCtr, ayCtr, azCtr
  real, dimension(:), allocatable :: inertialAcceleration, rr1, rr2, weightAccel
  real, dimension(:), allocatable :: xc, yc, zc

  real, external :: interpolate

  !==================================================
  
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkPtr(blockID, solnVec)
  
!! IF a variable index is explicitly specified, assume that as the potential
!! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else
     potVar=GPOT_VAR
  end if

  if (potVar == GPOL_VAR) then
     xCtr = sim_oxCtr
     yCtr = sim_oyCtr
     zCtr = sim_ozCtr
     axCtr = sim_oaxCtr
     ayCtr = sim_oayCtr
     azCtr = sim_oazCtr
  else
     xCtr = sim_xCtr
     yCtr = sim_yCtr
     zCtr = sim_zCtr
     axCtr = sim_axCtr
     ayCtr = sim_ayCtr
     azCtr = sim_azCtr
  endif

  iimin   = 1
  iimax   = numCells

  allocate(rr1(numCells))
  allocate(weightAccel(numCells))
  allocate(inertialAcceleration(numCells))
  allocate(xc(numCells))
  allocate(yc(numCells))
  allocate(zc(numCells))

  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zc, numCells)
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., yc, numCells)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xc, numCells)

  if (sweepDir == SWEEP_X) then                     ! x-direction

     rr1(:) = sqrt((xc(:)-xCtr(1))*(xc(:)-xCtr(1)) + &
          (yc(pos(1))-yCtr(1))*(yc(pos(1))-yCtr(1)) + &
          (zc(pos(2))-zCtr(1))*(zc(pos(2))-zCtr(1)))

     inertialAcceleration(:) = -axCtr(1)

     delxinv = real(NXB) / blockSize(IAXIS)
     
     gpot(:) = solnVec(potVar,:,pos(1),pos(2))

  elseif (sweepDir == SWEEP_Y) then                 ! y-direction

     rr1(:) = sqrt((xc(pos(1))-xCtr(1))*(xc(pos(1))-xCtr(1)) + &
          (yc(:)-yCtr(1))*(yc(:)-yCtr(1)) + &
          (zc(pos(2))-zCtr(1))*(zc(pos(2))-zCtr(1)))
     
     inertialAcceleration(:) = -ayCtr(1)

     delxinv = real(NYB) / blockSize(JAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),:,pos(2))

  else                                            ! z-direction

     rr1(:) = sqrt((xc(pos(1))-xCtr(1))*(xc(pos(1))-xCtr(1)) + &
          (yc(pos(2))-yCtr(1))*(yc(pos(2))-yCtr(1)) + &
          (zc(:)-zCtr(1))*(zc(:)-zCtr(1)))

     inertialAcceleration(:) = -azCtr(1)

     delxinv = real(NZB) / blockSize(KAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),pos(2),:)

  endif
  
  !-------------------------------------------------------------------------------
  
  !               Compute gravitational acceleration

  delxinv = 0.5e0 * delxinv
  
  do ii = iimin+1, iimax-1
     grav(ii) = delxinv * (gpot(ii-1) - gpot(ii+1))
  enddo
  
  grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
  grav(iimax) = grav(iimax-1)
  
  if (sim_mainClusterFixed .and. .not. sim_testSingleCluster) then

     where (rr1 <= sim_Rdamp)
        weightAccel = 1.0
     elsewhere
        weightAccel = exp(-(rr1-sim_Rdamp)/sim_Ldamp)
     endwhere

     grav = grav + inertialAcceleration*weightAccel

  endif

  call Grid_releaseBlkPtr(blockID, solnVec)
  
  deallocate(xc, yc, zc, rr1)
  deallocate(weightAccel)
  deallocate(inertialAcceleration)

  return
   
end subroutine Gravity_accelOneRow

