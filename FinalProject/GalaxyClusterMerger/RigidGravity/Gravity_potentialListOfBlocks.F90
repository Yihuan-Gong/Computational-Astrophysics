!!****if* source/Simulation/SimulationMain/magnetoHD/ClusterSloshingMagneticNew/Gravity_potentialListOfBlocks
!!
!!  NAME 
!!
!!     Gravity_potentialListOfBlocks
!!
!!  SYNOPSIS
!!
!!  Gravity_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                integer(IN) :: blockList(blockCount))
!!
!!  DESCRIPTION 
!!      This routine computes the gravitational potential for the gravity
!!      implementations (i.e., various Poisson implementations) which make
!!      use of it in computing the gravitational acceleration.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!
!!
!!
!!***

subroutine Gravity_potentialListOfBlocks(blockCount,blockList)

  use Gravity_data, ONLY : useGravity
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_solvePoisson, &
       GRID_PDE_BND_ISOLATED

  ! Hacks (shouldn't use this directly)

  use Simulation_data

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real, dimension(:), allocatable :: x, y, z

  integer       :: ierr

  integer       :: size(3)
  integer       :: lb, i, j, k
  real          :: rr1, rr2
  real          :: xc, yc, zc, atot1, atot2
  integer       :: bcTypes(6) = GRID_PDE_BND_ISOLATED
  real          :: bcValues(2,6) = 0.
  integer :: error
  character(len=124) :: str_buffer

  real, external :: interpolate

  !=========================================================================

  if(.not.useGravity) return
  
  call Timers_start("gravity Barrier")
  call MPI_Barrier (MPI_COMM_WORLD, ierr)
  call Timers_stop("gravity Barrier")

  call Timers_start("gravity")
 
  do lb = 1, blockCount

     call Grid_getBlkPtr(blocklist(lb), solnVec)

     solnVec(GPOL_VAR,:,:,:) = solnVec(GPOT_VAR,:,:,:)

     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     
     size(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     size(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     size(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

     allocate(x(size(1)), y(size(2)), z(size(3)))

     call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, .true., z, size(3))
     call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, .true., y, size(2))
     call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., x, size(1))

     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)

        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)

           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

              xc = x(i) - sim_xCtr(1)
              yc = y(j) - sim_yCtr(1)
              zc = z(k) - sim_zCtr(1)

              rr1 = sqrt(xc*xc + yc*yc + zc*zc)
              
              if (.not. sim_testSingleCluster) then

                 rr2 = sqrt((x(i)-sim_xCtr(2))**2 + &
                      (y(j)-sim_yCtr(2))**2 + &
                      (z(k)-sim_zCtr(2))**2)
                 
              endif
              
              if (rr1 <= sim_rMax1) then
#ifdef PDEN_VAR
                 solnVec(PDEN_VAR,i,j,k) = interpolate(pden1, numPoints1, r1, rr1)
#endif
                 solnVec(GPOT_VAR,i,j,k) = -interpolate(gpot1, numPoints1, r1, rr1)
              else
                 solnVec(GPOT_VAR,i,j,k) = -sim_Newton*sim_totMass1/rr1
              endif

              if (.not. sim_testSingleCluster) then
                 if (rr2 <= sim_rMax2) then
#ifdef PDEN_VAR
                    solnVec(PDEN_VAR,i,j,k) = solnVec(PDEN_VAR,i,j,k) + &
                         interpolate(pden2, numPoints2, r2, rr2)
#endif
                    solnVec(GPOT_VAR,i,j,k) = solnVec(GPOT_VAR,i,j,k) - &
                         interpolate(gpot2, numPoints2, r2, rr2)
                 else
                    solnVec(GPOT_VAR,i,j,k) = solnVec(GPOT_VAR,i,j,k) - &
                         sim_Newton*sim_totMass2/rr2
                 endif
              endif
                 
           enddo

        enddo

     enddo

     call Grid_releaseBlkPtr(blocklist(lb), solnVec)

     deallocate(x)
     deallocate(y)
     deallocate(z)
     
  enddo

#ifdef USEBARS
  call MPI_Barrier (MPI_Comm_World, ierr)
#endif  
  call Timers_stop ("gravity")
  
  return
end subroutine Gravity_potentialListOfBlocks
