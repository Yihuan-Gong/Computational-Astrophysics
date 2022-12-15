!!****f* source/Simulation/Simulation_initRestart
!!
!! NAME
!!  Simulation_initRestart
!!
!! SYNOPSIS
!!  Simulation_initRestart()
!!
!! DESCRIPTION
!!  This is where the user should place code for a setup that needs to adjust
!!  data on a restart, particularly if grid data, grid metadata or particle 
!!  data needs to be changed on restarting.
!!
!! ARGUMENTS
!!
!!***

subroutine Simulation_initRestart()

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords, &
       Grid_getListOfBlocks

  use Simulation_data, ONLY : sim_insertBubbles, sim_bubbleMethod, &
       sim_insertBubblesOnRestart, sim_bubbleX, sim_bubbleY, sim_bubbleZ, &
       sim_bubbleRadius, sim_BubbleEntropy, sim_numBubbles, &
       sim_useBubbleParticles, sim_meshMe

  use Eos_interface, ONLY : Eos_wrapped

  use Particles_data, ONLY : pt_numLocal, pt_typeInfo, pt_posInitialized

  use Particles_interface, ONLY : Particles_updateRefinement, Particles_updateAttributes

  use pt_interface, ONLY : pt_initPositions

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
#include "Particles.h"

  integer :: blockCount
  integer, dimension(MAXBLOCKS) :: blockList

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: i, j, k, lb
  integer :: size(3)
  
  logical :: partInit = .true.

  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real, dimension(:), allocatable :: xc, yc, zc

  real, parameter :: mue_twofifths = 1.052463228428472
  real, parameter :: mh = 1.6737352238051868e-24

  real :: bub_rad1, bub_rad2

  if (.not. sim_insertBubbles .and. .not. sim_insertBubblesOnRestart .or. &
       sim_bubbleMethod /= 1) then
     return
  endif

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do lb = 1, blockCount

     call Grid_getBlkPtr(blocklist(lb), solnVec)

     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     
     size(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     size(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     size(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

     allocate(xc(size(1)), yc(size(2)), zc(size(3)))

     call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, .true., zc, size(3))
     call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, .true., yc, size(2))
     call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., xc, size(1))

     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)

        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)

           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

              bub_rad1 = sqrt((xc(i)-sim_bubbleX(1))**2 + &
                   (yc(j)-sim_bubbleY(1))**2 + &
                   (zc(k)-sim_bubbleZ(1))**2)

              if (sim_numBubbles == 2) then
                 bub_rad2 = sqrt((xc(i)-sim_bubbleX(2))**2 + &
                      (yc(j)-sim_bubbleY(2))**2 + &
                      (zc(k)-sim_bubbleZ(2))**2)
              else
                 bub_rad2 = sim_bubbleRadius*100.0
              endif

              solnVec(BUBB_MSCALAR,i,j,k) = 0.0

              if (bub_rad1 <= sim_bubbleRadius .or. bub_rad2 <= sim_bubbleRadius) then
                 solnVec(DENS_VAR,i,j,k) = (solnVec(PRES_VAR,i,j,k) / sim_bubbleEntropy) ** 0.6 * mh * mue_twofifths
                 solnVec(BUBB_MSCALAR,i,j,k) = 1.0
              endif

           enddo

        enddo

     enddo

     call Eos_wrapped(MODE_DENS_PRES, blkLimits, blockList(lb))

     call Grid_releaseBlkPtr(blocklist(lb), solnVec)

     deallocate(xc)
     deallocate(yc)
     deallocate(zc)

     if (sim_useBubbleParticles) then
        call pt_initPositions(blockList(lb), partInit)
     endif
          
  enddo

  if (sim_meshMe == 0) then
     print *, "Done with particle positions"
  endif

  if (sim_useBubbleParticles) then
     pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_numLocal
     pt_posInitialized = .true.
     call Particles_updateRefinement(blockCount)
     call pt_createTag()
     call Particles_updateAttributes()
  endif

  return

end subroutine Simulation_initRestart
