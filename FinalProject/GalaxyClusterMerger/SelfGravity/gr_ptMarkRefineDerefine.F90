!!****if* source/Simulation/SimulationMain/GalaxyCluster/gr_ptMarkRefineDerefine
!!
!! NAME
!!  gr_ptMarkRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_ptMarkRefineDerefine()
!!
!! DESCRIPTION
!!
!! Routine to mark blocks for refinement or derefinement based upon
!! the count of particles in it.
!!
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptMarkRefineDerefine ()

#include "constants.h"
#include "Flash.h"

  use Particles_interface, ONLY : Particles_getCountPerBlk
  use Grid_data, ONLY : gr_refineOnParticleCount,gr_refineOnPdens,&
       gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  use tree, ONLY : refine,derefine,lrefine,lrefine_min,lrefine_max,stay
  use Grid_interface, ONLY : Grid_getListOfBlocks
  implicit none

  integer,dimension(MAXBLOCKS)::oneBlkCount,blkList
  integer :: i,blkCount,blockID
  real,dimension(:,:,:,:),pointer :: solnData
  logical :: deref

  if(.not.gr_refineOnParticleCount)return

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

  call Particles_getCountPerBlk(oneBlkCount)
  do i=1,blkCount
     blockID=blkList(i)
     refine(blockID)=(oneBlkCount(blockID)>gr_maxParticlesPerBlk) &
          .or. refine(blockID)
     refine(blockID)=refine(blockID).and.(lrefine(blockID)<lrefine_max)
     deref=(oneBlkCount(blockID)<gr_minParticlesPerBlk)
     derefine(blockID) = (.not. refine(blockID)) .and. &
          (lrefine(blockID)>lrefine_min) &
          .and. deref .and. (.not. stay(blockID))
  end do

  return
  
end subroutine gr_ptMarkRefineDerefine
