subroutine pt_initPositions(blockID, partPosInitialized)
    
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, &
       pt_posInitialized, pt_meshMe, pt_meshNumProcs
  use Simulation_data, ONLY: sim_numBubbleParticles, sim_useBubbleParticles, &
       pxbub1, pybub1, pzbub1, pxbub2, pybub2, pzbub2, sim_numBubbles
  use Grid_interface, ONLY: Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords

  use Particles_interface, ONLY : Particles_updateAttributes

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
  
  integer, INTENT(in) :: blockID
  logical, INTENT(inout) :: partPosInitialized

  ! Local variables.
  
  integer, save :: p

  integer       :: i
  logical       :: isInBlock

  real, dimension(MDIM) :: blockLower, blockUpper, blockSize, blockCenter

!==============================================================================

  if ((blockID == -1) .or. (.not. sim_useBubbleParticles)) return

  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  blockLower(:) = blockCenter(:) - 0.5*blockSize(:)
  blockUpper(:) = blockCenter(:) + 0.5*blockSize(:)

  ! Particle slot number (incremented and saved between calls)

  p = pt_numLocal
 
  partPosInitialized = .true.

  do i = 1, sim_numBubbleParticles

     ! Check to see if the particle lies within this block.
     isInBlock = (pxbub1(i) >= blockLower(1)) .and. (pxbub1(i) < blockUpper(1))
     isInBlock = isInBlock .and. ((pybub1(i) >= blockLower(2)) .and. (pybub1(i) < blockUpper(2)))
     isInBlock = isInBlock .and. ((pzbub1(i) >= blockLower(3)) .and. (pzbub1(i) < blockUpper(3)))
        
     if (isInBlock) then
        
        p = p + 1
        
        particles(POSX_PART_PROP,p) = pxbub1(i)
        particles(POSY_PART_PROP,p) = pybub1(i)
        particles(POSZ_PART_PROP,p) = pzbub1(i)
        ! For now velocities are zero                                                                                         
        particles(VELX_PART_PROP,p) = 0.0
        particles(VELY_PART_PROP,p) = 0.0
        particles(VELZ_PART_PROP,p) = 0.0
        particles(BLK_PART_PROP,p)  = real(blockID)
        particles(PROC_PART_PROP,p) = real(pt_meshMe)
        
     endif

  enddo

  if (sim_numBubbles == 2) then

     do i = 1, sim_numBubbleParticles

        ! Check to see if the particle lies within this block.  
                                                                 
        isInBlock = (pxbub2(i) >= blockLower(1)) .and. (pxbub2(i) < blockUpper(1))
        isInBlock = isInBlock .and. ((pybub2(i) >= blockLower(2)) .and. (pybub2(i) < blockUpper(2)))
        isInBlock = isInBlock .and. ((pzbub2(i) >= blockLower(3)) .and. (pzbub2(i) < blockUpper(3)))
     
        if (isInBlock) then
           
           p = p + 1

           particles(POSX_PART_PROP,p) = pxbub2(i)
           particles(POSY_PART_PROP,p) = pybub2(i)
           particles(POSZ_PART_PROP,p) = pzbub2(i)
           ! For now velocities are zero
           particles(VELX_PART_PROP,p) = 0.0
           particles(VELY_PART_PROP,p) = 0.0
           particles(VELZ_PART_PROP,p) = 0.0
           particles(BLK_PART_PROP,p)  = real(blockID)
           particles(PROC_PART_PROP,p) = real(pt_meshMe)

        endif

     enddo

  endif

!==============================================================================

! Set the particle database local number of particles.

  pt_numLocal = p

  return

!==============================================================================
  
end subroutine pt_initPositions
