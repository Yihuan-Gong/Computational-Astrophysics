!!****if* source/Simulation/SimulationMain/ClusterSloshing/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical, INTENT(inout) :: partPosInitialized,
!!                             logical, INTENT(out) :: updateRefine)
!!
!!
!! DESCRIPTION
!!
!!    Initialize particle locations. This routine calls pt_initPositions
!!    which a Particles unit's local routine to initialize the positions
!!    on leaf blocks. The routine also creates tags for all the particles
!!    This routine will initialize based on Lattice or with Density 
!!    distribution  depending upon which of the two is selected. 
!!
!! ARGUMENTS
!!
!!  partPosInitialized : boolean indicating whether particles positions were 
!!            successfully initialized. This is not really relevant
!!            for this version of the routine
!!
!! updateRefine : is true if the routine wished to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refine.
!!
!!***


subroutine Particles_initPositions (partPosInitialized,updateRefine)

  use Driver_interface, ONLY : Driver_abortFlash
  use Particles_interface, ONLY : Particles_updateRefinement
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getLocalNumBlks
  use pt_interface, ONLY : pt_initPositions,pt_createTag 
  use Particles_data, ONLY : pt_posInitialized,pt_numLocal, useParticles,&
       pt_typeInfo, particles, pt_maxPerProc, pt_meshMe

  use pt_interface, ONLY :  pt_initLocal

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "Flash_mpi.h"

  logical, INTENT(INOUT) :: partPosInitialized
  logical, INTENT(OUT) :: updateRefine

  integer       :: i, j, k, b
  integer       :: p, ierr
  integer       :: numLocalThisType, numNewLocalThisType, numLocalPreviousTypes
  integer       :: numPreviousLocal, localNumBlocks
  integer       :: blockID
  logical       :: IsInBlock
  real          :: xpos, ypos, zpos, bxl, byl, bzl, bxu, byu, bzu
  real          :: xvel, yvel, zvel

! NOTE dxParticle is particle spacing, not grid spacing
  real, dimension(MDIM) :: dxParticle = 0.0
  real, dimension(2,MDIM):: boundBox
  integer :: blkCount
  integer,dimension(MAXBLOCKS) :: blkList
!----------------------------------------------------------------------

  if(.not.useParticles) then
     partPosInitialized = .true.
  end if
  if(partPosInitialized) return

  !CD: We need to move the call to pt_initLocal to this level. 
  !Otherwise, if it is in pt_initPositions, we get a deadlock
  !when the number of blocks are not the same on all processors.
  call pt_initLocal()

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  call Grid_getLocalNumBlks(localNumBlocks)

  partPosInitialized=.true.
  updateRefine = .false.

  numLocalPreviousTypes = 0
  if(.not.updateRefine) then
     pt_numLocal = 0
  else
     pt_numLocal=sum(pt_typeInfo(PART_LOCAL,1:NPART_TYPES))
  end if
  do i = 1,NPART_TYPES
     if(.not.updateRefine) then
        numLocalThisType=0
     else
        numLocalThisType=pt_typeInfo(PART_LOCAL,i)
     end if

     b=0
     numNewLocalThisType = 0
     numPreviousLocal = pt_numLocal
     do while((b<blkCount).and.partPosInitialized)
        b=b+1
        blockID = blkList(b)
!!        print*,'pt_initPositions',blockID, pt_typeInfo(PART_INITMETHOD,i)
        select case(pt_typeInfo(PART_INITMETHOD,i))
        case(LATTICE)
           call pt_initPositionsLattice(blockID,partPosInitialized)
        case(WITH_DENSITY, CELLMASS, REJECTION)
           call pt_initPositionsWithDensity(blockID,partPosInitialized)
        case(CUSTOM)
           call pt_initPositions(blockID,partPosInitialized)
        case default
           call Driver_abortFlash("Particles_initPosition: no valid initialization method")
        end select
        numNewLocalThisType = pt_numLocal - numPreviousLocal
        pt_typeInfo(PART_LOCAL,i) = numNewLocalThisType + numLocalThisType
     enddo
#ifdef TYPE_PART_PROP
     particles(TYPE_PART_PROP, &
          pt_numLocal-numNewLocalThisType+1:pt_numLocal) = pt_typeInfo(PART_TYPE,i) 
#endif
     numLocalThisType=pt_typeInfo(PART_LOCAL,i)
     numLocalPreviousTypes = numLocalPreviousTypes + numLocalThisType
  end do
  pt_numLocal=sum(pt_typeInfo(PART_LOCAL,1:NPART_TYPES))
!!  print*,'Particles_initPositions: pt_numLocal now is',pt_numLocal

  pt_posInitialized = partPosInitialized

#ifdef DEBUG_PARTICLES
  if (pt_meshMe == MASTER_PE .OR. pt_meshNumProcs .LE. 4) then
     print*,'Particles_initPositions on processor', pt_meshMe, 'done, pt_numLocal=',pt_numLocal
  end if
#endif

  ! Now do the active particles

#ifdef MASS_PART_PROP

  call pt_initPositions(-1,partPosInitialized)

#ifdef PASSIVE_PART_TYPE
  pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE) = pt_numLocal - &
       pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE)
  pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE) = 1 + &
       pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE)
#else
  pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE) = pt_numLocal
  pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE) = 1  
#endif

#ifdef TYPE_PART_PROP
  particles(TYPE_PART_PROP, &
       pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE):pt_numLocal) = ACTIVE_PART_TYPE
#endif

  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (pt_meshMe == 0) print *, "Done with particles..."
  
  pt_posInitialized = .true.
  partPosInitialized = .true.
  
  call Particles_updateRefinement(localNumBlocks)

#endif

  call pt_createTag()

  pt_posInitialized = .true.
  partPosInitialized = .true.
  updateRefine = .false.

  return
  
end subroutine Particles_initPositions
