!!****if* source/Simulation/SimulationMain/magnetoHD/ClusterSloshingMagneticNew/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()
  
  use Simulation_data, ONLY : sim_refiningRadius1, sim_refiningRadius2, &
       sim_refinementDensityCutoff, sim_xCtr, sim_yCtr, sim_zCtr, &
       sim_refineOnSpheres, sim_testSingleCluster

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,&
                        gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_min, &
       lrefine_max, lrefine
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_markBlkRefine, Grid_getBlkIndexLimits, Grid_markRefineSpecialized, &
    Grid_markBlkDerefine, Grid_getDeltas
  use Driver_interface, ONLY : Driver_getNStep

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref, ierr
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  integer :: lb, blockCount, blockID
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blockList
  logical,dimension(maskSize) :: gcMask

  real, dimension(:,:,:,:), pointer :: solnData

  real :: specs(4)
  integer :: nstep

  logical :: deref

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  call Driver_getNStep(nstep)

  ! Restrict second-derivative refinement to 
  ! high-density regions
  
  do lb = 1, blockCount
     
     blockID = blockList(lb)
     
     call Grid_getBlkPtr(blockID, solnData)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     if (nstep > 1) then

        if (maxval(solnData(DENS_VAR, &
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) < &
             sim_refinementDensityCutoff) then
           
           refine(blockID) = .false.
        
        endif
     
     endif

     call Grid_releaseBlkPtr(blockID, solnData)
        
  enddo

  do lb = 1, blockCount

     blockID = blockList(lb)

     derefine(blockID) = (.not. refine(blockID)) .and. &
          (lrefine(blockID)>lrefine_min) &
          .and. (.not. stay(blockID))

  enddo

  if (sim_refineOnSpheres) then

    ! Refine a spherical volume around the centers of the clusters
    ! all the way to lrefine_max
    
     specs(1) = sim_xCtr(1)
     specs(2) = sim_yCtr(1)
     specs(3) = sim_zCtr(1)
     specs(4) = sim_refiningRadius1
     
     call Grid_markRefineSpecialized(INRADIUS, 4, specs, lrefine_max)
     
     if (.not. sim_testSingleCluster .and. sim_refiningRadius2 > 0.0) then

        specs(1) = sim_xCtr(2)
        specs(2) = sim_yCtr(2)
        specs(3) = sim_zCtr(2)
        specs(4) = sim_refiningRadius2
        
        call Grid_markRefineSpecialized(INRADIUS, 4, specs, lrefine_max)

     endif

  endif

  return
end subroutine Grid_markRefineDerefine

