!!****if* source/Simulation/SimulationMain/GalaxyCluster/Grid_markRefineDerefine
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
  
  use Simulation_data, ONLY : sim_refinementDensityCutoff, sim_pmass

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,&
                        gr_refineOnParticleCount, &
                        gr_minParticlesPerBlk, gr_maxParticlesPerBlk

  use tree, ONLY : newchild, refine, derefine, stay, lrefine_min, &
       lrefine_max, lrefine

  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_markBlkRefine, Grid_getBlkIndexLimits, Grid_markRefineSpecialized, &
    Grid_markBlkDerefine, Grid_getBlkBoundBox, Grid_getDeltas

  use Logfile_interface, ONLY: Logfile_stamp

  use Particles_interface, ONLY: Particles_updateGridVar

  use Driver_interface, ONLY: Driver_abortFlash

  ! Stuff we normally shouldn't use but I'm being sloppy

  use Particles_data, ONLY : pt_posInitialized
  use Driver_data, ONLY : dr_restart

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer   :: l,iref, ierr
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  integer :: lb, blockCount, blockID, oneBlkCount
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blockList
  logical,dimension(maskSize) :: gcMask

  real,dimension(MDIM) :: delta
  real,dimension(LOW:HIGH,MDIM) :: bndBox

  real, dimension(:,:,:,:), pointer :: solnData

  real :: xx, yy, zz
  integer :: i, j, k, imax, jmax, kmax, imin, jmin, kmin
  real :: dvol, delm

  real :: specs(4), lsum(4), asum(4), Mtot, Xcm, Ycm, Zcm

  integer :: nsum = 4

  logical :: deref

  character(len=124) :: str_buffer

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
  
  ! Turn off second-derivative refinement in low gas-density regions

#ifdef PRES_VAR
  do lb = 1, blockCount

     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID, solnData)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     if (maxval(solnData(DENS_VAR, &
          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))) < &
          sim_refinementDensityCutoff) then
        
        refine(blockID) = .false.

     endif

     call Grid_releaseBlkPtr(blockID, solnData)

  enddo
#endif

  ! Refine on particle count, or what would be the particle count

  if (pt_posInitialized .or. dr_restart) then

     call gr_ptMarkRefineDerefine()

  else

     do lb = 1, blockCount

        blockID = blockList(lb)

        call Grid_getBlkPtr(blockID, solnData)

        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

        call Grid_getDeltas(blockID, delta)

        dvol = delta(1)*delta(2)*delta(3)

        oneBlkCount = int(sum(solnData(PDEN_VAR, &
             blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))*dvol)/sim_pmass)

        refine(blockID) = (oneBlkCount > gr_maxParticlesPerBlk) &
             .or. refine(blockID)
        refine(blockID) = refine(blockID) .and. &
             (lrefine(blockID) < lrefine_max)
        deref = (oneBlkCount < gr_minParticlesPerBlk)
        derefine(blockID) = (.not. refine(blockID)) .and. &
          (lrefine(blockID)>lrefine_min) &
          .and. deref .and. (.not. stay(blockID))

        call Grid_releaseBlkPtr(blockID, solnData)

     enddo

  endif

  return
end subroutine Grid_markRefineDerefine

