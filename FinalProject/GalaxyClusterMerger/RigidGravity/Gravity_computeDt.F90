!!****if* source/Simulation/SimulationMain/magnetoHD/ClusterSloshingMagneticNew/Gravity_computeDt
!!
!! NAME
!!
!!  Gravity_computeDt
!!  
!! SYNOPSIS
!!
!!  Gravity_computeDt(integer(IN)        :: blockID,
!!                    real (OUT)         :: dt_grav,
!!                    integer(:)(INOUT)  :: dt_minloc(5))
!!
!! DESCRIPTION
!!
!!  Compute the timestep limiter due to the gravitational solver.
!!
!! ARGUMENTS
!!
!!  dt_grav:       Will Return the limiting timestep. Should be
!!                 set to a large value (1.D99) on input.
!!  dt_minloc(5):  An array to receive information about which
!!                 processor, block, and zone was responsible
!!                 for setting the limiting timestep.  The order
!!                 is i, j, k, b, p, where (i,j,k) = zone
!!                 indices, b = local block ID, and p = PE #.
!!                 This routine should only modify these values
!!                 if it changes dt_grav.
!!  blockID:       The local ID of the block to compute the
!!                 limiter on.
!!
!!***

subroutine Gravity_computeDt (blockID, dt_grav, dt_minloc)

!==============================================================================

  use Simulation_data, ONLY : sim_vxCtr, sim_vyCtr, sim_vzCtr, &
       sim_xCtr, sim_yCtr, sim_zCtr, sim_testSingleCluster, &
       sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
       sim_mainClusterFixed, sim_nBlockx, sim_nBlocky, sim_nBlockz, &
       sim_lrefineMax, sim_axCtr, sim_ayCtr, sim_azCtr

  use Grid_interface, ONLY: Grid_getBlkBoundBox, Grid_getDeltas

  use Gravity_data, ONLY : grv_meshMe

  implicit none
  
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer, intent(IN)    ::  blockID
  integer, intent(INOUT) ::  dt_minloc(5)
  real,intent(OUT)       ::  dt_grav

  real, parameter :: small = 1.0e-10
  real :: dtx, dty, dtz, dtnew, min_dx, min_dy, min_dz
  real :: delta(MDIM), lowerBound(MDIM), boundBox(2,MDIM)
  
  integer :: halo, icell, jcell, kcell

  if (.not. sim_testSingleCluster) then

     call Grid_getBlkBoundBox(blockID,boundBox)    ! physical bounding box of the block
     lowerBound = boundBox(1,:)
     call Grid_getDeltas(blockID,delta)

     if (sim_mainClusterFixed) then
        halo = 2
     else
        halo = 1
     endif

     min_dx = (sim_xMax - sim_xMin) / (2**(sim_lrefineMax-1) * NXB * sim_nBlockx)
     min_dy = (sim_yMax - sim_yMin) / (2**(sim_lrefineMax-1) * NYB * sim_nBlocky)
     min_dz = (sim_zMax - sim_zMin) / (2**(sim_lrefineMax-1) * NZB * sim_nBlockz)

     do while (halo <= 2)

        dtx = HUGE(1.0)
        dty = HUGE(1.0)
        dtz = HUGE(1.0)

        if (sim_vxCtr(halo) > small) dtx = min_dx / sim_vxCtr(halo)
        if (sim_vyCtr(halo) > small) dty = min_dy / sim_vyCtr(halo)
        if (sim_vzCtr(halo) > small) dtz = min_dz / sim_vzCtr(halo)

        dtnew = 0.5 * min(dtx, dty, dtz)
        
        if (dtnew < dt_grav) then

           dt_grav = dtnew

           icell = int((sim_xCtr(halo)-lowerBound(1))/delta(1)) + 1
           jcell = int((sim_yCtr(halo)-lowerBound(2))/delta(2)) + 1
           kcell = int((sim_zCtr(halo)-lowerBound(3))/delta(3)) + 1

           if (icell <= NXB .and. jcell <= NYB .and. kcell <= NZB) then
              !! information about where the minimum restriction took place
              dt_minloc(1) = icell + NGUARD 
              dt_minloc(2) = jcell + NGUARD
              dt_minloc(3) = kcell + NGUARD
              dt_minloc(4) = blockID
              dt_minloc(5) = grv_meshMe
           endif

        endif

	halo = halo + 1

     enddo

  else

     ! Don't really need to worry about it

     dt_grav = huge(1.)

  endif

  
  return

end subroutine Gravity_computeDt
