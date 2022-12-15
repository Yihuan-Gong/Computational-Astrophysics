subroutine pt_initPositions(blockID, partPosInitialized)
    
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, &
       pt_posInitialized, pt_meshMe, pt_meshNumProcs
  use Simulation_data, ONLY: sim_vxInit1, sim_vyInit1, sim_xCtr, sim_yCtr, &
       sim_zCtr, sim_numParticles1, sim_numParticles2, sim_vxInit2, sim_vyInit2, &
       sim_xMax, sim_xMin, sim_yMax, sim_particleFile1, sim_particleFile2, &
       sim_yMin, sim_zMax, sim_zMin, sim_testSingleCluster, sim_pmass1, sim_pmass2
  use Grid_data, ONLY: gr_refineOnParticleCount

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
  
  integer, INTENT(in) :: blockID
  logical, INTENT(inout) :: partPosInitialized

  ! Local variables.
  
  integer, save :: p

  integer       :: i, j, k, halo, ierr
  integer       :: localNumParticles1, localNumParticles2
  integer       :: particleStart1, particleStart2
  logical       :: isInDomain
  real          :: x,y,z,vx,vy,vz
  real, dimension(:), allocatable :: xpos1, ypos1, zpos1, xvel1, yvel1, zvel1
  real, dimension(:), allocatable :: xpos2, ypos2, zpos2, xvel2, yvel2, zvel2

!==============================================================================

#ifdef MASS_PART_PROP

  if (blockID > -1) return

  ! Particle slot number (incremented and saved between calls)

  p = pt_numLocal
 
  partPosInitialized = .true.

  allocate(xpos1(sim_numParticles1))
  allocate(ypos1(sim_numParticles1))
  allocate(zpos1(sim_numParticles1))
  allocate(xvel1(sim_numParticles1))
  allocate(yvel1(sim_numParticles1))
  allocate(zvel1(sim_numParticles1))

  if (.not. sim_testSingleCluster) then
     allocate(xpos2(sim_numParticles2))
     allocate(ypos2(sim_numParticles2))
     allocate(zpos2(sim_numParticles2))
     allocate(xvel2(sim_numParticles2))
     allocate(yvel2(sim_numParticles2))
     allocate(zvel2(sim_numParticles2))
  endif

  if (pt_meshMe == 0) then

     call read_particles(sim_particleFile1, &
          xpos1, ypos1, zpos1, xvel1, yvel1, zvel1)
     if (.not. sim_testSingleCluster) &
          call read_particles(sim_particleFile2, &
               xpos2, ypos2, zpos2, xvel2, yvel2, zvel2)

  endif

  call mpi_bcast(xpos1, sim_numParticles1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ypos1, sim_numParticles1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(zpos1, sim_numParticles1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(xvel1, sim_numParticles1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(yvel1, sim_numParticles1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(zvel1, sim_numParticles1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

  if (.not. sim_testSingleCluster) then

     call mpi_bcast(xpos2, sim_numParticles2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(ypos2, sim_numParticles2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(zpos2, sim_numParticles2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(xvel2, sim_numParticles2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(yvel2, sim_numParticles2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(zvel2, sim_numParticles2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     
  endif

  if (pt_meshMe == pt_meshNumProcs - 1) then
     localNumParticles1 = sim_numParticles1/pt_meshNumProcs + &
          mod(sim_numParticles1, pt_meshNumProcs)
     if (.not. sim_testSingleCluster) then
        localNumParticles2 = sim_numParticles2/pt_meshNumProcs + &
             mod(sim_numParticles2, pt_meshNumProcs)
     endif
  else
     localNumParticles1 = sim_numParticles1/pt_meshNumProcs
     if (.not. sim_testSingleCluster) then
        localNumParticles2 = sim_numParticles2/pt_meshNumProcs
     endif
  endif
  
  particleStart1 = pt_meshMe*(sim_numParticles1/pt_meshNumProcs) + 1
  if (.not. sim_testSingleCluster) &
       particleStart2 = pt_meshMe*(sim_numParticles2/pt_meshNumProcs) + 1

  do i = particleStart1, localNumParticles1+particleStart1-1
        
     x = xpos1(i) + sim_xCtr(1)
     y = ypos1(i) + sim_yCtr(1)
     z = zpos1(i) + sim_zCtr(1)

     !print *, x, y, z

     isInDomain = ((x > sim_xMin .and. x < sim_xMax) .and. &
          (y > sim_yMin .and. y < sim_yMax) .and. &
          (z > sim_zMin .and. z < sim_zMax))

     if (isInDomain) then

        vx = xvel1(i) + sim_vxInit1
        vy = yvel1(i) + sim_vyInit1
        vz = zvel1(i)
        
        p = p + 1
        
        particles(POSX_PART_PROP,p) = x
        particles(POSY_PART_PROP,p) = y
        particles(POSZ_PART_PROP,p) = z
        particles(VELX_PART_PROP,p) = vx
        particles(VELY_PART_PROP,p) = vy
        particles(VELZ_PART_PROP,p) = vz
        particles(MASS_PART_PROP,p) = sim_pmass1
        particles(BLK_PART_PROP,p)  = real(UNKNOWN)
        particles(PROC_PART_PROP,p)  = real(UNKNOWN)
        particles(HALO_PART_PROP,p)  = real(1.0)
        
     endif

  enddo
  
  if (.not. sim_testSingleCluster) then

       do i = particleStart2, localNumParticles2+particleStart2-1
        
          x = xpos2(i) + sim_xCtr(2)
          y = ypos2(i) + sim_yCtr(2)
          z = zpos2(i) + sim_zCtr(2)
          
          isInDomain = ((x > sim_xMin .and. x < sim_xMax) .and. &
               (y > sim_yMin .and. y < sim_yMax) .and. &
               (z > sim_zMin .and. z < sim_zMax))
          
          if (isInDomain) then
             
             vx = xvel2(i) + sim_vxInit2
             vy = yvel2(i) + sim_vyInit2
             vz = zvel2(i)
             
             p = p + 1
             
             particles(POSX_PART_PROP,p) = x
             particles(POSY_PART_PROP,p) = y
             particles(POSZ_PART_PROP,p) = z
             particles(VELX_PART_PROP,p) = vx
             particles(VELY_PART_PROP,p) = vy
             particles(VELZ_PART_PROP,p) = vz
             particles(MASS_PART_PROP,p) = sim_pmass2
             particles(BLK_PART_PROP,p)  = real(UNKNOWN)
             particles(PROC_PART_PROP,p)  = real(UNKNOWN)
             particles(HALO_PART_PROP,p)  = real(2.0)

          endif
          
       enddo

  endif

!==============================================================================

! Set the particle database local number of particles.

  pt_numLocal = p

  gr_refineOnParticleCount = .true.

  deallocate(xpos1,ypos1,zpos1)
  deallocate(xvel1,yvel1,zvel1)
  if (.not. sim_testSingleCluster) then
     deallocate(xpos2,ypos2,zpos2)
     deallocate(xvel2,yvel2,zvel2)
  endif

  return

#endif

!==============================================================================
  
end subroutine pt_initPositions
