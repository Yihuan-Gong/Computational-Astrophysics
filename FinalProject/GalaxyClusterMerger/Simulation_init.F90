subroutine Simulation_init()

#include "Flash.h"

  use Simulation_data

  use Grid_Data, ONLY : gr_refineOnParticleCount !! HACK

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use IO_interface, ONLY : IO_getScalar

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  ! Temporary variables.

  real, save :: kpc, Msun, keV, kms, Myr, erg_per_keV
  real       :: result1, result2, rmin, c_s
  integer :: ierr, i, numPEs, ii, jj, kk, totzones, halo, hydroLimit

  real, external :: sim_getMagField
  logical :: partPosInitialized, updateRefine

  logical, save :: activeParticles

  character(len=80) :: sim_prol_method_str
  character(len=1) :: dummy_char

  character(len=80) :: r_string = "/fields/radius"
  character(len=80) :: dens_string = "/fields/density"
  character(len=80) :: pden_string = "/fields/dark_matter_density"
  character(len=80) :: pres_string = "/fields/pressure"
  character(len=80) :: metl_string = "/fields/metallicity"
  character(len=80) :: vtan_string = "/fields/tangential_velocity"
  character(len=80) :: gpot_string = "/fields/gravitational_potential"
  character(len=80) :: grav_string = "/fields/gravitational_field"

  real, dimension(:), allocatable :: rbub, theta, phi

#ifdef MASS_PART_PROP
  activeParticles = .true.
#else
  activeParticles = .false.
#endif

!===============================================================================

  ! Initialization

  call mpi_comm_rank(MPI_COMM_WORLD, sim_meshMe, ierr)

  call PhysicalConstants_get("pi", sim_pi)

  kpc = 3.0856775807E21
  Myr = 31.5576e12
  Msun = 1.9889225E33
  k_B = 1.381e-16
  kms = 1.0e5
  N_a = 6.022e23
  mueinv = 0.875
  keV = 11604440.207109345
  m_e = 9.109e-28
  erg_per_keV = 1.602176562e-09

  call RuntimeParameters_get("profile1", sim_profile1)
  call RuntimeParameters_get("profile2", sim_profile2)
  call RuntimeParameters_get("magFile", sim_magFile)
  call RuntimeParameters_get("testSingleCluster", sim_testSingleCluster)
  call RuntimeParameters_get("mainClusterFixed", sim_mainClusterFixed)
  call RuntimeParameters_get("nsubzones", sim_subZones)
  call RuntimeParameters_get("isGas1", sim_isGas1)
!   call RuntimeParameters_get("isGas1", sim_isGas2)  ! should this be call RuntimeParameters_get("isGas2", sim_isGas2)
  call RuntimeParameters_get("isGas2", sim_isGas2)  
#ifdef MASS_PART_PROP
  call RuntimeParameters_get("particleFile1", sim_particleFile1)
  call RuntimeParameters_get("particleFile2", sim_particleFile2)
#endif
  call RuntimeParameters_get("refineOnSpheres", sim_refineOnSpheres)
  call RuntimeParameters_get("xmax", sim_xMax)
  call RuntimeParameters_get("xmin", sim_xMin)
  call RuntimeParameters_get("ymax", sim_yMax)
  call RuntimeParameters_get("ymin", sim_ymin)
  call RuntimeParameters_get("zmax", sim_zMax)
  call RuntimeParameters_get("zmin", sim_zMin)
  call RuntimeParameters_get("xInit1", sim_xInit1)
  call RuntimeParameters_get("yInit1", sim_yInit1)
  call RuntimeParameters_get("xInit2", sim_xInit2)
  call RuntimeParameters_get("yInit2", sim_yInit2)
  call RuntimeParameters_get("smlrho", sim_smlRho)
  call RuntimeParameters_get("smallx", sim_smallX)
  call RuntimeParameters_get("smallt", sim_smallT)
  call RuntimeParameters_get("smalle", sim_smalle)
  call RuntimeParameters_get("smallp", sim_smallP)
  call RuntimeParameters_get("RefinementDensityCutoff", &
       sim_refinementDensityCutoff)
  call RuntimeParameters_get("RefiningRadius1", sim_refiningRadius1)
  call RuntimeParameters_get("RefiningRadius2", sim_refiningRadius2)
  call RuntimeParameters_get("restart", sim_restart)
  call RuntimeParameters_get("lrefine_min", sim_lrefineMin)
  call RuntimeParameters_get("lrefine_max", sim_lrefineMax)
  call RuntimeParameters_get("setupSubclusterWithRestart", sim_setupSubclusterWithRestart)
  call RuntimeParameters_get("Ldamp", sim_Ldamp)
  call RuntimeParameters_get("Rdamp", sim_Rdamp)
  call RuntimeParameters_get("nblockx", sim_nblockx)
  call RuntimeParameters_get("nblocky", sim_nblocky)
  call RuntimeParameters_get("nblockz", sim_nblockz)
  call RuntimeParameters_get("vxInit1", sim_vxInit1)
  call RuntimeParameters_get("vxInit2", sim_vxInit2)
  call RuntimeParameters_get("vyInit1", sim_vyInit1)
  call RuntimeParameters_get("vyInit2", sim_vyInit2)
  call RuntimeParameters_get("readMetals", sim_readMetals)
  call RuntimeParameters_get("readVelocity", sim_readVelocity)
  call RuntimeParameters_get("perturbDensity", sim_perturbDensity)
  call RuntimeParameters_get("densPerturbScale", sim_densPerturbScale)
  call RuntimeParameters_get("bubbleMethod", sim_bubbleMethod)
  call RuntimeParameters_get("bubbleEnergy", sim_bubbleEnergy)
  call RuntimeParameters_get("bubbleStart", sim_bubbleStart)
  call RuntimeParameters_get("bubbleDuration", sim_bubbleDuration)
  call RuntimeParameters_get("bubbleXC", sim_bubbleXC)
  call RuntimeParameters_get("bubbleYC", sim_bubbleYC)
  call RuntimeParameters_get("bubbleZC", sim_bubbleZC)
  call RuntimeParameters_get("numBubbles", sim_numBubbles)
  call RuntimeParameters_get("turbVelocity", sim_turbVelocity)
#ifdef MAGX_VAR
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('forceHydroLimit',sim_forceHydroLimit)
#endif
  call RuntimeParameters_get("insertBubbles", sim_insertBubbles)
  call RuntimeParameters_get("insertBubblesOnRestart", sim_insertBubblesOnRestart)
  call RuntimeParameters_get("bubbleRadius", sim_bubbleRadius)
  call RuntimeParameters_get("bubbleHeight", sim_bubbleHeight)
  call RuntimeParameters_get("bubbleEntropy", sim_bubbleEntropy)
  call RuntimeParameters_get("useBubbleParticles", sim_useBubbleParticles)
  call RuntimeParameters_get("numBubbleParticles", sim_numBubbleParticles)
  call RuntimeParameters_get("bubbleAxis", sim_bubbleAxis)
  call RuntimeParameters_get("useJets", sim_useJets)
  call RuntimeParameters_get("jetAxis", sim_jetAxis)
  call RuntimeParameters_get("jetRadius", sim_jetRadius)
  call RuntimeParameters_get("jetVelocity", sim_jetVelocity)
  call RuntimeParameters_get("jetMach", sim_jetMach)
  call RuntimeParameters_get("jetDensity", sim_jetDensity)
  call RuntimeParameters_get("jetStart", sim_jetStart)
  call RuntimeParameters_get("jetDuration", sim_jetDuration)
  call RuntimeParameters_get("jetHeight", sim_jetHeight)
  call RuntimeParameters_get("jetAngle", sim_jetAngle)
  call RuntimeParameters_get("jetPeriod", sim_jetPeriod)
  call RuntimeParameters_get("rescaleMagEnergy", sim_rescaleMagEnergy)
  call RuntimeParameters_get("refine_on_particle_count", sim_refineOnParticleCount)
  call RuntimeParameters_get("smlrho", sim_smlrho)
  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("rClr1", sim_rClr1)
  call RuntimeParameters_get("rClr2", sim_rClr2)

#ifdef MAGX_VAR
  if (sim_forceHydroLimit .or. sim_restart) then
     hydroLimit = 1
  else
     hydroLimit = 0
  endif
#endif

  call mpi_comm_size(MPI_COMM_WORLD, numPEs, ierr)

  ! Compute derived quantities.

  call PhysicalConstants_get("Newton", sim_Newton)
  sim_poisfact = 4.*PI*sim_Newton

  sim_xInit1 = sim_xInit1*kpc
  sim_yInit1 = sim_yInit1*kpc
  sim_xInit2 = sim_xInit2*kpc
  sim_yInit2 = sim_yInit2*kpc
  sim_rClr1 = sim_rClr1*kpc
  sim_rClr2 = sim_rClr2*kpc
  sim_bubbleXC = sim_bubbleXC*kpc
  sim_bubbleYC = sim_bubbleYC*kpc
  sim_bubbleZC = sim_bubbleZC*kpc
  sim_refiningRadius1 = sim_refiningRadius1*kpc
  sim_refiningRadius2 = sim_refiningRadius2*kpc
  sim_bubbleRadius = sim_bubbleRadius*kpc
  sim_bubbleHeight = sim_bubbleHeight*kpc
  sim_jetRadius = sim_jetRadius*kpc
  sim_jetVelocity = sim_jetVelocity*1.0e5
  sim_jetHeight = sim_jetHeight*kpc
  sim_vxInit1 = sim_vxInit1 * 1.0e5
  sim_vxInit2 = sim_vxInit2 * 1.0e5
  sim_vyInit1 = sim_vyInit1 * 1.0e5
  sim_vyInit2 = sim_vyInit2 * 1.0e5
  sim_bubbleStart = sim_bubbleStart * Myr
  sim_bubbleDuration = sim_bubbleDuration * Myr
  sim_jetStart = sim_jetStart * Myr
  sim_jetDuration = sim_jetDuration * Myr
  sim_jetPeriod = sim_jetPeriod * Myr
  c_s = sim_jetVelocity/sim_jetMach
  sim_jetEnergy = (9./10.)*c_s*c_s
  sim_jetAngle = sim_jetAngle * PI/180.

  if (.not. sim_mainClusterFixed) then
     sim_Ldamp = 1.0e88
     sim_Rdamp = 1.0e88
  else
     sim_Ldamp = sim_Ldamp * kpc
     sim_Rdamp = sim_Rdamp * kpc
  endif

  if (sim_testSingleCluster) then
     sim_vxInit1 = 0.0
     sim_vxInit2 = 0.0
     sim_vyInit1 = 0.0
     sim_vyInit2 = 0.0
     sim_refiningRadius2 = 0.0
  endif

  if (sim_mainClusterFixed) then
     sim_vxInit1 = 0.0
     sim_vyInit1 = 0.0
  endif

  if (sim_meshMe == 0) then

     call read_num_points(sim_profile1, numPoints1)
     print *, "numPoints1 = ", numPoints1

  endif

  call mpi_bcast(numPoints1, 1, FLASH_INTEGER, 0, MPI_COMM_WORLD, ierr)

  allocate(r1(numPoints1))
  allocate(pden1(numPoints1))
  allocate(dens1(numPoints1))
  allocate(pres1(numPoints1))
  allocate(gpot1(numPoints1))
  allocate(grav1(numPoints1))
  allocate(metl1(numPoints1))
  allocate(vtan1(numPoints1))

  if (sim_meshMe == 0) then
     call read_profile(sim_profile1, r_string, r1)  ! Radius (see InitialCond.ipynb)
#ifdef PDEN_VAR
     call read_profile(sim_profile1, pden_string, pden1)
#endif
     call read_profile(sim_profile1, gpot_string, gpot1)
     call read_profile(sim_profile1, grav_string, grav1)

     if (sim_isGas1) then
        call read_profile(sim_profile1, dens_string, dens1)
        call read_profile(sim_profile1, pres_string, pres1)
        if (sim_readMetals) then
            call read_profile(sim_profile1, metl_string, metl1)
        endif
        if (sim_readVelocity) then
            call read_profile(sim_profile1, vtan_string, vtan1)
        endif
     endif

     gpot1(:) = -gpot1(:)
     grav1(:) = -grav1(:)

  endif

  call mpi_bcast(r1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(pden1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dens1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(pres1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(gpot1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(grav1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(metl1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(vtan1, numPoints1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

  sim_rMax1 = r1(numPoints1)
  sim_totMass1 = grav1(numPoints1)*sim_rMax1*sim_rMax1/sim_Newton

  if (.not. sim_testSingleCluster) then

     if (sim_meshMe == 0) then

        call read_num_points(sim_profile2, numPoints2)
        print *, "numPoints2 = ", numPoints2

     endif

     call mpi_bcast(numPoints2, 1, FLASH_INTEGER, 0, MPI_COMM_WORLD, ierr)

     allocate(r2(numPoints2))
     allocate(pden2(numPoints2))
     allocate(gpot2(numPoints2))
     allocate(grav2(numPoints2))
     allocate(dens2(numPoints2))
     allocate(pres2(numPoints2))
     allocate(metl2(numPoints2))

     if (sim_meshMe == 0) then

       call read_profile(sim_profile2, r_string, r2)
#ifdef PDEN_VAR
       call read_profile(sim_profile2, pden_string, pden2)  ! Read DM density for cluster 2
#endif
       if (sim_isGas2) then
          call read_profile(sim_profile2, dens_string, dens2)  ! Read gas density for cluster 2
          call read_profile(sim_profile2, pres_string, pres2)
          if (sim_readMetals) then
              call read_profile(sim_profile2, metl_string, metl2)
          endif
       endif

       call read_profile(sim_profile2, gpot_string, gpot2)
       call read_profile(sim_profile2, grav_string, grav2)

       gpot2(:) = -gpot2(:)
       grav2(:) = -grav2(:)

     endif

     call mpi_bcast(r2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(pden2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(dens2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(pres2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(gpot2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(grav2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(metl2, numPoints2, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

     sim_rMax2 = r2(numPoints2)
     sim_totMass2 = grav2(numPoints2)*sim_rMax2*sim_rMax2/sim_Newton

  endif

  if (.not. sim_restart) then

#ifdef MASS_PART_PROP
     if (sim_meshMe == 0) then

        call read_particle_number(sim_particleFile1, sim_numParticles1)
        print *, "Number of particles for halo 1 = ", sim_numParticles1
        if (.not. sim_testSingleCluster) then
           call read_particle_number(sim_particleFile2, sim_numParticles2)
           print *, "Number of particles for halo 2 = ", sim_numParticles2
        endif
        call read_particle_mass(sim_particleFile1, sim_pmass1)
        if (.not. sim_testSingleCluster) then
           call read_particle_mass(sim_particleFile2, sim_pmass2)
        endif
     endif

     call mpi_bcast(sim_numParticles1, 1, FLASH_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(sim_numParticles2, 1, FLASH_INTEGER, 0, MPI_COMM_WORLD, ierr)

     call mpi_bcast(sim_pmass1, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(sim_pmass2, 1, FLASH_REAL, 0, MPI_COMM_WORLD, ierr)

     sim_pmass = 0.5*(sim_pmass1+sim_pmass2)

     if (sim_meshMe == 0) then
        print *, "Particle mass = ", sim_pmass, sim_pmass1, sim_pmass2
     endif

#endif

     if (sim_refineOnParticleCount) gr_refineOnParticleCount = .false.

     sim_xCtr(1) = sim_xInit1
     sim_yCtr(1) = sim_yInit1
     sim_zCtr(1) = 0.5*(sim_zMin+sim_zMax)

     sim_xCtr(2) = sim_xInit2
     sim_yCtr(2) = sim_yInit2
     sim_zCtr(2) = 0.5*(sim_zMin+sim_zMax)

     nsubinv = 1./real(sim_subZones)
     nsubvolinv = nsubinv**3

     call mpi_barrier(MPI_COMM_WORLD, ierr)

#ifdef MAGX_VAR     

      ! Read the magnetic_file.h5 if sim_forceHydroLimit==False
     if (.not. sim_forceHydroLimit .or. sim_perturbDensity .or. sim_turbVelocity) then

        call open_mag_file(sim_magFile, sim_nx, sim_ny, sim_nz, &
             sim_Adx, sim_Ady, sim_Adz, &
             sim_Axmin, sim_Aymin, sim_Azmin)

        allocate(sim_Axcoord(sim_nx))
        allocate(sim_Aycoord(sim_ny))
        allocate(sim_Azcoord(sim_nz))

        do ii = 1, sim_nx
           sim_Axcoord(ii) = (real(ii)-0.501)*sim_Adx + sim_Axmin
        enddo

        do jj = 1, sim_ny
           sim_Aycoord(jj) = (real(jj)-0.501)*sim_Ady + sim_Aymin
        enddo

        do kk = 1, sim_nz
           sim_Azcoord(kk) = (real(kk)-0.501)*sim_Adz + sim_Azmin
        enddo

     endif

#endif

     if (.not. sim_mainClusterFixed) then

        sim_vxCtr(1) = sim_vxInit1
        sim_vyCtr(1) = sim_vyInit1
        sim_vzCtr(1) = 0.0

        sim_oaxCtr(1) = 0.
        sim_oayCtr(1) = 0.
        sim_oazCtr(1) = 0.

     endif

     sim_vxCtr(2) = sim_vxInit2
     sim_vyCtr(2) = sim_vyInit2
     sim_vzCtr(2) = 0.0

     sim_oaxCtr(2) = 0.
     sim_oayCtr(2) = 0.
     sim_oazCtr(2) = 0.

  else if (.not. activeParticles) then

     if (sim_mainClusterFixed) then

        sim_xCtr(1) = sim_xInit1
        sim_yCtr(1) = sim_yInit1
        sim_zCtr(1) = 0.5*(sim_zMin+sim_zMax)

        sim_vxCtr(1) = 0.0
        sim_vyCtr(1) = 0.0
        sim_vzCtr(1) = 0.0

        sim_oaxCtr(1) = 0.0
        sim_oayCtr(1) = 0.0
        sim_oazCtr(1) = 0.0

     else

        call IO_getScalar("cluster 1 x", sim_xCtr(1))
        call IO_getScalar("cluster 1 y", sim_yCtr(1))
        call IO_getScalar("cluster 1 z", sim_zCtr(1))
        call IO_getScalar("cluster 1 vx", sim_vxCtr(1))
        call IO_getScalar("cluster 1 vy", sim_vyCtr(1))
        call IO_getScalar("cluster 1 vz", sim_vzCtr(1))
        call IO_getScalar("cluster 1 ax", sim_oaxCtr(1))
        call IO_getScalar("cluster 1 ay", sim_oayCtr(1))
        call IO_getScalar("cluster 1 az", sim_oazCtr(1))

     endif

     if (.not. sim_testSingleCluster) then

        if (sim_setupSubclusterWithRestart) then

           sim_xCtr(2) = sim_xInit2
           sim_yCtr(2) = sim_yInit2
           sim_zCtr(2) = 0.5*(sim_zMin+sim_zMax)
           sim_vxCtr(2) = sim_vxInit2
           sim_vyCtr(2) = sim_vyInit2
           sim_vzCtr(2) = 0.0

           sim_oaxCtr(2) = 0.
           sim_oayCtr(2) = 0.
           sim_oazCtr(2) = 0.

        else

           call IO_getScalar("cluster 2 x", sim_xCtr(2))
           call IO_getScalar("cluster 2 y", sim_yCtr(2))
           call IO_getScalar("cluster 2 z", sim_zCtr(2))
           call IO_getScalar("cluster 2 vx", sim_vxCtr(2))
           call IO_getScalar("cluster 2 vy", sim_vyCtr(2))
           call IO_getScalar("cluster 2 vz", sim_vzCtr(2))
           call IO_getScalar("cluster 2 ax", sim_oaxCtr(2))
           call IO_getScalar("cluster 2 ay", sim_oayCtr(2))
           call IO_getScalar("cluster 2 az", sim_oazCtr(2))

        endif

     endif

  endif

  if (sim_insertBubbles) then
     
     sim_bubbleVolume = 4.*PI/3.*sim_bubbleRadius**3
     
     sim_bubbleX(:) = sim_bubbleXC
     sim_bubbleY(:) = sim_bubbleYC
     sim_bubbleZ(:) = sim_bubbleZC
    
     if (sim_bubbleAxis == 1) then
        sim_bubbleX(1) = sim_bubbleX(1) + sim_bubbleHeight
        sim_bubbleX(2) = sim_bubbleX(2) - sim_bubbleHeight
     else if (sim_bubbleAxis == 2) then
        sim_bubbleY(1) = sim_bubbleY(1) + sim_bubbleHeight
        sim_bubbleY(2) = sim_bubbleY(2) - sim_bubbleHeight
     else if (sim_bubbleAxis == 3) then
        sim_bubbleZ(1) = sim_bubbleZ(1) + sim_bubbleHeight
        sim_bubbleZ(2) = sim_bubbleZ(2) - sim_bubbleHeight
     endif
     
     if (sim_bubbleMethod == 2) then
        sim_bubblePower = sim_bubbleEnergy / sim_bubbleDuration
     else
        sim_bubbleEntropy = sim_bubbleEntropy * erg_per_keV
     endif

     if (sim_useBubbleParticles .and. sim_bubbleMethod == 1) then

        allocate(pxbub1(sim_numBubbleParticles))
        allocate(pybub1(sim_numBubbleParticles))
        allocate(pzbub1(sim_numBubbleParticles))
        allocate(pxbub2(sim_numBubbleParticles))
        allocate(pybub2(sim_numBubbleParticles))
        allocate(pzbub2(sim_numBubbleParticles))
        
        allocate(rbub(sim_numBubbleParticles))
        allocate(theta(sim_numBubbleParticles))
        allocate(phi(sim_numBubbleParticles))
        
        call random_number(rbub)
        call random_number(theta)
        call random_number(phi)
        
        rbub = sim_bubbleRadius*rbub**(1./3.)
        theta = acos(1.0-2.0*theta)
        phi = 2.0*PI*phi

        do i = 1, sim_numBubbleParticles
           
           pxbub1(i) = rbub(i)*cos(phi(i))*sin(theta(i)) + sim_bubbleX(1)
           pybub1(i) = rbub(i)*sin(phi(i))*sin(theta(i)) + sim_bubbleY(1)
           pzbub1(i) = rbub(i)*cos(theta(i)) + sim_bubbleZ(1)
           
        enddo
        
        if (sim_numBubbles == 2) then
           
           call random_number(rbub)
           call random_number(theta)
           call random_number(phi)

           rbub = sim_bubbleRadius*rbub**(1./3.)
           theta = acos(1.0-2.0*theta)
           phi = 2.0*PI*phi
           
           do i = 1, sim_numBubbleParticles
              
              pxbub2(i) = rbub(i)*cos(phi(i))*sin(theta(i)) + sim_bubbleX(2)
              pybub2(i) = rbub(i)*sin(phi(i))*sin(theta(i)) + sim_bubbleY(2)
              pzbub2(i) = rbub(i)*cos(theta(i)) + sim_bubbleZ(2)

           enddo

        endif

        deallocate(rbub)
        deallocate(theta)
        deallocate(phi)
        
     endif
     
  endif
  
  sim_dvolMin = (sim_xMax-sim_xMin)*(sim_zMax-sim_zMin)*(sim_zMax-sim_zMin)
  sim_dvolMin = sim_dvolMin / (NXB*NYB*NZB*sim_nblockx*sim_nblocky*sim_nblockz*2**(3*sim_lrefineMax-3))

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  return

end subroutine Simulation_init
