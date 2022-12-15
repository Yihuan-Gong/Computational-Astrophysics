module Simulation_data

  ! Parameters of the problem.

  character(len=80), save :: sim_profile1, sim_profile2
  character(len=80), save :: sim_particleFile1, sim_particleFile2
  character(len=80), save :: sim_magFile
  real, save    :: sim_refiningRadius1, sim_refiningRadius2
  integer, save :: sim_subZones, sim_lrefineMin, sim_lrefineMax
  integer, save :: sim_nblockx, sim_nblocky, sim_nblockz, sim_bubbleMethod
  integer, save :: sim_numParticles1, sim_numParticles2, sim_bubbleAxis
  real, save    :: sim_smalle, sim_smallP, sim_smallT, sim_pi, sim_smlRho
  real, save    :: sim_smallX, sim_bubbleEnergy
  real, save    :: sim_xInit1, sim_yInit1
  real, save  	:: sim_xInit2, sim_yInit2
  real, save    :: sim_refinementDensityCutoff
  real, save    :: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin
  logical, save :: sim_testSingleCluster, sim_isGas1, sim_restart, sim_insertBubbles
  logical, save :: sim_forceHydroLimit, sim_readMetals, sim_refineOnSpheres, sim_isGas2
  logical, save :: sim_mainClusterFixed, sim_perturbDensity
  logical, save :: sim_readVelocity, sim_setupSubclusterWithRestart, sim_turbVelocity
  real, save    :: sim_Ldamp, sim_Rdamp, sim_rClr1, sim_rClr2
  real, save    :: sim_densPerturbScale, sim_bubbleXC, sim_bubbleYC, sim_bubbleZC
  real, save    :: sim_vxInit1, sim_vxInit2, sim_vyInit1, sim_vyInit2, sim_rescaleMagEnergy
  integer, save :: sim_numBubbles
  real, save    :: sim_bubbleRadius, sim_bubbleHeight, sim_bubbleEntropy
  real, save    :: sim_bubbleStart, sim_bubbleDuration, sim_jetAngle, sim_jetPeriod
  real, save    :: sim_jetStart, sim_jetDuration, sim_jetMach, sim_jetHeight
  real, save    :: sim_jetRadius, sim_jetVelocity, sim_jetDensity, sim_jetEnergy
  integer, save :: sim_jetAxis, sim_numBubbleParticles
  logical, save :: sim_useJets, sim_refineOnParticleCount, sim_useBubbleParticles
  logical, save :: sim_insertBubblesOnRestart

  ! Other variables

  real, save    :: nsubinv, nsubvolinv, sim_poisfact, sim_dvolMin
  real, save    :: sim_xCtr(2), sim_yCtr(2), sim_zCtr(2)
  real, save    :: sim_oxCtr(2), sim_oyCtr(2), sim_ozCtr(2)
  real, save    :: sim_Newton, sim_pmass, sim_bubbleVolume, sim_bubblePower
  real, save    :: sim_vxCtr(2), sim_vyCtr(2), sim_vzCtr(2), sim_pmass1
  real, save    :: sim_axCtr(2), sim_ayCtr(2), sim_azCtr(2), sim_pmass2
  real, save    :: sim_oaxCtr(2), sim_oayCtr(2), sim_oazCtr(2)
  integer, save :: sim_nx, sim_ny, sim_nz, sim_nAx, sim_nAy, sim_nAz
  real, save :: sim_Axmin, sim_Aymin, sim_Azmin, sim_Adx, sim_Ady, sim_Adz
  real, save :: sim_Axmax, sim_Aymax, sim_Azmax
  real, allocatable, dimension(:), save :: sim_Axcoord, sim_Aycoord, sim_Azcoord
  real, save :: N_a, k_B, m_e, mueinv
  real, save :: sim_bubbleX(2), sim_bubbleY(2), sim_bubbleZ(2)

  logical, save :: sim_killdivb

  real, save, dimension(:,:,:), allocatable :: sim_Ax
  real, save, dimension(:,:,:), allocatable :: sim_Ay
  real, save, dimension(:,:,:), allocatable :: sim_Az
  real, save, dimension(:,:,:), allocatable :: sim_Dp

  integer, save :: sim_meshMe, numPoints1, numPoints2

  real, allocatable, save, dimension(:) :: r1, pden1, dens1, pres1, gpot1
  real, allocatable, save, dimension(:) :: grav1, metl1, vtan1
  real, allocatable, save, dimension(:) :: r2, pden2, dens2, pres2, gpot2
  real, allocatable, save, dimension(:) :: grav2, metl2
  real, allocatable, save, dimension(:) :: pxbub1, pybub1, pzbub1
  real, allocatable, save, dimension(:) :: pxbub2, pybub2, pzbub2

  real, save :: sim_totMass1, sim_totMass2, sim_rMax1, sim_rMax2

end module Simulation_data
