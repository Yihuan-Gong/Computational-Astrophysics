subroutine Simulation_initBlock(blockId)

  use Simulation_data

  use Eos_interface, ONLY : Eos

  use ut_interpolationInterface, ONLY : ut_hunt

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getDeltas, &
      Grid_putPointData, Grid_getCellCoords, Grid_putBlkData, &
      Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkRefineLevel, &
      Grid_getBlkBoundBox
 
  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
  
  integer, INTENT(in) :: blockId

  ! Temporary variables.
   
  integer :: i, j, k, n, ii, jj, kk
  integer :: ir1, ir2, imin, imax, jmin, jmax, kmin, kmax
  integer :: blockRefineLevel

  real :: temp_zone, dens_zone, pres_zone, velx_zone
  real :: vely_zone, ener_zone, clr1_zone, clr2_zone
  real :: gamc_zone, metl_zone, pden_zone, vtan_zone
  real :: bubb_zone

  real :: sum_dens, sum_pres, sum_mass1, sum_mass2, sum_metl, sum_vtan
  real :: pres_sample, dens_sample, Bmag, metl_sample, sum_pden
  real :: metl_sample1, metl_sample2, pden_sample1, pden_sample2
  real :: pres_sample1, dens_sample1, pres_center, dx_min, pden_sample
  real :: pres_sample2, dens_sample2, Btheta, Bphi, theta
  real :: vtan_sample, vtan_sample1

  real, save :: five_thirds = 5./3.
  real, save :: four_thirds = 4./3.
  real, save :: two_thirds = 2./3.

  real, parameter :: Mpc = 3.0856e24
  real, parameter :: mue_twofifths = 1.052463228428472
  real, parameter :: mh = 1.6737352238051868e-24

  real :: dr1, dr2, rho_aC, B_aC, sample_fact
  real :: bub_rad1, bub_rad2
  
  real, dimension(EOS_NUM) :: eosData
  logical,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: vecLen

  real, external :: interpolate

  ! Block coordinate info etc.

  real, dimension(MDIM)             :: del
  real, dimension(2, MDIM)          :: boundBox
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(LOW:HIGH,MDIM) :: eosRange
  integer :: pos(3), size(3), ipos(3), sample_res
  integer :: nAx, nAy, nAz, ibegin, iend, jbegin, jend, kbegin, kend

  real, dimension(:), allocatable :: xl, yl, zl, xc, yc, zc, xr, yr, zr
  real :: bxmin, bymin, bzmin, bxmax, bymax, bzmax

  real :: xx, yy, zz, dxx, dyy, dzz, dvol, rr, xx2, yy2, zz2, dx, dy, dz, rr1, rr2
  real :: x1, x2, x3, r_cyl, phi
  real, dimension(2) :: xh, yh, zh, dh, db
 
  real, dimension(:,:,:), allocatable :: rhog, t, vx, vy, vz, metl, rhod, &
       p, game, gamc, e, ei, clr1, clr2, bp, bx, by, bz, divb, Ax, Ay, Az, &
       Dp, bubb

#if NSPECIES > 0
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn
#endif

  logical, save :: gcell = .true.

  real, external :: vecPot, scalarPot
  real, pointer, dimension(:,:,:,:) :: facexData, faceyData, facezData

!===============================================================================

  pos(:) = 1

  eosMask = .false.

  call Grid_getBlkRefineLevel(blockID,blockRefineLevel)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockID,del)
  
  size(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  size(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  size(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

#ifdef PDEN_VAR
  allocate(rhod(size(1), size(2), size(3)))
#endif

  if (sim_isGas1 .or. sim_isGas2) then

      allocate(clr1(size(1), size(2), size(3)))
      allocate(clr2(size(1), size(2), size(3)))
      allocate(bubb(size(1), size(2), size(3)))
      allocate(rhog(size(1), size(2), size(3)))
      allocate(t(size(1), size(2), size(3)))
      allocate(vx(size(1), size(2), size(3)))
      allocate(vy(size(1), size(2), size(3)))
      allocate(vz(size(1), size(2), size(3)))
      allocate(e(size(1), size(2), size(3)))
      allocate(ei(size(1), size(2), size(3)))
      allocate(p(size(1), size(2), size(3)))
      allocate(game(size(1), size(2), size(3)))
      allocate(gamc(size(1), size(2), size(3)))
      allocate(bx(size(1), size(2), size(3)))
      allocate(by(size(1), size(2), size(3)))
      allocate(bz(size(1), size(2), size(3)))
      allocate(bp(size(1), size(2), size(3)))
      allocate(divb(size(1), size(2), size(3)))
      allocate(Ax(size(1)+1, size(2)+1, size(3)+1))
      allocate(Ay(size(1)+1, size(2)+1, size(3)+1))
      allocate(Az(size(1)+1, size(2)+1, size(3)+1))
      allocate(Dp(size(1), size(2), size(3)))
      allocate(metl(size(1), size(2), size(3)))

  endif
  
  allocate(xl(size(1)), yl(size(2)), zl(size(3)))
  allocate(xc(size(1)), yc(size(2)), zc(size(3)))
  allocate(xr(size(1)), yr(size(2)), zr(size(3)))

  call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, gcell, zl, size(3))
  call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, gcell, yl, size(2))
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xl, size(1))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, zc, size(3))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, yc, size(2))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xc, size(1))
  call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE, gcell, zr, size(3))
  call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, gcell, yr, size(2))
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xr, size(1))

  dx = del(1)
  dy = del(2)
  dz = del(3)
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                    !
  !   Expand initial conditions (e.g. density, temperature)            !
  !   written as f(r) into f(x,y,z) by the initial core position       !
  !   of cluster1 (sim_Xctr(1), ...) and cluster2 (sim_Xctr(2), ...)   !
  !                                                                    !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)

     dzz = del(3) * nsubinv

     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)

        dyy = del(2) * nsubinv

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

           dxx = del(1) * nsubinv

           sum_dens = 0.
           sum_pres = 0.
           sum_mass1 = 0.
           sum_mass2 = 0.
           sum_metl = 0.
           sum_pden = 0.
           sum_vtan = 0.

           dvol = dxx*dyy*dzz

           do kk = 1, sim_subZones
              
              zz = zl(k) + (real(kk) - 0.5)*dzz
              zh(:) = zz - sim_zCtr(:)
           
              do jj = 1, sim_subZones
                 
                 yy = yl(j) + (real(jj) - 0.5)*dyy
                 yh(:) = yy - sim_yCtr(:)
           
                 do ii = 1, sim_subZones

                    xx = xl(i) + (real(ii) - 0.5)*dxx
                    xh(:) = xx - sim_xCtr(:)
          
                    dh(:) = sqrt(xh(:)**2 + yh(:)**2 + zh(:)**2)  ! The distance from cluster center
           
                    pres_sample1 = 0.
                    pres_sample2 = 0.
                    dens_sample1 = 0.
                    dens_sample2 = 0.
                    metl_sample1 = 0.
                    metl_sample2 = 0.
                    pden_sample1 = 0.
                    pden_sample2 = 0.
		                vtan_sample1 = 0.

                    rr = dh(1)  ! The distance from cluster center

#ifdef PDEN_VAR
                    ! Interplotate pden1 (dark matter density of 1) with r1 (radius list) 
                    ! and rr (distance from the cluster core at current position(x,y,z))
                    pden_sample1 = interpolate(pden1, numPoints1, r1, rr)
#endif

                    if (sim_isGas1) then

                        dens_sample1 = interpolate(dens1, numPoints1, r1, rr)
                        pres_sample1 = interpolate(pres1, numPoints1, r1, rr)
                        if (sim_readMetals) &
                            metl_sample1 = interpolate(metl1, numPoints1, r1, rr)

                        if (sim_readVelocity) then
                          r_cyl = sqrt(xh(1)**2 + yh(1)**2)
                          vtan_sample1 = interpolate(vtan1, numPoints1, r1, r_cyl)
                        endif
                    
                    endif

                    if (.not. sim_testSingleCluster) then

                       rr = dh(2)
#ifdef PDEN_VAR
                       pden_sample2 = interpolate(pden2, numPoints2, r2, rr)
#endif

                       if (sim_isGas2) then 
                          
                          dens_sample2 = interpolate(dens2, numPoints2, r2, rr)
                          pres_sample2 = interpolate(pres2, numPoints2, r2, rr)
                          metl_sample2 = interpolate(metl2, numPoints2, r2, rr)
                          
                       endif
                                              
                    endif

                    if (sim_isGas1 .or. sim_isGas2) then

                        dens_sample = dens_sample1 + dens_sample2
                        pres_sample = pres_sample1 + pres_sample2

                        sum_dens = sum_dens + dens_sample
                        sum_pres = sum_pres + pres_sample

                        if (sim_readMetals) then
                          metl_sample = metl_sample1*dens_sample1 + &
                              metl_sample2*dens_sample2
                          sum_metl = sum_metl + metl_sample
                        endif

                        if (sim_readVelocity) then
                          vtan_sample = vtan_sample1*dens_sample1
                          sum_vtan = sum_vtan + vtan_sample
                        endif

                        if (.not. sim_testSingleCluster) then
                          sum_mass1 = sum_mass1 + dens_sample1
                          sum_mass2 = sum_mass2 + dens_sample2
                        endif

                    endif

                    pden_sample = pden_sample1 + pden_sample2
                    sum_pden = sum_pden + pden_sample

                 enddo
              enddo
           enddo

           pden_zone = sum_pden*nsubvolinv

           if (sim_isGas1 .or. sim_isGas2) then

              dens_zone = sum_dens*nsubvolinv
              pres_zone = sum_pres*nsubvolinv
              bubb_zone = 0.0

              if (sim_insertBubbles .and. sim_bubbleMethod == 1) then
                  
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
                  
                  if (bub_rad1 <= sim_bubbleRadius .or. bub_rad2 <= sim_bubbleRadius) then
                    dens_zone = (pres_zone / sim_bubbleEntropy) ** 0.6 * mh * mue_twofifths
                    bubb_zone = 1.0
                  endif

              endif
            
              if (sim_readMetals) &
                    metl_zone = sum_metl/sum_dens
              if (sim_readVelocity) &
                    vtan_zone = sum_vtan/sum_dens

              clr1_zone = 0.0
              clr2_zone = 0.0

              if (sim_testSingleCluster .or. .not. sim_isGas2) then
                  clr1_zone = 1.0
              else
                  rr1 = sqrt((xc(i)-sim_xCtr(1))**2 + &
                      (yc(j)-sim_yCtr(1))**2 + &
                      (zc(k)-sim_zCtr(1))**2)
                  rr2 = sqrt((xc(i)-sim_xCtr(2))**2 + &
                      (yc(j)-sim_yCtr(2))**2 + &
                      (zc(k)-sim_zCtr(2))**2)
                  if (rr1 <= sim_rClr1) clr1_zone = 1.0
                  if (rr2 <= sim_rClr2) clr2_zone = 1.0
              endif

              if (sim_testSingleCluster) then
                  velx_zone = 0.0
                  vely_zone = 0.0
              else
                  if (sim_isGas2) then
                    velx_zone = (sum_mass1*sim_vxInit1+sum_mass2*sim_vxInit2) / &
                        (sum_mass1+sum_mass2)
                    vely_zone = (sum_mass1*sim_vyInit1+sum_mass2*sim_vyInit2) / &
                        (sum_mass1+sum_mass2)
                  else
                    velx_zone = 0.0
                    vely_zone = 0.0
                  endif
              endif

              if (sim_readVelocity) then
                  phi = atan2(yc(j)-sim_yCtr(1), xc(i)-sim_xCtr(1))
              if (phi < 0) phi = phi + 2.*PI
                  velx_zone = velx_zone - vtan_zone*sin(phi)
                  vely_zone = vely_zone + vtan_zone*cos(phi)
              endif

              gamc_zone = five_thirds

              vx(i,j,k)   = velx_zone
              vy(i,j,k)   = vely_zone
              vz(i,j,k)   = 0.
              rhog(i,j,k) = dens_zone
              clr1(i,j,k) = clr1_zone
              clr2(i,j,k) = clr2_zone
              p(i,j,k)    = pres_zone
              game(i,j,k) = gamc_zone
              gamc(i,j,k) = gamc_zone
              metl(i,j,k) = metl_zone
              bubb(i,j,k) = bubb_zone

              eosData(EOS_PRES) = pres_zone
              eosData(EOS_DENS) = dens_zone

              vecLen = 1

              call Eos(MODE_DENS_PRES, vecLen, eosData)
              
              temp_zone   = eosData(EOS_TEMP)
              ener_zone   = eosData(EOS_EINT)

              t(i,j,k)    = temp_zone
              ei(i,j,k)   = ener_zone
          
          endif

#ifdef PDEN_VAR
          rhod(i,j,k) = pden_zone
#endif

        enddo
     enddo
  enddo

  if (sim_isGas1 .or. sim_isGas2) then

      if (sim_perturbDensity) then

        call Grid_getBlkBoundBox(blockID, boundBox)

        bxmin = boundBox(LOW,IAXIS)-NGUARD*del(IAXIS)
        bymin = boundBox(LOW,JAXIS)-NGUARD*del(JAXIS)
        bzmin = boundBox(LOW,KAXIS)-NGUARD*del(KAXIS)
        
        bxmax = boundBox(HIGH,IAXIS)+NGUARD*del(IAXIS)
        bymax = boundBox(HIGH,JAXIS)+NGUARD*del(JAXIS)
        bzmax = boundBox(HIGH,KAXIS)+NGUARD*del(KAXIS)
        
        ibegin = max(int((bxmin-sim_Axmin)/sim_Adx), 1)
        jbegin = max(int((bymin-sim_Aymin)/sim_Ady), 1)
        kbegin = max(int((bzmin-sim_Azmin)/sim_Adz), 1)

        iend = min(int((bxmax-sim_Axmin)/sim_Adx) + 2, sim_nx)
        jend = min(int((bymax-sim_Aymin)/sim_Ady) + 2, sim_ny)
        kend = min(int((bzmax-sim_Azmin)/sim_Adz) + 2, sim_nz)

        sim_nAx = iend-ibegin+1
        sim_nAy = jend-jbegin+1
        sim_nAz = kend-kbegin+1

        Dp = 0.0

        allocate(sim_Dp(sim_nAz,sim_nAy,sim_nAx))

        call read_density_field(ibegin, jbegin, kbegin, &
              iend, jend, kend, sim_Dp)

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

                  Dp(i,j,k) = scalarPot(ibegin,jbegin,kbegin, &
                      xc(i),yc(j),zc(k))

              enddo
            enddo
        enddo

        deallocate(sim_Dp)

        rhog = rhog*(1.+sim_densPerturbScale*Dp)

      endif

      if (sim_turbVelocity) then
        
        call Grid_getBlkBoundBox(blockID, boundBox)

        bxmin = boundBox(LOW,IAXIS)-NGUARD*del(IAXIS)
        bymin = boundBox(LOW,JAXIS)-NGUARD*del(JAXIS)
        bzmin = boundBox(LOW,KAXIS)-NGUARD*del(KAXIS)
        
        bxmax = boundBox(HIGH,IAXIS)+NGUARD*del(IAXIS)
        bymax = boundBox(HIGH,JAXIS)+NGUARD*del(JAXIS)
        bzmax = boundBox(HIGH,KAXIS)+NGUARD*del(KAXIS)
        
        ibegin = max(int((bxmin-sim_Axmin)/sim_Adx), 1)
        jbegin = max(int((bymin-sim_Aymin)/sim_Ady), 1)
        kbegin = max(int((bzmin-sim_Azmin)/sim_Adz), 1)

        iend = min(int((bxmax-sim_Axmin)/sim_Adx) + 2, sim_nx)
        jend = min(int((bymax-sim_Aymin)/sim_Ady) + 2, sim_ny)
        kend = min(int((bzmax-sim_Azmin)/sim_Adz) + 2, sim_nz)

        sim_nAx = iend-ibegin+1
        sim_nAy = jend-jbegin+1
        sim_nAz = kend-kbegin+1
        
        allocate(sim_Ax(sim_nAz,sim_nAy,sim_nAx))
        allocate(sim_Ay(sim_nAz,sim_nAy,sim_nAx))
        allocate(sim_Az(sim_nAz,sim_nAy,sim_nAx))

        call read_velocity_field(ibegin, jbegin, kbegin, &
              iend, jend, kend, sim_Ax, sim_Ay, sim_Az)
        
        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                  
                  vx(i,j,k) = vx(i,j,k) + vecPot(1,ibegin,jbegin,kbegin,xc(i),yc(j),zc(k)) 
                  vy(i,j,k) = vy(i,j,k) + vecPot(2,ibegin,jbegin,kbegin,xc(i),yc(j),zc(k)) 
                  vz(i,j,k) = vz(i,j,k) + vecPot(3,ibegin,jbegin,kbegin,xc(i),yc(j),zc(k)) 
            
              enddo
            enddo
        enddo
        
        deallocate(sim_Ax)
        deallocate(sim_Ay)
        deallocate(sim_Az)
        
      endif

#ifdef MAGX_VAR
      
      if (sim_killdivb .and. .not. sim_forceHydroLimit .and. &
          blockRefineLevel >= sim_lrefineMin) then
        
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        sample_res = 2**(sim_lrefineMax-blockRefineLevel)
        sample_fact = 1.0/real(sample_res)

        Ax = 0.0
        Ay = 0.0
        Az = 0.0

        call Grid_getBlkBoundBox(blockID, boundBox)
        
        bxmin = boundBox(LOW,IAXIS)-NGUARD*del(IAXIS)
        bymin = boundBox(LOW,JAXIS)-NGUARD*del(JAXIS)
        bzmin = boundBox(LOW,KAXIS)-NGUARD*del(KAXIS)

        bxmax = boundBox(HIGH,IAXIS)+NGUARD*del(IAXIS)
        bymax = boundBox(HIGH,JAXIS)+NGUARD*del(JAXIS)
        bzmax = boundBox(HIGH,KAXIS)+NGUARD*del(KAXIS)

        ibegin = max(int((bxmin-sim_Axmin)/sim_Adx), 1)
        jbegin = max(int((bymin-sim_Aymin)/sim_Ady), 1)
        kbegin = max(int((bzmin-sim_Azmin)/sim_Adz), 1)
        
        iend = min(int((bxmax-sim_Axmin)/sim_Adx) + 2, sim_nx)
        jend = min(int((bymax-sim_Aymin)/sim_Ady) + 2, sim_ny)
        kend = min(int((bzmax-sim_Azmin)/sim_Adz) + 2, sim_nz)
        
        sim_nAx = iend-ibegin+1
        sim_nAy = jend-jbegin+1
        sim_nAz = kend-kbegin+1

        allocate(sim_Ax(sim_nAz,sim_nAy,sim_nAx))
        allocate(sim_Ay(sim_nAz,sim_nAy,sim_nAx))
        allocate(sim_Az(sim_nAz,sim_nAy,sim_nAx))

        call read_magnetic_field(ibegin, jbegin, kbegin, &
              iend, jend, kend, sim_Ax, sim_Ay, sim_Az)

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
                  
                  if (i <= blkLimitsGC(HIGH,IAXIS)) then
                    x1 = xl(i)
                  else
                    x1 = xr(i-1) 
                  endif
                  
                  if (j <= blkLimitsGC(HIGH,JAXIS)) then
                    x2 = yl(j) 
                  else
                    x2 = yr(j-1)
                  endif
                  
                  if (k <= blkLimitsGC(HIGH,KAXIS)) then
                    x3 = zl(k)
                  else
                    x3 = zr(k-1)
                  endif

                  do ii = 1, sample_res
                    dxx = (real(ii)-0.5)*dx*sample_fact
                    Ax(i,j,k) = Ax(i,j,k) + vecPot(1,ibegin,jbegin,kbegin, &
                          x1+dxx,x2,x3) 
                  enddo

                  do jj = 1, sample_res
                    dyy = (real(jj)-0.5)*dy*sample_fact
                    Ay(i,j,k) = Ay(i,j,k) + vecPot(2,ibegin,jbegin,kbegin, &
                          x1,x2+dyy,x3)
                  enddo

                  do kk = 1, sample_res
                    dzz = (real(kk)-0.5)*dz*sample_fact
                    Az(i,j,k) = Az(i,j,k) + vecPot(3,ibegin,jbegin,kbegin, &
                          x1,x2,x3+dzz)
                  enddo

                  Ax(i,j,k) = Ax(i,j,k)*sample_fact
                  Ay(i,j,k) = Ay(i,j,k)*sample_fact
                  Az(i,j,k) = Az(i,j,k)*sample_fact

              enddo
            enddo
        enddo

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
                  facexData(MAG_FACE_VAR,i,j,k) = (Az(i,j+1,k)-Az(i,j,k))/dy - &
                      (Ay(i,j,k+1)-Ay(i,j,k))/dz
              enddo
            enddo
        enddo

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                  faceyData(MAG_FACE_VAR,i,j,k) = (Ax(i,j,k+1)-Ax(i,j,k))/dz - &
                      (Az(i+1,j,k)-Az(i,j,k))/dx
              enddo
            enddo
        enddo     

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                  facezData(MAG_FACE_VAR,i,j,k) = (Ay(i+1,j,k)-Ay(i,j,k))/dx - &
                      (Ax(i,j+1,k)-Ax(i,j,k))/dy
              enddo
            enddo
        enddo
        
        facexData(MAG_FACE_VAR,:,:,:) = facexData(MAG_FACE_VAR,:,:,:) * sim_rescaleMagEnergy / sqrt(4.*PI)
        faceyData(MAG_FACE_VAR,:,:,:) = faceyData(MAG_FACE_VAR,:,:,:) * sim_rescaleMagEnergy / sqrt(4.*PI)
        facezData(MAG_FACE_VAR,:,:,:) = facezData(MAG_FACE_VAR,:,:,:) * sim_rescaleMagEnergy / sqrt(4.*PI)

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)-1
            
            do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)-1
              
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)-1

                  bx(i,j,k)   = 0.5*(facexData(MAG_FACE_VAR,i+1,j,k) + &
                      facexData(MAG_FACE_VAR,i,j,k))
                  by(i,j,k)   = 0.5*(faceyData(MAG_FACE_VAR,i,j+1,k) + &
                      faceyData(MAG_FACE_VAR,i,j,k))
                  bz(i,j,k)   = 0.5*(facezData(MAG_FACE_VAR,i,j,k+1) + &
                      facezData(MAG_FACE_VAR,i,j,k))
                  
                  divb(i,j,k) = (facexData(MAG_FACE_VAR,i+1,j,k) - &
                      facexData(MAG_FACE_VAR,i,j,k))/del(1) &
                      + (faceyData(MAG_FACE_VAR,i,j+1,k) - &
                      faceyData(MAG_FACE_VAR,i,j,k))/del(2) &
                      + (facezData(MAG_FACE_VAR,i,j,k+1) - &
                      facezData(MAG_FACE_VAR,i,j,k))/del(3)
                  
                  bp(i,j,k)   = 0.5*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)

              enddo
            enddo
        enddo

        deallocate(sim_Ax)
        deallocate(sim_Ay)
        deallocate(sim_Az)

        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

      endif

#endif

      e = ei + 0.5*(vx*vx+vy*vy+vz*vz)

      call Grid_putBlkData(blockID, CENTER, DENS_VAR, EXTERIOR, pos,rhog, size)
#ifdef PRES_VAR
      call Grid_putBlkData(blockID, CENTER, PRES_VAR, EXTERIOR, pos,p,    size)
      call Grid_putBlkData(blockID, CENTER, GAME_VAR, EXTERIOR, pos,game, size)
      call Grid_putBlkData(blockID, CENTER, GAMC_VAR, EXTERIOR, pos,gamc, size)
      call Grid_putBlkData(blockID, CENTER, CLR1_MSCALAR, EXTERIOR, pos, clr1, size)
      call Grid_putBlkData(blockID, CENTER, CLR2_MSCALAR, EXTERIOR, pos, clr2, size)
      call Grid_putBlkData(blockID, CENTER, METL_MSCALAR, EXTERIOR, pos, metl, size)
      call Grid_putBlkData(blockID, CENTER, BUBB_MSCALAR, EXTERIOR, pos, bubb, size)
      call Grid_putBlkData(blockID, CENTER, TEMP_VAR,     EXTERIOR, pos,t,    size)
      call Grid_putBlkData(blockID, CENTER, VELX_VAR,     EXTERIOR, pos,vx,   size)
      call Grid_putBlkData(blockID, CENTER, VELY_VAR,     EXTERIOR, pos,vy,   size)
      call Grid_putBlkData(blockID, CENTER, VELZ_VAR,     EXTERIOR, pos,vz,   size)
      call Grid_putBlkData(blockID, CENTER, ENER_VAR,     EXTERIOR, pos,e,    size)
      call Grid_putBlkData(blockID, CENTER, EINT_VAR,     EXTERIOR, pos,ei,   size)
#endif

#ifdef MAGX_VAR
      if (sim_killdivb .and. .not. sim_forceHydroLimit) then
        call Grid_putBlkData(blockID, CENTER, MAGP_VAR, EXTERIOR, pos,bp,   size)
        call Grid_putBlkData(blockID, CENTER, MAGX_VAR, EXTERIOR, pos,bx,   size)
        call Grid_putBlkData(blockID, CENTER, MAGY_VAR, EXTERIOR, pos,by,   size)
        call Grid_putBlkData(blockID, CENTER, MAGZ_VAR, EXTERIOR, pos,bz,   size)
        call Grid_putBlkData(blockID, CENTER, DIVB_VAR, EXTERIOR, pos,divb,   size)
      endif
#endif

      deallocate(clr1)
      deallocate(clr2)
      deallocate(rhog)
      deallocate(bubb)
      deallocate(metl)
      deallocate(t)
      deallocate(vx)
      deallocate(vy)
      deallocate(vz)
      deallocate(e)
      deallocate(ei)
      deallocate(p)
      deallocate(gamc)
      deallocate(game)
      deallocate(bp)
      deallocate(bx)
      deallocate(by)
      deallocate(bz)
      deallocate(divb)
      deallocate(Ax)
      deallocate(Ay)
      deallocate(Az)
      deallocate(Dp)

  endif

#ifdef PDEN_VAR
  call Grid_putBlkData(blockID, CENTER, PDEN_VAR,     EXTERIOR, pos,rhod,   size)
  deallocate(rhod)
#endif

  deallocate(xl)
  deallocate(yl)
  deallocate(zl)
  deallocate(xc)
  deallocate(yc)
  deallocate(zc)
  deallocate(xr)
  deallocate(yr)
  deallocate(zr)
  
  return

end subroutine Simulation_initBlock

function scalarPot(ibegin, jbegin, kbegin, xx, yy, zz)

  use Simulation_data

  implicit none

  integer, intent(IN) :: ibegin, jbegin, kbegin
  real, intent(IN) :: xx, yy, zz

  real, external :: TSC_weight

  real :: scalarPot

  real :: dx, dy, dz

  integer :: i, j, k, ii, jj, kk, ib, jb, kb

  scalarPot = 0.0

  do i = -1, 1
     dx = (xx - sim_Axcoord(ii+i))/sim_Adx
     do j = -1, 1
        dy = (yy - sim_Aycoord(jj+j))/sim_Ady
        do k = -1, 1
           dz = (zz - sim_Azcoord(kk+k))/sim_Adz
           if (ib > 1 .and. ib < sim_nAx .and. &
               jb > 1 .and. jb < sim_nAy .and. &
               kb > 1 .and. kb < sim_nAz) then
              scalarPot = scalarPot + sim_Dp(kb+k,jb+j,ib+i) * &
                   TSC_weight(dx)*TSC_weight(dy)*TSC_weight(dz)
           endif
        enddo
     enddo
  enddo

  return

end function scalarPot

function vecPot(axis, ibegin, jbegin, kbegin, xx, yy, zz)
 
  use Simulation_data

  implicit none

  integer, intent(IN) :: axis, ibegin, jbegin, kbegin
  real, intent(IN) :: xx, yy, zz

  real, external :: TSC_weight

  real :: vecPot

  real :: dx, dy, dz

  integer :: i, j, k, ii, jj, kk, ib, jb, kb

  vecPot = 0.0

  ii = int((xx-sim_Axmin)/sim_Adx) + 1
  jj = int((yy-sim_Aymin)/sim_Ady) + 1
  kk = int((zz-sim_Azmin)/sim_Adz) + 1

  ib = ii - ibegin + 1
  jb = jj - jbegin + 1
  kb = kk - kbegin + 1

  do i = -1, 1
     dx = (xx - sim_Axcoord(ii+i))/sim_Adx
     do j = -1, 1
        dy = (yy - sim_Aycoord(jj+j))/sim_Ady
        do k = -1, 1
           dz = (zz - sim_Azcoord(kk+k))/sim_Adz
           if (ib > 1 .and. ib < sim_nAx .and. &
               jb > 1 .and. jb < sim_nAy .and. & 
               kb > 1 .and. kb < sim_nAz) then
              if (axis == 1) then
                 vecPot = vecPot + sim_Ax(kb+k,jb+j,ib+i) * &
                      TSC_weight(dx)*TSC_weight(dy)*TSC_weight(dz)
              else if (axis == 2) then
                 vecPot = vecPot + sim_Ay(kb+k,jb+j,ib+i) * &
                      TSC_weight(dx)*TSC_weight(dy)*TSC_weight(dz)
              else if (axis == 3) then
                 vecPot = vecPot + sim_Az(kb+k,jb+j,ib+i) * &
                      TSC_weight(dx)*TSC_weight(dy)*TSC_weight(dz)
              endif
           endif
        enddo
     enddo
  enddo

  return

end function vecPot

function TSC_weight(x)

  implicit none

  real :: TSC_weight, xx
  real, intent(IN) :: x

  xx = abs(x)

  if (xx <= 0.5) then

     TSC_weight = 0.75 - xx*xx
     
  else if (xx >= 0.5 .and. xx <= 1.5) then

     TSC_weight = 0.5*(1.5-xx)*(1.5-xx)

  else
     
     TSC_weight = 0.0

  endif

  return

end function TSC_weight
