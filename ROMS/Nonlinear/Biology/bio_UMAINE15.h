      SUBROUTINE biology (ng,tile)
!
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  UMaine CoSiNE ecosystem model 13                                    !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields. The model is based on the        !
!  Carbon, Silicon, Nitrogen Ecosystem (CoSiNE) model (Chai et al.,    !
!  2002). The model state variables are:                               !
!                                                                      !
!    iNO3_      Nitrate                                                !  
!    iSiOH      Silicate                                               !
!    iNH4_      Ammonium                                               !
!    iSphy      Small Phytoplankton                                    !
!    iLphy      Diatoms                                                !
!    iSzoo      Micro Zooplankton                                      !
!    iLphy      Meso Zooplankton                                       !
!    iSDet      Detritus-nitrogen                                      !
!    iopal      Detritus-silicate                                      !
!    iPH4_      Phosphate                                              !
!    iOxyg      Dissolved Oxygen                                       !
!    ITIC_      Total CO2                                              !
!    iTAlk      Total Alkalinity                                       !
!    iChl1      Chlorophyll for Sphy                                   !
!    iChl2      Chlorophyll for diatom                                 !
#ifdef HAB
!    iHphy      HAB Phytoplankton                                      !
!    iChl3      Chlorophyll in HAB Phytoplankton                       !
#endif
!                                                                      !
!  Please cite:                                                        !
!                                                                      !
!    Chai, F., R.C. Dugdale, T-H Peng, F.P.Wilkerson, and R.T. Barber  !
!      (2002): One dimensional Ecosystem Model of the Equatorial       !
!      Pacific Upwelling System, Part I: Model Development and Silicon !
!      and Nitorgen Cycle. Deep-Sea Res. II, Vol. 49, No. 13-14,       !
!      2713-2745.                                                      !
!                                                                      !
!    Xiu, P., and F. Chai, 2011. Modeled biogeochemical responses to   !
!       mesoscale eddies in the South China Sea, Journal of Geophysical! 
!       Research, 116, C10006, doi:10.1029/2010JC006800.               !
!                                                                      !
!  Adapted from 1D code develped by Fei Chai of UMaine, and further    !
!  modified and developed by Lei Shi and adding 4 more biological      !
!  statevariables                                                      !
!                                                                      !
!  Last: PENG XIU 10/2016                                                !
!                                                                      !
!  Options you can use: OXYGEN, CARBON, SINK_OP1, SINK_OP2             !
!                       TALK_NONCONSERV,DIAGNOSTICS_BIO
!  DU: 8/9/2021: removed SINK_OP selection. using SINK_OP2 only         !
!  DU: 6/9/2023: adding 3rd phytoplankton group representing HAB phyto !
!***********************************************************************                                
!
      USE mod_param
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
#ifdef SEDBIO
      USE mod_sedbio
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                     &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),            &
     &                   IminS, ImaxS, JminS, JmaxS,                   &
     &                   nstp(ng), nnew(ng),                           &
#ifdef MASKING
     &                   GRID(ng) % rmask,                             &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                   GRID(ng) % rmask_io,                          &
# endif
#endif
     &                   GRID(ng) % Hz,                                &
     &                   GRID(ng) % H,                                 &
     &                   GRID(ng) % z_r,                               &
     &                   GRID(ng) % z_w,                               &
     &                   FORCES(ng) % srflx,                           &
#if defined OXYGEN || defined CARBON
# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                           &
     &                   FORCES(ng) % Vwind,                           &
# else
     &                   FORCES(ng) % sustr,                           &
     &                   FORCES(ng) % svstr,                           &
# endif
#endif
#ifdef CARBON
     &                   OCEAN(ng) % pH,                               &
     &                   OCEAN(ng) % rho,                              &
! AKB 6/22/23
     &                   OCEAN(ng) % fgCO2,                            &
     &                   OCEAN(ng) % pCO2s,                            &
#endif
#ifdef OPTICS_OP1
     &                   OCEAN(ng) % kdpar,                            &
# ifdef READ_ZENITH
     &                   FORCES(ng) % zenith_ang,                      &
# endif
#endif
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                         &
     &                   DIAGS(ng) % DiaBio3d,                         &
#endif
#ifdef SEDBIO
     &                   SEDBIOL(ng) % sedPOM,                          &
     &                   SEDBIOL(ng) % sedPoreWaterCon,                 &
     &                   SEDBIOL(ng) % sedDecayRate,                    &
     &                   SEDBIOL(ng) % sedFlux,                         &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                               &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,           &
     &                         IminS, ImaxS, JminS, JmaxS,             &
     &                         nstp, nnew,                             &
#ifdef MASKING
     &                         rmask,                                  &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                         rmask_io,                               &
# endif
#endif
     &                         Hz,hh, z_r, z_w, srflx,                 &
#if defined OXYGEN || defined CARBON
# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                           &
# else
     &                         sustr, svstr,                           &
# endif
#endif
#ifdef CARBON
     &                         pH, rho,                                &
! AKB 6/22/23
     &                         fgCO2,                                  &
     &                         pCO2s,                                  &
#endif
#ifdef OPTICS_OP1
     &                         kdpar,                                  &
# ifdef READ_ZENITH
     &                         zenith_ang,                             &
# endif
#endif
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                     &
#endif
#ifdef SEDBIO
     &                         sedPOM, sedPoreWaterCon,                &
     &                         sedDecayRate, sedFlux,                  &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mod_sources
      USE mod_sedbio
#ifdef SEDBIO
      USE sediment_biology_mod,  ONLY : sediment_bio
# ifdef DISTRIBUTE
	USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
#endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_io(LBi:,LBj:)
#  endif      
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: hh(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# if defined OXYGEN || defined CARBON
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
! AKB 6/22/23
      real(r8), intent(inout) :: fgCO2(LBi:,LBj:)
      real(r8), intent(inout) :: pCO2s(LBi:,LBj:)
# endif
# ifdef OPTICS_OP1
      real(r8), intent(inout) :: kdpar(LBi:,LBj:,:)
#  ifdef READ_ZENITH
      real(r8), intent(in) :: zenith_ang(LBi:,LBj:)
#  endif
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
# ifdef SEDBIO
      real(r8), intent(inout) :: sedPOM(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: sedPoreWaterCon(LBi:,LBj:,:)
      real(r8), intent(inout) :: sedDecayRate(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: sedFlux(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      
#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_io(LBi:UBi,LBj:UBj)
#  endif      
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: hh(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# if defined OXYGEN || defined CARBON
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,UBk)
! AKB 6/22/23
      real(r8), intent(inout) :: fgCO2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: pCO2s(LBi:UBi,LBj:UBj)
# endif
# ifdef OPTICS_OP1
      real(r8), intent(inout) :: kdpar(LBi:UBi,LBj:UBj,UBk)
#  ifdef READ_ZENITH
      real(r8), intent(in) :: zenith_ang(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
# ifdef SEDBIO
      real(r8), intent(inout) :: sedPOM(LBi:UBi,LBj:UBj,nspc,NPOM)
      real(r8), intent(inout) :: sedPoreWaterCon(LBi:UBi,LBj:UBj,NPWC)
      real(r8), intent(inout) :: sedDecayRate(LBi:UBi,LBj:UBj,nspc,NDR)
      real(r8), intent(inout) :: sedFlux(LBi:UBi,LBj:UBj,NSF)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
      
#endif
!
!  Local variable declarations.
!
#ifdef HAB
      integer, parameter :: Nsink = 8
#else
      integer, parameter :: Nsink = 6
#endif
#if defined OXYGEN || defined CARBON
      real(r8) :: u10squ, u10spd
#endif

      integer :: Iter, i, indx, isink, ibio, ivar, j, k, ks,is,ii,jj

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: Minval = 1.0e-6_r8
      real(r8), parameter :: zeptic = 100.0_r8

      real(r8) :: dtdays

      real(r8) :: cff, cff0, cff1

      real(r8), dimension(Nsink) :: Wbio

      real(r8), dimension(LBi:UBi) :: PARsur

#if defined OXYGEN || defined CARBON
      real(r8), dimension(LBi:UBi) :: kw660
#endif

#ifdef OXYGEN
      real(r8), dimension(LBi:UBi) :: o2sat
      real(r8), dimension(LBi:UBi) :: o2flx
#endif

#ifdef CARBON
      real(r8), dimension(LBi:UBi) :: co2flx
      real(r8), dimension(LBi:UBi) :: pco2sf
#endif

      integer, parameter :: mmax = 31
      real(r8), dimension(LBi:UBi,N(ng),NT(ng)) :: Bio
      real(r8), dimension(LBi:UBi,N(ng),NT(ng)) :: Bio_bak
#ifdef OPTICS_OP1
      real(r8), dimension(LBi:UBi,N(ng)) :: kd
#endif
! AKB Added new SEDBIO arguments on 3/20/23--------------------------

#ifdef SEDBIO
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bPOC
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bPON
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bPOP
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bPSi
      real(r8), dimension(LBi:UBi,LBj:UBj) :: bNO3
      real(r8), dimension(LBi:UBi,LBj:UBj) :: bNH4
      real(r8), dimension(LBi:UBi,LBj:UBj) :: bPO4
      real(r8), dimension(LBi:UBi,LBj:UBj) :: bSi
# ifdef CARBON
      real(r8), dimension(LBi:UBi,LBj:UBj) :: bTIC
#  ifdef TALK_NONCONSERV
      real(r8), dimension(LBi:UBi,LBj:UBj) :: bAlk
#  endif
# endif
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bUC
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bUN
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bUP
      real(r8), dimension(LBi:UBi,LBj:UBj,nspc) :: bUS
      real(r8), dimension(LBi:UBi,LBj:UBj) :: JNO3
      real(r8), dimension(LBi:UBi,LBj:UBj) :: JNH4
      real(r8), dimension(LBi:UBi,LBj:UBj) :: JPO4
      real(r8), dimension(LBi:UBi,LBj:UBj) :: JSi
      real(r8), dimension(LBi:UBi,LBj:UBj) :: SOD
# ifdef CARBON
      real(r8), dimension(LBi:UBi,LBj:UBj) :: JTIC
#  ifdef TALK_NONCONSERV
      real(r8), dimension(LBi:UBi,LBj:UBj) :: JAlk
#  endif
# endif
#endif
! AKB Added new SEDBIO arguments on 3/20/23--------------------------

      real(r8), dimension(LBi:UBi,N(ng)) :: Hz_inv
      real(r8), dimension(LBi:UBi,N(ng)) :: Hz_inv2
      real(r8), dimension(LBi:UBi,N(ng)) :: Hz_inv3
      real(r8), dimension(LBi:UBi,N(ng)) :: hzl

      real(r8), dimension(LBi:UBi,N(ng)+1) :: PIO
      real(r8), dimension(LBi:UBi,N(ng)) :: PAR
      real(r8), dimension(LBi:UBi,N(ng)) :: ADPT
#ifdef OPTICS_OP1
      real(r8) :: dt_optics
      real(r8), parameter :: kdIniVal = 0.05_r8
#endif
! parameters controling light inhibition of nitrification (W/m^2)
      real(r8), parameter :: II00 = 0.0095_r8
      real(r8), parameter :: KI_nitrif = 0.10_r8

      real(r8), dimension(N(ng)) :: sinkindx

      real(r8) :: cents1,cens1,cens2,cents2,cent1
      real(r8) :: thick,upo4s1,uco2s1,pnh4s2,OXR
	real(r8) :: q10pp,q10pr,q10br,q10od,q10zr,q10zg,q10ni
      real(r8) :: alts1, alts2, grows1, grows2,unh4s1,pnh4s1
      real(r8) :: uno3s2, unh4s2, uno3s1,usio4s2,upo4s2,uco2s2
      real(r8) :: xco2_in,gno3s1,gnh4s1,gno3s2,gnh4s2,gsio4s2
      real(r8) :: nps1,rps1,nps2,rps2,morts1,morts2,excrz1
      real(r8) :: excrz2,remvz2,nitrif,midd,middsi
      real(r8) :: gs1zz1,ro8,ro9,gs2zz2,gzz1zz2,gddzz2,gtzz2
      real(r8) :: Qsms1,Qsms2,Qsms3,Qsms4,Qsms5,Qsms6,Qsms7
      real(r8) :: Qsms8,Qsms9,Qsms10,Qsms11,Qsms12,Qsms13
      real(r8) :: NQsms1,NQsms2,NQsms3,NQsms4,NQsms5,NQsms6,NQsms7
      real(r8) :: NQsms8,NQsms9,NQsms10,NQsms11,NQsms12,NQsms13
      real(r8) :: sms1,sms2,sms3,sms4,sms5,sms6,sms7,SOC
      real(r8) :: sms8,sms9,sms10,sms11,sms12,sms13
      
      real(r8) :: npchl1,npchl2,gchl1zz1,gchl2zz2,morchl1,morchl2
      real(r8) :: Qsms14,Qsms15,NQsms14,NQsms15,sms14,sms15
      real(r8) :: Chl2C_s1,Chl2C_s2,excrz1_2,excrz2_2
      real(r8) :: aggregs1,aggregs2,aggregchl1,aggregchl2,sumphy
#ifdef PHYTO_RESP
      real(r8) :: resps1,resps2,respchl1,respchl2
#endif
#ifdef HAB
      real(r8) :: alts3,grows3,unh4s3,uno3s3,usio4s3,upo4s3,uco2s3
      real(r8) :: gno3s3,gnh4s3,gsio4s3,pnh4s3
      real(r8) :: nps3,rps3,morts3,cents3
      real(r8) :: npchl3,gchl3zz3,morchl3
      real(r8) :: Qsms16,Qsms17,NQsms16,NQsms17,sms16,sms17
      real(r8) :: Chl2C_s3
      real(r8) :: aggregs3,aggregchl3
      real(r8) :: gs3zz2,gchl3zz2
# ifdef PHYTO_RESP
      real(r8) :: resps3,respchl3
# endif
#endif
      
      real(r8) :: cffL, cffR, cu, dltL, dltR
      integer, dimension(LBi:UBi,N(ng)) :: ksource
      real(r8), dimension(LBi:UBi,0:N(ng)) :: FC
      real(r8), dimension(LBi:UBi,N(ng)) :: WL
      real(r8), dimension(LBi:UBi,N(ng)) :: WR
      real(r8), dimension(LBi:UBi,N(ng)) :: bL
      real(r8), dimension(LBi:UBi,N(ng)) :: bR
      real(r8), dimension(LBi:UBi,N(ng)) :: qc
#ifdef SEDBIO
      real(r8), dimension(LBi:UBi) :: FPOC
      real(r8), dimension(LBi:UBi) :: FPON
      real(r8), dimension(LBi:UBi) :: FPOP
      real(r8), dimension(LBi:UBi) :: FPSi
#endif

#include "set_bounds.h"

#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!

      IF (((iic(ng).gt.ntsDIA(ng)).and.                                &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                           &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.            &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN

        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
#ifdef OPTICS_OP1
!  Set time step for updating kdpar
      dt_optics=dt(ng)*optic_upd_fac(ng)/REAL(BioIter(ng),r8)
#endif
!
!  Set vertical sinking identification vector.
!
      idsink(1)=iLphy
      idsink(2)=iSDet
      idsink(3)=iopal
      idsink(4)=iSphy
      idsink(5)=iChl1
      idsink(6)=iChl2
#ifdef HAB
      idsink(7)=iHphy
      idsink(8)=iChl3
#endif
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wsp(ng)                ! iLphy
      Wbio(2)=wsd(ng)                ! iSDet
      Wbio(3)=wsdsi(ng)              ! iopal
      Wbio(4)=0.0_r8                ! iSphy
      Wbio(5)=0.0_r8                ! iSphy ichl1
      Wbio(6)=wsp(ng)                ! iLphy ichl2
#ifdef HAB
      Wbio(7)=wsp3(ng)              ! iHphy
      Wbio(8)=wsp3(ng)              ! iChl3
#endif
#ifdef SEDBIO
!  Extract sediment biology variables from full arrays
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,nspc                           ! nspc = 3: number of sed particulate classes
            bPOC(i,j,k)=sedPOM(i,j,k,ibPOC)     ! NPOM = ibPOC = 1: benthic POC conc
            bPON(i,j,k)=sedPOM(i,j,k,ibPON)     ! NPOM = ibPON = 2: benthic PON conc
            bPOP(i,j,k)=sedPOM(i,j,k,ibPOP)     ! NPOM = ibPOP = 3: benthic POP conc
            bPSi(i,j,k)=sedPOM(i,j,k,ibPSi)     ! NPOM = ibPSi = 4: benthic Si conc

            bUC(i,j,k)=sedDecayRate(i,j,k,ibUC) ! NDR = ibUC = 1: decay rate of POC
            bUN(i,j,k)=sedDecayRate(i,j,k,ibUN) ! NDR = ibUN = 2: decay rate PON
            bUP(i,j,k)=sedDecayRate(i,j,k,ibUP) ! NDR = ibUP = 3: decay rate POP
            bUS(i,j,k)=sedDecayRate(i,j,k,ibUS) ! NDR = ibUS = 4: decay rate PSi
          END DO

          bNO3(i,j)=sedPoreWaterCon(i,j,ibNO3)      ! NPWC = ibNO3 = 1
          bNH4(i,j)=sedPoreWaterCon(i,j,ibNH4)      ! NPWC = ibNH4 = 2
          bPO4(i,j)=sedPoreWaterCon(i,j,ibPO4)      ! NPWC = ibPO4 = 3
          bSi(i,j)=sedPoreWaterCon(i,j,ibSi)        ! NPWC = ibSi = 4
# ifdef CARBON
          bTIC(i,j)=sedPoreWaterCon(i,j,ibTIC)      ! NPWC = ibTIC = 5
#  ifdef TALK_NONCONSERV
          bAlk(i,j)=sedPoreWaterCon(i,j,ibAlk)      ! NPWC = ibAlk = 6
#  endif
# endif

          JNO3(i,j)=sedFlux(i,j,ibJNO3)         ! NSF = ibJNO3 = 1: benthic flux of NO3
          JNH4(i,j)=sedFlux(i,j,ibJNH4)         ! NSF = ibJNH4 = 2: benthic flux of NH4
          JPO4(i,j)=sedFlux(i,j,ibJPO4)         ! NSF = ibJPO4 = 3: benthic flux of PO4
          JSi(i,j)=sedFlux(i,j,ibJSi)           ! NSF = ibJSi  = 4: benthic flux of Si
          SOD(i,j)=sedFlux(i,j,ibSOD)           ! NSF = ibSOD  = 5: benthic oxygen demand
# ifdef CARBON
          JTIC(i,j)=sedFlux(i,j,ibJTIC)         ! NSF = ibJTIC = 6: benthic flux of TIC
#  ifdef TALK_NONCONSERV
          JAlk(i,j)=sedFlux(i,j,ibJAlk)         ! NSF - ibJAlk = 7: benthic flux of alkalinity
#  endif
# endif
        END DO
      END DO
#endif

      J_LOOP : DO j=Jstr,Jend
   
#ifdef SEDBIO
!  Initialize benthic flux vectors to zero
        DO i=Istr,Iend
            FPOC(i)=0.0
            FPON(i)=0.0
            FPOP(i)=0.0
            FPSi(i)=0.0
        END DO
#endif
!
!  Compute inverse thickness to avoid repeated divisions.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            hzl(i,k)=Hz(i,j,k)
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO ibio=1,NBT
          indx=idbio(ibio)
            
          DO k=1,N(ng)
            DO i=Istr,Iend
!        if( (indx.eq.iopal.or.indx.eq.iSiOH).and.t(i,j,k,nstp,indx).lt.0.0_r8) then
!            write(*,*)'< 0: indx = ',indx,'i=',i,'j=',j,'k=',k
!        end if
!              Bio_bak(i,k,indx)=MAX(t(i,j,k,nstp,indx),0.000001_r8)
		  Bio_bak(i,k,indx)=MAX(t(i,j,k,nstp,indx),0.000000_r8)
              Bio(i,k,indx)=Bio_bak(i,k,indx)
            END DO
          END DO
        END DO

#ifdef CARBON
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,iTIC_)=Bio_bak(i,k,iTIC_)
          END DO
        END DO
#endif
#ifdef OXYGEN
        DO k=1,N(ng)
          DO i=Istr,Iend
              Bio_bak(i,k,iOxyg)=MIN(Bio(i,k,iOxyg),800.0_r8)
              Bio(i,k,iOxyg)=Bio_bak(i,k,iOxyg)
           END DO
         END DO  
#endif
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
          END DO
        END DO

!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!=======================================================================
!
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  microphytoplankotn and diatom concentration within each grid box. 
!  Then, photosynthetically available radiation (PAR) are calculated
!  as the depth-averaged PAR within the vertical grid. Since the light
!  penetreated from surface to bottom, so the PAR are calculated in the 
!  same order. PARsur is surface PAR value. PIO is the PAR value at 
!  surface or bottom of a vertical grid, or at w location vertically.
!
#ifdef OPTICS_OP1
      IF( mod(time(ng),dt_optics).lt.0.000001_r8) then
!        write(*,*)'calling kdPAR_compute'
!        call optic_property(Istr, Iend, ng,                             &
        call kdPAR_compute(Istr, Iend, ng,                              &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, j,                           &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
#  ifdef READ_ZENITH
     &                       zenith_ang,                                &
#  endif
     &                       Bio(LBi:,:,isalt), hzl, z_r,               &
     &                       Bio(LBi:,:,iChl1), Bio(LBi:,:,iChl2),      &
     &                       Bio(LBi:,:,iSphy), Bio(LBi:,:,iLphy),      &
     &                       Bio(LBi:,:,iSDet),                         &
#  ifdef HAB
     &                       Bio(LBi:,:,iHphy), Bio(LBi:,:,iChl3),      &
#  endif
     &                       kd)

! store kdpar
        DO k=1,N(ng)
          DO i=Istr,Iend
            kdpar(i,j,k)=kd(i,k)
          END DO
        END DO

      ENDIF
#endif

    DO i=Istr,Iend
      PIO(i,N(ng)+1)=PARsur(i)
      IF (PIO(i,N(ng)+1) .lt. 0.0_r8) PIO(i,N(ng)+1)=0.0_r8
      DO k=N(ng),1,-1

#ifdef OPTICS_OP1
        cff1=kdpar(i,j,k)*HZ(i,j,k)
#else
#   ifdef HAB
        cff1=(AK1(ng)+(Bio(i,k,iChl1)+Bio(i,k,iChl2)+Bio(i,k,iChl3)*AK2(ng))*HZ(i,j,k)
#   else
        cff1=(AK1(ng)+(Bio(i,k,iChl1)+Bio(i,k,iChl2))*AK2(ng))*HZ(i,j,k)
#   endif
#endif
       
       PIO(i,K)=PIO(i,K+1)*EXP(-cff1)
! average PAR over layer:
       PAR(i,K)=(PIO(i,K+1)-PIO(i,K))/cff1
! variable ADPT only used if no chlorophyll variable used
       ADPT(i,K) = 1.0_r8-4.0_r8*z_r(i,j,k)/zeptic
      END DO
    END DO

       DO k=1,N(ng)
          DO i=Istr,Iend
!
!-----------------------------------------------------------------------
!     CALCULATING the temperature dependence of biology processes
!-----------------------------------------------------------------------
! compute q10 factors
!	Q10=exp(0.069_r8*(Bio(i,k,itemp)-15.0_r8))
      q10pp=exp(a1(ng)*(Bio(i,k,itemp)-15.0_r8))	! phyto production
      q10pr=exp(a2(ng)*(Bio(i,k,itemp)-15.0_r8))	! phyto respiration
      q10br=exp(a3(ng)*(Bio(i,k,itemp)-15.0_r8))	! bacterial remineralization
      q10od=exp(a4(ng)*(Bio(i,k,itemp)-15.0_r8))	! opal dissolution
      q10zr=exp(a5(ng)*(Bio(i,k,itemp)-15.0_r8))	! zoopl respiration
      q10zg=exp(a6(ng)*(Bio(i,k,itemp)-15.0_r8))	! zoopl grazing
      q10ni=exp(a7(ng)*(Bio(i,k,itemp)-15.0_r8))	! nitrification
!
!-----------------------------------------------------------------------
!     CALCULATING THE OXIDATION RATE OF ORGANIC MATTER
!-----------------------------------------------------------------------
!
!  Any biology processes that consume oxygen will be limited by the
!  availability of dissolved oxygen.
!
#ifdef OXYGEN
      OXR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg)+AKOX(ng))
#else
      OXR = 1.0_r8
#endif
!
!-----------------------------------------------------------------------
!     CALCULATING THE GROWTH RATE AS NO3,NH4, AND LIGHT;
!     GRAZING, PARTICLE SINKING AND REGENERATION
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  small phytoplankton nutrient uptake, and growth
!-----------------------------------------------------------------------

!  nitrogen limited growth
	pnh4s1=exp(-pis1(ng)*Bio(i,k,iNH4_))

!another option
!    pnh4s1=	1.0_r8/(1.0_r8+Bio(i,k,iNH4_)/aknh4s1(ng))
	
    uno3s1 = pnh4s1*Bio(i,k,iNO3_)/(akno3s1(ng)+Bio(i,k,iNO3_))
    unh4s1 = Bio(i,k,iNH4_)/(aknh4s1(ng)+Bio(i,k,iNH4_))

!  phosphate limited growth
	upo4s1=Bio(i,k,iPO4_)/(akpo4s1(ng)+Bio(i,k,iPO4_))
	
!  tco2 limited growth
#ifdef CARBON
	uco2s1=Bio(i,k,iTIC_)/(akco2s1(ng)+Bio(i,k,iTIC_))
#else
	uco2s1=1.0_r8
#endif

! light adaption if you dont have chloropyll variable
!    alts1= 1.0_r8 - exp(-PAR(i,k)*(ADPT(i,k)/40.0_r8)/gmaxs1(ng))
!

! if you have chloropyll variable

    alts1=1.0_r8 - exp(-PAR(i,k)*amaxs1(ng)/gmaxs1(ng))   

!    alts1=amaxs1(ng)*PAR(i,k)/sqrt(gmaxs1(ng)*gmaxs1(ng)+   &
!   &             amaxs1(ng)*PAR(i,k)*amaxs1(ng)*PAR(i,k))
   
      
! growth
    
    grows1=min(uno3s1+unh4s1,upo4s1)*alts1

!! adjustments
!
    cents1=grows1/(uno3s1+unh4s1+Minval)
    uno3s1=cents1*uno3s1
    unh4s1=cents1*unh4s1
    upo4s1=grows1
    uco2s1=grows1
!
!-----------------------------------------------------------------------
!  diatom nutrient uptake, and growth
!-----------------------------------------------------------------------
!
    pnh4s2= exp(-pis2(ng)*Bio(i,k,iNH4_))
!another option
!    pnh4s2=	1.0_r8/(1.0_r8+Bio(i,k,iNH4_)/aknh4s2(ng))
    
    uno3s2 = pnh4s2*Bio(i,k,iNO3_)/(akno3s2(ng)+Bio(i,k,iNO3_))
    unh4s2 = Bio(i,k,iNH4_)/(aknh4s2(ng)+Bio(i,k,iNH4_))
    usio4s2 = Bio(i,k,iSiOH)/(aksio4s2(ng)+Bio(i,k,iSiOH))
    upo4s2 = Bio(i,k,iPO4_)/(akpo4s2(ng)+Bio(i,k,iPO4_))
!  tco2 limited growth
#ifdef CARBON
    uco2s2 = Bio(i,k,iTIC_)/(akco2s2(ng)+Bio(i,k,iTIC_))
#else
	uco2s2=1.0_r8
#endif

! light adaption if you dont have chloropyll variable
!    alts2= 1.0_r8 - exp(-PAR(i,k)*(ADPT(i,k)/40.0_r8)/gmaxs2(ng))
!if you have chloropyll variable
    alts2=1.0_r8 - exp(-PAR(i,k)*amaxs2(ng)/gmaxs2(ng))
    
!    alts2=amaxs2(ng)*PAR(i,k)/sqrt(gmaxs2(ng)*gmaxs2(ng)+   &
!   &             amaxs2(ng)*PAR(i,k)*amaxs2(ng)*PAR(i,k))

! careful with this limitation term for diatoms

    grows2=min(uno3s2+unh4s2,usio4s2,upo4s2)*alts2
!    grows2=usio4s2*alts2

    cents2=grows2/(uno3s2+unh4s2+Minval)
    uno3s2=cents2*uno3s2
    unh4s2=cents2*unh4s2
    
    usio4s2=grows2
    upo4s2=grows2
    uco2s2=grows2
#ifdef HAB
!
!-----------------------------------------------------------------------
!  HAB phyto nutrient uptake, and growth
!-----------------------------------------------------------------------
!
    pnh4s3= exp(-pis3(ng)*Bio(i,k,iNH4_))
!another option
!    pnh4s3=    1.0_r8/(1.0_r8+Bio(i,k,iNH4_)/aknh4s3(ng))
                  
    uno3s3 = pnh4s3*Bio(i,k,iNO3_)/(akno3s3(ng)+Bio(i,k,iNO3_))
    unh4s3 = Bio(i,k,iNH4_)/(aknh4s3(ng)+Bio(i,k,iNH4_))
    usio4s3 = Bio(i,k,iSiOH)/(aksio4s3(ng)+Bio(i,k,iSiOH))
    upo4s3 = Bio(i,k,iPO4_)/(akpo4s3(ng)+Bio(i,k,iPO4_))
!  tco2 limited growth
# ifdef CARBON
    uco2s3 = Bio(i,k,iTIC_)/(akco2s3(ng)+Bio(i,k,iTIC_))
# else
    uco2s3=1.0_r8
# endif

! light adaption if you dont have chloropyll variable
!    alts2= 1.0_r8 - exp(-PAR(i,k)*(ADPT(i,k)/40.0_r8)/gmaxs2(ng))
!if you have chloropyll variable
    alts3=1.0_r8 - exp(-PAR(i,k)*amaxs3(ng)/gmaxs3(ng))
                  
!    alts3=amaxs3(ng)*PAR(i,k)/sqrt(gmaxs3(ng)*gmaxs3(ng)+   &
!   &             amaxs3(ng)*PAR(i,k)*amaxs3(ng)*PAR(i,k))

! careful with this limitation term for diatoms

    grows3=min(uno3s3+unh4s3,usio4s3,upo4s3)*alts3
!    grows3=usio4s3*alts3

    cents3=grows3/(uno3s3+unh4s3+Minval)
    uno3s3=cents3*uno3s3
    unh4s3=cents3*unh4s3
                  
    usio4s3=grows3
    upo4s3=grows3
    uco2s3=grows3
#endif
    
! uptake     
    gno3s1  = gmaxs1(ng)*uno3s1   
    gnh4s1  = gmaxs1(ng)*unh4s1
    gno3s2  = gmaxs2(ng)*uno3s2
    gnh4s2  = gmaxs2(ng)*unh4s2
    gsio4s2 = gmaxs2(ng)*usio4s2
#ifdef HAB
    gno3s3  = gmaxs3(ng)*uno3s3
    gnh4s3  = gmaxs3(ng)*unh4s3
    gsio4s3 = gmaxs3(ng)*usio4s3
#endif
              
!     -------------------------------------------------------
!     CALCULATING THE NEW,REGENERATED,AND PRIMARY PRODUCTION 
!     -------------------------------------------------------
      nps1 =  gno3s1 * Bio(i,k,iSphy) * q10pp
      rps1 =  gnh4s1 * Bio(i,k,iSphy) * q10pp
      nps2 =  gno3s2 * Bio(i,k,iLphy)  * q10pp
      rps2 =  gnh4s2 * Bio(i,k,iLphy) * q10pp
#ifdef PHYTO_RESP
!     phytoplankton respiration (basal plus growth related)
      resps1=(rrb1(ng)*Bio(i,k,iSphy) + rrg1(ng)*(nps1+rps1))*q10pr
      resps2=(rrb2(ng)*Bio(i,k,iLphy) + rrg2(ng)*(nps2+rps2))*q10pr
#endif
#ifdef HAB
      nps3 =  gno3s3 * Bio(i,k,iHphy) * q10pp
      rps3 =  gnh4s3 * Bio(i,k,iHphy) * q10pp
# ifdef PHYTO_RESP
!     phytoplankton respiration (basal plus growth related)
      resps3=(rrb3(ng)*Bio(i,k,iHphy) + rrg3(ng)*(nps3+rps3))*q10pr
# endif
#endif
!      
!     *** Chla for s1: Xiu and Geng ***  npchl1; npchl2; Chl2cs1_m(ng); Chl2cs2_m(ng)
!

      cff=c2n(ng)*12.0_r8
      Chl2C_s1=MIN(Bio(i,k,iChl1)/(Bio(i,k,iSphy)*cff+Minval),  &
     &                    Chl2cs1_m(ng))
     
      Chl2C_s2=MIN(Bio(i,k,iChl2)/(Bio(i,k,iLphy)*cff+Minval),  &
     &                    Chl2cs2_m(ng))
#ifdef HAB
      Chl2C_s3=MIN(Bio(i,k,iChl3)/(Bio(i,k,iHphy)*cff+Minval),  &
     &                    Chl2cs3_m(ng))
#endif
      npchl1=(gno3s1+gnh4s1)*(gno3s1+gnh4s1)*Chl2cs1_m(ng)      &
     &   *Bio(i,k,iChl1)/                             &
     &   (amaxs1(ng)*MAX(Chl2C_s1,Minval)*PAR(i,k)+Minval)

      npchl2=(gno3s2+gnh4s2)*(gno3s2+gnh4s2)*Chl2cs2_m(ng)     &
     &   *Bio(i,k,iChl2)/                             &
     &   (amaxs2(ng)*MAX(Chl2C_s2,Minval)*PAR(i,k)+Minval)
		  
      npchl1=npchl1*q10pp
      npchl2=npchl2*q10pp
#ifdef PHYTO_RESP
	respchl1=resps1*c2n(ng)*Chl2C_s1*12.0_r8*q10pr
      respchl2=resps2*c2n(ng)*Chl2C_s2*12.0_r8*q10pr	  
#endif
#ifdef HAB
      npchl3=(gno3s3+gnh4s3)*(gno3s3+gnh4s3)*Chl2cs3_m(ng)     &
     &   *Bio(i,k,iChl3)/                             &
     &   (amaxs3(ng)*MAX(Chl2C_s3,Minval)*PAR(i,k)+Minval)
      npchl3=npchl3*q10pp
# ifdef PHYTO_RESP
      respchl3=resps3*c2n(ng)*Chl2C_s3*12.0_r8*q10pr
# endif
#endif

!     -------------------------------------------------------
!     CALCULATING THE mortality of phyto and zoopl
!     -------------------------------------------------------
! use respiration q10s
      morts1=bgamma3(ng)*Bio(i,k,iSphy)*q10pr
      morts2=bgamma4(ng)*Bio(i,k,iLphy)*q10pr
#ifdef HAB
      morts3=bgamma4s3(ng)*Bio(i,k,iHphy)*q10pr
#endif
      remvz2 = bgamma0(ng)*Bio(i,k,iLzoo)*Bio(i,k,iLzoo)*q10zr

!     *** Chla for s1 and s2: Xiu and Geng
!
       morchl1=bgamma3(ng)*Bio(i,k,iChl1)*q10pr
       morchl2=bgamma4(ng)*Bio(i,k,iChl2)*q10pr
#ifdef HAB
       morchl3=bgamma4s3(ng)*Bio(i,k,iChl3)*q10pr
#endif
      
! aggregates
!  use respiration q10 
#ifdef HAB
      sumphy=Bio(i,k,iSphy)+Bio(i,k,iLphy)+Bio(i,k,iHphy)
      aggregs1=bgamma6(ng)*sumphy*Bio(i,k,iSphy)*q10pr
      aggregs2=bgamma6(ng)*sumphy*Bio(i,k,iLphy)*q10pr
      aggregs3=bgamma6(ng)*sumphy*Bio(i,k,iHphy)*q10pr
      aggregchl1=bgamma6(ng)*sumphy*Bio(i,k,iChl1)*q10pr
      aggregchl2=bgamma6(ng)*sumphy*Bio(i,k,iChl2)*q10pr
      aggregchl3=bgamma6(ng)*sumphy*Bio(i,k,iChl3)*q10pr
#else
      sumphy=Bio(i,k,iSphy)+Bio(i,k,iLphy)
      aggregs1=bgamma6(ng)*sumphy*Bio(i,k,iSphy)*q10pr
      aggregs2=bgamma6(ng)*sumphy*Bio(i,k,iLphy)*q10pr
      aggregchl1=bgamma6(ng)*sumphy*Bio(i,k,iChl1)*q10pr
      aggregchl2=bgamma6(ng)*sumphy*Bio(i,k,iChl2)*q10pr
#endif
!      -------------------------------------------------------
!     CALCULATING THE nitrification and reminalization
!     -------------------------------------------------------
! formulation in which light inhibits nitrification
        nitrif = bgamma7(ng)*(1.0_r8-max(0.0_r8,(PAR(i,k)-II00)/   &
     &       (KI_nitrif+PAR(i,k)-II00)))*Bio(i,k,iNH4_)
	nitrif=nitrif*q10ni
		  
! remineralization (nitrogenous detritus)
      cent1=bgamma5(ng)
      midd = cent1*Bio(i,k,iSDet)*q10br
      
! remineralization (silicious detritus)
      cent1=bgamma5s(ng)
      middsi = cent1*Bio(i,k,iopal)*q10od
           
!     -------------------------------------------------------
!     CALCULATING THE GRAZING RATE
!     -------------------------------------------------------
! type II, michaelis-menten:
!      gs1zz1 = beta1(ng)*Bio(i,k,iSzoo)*Bio(i,k,iSphy)                 &
!              &    /(akz1(ng)+Bio(i,k,iSphy))
! type III, michaelis-menten:
        gs1zz1 = beta1(ng)*Bio(i,k,iSzoo)*Bio(i,k,iSphy)*Bio(i,k,iSphy) &
     &    /(akz1(ng)*akz1(ng)+Bio(i,k,iSphy)*Bio(i,k,iSphy))

!     *** Chla for s1: Xiu and Geng
!              
!      gchl1zz1 = gs1zz1*Bio(i,k,iChl1)/Bio(i,k,iSphy)
	gchl1zz1 = gs1zz1*Bio(i,k,iChl1)/(Bio(i,k,iSphy)+Minval)
        
      gs1zz1=gs1zz1*q10zg
      gchl1zz1=gchl1zz1*q10zg
		  
! type III, michaelis-menten
#ifdef HAB
!      ro8=ro5(ng)*Bio(i,k,iLphy)+ro6(ng)*Bio(i,k,iSzoo)                &
!     & +ro7(ng)*Bio(i,k,iSDet) +ro5H(ng)*Bio(i,k,iHphy)
      ro9=ro5(ng)*Bio(i,k,iLphy)*Bio(i,k,iLphy)                        &
     & +ro6(ng)*Bio(i,k,iSzoo)*Bio(i,k,iSzoo)                          &
     & +ro7(ng)*Bio(i,k,iSDet)*Bio(i,k,iSDet)                          &
     & +ro5H(ng)*Bio(i,k,iHphy)*Bio(i,k,iHphy)

!      if((ro8 .le. 0.0_r8) .and. (ro9 .le. 0.0_r8))then
      IF( ro9 .le. 0.0_r8 )then
        gs2zz2  = 0.0_r8
        gs3zz2  = 0.0_r8
        gddzz2  = 0.0_r8
        gzz1zz2 = 0.0_r8
        gchl2zz2 = 0.0_r8    ! Xiu and Geng
        gchl3zz2 = 0.0_r8
      ELSE
        gs2zz2= beta2(ng)*ro5(ng)*Bio(i,k,iLphy)*Bio(i,k,iLphy)        &
     & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
        gs3zz2= beta2(ng)*ro5H(ng)*Bio(i,k,iHphy)*Bio(i,k,iHphy)       &
     & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
        gzz1zz2=beta2(ng)*ro6(ng)*Bio(i,k,iSzoo)*Bio(i,k,iSzoo)        &
     & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
        gddzz2=beta2(ng)*ro7(ng)*Bio(i,k,iSDet)*Bio(i,k,iSDet)         &
     & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
            
        !     *** Chla for s2 and s3: Xiu and Geng
        !
!        gchl2zz2 = gs2zz2*Bio(i,k,iChl2)/Bio(i,k,iLphy)
!        gchl3zz2 = gs3zz2*Bio(i,k,iChl3)/Bio(i,k,iHphy)
	  gchl2zz2 = gs2zz2*Bio(i,k,iChl2)/(Bio(i,k,iLphy)+Minval)
        gchl3zz2 = gs3zz2*Bio(i,k,iChl3)/(Bio(i,k,iHphy)+Minval)
      ENDIF
      
      gs2zz2=gs2zz2*q10zg
      gs3zz2=gs3zz2*q10zg
      gzz1zz2=gzz1zz2*q10zg
      gddzz2=gddzz2*q10zg
      gchl2zz2=gchl2zz2*q10zg
      gchl3zz2=gchl3zz2*q10zg
#else
!      ro8=ro5(ng)*Bio(i,k,iLphy)+ro6(ng)*Bio(i,k,iSzoo)                &
!              & +ro7(ng)*Bio(i,k,iSDet)
      ro9=ro5(ng)*Bio(i,k,iLphy)*Bio(i,k,iLphy)+ro6(ng)                &
            & *Bio(i,k,iSzoo)*Bio(i,k,iSzoo)                           &
            & +ro7(ng)*Bio(i,k,iSDet)*Bio(i,k,iSDet)

 !     if((ro8 .le. 0.0_r8) .and. (ro9 .le. 0.0_r8))then
      IF( ro9 .le. 0.0_r8 )then
           gs2zz2  = 0.0_r8
           gddzz2  = 0.0_r8
           gzz1zz2 = 0.0_r8
           gchl2zz2 = 0.0_r8    ! Xiu and Geng
      ELSE
      gs2zz2= beta2(ng)*ro5(ng)*Bio(i,k,iLphy)*Bio(i,k,iLphy)          &
           & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
      gzz1zz2=beta2(ng)*ro6(ng)*Bio(i,k,iSzoo)*Bio(i,k,iSzoo)          &
            & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
      gddzz2=beta2(ng)*ro7(ng)*Bio(i,k,iSDet)*Bio(i,k,iSDet)           &
            & *Bio(i,k,iLzoo)/(akz2(ng)*akz2(ng)+ro9)
            
            !
!     *** Chla for s2: Xiu and Geng
!
!      gchl2zz2 = gs2zz2*Bio(i,k,iChl2)/Bio(i,k,iLphy)
	  gchl2zz2 = gs2zz2*Bio(i,k,iChl2)/(Bio(i,k,iLphy)+Minval)
      ENDIF
		  
      gs2zz2=gs2zz2*q10zg
      gzz1zz2=gzz1zz2*q10zg
      gddzz2=gddzz2*q10zg
      gchl2zz2=gchl2zz2*q10zg
#endif

! critical values to stop grazing (shouldn't be needed for type III grazing)
!    IF(Bio(i,k,iSphy) .le. 0.005_r8)then
!        gs1zz1=0.0_r8
!        gchl1zz1=0.0_r8    ! Xiu and Geng
!    ENDIF
!    IF(Bio(i,k,iLphy) .le. 0.005_r8)then
!        gs2zz2=0.0_r8
!        gchl2zz2=0.0_r8    ! Xiu and Geng
!    ENDIF
          
! total grazing for large zoo

   gtzz2 =  gddzz2 + gzz1zz2 + gs2zz2
#ifdef HAB
!       IF(Bio(i,k,iHphy) .le. 0.005_r8) then
!              gs3zz2=0.0_r8
!              gchl3zz2=0.0_r8
!       ENDIF
       gtzz2=gtzz2 + gs3zz2
#endif

! Excretion by zooplankton
!  use q10 for respiration
!   basal excretion
       excrz1=reg1(ng)*Bio(i,k,iSzoo)*q10zr
       excrz2=reg2(ng)*Bio(i,k,iLzoo)*q10zr
               
! ingestion related excretion (from Fennel 2006)
!  this is a fraction of the grazing rate
       excrz1_2=bgamma1(ng)*reg1(ng)*gs1zz1/beta1(ng)
       excrz2_2=bgamma2(ng)*reg2(ng)*gtzz2/beta2(ng)
					     
       excrz1_2=excrz1_2*q10zr
       excrz2_2=excrz2_2*q10zr
              
!sediment o2 consumption
! moved below calculation of sinking flux and, if SEDBIO defined, calculation
!  of sediment oxygen demand.

!  this next expression is the one used in the provided code (DU 8/9/21)
!  SOC is now computed in subroutine sediment_bio (SOD)
!     SOC=0.614_r8*(2.0_r8**(Bio(i,k,itemp)/5.474_r8))
     
        Qsms1 = - nps1 - nps2 + OXR*nitrif
        Qsms3 = - rps1 - rps2 + OXR*excrz1 + OXR*excrz2             &
     &          - OXR*nitrif+OXR*midd+OXR*(excrz1_2+excrz2_2)
     
        Qsms4 = + nps1 + rps1 - gs1zz1 - morts1 - aggregs1 
        Qsms5 = + nps2 + rps2 - gs2zz2 - morts2 - aggregs2
        Qsms6 = + bgamma1(ng)*gs1zz1 - OXR*excrz1 - gzz1zz2 -OXR*excrz1_2
        Qsms7 = + bgamma2(ng)*gtzz2 - OXR*excrz2 - remvz2   -OXR*excrz2_2
        Qsms8 = + (1-bgamma1(ng))*gs1zz1 + (1-bgamma2(ng))*gtzz2       &
     &             - gddzz2+ morts1 + morts2- OXR*midd+ remvz2   &
     &             + aggregs1 + aggregs2
     
        Qsms2 = - (nps2 + rps2)*si2n(ng) + middsi
        Qsms9 = (gs2zz2 + morts2+aggregs2)*si2n(ng) - middsi
!        Qsms10= - (nps1+rps1+nps2+rps2)*p2n(ng)                        &
!     &          +OXR*(excrz1+excrz2)*p2n(ng)+ OXR*midd*p2n(ng)
! added effect of ingestion-related excretion
        Qsms10= - (nps1+rps1+nps2+rps2)*p2n(ng)                        &
     &          +OXR*(excrz1+excrz2+excrz1_2+excrz2_2)*p2n(ng)         &
     &          + OXR*midd*p2n(ng)
#ifdef PHYTO_RESP
        Qsms4 = Qsms4 - OXR*resps1
        Qsms5 = Qsms5 - OXR*resps2
        Qsms3 = Qsms3 + OXR*(resps1 + resps2)
        Qsms10 = Qsms10 + OXR*(resps1 + resps2)*p2n(ng)
!        Qsms2 = Qsms2 + OXR*(resps1+resps2)*si2n(ng)
        Qsms2 = Qsms2 + OXR*resps2*si2n(ng)
#endif
#ifdef HAB
!  note that detritus here does not get the sloppy-feeding part as this
!  is included in gtzz2 above. likewise excretion above includes the HAB part,
!  so no need to include here:
        Qsms1 = Qsms1 - nps3
        Qsms3 = Qsms3 - rps3
        Qsms8 = Qsms8 + morts3 + aggregs3
        Qsms2 = Qsms2 - (nps3 + rps3)*si2n(ng)
        Qsms9 = Qsms9 + (gs3zz2 + morts3 + aggregs3)*si2n(ng)
        Qsms10 = Qsms10 - (nps3 + rps3)*p2n(ng)
        Qsms16 = + nps3 + rps3 - gs3zz2 - morts3 - aggregs3     ! HAB
# ifdef PHYTO_RESP
        Qsms16 = Qsms16 - OXR*resps3
        Qsms3 = Qsms3 + OXR*resps3
        Qsms10 = Qsms10 + OXR*resps3*p2n(ng)
        Qsms2 = Qsms2 + OXR*resps3*si2n(ng)
# endif
#endif
#ifdef OXYGEN
!        Qsms11 = (nps1+nps2)*o2no(ng)+(rps1+rps2)*o2nh(ng)             &
!     &         -2.0_r8*OXR*nitrif - OXR*(excrz1 + excrz2)*o2nh(ng)     &
!     &         - OXR*midd*o2nh(ng)
! added effect of ingestion-related excretion
        Qsms11 = (nps1+nps2)*o2no(ng)+(rps1+rps2)*o2nh(ng)             &
     &         -2.0_r8*OXR*nitrif - OXR*(excrz1 + excrz2)*o2nh(ng)     &
     &         -OXR*(excrz1_2+excrz2_2)*o2nh(ng)- OXR*midd*o2nh(ng)
# ifdef PHYTO_RESP
        Qsms11 = Qsms11 - OXR*( resps1 + resps2 )*o2nh(ng)                          
# endif
#  ifdef HAB
        Qsms11 = Qsms11 + nps3*o2no(ng) + rps3*o2nh(ng)
#   ifdef PHYTO_RESP
        Qsms11 = Qsms11 - OXR*resps3*o2nh(ng)                          
#   endif
#  endif
#endif

#ifdef CARBON
!        Qsms12= - (nps1+rps1+nps2+rps2)*c2n(ng)+OXR*(excrz1+excrz2)    &
!     &          *c2n(ng)+ OXR*midd*c2n(ng)
! added effect of ingestion-related excretion
         Qsms12= - (nps1+rps1+nps2+rps2)*c2n(ng)+ OXR*midd*c2n(ng)    &
     &           +OXR*(excrz1+excrz2+excrz1_2+excrz2_2)*c2n(ng)
# ifdef PHYTO_RESP
         Qsms12 = Qsms12 + OXR*(resps1 +resps2)*c2n(ng)                            
# endif
        Qsms13= - Qsms1 + Qsms3
#  ifdef HAB
        Qsms12 = Qsms12 - (nps3 + rps3)*c2n(ng)
#    ifdef PHYTO_RESP
        Qsms12 = Qsms12 + OXR*resps3*c2n(ng)                            
#    endif
#  endif
#endif
!
!       Chla: Xiu and Geng
!

        if(Bio(i,k,iChl1) .ge. 0.01_r8) then
          Qsms14 = + npchl1 - morchl1 - gchl1zz1 - aggregchl1
#ifdef PHYTO_RESP
	    Qsms14 = Qsms14 - respchl1
#endif
        else
          Qsms14 = npchl1
        endif
        
        if(Bio(i,k,iChl2) .ge. 0.01_r8) then
          Qsms15 = + npchl2 - morchl2 - gchl2zz2 - aggregchl2
#ifdef PHYTO_RESP
	    Qsms15 = Qsms15 - respchl2
#endif
        else
          Qsms15 = npchl2
        endif
#ifdef HAB
        IF(Bio(i,k,iChl3) .ge. 0.01_r8) then
            Qsms17 = + npchl3 - morchl3 - gchl3zz2 - aggregchl3
# ifdef PHYTO_RESP
		Qsms17 = Qsms17 - respchl3
# endif
        ELSE
            Qsms17 = npchl3
        ENDIF
#endif
       
#ifdef DIAGNOSTICS_BIO
        DiaBio3d(i,j,k,iPPro1)=DiaBio3d(i,j,k,iPPro1)+                 &
# ifdef WET_DRY
     &            rmask_io(i,j)*                                       &
# endif
     &         (nps1 + rps1)*dtdays
        
        DiaBio3d(i,j,k,iPPro2)=DiaBio3d(i,j,k,iPPro2)+                 &
# ifdef WET_DRY
     &            rmask_io(i,j)*                                       &
# endif        
     &            (nps2 + rps2)*dtdays
        
        DiaBio3d(i,j,k,iNO3u)=DiaBio3d(i,j,k,iNO3u)+                   &
# ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
# endif        
     &              (nps1+nps2)*dtdays
              
# ifdef OXYGEN
        DiaBio3d(i,j,k,iO2pr)=DiaBio3d(i,j,k,iO2pr)+                   &
#  ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
#  endif
     &   ((nps1+nps2)*o2no(ng)+(rps1+rps2)*o2nh(ng))*dtdays
#  ifdef HAB
        DiaBio3d(i,j,k,iO2pr)=DiaBio3d(i,j,k,iO2pr)+                   &
#   ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
#   endif
     &   (nps3*o2no(ng)+rps3*o2nh(ng))*dtdays
#  endif
              
        DiaBio3d(i,j,k,initri)=DiaBio3d(i,j,k,initri)+                 &
#  ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
#  endif
     &       (-2.0_r8*OXR*nitrif)*dtdays
              
        DiaBio3d(i,j,k,iremin)=DiaBio3d(i,j,k,iremin)+                 &
#  ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
#  endif
     &       (-OXR*midd*o2nh(ng))*dtdays
        
        DiaBio3d(i,j,k,izoopl)=DiaBio3d(i,j,k,izoopl)+                 &
#  ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
#  endif
     &       (-OXR*(excrz1 + excrz2                                    &
     &       +excrz1_2 + excrz2_2)*o2nh(ng)*dtdays)
#  ifdef PHYTO_RESP    
        DiaBio3d(i,j,k,iphyres)=DiaBio3d(i,j,k,iphyres)+               &
#   ifdef WET_DRY
     &              rmask_io(i,j)*                                     &
#   endif 
#   ifdef HAB
     &       (-OXR*(resps1 + resps2 +resps3)*o2nh(ng)*dtdays)                        
#   else
     &       (-OXR*(resps1 + resps2)*o2nh(ng)*dtdays)
#   endif
#  endif
# endif
              
        IF ( (MOD(iic(ng),nDIA(ng)).eq.0) ) THEN
! divide diagnostic sums by the number of bio iterations.
              DO ivar=1,NDbio3d
                DiaBio3d(i,j,k,ivar)=DiaBio3d(i,j,k,ivar)/             &
     &            (REAL(BioIter(ng),r8))
              END DO
        END IF
#endif
              
! update values 
      
        Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+dtdays*Qsms1
        Bio(i,k,iSiOH)=Bio(i,k,iSiOH)+dtdays*Qsms2
        Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+dtdays*Qsms3
        Bio(i,k,iSphy)=Bio(i,k,iSphy)+dtdays*Qsms4
        Bio(i,k,iLphy)=Bio(i,k,iLphy)+dtdays*Qsms5
        Bio(i,k,iSzoo)=Bio(i,k,iSzoo)+dtdays*Qsms6
        Bio(i,k,iLzoo)=Bio(i,k,iLzoo)+dtdays*Qsms7
        Bio(i,k,iSDet)=Bio(i,k,iSDet)+dtdays*Qsms8
        Bio(i,k,iopal)=Bio(i,k,iopal)+dtdays*Qsms9
        Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+dtdays*Qsms10

#ifdef OXYGEN
        Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+dtdays*Qsms11
#endif
#ifdef CARBON
        Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+dtdays*Qsms12
# ifdef TALK_NONCONSERV
        Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+dtdays*Qsms13
# endif	
#endif
!
!      Chla : Xiu and Geng
!
        Bio(i,k,iChl1)=Bio(i,k,iChl1)+dtdays*Qsms14
        Bio(i,k,iChl2)=Bio(i,k,iChl2)+dtdays*Qsms15
#ifdef HAB
        Bio(i,k,iHphy)=Bio(i,k,iHphy)+dtdays*Qsms16
        Bio(i,k,iChl3)=Bio(i,k,iChl3)+dtdays*Qsms17
#endif
        END DO  !k loop
      END DO  !i loop

!other flux

#if defined OXYGEN || defined CARBON
!
!-----------------------------------------------------------------------
!     CALCULATING gas transfer velocity at a Schmidt number 
!-----------------------------------------------------------------------
!
      k=N(ng)
        DO i=Istr,Iend
!
!  Compute wind speed.
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
!
!  drag coefficient is 0.001, and air density is 1.2 kg per cube meter
!
            cff1=rho0/(0.001_r8*1.2_r8)
!convert wind stress to wind speed square

            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+          &
     &                 (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
            u10spd=sqrt(u10squ)
!
! Compute gas transfer velocity at a Schmidt number of 660.
! climatology wind speed (Wanninkhof & Mcgillis, 1999).
!      kw660=1.09*u10spd-0.333*u10squ+0.078*u10spd*u10squ      
! (in units of cm/hr), the one is too pronounced with large speed
!  short-term (<1 day) winds (Wanninkhof & Mcgillis, 1999).
!      kw660=0.0283*u10spd*u10squ                     
!( in units of cm/hr)
!
! the following is from Wanninkhof 1992:
!            kw660(i)=0.31*u10squ
! the following is from Wanninkhof 2014:
            kw660(i)=0.251_r8*u10squ
       END DO
#endif

#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!     CALCULATING THE O2 SURFACE saturation concentration and FLUX  
!-----------------------------------------------------------------------
!
      k=N(ng)
      CALL O2_flux (Istr, Iend, LBi, UBi, LBj, UBj,                      &
     &                     IminS, ImaxS, j,                            &
# ifdef MASKING
     &                     rmask,                                      &
# endif
     &                     Bio(LBi:,k,itemp), Bio(LBi:,k,isalt),   &
     &                     Bio(LBi:,k,iOxyg), kw660,                 &
     &                     1.0_r8, o2sat, o2flx)     
      DO i=Istr,Iend
        Bio(i,k,iOxyg)=bio(i,k,iOxyg)+dtdays*o2flx(i)*Hz_inv(i,k)
    
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iO2fx)=DiaBio2d(i,j,iO2fx)+                   &
#  ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#  endif
     &                          o2flx(i)*dtdays
# endif

      END DO
    
#endif

#ifdef CARBON
!
!-----------------------------------------------------------------------
!  CALCULATING
!  Surface equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface, and CO2 gas exchange.
!  CO2 flux is downward/positive into the water.
!-----------------------------------------------------------------------
!
      k=N(ng)

! *************************************************************************
! AIR-SEA CO2 FLUX
! Returns:
!     - CO2 downward surface flux in mmol/m3/day
!     - surface pCO2 in ppmv
! *************************************************************************

      CALL CO2_flux (Istr, Iend, LBi, UBi, LBj, UBj,                   &
     &                     IminS, ImaxS, j,                            &
# ifdef MASKING
     &                     rmask,                                      &
# endif
     &                     Bio(LBi:,k,itemp), Bio(LBi:,k,isalt),   &
     &                     Bio(LBi:,k,iTIC_), Bio(LBi:,k,iTAlk),   &
     &                     Bio(LBi:,k,iPO4_), Bio(LBi:,k,iSiOH),   &
     &                     rho(LBi:,j,k),kw660, 1.0_r8,pco2a(ng),  &
     &                     co2flx,pco2sf)
      DO i=Istr,Iend
        Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+dtdays*co2flx(i)*Hz_inv(i,k)
! AKB 6/22/23
! Added pCO2 and pH outputs on 7/23/23 - AKB
        fgCO2(i,j)=co2flx(i)
        pCO2s(i,j)=pco2sf(i)
	  
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iCOfx)=DiaBio2d(i,j,iCOfx)+                   &
#  ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#  endif
     &                          co2flx(i)*dtdays
     
!            DiaBio2d(i,j,ipCO2)=pco2s(i)
              DiaBio2d(i,j,ipCO2)=pco2sf(i)
#  ifdef WET_DRY
            DiaBio2d(i,j,ipCO2)=DiaBio2d(i,j,ipCO2)*rmask_io(i,j)
#  endif
# endif
	  
      END DO
#endif
! 
!-----------------------------------------------------------------------
!     CALCULATING THE SINKING FLUX
!-----------------------------------------------------------------------

!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
    SINK_LOOP: DO isink=1,Nsink
            indx=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,indx)
!                if( i.eq.309.and.j.eq.387.and.indx.eq.11 ) then
!		    	write(*,*)'prior: i=',i,'j=',j,'k=',k,'bio=',Bio(i,k,indx)
!		    end if
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                   &
     &                    Hz(i,j,ks)*cu*                               &
     &                    (bL(i,ks)+                                   &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-             &
     &                         (1.5_r8-cu)*                            &
     &                         (bR(i,ks)+bL(i,ks)-                     &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,indx)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
!                if( i.eq.309.and.j.eq.387.and.indx.eq.11 ) then
!		    	write(*,*)'after: i=',i,'j=',j,'k=',k,'bio=',Bio(i,k,indx)
!		    end if
              END DO
            END DO
#  ifdef SEDBIO
!
!   Get fluxes to sediment
          DO i=Istr,Iend
!            IF( isink.eq.1 ) THEN
            IF( idsink(isink).eq.iLphy ) then
                ! diatoms
                FPON(i)=FPON(i)+FC(i,0)/dtdays
                FPSi(i)=FPSi(i)+FC(i,0)*si2n(ng)/dtdays
!            ELSE IF( isink.eq.2 ) THEN
            ELSE IF( idsink(isink).eq.iSDet ) THEN
               ! nitrogenous detritus
                FPON(i)=FPON(i)+FC(i,0)/dtdays

!            ELSE IF( isink.eq.3 ) THEN
            ELSE IF( idsink(isink).eq.iopal ) THEN
               ! siliceous detritus
               FPSi(i)=FPSi(i)+FC(i,0)/dtdays
            ELSE IF( idsink(isink).eq.iSphy ) THEN
               !small phyto
               FPON(i)=FPON(i)+FC(i,0)/dtdays
#    ifdef HAB
            ELSE IF( idsink(isink).eq.iHphy ) THEN
              ! HAB phyto (silica-containing)
              FPON(i)=FPON(i)+FC(i,0)/dtdays
              FPSi(i)=FPSi(i)+FC(i,0)*si2n(ng)/dtdays
#    endif
            END IF
            
          END DO
#  endif
    END DO SINK_LOOP
      
        DO i=Istr,Iend
            ! compute carbon, phosporous flux from nitrogen flux
            FPOC(i)= FPON(i)*c2n(ng)
            FPOP(i)= FPON(i)*p2n(ng)
                                   
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,ibPONfx)=DiaBio2d(i,j,ibPONfx)+               &
#  ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#  endif
     &                          FPON(i)*dtdays
                                   
            DiaBio2d(i,j,ibPSifx)=DiaBio2d(i,j,ibPSifx)+               &
#  ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#  endif
     &                          FPSi(i)*dtdays
# endif
                                   
        END DO
              
# ifdef SEDBIO
	  
        call sediment_bio(ng, Istr, Iend, LBi, UBi,                    &
     &                      dtdays,                                    &
#   ifdef MASKING
     &                      rmask(LBi:UBi,j),                          &
#   endif
     &                      Bio(LBi:,1,itemp),Bio(LBi:,1,iOxyg),       &
     &                      Bio(LBi:,1,iNO3_),Bio(LBi:,1,iNH4_),       &
     &                      Bio(LBi:,1,iPO4_),Bio(LBi:,1,iSiOH),       &
     &                      FPOC, FPON, FPOP, FPSi,                    &
     &                      bPOC(LBi:,j,:), bPON(LBi:,j,:),            &
     &                      bPOP(LBi:,j,:), bPSi(LBi:,j,:),            &
     &                      bNO3(LBi:,j), bNH4(LBi:,j),                &
     &                      bPO4(LBi:,j), bSi(LBi:,j),                 &
     &                      bUC(LBi:,j,:), bUN(LBi:,j,:),              &
     &                      bUP(LBi:,j,:), bUS(LBi:,j,:),              &
     &                      JNO3(LBi:,j), JNH4(LBi:,j),                &
     &                      JPO4(LBi:,j), JSi(LBi:,j), SOD(LBi:,j),    &
#   ifdef CARBON
    &                      Bio(LBi:,1,iTIC_), bTIC(LBi:,j),            &
    &                      JTIC(LBi:,j),                               &
#    ifdef TALK_NONCONSERV
    &                      Bio(LBi:,1,iTAlk), bAlk(LBi:,j),            &
    &                      JAlk(LBi:,j))
#    endif
#   endif

! WRITTEN IN SEDIMENT_BIOLOGY.F AS FOLLOWS-------------------------------
! SUBROUTINE sediment_bio(ng, Istr, Iend, LBi, UBi,                &    
! &                      dtdays,                                   &    
!   ifdef MASKING
!     &                      rmask,                                     &    
!   endif
!     &                      temp, DOX, dNO3, dNH4, dPO4, dSiO4,        &    
!     &                      FPOC, FPON, FPOP, FPSi,                    &    
!     &                      bPOC, bPON, bPOP, bPSi,                    &    
!     &                      bNO3, bNH4, bPO4, bSi,                     &    
!     &                      bUC, bUN, bUP, bUS,                        &    
!     &                      JNO3, JNH4, JPO4, JSi, SOD,                &    
!ifdef CARBON
!     &                      dTIC, bTIC, JTIC,                          &    
!ifdef TALK_NONCONSERV
!     &                      dAlk, bAlk, JAlk)
!endif
!endif
! -----------------------------------------------------------------------
! ORDER OF ARGUMENTS:           temp, oxygen, dNO3, dNH4, dPO4, dSiO4, FPOC, FPON,
!                               FPOP, FPSi, bPOC, bPON, bPOP, bPSi, bNO3, bNH4,
!                               bPO4, bSi, bUC, bUN, bUP, bUS, JNO3, JNH4, JPO4 
!                               JSi, SOD
! ADDED ARGUMENTS (3/20/23):    dTIC, bTIC, JTIC, dAlk, bAlk, JAlK
 
!
! update near-bottom concentrations in response to fluxes
!
! Near bottom variables: iOxyg, iNO3_, iNH4_, iPO4_, iSiOh, iTIC_, iTAlk
!                        where k = 1 for bottom sigma level
! Form: [conc]_new = [conc]_old + dt * flux / H_sed
! Negative for SOD since SOD is a comsumption term

        k=1
        DO i=Istr,Iend
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-dtdays*SOD(i,j)*Hz_inv(i,k)
            Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+dtdays*JNO3(i,j)*Hz_inv(i,k)
            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+dtdays*JNH4(i,j)*Hz_inv(i,k)
            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+dtdays*JPO4(i,j)*Hz_inv(i,k)
            Bio(i,k,iSiOH)=Bio(i,k,iSiOH)+dtdays*JSi(i,j)*Hz_inv(i,k)
            ! TIC and Alk added on 2/28/23 by A. Baskind
            Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+dtdays*JTIC(i,j)*Hz_inv(i,k)
            Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+dtdays*JAlk(i,j)*Hz_inv(i,k)
        
! Diagnostics
#  ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,ibNO3fx)=DiaBio2d(i,j,ibNO3fx)+               &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          JNO3(i,j)*dtdays
            DiaBio2d(i,j,ibNH4fx)=DiaBio2d(i,j,ibNH4fx)+               &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          JNH4(i,j)*dtdays
            DiaBio2d(i,j,ibPO4fx)=DiaBio2d(i,j,ibPO4fx)+               &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          JPO4(i,j)*dtdays
            DiaBio2d(i,j,ibSiO2fx)=DiaBio2d(i,j,ibSiO2fx)+             &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          JSi(i,j)*dtdays
            DiaBio2d(i,j,ibTICfx)=DiaBio2d(i,j,ibTICfx)+               &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          JTIC(i,j)*dtdays
            DiaBio2d(i,j,ibAlkfx)=DiaBio2d(i,j,ibAlkfx)+               &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          JAlk(i,j)*dtdays
! minus sign for SOD:
            DiaBio2d(i,j,ibO2fx)=DiaBio2d(i,j,ibO2fx)-                 &
#   ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#   endif
     &                          SOD(i,j)*dtdays

            IF ( (MOD(iic(ng),nDIA(ng)).eq.0) ) THEN
            ! divide 2d diagnostic sums by the number of bio iterations.
                DO ivar=1,NDbio2d
                    DiaBio2d(i,j,ivar)=DiaBio2d(i,j,ivar)/             &
     &            (REAL(BioIter(ng),r8))
                END DO
            END IF
#  endif                                   
        END DO
# else
! for non sediment bio case, estimate sediment oxygen consumption as !
! function of temperature (as was done in v1 of bio_UMAINE15)
! included multiplication by Q10 factor to be consistent with v1.
! use q10 for phyto respiration
        k=1
        DO i=Istr,Iend
!            Q10=exp(0.069_r8*(Bio(i,k,itemp)-25.0_r8))
            SOC=0.614_r8*(2.0_r8**(Bio(i,k,itemp)/5.474_r8))
            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+dtdays*0.036_r8*SOC*q10pr
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-dtdays*SOC*q10pr
#  ifdef CARBON
            Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+dtdays*0.9475_r8*SOC*q10pr
#  endif
        END DO
# endif
	
      END DO ITER_LOOP

!
!-----------------------------------------------------------------------
!  Update global tracer variables.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012)
!-----------------------------------------------------------------------

        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
!              t(i,j,k,nnew,indx)=max(t(i,j,k,nnew,indx)+                   &
!     &            (Bio(i,k,indx)-Bio_bak(i,k,indx))*  Hz(i,j,k),0.00001_r8)
              t(i,j,k,nnew,indx)=t(i,j,k,nnew,indx)+                   &
     &            (Bio(i,k,indx)-Bio_bak(i,k,indx))*  Hz(i,j,k)

! Note that following test for NaN only works if debugging is on
!              if( isnan(t(i,j,k,nnew,indx))) then
!                write(*,*)' time = ',time(ng),' indx = ',indx,'NAN found'
!              end if

              
            END DO
          END DO
        END DO

      END DO J_LOOP

#ifdef SEDBIO
!
! Update sediment biology variables
      DO j=Jstr,Jend
        DO i=Istr,Iend
            DO k=1,nspc
              sedPOM(i,j,k,ibPOC)=bPOC(i,j,k)
              sedPOM(i,j,k,ibPON)=bPON(i,j,k)
              sedPOM(i,j,k,ibPOP)=bPOP(i,j,k)
              sedPOM(i,j,k,ibPSi)=bPSi(i,j,k)
              sedDecayRate(i,j,k,ibUC)=bUC(i,j,k)
              sedDecayRate(i,j,k,ibUN)=bUN(i,j,k)
              sedDecayRate(i,j,k,ibUP)=bUP(i,j,k)
              sedDecayRate(i,j,k,ibUS)=bUS(i,j,k)
            END DO
            sedPoreWaterCon(i,j,ibNO3)=bNO3(i,j)
            sedPoreWaterCon(i,j,ibNH4)=bNH4(i,j)
            sedPoreWaterCon(i,j,ibPO4)=bPO4(i,j)
            sedPoreWaterCon(i,j,ibSi)=bSi(i,j)
            sedFlux(i,j,ibJNO3)=JNO3(i,j)
            sedFlux(i,j,ibJNH4)=JNH4(i,j)
            sedFlux(i,j,ibJPO4)=JPO4(i,j)
            sedFlux(i,j,ibJSi)=JSi(i,j)
            sedFlux(i,j,ibSOD)=SOD(i,j)
# ifdef CARBON
            sedPoreWaterCon(i,j,ibTIC)=bTIC(i,j)    ! Added 2/28/23 by A.Baskind
            sedFlux(i,j,ibJTIC)=JTIC(i,j)           ! Added 2/28/23 by A.Baskind
#  ifdef TALK_NONCONSERV
            sedPoreWaterCon(i,j,ibAlk)=bAlk(i,j)    ! Added 2/28/23 by A.Baskind
            sedFlux(i,j,ibJAlk)=JAlk(i,j)           ! Added 2/28/23 by A.Baskind
#  endif
# endif
        END DO
      END DO
# ifdef DISTRIBUTE
!
! exchange sediment biology variables between tiles to fill in ghost points
	call mp_exchange3d(ng, tile, iNLM, 1,					&
     &			 LBi, UBi, LBj, UBj, 1, NPWC,				&
     &			 NghostPoints, EWperiodic(ng), NSperiodic(ng),  &
     &			 sedPoreWaterCon)

	call mp_exchange3d(ng, tile, iNLM, 1,					&
     &			 LBi, UBi, LBj, UBj, 1, NSF,				&
     &			 NghostPoints, EWperiodic(ng), NSperiodic(ng),  &
     &			 sedFlux)
     
      call mp_exchange4d(ng, tile, iNLM, 1,					&
     &			 LBi, UBi, LBj, UBj, 1, nspc, 1, NDR,		&
     &			 NghostPoints, EWperiodic(ng), NSperiodic(ng),  &
     &			 sedDecayRate)
     
      call mp_exchange4d(ng, tile, iNLM, 1,					&
     &			 LBi, UBi, LBj, UBj, 1, nspc, 1, NPOM,		&
     &			 NghostPoints, EWperiodic(ng), NSperiodic(ng),  &
     &			 sedPOM)
     
     call mp_exchange3d(ng, tile, iNLM, 1,					&
     &			 LBi, UBi, LBj, UBj, 1, N(ng),	   	      &
     &			 NghostPoints, EWperiodic(ng), NSperiodic(ng),  &
     &			 kdpar)
# endif   
#endif
      RETURN
      
    END SUBROUTINE biology_tile
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Subroutines used in the model
!

#ifdef CARBON
        subroutine   CO2_flux (Istr, Iend, LBi, UBi, LBj, UBj,         &
     &                     IminS, ImaxS, j,                            &
# ifdef MASKING
     &                     rmask,                                      &
# endif
     &                     t, ss,dic, alk,po4,si,rho,kw660,            &
     &                     ppo, xco2,                                  &
     &                     co2ex,pco2sf)
!
!**********************************************************************
!
!  Computes the time rate of change of DIC in the surface
!  layer due to air-sea gas exchange in mmol/m^3/day.
!
!  Inputs:
!    t        model surface temperature (deg C)
!    ss        model surface salinity (permil)
!    kw660    gas transfer velocity at a Schmidt number of 660,
!               accounting for sea ice fraction (cm/hr)
!    rho      surface water density anomaly (kg/m^3)
!    ppo      surface pressure divided by 1 atm
!    dic      surface DIC concentration (mol/m^3)
!    alk      surface alkalinity (eq/m^3)
!    po4      surface phosphate concentration (mol/m^3)
!    si       surface silicate concentration (mol/m^3)
!    xco2     atmospheric CO2 mixing ratio (ppm)
!  Output:
!    co2ex    time rate of change of DIC in the surface layer due
!               to air-sea exchange (mmol/m^3/day)
!    pco2sf   surface pCO2 [ppmv]
!c**********************************************************************
      
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
      real(r8),  intent(in) :: ppo,xco2
      integer :: i
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: t(LBi:)
      real(r8), intent(in) :: ss(LBi:)
      real(r8), intent(in) :: dic(LBi:)
      real(r8), intent(in) :: alk(LBi:)
      real(r8), intent(in) :: po4(LBi:)
      real(r8), intent(in) :: si(LBi:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: t(LBi:UBi)
      real(r8), intent(in) :: ss(LBi:UBi)
      real(r8), intent(in) :: dic(LBi:UBi)
      real(r8), intent(in) :: alk(LBi:UBi)
      real(r8), intent(in) :: po4(LBi:UBi)
      real(r8), intent(in) :: si(LBi:UBi)
#  endif

      real(r8), intent(in) :: kw660(LBi:UBi)
      real(r8), intent(in) :: rho(LBi:UBi)
      real(r8), intent(out) :: co2ex(LBi:UBi)
      real(r8), intent(out) :: pco2sf(LBi:UBi)
      real(r8) :: s(LBi:UBi)
      real(r8) :: scco2,kwco2,phlo,phhi
      real(r8) :: co2star,dco2star,pCO2surf,dpco2,ph
      real(r8) :: dic2,alk2,po42,si2,rhos
      real(r8) :: A,B,C,D,E,tsq
! schmidt number polynomial coefficients (freshwater):
      real(r8), parameter :: A0 = 1923.6_r8
      real(r8), parameter :: B0 = -125.06_r8
      real(r8), parameter :: C0 = 4.3773_r8
      real(r8), parameter :: D0 = -0.085681_r8
      real(r8), parameter :: E0 = 0.00070284_r8
! schmidt number polynomial coefficients (35 psu):
      real(r8), parameter :: A35 = 2116.8_r8
      real(r8), parameter :: B35 = -136.25_r8
      real(r8), parameter :: C35 = 4.7353_r8
      real(r8), parameter :: D35 = -0.092307_r8
      real(r8), parameter :: E35 = 0.0007555_r8
      
       I_LOOP: DO i=Istr,Iend

         dic2=dic(i)/1000.0_r8
         alk2=alk(i)/1000.0_r8
         po42=po4(i)/1000.0_r8
         si2 = si(i)/1000.0_r8
         
!         !not work well for low s
!         IF (ss(i) .lt. 5.0_r8) then
!            s(i)=5.0_r8
!         ELSE
            s(i)=ss(i)
!         ENDIF
         
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
! third order polynomial for Schmidt No. from Wanninkhof 1992:
!      scco2 = 2073.1_r8 - 125.62_r8*t(i) + 3.6276_r8*t(i)*t(i)         &
!     & - 0.043219_r8*t(i)*t(i)*t(i)
! 4th order polynomial from Wanninkhof 2014 (for 35 psu):
!      scco2 = 2116.8_r8 - 136.25_r8*t(i) + 4.7353_r8*t(i)*t(i)         &
!     & - 0.092307_r8*t(i)*t(i)*t(i) + 0.0007555_r8*t(i)*t(i)*t(i)*t(i)
! interpolate schmidt number polynomial coefficients using surface salinity:
      A=A0+ss(i)*(A35-A0)/35.0_r8
      B=B0+ss(i)*(B35-B0)/35.0_r8
      C=C0+ss(i)*(C35-C0)/35.0_r8
      D=D0+ss(i)*(D35-D0)/35.0_r8
      E=E0+ss(i)*(E35-E0)/35.0_r8
      tsq=t(i)*t(i)
! schmidt number (4th order polynomial)
      scco2= A + B*t(i) + C*tsq + D*tsq*t(i) + E*tsq*tsq

      kwco2 = Kw660(i) * (660.0_r8/scco2)**0.5_r8             
!     (in units of cm/hr)
!  Compute the transfer velocity for CO2 in m/day

      kwco2=kwco2*0.01_r8*24.0_r8                        

      phlo = 6.0_r8
!      phhi = 9.0_r8
      phhi = 11.0_r8
! surface density (kg/m^3)
      rhos=1000_r8+rho(i)

      CALL co2calc(t(i),s(i),dic2,alk2,po42,si2,phlo,phhi,             &
     & 	   xco2,ppo,co2star,dco2star,pCO2surf,dpco2,ph,rhos)

      co2ex(i) = kwco2*dco2star

!  Compute time rate of change of CO2 due to gas exchange [1] in 
!  mmol/m^3/day.
!  Return surface pCO2 in ppmv.

      co2ex(i) = 1000.0_r8*co2ex(i)
      pco2sf(i) = pCO2surf
 
!    write(*,*) '++++++++++++++++++++++'
!    write(*,*)t(i),s(i),dic(i),alk(i),po4(i),si(i),xco2,pCO2surf,     &
!   & ph,co2ex(i)

#  ifdef MASKING
      ELSE
        co2ex(i)=0.0_r8
        pco2sf(i) = 0.0_r8
      END IF
#  endif

      END DO I_LOOP
      RETURN
      END SUBROUTINE CO2_flux
      
!------------------------------------------------   
      
       subroutine co2calc(t,s,dic_in,ta_in,pt_in,sit_in                &
     &                  ,phlo,phhi,xco2_in,atmpres                     &
     &                  ,co2star,dco2star,pCO2surf,dpco2,ph,rhos)
      USE mod_kinds
      implicit none
!**********************************************************************
!
! SUBROUTINE CO2CALC
!
! PURPOSE
!      Calculate delta co2* from total alkalinity and total CO2 at
! temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
!
! USAGE
!       call co2calc(t,s,dic_in,ta_in,pt_in,sit_in
!    &                  ,phlo,phhi,ph,xco2_in,atmpres
!    &                  ,co2star,dco2star,pCO2surf,dpco2,ph,rhos)
!
! INPUT
!      dic_in = total inorganic carbon (mol/m^3) 
!                where 1 T = 1 metric ton = 1000 kg
!      ta_in  = total alkalinity (eq/m^3) 
!      pt_in  = inorganic phosphate (mol/m^3) 
!      sit_in = inorganic silicate (mol/m^3) 
!      t      = temperature (degrees C)
!      s      = salinity (PSU)
!      phlo   = lower limit of pH range
!      phhi   = upper limit of pH range
!      xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
!      atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
!      ph = pH
!      rhos = in situ surface density (kg/m^3)
!
!       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
!             used to initialize variables dic, ta, pt, sit, and xco2.
!             * Variables dic, ta, pt, and sit are in the common block 
!               "species".
!             * Variable xco2 is a local variable.
!             * Variables with "_in" suffix have different units 
!               than those without.
! OUTPUT
!      co2star  = CO2*water (mol/m^3)
!      dco2star = delta CO2 (mol/m^3)
!       pco2surf = oceanic pCO2 (ppmv)
!       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
!
! IMPORTANT: Some words about units - (JCO, 4/4/1999)
!  - Models carry tracers in mol/m^3 (on a per volume basis)
!  - Conversely, this routine, which was written by observationalists 
!       (C. Sabine and R. Key), passes input arguments in umol/kg  
!       (i.e., on a per mass basis)
!  - I have changed things slightly so that input arguments are in 
!     mol/m^3,
!  - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
!    should be given in mol/m^3; output arguments "co2star" and 
!    "dco2star" are likewise in mol/m^3.
!**********************************************************************
       real(r8),intent(in)  :: t,s,dic_in,ta_in,pt_in,sit_in
       real(r8),intent(in)  :: phlo,phhi,xco2_in,atmpres,rhos
       real(r8),intent(out) :: co2star,dco2star,pCO2surf,dpco2,ph
!
!  Local variable declarations.
!       
       real(r8) :: invtk,is,is2,bt,st,ft,sit,pt,dic,ta
       real(r8) :: k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,htotal
       real(r8) :: permil,permeg,xco2,tk,tk100,tk1002,dlogtk,sqrtis,s2
       real(r8) :: sqrts,s15,scl,x1,x2,xacc,htotal2,co2starair

!       Change units from the input of mol/m^3 -> mol/kg:
!       (1 mol/m^3)  x (1 m^3/1024.5 kg)
!       where the ocean''s mean surface density is 1024.5 kg/m^3
!       Note: mol/kg are actually what the body of this routine uses 
!       for calculations.  
!       DU (7/10/2023): changed to use model surface density (rhos)
!       permil = 1.0_r8 / 1024.5_r8
       permil = 1.0_r8 / rhos
       pt=pt_in*permil
       sit=sit_in*permil
       ta=ta_in*permil
       dic=dic_in*permil
       permeg=0.000001_r8
!       To convert input in uatm -> atm
      xco2=xco2_in*permeg
!
! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code. 
! Once calculated, the constants are
! stored and passed in the common block "const". The original version 
! of this code was based on the code by Dickson in Version 2 of 
! "Handbook of Methods for the Analysis of the Various Parameters of 
! the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
!
! Derive simple terms used more than once
!
      tk = 273.15_r8 + t
      tk100 = tk/100.0_r8
      tk1002=tk100*tk100
      invtk=1.0_r8/tk
      dlogtk=log(tk)
      is=19.924_r8*s/(1000.0_r8-1.005_r8*s)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5_r8
      scl=s/1.80655_r8
!
! f = k0(1-pH2O)*correction term for non-ideality
!
! Weiss & Price (1980, Mar. Chem., 8,347-359; Eq 13 with table 6 values)
!
      ff = exp(-162.8301_r8 + 218.2968_r8/tk100  +                     &
     &            90.9241_r8*log(tk100) - 1.47696_r8*tk1002 +          &
     &            s * (0.025695_r8 - 0.025225_r8*tk100 +               &
     &            0.0049867_r8*tk1002))
!
! K0 (Weiss 1974) IS THE CO2 SOLUBILITY IN SEAWATER (IN MMOL M-3 UATM-1)
!
      k0 = exp(93.4517_r8/tk100 - 60.2409_r8 + 23.3585_r8 * log(tk100) &
     & +s*(0.023517_r8 - 0.023656_r8 * tk100 + 0.0047036_r8 * tk1002))

!
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
!
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
!
      k1=10.0_r8**(-1.0_r8*(3670.7_r8*invtk - 62.008_r8 +              &
     & 	      9.7944_r8*dlogtk -0.0118_r8 * s + 0.000116_r8*s2))
!
      k2=10.0_r8**(-1.0_r8*(1394.7_r8*invtk + 4.777_r8 -               &
     &            0.0184_r8*s + 0.000118_r8*s2))
!
! kb = [H][BO2]/[HBO2]
!
! Millero p.669 (1995) using data from Dickson (1990)
!
      kb=exp((-8966.90_r8 - 2890.53_r8*sqrts - 77.942_r8*s +           &
     &            1.728_r8*s15 - 0.0996_r8*s2)*invtk +                 &
     &            (148.0248_r8 + 137.1942_r8*sqrts + 1.62142_r8*s) +   &
     &            (-24.4344_r8 - 25.085_r8*sqrts - 0.2474_r8*s) *      &
     &            dlogtk + 0.053105_r8*sqrts*tk)
!
! k1p = [H][H2PO4]/[H3PO4]
!
! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!
      k1p = exp(-4576.752_r8*invtk + 115.525_r8 - 18.453_r8 * dlogtk + &
     &            (-106.736_r8*invtk + 0.69171_r8) * sqrts +           &
     &            (-0.65643_r8*invtk - 0.01844_r8) * s)      
!
! k2p = [H][HPO4]/[H2PO4]
!
! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!
      k2p = exp(-8814.715_r8*invtk + 172.0883_r8 - 27.927_r8 * dlogtk+ &
     &            (-160.340_r8*invtk + 1.3566_r8) * sqrts +            &
     &            (0.37335_r8*invtk - 0.05778_r8) * s)
!
! k3p = [H][PO4]/[HPO4]
!
! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!
      k3p = exp(-3070.75_r8*invtk - 18.141_r8 +                        &
     &            (17.27039_r8*invtk + 2.81197_r8) *                   &
     &            sqrts + (-44.99486_r8*invtk - 0.09984_r8) * s)
!
! ksi = [H][SiO(OH)3]/[Si(OH)4]
!
! Millero p.671 (1995) using data from Yao and Millero (1995)
!
      ksi = exp(-8904.2_r8*invtk + 117.385_r8 - 19.334_r8 * dlogtk +   &
     &            (-458.79_r8*invtk + 3.5913_r8) * sqrtis +            &
     &            (188.74_r8*invtk - 1.5998_r8) * is +                 &
     &            (-12.1652_r8*invtk + 0.07871_r8) * is2 +             &
     &            log(1.0_r8-0.001005_r8*s))
!
! kw = [H][OH]
!
! Millero p.670 (1995) using composite data
!
      kw = exp(-13847.26_r8*invtk + 148.9652_r8 - 23.6521_r8 *dlogtk+  &
     &            (118.67_r8*invtk - 5.977_r8 + 1.0495_r8 * dlogtk) *  &
     &            sqrts - 0.01615_r8 * s)
!
! ks = [H][SO4]/[HSO4]
!
! Dickson (1990, J. chem. Thermodynamics 22, 113)
!
      ks=exp(-4276.1_r8*invtk + 141.328_r8 - 23.093_r8*dlogtk +        &
     &    (-13856.0_r8*invtk +324.57_r8-47.986_r8*dlogtk)*sqrtis+      &
     &    (35474.0_r8*invtk - 771.54_r8 + 114.723_r8*dlogtk) * is -    &
     &     2698.0_r8*invtk*is**1.5_r8 + 1776.0_r8*invtk*is2 +          &
     &     log(1.0_r8 - 0.001005_r8*s))
!
! kf = [H][F]/[HF]
!
! Dickson and Riley (1979) -- change pH scale to total
!
      kf=exp(1590.2_r8*invtk - 12.641_r8 + 1.525_r8*sqrtis +           &
     &            log(1.0_r8 - 0.001005_r8*s) +                        &
     &            log(1.0_r8 + (0.1400_r8/96.062_r8)*(scl)/ks))
!
! Calculate concentrations for borate, sulfate, and fluoride
!
! Uppstrom (1974)
      bt = 0.000232_r8 * scl/10.811_r8
! Morris & Riley (1966)
      st = 0.14_r8 * scl/96.062_r8
! Riley (1965)
      ft = 0.000067_r8 * scl/18.9984_r8
!
!
! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.
!
! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate htotal. In this case we use the 
! Newton-Raphson "safe" method taken from "Numerical Recipes" (function 
! "rtsafe.f" with error trapping removed).
!
! As currently set, this procedure iterates about 12 times. The x1 and 
! x2 values set below will accomodate ANY oceanographic values. If an 
! initial guess of the pH is known, then the number of iterations can 
! be reduced to about 5 by narrowing the gap between x1 and x2. 
! It is recommended that the first few time steps be run with x1 and x2 
! set as below. After that, set x1 and x2 to the previous value of the 
! pH +/- ~0.5. The current setting of xacc will result in co2star 
! accurate to 3 significant figures (xx.y). Making xacc bigger will 
! result in faster convergence also, but this is not recommended 
! (xacc of 10**-9 drops precision to 2 significant figures).
!
! Parentheses added around negative exponents (Keith Lindsay)
!
      x1 = 10.0_r8**(-phhi)
      x2 = 10.0_r8**(-phlo)
!      xacc = 0.0000000001_r8
      xacc = 0.00000000001_r8
      call drtsafe(x1,x2,xacc,                                         &
    & k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,kw,pt,sit,ksi,ft,kf,ta,ff,  &
    & htotal )
!
! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
!
      htotal2=htotal*htotal
      co2star=dic*htotal2/(htotal2 + k1*htotal + k1*k2)
      co2starair=xco2*ff*atmpres
      dco2star=co2starair-co2star
      ph=-log10(htotal)

        pCO2surf = co2star / ff
        dpCO2    = pCO2surf - xco2*atmpres
!
!  Convert units of output arguments
!      Note: co2star and dco2star are calculated in mol/kg within this  
!      routine. Thus Convert now from mol/kg -> mol/m^3

       co2star  = co2star / permil
       dco2star = dco2star / permil

       pCO2surf = pCO2surf / permeg
       dpCO2    = dpCO2 / permeg
!       write(*,*) '++++++++',pCO2surf,dpCO2,co2star,ff,ph,htotal,  &
!     &k1,k2,permil,permeg      
      RETURN
      END SUBROUTINE co2calc
 
!-----------------------------------------------------------------------     
     SUBROUTINE DRTSAFE(X1,X2,XACC,                                    &
   & k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,kw,pt,sit,ksi,ft,kf,ta,ff,   &
   & DRTSAFE2)
      USE mod_kinds
      implicit none
      real(r8),intent(in)  :: X1,X2,XACC
      real(r8),intent(in)  :: k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt
      real(r8),intent(in)  :: kb,kw,pt,sit,ksi,ft,kf,ta,ff
      real(r8),intent(out) :: DRTSAFE2
      
      integer  :: j,MAXIT
      real(r8) :: FL,DF,FH,XL,XH,SWAP,DXOLD,DX,TEMP,F
!
!      File taken from Numerical Recipes. Modified  R.M.Key 4/94
!
      MAXIT=100
      CALL ta_iter_1(X1,FL,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,   &
    &  	      kw,pt,sit,ksi,ft,kf,ta,ff)
      CALL ta_iter_1(X2,FH,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,   &
    &  	      kw,pt,sit,ksi,ft,kf,ta,ff)
      IF(FL .LT. 0.0_r8) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE2=0.5_r8*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL ta_iter_1(DRTSAFE2,F,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt, &
     & 	      kb,kw,pt,sit,ksi,ft,kf,ta,ff)
      DO J=1,MAXIT
        IF(((DRTSAFE2-XH)*DF-F)*((DRTSAFE2-XL)*DF-F) .GE. 0.0_r8 .OR. &
     &            ABS(2.0_r8*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5_r8*(XH-XL)
          DRTSAFE2=XL+DX
          IF(XL .EQ. DRTSAFE2) RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE2
          DRTSAFE2=DRTSAFE2-DX
          IF(TEMP .EQ. DRTSAFE2) RETURN
      END IF
        IF(ABS(DX) .LT. XACC) RETURN
      CALL ta_iter_1(DRTSAFE2,F,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,  &
     &   	kb,kw,pt,sit,ksi,ft,kf,ta,ff)
        IF(F .LT. 0.0_r8) THEN
          XL=DRTSAFE2
          FL=F
        ELSE
          XH=DRTSAFE2
          FH=F
        END IF
       END DO
      RETURN
      END SUBROUTINE DRTSAFE
!-----------------------------------------------------------------------  
  SUBROUTINE ta_iter_1(x,fn,df,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,&
    &  	      kw,pt,sit,ksi,ft,kf,ta,ff)
      USE mod_kinds
      implicit none
      real(r8),intent(in)  :: x,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb
      real(r8),intent(in)  :: kw,pt,sit,ksi,ft,kf,ta,ff
      real(r8),intent(out) :: fn,df
      
      real(r8) ::  k12,k12p,k123p,x2,x3,c,a,a2,da,b,b2,db
  
!
! This routine expresses TA as a function of DIC, htotal and constants.
! It also calculates the derivative of this function with respect to 
! htotal. It is used in the iterative solution for htotal. In the call
! "x" is the input value for htotal, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhtotal
!
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0_r8 + st/ks
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0_r8*x2 + 2.0_r8*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0_r8*x + k1
!
!   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
!
      fn = k1*x*dic/b +                                   &
     &           2.0_r8*dic*k12/b +                       &
     &           bt/(1.0_r8 + x/kb) +                     &
     &           kw/x +                                   &
     &           pt*k12p*x/a +                            &
     &           2.0_r8*pt*k123p/a +                      &
     &           sit/(1.0_r8 + x/ksi) -                   &
     &           x/c -                                    &
     &           st/(1.0_r8 + ks/x/c) -                   &
     &           ft/(1.0_r8 + kf/x) -                     &
     &           pt*x3/a -                                &
     &           ta
!
!      df = dfn/dx
!
!      df = ((k1*dic*b) - k1*x*dic*db)/b2 -                             &
!     &           2.0_r8*dic*k12*db/b2 -                                &
!     &           bt/kb/(1.0_r8+x/kb)**2.0_r8 -                         &
!     &           kw/x2 +                                               &
!     &           (pt*k12p*(a - x*da))/a2 -                             &
!     &           2.0_r8*pt*k123p*da/a2 -                               &
!     &           sit/ksi/(1.0_r8+x/ksi)**2.0_r8 -                      &
!     &           1.0_r8/c +                                            &
!     &           st*(1.0_r8 + ks/x/c)**(-2.0_r8)*(ks/c/x2) +           &
!     &           ft*(1.0_r8 + kf/x)**(-2.0_r8)*kf/x2 -                 &
!     &           pt*x2*(3.0_r8*a-x*da)/a2

! DU: slightly modified to correct two sign errors (small terms)
      df = ((k1*dic*b) - k1*x*dic*db)/b2 -                             &
     &           2.0_r8*dic*k12*db/b2 -                                &
     &           bt/kb/(1.0_r8+x/kb)**2.0_r8 -                         &
     &           kw/x2 +                                               &
     &           (pt*k12p*(a - x*da))/a2 -                             &
     &           2.0_r8*pt*k123p*da/a2 -                               &
     &           sit/ksi/(1.0_r8+x/ksi)**2.0_r8 -                      &
     &           1.0_r8/c -                                            &
     &           st*(1.0_r8 + ks/x/c)**(-2.0_r8)*(ks/c/x2) -           &
     &           ft*(1.0_r8 + kf/x)**(-2.0_r8)*kf/x2 -                 &
     &           pt*x2*(3.0_r8*a-x*da)/a2
     
      return
      END SUBROUTINE ta_iter_1
!--------------------
#endif    
      

#ifdef OXYGEN
      SUBROUTINE O2_flux (Istr, Iend,                                  &
     &                       LBi, UBi, LBj, UBj,                       &
     &                       IminS, ImaxS, j,                          &
#  ifdef MASKING
     &                       rmask,                                    &
#  endif
     &                       T, S, O2, kw660, ppo, o2sat, O2flx)
!
!***********************************************************************
!                                                                      !
!  Computes the time rate of change of oxygen in the surface           !
!  layer due to air-sea gas exchange in mol/m^3/day                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     O2         Dissolevd oxygen concentration (micromole O2/m^3)     !
!     kw660      gas transfer velocity at a Schmidt number of 660,     !
!                  accounting for sea ice fraction (cm/hr)             !
!     ppo        surface pressure divided by 1 atm.                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     o2sat      dissolved oxygen saturation concentration (mmol/m^3)  !
!                  due to air-sea exchange (mmol/m^2/day)              !
!     o2flx      time rate of oxygen O2 flux in the sea surface        !
!                  due to air-sea exchange (mmol/m^2/day)              !
!                                                                      !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
!
      real(r8),  intent(in) :: ppo
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(LBi:)
      real(r8), intent(in) :: S(LBi:)
      real(r8), intent(in) :: O2(LBi:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(LBi:UBi)
      real(r8), intent(in) :: S(LBi:UBi)
      real(r8), intent(in) :: O2(LBi:UBi)
#  endif

      real(r8), intent(in) :: kw660(LBi:UBi)

      real(r8), intent(out) :: o2sat(LBi:UBi)
      real(r8), intent(out) :: o2flx(LBi:UBi)
!
!  Local variable declarations.
!
      integer :: i

      real(r8) :: sco2, kwo2
      real(r8) :: TT, TK, TS, TS2, TS3, TS4, TS5, CO

      real(r8), parameter :: A0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: A1 = 3.22014_r8       ! saturation
      real(r8), parameter :: A2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: A3 = 4.94457_r8
      real(r8), parameter :: A4 =-0.256847_r8
      real(r8), parameter :: A5 = 3.88767_r8
      real(r8), parameter :: B0 =-0.00624523_r8
      real(r8), parameter :: B1 =-0.00737614_r8
      real(r8), parameter :: B2 =-0.0103410_r8
      real(r8), parameter :: B3 =-0.00817083_r8
      real(r8), parameter :: C0 =-0.000000488682_r8
      real(r8) :: A,B,C,D,E,tsq
! schmidt number polynomial coefficients (freshwater):
      real(r8), parameter :: A0s = 1745.1_r8
      real(r8), parameter :: B0s = -124.34_r8
      real(r8), parameter :: C0s = 4.8055_r8
      real(r8), parameter :: D0s = -0.10115_r8
      real(r8), parameter :: E0s = 0.00086842_r8
! schmidt number polynomial coefficients (35 psu):
      real(r8), parameter :: A35s = 1920.4_r8
      real(r8), parameter :: B35s = -135.6_r8
      real(r8), parameter :: C35s = 5.2122_r8
      real(r8), parameter :: D35s = -0.10939_r8
      real(r8), parameter :: E35s = 0.00093777_r8
!
!=======================================================================
!  Determine coefficients.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
!
! ********************************************************************
!                                     
! Computes the oxygen saturation concentration at 1 atm total pressure
! in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
! in permil). 
!
! FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
! THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
! o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
! 0 permil <= S <= 42 permil
! C
! CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
! o2sato = 282.015 mmol/m^3
!
! ********************************************************************
!
      TT  = 298.15_r8-T(i)
      TK  = 273.15_r8+T(i)
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5             &
     &     + S(i)*(B0 + B1*TS + B2*TS2 + B3*TS3)                       &
     &     + C0*(S(i)*S(i))
      o2sat(i) = EXP(CO)
!
!  Convert from ml/l to mol/m^3
!
      o2sat(i) = (o2sat(i)/22391.6_r8)*1000.0_r8
!
!  Convert from mol/m^3 to mmol/m^3
!
      o2sat(i) = o2sat(i)*1000.0_r8
!
!
!*********************************************************************
!
!  Computes the Schmidt number of oxygen in seawater using the
!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!  Cycles, 12, 141-163).  Input is temperature in deg C.
!  DU (07/2023): changed to use Wanninkhof (2014) 4th order polynomial,
!  using coefficients interpolated to the in situ salinity from the
!  freshwater and 35 psu values given in that paper.
!
!*********************************************************************
!
!      sco2 = 1638.0_r8 - 81.83_r8*t(i) +                               &
!     &       1.483_r8*t(i)**2 - 0.008004_r8*t(i)**3
! DU (7/11/23): changed to use Wanninkhof (2014) relation (for 35 psu)
!       sco2 = 1920.4_r8 - 135.6_r8*t(i) + 5.2122_r8*t(i)**2             &
!      &       - 0.10939_r8*t(i)**3 + 0.00093777_r8*t(i)**4
! interpolate schmidt number polynomial coefficients using surface salinity
! using wanninkhof (2014) polynomials for 35 psu and 0 psu:
      A=A0s+S(i)*(A35s-A0s)/35.0_r8
      B=B0s+S(i)*(B35s-B0s)/35.0_r8
      C=C0s+S(i)*(C35s-C0s)/35.0_r8
      D=D0s+S(i)*(D35s-D0s)/35.0_r8
      E=E0s+S(i)*(E35s-E0s)/35.0_r8
      tsq=T(i)*T(i)
! schmidt number (4th order polynomial)
      sco2= A + B*T(i) + C*tsq + D*tsq*T(i) + E*tsq*tsq
!
!  Compute the transfer velocity for O2 in m/s
!
!    kwo2 = Kw660 * (sco2(t)/660)**-0.5*0.01/3600.0    !(in  m/sec)
!
      kwo2 = Kw660(i) * sqrt(660.0_r8/sco2)   !(in units of cm/hr)
!
!  Compute the transfer velocity for O2 in m/day
!
      KWO2=KWO2*0.01_r8*24.0_r8
!
!  (in units of m/day)
!
!  Compute the saturation concentrations for O2
!
!      o2sat(i) = o2sato(t,s)*ppo       
!  OCMIP
!      o2sat = dosat(t+273.15,s)     
!  Weiss
!
!  Compute time rate of O2 gas exchange
!
      o2flx(i) = kwo2*(o2sat(i)*ppo-o2(i))
#  ifdef MASKING
      ELSE
        o2sat(i)=0.0_r8
        o2flx(i)=0.0_r8
      END IF
#  endif
!        if( (i.eq.535.and.j.eq.690) ) THEN
!            write(*,*)'KWO2=',KWO2,'sco2 =',sco2,'o2sat=',o2sat(i),    &
!     &        'o2flx=',o2flx(i),'T=',T(i),'S=',S(i)
!        end if

      END DO I_LOOP
      RETURN
      END SUBROUTINE O2_flux
#endif

#ifdef OPTICS_OP1
      subroutine kdPAR_compute(Istr, Iend, ng,                          &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, j,                           &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
#  ifdef READ_ZENITH
     &                       zenith_ang,                                &
#  endif
     &                       salt,                                      &
     &                       hzl,                                       &
     &                       z_r,                                       &
     &                       chl1, chl2,                                &
     &                       s1, s2,                                    &
     &                       ddn,                                       &
#  ifdef HAB
     &                       s3, chl3,                                  &
#  endif
     &                       kd)
!
!***********************************************************************
!                                                                      !
!  This routine computes water extinction coefficient kd               !
!                                                                      !
!  Dave Ullman (8/2024), based on code by Peng Xiu                      !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
      USE mod_param
      USE mod_biology

      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: ng
      integer,  intent(in) :: LBi, UBi, LBj, UBj, UBk, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef READ_ZENITH
      real(r8), intent(in) :: zenith_ang(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: salt(LBi:,:)
      real(r8), intent(in) :: hzl(LBi:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: chl1(LBi:,:)
      real(r8), intent(in) :: chl2(LBi:,:)
      real(r8), intent(in) :: s1(LBi:,:)
      real(r8), intent(in) :: s2(LBi:,:)
#   ifdef HAB
      real(r8), intent(in) :: s3(LBi:,:)
      real(r8), intent(in) :: chl3(LBi:,:)
#   endif
      real(r8), intent(in) :: ddn(LBi:,:)
      real(r8), intent(out) :: kd(LBi:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef READ_ZENITH
      real(r8), intent(in) :: zenith_ang(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: salt(LBi:UBi,N(ng))
      real(r8), intent(in) :: hzl(LBi:UBi,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: chl1(LBi:UBi,N(ng))
      real(r8), intent(in) :: chl2(LBi:UBi,N(ng))
      real(r8), intent(in) :: s1(LBi:UBi,N(ng))
      real(r8), intent(in) :: s2(LBi:UBi,N(ng))
#   ifdef HAB
      real(r8), intent(in) :: s3(LBi:UBi,N(ng))
      real(r8), intent(in) :: chl3(LBi:UBi,N(ng))
#   endif
      real(r8), intent(in) :: ddn(LBi:UBi,N(ng))
      real(r8), intent(out) :: kd(LBi:UBi,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, k

      real(r8), dimension(N(ng)) :: bbw
      real(r8), dimension(N(ng)) :: achl1
      real(r8), dimension(N(ng)) :: achl2
      real(r8), dimension(N(ng)) :: ap
      real(r8), dimension(N(ng)) :: adet
      real(r8), dimension(N(ng)) :: acdom
      real(r8), dimension(N(ng)) :: bbp1
      real(r8), dimension(N(ng)) :: bbp2
#   ifdef HAB
      real(r8), dimension(N(ng)) :: achl3
      real(r8), dimension(N(ng)) :: bbp3
#   endif
      real(r8), dimension(N(ng)) :: a_abs
      real(r8), dimension(N(ng)) :: bbp
      real(r8), dimension(N(ng)) :: bb

      real(r8), parameter :: kapa0=-0.057_r8
      real(r8), parameter :: kapa1=0.482_r8
      real(r8), parameter :: kapa2=4.221_r8
      real(r8), parameter :: ksai0=0.183_r8
      real(r8), parameter :: ksai1=0.702_r8
      real(r8), parameter :: ksai2=-2.567_r8
      real(r8), parameter :: alpha0=0.090_r8
      real(r8), parameter :: alpha1=1.465_r8
      real(r8), parameter :: alpha2=-0.667_r8
!thetaa: solar zenith angle radians
!      real(r8), parameter :: thetaa=5.0_r8*3.1416_r8/180.0_r8 ! original
!      real(r8), parameter :: thetaa=41.5_r8*3.1416_r8/180.0_r8)
# if !defined READ_ZENITH
      real(r8), parameter :: zenith_ang=41.5_r8     !constant solar zenith angle
# endif
      real(r8), parameter :: deg2rad=3.141592654_r8/180.0_r8
      real(r8), parameter :: massPOC=12.0_r8       !atomic weight of carbon
      real(r8), parameter :: bbg=0.00035_r8
      real(r8), parameter :: r_phy_POC=0.3_r8      ! ratio of phyto. carbon to total POC
      real(r8), parameter :: mmol2mol=0.001_r8     ! mmole to mol conversion
      real(r8), parameter :: adet_ss=0.011_r8      ! spectral slope for adet
      real(r8), parameter :: acdom_ss=0.020_r8     ! spectral slope for aCDOM
      real(r8), parameter :: a_cdom_sal=4.372_r8   ! constant term for aCDOM(412) from salt
      real(r8), parameter :: b_cdom_sal=-0.129_r8   ! linear term for aCDOM(412) from salt
      real(r8), parameter :: chl_exp=0.65_r8       ! exponent for chl absorption
      real(r8), parameter :: bbrat_sw=0.5_r8       ! backscattering ratio for seawater
      real(r8), parameter :: a_bbw_salt=1.0_r8     ! constant term for bbw from salt
      real(r8), parameter :: b_bbw_salt=0.00811_r8 ! linear term for bbw from salt

! values at 490 nm:
      real(r8), parameter :: aw_abs=0.0150_r8
      real(r8), parameter :: bw_abs=0.0031_r8
      real(r8), parameter :: achlstar=0.7500_r8
      real(r8), parameter :: adetstar=0.1000_r8

      real(r8) :: cff, adet_fac1, adet_fac2, acdom_fac
      real(r8) :: bp1_fac1,bp1_exp,bp2_fac1,bp2_exp,bp_fac
      real(r8) :: aa_int,bb_int,z_int,aa_avg,bb_avg
      real(r8) :: kpar1,kpar2,thetaa,acdtmp

! factors for detrital absorption calculation
      adet_fac1=mmol2mol*massPOC*c2n(ng)
      adet_fac2=exp(-adet_ss*(490-440))     ! exponential correction from 440 to 490 nm
! factors for CDOM absorption calculation
      acdom_fac=exp(-acdom_ss*(490-412))    ! exponential correction from 412 to 490 nm
! factors for backscatter due to small phyto. (s1)
      bp1_fac1=c2n(ng)*massPOC/(r_phy_POC*476935.8_r8)
      bp1_exp=1.0_r8/1.277_r8;
! factors for backscatter due to diatoms (s2)
      bp2_fac1=c2n(ng)*massPOC/(r_phy_POC*17069.0_r8)
      bp2_exp=1.0_r8/0.859_r8;
      bp_fac=1/sqrt(490.0_r8/512.0_r8)        ! factor to convert from 512 nm to 490

      I_LOOP: DO i=Istr,Iend
        
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
#  ifdef READ_ZENITH
          thetaa=deg2rad*zenith_ang(i,j)
#  else
          thetaa=deg2rad*zenith_ang
#  endif
          do k=1,N(ng)

          ! absorption due to chlorophyll
             achl1(k)=0.06_r8*achlstar*(chl1(i,k))**chl_exp
             achl2(k)=0.06_r8*achlstar*(chl2(i,k))**chl_exp
#   ifdef HAB
             achl3(k)=0.06_r8*achlstar*(chl3(i,k))**chl_exp
             ap(k)=achl1(k)+achl2(k)+achl3(k)
#   else
             ap(k)=achl1(k)+achl2(k)
#   endif

          ! absorption due to detritus
             adet(k)=adetstar*ddn(i,k)*adet_fac1*adet_fac2
          ! absorption due to CDOM
!             aCDOM(k)=(a_cdom_sal+b_cdom_sal*salt(i,k))*acdom_fac
          ! make sure aCDOM does not go negative (at very high salinity)
             acdtmp=(a_cdom_sal+b_cdom_sal*salt(i,k))*acdom_fac
             aCDOM(k)=MAX(acdtmp,0.0_r8)
          ! total absorption coefficient
             a_abs(k)=ap(k)+adet(k)+aCDOM(k)+aw_abs

          ! backscatter due to small phytoplankton
             bbp1(k)=( (s1(i,k)*bp1_fac1)**bp1_exp )*bp_fac
          ! backscatter due to diatoms
             bbp2(k)=( (s2(i,k)*bp2_fac1)**bp2_exp )*bp_fac
#   ifdef HAB
          ! backscatter due to HAB phytoplankton (P-n diatoms)
             bbp3(k)=( (s3(i,k)*bp2_fac1)**bp2_exp )*bp_fac
          ! total backscatter due to phytoplankton
             bbp(k)=bbp1(k)+bbp2(k)+bbp3(k)+bbg
#   else
             bbp(k)=bbp1(k)+bbp2(k)+bbg
#   endif
          ! backscatter due to seawater
             bbw(k)=bbrat_sw*bw_abs*( a_bbw_salt + b_bbw_salt*salt(i,k) )
          ! total backscattering coefficient
             bb(k)=bbw(k)+bbp(k)
          enddo   ! end k loop

         !Kpar method from Lee et al., (2005)
         ! use average a_abs and bb over portion of the water column above the kth level
         aa_int=0.0_r8
         bb_int=0.0_r8
         z_int=0.0_r8
         do k=N(ng),1,-1
           aa_int=aa_int+a_abs(k)*hzl(i,k)
           bb_int=bb_int+bb(k)*hzl(i,k)
           z_int=z_int+hzl(i,k)
           aa_avg=aa_int/z_int
           bb_avg=bb_int/z_int
 
           kpar1=(kapa0+kapa1*sqrt(aa_avg)+kapa2*bb_avg)*(1.0_r8+alpha0*sin(thetaa))
           kpar2=(ksai0+ksai1*aa_avg+ksai2*bb_avg)*(alpha1+alpha2*(cos(thetaa)))

           kd(i,k)= kpar1+kpar2/sqrt(1.0_r8+abs(z_r(i,j,k)))
         enddo

#  ifdef MASKING
       ELSE
         do k=1,N(ng)
           kd(i,k)=1000000.00_r8
         enddo
       END IF
#  endif
      END DO I_LOOP

      RETURN
      END SUBROUTINE kdPAR_compute
#endif







