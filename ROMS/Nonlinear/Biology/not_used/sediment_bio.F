#include "cppdefs.h"

      MODULE sediment_bio_mod

#if defined NONLINEAR && defined BIOLOGY && defined SEDBIO
!
!
      implicit none

      PRIVATE
      PUBLIC  :: sediment_bio

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment_bio (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbio
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
! need to get FPOC, FPON, etc into this subroutine!!!
      CALL sediment_bio_tile (ng, tile,                                 &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
# ifdef WET_DRY
     &                   GRID(ng) % rmask_wet,                          &
# endif
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
     &                   OCEAN(ng) % t,                                 &
     &                   SEDBIO(ng) % bPOC,                             &
     &                   SEDBIO(ng) % bPON,                             &
     &                   SEDBIO(ng) % bPOP,                             &
     &                   SEDBIO(ng) % bPOS,                             &
     &                   SEDBIO(ng) % bNO3,                             &
     &                   SEDBIO(ng) % bNH4,                             &
     &                   SEDBIO(ng) % bPO4,                             &
     &                   SEDBIO(ng) % bSi,                              &
     &                   SEDBIO(ng) % bUC,                              &
     &                   SEDBIO(ng) % bUN,                              &
     &                   SEDBIO(ng) % bUP,                              &
     &                   SEDBIO(ng) % bUS,                              &
     &                   SEDBIO(ng) % JNO3,                             &
     &                   SEDBIO(ng) % JNH4,                             &
     &                   SEDBIO(ng) % JPO4,                             &
     &                   SEDBIO(ng) % JSi,                              &
     &                   SEDBIO(ng) % SOD)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sediment_bio
!
!***********************************************************************
      SUBROUTINE sediment_bio_tile (ng, tile,                           &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
# ifdef WET_DRY
     &                         rmask_wet,                               &
# endif
     &                         bustr, bvstr,                            &
     &                         t,                                       &
     &                         bPOC, bPON, bPOP, bPSi,                  &
     &                         bNO3, bNH4, bPO4, bSi,                   &
     &                         bUC, bUN, bUP, bUS,                      &
     &                         JNO3, JNH4, JPO4, JSi, SOD)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sedbio
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bPOC(LBi:,LBj:,:)
      real(r8), intent(inout) :: bPON(LBi:,LBj:,:)
      real(r8), intent(inout) :: bPOP(LBi:,LBj:,:)
      real(r8), intent(inout) :: bPSi(LBi:,LBj:,:)
      real(r8), intent(inout) :: bNO3(LBi:,LBj:)
      real(r8), intent(inout) :: bNH4(LBi:,LBj:)
      real(r8), intent(inout) :: bPO4(LBi:,LBj:)
      real(r8), intent(inout) :: bSi(LBi:,LBj:)
      real(r8), intent(inout) :: bUC(LBi:,LBj:,:)
      real(r8), intent(inout) :: bUN(LBi:,LBj:,:)
      real(r8), intent(inout) :: bUP(LBi:,LBj:,:)
      real(r8), intent(inout) :: bUS(LBi:,LBj:,:)
      real(r8), intent(inout) :: JNO3(LBi:,LBj:)
      real(r8), intent(inout) :: JNH4(LBi:,LBj:)
      real(r8), intent(inout) :: JPO4(LBi:,LBj:)
      real(r8), intent(inout) :: JSi(LBi:,LBj:)
      real(r8), intent(inout) :: SOD(LBi:,LBj:)
# else
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: bPOC(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bPON(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bPOP(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bPSi(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bNO3(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bNH4(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bPO4(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bSi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bUC(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bUN(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bUP(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: bUS(LBi:UBi,LBj:UBj,nspc)
      real(r8), intent(inout) :: JNO3(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: JNH4(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: JPO4(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: JSi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: SOD(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations. (NEED TO MODIFY or delete)
!
      integer :: Ksed, i, ised, j, k, ks
      integer :: bnew

! new ones:
! near-bed water column dissolved concentrations
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dOX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dNO3
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dNH4
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dPO4
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dSiO4
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: temp

      integer :: i, j, k
      real(r8) :: Qt
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,nspc) :: JPOC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,nspc) :: JPON
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,nspc) :: JPOP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,nspc) :: JPSi
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: O2
      real(r8) :: stc, CSOD, NSOD, Jnit, Jdenit, JCO2

      real(r8) :: w2, cff, cff1, cff2, cff3
      real(r8), parameter :: TR = 20.0_r8
      real(r8), parameter :: eps = 1.0e-10_r8

# include "set_bounds.h"

      w2=bw*2.74e-5_r8   ! conversion from cm/yr to m/day
! copy near-bottom water column concentrations to local arrays
      DO j=Jstr,Jend
        DO i=Istr,Iend
            dOX(i,j)=t(i,j,1,nstp,iOxyg)
            O2(i,j)=dOX(i,j)*32e-3_r8 ! conversion from mmol/m3 to mg/l
            dNO3(i,j)=t(i,j,1,nstp,iNO3_)
            dNH4(i,j)=t(i,j,1,nstp,iNH4_)
            dPO4(i,j)=t(i,j,1,nstp,iPO4_)
            dSiO4(i,j)=t(i,j,1,nstp,iSiOH)
            temp(i,j)=t(i,j,1,nstp,itemp)
#   ifdef WET_DRY
            dOX(i,j)=dOX(i,j)*rmask_wet(i,j)
            O2(i,j)=O2(i,j)*rmask_wet(i,j)
            dNO3(i,j)=dNO3(i,j)*rmask_wet(i,j)
            dNH4(i,j)=dNH4(i,j)*rmask_wet(i,j)
            dPO4(i,j)=dPO4(i,j)*rmask_wet(i,j)
            dSiO4(i,j)=dSiO4(i,j)*rmask_wet(i,j)
            temp(i,j)=temp(i,j)*rmask_wet(i,j)
#   endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Compute diagenesis flux (mmol/m3/day) and update POM concentration
!-----------------------------------------------------------------------
!
      cff=dtdays/bdep
      DO j=Jstr,Jend
        DO i=Istr,Iend
            Qt=btheta_diag**(temp(i,j)-TR)
            DO k=1,nspc
!        Diagenesis flux (mmol/m3/day)
                JPOC(i,j,k)=bUC(i,j,k)*bPOC(i,j,k)*Qt
                JPON(i,j,k)=bUN(i,j,k)*bPON(i,j,k)*Qt
                JPOP(i,j,k)=bUP(i,j,k)*bPOP(i,j,k)*Qt
                JPSi(i,j,k)=bUS(i,j,k)*bPSi(i,j,k)*Qt
!        POM concentration (mmol/m3)
                bPOC(i,j,k)=bPOC(i,j,k)+FPOC(i,j)*bfc(k)*cff            &
     &                      -JPOC(i,j,k)*dtdays-w2*bPOC(i,j,k)*cff
                bPON(i,j,k)=bPON(i,j,k)+FPON(i,j)*bfn(k)*cff            &
     &                      -JPON(i,j,k)*dtdays-w2*bPON(i,j,k)*cff
                bPOP(i,j,k)=bPOP(i,j,k)+FPOP(i,j)*bfp(k)*cff            &
     &                      -JPOP(i,j,k)*dtdays-w2*bPOP(i,j,k)*cff
                bPSi(i,j,k)=bPSi(i,j,k)+FPSi(i,j)*bfs(k)*cff            &
     &                      -JPSi(i,j,k)*dtdays-w2*bPSi(i,j,k)*cff
            END DO
!        Nutrient fluxes and SOD (mmol/m2/day)
            Jnit=bNH4(i,j)*bnit*(O2(i,j)/(O2(i,j)+bdo_nit))*            &
     &           btheta_nit**(temp(i,j)-TR)
            CSOD=0.0_r8
            DO k=1,nspc
                CSOD=CSOD+JPOC(i,j,k)
            END DO
            JCO2=CSOD   ! computing CO2 flux, but at present not using
            CSOD=CSOD*(O2(i,j)/(O2(i,j)+bdo_c))*bdep
            NSOD=bao2*Jnit*bdep
            SOD(i,j)=CSOD+NSOD
!        Nutrient fluxes
            stc=SOD(i,j)/max(DOX(i,j),1.0e-6_r8)
            JNO3(i,j)=stc*(bNO3(i,j)-dNO3(i,j))
            JNH4(i,j)=stc*(bNH4(i,j)-dNH4(i,j))
            cff=bpi*exp(-O2(i,j)/bpsi_p)
            JPO4(i,j)=stc*(bPO4(i,j)*cff-dPO4(i,j))
            JSi(i,j)=stc*(bSi(i,j)-dSiO4(i,j))
!        Update pore water nutrient concentrations
            cff=dtdays/bdep
            bNO3(i,j)=bNO3(i,j)+Jnit*dtdays-JNO3(i,j)*cff
            IF (bNO3(i,j).gt.0.0_r8) THEN
                Jdenit=min(bdenit*btheta_denit**(temp(i,j)-TR)*         &
     &                  exp(-O2(i,j)/bpsi_n)*dNO3(i,j),0.8_r8*CSOD/bdep)
                bNO3(i,j)=bNO3(i,j)-Jdenit*dtdays
            END IF
            cff1=0.0_r8
            cff2=0.0_r8
            cff3=0.0_r8
            DO k=1,nspc
                cff1=cff1+JPON(i,j,k)
                cff2=cff2+JPOP(i,j,k)
                cff3=cff3+JPSi(i,j,k)
            END DO
            bNH4(i,j)=bNH4(i,j)+cff1*dtdays-Jnit*dtdays-JNH4(i,j)*cff
            bPO4(i,j)=bPO4(i,j)+cff2*dtdays-JPO4(i,j)*cff
            bSi(i,j)=bSi(i,j)+cff3*dtdays-JSi(i,j)*cff
!         ensure that concentrations are >= 0
            bNO3(i,j)=max(bNO3(i,j),0.0_r8)
            bNH4(i,j)=max(bNH4(i,j),0.0_r8)
            bPO4(i,j)=max(bPO4(i,j),0.0_r8)
            bSi(i,j)=max(bSi(i,j),0.0_r8)
            DO k=1,nspc
                bPOC(i,j,k)=max(bPOC(i,j,k),0.0_r8)
                bPON(i,j,k)=max(bPON(i,j,k),0.0_r8)
                bPOP(i,j,k)=max(bPOP(i,j,k),0.0_r8)
                bPSi(i,j,k)=max(bPSi(i,j,k),0.0_r8)
            END DO
!         compute decay rates
            DO k=1,nspc
                bUC(i,j,k)=bUC(i,j,k)+balpha*bUmax(k)*bUC(i,j,k)*dtdays &
     &              -FPOC(i,j)*bfc(k)*dtdays*bUC(i,j,k)/bPOC(i,j,k)/bdep
                bUN(i,j,k)=bUN(i,j,k)+balpha*bUmax(k)*bUN(i,j,k)*dtdays &
     &              -FPON(i,j)*bfn(k)*dtdays*bUN(i,j,k)/bPON(i,j,k)/bdep
                bUP(i,j,k)=bUP(i,j,k)+balpha*bUmax(k)*bUP(i,j,k)*dtdays &
     &              -FPOP(i,j)*bfp(k)*dtdays*bUP(i,j,k)/bPOP(i,j,k)/bdep
                bUS(i,j,k)=bUS(i,j,k)+balpha*bUmax(k)*bUS(i,j,k)*dtdays &
     &              -FPSi(i,j)*bfs(k)*dtdays*bUS(i,j,k)/bPSi(i,j,k)/bdep
!         apply limits to decay rates and ensure positivity
                bUC(i,j,k)=max( min(bUC(i,j,k),bUmax(k)),0.0_r8 )
                bUN(i,j,k)=max( min(bUN(i,j,k),bUmax(k)),0.0_r8 )
                bUP(i,j,k)=max( min(bUP(i,j,k),bUmax(k)),0.0_r8 )
                bUS(i,j,k)=max( min(bUS(i,j,k),bUmax(k)),0.0_r8 )
!         if particulate concentration < eps, then set decay rate to max value
                if( bPOC(i,j,k).lt.eps) bUC(i,j,k)=bUmax(k)
                if( bPON(i,j,k).lt.eps) bUN(i,j,k)=bUmax(k)
                if( bPOP(i,j,k).lt.eps) bUP(i,j,k)=bUmax(k)
                if( bPSi(i,j,k).lt.eps) bUS(i,j,k)=bUmax(k)
            END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE sediment_bio_tile
#endif
      END MODULE sediment_bio_mod
