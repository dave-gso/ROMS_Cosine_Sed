#include "cppdefs.h"
      MODULE sediment_biology_mod

#ifdef SEDBIO
      implicit none
      PRIVATE
      PUBLIC :: sediment_bio

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment_bio(ng, Istr, Iend, LBi, UBi,                &
     &                      dtdays,                                   &
#   ifdef MASKING
     &                      rmask,                                     &
#   endif
     &                      temp, DOX, dNO3, dNH4, dPO4, dSiO4,        &
     &                      FPOC, FPON, FPOP, FPSi,                    &
     &                      bPOC, bPON, bPOP, bPSi,                    &
     &                      bNO3, bNH4, bPO4, bSi,                     &
     &                      bUC, bUN, bUP, bUS,                        &
     &                      JNO3, JNH4, JPO4, JSi, SOD)
!***********************************************************************
!
      USE mod_kinds
      USE mod_biology, ONLY : bUmax, bdep, balpha, bw, btheta_diag,     &
     &                       bnit, btheta_nit, bdo_nit, bdenit,        &
     &                       btheta_denit, bpsi_n, bdo_c, bao2, bpi,   &
     &                       bpsi_p, bfc, bfn, bfp, bfs
      USE mod_sedbio, ONLY : nspc

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Istr, Iend
      integer, intent(in) :: LBi, UBi
      real(r8), intent(in) :: dtdays
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:)
#  endif
      real(r8), intent(in) :: temp(LBi:)
      real(r8), intent(in) :: DOX(LBi:)
      real(r8), intent(in) :: dNO3(LBi:)
      real(r8), intent(in) :: dNH4(LBi:)
      real(r8), intent(in) :: dPO4(LBi:)
      real(r8), intent(in) :: dSiO4(LBi:)
      real(r8), intent(in) :: FPOC(LBi:)
      real(r8), intent(in) :: FPON(LBi:)
      real(r8), intent(in) :: FPOP(LBi:)
      real(r8), intent(in) :: FPSi(LBi:)
      real(r8), intent(inout) :: bPOC(LBi:,:)
      real(r8), intent(inout) :: bPON(LBi:,:)
      real(r8), intent(inout) :: bPOP(LBi:,:)
      real(r8), intent(inout) :: bPSi(LBi:,:)
      real(r8), intent(inout) :: bNO3(LBi:)
      real(r8), intent(inout) :: bNH4(LBi:)
      real(r8), intent(inout) :: bPO4(LBi:)
      real(r8), intent(inout) :: bSi(LBi:)
      real(r8), intent(inout) :: bUC(LBi:,:)
      real(r8), intent(inout) :: bUN(LBi:,:)
      real(r8), intent(inout) :: bUP(LBi:,:)
      real(r8), intent(inout) :: bUS(LBi:,:)
      real(r8), intent(inout) :: JNO3(LBi:)
      real(r8), intent(inout) :: JNH4(LBi:)
      real(r8), intent(inout) :: JPO4(LBi:)
      real(r8), intent(inout) :: JSi(LBi:)
      real(r8), intent(inout) :: SOD(LBi:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi)
#  endif
      real(r8), intent(in) :: temp(LBi:UBi)
      real(r8), intent(in) :: DOX(LBi:UBi)
      real(r8), intent(in) :: dNO3(LBi:UBi)
      real(r8), intent(in) :: dNH4(LBi:UBi)
      real(r8), intent(in) :: dPO4(LBi:UBi)
      real(r8), intent(in) :: dSiO4(LBi:UBi)
      real(r8), intent(in) :: FPOC(LBi:UBi)
      real(r8), intent(in) :: FPON(LBi:UBi)
      real(r8), intent(in) :: FPOP(LBi:UBi)
      real(r8), intent(in) :: FPSi(LBi:UBi)
      real(r8), intent(inout) :: bPOC(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPON(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPOP(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPSi(LBi:UBi,nspc)
      real(r8), intent(inout) :: bNO3(LBi:UBi)
      real(r8), intent(inout) :: bNH4(LBi:UBi)
      real(r8), intent(inout) :: bPO4(LBi:UBi)
      real(r8), intent(inout) :: bSi(LBi:UBi)
      real(r8), intent(inout) :: bUC(LBi:UBi,nspc)
      real(r8), intent(inout) :: bUN(LBi:UBi,nspc)
      real(r8), intent(inout) :: bUP(LBi:UBi,nspc)
      real(r8), intent(inout) :: bUS(LBi:UBi,nspc)
      real(r8), intent(inout) :: JNO3(LBi:UBi)
      real(r8), intent(inout) :: JNH4(LBi:UBi)
      real(r8), intent(inout) :: JPO4(LBi:UBi)
      real(r8), intent(inout) :: JSi(LBi:UBi)
      real(r8), intent(inout) :: SOD(LBi:UBi)
# endif

! local variables:

      integer :: i, k
      real(r8) :: Qt
      real(r8), dimension(nspc) :: JPOC
      real(r8), dimension(nspc) :: JPON
      real(r8), dimension(nspc) :: JPOP
      real(r8), dimension(nspc) :: JPSi
      real(r8) :: stc, CSOD, NSOD, Jnit, Jdenit, JCO2, O2, DOxyg

      real(r8) :: w2, cff, cff1, cff2, cff3
      real(r8), parameter :: TR = 20.0_r8
      real(r8), parameter :: eps = 1.0e-10_r8
      real(r8), parameter :: o2conv = 32.0e-03_r8

      w2=bw(ng)*2.74e-5_r8   ! conversion from cm/yr to m/day
!
!-----------------------------------------------------------------------
! Compute diagenesis flux (mmol/m3/day) and update POM concentration
!-----------------------------------------------------------------------
!
        cff=dtdays/bdep(ng)
        DO i=Istr,Iend
# ifdef MASKING
          IF( rmask(i).gt.0) THEN
# endif
! if near-bottom DOX is less than zero, make it zero:
            DOxyg=max(DOX(i),0.0_r8)
!            O2=DOX(i)*o2conv !mmol.m-3 to mg/L
            O2=DOxyg*o2conv !mmol.m-3 to mg/L
            Qt=btheta_diag(ng)**(temp(i)-TR)
            DO k=1,nspc
!        Diagenesis flux (mmol/m3/day)
                JPOC(k)=bUC(i,k)*bPOC(i,k)*Qt
                JPON(k)=bUN(i,k)*bPON(i,k)*Qt
                JPOP(k)=bUP(i,k)*bPOP(i,k)*Qt
                JPSi(k)=bUS(i,k)*bPSi(i,k)*Qt
!        POM concentration (mmol/m3)
                bPOC(i,k)=bPOC(i,k)+FPOC(i)*bfc(ng,k)*cff            &
     &                      -JPOC(k)*dtdays-w2*bPOC(i,k)*cff
                bPON(i,k)=bPON(i,k)+FPON(i)*bfn(ng,k)*cff            &
     &                      -JPON(k)*dtdays-w2*bPON(i,k)*cff
                bPOP(i,k)=bPOP(i,k)+FPOP(i)*bfp(ng,k)*cff            &
     &                      -JPOP(k)*dtdays-w2*bPOP(i,k)*cff
                bPSi(i,k)=bPSi(i,k)+FPSi(i)*bfs(ng,k)*cff            &
     &                      -JPSi(k)*dtdays-w2*bPSi(i,k)*cff
            END DO
!        Nutrient fluxes and SOD (mmol/m2/day)
            Jnit=bNH4(i)*bnit(ng)*(O2/(O2+bdo_nit(ng)))*               &
     &           btheta_nit(ng)**(temp(i)-TR)
            CSOD=0.0_r8
            DO k=1,nspc
                CSOD=CSOD+JPOC(k)
            END DO
            JCO2=CSOD   ! computing CO2 flux, but at present not using
            CSOD=CSOD*(O2/(O2+bdo_c(ng)))*bdep(ng)
            NSOD=bao2(ng)*Jnit*bdep(ng)
            if( CSOD.lt.0.0_r8.or.NSOD.lt.0.0_r8) then
                write(*,*)'i=',i,'CSOD=',CSOD,'NSOD=',NSOD
                write(*,*)'DOX(i)=',DOX(i),'JPOC(1:3)=',JPOC(1:3)
            end if
            SOD(i)=CSOD+NSOD
!        Nutrient fluxes
!            stc=SOD(i)/max(DOX(i),1.0e-6_r8)
            stc=SOD(i)/max(DOxyg,1.0e-6_r8)
            JNO3(i)=stc*(bNO3(i)-dNO3(i))
            JNH4(i)=stc*(bNH4(i)-dNH4(i))
            cff=bpi(ng)*exp(-O2/bpsi_p(ng))
            JPO4(i)=stc*(bPO4(i)*cff-dPO4(i))
            JSi(i)=stc*(bSi(i)-dSiO4(i))
!        Update pore water nutrient concentrations
            cff=dtdays/bdep(ng)
            bNO3(i)=bNO3(i)+Jnit*dtdays-JNO3(i)*cff
            IF (bNO3(i).gt.0.0_r8) THEN
                Jdenit=min(bdenit(ng)*btheta_denit(ng)**(temp(i)-TR)*    &
     &                  exp(-O2/bpsi_n(ng))*dNO3(i),0.8_r8*CSOD/bdep(ng))
                bNO3(i)=bNO3(i)-Jdenit*dtdays
            END IF
 !           if( (i.eq.531.and.j.eq.693) ) then
 !               write(*,*)'sediment_bio: i=',i,'j=',j,'SOD=',SOD(i)
 !               write(*,*)'sediment_bio: i=',i,'j=',j,'BPSi(i,1)=',BPSi(i,1)
 !               write(*,*)'                            BPSi(i,2)=',BPSi(i,2)
 !           end if
            cff1=0.0_r8
            cff2=0.0_r8
            cff3=0.0_r8
            DO k=1,nspc
                cff1=cff1+JPON(k)
                cff2=cff2+JPOP(k)
                cff3=cff3+JPSi(k)
            END DO
            bNH4(i)=bNH4(i)+cff1*dtdays-Jnit*dtdays-JNH4(i)*cff
            bPO4(i)=bPO4(i)+cff2*dtdays-JPO4(i)*cff
            bSi(i)=bSi(i)+cff3*dtdays-JSi(i)*cff
!         ensure that concentrations are >= 0
            bNO3(i)=max(bNO3(i),0.0_r8)
            bNH4(i)=max(bNH4(i),0.0_r8)
            bPO4(i)=max(bPO4(i),0.0_r8)
            bSi(i)=max(bSi(i),0.0_r8)
            DO k=1,nspc
                bPOC(i,k)=max(bPOC(i,k),0.0_r8)
                bPON(i,k)=max(bPON(i,k),0.0_r8)
                bPOP(i,k)=max(bPOP(i,k),0.0_r8)
                bPSi(i,k)=max(bPSi(i,k),0.0_r8)
            END DO
!         compute decay rates
            DO k=1,nspc
                bUC(i,k)=bUC(i,k)+balpha(ng)*bUmax(ng,k)*bUC(i,k)*dtdays &
     &              -FPOC(i)*bfc(ng,k)*dtdays*bUC(i,k)/bPOC(i,k)/bdep(ng)
                bUN(i,k)=bUN(i,k)+balpha(ng)*bUmax(ng,k)*bUN(i,k)*dtdays &
     &              -FPON(i)*bfn(ng,k)*dtdays*bUN(i,k)/bPON(i,k)/bdep(ng)
                bUP(i,k)=bUP(i,k)+balpha(ng)*bUmax(ng,k)*bUP(i,k)*dtdays &
     &              -FPOP(i)*bfp(ng,k)*dtdays*bUP(i,k)/bPOP(i,k)/bdep(ng)
!                bUS(i,k)=bUS(i,k)+balpha(ng)*bUmax(ng,k)*bUS(i,k)*dtdays &
!     &              -FPSi(i)*bfs(ng,k)*dtdays*bUS(i,k)/bPSi(i,k)/bdep(ng)
!         for silicate, decay rates are not bacterially-mediated, so
!           set to max values
                bUS(i,k)=bUmax(ng,k)
!         apply limits to decay rates and ensure positivity
                bUC(i,k)=max( min(bUC(i,k),bUmax(ng,k)),0.0_r8 )
                bUN(i,k)=max( min(bUN(i,k),bUmax(ng,k)),0.0_r8 )
                bUP(i,k)=max( min(bUP(i,k),bUmax(ng,k)),0.0_r8 )
!                bUS(i,k)=max( min(bUS(i,k),bUmax(ng,k)),0.0_r8 )
!         if particulate concentration < eps, then set decay rate to max value
                if( bPOC(i,k).lt.eps) bUC(i,k)=bUmax(ng,k)
                if( bPON(i,k).lt.eps) bUN(i,k)=bUmax(ng,k)
                if( bPOP(i,k).lt.eps) bUP(i,k)=bUmax(ng,k)
!                if( bPSi(i,k).lt.eps) bUS(i,k)=bUmax(ng,k)
            END DO
!# ifdef MASKING
!            bNO3(i)=bNO3(i)*rmask(i)
!            bNH4(i)=bNH4(i)*rmask(i)
!            bPO4(i)=bPO4(i)*rmask(i)
!            bSi(i)=bSi(i)*rmask(i)
!            JNO3(i)=JNO3(i)*rmask(i)
!            JNH4(i)=JNH4(i)*rmask(i)
!            JPO4(i)=JPO4(i)*rmask(i)
!            JSi(i)=JSi(i)*rmask(i)
!            SOD(i)=SOD(i)*rmask(i)
!           DO k=1,nspc
!                bPOC(i,k)=bPOC(i,k)*rmask(i)
!                bPON(i,k)=bPON(i,k)*rmask(i)
!                bPOP(i,k)=bPOP(i,k)*rmask(i)
!                bPSi(i,k)=bPSi(i,k)*rmask(i)
!                bUC(i,k)=bUC(i,k)*rmask(i)
!                bUN(i,k)=bUN(i,k)*rmask(i)
!                bUP(i,k)=bUP(i,k)*rmask(i)
!                bUS(i,k)=bUS(i,k)*rmask(i)
!            END DO
!# endif
!            if(JNO3(i).gt.1.0e10_r8) then
!                write(*,*)'i=',i,'JNO3=',JNO3(i)
!            end if
# ifdef MASKING
          END IF
# endif
        END DO


      RETURN
      END SUBROUTINE sediment_bio
#endif
      END MODULE sediment_biology_mod
