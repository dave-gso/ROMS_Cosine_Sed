
!
      USE mod_kinds
!
      implicit none
      integer, parameter :: nspc = 3   ! number of sediment particulate classes
      integer, parameter :: NPOM = 4   ! number of POM variables
                                       ! (ibPOC, ibPON, ibPOP, ibPSi)
      integer, parameter :: NPWC = 6   ! number of pore water conc. variables
                                       ! (ibNO3, ibNH4, ibPO4, ibSi, ibTIC, ibAlk)
      integer, parameter :: NDR = 4    ! number of decay rate variables
                                       ! (ibUC, ibUN, ibUP, ibUS)
      integer, parameter :: NSF = 7    ! number of sediment-water flux variables
                                       ! (ibJNO3, ibJNH4, ibJPO4, ibJSi, ibSOD, ibJTIC, ibJALK)
      integer :: NBBT                  ! total number of sediment bio variables

! benthic variables indices
      integer :: idbPM(NPOM)      ! Benthic particulate matter tracers (indices for i/o?)
      integer :: ibbPM(NPOM)      ! indices for SEDBIOL arrays
      integer, parameter  :: ibPOC = 1        ! benthic POC concentration
      integer, parameter  :: ibPON = 2        ! benthic PON concentration
      integer, parameter  :: ibPOP = 3        ! benthic POP concentration
      integer, parameter  :: ibPSi = 4        ! benthic Silicate concentration
      integer :: idbPW(NPWC)      ! Benthic pore water concentration tracers
      integer :: ibbPW(NPWC)      ! indices for SEDBIOL arrays
      integer, parameter  :: ibNO3 = 1        ! pore water NO3 concentration
      integer, parameter  :: ibNH4 = 2        ! pore water NH4 concentration
      integer, parameter  :: ibPO4 = 3        ! pore water PO4 concentration
      integer, parameter  :: ibSi  = 4        ! pore water Silicate concentration
      integer, parameter  :: ibTIC = 5        ! pore water TIC concentration
      integer, parameter  :: ibAlk = 6        ! pore water alkalinity concentration
      integer :: idbDR(NDR)       ! Benthic POM decay rates
      integer :: ibbDR(NDR)       ! indices for SEDBIOL arrays
      integer, parameter   :: ibUC = 1        ! decay rate of particulate carbon
      integer, parameter   :: ibUN = 2        ! decay rate of particulate nitrogen
      integer, parameter   :: ibUP = 3        ! decay rate of particulate Phosphorus
      integer, parameter   :: ibUS = 4        ! decay rate of particulate Silicate
      integer :: idbSF(NSF)       ! Benthic fluxes
      integer :: ibbSF(NSF)       ! indices for SEDBIOL arrays
      integer, parameter  :: ibJNO3 = 1       ! benthic flux of NO3
      integer, parameter  :: ibJNH4 = 2       ! benthic flux of NH4
      integer, parameter  :: ibJPO4 = 3       ! benthic flux of PO4
      integer, parameter  :: ibJSi  = 4       ! benthic flux of Silicate
      integer, parameter  :: ibSOD  = 5       ! benthic oxygen demand
      integer, parameter  :: ibJTIC = 6       ! benthic flux of TIC
      integer, parameter  :: ibJAlk = 7       ! benthic flux of alkalinity
    
      TYPE T_SEDBIOL
!
! sediment model variables
!
!   sedPOM(:,:,:,ibPOC) => Particulate organic carbon concentration
!   sedPOM(:,:,:,ibPON) => Particulate organic nitrogen concentration
!   sedPOM(:,:,:,ibPOP) => Particulate organic phosphorous concentration
!   sedPOM(:,:,:,ibPSi) => Particulate silicon concentration

!   sedPoreWaterCon(:,:,ibNO3) => pore water nitrate concentration
!   sedPoreWaterCon(:,:,ibNH4) => pore water ammonium concentration
!   sedPoreWaterCon(:,:,ibPO4) => pore water phosphate concentration
!   sedPoreWaterCon(:,:,ibSi)  => pore water silicate concentration
!   sedPoreWaterCon(:,:,ibTIC)  => pore water TIC concentration
!   sedPoreWaterCon(:,:,ibAlk)  => pore water alkalinity concentration

!   sedDecayRate(:,:,:,ibUC)   => decay rate of particulate carbon
!   sedDecayRate(:,:,:,ibUN)   => decay rate of particulate nitrogen
!   sedDecayRate(:,:,:,ibUP)   => decay rate of particulate phosphorous
!   sedDecayRate(:,:,:,ibUS)   => decay rate of particulate silicon

!   sedFlux(:,:,ibJNO3)  => flux of nitrate between sediment and water
!   sedFlux(:,:,ibJNH4)  => flux of ammonium between sediment and water
!   sedFlux(:,:,ibJPO4)  => flux of phosphate between sediment and water
!   sedFlux(:,:,ibJSi)  => flux of silicate between sediment and water
!   sedFlux(:,:,ibSOD)  => flux of oxygen between sediment and water
!   sedFlux(:,:,ibJTIC)  => flux of TIC between sediment and water
!   sedFlux(:,:,ibJAlk)  => flux of alkalinity between sediment and water

        real(r8), pointer :: sedPOM(:,:,:,:)
        real(r8), pointer :: sedPoreWaterCon(:,:,:)
        real(r8), pointer :: sedDecayRate(:,:,:,:)
        real(r8), pointer :: sedFlux(:,:,:)
        
      END TYPE T_SEDBIOL

      TYPE (T_SEDBIOL), allocatable :: SEDBIOL(:)

      CONTAINS

      SUBROUTINE allocate_sedbio_arrays (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
!      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Allocate Structure arrays
!
        IF (ng.eq.1) allocate ( SEDBIOL(Ngrids) )

!        write(*,*)  'allocate_sedbio_arrays: ng =',ng
!        write(*,*) 'allocate_sedbio_arrays:',LBi,UBi,LBj,UBj,nspc,NPOM

        allocate ( SEDBIOL(ng) % sedPOM(LBi:UBi,LBj:UBj,nspc,NPOM) )
        allocate ( SEDBIOL(ng) % sedPoreWaterCon(LBi:UBi,LBj:UBj,NPWC) )
        allocate ( SEDBIOL(ng) % sedDecayRate(LBi:UBi,LBj:UBj,nspc,NDR) )
        allocate ( SEDBIOL(ng) % sedFlux(LBi:UBi,LBj:UBj,NSF) )

      RETURN
      END SUBROUTINE allocate_sedbio_arrays

      SUBROUTINE initialize_sedbio_arrays (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the  shared     !
!  arrays  across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
!      USE mod_ncparam
!      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k, m

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize sediment biology structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN


        DO j=Jmin,Jmax
            DO i=Imin,Imax
                DO k=1,nspc
                    DO m=1,NPOM
                        SEDBIOL(ng) % sedPOM(i,j,k,m) = IniVal
                    END DO
                    DO m=1,NDR
                        SEDBIOL(ng) % sedDecayRate(i,j,k,m) = IniVal
                    END DO
                END DO
                DO m=1,NPWC
                    SEDBIOL(ng) % sedPoreWaterCon(i,j,m) = IniVal
                END DO
                DO m=1,NSF
                    SEDBIOL(ng) % sedFlux(i,j,m) = IniVal
                END DO
            END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE initialize_sedbio_arrays
