!                                                                      !
!  Parameters for UMaine CoSiNE model with added HAB phyto class:      !
!                                                                      !
!   reg1     Microzooplankton excretion rate to ammonium [1/day].      !
!   reg2     Mesozooplankton excretion rate to ammonium [1/day].       !
!   gmaxs1   Maximum specific growth rate of small phytoplankton       !
!              [1/day]                                                 !
!   gmaxs2   Maximum specific growth rate of diatom [1/day]            !
!   gmaxs3   Maximum specific growth rate of HAB phyto [1/day] 	     !
!! only if PHYTO_RESP defined:
!   rrb1	 Basal respiration rate of small phytoplankton	[1/day]    !
!   rrb2	 Basal respiration rate of diatom	[1/day]    		     !
!   rrb3	 Basal respiration rate of HAB phytoplankton	[1/day]    !
!   rrg1	 Growth-dependent respiration rate of small phyto [no unit]!
!   rrg2	 Growth-dependent respiration rate of diatom [no unit]     !
!   rrg3	 Growth-dependent respiration rate of HAB phyto [no unit]  !
!!
!   beta1    Microzooplankton maximum grazing rate [1/day]             !
!   beta2    Mesozooplankton maximum grazing rate [1/day]              !
!   akz1     Half saturation constant for microzooplankton grazing     !
!              [mmol_N/m3]                                             !
!   akz2     Half saturation constant for mesozooplankton grazing      !
!              [mmol_N/m3]                                             !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   amaxs1   Initial slope of P-I curve of small phytoplankton         !
!              [1/(Watts/m2)/day]                                      !
!   amaxs2   Initial slope of P-I curve of large phytoplankton         !
!              [1/(Watts/m2)/day]                                      !
!   amaxs3   Initial slope of P-I curve of HAB phytoplankton           !
!              [1/(Watts/m2)/day]                                      !
!   parsats1 PAR saturation onset parameter of small phytoplankton     !
!              [Watts/m2]                                              !
!   parsats2 PAR saturation onset parameter of diatom [Watts/m2]       !
!              [Watts/m2]                                              !
!   parsats3 PAR saturation onset parameter of HAB phyto [Watts/m2]    !
!              [Watts/m2]                                              !
!   pis1     Ammonium inhibition parameter for small phytoplankton     !
!              [mmol_N/m3]                                             !
!   pis2     Ammonium inhibition parameter for diatom [mmol_N/m3]      !
!   pis3     Ammonium inhibition parameter for HAB phyto [mmol_N/m3]   !
!   akno3s1  Half saturation concentration for nitrate uptake by       !
!              small phytoplankton [mmol_N/m3].                        !
!   akno3s2  Half saturation concentration for nitrate uptake by       !
!              diatom [mmol_N/m3].                                     !
!   akno3s3  Half saturation concentration for nitrate uptake by       !
!              HAB phyto [mmol_N/m3].                                  !
!   aknh4s1  Half saturation concentration for ammonium uptake by      !
!              small phytoplankton [mmol_N/m3].                        !
!   aknh4s2  Half saturation concentration for ammonium uptake by      !
!              diatom [mmol_N/m3].                                     !
!   aknh4s3  Half saturation concentration for ammonium uptake by      !
!              HAB phyto [mmol_N/m3].                                  !
!   akpo4s1  Half saturation concentration for phosphate uptake by     !
!              small phytoplankton [mmol_P/m3].                        !
!   akpo4s2  Half saturation concentration for phosphate uptake by     !
!              diatom [mmol_P/m3].                                     !
!   akpo4s3  Half saturation concentration for phosphate uptake by     !
!              HAB phyto [mmol_P/m3].                                  !
!   akco2s1  Half saturation concentration for co2 uptake by           !
!              small phytoplankton [mmol_C/m3].                        !
!   akco2s2  Half saturation concentration for co2 uptake by           !
!              diatom [mmol_C/m3].                                     !
!   akco2s3  Half saturation concentration for co2 uptake by           !
!              HAB phyto [mmol_C/m3].                                  !
!   aksio4s2 Half saturation constant for silicate uptake by           !
!              diatom [mmol_Si/m3].                                    !
!   aksio4s3 Half saturation constant for silicate uptake by           !
!              HAB phyto [mmol_Si/m3].                                 !
!   akox     Half saturation constant for oxidation                    !
!                 [mmol_O/m3].                                         !
!   ak1      Light attenuation coefficient of water [1/m]              !
!   ak2      Specific light attenuation coefficient for                !
!              phytoplankton [1/m/(mmol_N/m3)].                        !
!   bgamma0   Mesozooplankton specific mortality rate [1/day].         !
!   bgamma1   Grazing efficiency of microzooplankton [nondimensional]. !
!   bgamma2   Grazing efficiency of mesozooplankton [nondimensional].  !
!   bgamma3   Death rate of small phytoplankton [1/day].               !
!   bgamma4   Death rate of large phytoplankton [1/day].               !
!   bgamma4s3   Death rate of HAB phytoplankton [1/day].               !
!   bgamma5   Decay rate of detritus (nitrogenous) [1/day].            !
!   bgamma5s   Decay rate of detritus (silicious) [1/day].             !
!   bgamma6                                                            !
!   bgamma7   Nitrafication rate [1/day].                              !
!   wsd      Sinking velocity of detritus [m/day].                     !
!   wsdsi    Sinking velocity of detritus silicate [m/day].            !
!   wsp      Sinking velocity of large phytoplankton [m/day].          !
!   wsps3    Sinking velocity of HAB phytoplankton [m/day].            !
!   pco2a    Air pCO2 [ppmv].                                          !
!   si2n     Silicate to nitrogen ratio [mol_Si/mol_N].                !
!   p2n      Phosphorus to nitrogen ratio [mol_P/mol_N].               !
!   o2no     Oxygen to nitrate ratio [mol_O2/mol_NO3].                 !
!   o2nh     Oxygen to ammonium ratio [mol_O2/mol_NH4].                !
!   c2n      Carbon to nitrogen ratio [mol_C/mol_N].                   !
!   ro5      Grazing preference for diatom [nondimensional].           !
!   ro6      Grazing preference for mesozooplankton [nondimensional]   !
!   ro7      Grazing preference for detritus [nondimensional].         !
!   ro5H      Grazing preference for HAB phyto [nondimensional].         !
!   By: PENG XIU 12/2013        
!   Chl2cs1_m  Maximum chlorophyll to carbon ratio for sphytoplankton    
!                 [mg_Chl/mg_C]                                        !
!   Chl2cs2_m  Maximum chlorophyll to carbon ratio for Diatom
!                 [mg_Chl/mg_C]
!   Chl2cs3_m  Maximum chlorophyll to carbon ratio for HAB phyto
!                 [mg_Chl/mg_C]
#ifdef OPTICS_OP1
!   optic_upd_fac  number of timesteps between updates of kdPAR       !
#endif
!   q10_phyto_prod	Q10 factor for phytoplankton production
!   q10_phyto_resp	Q10 factor for phytoplankton respiration
!   q10_bact_remin	Q10 factor for bacterial remineralization
!   q10_opal_dissol	Q10 factor for opal dissolution
!   q10_zoopl_resp	Q10 factor for zooplankton respiration
!   q10_zoopl_graz	Q10 factor for zooplankton grazing
!   q10_nitrif		Q10 factor for nitrification
!   By: Xiu and Geng 9/2016                                    !
!=======================================================================
!
      USE mod_param
#ifdef SEDBIO
      USE mod_sedbio
#endif
      implicit none

!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:) ! Biological tracers
      integer :: iNO3_                 ! Nitrate concentration
      integer :: iNH4_                 ! Ammonium concentration
      integer :: iSiOH                 ! Silicate concentration
      integer :: iSphy                 ! Samll phytoplankton
      integer :: iLphy                 ! Diatom concentration
      integer :: iSzoo                 ! Small zooplankotn concentration
      integer :: iLzoo                 ! Mesozooplankotn concentration
      integer :: iSDet                 ! Detritus notrogen concentration
      integer :: iopal                 ! Biogenic silicate concentration
      integer :: iPO4_                 ! Phosphate concentration
#ifdef OXYGEN
      integer :: iOxyg                 ! Dissolved oxygen concentration
#endif
#ifdef CARBON
      integer :: iTIC_                 ! Total inorganic carbon
      integer :: iTAlk                 ! Total alkalinity
#endif
      integer :: iChl1                 ! Chlorophyll concentration for small phytoplankton
      integer :: iChl2                 ! Chlorophyll concentration for diatom
#ifdef HAB
      integer :: iHphy                 ! HAB phytoplankton
      integer :: iChl3                 ! Chlorophyll concentration for HAB phytoplankton
#endif
!
#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Biological 2D diagnostic variable IDs.
!
      integer, allocatable :: iDbio2(:)       ! 2D biological terms

      integer  :: iCOfx                       ! air-sea CO2 flux
      integer  :: ipCO2                       ! partial pressure of CO2
      integer  :: iO2fx                       ! air-sea O2 flux
# ifdef SEDBIO
      integer  :: ibNO3fx                     ! benthic NO3 flux
      integer  :: ibNH4fx                     ! benthic NH4 flux
      integer  :: ibPO4fx                     ! benthic PO4 flux
      integer  :: ibSiO2fx                    ! benthic SiO2 flux
      integer  :: ibPONfx                     ! flux of PON to sediment
      integer  :: ibPSifx                     ! flux of PSi to sediment
#  ifdef CARBON
      integer  :: ibTICfx                     ! benthic TIC flux
      integer  :: ibAlkfx                     ! benthic Alk flux
#  endif
#  ifdef OXYGEN
      integer  :: ibO2fx                     ! benthic O2 flux
#  endif
# endif
!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms

      integer  :: iPPro1 = 1                   ! primary productivity for s1
      integer  :: iPPro2 = 2                   ! primary productivity for s2
      integer  :: iNO3u = 3                   ! NO3 uptake
# ifdef OXYGEN
      integer  :: iO2pr = 4                    ! O2 production
      integer  :: initri = 5                   ! O2 loss, nitrification
      integer  :: iremin = 6                   ! O2 loss, remineralization
      integer  :: izoopl = 7                   ! O2 loss, zooplankton
#  ifdef PHYTO_RESP
	integer  :: iphyres = 8			     ! O2 loss, phytoplankton respiration
#  endif
# endif
#endif


      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: reg1(:)            ! 1/day
      real(r8), allocatable :: reg2(:)            ! 1/day
      real(r8), allocatable :: gmaxs1(:)          ! 1/day
      real(r8), allocatable :: gmaxs2(:)          ! 1/day
#ifdef PHYTO_RESP
	real(r8), allocatable :: rrb1(:)		  ! 1/day
	real(r8), allocatable :: rrb2(:)		  ! 1/day
	real(r8), allocatable :: rrg1(:)		  ! nondimensional
	real(r8), allocatable :: rrg2(:)		  ! nondimensional
#endif
      real(r8), allocatable :: beta1(:)           ! 1/day
      real(r8), allocatable :: beta2(:)           ! 1/day
      real(r8), allocatable :: akz1(:)            ! mmol_N/m3
      real(r8), allocatable :: akz2(:)            ! mmol_N/m3
      real(r8), allocatable :: PARfrac(:)         ! nondimensional
      real(r8), allocatable :: amaxs1(:)          ! 1/(Watts/m2)/day
      real(r8), allocatable :: amaxs2(:)          ! 1/(Watts/m2)/day
      real(r8), allocatable :: pis1(:)            ! m3/mmol_N
      real(r8), allocatable :: pis2(:)            ! m3/mmol_N
      real(r8), allocatable :: akno3s1(:)         ! mmol_N/m3
      real(r8), allocatable :: aknh4s1(:)         ! mmol_N/m3
      real(r8), allocatable :: akpo4s1(:)         ! mmol_P/m3
      real(r8), allocatable :: akco2s1(:)         ! mmol_C/m3
      real(r8), allocatable :: akno3s2(:)         ! mmol_N/m3
      real(r8), allocatable :: aknh4s2(:)         ! mmol_N/m3
      real(r8), allocatable :: aksio4s2(:)        ! mmol_Si/m3
      real(r8), allocatable :: akpo4s2(:)         ! mmol_P/m3
      real(r8), allocatable :: akco2s2(:)         ! mmol_C/m3
      real(r8), allocatable :: akox(:)            ! mmol_O/m3
      real(r8), allocatable :: ak1(:)             ! 1/m
      real(r8), allocatable :: ak2(:)             ! 1/m/(mmol_N/m3)
      real(r8), allocatable :: parsats1(:)        ! Watts/m2
      real(r8), allocatable :: parsats2(:)        ! Watts/m2
      real(r8), allocatable :: bgamma0(:)         ! 1/day
      real(r8), allocatable :: bgamma1(:)         ! [nondimensional]
      real(r8), allocatable :: bgamma2(:)         ! [nondimensional]
      real(r8), allocatable :: bgamma3(:)         ! 1/day
      real(r8), allocatable :: bgamma4(:)         ! 1/day
      real(r8), allocatable :: bgamma5(:)         ! 1/day
      real(r8), allocatable :: bgamma5s(:)        ! 1/day
      real(r8), allocatable :: bgamma6(:)         !
      real(r8), allocatable :: bgamma7(:)         ! 1/day
      real(r8), allocatable :: wsd(:)             ! m/day
      real(r8), allocatable :: wsdsi(:)           ! m/day
      real(r8), allocatable :: wsp(:)             ! m/day
      real(r8), allocatable :: si2n(:)            ! mol_Si/mol_N
      real(r8), allocatable :: pco2a(:)           ! ppmv
      real(r8), allocatable :: p2n(:)             ! mol_P/mol_N
      real(r8), allocatable :: o2no(:)            ! mol_O2/mol_NO3
      real(r8), allocatable :: o2nh(:)            ! mol_O2/mol_NH4
      real(r8), allocatable :: c2n(:)             ! mol_C/mol_N
      real(r8), allocatable :: ro5(:)             ! nondimensional
      real(r8), allocatable :: ro6(:)             ! nondimensional
      real(r8), allocatable :: ro7(:)             ! nondimensional
      real(r8), allocatable :: pCO2air(:)         ! ppmv
      real(r8), allocatable :: Chl2cs1_m(:)         ! mg_Chl/mg_C
      real(r8), allocatable :: Chl2cs2_m(:)         ! mg_Chl/mg_C
#ifdef HAB
      real(r8), allocatable :: gmaxs3(:)          ! 1/day
      real(r8), allocatable :: amaxs3(:)          ! 1/(Watts/m2)/day
#  ifdef PHYTO_RESP
	real(r8), allocatable :: rrb3(:)		  ! 1/day
	real(r8), allocatable :: rrg3(:)		  ! nondimensional
#  endif
      real(r8), allocatable :: parsats3(:)        ! Watts/m2
      real(r8), allocatable :: pis3(:)            ! m3/mmol_N
      real(r8), allocatable :: akno3s3(:)         ! mmol_N/m3
      real(r8), allocatable :: aknh4s3(:)         ! mmol_N/m3
      real(r8), allocatable :: akpo4s3(:)         ! mmol_P/m3
      real(r8), allocatable :: akco2s3(:)         ! mmol_C/m3
      real(r8), allocatable :: aksio4s3(:)        ! mmol_Si/m3
      real(r8), allocatable :: bgamma4s3(:)         ! 1/day
      real(r8), allocatable :: wsp3(:)             ! m/day
      real(r8), allocatable :: ro5H(:)             ! nondimensional
      real(r8), allocatable :: Chl2cs3_m(:)         ! mg_Chl/mg_C
#endif
#ifdef OPTICS_OP1
      real(r8), allocatable :: optic_upd_fac(:)    ! nondimensional
#endif
	real(r8), allocatable :: q10_phyto_prod(:)	   ! nondimensional
	real(r8), allocatable :: q10_phyto_resp(:)	   ! nondimensional
	real(r8), allocatable :: q10_bact_remin(:)	   ! nondimensional
	real(r8), allocatable :: q10_opal_dissol(:)	   ! nondimensional
	real(r8), allocatable :: q10_zoopl_resp(:)	   ! nondimensional
	real(r8), allocatable :: q10_zoopl_graz(:)	   ! nondimensional
	real(r8), allocatable :: q10_nitrif(:)		   ! nondimensional
! a1-a7 are the constants in Q10=exp(a*(T-Tref)) which are computed
!	  from the input q10 values
	real(r8), allocatable :: a1(:)			   ! 1/degC
	real(r8), allocatable :: a2(:)			   ! 1/degC
	real(r8), allocatable :: a3(:)			   ! 1/degC
	real(r8), allocatable :: a4(:)			   ! 1/degC
	real(r8), allocatable :: a5(:)			   ! 1/degC
      real(r8), allocatable :: a6(:)			   ! 1/degC
	real(r8), allocatable :: a7(:)			   ! 1/degC

#ifdef SEDBIO
! benthic biology parameters
      real(r8), allocatable :: bUmax(:,:)
      real(r8), allocatable :: bUmaxSi(:,:)
      real(r8), allocatable :: bdep(:)
      real(r8), allocatable :: balpha(:)
      real(r8), allocatable :: bw(:)
      real(r8), allocatable :: btheta_diag(:)
      real(r8), allocatable :: bnit(:)
      real(r8), allocatable :: btheta_nit(:)
      real(r8), allocatable :: bdo_nit(:)
      real(r8), allocatable :: bdenit(:)
      real(r8), allocatable :: btheta_denit(:)
      real(r8), allocatable :: bpsi_n(:)
      real(r8), allocatable :: bdo_c(:)
      real(r8), allocatable :: bao2(:)
      real(r8), allocatable :: bpi(:)
      real(r8), allocatable :: bpsi_p(:)
      real(r8), allocatable :: bfc(:,:)
      real(r8), allocatable :: bfn(:,:)
      real(r8), allocatable :: bfp(:,:)
      real(r8), allocatable :: bfs(:,:)
#endif
      CONTAINS

      SUBROUTINE initialize_biology
#ifdef SEDBIO
      USE mod_sedbio, ONLY: nspc, NPOM, NPWC, NDR, NSF, NBBT
#endif
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
#ifdef CARBON
# ifdef OXYGEN
      NBT=15
# else
      NBT=14
# endif
#else
# ifdef OXYGEN
      NBT=13
# else
      NBT=12
# endif
#endif
#ifdef HAB
      NBT=NBT+2
#endif
#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
      NDbio3d=3
# ifdef OXYGEN
        NDbio3d=NDbio3d+4
#  ifdef PHYTO_RESP
	NDbio3d=NDbio3d+1
#  endif
# endif

      NDbio2d=0
! benthic flux diagnostics:
# ifdef SEDBIO
      NDbio2d=NDbio2d+6
#  ifdef CARBON
      NDbio2d=NDbio2d+2
#  endif
#  ifdef OXYGEN
      NDbio2d=NDbio2d+1
#  endif
# endif

# ifdef CARBON
      NDbio2d=NDbio2d+2
# endif                                                                                               
# ifdef OXYGEN
      NDbio2d=NDbio2d+1
# endif

     ic=0
# ifdef SEDBIO
      ibNO3fx=ic+1
      ibNH4fx=ic+2
      ibPO4fx=ic+3
      ibSiO2fx=ic+4
      ibPONfx=ic+5
      ibPSifx=ic+6
      ic=ic+6
# endif
# ifdef CARBON
      iCOfx=ic+1
      ipCO2=ic+2
      ic=ic+2
#  ifdef SEDBIO
      ibTICfx=ic+1
      ibAlkfx=ic+2
      ic=ic+2
#  endif
# endif
# ifdef OXYGEN
      iO2fx=ic+1
      ic=ic+1
#  ifdef SEDBIO
      ibO2fx=ic+1
      ic=ic+1
#  endif
# endif
#endif
#ifdef SEDBIO
! total number of sediment biology variables
      NBBT=NPOM+NPWC+NDR+NSF
#endif
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(reg1)) THEN
        allocate ( reg1(Ngrids) )
      END IF
      IF (.not.allocated(reg2)) THEN
        allocate ( reg2(Ngrids) )
      END IF
      IF (.not.allocated(gmaxs1)) THEN
        allocate ( gmaxs1(Ngrids) )
      END IF
      IF (.not.allocated(gmaxs2)) THEN
        allocate ( gmaxs2(Ngrids) )
      END IF
#ifdef PHYTO_RESP
	IF (.not.allocated(rrb1)) THEN
	  allocate ( rrb1(Ngrids) )
	END IF
	IF (.not.allocated(rrb2)) THEN
	  allocate ( rrb2(Ngrids) )
	END IF
	IF (.not.allocated(rrg1)) THEN
	  allocate ( rrg1(Ngrids) )
	END IF
	IF (.not.allocated(rrg2)) THEN
	  allocate ( rrg2(Ngrids) )
	END IF
#endif
      IF (.not.allocated(beta1)) THEN
        allocate ( beta1(Ngrids) )
      END IF
      IF (.not.allocated(beta2)) THEN
        allocate ( beta2(Ngrids) )
      END IF
      IF (.not.allocated(akz1)) THEN
        allocate ( akz1(Ngrids) )
      END IF
      IF (.not.allocated(akz2)) THEN
        allocate ( akz2(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(amaxs1)) THEN
        allocate ( amaxs1(Ngrids) )
      END IF
      IF (.not.allocated(amaxs2)) THEN
        allocate ( amaxs2(Ngrids) )
      END IF
      IF (.not.allocated(pis1)) THEN
        allocate ( pis1(Ngrids) )
      END IF
      IF (.not.allocated(pis2)) THEN
        allocate ( pis2(Ngrids) )
      END IF
      IF (.not.allocated(akno3s1)) THEN
        allocate ( akno3s1(Ngrids) )
      END IF
      IF (.not.allocated(aknh4s1)) THEN
        allocate ( aknh4s1(Ngrids) )
      END IF
      IF (.not.allocated(akpo4s1)) THEN
        allocate ( akpo4s1(Ngrids) )
      END IF
      IF (.not.allocated(akco2s1)) THEN
        allocate ( akco2s1(Ngrids) )
      END IF
      IF (.not.allocated(akno3s2)) THEN
        allocate ( akno3s2(Ngrids) )
      END IF
      IF (.not.allocated(aknh4s2)) THEN
        allocate ( aknh4s2(Ngrids) )
      END IF
      IF (.not.allocated(aksio4s2)) THEN
        allocate ( aksio4s2(Ngrids) )
      END IF
      IF (.not.allocated(akpo4s2)) THEN
        allocate ( akpo4s2(Ngrids) )
      END IF
      IF (.not.allocated(akco2s2)) THEN
        allocate ( akco2s2(Ngrids) )
      END IF
      IF (.not.allocated(akox)) THEN
        allocate ( akox(Ngrids) )
      END IF
      IF (.not.allocated(ak1)) THEN
        allocate ( ak1(Ngrids) )
      END IF
      IF (.not.allocated(ak2)) THEN
        allocate ( ak2(Ngrids) )
      END IF
      IF (.not.allocated(parsats1)) THEN
        allocate ( parsats1(Ngrids) )
      END IF
      IF (.not.allocated(parsats2)) THEN
        allocate ( parsats2(Ngrids) )
      END IF
      IF (.not.allocated(bgamma0)) THEN
        allocate ( bgamma0(Ngrids) )
      END IF
      IF (.not.allocated(bgamma1)) THEN
        allocate ( bgamma1(Ngrids) )
      END IF
      IF (.not.allocated(bgamma2)) THEN
        allocate ( bgamma2(Ngrids) )
      END IF
      IF (.not.allocated(bgamma3)) THEN
        allocate ( bgamma3(Ngrids) )
      END IF
      IF (.not.allocated(bgamma4)) THEN
        allocate ( bgamma4(Ngrids) )
      END IF
      IF (.not.allocated(bgamma5)) THEN
        allocate ( bgamma5(Ngrids) )
      END IF
      IF (.not.allocated(bgamma5s)) THEN
        allocate ( bgamma5s(Ngrids) )
      END IF
      IF (.not.allocated(bgamma6)) THEN
        allocate ( bgamma6(Ngrids) )
      END IF
      IF (.not.allocated(bgamma7)) THEN
        allocate ( bgamma7(Ngrids) )
      END IF
      IF (.not.allocated(wsd)) THEN
        allocate ( wsd(Ngrids) )
      END IF
      IF (.not.allocated(wsdsi)) THEN
        allocate ( wsdsi(Ngrids) )
      END IF
      IF (.not.allocated(wsp)) THEN
        allocate ( wsp(Ngrids) )
      END IF
      IF (.not.allocated(si2n)) THEN
        allocate ( si2n(Ngrids) )
      END IF
      IF (.not.allocated(pco2a)) THEN
        allocate ( pco2a(Ngrids) )
      END IF
      IF (.not.allocated(p2n)) THEN
        allocate ( p2n(Ngrids) )
      END IF
      IF (.not.allocated(o2no)) THEN
        allocate ( o2no(Ngrids) )
      END IF
      IF (.not.allocated(o2nh)) THEN
        allocate ( o2nh(Ngrids) )
      END IF
      IF (.not.allocated(c2n)) THEN
        allocate ( c2n(Ngrids) )
      END IF
      IF (.not.allocated(ro5)) THEN
        allocate ( ro5(Ngrids) )
      END IF
      IF (.not.allocated(ro6)) THEN
        allocate ( ro6(Ngrids) )
      END IF
      IF (.not.allocated(ro7)) THEN
        allocate ( ro7(Ngrids) )
      END IF
      IF (.not.allocated(pCO2air)) THEN
        allocate ( pCO2air(Ngrids) )
      END if
      if (.not.allocated(Chl2cs1_m)) THEN
        allocate ( Chl2cs1_m(Ngrids) )
      END if
      if (.not.allocated(Chl2cs2_m)) THEN
        allocate ( Chl2cs2_m(Ngrids) )
      END IF
#ifdef HAB
      IF (.not.allocated(gmaxs3)) THEN
        allocate ( gmaxs3(Ngrids) )
      END IF
      IF (.not.allocated(amaxs3)) THEN
        allocate ( amaxs3(Ngrids) )
      END IF
#  ifdef PHYTO_RESP
	IF (.not.allocated(rrb3)) THEN
	  allocate ( rrb3(Ngrids) )
	END IF
	IF (.not.allocated(rrg3)) THEN
	  allocate ( rrg3(Ngrids) )
	END IF
#  endif
      IF (.not.allocated(parsats3)) THEN
        allocate ( parsats3(Ngrids) )
      END IF
      IF (.not.allocated(pis3)) THEN
        allocate ( pis3(Ngrids) )
      END IF
      IF (.not.allocated(akno3s3)) THEN
        allocate ( akno3s3(Ngrids) )
      END IF
      IF (.not.allocated(aknh4s3)) THEN
        allocate ( aknh4s3(Ngrids) )
      END IF
      IF (.not.allocated(akpo4s3)) THEN
        allocate ( akpo4s3(Ngrids) )
      END IF
      IF (.not.allocated(akco2s3)) THEN
        allocate ( akco2s3(Ngrids) )
      END IF
      IF (.not.allocated(aksio4s3)) THEN
        allocate ( aksio4s3(Ngrids) )
      END IF
      IF (.not.allocated(Chl2cs3_m)) THEN
        allocate ( Chl2cs3_m(Ngrids) )
      END IF
      IF (.not.allocated(ro5H)) THEN
        allocate ( ro5H(Ngrids) )
      END IF
      IF (.not.allocated(wsp3)) THEN
        allocate ( wsp3(Ngrids) )
      END IF
      IF (.not.allocated(bgamma4s3)) THEN
        allocate ( bgamma4s3(Ngrids) )
      END IF
#endif
#ifdef OPTICS_OP1
      if (.not.allocated(optic_upd_fac)) THEN
        allocate ( optic_upd_fac(Ngrids) )
      END IF
#endif
#ifdef SEDBIO
      if (.not.allocated(bUmax)) THEN
        allocate ( bUmax(Ngrids,nspc) )
      END IF
      if (.not.allocated(bUmaxSi)) THEN
        allocate ( bUmaxSi(Ngrids,nspc) )
      END IF
      if (.not.allocated(bdep)) THEN
        allocate ( bdep(Ngrids) )
      END IF
      if (.not.allocated(balpha)) THEN
        allocate ( balpha(Ngrids) )
      END IF
      if (.not.allocated(bw)) THEN
        allocate ( bw(Ngrids) )
      END IF
      if (.not.allocated(btheta_diag)) THEN
        allocate ( btheta_diag(Ngrids) )
      END IF
      if (.not.allocated(bnit)) THEN
        allocate ( bnit(Ngrids) )
      END IF
      if (.not.allocated(btheta_nit)) THEN
        allocate ( btheta_nit(Ngrids) )
      END IF
      if (.not.allocated(bdo_nit)) THEN
        allocate ( bdo_nit(Ngrids) )
      END IF
      if (.not.allocated(bdenit)) THEN
        allocate ( bdenit(Ngrids) )
      END IF
      if (.not.allocated(btheta_denit)) THEN
        allocate ( btheta_denit(Ngrids) )
      END IF
      if (.not.allocated(bpsi_n)) THEN
        allocate ( bpsi_n(Ngrids) )
      END IF
      if (.not.allocated(bdo_c)) THEN
        allocate ( bdo_c(Ngrids) )
      END IF
      if (.not.allocated(bao2)) THEN
        allocate ( bao2(Ngrids) )
      END IF
      if (.not.allocated(bpi)) THEN
        allocate ( bpi(Ngrids) )
      END IF
      if (.not.allocated(bpsi_p)) THEN
        allocate ( bpsi_p(Ngrids) )
      END IF
      if (.not.allocated(bfc)) THEN
        allocate ( bfc(Ngrids,nspc) )
      END IF
      if (.not.allocated(bfn)) THEN
        allocate ( bfn(Ngrids,nspc) )
      END IF
      if (.not.allocated(bfp)) THEN
        allocate ( bfp(Ngrids,nspc) )
      END IF
      if (.not.allocated(bfs)) THEN
        allocate ( bfs(Ngrids,nspc) )
      END IF
#endif
	IF (.not.allocated(q10_phyto_prod)) THEN
	  allocate ( q10_phyto_prod(Ngrids) )
	END IF
	IF (.not.allocated(q10_phyto_resp)) THEN
	  allocate ( q10_phyto_resp(Ngrids) )
	END IF
	IF (.not.allocated(q10_bact_remin)) THEN
	  allocate ( q10_bact_remin(Ngrids) )
	END IF
	IF (.not.allocated(q10_opal_dissol)) THEN
	  allocate ( q10_opal_dissol(Ngrids) )
	END IF
	IF (.not.allocated(q10_zoopl_resp)) THEN
	  allocate ( q10_zoopl_resp(Ngrids) )
	END IF
	IF (.not.allocated(q10_zoopl_graz)) THEN
	  allocate ( q10_zoopl_graz(Ngrids) )
	END IF
	IF (.not.allocated(q10_nitrif)) THEN
	  allocate ( q10_nitrif(Ngrids) )
	END IF
	IF (.not.allocated(a1)) THEN
	  allocate ( a1(Ngrids) )
	END IF
	IF (.not.allocated(a2)) THEN
	  allocate ( a2(Ngrids) )
	END IF
	IF (.not.allocated(a3)) THEN
	  allocate ( a3(Ngrids) )
	END IF
	IF (.not.allocated(a4)) THEN
	  allocate ( a4(Ngrids) )
	END IF
	IF (.not.allocated(a5)) THEN
	  allocate ( a5(Ngrids) )
	END IF
	IF (.not.allocated(a6)) THEN
	  allocate ( a6(Ngrids) )
	END IF
	IF (.not.allocated(a7)) THEN
	  allocate ( a7(Ngrids) )
	END IF
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!     
!  Allocate biological diagnostics vectors
!       
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
      END IF
#endif
!
!  Set identification indices.
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iNO3_=ic+1
      iNH4_=ic+2
      iSiOH=ic+3
      iSphy=ic+4
      iLphy=ic+5
      iSzoo=ic+6
      iLzoo=ic+7
      iSDet=ic+8
      iopal=ic+9
      iPO4_=ic+10
      iChl1=ic+11
      iChl2=ic+12
      ic=ic+12
# ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
# endif
# ifdef CARBON
      iTIC_=ic+1
      iTAlk=ic+2
      ic=ic+2
# endif
# ifdef HAB
      iHphy=ic+1
      iChl3=ic+2
      ic=ic+2
# endif

      RETURN
      END SUBROUTINE initialize_biology
