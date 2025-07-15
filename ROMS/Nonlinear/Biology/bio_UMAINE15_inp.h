      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!  By: PENG XIU 12/2013                                                !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
#ifdef SEDBIO
      USE mod_iounits
#endif
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: find_file
      integer :: Npts, Nval, i, itrc, ng, status

      integer :: igrid, itracer,iTrcStr, iTrcEnd,nline,ifield
      
      integer :: decode_line, load_i, load_l, load_r,load_lbc

      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=1  ), parameter :: blank = ' '
      character (len=40 ) :: KeyWord
      character (len=256) :: line
!     character (len=256), dimension(100) :: Cval
      character (len=256), dimension(200) :: Cval
      character (len=256) :: fname
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!
!-----------------------------------------------------------------------
!  Read in UMaine CoSiNE biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
        SELECT CASE (TRIM(KeyWord))
          CASE('Lbiology')
             Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          CASE('BioIter') 
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          CASE('reg1') 
            Npts=load_r(Nval, Rval, Ngrids, reg1)
          CASE('reg2') 
            Npts=load_r(Nval, Rval, Ngrids, reg2)
          CASE('gmaxs1') 
            Npts=load_r(Nval, Rval, Ngrids, gmaxs1)
          CASE('gmaxs2') 
            Npts=load_r(Nval, Rval, Ngrids, gmaxs2)
          CASE('beta1') 
            Npts=load_r(Nval, Rval, Ngrids, beta1)
          CASE('beta2') 
            Npts=load_r(Nval, Rval, Ngrids, beta2)
          CASE('akz1') 
            Npts=load_r(Nval, Rval, Ngrids, akz1)
          CASE('akz2') 
            Npts=load_r(Nval, Rval, Ngrids, akz2)
          CASE('PARfrac') 
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          CASE('amaxs1') 
            Npts=load_r(Nval, Rval, Ngrids, amaxs1)
          CASE('amaxs2') 
            Npts=load_r(Nval, Rval, Ngrids, amaxs2)
          CASE('parsats1') 
            Npts=load_r(Nval, Rval, Ngrids, parsats1)
          CASE('parsats2') 
            Npts=load_r(Nval, Rval, Ngrids, parsats2)
          CASE('pis1') 
            Npts=load_r(Nval, Rval, Ngrids, pis1)
          CASE('pis2') 
            Npts=load_r(Nval, Rval, Ngrids, pis2)
          CASE('akno3s1') 
            Npts=load_r(Nval, Rval, Ngrids, akno3s1)
          CASE('akno3s2') 
            Npts=load_r(Nval, Rval, Ngrids, akno3s2)
          CASE('aknh4s1') 
            Npts=load_r(Nval, Rval, Ngrids, aknh4s1)
          CASE('aknh4s2') 
            Npts=load_r(Nval, Rval, Ngrids, aknh4s2)
          CASE('akpo4s1') 
            Npts=load_r(Nval, Rval, Ngrids, akpo4s1)
          CASE('akpo4s2') 
            Npts=load_r(Nval, Rval, Ngrids, akpo4s2)
          CASE('akco2s1') 
            Npts=load_r(Nval, Rval, Ngrids, akco2s1)
          CASE('akco2s2') 
            Npts=load_r(Nval, Rval, Ngrids, akco2s2)
          CASE('aksio4s2') 
            Npts=load_r(Nval, Rval, Ngrids, aksio4s2)
          CASE('akox') 
            Npts=load_r(Nval, Rval, Ngrids, akox)
          CASE('ak1') 
            Npts=load_r(Nval, Rval, Ngrids, ak1)
          CASE('ak2') 
            Npts=load_r(Nval, Rval, Ngrids, ak2)
          CASE('bgamma0') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma0)
          CASE('bgamma1') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma1)
          CASE('bgamma2') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma2)
          CASE('bgamma3') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma3)
          CASE('bgamma4') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma4)
          CASE('bgamma5') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma5)
          CASE('bgamma5s')
            Npts=load_r(Nval, Rval, Ngrids, bgamma5s)
          CASE('bgamma6')
            Npts=load_r(Nval, Rval, Ngrids, bgamma6)
          CASE('bgamma7') 
            Npts=load_r(Nval, Rval, Ngrids, bgamma7)
          CASE('wsd') 
            Npts=load_r(Nval, Rval, Ngrids, wsd)
          CASE('wsdsi') 
            Npts=load_r(Nval, Rval, Ngrids, wsdsi)
          CASE('wsp') 
            Npts=load_r(Nval, Rval, Ngrids, wsp)
          CASE('pco2a') 
            Npts=load_r(Nval, Rval, Ngrids, pco2a)
          CASE('si2n') 
            Npts=load_r(Nval, Rval, Ngrids, si2n)
          CASE('p2n') 
            Npts=load_r(Nval, Rval, Ngrids, p2n)
          CASE('o2no') 
            Npts=load_r(Nval, Rval, Ngrids, o2no)
          CASE('o2nh') 
            Npts=load_r(Nval, Rval, Ngrids, o2nh)
          CASE('c2n') 
            Npts=load_r(Nval, Rval, Ngrids, c2n)
          CASE('ro5') 
            Npts=load_r(Nval, Rval, Ngrids, ro5)
          CASE('ro6') 
            Npts=load_r(Nval, Rval, Ngrids, ro6)
          CASE('ro7') 
            Npts=load_r(Nval, Rval, Ngrids, ro7)
          CASE('Chl2cs1_m') 
            Npts=load_r(Nval, Rval, Ngrids, Chl2cs1_m)
          CASE('Chl2cs2_m') 
            Npts=load_r(Nval, Rval, Ngrids, Chl2cs2_m)
#ifdef OPTICS_OP1
          CASE('optic_upd_fac')
            Npts=load_r(Nval, Rval, Ngrids, optic_upd_fac)
#endif
#ifdef PHYTO_RESP
          CASE('rrb1')
            Npts=load_r(Nval, Rval, Ngrids, rrb1)
          CASE('rrb2')
            Npts=load_r(Nval, Rval, Ngrids, rrb2)
          CASE('rrg1')
            Npts=load_r(Nval, Rval, Ngrids, rrg1)
          CASE('rrg2')
            Npts=load_r(Nval, Rval, Ngrids, rrg2)
#endif
#ifdef HAB
          CASE('gmaxs3')
            Npts=load_r(Nval, Rval, Ngrids, gmaxs3)
          CASE('amaxs3')
            Npts=load_r(Nval, Rval, Ngrids, amaxs3)
          CASE('parsats3')
            Npts=load_r(Nval, Rval, Ngrids, parsats3)
          CASE('pis3')
            Npts=load_r(Nval, Rval, Ngrids, pis3)
          CASE('akno3s3')
            Npts=load_r(Nval, Rval, Ngrids, akno3s3)
          CASE('aknh4s3')
            Npts=load_r(Nval, Rval, Ngrids, aknh4s3)
          CASE('akpo4s3')
            Npts=load_r(Nval, Rval, Ngrids, akpo4s3)
          CASE('akco2s3')
            Npts=load_r(Nval, Rval, Ngrids, akco2s3)
          CASE('aksio4s3')
            Npts=load_r(Nval, Rval, Ngrids, aksio4s3)
          CASE('bgamma4s3')
            Npts=load_r(Nval, Rval, Ngrids, bgamma4s3)
          CASE('wsp3')
            Npts=load_r(Nval, Rval, Ngrids, wsp3)
          CASE('ro5H')
            Npts=load_r(Nval, Rval, Ngrids, ro5H)
          CASE('Chl2cs3_m')
            Npts=load_r(Nval, Rval, Ngrids, Chl2cs3_m)
#ifdef PHYTO_RESP
          CASE('rrb3')
            Npts=load_r(Nval, Rval, Ngrids, rrb3)
          CASE('rrg3')
            Npts=load_r(Nval, Rval, Ngrids, rrg3)
#endif
#endif
	    CASE('q10_phyto_prod')
		Npts=load_r(Nval, Rval, Ngrids, q10_phyto_prod)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a1(ng)=log(q10_phyto_prod(ng))/10.0_r8
		END DO
	    CASE('q10_phyto_resp')
		Npts=load_r(Nval, Rval, Ngrids, q10_phyto_resp)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a2(ng)=log(q10_phyto_resp(ng))/10.0_r8
		END DO
	    CASE('q10_bact_remin')
		Npts=load_r(Nval, Rval, Ngrids, q10_bact_remin)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a3(ng)=log(q10_bact_remin(ng))/10.0_r8
		END DO
	    CASE('q10_opal_dissol')
		Npts=load_r(Nval, Rval, Ngrids, q10_opal_dissol)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a4(ng)=log(q10_opal_dissol(ng))/10.0_r8
		END DO
	    CASE('q10_zoopl_resp')
		Npts=load_r(Nval, Rval, Ngrids, q10_zoopl_resp)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a5(ng)=log(q10_zoopl_resp(ng))/10.0_r8
		END DO
	    CASE('q10_zoopl_graz')
		Npts=load_r(Nval, Rval, Ngrids, q10_zoopl_graz)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a6(ng)=log(q10_zoopl_graz(ng))/10.0_r8
		END DO
	    CASE('q10_nitrif')
		Npts=load_r(Nval, Rval, Ngrids, q10_nitrif)
!	compute exponential constant to get this q10
		DO ng=1,Ngrids
		  a7(ng)=log(q10_nitrif(ng))/10.0_r8
		END DO

#ifdef SEDBIO
          CASE ('BSEDPARNAM')
            DO i=1,LEN(bsedparnam)
                bsedparnam(i:i)=blank
            END DO
            bsedparnam=TRIM(ADJUSTL(Cval(Nval)))
#endif
          CASE('TNU2') 
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                nl_tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          CASE('TNU4') 
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                nl_tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu2(i,ng)=Rbio(itrc,ng)
                  tl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu4(i,ng)=Rbio(itrc,ng)
                  tl_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
          CASE('AKT_BAK') 
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          CASE('TNUDG') 
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
            
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,    &
     &                      iTrcStr, iTrcEnd,                          &
     &                      Vname(1,idTvar(idbio(itracer))), LBC)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,    &
     &                      iTrcStr, iTrcEnd,                          &
     &                      Vname(1,idTvar(idbio(itracer))), ad_LBC)
#endif
#ifdef TCLIMATOLOGY
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTvar)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                         &
     &                                'idTvar(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTsur)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                         &
     &                                'idTsur(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#ifdef OPTICS_OP1
            CASE ('Hout(idkdPAR)')
              IF (idkdPAR.eq.0) THEN
                IF (Master) WRITE (out,30) 'idkdPAR'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idkdPAR,:))
#endif
! AKB 6/22/23
#ifdef CARBON
            CASE ('Hout(idfgCO2)')
              IF (idfgCO2.eq.0) THEN
                IF (Master) WRITE (out,30) 'idfgCO2'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idfgCO2,:))
            CASE ('Hout(idpCO2s)')
              IF (idpCO2s.eq.0) THEN
                IF (Master) WRITE (out,30) 'idpCO2s'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idpCO2s,:))
#endif

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            CASE ('Aout(idTvar)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idTTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idUTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idVTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idVTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHUTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=iHUTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHVTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=iHVTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
#ifdef DIAGNOSTICS_TS
            CASE ('Dout(iTrate)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTrate),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iThadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iThadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTxadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTyadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTvadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltrc(i,ng)
                END DO
              END DO
# if defined TS_DIF2 || defined TS_DIF4
            CASE ('Dout(iThdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iThdif),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTxdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTydif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTydif),ng)=Ltrc(i,ng)
                END DO
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            CASE ('Dout(iTsdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltrc(i,ng)
                END DO
              END DO
#  endif
# endif
            CASE ('Dout(iTvdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltrc(i,ng)
                END DO
              END DO
#endif
#ifdef DIAGNOSTICS_BIO
# ifdef SEDBIO
            CASE ('Dout(ibNO3fx)')
              IF (iDbio2(ibNO3fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibNO3fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibNO3fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ibNH4fx)')
              IF (iDbio2(ibNH4fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibNH4fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibNH4fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ibPO4fx)')
              IF (iDbio2(ibPO4fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibPO4fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibPO4fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ibSiO2fx)')
              IF (iDbio2(ibSiO2fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibSiO2fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibSiO2fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ibPONfx)')
              IF (iDbio2(ibPONfx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibPONfx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibPONfx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ibPSifx)')
              IF (iDbio2(ibPSifx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibPSifx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibPSifx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
#  ifdef CARBON
            CASE ('Dout(ibTICfx)')
              IF (iDbio2(ibTICfx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibTICfx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibTICfx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ibAlkfx)')
              IF (iDbio2(ibAlkfx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibAlkfx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibAlkfx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
#  endif
#  ifdef OXYGEN
            CASE ('Dout(ibO2fx)')
              IF (iDbio2(ibO2fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ibO2fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ibO2fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
#  endif
# endif

# ifdef CARBON
            CASE ('Dout(iCOfx)')
              IF (iDbio2(iCOfx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(iCOfx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(iCOfx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif

# ifdef CARBON
            CASE ('Dout(ipCO2)')
              IF (iDbio2(ipCO2).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ipCO2)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ipCO2)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
# ifdef OXYGEN
            CASE ('Dout(iO2fx)')
              IF (iDbio2(iO2fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(iO2fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(iO2fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
            CASE ('Dout(iPPro1)')
              IF (iDbio3(iPPro1).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iPPro1)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iPPro1)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
           
            CASE ('Dout(iPPro2)')
              IF (iDbio3(iPPro2).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iPPro2)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iPPro2)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO

            CASE ('Dout(iNO3u)')
              IF (iDbio3(iNO3u).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iNO3u)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iNO3u)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# ifdef OXYGEN
            CASE ('Dout(iO2pr)')
              IF (iDbio3(iO2pr).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iO2pr)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iO2pr)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            
            CASE ('Dout(initri)')
              IF (iDbio3(initri).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(initri)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(initri)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            
            CASE ('Dout(iremin)')
              IF (iDbio3(iremin).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iremin)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iremin)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO

            CASE ('Dout(izoopl)')
              IF (iDbio3(izoopl).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(izoopl)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(izoopl)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
#  ifdef PHYTO_RESP
            CASE ('Dout(iphyres)')
              IF (iDbio3(iphyres).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iphyres)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iphyres)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
#  endif
# endif
#endif
 
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                     &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,110) reg1(ng), 'reg1',                          &
     &            'Microzooplankton excretion rate to ammonium',       &
     &            '[1/day].'
            WRITE (out,100) reg2(ng), 'reg2',                          &
     &            'Mesozooplankton excretion rate to ammonium [1/day].'
            WRITE (out,110) gmaxs1(ng), 'gmaxs1',                      &
     &           'Maximum specific growth rate of small phytoplankton',&
     &            '[1/day].'
            WRITE (out,100) gmaxs2(ng), 'gmaxs2',                      &
     &            'Maximum specific growth rate of diatom [1/day].'
#ifdef PHYTO_RESP
            WRITE (out,100) rrb1(ng), 'rrb1',                          &
     &            'basal respiration rate of small phyto. [1/day].'
            WRITE (out,100) rrg1(ng), 'rrg1',                          &
     &            'production-based resp. rate of small phyto. []'     
		WRITE (out,100) rrb2(ng), 'rrb2',                          &
     &            'basal respiration rate of diatom [1/day].'
		WRITE (out,100) rrg2(ng), 'rrg2',                          &
     &            'production-based resp. rate of diatom []'
#endif
            WRITE (out,100) beta1(ng), 'beta1',                        &
     &            'Microzooplankton maximum grazing rate [1/day].'
            WRITE (out,100) beta2(ng), 'beta2',                        &
     &            'Mesozooplankton maximum grazing rate [1/day].'
            WRITE (out,110) akz1(ng), 'akz1',                          &
     &            'Half saturation constant for microzooplankton',     &
     &            'grazing [mmol_N/m3].'
            WRITE (out,110) akz2(ng), 'akz2',                          &
     &            'Half saturation constant for mesozooplankton',      &
     &            'grazing [mmol_N/m3].'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                    &
     &            'Fraction of shortwave radiation that is',           &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) amaxs1(ng), 'amaxs1',                      &
     &            'Initial slope of P-I curve of small',               &
     &            'phytoplankton [1/(Watts/m2)/day].'
            WRITE (out,110) amaxs2(ng), 'amaxs2',                      &
     &            'Initial slope of P-I curve of large',               &
     &            'phytoplankton [1/(Watts/m2)/day].'
            WRITE (out,110) parsats1(ng), 'parsats1',                  &
     &            'PAR saturation onset parameter of',                 &
     &            'small phytoplankton [Watts/m2].'
            WRITE (out,110) parsats2(ng), 'parsats2',                  &
     &            'PAR saturation onset parameter of diatom',          &
     &            '[Watts/m2].'
            WRITE (out,110) pis1(ng), 'pis1',                          &
     &            'Ammonium inhibition parameter for small',           &
     &            'phytoplankton [mmol_N/m3].'
            WRITE (out,110) pis2(ng), 'pis2',                          &
     &            'Ammonium inhibition parameter for diatom',          &
     &            '[mmol_N/m3].'
            WRITE (out,110) akno3s1(ng), 'akno3s1',                    &
     &            'Half saturation concentration for nitrate',         &
     &            'uptake by small phytoplankton [mmol_N/m3].'
            WRITE (out,110) akno3s2(ng), 'akno3s2',                    &
     &            'Half saturation concentration for nitrate',         &
     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) aknh4s1(ng), 'aknh4s1',                    &
     &            'Half saturation concentration for ammonium',        &
     &            'uptake by small phytoplankton [mmol_N/m3].'
            WRITE (out,110) aknh4s2(ng), 'aknh4s2',                    &
     &            'Half saturation concentration for ammonium',        &
     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) akpo4s1(ng), 'akpo4s1',                    &
     &            'Half saturation concentration for phosphate',       &
     &            'uptake by small phytoplankton [mmol_P/m3].'
            WRITE (out,110) akpo4s2(ng), 'akpo4s2',                    &
     &            'Half saturation concentration for phosphate',       &
     &            'uptake by diatom [mmol_P/m3].'
            WRITE (out,110) akco2s1(ng), 'akco2s1',                    &
     &            'Half saturation concentration for co2',             &
     &            'uptake by small phytoplankton [mmol_C/m3].'
            WRITE (out,110) akco2s2(ng), 'akco2s2',                    &
     &            'Half saturation concentration for co2',             &
     &            'uptake by diatom [mmol_C/m3].'
            WRITE (out,110) aksio4s2(ng), 'aksio4s2',                  &
     &            'Half saturation constant for silicate',             &
     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) akox(ng), 'akox',                          &
     &            'Half saturation constant for oxidation',            &
     &            ' [mmol_O/m3].'
            WRITE (out,100) ak1(ng), 'ak1',                            &
     &            'Light attenuation coefficient of water [1/m].'
            WRITE (out,110) ak2(ng), 'ak2',                            &
     &            'Specific light attenuation coefficient for',        &
     &            'phytoplankton [1/m/(mmol_N/m3)].'
            WRITE (out,100) bgamma0(ng), 'bgamma0',                    &
     &            'Mesozooplankton specific mortality rate [1/day].'
            WRITE (out,110) bgamma1(ng), 'bgamma1',                    &
     &            'Grazing efficiency of microzooplankton',            &
     &            '[nondimensional].'
            WRITE (out,110) bgamma2(ng), 'bgamma2',                    &
     &            ' Grazing efficiency of mesozooplankton',            &
     &            '[nondimensional].'
            WRITE (out,100) bgamma3(ng), 'bgamma3',                    &
     &            'Death rate of small phytoplankton [1/day].'
            WRITE (out,100) bgamma4(ng), 'bgamma4',                    &
     &            'Death rate of large phytoplankton [1/day].'
            WRITE (out,100) bgamma5(ng), 'bgamma5',                    &
     &            'Decay rate of nitrogenous detritus [1/day].'
            WRITE (out,100) bgamma5s(ng), 'bgamma5s',                    &
     &            'Decay rate of silicious detritus [1/day].'
            WRITE (out,100) bgamma6(ng), 'bgamma6',                    &
     &            ' '
            WRITE (out,100) bgamma7(ng), 'bgamma7',                    &
     &            'Nitrafication rate [1/day].'
            WRITE (out,100) wsd(ng), 'wsd',                            &
     &            'Sinking velocity of detritus [m/day].'
            WRITE (out,100) wsdsi(ng), 'wsdsi',                        &
     &            'Sinking velocity of detritus silicate [m/day].'
            WRITE (out,100) wsp(ng), 'wsp',                            &
     &            'Sinking velocity of large phytoplankton [m/day].'
            WRITE (out,100) pco2a(ng), 'pco2a',                        &
     &            'Air pCO2 [ppmv].'
            WRITE (out,100) si2n(ng), 'si2n',                          &
     &            'Silicate to nitrogen ratio [mol_Si/mol_N].'
            WRITE (out,100) p2n(ng), 'p2n',                            &
     &            'Phosphorus to nitrogen ratio [mol_P/mol_N].'
            WRITE (out,100) o2no(ng), 'o2no',                          &
     &            'Oxygen to nitrate ratio [mol_O2/mol_NO3].'
            WRITE (out,100) o2nh(ng), 'o2nh',                          &
     &            'Oxygen to ammonium ratio [mol_O2/mol_NH4].'
            WRITE (out,100) c2n(ng), 'c2n',                            &
     &            'Carbon to nitrogen ratio [mol_C/mol_N].'
            WRITE (out,100) ro5(ng), 'ro5',                            &
     &            'Grazing preference for diatom [nondimensional].'
            WRITE (out,110) ro6(ng), 'ro6',                            &
     &            'Grazing preference for mesozooplankton',            &
     &            '[nondimensional].'
            WRITE (out,100) ro7(ng), 'ro7',                            &
     &            'Grazing preference for detritus [nondimensional].'
            WRITE (out,110) Chl2cs1_m(ng), 'Chl2cs1_m',                &
     &            'Maximum chlorophyll to carbon ratio for',           &
     &            'sphytoplankton [mg_Chl/mg_C)].'
            WRITE (out,110) Chl2cs2_m(ng), 'Chl2cs2_m',                &
     &            'Maximum chlorophyll to carbon ratio for',           &
     &            'Diatom [mg_Chl/mg_C)].'
#ifdef HAB
            WRITE (out,100) gmaxs3(ng), 'gmaxs3',                      &
     &            'Maximum specific growth rate of HAB phyto [1/day].'
            WRITE (out,110) amaxs3(ng), 'amaxs3',                      &
     &            'Initial slope of P-I curve of HAB',                 &
     &            'phytoplankton [1/(Watts/m2)/day].'
#  ifdef PHYTO_RESP
		WRITE (out,100) rrb3(ng), 'rrb3',                          &
     &            'basal respiration rate of HAB phyto. [1/day].'
		WRITE (out,100) rrg3(ng), 'rrg3',                          &
     &            'production-based resp. rate of HAB phyto. []'
#  endif
            WRITE (out,110) parsats3(ng), 'parsats3',                  &
     &            'PAR saturation onset parameter of HAB phyto',       &
     &            '[Watts/m2].'
            WRITE (out,110) pis3(ng), 'pis3',                          &
     &            'Ammonium inhibition parameter for HAB phyto',       &
     &            '[mmol_N/m3].'
            WRITE (out,110) akno3s3(ng), 'akno3s3',                    &
     &            'Half saturation concentration for nitrate',         &
     &            'uptake by HAB phyto [mmol_N/m3].'
            WRITE (out,110) aknh4s3(ng), 'aknh4s3',                    &
     &            'Half saturation concentration for ammonium',        &
     &            'uptake by HAB phyto [mmol_N/m3].'
            WRITE (out,110) akpo4s3(ng), 'akpo4s3',                    &
     &            'Half saturation concentration for phosphate',       &
     &            'uptake by HAB phyto [mmol_P/m3].'
            WRITE (out,110) akco2s3(ng), 'akco2s3',                    &
     &            'Half saturation concentration for co2',             &
     &            'uptake by HAB phyto [mmol_C/m3].'
            WRITE (out,110) aksio4s3(ng), 'aksio4s3',                  &
     &            'Half saturation constant for silicate',             &
     &            'uptake by HAB phyto [mmol_N/m3].'
            WRITE (out,100) bgamma4s3(ng), 'bgamma4s3',                    &
     &            'Death rate of HAB phytoplankton [1/day].'
            WRITE (out,100) wsp3(ng), 'wsp3',                            &
     &            'Sinking velocity of HAB phytoplankton [m/day].'
            WRITE (out,100) ro5H(ng), 'ro5H',                            &
     &            'Grazing preference for HAB phyto [nondimensional].'
            WRITE (out,110) Chl2cs3_m(ng), 'Chl2cs3_m',                &
     &            'Maximum chlorophyll to carbon ratio for',           &
     &            'HAB phyto [mg_Chl/mg_C)].'
#endif
		WRITE (out,110) q10_phyto_prod(ng), 'q10_phyto_prod',      &
     &            'Q10 for phytoplankton production [nondimensional]'   
		WRITE (out,110) q10_phyto_resp(ng), 'q10_phyto_resp',      &
     &            'Q10 for phytoplankton respiration [nondimensional]'  
		WRITE (out,110) q10_bact_remin(ng), 'q10_bact_remin',      &
     &            'Q10 for bacterial remineralization [nondimensional]' 
		WRITE (out,110) q10_opal_dissol(ng), 'q10_opal_dissol',    &
     &            'Q10 for opal dissolution [nondimensional]'  
		WRITE (out,110) q10_zoopl_resp(ng), 'q10_zoopl_resp',      &
     &            'Q10 for zooplankton respiration [nondimensional]' 
		WRITE (out,110) q10_zoopl_graz(ng), 'q10_zoopl_graz',      &
     &            'Q10 for zooplankton grazing [nondimensional]' 
		WRITE (out,110) q10_nitrif(ng), 'q10_nitrif',      	     &
     &            'Q10 for nitrification [nondimensional]' 

#ifdef SEDBIO
            fname=bsedparnam
            IF (.not.find_file(ng, fname, 'BSEDPOSNAM')) THEN
                IF (MASTER) WRITE(out,270) TRIM(fname)
                exit_flag=4
                RETURN
            ELSE
                WRITE (out,230) ' Sediment Biology Parameters File: ', &
     &                    TRIM(fname)
            END IF
#endif
#ifdef OPTICS_OP1
            WRITE (out,110) optic_upd_fac(ng), 'optic_upd_fac',        &
     &            'Compute kdPAR every optic_upd_fac timesteps',       &
     &            '[no dimension].'
#endif

#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu2(i,ng), 'tnu2', i,                 &
     &              'NLM Horizontal, harmonic mixing coefficient',     &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu2(i,ng), 'ad_tnu2', i,              &
     &              'ADM Horizontal, harmonic mixing coefficient',     &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu2(i,ng), 'tl_tnu2', i,              &
     &              'TLM Horizontal, harmonic mixing coefficient',     &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu4(i,ng), 'tnu4', i,                 &
     &              'NLM Horizontal, biharmonic mixing coefficient',   &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu4(i,ng), 'ad_tnu4', i,              &
     &              'ADM Horizontal, biharmonic mixing coefficient',   &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu4(i,ng), 'tl_tnu4', i,              &
     &              'TLM Horizontal, biharmonic mixing coefficient',   &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,               &
     &             'Background vertical mixing coefficient (m2/s)',    &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef FORWARD_MIXING
            DO itrc=1,NBT
              i=idbio(itrc)
# ifdef ADJOINT
              WRITE (out,90) ad_Akt_fac(i,ng), 'ad_Akt_fac', i,        &
     &              'ADM basic state vertical mixing scale factor',    &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_Akt_fac(i,ng), 'tl_Akt_fac', i,        &
     &              'TLM basic state vertical mixing scale factor',    &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                  &
     &              'Nudging/relaxation time scale (days)',            &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',        &
     &              i, 'Turning ON  point sources/Sink on tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',        &
     &              i, 'Turning OFF point sources/Sink on tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
#ifdef TCLIMATOLOGY
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM',          &
     &              i, 'Processing climatology on tracer ', i,         &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                   &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                  &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                   &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                  &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef OPTICS_OP1
            IF (Hout(idkdPAR,ng)) WRITE (out,65) Hout(idkdPAR,ng),         &
     &       'Hout(idkdPAR)',                                            &
     &       'Write out light attenuation coefficient ',TRIM(Vname(1,idkdPAR))
#endif
! AKB 6/22/23
#ifdef CARBON
            IF (Hout(idfgCO2,ng)) WRITE (out,65) Hout(idfgCO2,ng),         &
     &       'Hout(idfgCO2)',                                            &
     &       'Write out CO2 downward surface flux ',TRIM(Vname(1,idfgCO2))
            IF (Hout(idpCO2s,ng)) WRITE (out,65) Hout(idpCO2s,ng),         &
     &       'Hout(idpCO2s)',                                            &
     &       'Write out surface partial pressure of CO2. ',TRIM(Vname(1,idpCO2s))
#endif

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            WRITE (out,'(1x)')
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTvar(i),ng)) WRITE (out,60)                   &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                  &
     &            'Write out averaged tracer ', i,                     &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTTav(i),ng)) WRITE (out,60)                   &
     &            Aout(idTTav(i),ng), 'Aout(idTTav)',                  &
     &            'Write out averaged <t*t> for tracer ', i,           &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idUTav(i),ng)) WRITE (out,60)                   &
     &            Aout(idUTav(i),ng), 'Aout(idUTav)',                  &
     &            'Write out averaged <u*t> for tracer ', i,           &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idVTav(i),ng)) WRITE (out,60)                   &
     &            Aout(idVTav(i),ng), 'Aout(idVTav)',                  &
     &            'Write out averaged <v*t> for tracer ', i,           &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHUTav(i),ng)) WRITE (out,60)                   &
     &            Aout(iHUTav(i),ng), 'Aout(iHUTav)',                  &
     &            'Write out averaged <Huon*t> for tracer ', i,        &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHVTav(i),ng)) WRITE (out,60)                   &
     &            Aout(iHVTav(i),ng), 'Aout(iHVTav)',                  &
     &            'Write out averaged <Hvom*t> for tracer ', i,        &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
#endif
#ifdef DIAGNOSTICS_TS
            WRITE (out,'(1x)')
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTrate),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTrate)',                 &
     &              'Write out rate of change of tracer ', itrc,       &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThadv),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iThadv)',                 &
     &              'Write out horizontal advection, tracer ', itrc,   &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTxadv),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTxadv)',                 &
     &              'Write out horizontal X-advection, tracer ', itrc, &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTyadv),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTyadv)',                 &
     &              'Write out horizontal Y-advection, tracer ', itrc, &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvadv),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTvadv)',                 &
     &              'Write out vertical advection, tracer ', itrc,     &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
# if defined TS_DIF2 || defined TS_DIF4
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThdif),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iThdif)',                 &
     &              'Write out horizontal diffusion, tracer ', itrc,   &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(i,iTxdif),ng))                           &
     &          WRITE (out,60) .TRUE., 'Dout(iTxdif)',                 &
     &              'Write out horizontal X-diffusion, tracer ', itrc, &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTydif),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTydif)',                 &
     &              'Write out horizontal Y-diffusion, tracer ', itrc, &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTsdif),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTsdif)',                 &
     &              'Write out horizontal S-diffusion, tracer ', itrc, &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
#  endif
# endif
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvdif),ng))                        &
     &          WRITE (out,60) .TRUE., 'Dout(iTvdif)',                 &
     &              'Write out vertical diffusion, tracer ', itrc,     &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
#endif
#ifdef DIAGNOSTICS_BIO
            IF (NDbio2d.gt.0) THEN
              DO itrc=1,NDbio2d
                i=iDbio2(itrc)
                IF (Dout(i,ng)) WRITE (out,140)                        &
     &              Dout(i,ng), 'Hout(iDbio2)',                        &
     &              'Write out diagnostics for', TRIM(Vname(1,i))
              END DO
            END IF
            DO itrc=1,NDbio3d
              i=iDbio3(itrc)
              IF (Dout(i,ng)) WRITE (out,140)                          &
     &            Dout(i,ng), 'Hout(iDbio3)',                          &
     &            'Write out diagnostics for', TRIM(Vname(1,i))
            END DO
#endif

          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          nl_tnu4(i,ng)=SQRT(ABS(nl_tnu4(i,ng)))
#ifdef ADJOINT
          ad_tnu4(i,ng)=SQRT(ABS(ad_tnu4(i,ng)))
#endif
#if defined TANGENT || defined TL_IOMS
          tl_tnu4(i,ng)=SQRT(ABS(tl_tnu4(i,ng)))
#endif
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' UMaine CoSiNE Model Parameters, Grid: ',i2.2,      &
     &        /,  ' =================================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  65  FORMAT (10x,l1,2x,a,t30,a,':',1x,a)
  70  FORMAT (f11.3,2x,a,t30,a)
  80  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',       &
     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t30,a,1x,a) 
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 230  FORMAT (2x,a,a)
 270  FORMAT (/,' READ_BIOPAR - could not find input file:  ',a)

      RETURN
      END SUBROUTINE read_BioPar
