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
      integer :: Npts, Nval, i, itrc, ng, status
      integer :: igrid, itracer,iTrcStr, iTrcEnd,nline,ifield
      integer :: decode_line, load_i, load_l, load_r,load_lbc
      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc
      real(r8), dimension(NBT,Ngrids) :: Rbio
      real(r8), dimension(100) :: Rval
      character (len=40 ) :: KeyWord
      character (len=256) :: line
!     character (len=256), dimension(100) :: Cval
      character (len=256), dimension(200) :: Cval
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
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &

     &                      iTrcStr, iTrcEnd,                           &

     &                      Vname(1,idTvar(idbio(itracer))), LBC)
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
                    IF (Master) WRITE (out,30)                          &

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
                    IF (Master) WRITE (out,30)                          &

     &                                'idTsur(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
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
            WRITE (out,50) BioIter(ng), 'BioIter',                      &

     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,110) reg1(ng), 'reg1',                           &

     &            'Microzooplankton excretion rate to ammonium',        &

     &            '[1/day].'
            WRITE (out,100) reg2(ng), 'reg2',                           &

     &            'Mesozooplankton excretion rate to ammonium [1/day].'
            WRITE (out,110) gmaxs1(ng), 'gmaxs1',                       &

     &           'Maximum specific growth rate of small phytoplankton', &

     &            '[1/day].'
            WRITE (out,100) gmaxs2(ng), 'gmaxs2',                       &

     &            'Maximum specific growth rate of diatom [1/day].'
            WRITE (out,100) beta1(ng), 'beta1',                         &

     &            'Microzooplankton maximum grazing rate [1/day].'
            WRITE (out,100) beta2(ng), 'beta2',                         &

     &            'Mesozooplankton maximum grazing rate [1/day].'
            WRITE (out,110) akz1(ng), 'akz1',                           &

     &            'Half saturation constant for microzooplankton',      &

     &            'grazing [mmol_N/m3].'
            WRITE (out,110) akz2(ng), 'akz2',                           &

     &            'Half saturation constant for mesozooplankton',       &

     &            'grazing [mmol_N/m3].'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &

     &            'Fraction of shortwave radiation that is',            &

     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) amaxs1(ng), 'amaxs1',                       &

     &            'Initial slope of P-I curve of small',                &

     &            'phytoplankton [1/(Watts/m2)/day].'
            WRITE (out,110) amaxs2(ng), 'amaxs2',                       &

     &            'Initial slope of P-I curve of large',                &

     &            'phytoplankton [1/(Watts/m2)/day].'
            WRITE (out,110) parsats1(ng), 'parsats1',                   &

     &            'PAR saturation onset parameter of',                  &

     &            'small phytoplankton [Watts/m2].'
            WRITE (out,110) parsats2(ng), 'parsats2',                   &

     &            'PAR saturation onset parameter of diatom',           &

     &            '[Watts/m2].'
            WRITE (out,110) pis1(ng), 'pis1',                           &

     &            'Ammonium inhibition parameter for small',            &

     &            'phytoplankton [mmol_N/m3].'
            WRITE (out,110) pis2(ng), 'pis2',                           &

     &            'Ammonium inhibition parameter for diatom',           &

     &            '[mmol_N/m3].'
            WRITE (out,110) akno3s1(ng), 'akno3s1',                     &

     &            'Half saturation concentration for nitrate',          &

     &            'uptake by small phytoplankton [mmol_N/m3].'
            WRITE (out,110) akno3s2(ng), 'akno3s2',                     &

     &            'Half saturation concentration for nitrate',          &

     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) aknh4s1(ng), 'aknh4s1',                     &

     &            'Half saturation concentration for ammonium',         &

     &            'uptake by small phytoplankton [mmol_N/m3].'
            WRITE (out,110) aknh4s2(ng), 'aknh4s2',                     &

     &            'Half saturation concentration for ammonium',         &

     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) akpo4s1(ng), 'akpo4s1',                     &

     &            'Half saturation concentration for phosphate',        &

     &            'uptake by small phytoplankton [mmol_P/m3].'
            WRITE (out,110) akpo4s2(ng), 'akpo4s2',                     &

     &            'Half saturation concentration for phosphate',        &

     &            'uptake by diatom [mmol_P/m3].'
            WRITE (out,110) akco2s1(ng), 'akco2s1',                     &

     &            'Half saturation concentration for co2',              &

     &            'uptake by small phytoplankton [mmol_C/m3].'
            WRITE (out,110) akco2s2(ng), 'akco2s2',                     &

     &            'Half saturation concentration for co2',              &

     &            'uptake by diatom [mmol_C/m3].'
            WRITE (out,110) aksio4s2(ng), 'aksio4s2',                   &

     &            'Half saturation constant for silicate',              &

     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) akox(ng), 'akox',                           &

     &            'Half saturation constant for oxidation',             &

     &            ' [mmol_O/m3].'
            WRITE (out,100) ak1(ng), 'ak1',                             &

     &            'Light attenuation coefficient of water [1/m].'
            WRITE (out,110) ak2(ng), 'ak2',                             &

     &            'Specific light attenuation coefficient for',         &

     &            'phytoplankton [1/m/(mmol_N/m3)].'
            WRITE (out,100) bgamma0(ng), 'bgamma0',                     &

     &            'Mesozooplankton specific mortality rate [1/day].'
            WRITE (out,110) bgamma1(ng), 'bgamma1',                     &

     &            'Grazing efficiency of microzooplankton',             &

     &            '[nondimensional].'
            WRITE (out,110) bgamma2(ng), 'bgamma2',                     &

     &            ' Grazing efficiency of mesozooplankton',             &

     &            '[nondimensional].'
            WRITE (out,100) bgamma3(ng), 'bgamma3',                     &

     &            'Death rate of small phytoplankton [1/day].'
            WRITE (out,100) bgamma4(ng), 'bgamma4',                     &

     &            'Death rate of large phytoplankton [1/day].'
            WRITE (out,100) bgamma5(ng), 'bgamma5',                     &

     &            'Decay rate of detritus [1/day].'
            WRITE (out,100) bgamma6(ng), 'bgamma6',                     &

     &            ' '
            WRITE (out,100) bgamma7(ng), 'bgamma7',                     &

     &            'Nitrafication rate [1/day].'
            WRITE (out,100) wsd(ng), 'wsd',                             &

     &            'Sinking velocity of detritus [m/day].'
            WRITE (out,100) wsdsi(ng), 'wsdsi',                         &

     &            'Sinking velocity of detritus silicate [m/day].'
            WRITE (out,100) wsp(ng), 'wsp',                             &

     &            'Sinking velocity of large phytoplankton [m/day].'
            WRITE (out,100) pco2a(ng), 'pco2a',                         &

     &            'Air pCO2 [ppmv].'
            WRITE (out,100) si2n(ng), 'si2n',                           &

     &            'Silicate to nitrogen ratio [mol_Si/mol_N].'
            WRITE (out,100) p2n(ng), 'p2n',                             &

     &            'Phosphorus to nitrogen ratio [mol_P/mol_N].'
            WRITE (out,100) o2no(ng), 'o2no',                           &

     &            'Oxygen to nitrate ratio [mol_O2/mol_NO3].'
            WRITE (out,100) o2nh(ng), 'o2nh',                           &

     &            'Oxygen to ammonium ratio [mol_O2/mol_NH4].'
            WRITE (out,100) c2n(ng), 'c2n',                             &

     &            'Carbon to nitrogen ratio [mol_C/mol_N].'
            WRITE (out,100) ro5(ng), 'ro5',                             &

     &            'Grazing preference for diatom [nondimensional].'
            WRITE (out,110) ro6(ng), 'ro6',                             &

     &            'Grazing preference for mesozooplankton',             &

     &            '[nondimensional].'
            WRITE (out,100) ro7(ng), 'ro7',                             &

     &            'Grazing preference for detritus [nondimensional].'
            WRITE (out,110) Chl2cs1_m(ng), 'Chl2cs1_m',                 &

     &            'Maximum chlorophyll to carbon ratio for',            &

     &            'sphytoplankton [mg_Chl/mg_C)].'
            WRITE (out,110) Chl2cs2_m(ng), 'Chl2cs1_m',                 &

     &            'Maximum chlorophyll to carbon ratio for',            &

     &            'Diatom [mg_Chl/mg_C)].'               
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu2(i,ng), 'tnu2', i,                  &

     &              'NLM Horizontal, harmonic mixing coefficient',      &

     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &

     &             'Background vertical mixing coefficient (m2/s)',     &

     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &

     &              'Nudging/relaxation time scale (days)',             &

     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',         &

     &              i, 'Turning ON  point sources/Sink on tracer ', i,  &

     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',         &

     &              i, 'Turning OFF point sources/Sink on tracer ', i,  &

     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &

     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &

     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &

     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &

     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
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
  40  FORMAT (/,/,' UMaine CoSiNE Model Parameters, Grid: ',i2.2,       &

     &        /,  ' =================================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t30,a)
  80  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &

     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t30,a,1x,a) 
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
      RETURN
      END SUBROUTINE read_BioPar
