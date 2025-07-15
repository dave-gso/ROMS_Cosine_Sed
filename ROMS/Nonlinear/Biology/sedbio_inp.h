      SUBROUTINE read_SedBioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in sediment biological model input parameters.   !
!                                                                      !
!  By: Dave Ullman 08/2021                                                !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_sedbio
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

! Changed dimensions of Ltrc from 5 to NSF since NSF is 7 and larger than 5 (number of sediment fluxes)
! The robust way to do this would be to have the dimension be the largest of NSF, NPWC, NPOM, N-decay rate
! AKB 5/18/23
!
      logical, dimension(NSF,Ngrids) :: Ltrc
      real(r8), dimension(NSF,Ngrids) :: Rbio
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
      iTrcStr=1                          ! first SEDBIOL tracer to process
      iTrcEnd=NBBT                       ! last SEDBIOL tracer to process
      nline=0                            ! LBC multi-line counter
!
!
!-----------------------------------------------------------------------
!  Read in UMaine CoSiNE sediment biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
        SELECT CASE (TRIM(KeyWord))
          CASE('Lsedbiology')
             Npts=load_l(Nval, Cval, Ngrids, Lsedbiology)
          CASE('bUmax')
             Npts=load_r(Nval, Rval, nspc*Ngrids, Rbio)
             DO ng=1,Ngrids
                DO itrc=1,nspc
                    bUmax(ng,itrc)=Rbio(itrc,ng)
                END DO
             END DO
          CASE('bUmaxSi')
            Npts=load_r(Nval, Rval, nspc*Ngrids, Rbio)
            DO ng=1,Ngrids
                DO itrc=1,nspc
                    bUmaxSi(ng,itrc)=Rbio(itrc,ng)
                END DO
            END DO
          CASE('bdep')
            Npts=load_r(Nval, Rval, Ngrids, bdep)
          CASE('balpha')
            Npts=load_r(Nval, Rval, Ngrids, balpha)
          CASE('bw')
            Npts=load_r(Nval, Rval, Ngrids, bw)
          CASE('btheta_diag')
            Npts=load_r(Nval, Rval, Ngrids, btheta_diag)
          CASE('bnit')
            Npts=load_r(Nval, Rval, Ngrids, bnit)
          CASE('btheta_nit')
            Npts=load_r(Nval, Rval, Ngrids, btheta_nit)
          CASE('bdo_nit')
            Npts=load_r(Nval, Rval, Ngrids, bdo_nit)
          CASE('bdenit')
            Npts=load_r(Nval, Rval, Ngrids, bdenit)
          CASE('btheta_denit')
            Npts=load_r(Nval, Rval, Ngrids, btheta_denit)
          CASE('bpsi_n')
            Npts=load_r(Nval, Rval, Ngrids, bpsi_n)
          CASE('bdo_c')
            Npts=load_r(Nval, Rval, Ngrids, bdo_c)
          CASE('bao2')
            Npts=load_r(Nval, Rval, Ngrids, bao2)
          CASE('bpi')
            Npts=load_r(Nval, Rval, Ngrids, bpi)
          CASE('bpsi_p')
            Npts=load_r(Nval, Rval, Ngrids, bpsi_p)
          CASE('bfc')
            Npts=load_r(Nval, Rval, nspc*Ngrids, Rbio)
            DO ng=1,Ngrids
                DO itrc=1,nspc
                    bfc(ng,itrc)=Rbio(itrc,ng)
                END DO
            END DO
          CASE('bfn')
            Npts=load_r(Nval, Rval, nspc*Ngrids, Rbio)
            DO ng=1,Ngrids
                DO itrc=1,nspc
                    bfn(ng,itrc)=Rbio(itrc,ng)
                END DO
            END DO
          CASE('bfp')
            Npts=load_r(Nval, Rval, nspc*Ngrids, Rbio)
            DO ng=1,Ngrids
                DO itrc=1,nspc
                    bfp(ng,itrc)=Rbio(itrc,ng)
                END DO
            END DO
          CASE('bfs')
            Npts=load_r(Nval, Rval, nspc*Ngrids, Rbio)
            DO ng=1,Ngrids
                DO itrc=1,nspc
                    bfs(ng,itrc)=Rbio(itrc,ng)
                END DO
            END DO
          CASE ('Hout(idbPM)')
            Npts=load_l(Nval, Cval, NPOM*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NPOM
                i=idbPM(itrc)
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                         &
     &                                'idbPM(', itrc, ')'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          CASE ('Hout(idbPW)')
            Npts=load_l(Nval, Cval, NPWC*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NPWC
                i=idbPW(itrc)
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                         &
     &                                'idbPW(', itrc, ')'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          CASE ('Hout(idbDR)')
            Npts=load_l(Nval, Cval, NDR*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NDR
                i=idbDR(itrc)
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                         &
     &                                'idbDR(', itrc, ')'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          CASE ('Hout(idbSF)')
            Npts=load_l(Nval, Cval, NSF*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NSF
                i=idbSF(itrc)
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                         &
     &                                'idbSF(', itrc, ')'
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
          IF (Lsedbiology(ng)) THEN
            WRITE (out,40) ng
            DO itrc=1,nspc
              WRITE (out,90) bUmax(ng,itrc), 'bUmax', itrc,                  &
     &              'Maximum particulate organic matter decay rate',   &
     &              'for POM class [1/day] '
            ENDDO
            DO itrc=1,nspc
              WRITE (out,90) bUmaxSi(ng,itrc), 'bUmaxSi', itrc,              &
     &              'Maximum particulate silicate decay rate',         &
     &              'for POM class [1/day] '
            ENDDO
            WRITE (out,100) bdep(ng), 'bdep',                          &
     &            'Thickness of benthic sediment layer (m).'
            WRITE (out,110) balpha(ng), 'balpha',                      &
     &            'Ratio of bacterial growth rate to max. decay rate', &
     &            '[nondimensional].'
            WRITE (out,100) bw(ng), 'bw',                              &
     &            'Sedimentation rate [cm/year].'
            WRITE (out,110) btheta_diag(ng), 'btheta_diag',            &
     &           'Base for temperature dependence of diagenesis',      &
     &            '[nondimensional].'
            WRITE (out,110) bnit(ng), 'bnit',                          &
     &            'Nitrification rate at reference temperature and',   &
     &             'unlimited O2 [1/day]'
            WRITE (out,110) btheta_nit(ng), 'btheta_nit',               &
     &            'Base for temperature dependence of nitrification',   &
     &            '[nondimensional]'
            WRITE (out,110) bdo_nit(ng), 'bdo_nit',                    &
     &            'Half saturation constant for O2 dependence of',      &
     &            'nitrification [mg/L]'
            WRITE (out,110) bdenit(ng), 'bdenit',                      &
     &            'Denitrification rate at reference temperature and', &
     &            'anoxic conditions [1/day].'
            WRITE (out,110) btheta_denit(ng), 'btheta_denit',          &
     &            'Base for temperature dependence of denitrification', &
     &            '[nondimensional].'
            WRITE (out,110) bpsi_n(ng), 'bpsi_n',                      &
     &            'Constant for oxygen dependence of denitrification', &
     &            '[mg/L]'
            WRITE (out,110) bdo_c(ng), 'bdo_c',                        &
     &            'Half saturation constant for carbon oxidation',     &
     &            '[mg/L].'
            WRITE (out,110) bao2(ng), 'bao2',                          &
     &            'Stoichiometric coefficient for oxygen consumed by', &
     &            'nitrification [nondimensional].'
            WRITE (out,110) bpi(ng), 'bpi',                            &
     &            'fraction of PO4 dissolved in pore water under',     &
     &            'anoxic conditions [nondimensional].'
            WRITE (out,110) bpsi_p(ng), 'bpsi_p',                      &
     &            'Constant for oxygen dependence of dissolved PO4',   &
     &            '[mg/L].'
      
            DO itrc=1,nspc
              WRITE (out,90) bfc(ng,itrc), 'bfc', itrc,                      &
     &              'Particulate carbon partitioning coefficient',     &
     &              '[nondimensional]'
            ENDDO
            DO itrc=1,nspc
              WRITE (out,90) bfn(ng,itrc), 'bfn', itrc,                      &
     &              'Particulate nitrogen partitioning coefficient',   &
     &              '[nondimensional]'
            ENDDO
            DO itrc=1,nspc
              WRITE (out,90) bfp(ng,itrc), 'bfp', itrc,                      &
     &              'Particulate phosphorous partitioning coefficient',&
     &              '[nondimensional]'
            ENDDO
            DO itrc=1,nspc
              WRITE (out,90) bfs(ng,itrc), 'bfs', itrc,                      &
&              'Particulate silicon partitioning coefficient',         &
&              '[nondimensional]'
            ENDDO

          END IF
        END DO
      END IF
!

  30  FORMAT (/,' read_SedBioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' UMaine CoSiNE Sediment Model Parameters, Grid: ',  &
     &         i2.2,/,  ' =================================',/)
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t30,a)
  80  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
!90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',       &
     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t30,a,1x,a) 
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_SedBioPar
