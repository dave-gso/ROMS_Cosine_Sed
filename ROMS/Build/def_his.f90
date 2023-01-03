      SUBROUTINE def_his (ng, ldef)
!
!svn $Id: def_his.F 764 2015-06-23 22:40:38Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates history NetCDF file, it defines its            !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
      integer, parameter :: Natt = 25
      integer :: i, j, ifield, itrc, nvd3, nvd4, varid
      integer :: recdim, status
      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)
      integer :: def_dim
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
      real(r8) :: Aval(6)
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_his.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=HIS(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new history file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), HIS(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'N',             &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, HIS(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
!
!  Define dimension vector for staggered w-momentum type variables.
!
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        HIS(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, HIS(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idFsur),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D U-momentum component.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbar),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D V-momentum component.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbar),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Hout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu2dE),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Hout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv2dN),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D U-momentum component.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUvel),   &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D V-momentum component.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvel),   &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Hout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu3dE),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Hout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv3dN),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWvel),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idOvel),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Tid(itrc),   &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idDano),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvis),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTdif),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSdif),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idMtke),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idMtls),   &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface air pressure.
!
        IF (Hout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface winds.
!
        IF (Hout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (Hout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            status=def_var(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idTsur(itrc)), NF_FOUT,          &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLhea),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idShea),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLrad),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define E-P flux (as computed by bulk_flux.F).
!
        IF (Hout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idEmPf),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idevap),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idrain),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSrad),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUsms),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVsms),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbms),   &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbms),   &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, HIS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, HIS(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing history file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=HIS(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open history file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, HIS(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  history variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            HIS(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            HIS(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            HIS(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            HIS(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            HIS(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            HIS(ng)%Vid(idv2dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            HIS(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            HIS(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            HIS(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            HIS(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            HIS(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            HIS(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            HIS(ng)%Vid(idDano)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            HIS(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            HIS(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            HIS(ng)%Vid(idSdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            HIS(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            HIS(ng)%Vid(idMtls)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            HIS(ng)%Vid(idPair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            HIS(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            HIS(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            HIS(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            HIS(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            HIS(ng)%Vid(idLrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            HIS(ng)%Vid(idTair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            HIS(ng)%Vid(idEmPf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            HIS(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            HIS(ng)%Vid(idrain)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            HIS(ng)%Vid(idSrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            HIS(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            HIS(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            HIS(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            HIS(ng)%Vid(idVbms)=var_id(i)
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              HIS(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              HIS(ng)%Vid(idTsur(itrc))=var_id(i)
            END IF
          END DO
        END DO
!
!  Check if history variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Hout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Hout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Hout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Hout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPair).and.Hout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idEmPf).and.Hout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefHIS(ng).gt.0) THEN
          HIS(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefHIS(ng)*((ntstart(ng)-1)/ndefHIS(ng)))/   &
     &                   nHIS(ng)
        ELSE
          HIS(ng)%Rindex=(ntstart(ng)-1)/nHIS(ng)
        END IF
        HIS(ng)%Rindex=MIN(HIS(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (6x,'DEF_HIS   - creating history file, Grid ',i2.2,': ',  &
     &        a)
  20  FORMAT (6x,'DEF_HIS   - inquiring history file, Grid ',i2.2,': ', &
     &        a)
  30  FORMAT (/,' DEF_HIS - unable to create history NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_HIS - unable to open history NetCDF file: ',a)
  60  FORMAT (/,' DEF_HIS - unable to find variable: ',a,2x,            &
     &        ' in history NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_his
