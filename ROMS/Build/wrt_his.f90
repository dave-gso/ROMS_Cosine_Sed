      SUBROUTINE wrt_his (ng)
!
!svn $Id: wrt_his.F 764 2015-06-23 22:40:38Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested model fields at requested levels      !
!  into history NetCDF file.                                           !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod, ONLY : scale_omega
      USE uv_rotate_mod, ONLY : uv_rotate2d
      USE uv_rotate_mod, ONLY : uv_rotate3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, status, tile
      integer :: i, itrc, j, k
      real(r8) :: scale
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
!
      SourceFile='wrt_his.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%Fcount
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, HIS(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/HIS(ng)%Rindex/), (/1/),                  &
     &                      ncid = HIS(ng)%ncid,                        &
     &                      varid = HIS(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idFsur), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % zeta(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbar), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % ubar(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbar), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % vbar(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % ubar(:,:,kstp(ng)),               &
     &                    OCEAN(ng) % vbar(:,:,kstp(ng)),               &
     &                    Ur2d, Vr2d)
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu2dE), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv2dN), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % u(:,:,:,nrhs(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % v(:,:,:,nrhs(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % u(:,:,:,nrhs(ng)),                &
     &                    OCEAN(ng) % v(:,:,:,nrhs(ng)),                &
     &                    Ur3d, Vr3d)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu3dE), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     Ur3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv3dN), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     Vr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_r8
        gtype=gfactor*w3dvar
        tile=MyRank
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idOvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     Wr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % wvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Tid(itrc), &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idDano), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % rho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvis), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTdif), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSdif), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idMtke), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % tke(:,:,:,nrhs(ng)),            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idMtls), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % gls(:,:,:,nrhs(ng)),            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPair), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Pair)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUair), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Uwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Hout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVair), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Vwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface active traces fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                    &
     &                       HIS(ng)%Vid(idTsur(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLhea), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % lhflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idShea), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % shflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLrad), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % lrflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idEmPf), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % EminusP)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idevap), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % evap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idrain), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % rain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSrad), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % srflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUsms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVsms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % bustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % bvstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HIS(ng)%name, HIS(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) kstp(ng), nrhs(ng), HIS(ng)%Rindex
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
      RETURN
      END SUBROUTINE wrt_his
