      MODULE nf_fwrite3d_mod
!
!svn $Id: nf_fwrite3d.F 751 2015-01-07 22:56:36Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 3D array into an  !
!  output NetCDF file.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvarid      NetCDF variable ID.                                 !
!     tindex       NetCDF time record index to write.                  !
!     gtype        Grid type. If negative, only write water points.    !
!     LBi          I-dimension Lower bound.                            !
!     UBi          I-dimension Upper bound.                            !
!     LBj          J-dimension Lower bound.                            !
!     UBj          J-dimension Upper bound.                            !
!     LBk          K-dimension Lower bound.                            !
!     UBk          K-dimension Upper bound.                            !
!     Amask        land/Sea mask, if any (real).                       !
!     Ascl         Factor to scale field before writing (real).        !
!     A            Field to write out (real).                          !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional).                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite3d  Error flag (integer).                               !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      FUNCTION nf_fwrite3d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, LBk, UBk, Ascl,         &
     &                      Amask,                                      &
     &                      A,  SetFillVal)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_bcasti
      USE distribute_mod, ONLY : mp_gather3d
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      real(r8), intent(in) :: Ascl
      real(r8), intent(in) :: Amask(LBi:,LBj:)
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: LandFill
      integer :: i, j, k, ic, Npts
      integer :: Imin, Imax, Jmin, Jmax, Koff
      integer :: Ilen, Jlen, Klen, IJlen, MyType, status
      integer, dimension(4) :: start, total
      integer :: nf_fwrite3d
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Awrk
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      MyType=gtype
      SELECT CASE (ABS(MyType))
        CASE (p3dvar)
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
        CASE (r3dvar)
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
        CASE (u3dvar)
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
        CASE (v3dvar)
          Imin=IOBOUNDS(ng)%ILB_v
          Imax=IOBOUNDS(ng)%IUB_v
          Jmin=IOBOUNDS(ng)%JLB_v
          Jmax=IOBOUNDS(ng)%JUB_v
        CASE DEFAULT
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
      END SELECT
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
!
!  Set switch to replace land areas with fill value, spval.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Awrk=0.0_r8
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order.
!  Overwrite masked points with special value.
!-----------------------------------------------------------------------
!
        CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                    tindex, gtype, Ascl,                          &
     &                    Amask,                                        &
     &                    A, Npts, Awrk, SetFillVal)
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
        nf_fwrite3d=nf90_noerr
        IF (OutThread) THEN
          IF (gtype.gt.0) THEN
            start(1)=1
            total(1)=Ilen
            start(2)=1
            total(2)=Jlen
            start(3)=1
            total(3)=Klen
            start(4)=tindex
            total(4)=1
          ELSE
            start(1)=1
            total(1)=Npts
            start(2)=tindex
            total(2)=1
          END IF
          status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
          nf_fwrite3d=status
        END IF
!
!-----------------------------------------------------------------------
!  Broadcast IO error flag to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcasti (ng, model, status)
      nf_fwrite3d=status
      RETURN
      END FUNCTION nf_fwrite3d
      END MODULE nf_fwrite3d_mod
