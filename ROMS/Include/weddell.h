/*
** svn $Id: weddell.h 778 2015-08-22 01:20:04Z arango $
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Idealized Weddell Sea Application: Tides and Ice Shelf Test.
**
** Application flag:   WEDDELL
** Input script:       ocean_wedell.in
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#undef  UV_VIS4
#undef  MIX_S_UV
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define ICESHELF
#define AVERAGES
#define RADIATION_2D
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_SRFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
