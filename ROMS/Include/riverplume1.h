/*
** svn $Id: riverplume1.h 778 2015-08-22 01:20:04Z arango $
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for River Plume Test (original version).
**
** Application flag:   RIVERPLUME1
** Input script:       ocean_riverplume1.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define AVERAGES

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#define ANA_GRID
#define ANA_INITIAL
#define ANA_PSOURCE
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
