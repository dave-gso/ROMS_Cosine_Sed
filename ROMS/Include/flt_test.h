/*
** svn $Id: flt_test.h 779 2015-08-22 03:11:58Z arango $
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Floats Tracking Test.
**
** Application flag:   FLT_TEST
** Input script:       ocean_flt_test2d.in,  ocean_flt_test3d.in
**                     floats_flt_test2d.in, floats_flt_test3d.in
*/

#define UV_ADV
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define FLOATS
#define MASKING
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

#ifdef SOLVE3D
# define DJ_GRADPS
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define TS_A4HADVECTION
# define TS_A4VADVECTION
# define BODYFORCE
# define ANA_BTFLUX
# define ANA_STFLUX
#endif

