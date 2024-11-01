/*
** svn $Id: channel_neck.h 779 2015-08-22 03:11:58Z arango $
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Channel with a Constriction Test.
**
** Application flag:   CHANNEL_NECK
** Input script:       ocean_channel_neck1.in, ocean_channel_neck3.in
*/

#undef  AD_SENSITIVITY
#define CORRELATION
#undef  MULTI_DOMAIN
#undef  OPT_OBSERVATIONS
#undef  OPT_PERTURBATION
#undef  SANITY_CHECK
#undef  TLM_DRIVER

#define UV_ADV
#define UV_COR
#define UV_LDRAG
#define UV_VIS2
#define MIX_GEO_UV
#define DJ_GRADPS
#undef  TS_A4HADVECTION
#undef  TS_A4VADVECTION
#undef  TS_DIF2
#undef  MIX_GEO_TS
#define SOLVE3D
#define MASKING
#undef  ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#undef  ANA_VMIX

#ifdef CORRELATION
# define VCONVOLUTION
# define IMPLICIT_VCONV
#endif

#ifdef TLM_DRIVER
# define FORWARD_READ
# undef  FORWARD_MIXING
#endif

#ifdef OPT_PERTURBATION
# define FORWARD_READ
# undef  FORWARD_MIXING
#endif

#ifdef SANITY_CHECK
# define ANA_PERTURB
# define FORWARD_READ
# undef  FORWARD_MIXING
#endif

#ifdef AD_SENSITIVITY
# define FORWARD_READ
# undef  FORWARD_MIXING
#endif

#ifdef OPT_OBSERVATIONS
# define VCONVOLUTION
# define IMPLICIT_VCONV
# define FORWARD_READ
# undef  FORWARD_MIXING
#endif

#ifdef MULTI_DOMAIN
# undef MASKING
#endif
