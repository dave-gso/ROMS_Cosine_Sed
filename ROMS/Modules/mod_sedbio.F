#include "cppdefs.h"
      MODULE mod_sedbio

/*
** svn $Id: mod_sedbed.F 751 2015-01-07 22:56:36Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2015 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module declares  Sediment Biology Model Kernel Structure
** containing variables describing benthic biology                    **
**                                                                    **
** Note that the  *.h file is located in ROMS/Nonlinear/Biology and  **
** it is included within <...> to allow the user customize it in the  **
** project directory, while keeping the distributed file intact       **
** (check build scripts for details).                                 **
**                                                                    **
************************************************************************
*/

#if defined SEDBIO
# include <bio_UMAINE15_sedbio_mod.h>
#endif

      END MODULE mod_sedbio
