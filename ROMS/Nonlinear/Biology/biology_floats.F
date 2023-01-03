#include "cppdefs.h"
      MODULE biology_floats_mod
#if defined NONLINEAR && defined FLOATS && defined FLOAT_BIOLOGY
!
!svn $Id: biology_floats.F 751 2015-01-07 22:56:36Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routine computes the behavior for biological floats.          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: biology_floats

      CONTAINS

# ifdef FLOAT_OYSTER
#  include <oyster_floats.h>
# endif

#endif

      END MODULE biology_floats_mod
