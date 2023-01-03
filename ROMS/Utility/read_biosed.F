#include "cppdefs.h"
#ifdef BIOLOGY
!!                                                                     !
!!  This routine reads and reports sediment biological model input     !
!!  parameters.                                                        !
!!                                                                     !
!!======================================================================
!!
# if defined BIO_UMAINE15 & defined SEDBIO
#  include <sedbio_inp.h>
# else
      SUBROUTINE read_SedBioPar
      END SUBROUTINE read_SedBioPar
# endif
#else
      SUBROUTINE read_SedBioPar
      END SUBROUTINE read_SedBioPar
#endif
