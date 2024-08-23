#include <misc.h>
#include <params.h>
      real function esat(t)
C-----------------------------------------------------------------------
C
C Utility procedure to look up saturation vapor pressure from 
C precomputed table 
C
C----------------------------Code History-------------------------------
C
C Original version:  J. Hack
C Standardized:      J. Rosinski, Jun 1992
C                    L. Buja,     Feb 1996
C Reviewed:          J. Hack,     Aug 1992, Feb 1996
C
C-----------------------------------------------------------------------
c
c $Id: esat.F,v 1.1 1998/04/01 07:21:41 ccm Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real t                  ! temperature
C
C-----------------------------------------------------------------------
#include <eslookup.h>
C-----------------------------------------------------------------------
C
      esat = estblf(t)
C
      return
      end
 
