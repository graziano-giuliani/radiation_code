#include <misc.h>
#include <params.h>
      real function qsat(t,p)
C-----------------------------------------------------------------------
C
C Utility procedure to look up saturation vapor pressure from 
C precomputed table and calculate saturation specific humidity (g/g)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Hack
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, G. Taylor, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: qsat.F,v 1.1 1998/04/01 07:22:11 ccm Exp $
c
#include <implicit.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real t                  ! temperature
      real p                  ! pressure
C
C---------------------------Local workspace-----------------------------
C
      real es                 ! saturation vapor pressure
C-----------------------------------------------------------------------
#include <eslookup.h>
C-----------------------------------------------------------------------
      es   = estblf(t)
      qsat = epsqs*es/(p - (1. - epsqs)*es)
C
C The following check is to avoid the generation of negative values
C that can occur in the upper stratosphere and mesosphere
C
      qsat = min(1.0,qsat)
      if (qsat .lt. 0.0) qsat = 1.0
C
      return
      end
 
