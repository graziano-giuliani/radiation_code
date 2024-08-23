#include <misc.h>
#include <params.h>
      subroutine tsinti(tmeltx  ,latvapx ,rairx   ,stebolx ,dtimex  )
C-----------------------------------------------------------------------
C
C Initialize surface temperature calculation constants
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          B. Briegleb, J. Hack, Aug 1992, March 1996
C
C-----------------------------------------------------------------------
c
c $Id: tsinti.F,v 1.1 1998/04/01 07:22:49 ccm Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <comtsc.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real tmeltx         ! Melting temperature of snow and ice
      real latvapx        ! Latent heat of vaporization
      real rairx          ! Gas constant for dry air
      real stebolx        ! Stefan-Boltzmann constant
      real dtimex         ! Model time step (delta-t)
C
C-----------------------------------------------------------------------
C
      latice = 3.336e5    ! Latent heat of fusion at 0'C = 3.336e5 J/Kg
      tmelt  = tmeltx
      latvap = latvapx
      rair   = rairx
      stebol = stebolx
      dtime  = dtimex
      snwedp = 10.0       ! 10:1 Snow:water equivalent depth factor
C
      return
      end
 
