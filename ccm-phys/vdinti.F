#include <misc.h>
#include <params.h>
      subroutine vdinti(cpwvx   ,cpairx  ,gravx   ,rairx   )
C-----------------------------------------------------------------------
C
C Initialization of time independent fields for vertical diffusion.
C Call initialization routine for boundary layer scheme.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: vdinti.F,v 1.1 1998/04/01 07:22:52 ccm Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C-----------------------------------------------------------------------
#include <comhyb.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real cpwvx        ! spec. heat of water vapor at const. pressure
      real cpairx       ! specific heat of dry air
      real gravx        ! acceleration due to gravity
      real rairx        ! gas constant for dry air
C
C---------------------------Local workspace-----------------------------
C
      integer k         ! vertical loop index
C
C-----------------------------------------------------------------------
C
C Hard-wired numbers.
C zkmin = minimum k = kneutral*f(ri)
C
      zkmin = 0.01
C
C Set physical constants for vertical diffusion and pbl:
C
      cpair  = cpairx
      gravit = gravx
      rair   = rairx
C
C Derived constants
C ntopfl = top level to which v-diff is applied
C npbl = max number of levels (from bottom) in pbl
C
c
c Limit pbl height to regions below 400 mb
c
      npbl = 0
      do k=plev,1,-1
        if (hypm(k).ge.4.e4) then
          npbl = npbl + 1
        end if
      end do
      npbl = max(npbl,1)
      if (masterproc) then
        write(6,*)'VDINTI: PBL height will be limited to bottom ',npbl,
     $' model levels. Top is ',hypm(plevp-npbl),' pascals'
      end if

      ntopfl = 1
      if (plev.eq.1) ntopfl = 0
      cpvir = cpwvx/cpairx - 1.
C
C Set the square of the mixing lengths.
C
      ml2(1) = 0.
      do k=2,plev
        ml2(k) = 30.0**2
      end do
      ml2(plevp) = 0.
C
C Initialize pbl variables
C
      call pbinti(gravx)
C
      return
      end
 
