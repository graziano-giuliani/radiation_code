#include <misc.h>
#include <params.h>
      subroutine cldsav(cld     ,pmid    ,cldtot  ,cldlow  ,cldmed  ,
     $                  cldhgh  )
C-----------------------------------------------------------------------
C

C Compute total & 3 levels of cloud fraction assuming random overlap.
C These diagnostics must be done on the fly to correctly deal with
C random overlap. Pressure ranges for the 3 cloud levels are specified.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      J. Rosinski,       Jun 1992
C                    L. Buja,           Feb 1996
C Reviewed:          J. Hack, J. Kiehl, Aug 1992, Apr 1996
C
C-----------------------------------------------------------------------
c
c $Id: cldsav.F,v 1.1 1998/04/01 07:21:07 ccm Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
      real plowmax             ! Max prs for low cloud cover range
      real plowmin             ! Min prs for low cloud cover range
      real pmedmax             ! Max prs for mid cloud cover range
      real pmedmin             ! Min prs for mid cloud cover range
      real phghmax             ! Max prs for hgh cloud cover range
      real phghmin             ! Min prs for hgh cloud cover range
C
      parameter (plowmax = 120000.,plowmin = 70000.,
     $           pmedmax =  70000.,pmedmin = 40000.,
     $           phghmax =  40000.,phghmin =  5000.)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real cld(plond,plevp)    ! Cloud fraction
      real pmid(plond,plev)    ! Level pressures
C
C Output arguments
C
      real cldtot(plond)       ! Total random overlap cloud cover
      real cldlow(plond)       ! Low random overlap cloud cover
      real cldmed(plond)       ! Middle random overlap cloud cover
      real cldhgh(plond)       ! High random overlap cloud cover
C
C---------------------------Local workspace-----------------------------
C
      integer i,k              ! Longitude,level indices
      real clrsky(plond)       ! Random overlap clear sky fraction
C
C-----------------------------------------------------------------------
C
C Total cloud
C
      do i=1,plon
        clrsky(i) = 1.0
      end do
      do k=plev,1,-1
        do i=1,plon
          clrsky(i) = clrsky(i)*(1. - cld(i,k))
        end do
      end do
      do i=1,plon
        cldtot(i) = 1. - clrsky(i)
      end do
C
C High cloud (400mb -> 50mb)
C
      do i=1,plon
        clrsky(i) = 1.0
      end do
      do k=plev,1,-1
        do i=1,plon
          if (phghmin.le.pmid(i,k) .and. pmid(i,k).lt.phghmax) then
            clrsky(i) = clrsky(i)*(1. - cld(i,k))
          end if
        end do
      end do
      do i=1,plon
        cldhgh(i) = 1. - clrsky(i)
      end do
C
C Medium cloud (700mb -> 400mb)
C
      do i=1,plon
        clrsky(i) = 1.0
      end do
      do k=plev,1,-1
        do i=1,plon
          if (pmedmin.le.pmid(i,k) .and. pmid(i,k).lt.pmedmax) then
            clrsky(i) = clrsky(i)*(1. - cld(i,k))
          end if
        end do
      end do
      do i=1,plon
        cldmed(i) = 1. - clrsky(i)
      end do
C
C Low cloud (sfc -> 700mb)
C
      do i=1,plon
        clrsky(i) = 1.0
      end do
      do k=plev,1,-1
        do i=1,plon
          if (plowmin.le.pmid(i,k) .and. pmid(i,k).lt.plowmax) then
            clrsky(i) = clrsky(i)*(1. - cld(i,k))
          end if
        end do
      end do
      do i=1,plon
        cldlow(i) = 1. - clrsky(i)
      end do
C
      return
      end
 
