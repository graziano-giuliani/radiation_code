#include <misc.h>
#include <params.h>
      subroutine torque(um1    ,psm1   ,pdelm1   ,phis     ,clatitude,
     $                  taux   ,taugsx ,taugtx   ,tor      ,mom      )
C-----------------------------------------------------------------------
C
C Calculate total torque and relative angular momementum.
C
C---------------------------Code history--------------------------------
C
C Original version:  E. Kluzek, December 1998
C
C-----------------------------------------------------------------------
c
c $Id: torque.F,v 1.1.2.1 1999/01/05 22:39:04 erik Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C-----------------------------------------------------------------------
#include <pspect.h>
      real pie                    ! Mathematical constant
      parameter( pie = 3.14159265358979323846 )
C------------------------------Commons----------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <commap.h>
C------------------------------Input Arguments--------------------------
      real um1(plond,plev)        ! U horizontal wind component
      real psm1(plond)            ! Surface Pressure
      real pdelm1(plond,plev)     ! Pressure depth in layer
      real phis(plond)            ! Surface geopotential
      real clatitude              ! latitude (radians)
      real taux(plond)            ! stress at surface (zonal)
      real taugtx(plond)          ! Gravity wave c=0 stress at model top (zonal)
      real taugsx(plond)          ! Gravity wave c=0 stress at surface (zonal)
C------------------------------Output Arguments-------------------------
      real tor(plond)             ! Torque
      real mom(plond)             ! Angular momementum
C------------------------------Local Variables--------------------------
      integer i, k                ! Loop variables
      real coslat                 ! cosine of the latitude
      real dlambda                ! delta longitude (for pressure difference across mountains)
      real dphis(plond)           ! delta in geopotential height divide by dlambda (dphis/dlambda)
      real sum_u                  ! mass weighted sum of east-west wind
C
      coslat = cos(clatitude)
      do i = 1, plon
        if ( i .eq. 1 )then
          dlambda  = clon(i+1) - clon(plon) + 2.0*pie
          dphis(i) = (phis(i+1) - phis(plon))/dlambda
        else if ( i .eq. plon )then
          dlambda  = clon(1) - clon(i-1) + 2.0*pie
          dphis(i) = (phis(1) - phis(i-1))/dlambda
        else
          dlambda  = clon(i+1) - clon(i-1)
          dphis(i) = (phis(i+1) - phis(i-1))/dlambda
        end if
      end do
      do i = 1, plon
C Mountain torque
        tor(i) = rga*psm1(i)*dphis(i)
C Friction torque
        tor(i) = tor(i) + rearth*coslat*taux(i)
C Gravity wave drag torque
        tor(i) = tor(i) + rearth*coslat*(taugsx(i) - taugtx(i))
C Angular momentum
        sum_u = 0.0
        do k = 1, plev
          sum_u = sum_u + um1(i,k)*pdelm1(i,k)
        end do
        mom(i) = rearth*rga*coslat*sum_u
      end do
C
      return
      end
 
