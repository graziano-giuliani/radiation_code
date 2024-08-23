#include <misc.h>
#include <params.h>
      subroutine mgwinti (cpx, cpwv, gx, rx)
C-----------------------------------------------------------------------
C
C Time independent initialization for multiple gravity wave 
C parameterization.
C
C-------------------------Code History----------------------------------
C
C Original version:  B. Boville
C Standardized:      T. Acker, Feb 1996
C Reviewed:          B. Boville, April 1996
C 
C-----------------------------------------------------------------------
c
c $Id: mgwinti.F,v 1.1 1998/04/01 07:21:58 ccm Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <mgw.h>
C-----------------------------------------------------------------------
C
C Input variables
C
      real cpx                  ! specific heat of dry air (constant p)
      real cpwv                 ! specific heat of water vapor (const. p)
      real gx                   ! acceleration of gravity
      real rx                   ! gas constant for dry air
C
C Local variables
C
      integer k
      real oroe                 ! efficiency factor for orographic term
C
C-----------------------------------------------------------------------
C
c
c Copy model constants
c
      cp     = cpx
      cpvir  = cpwv / cp - 1
      g      = gx
      r      = rx
c
c Set MGWD constants
c
      kwv    = 6.28e-5          ! 100 km wave length
      dback  = 0.05             ! background diffusivity
      orofc2 = OROFC2           ! critical froude number squared
      oroe   = OROE             ! efficiency for orographic generation
      tauscal= TAUSCAL          ! scale factor for background stress
      taunor = TAUNOR           ! northern hemisphere amplitude
      tausou = TAUSOU           ! southern hemisphere amplitude
      taueq  = TAUEQ            ! equatorial amplitude
      taueqw = TAUEQW*3.14159/180. ! equatorial width, exp{-(y/L)^2}
c
c Set phase speeds 
c
      do k = -pgwv, pgwv
         c(k)   = 10. * k       ! 0, +/- 10, +/- 20, ... m/s
      End Do
c
c Set radiative damping times
c
      do k = 0, plev
         alpha(k) = 1.e-6       ! about 10 days.
      end do
c
c Min and max values to keep things reasonable
c
      mxasym = MXASYM       ! max factor of 10 from |tau(c)| to |tau(-c)|
      mxrange= MXRANGE      ! factor of 100 from max to min |tau(c)|
      n2min  = 1.e-8        ! min value of Brunt-Vaisalla freq squared
      orohmin= 10.              ! min sdv for orographic wave drag
      orovmin=  2.              ! min wind speed for orographic wave drag
      taumin = 1.e-10           ! min stress considered > 0
      tndmax = 500. / 86400.    ! 250 m/s/day max tendency
      ubmc2mn= 0.01             ! min value of (u-c)^2
c
c Determine efficiency factor from gravity wave parameters
c
      oroeko2= 0.5 * oroe * kwv
      efcncy = oroe * orofc2
c
c Determine other derived constants
c
      rog    = r/g
      efkw   = efcncy * kwv
c
      return
      end
 
