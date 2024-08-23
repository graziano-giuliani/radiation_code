#include <misc.h>
#include <params.h>
      subroutine initesttr( q3 )
C-----------------------------------------------------------------------
C
C Initialize test tracers.  The test tracers are:
C
C    1) Radon, init to zero, surface fluxes from WCRP95, 5.5 
C       day e-folding decay.
C    2) conserved unit tracer
C    3) ozone-like tracer, init to 1.e-9 above ~100mb, zero 
C       elsewhere, re-zero the bottom level at each timestep.
C Note that:
C    o ixtrct   = index of radon advected tracer
C    o ixtrct+1 = index of conserved unit tracer
C    o ixtrct+2 = index of ozone-like tracer
C
C-------------------------Code History----------------------------------
C
C Original version:  B. Eaton, 1995
C Standardized:      T. Acker, Feb 1996
C Reviewed:          
C 
C-----------------------------------------------------------------------
#include <implicit.h>
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <comctl.h>
c-----------------------------------------------------------------------
#include <comindx.h>
C-----------------------------Arguments---------------------------------
C
C Output arguments:
C
      real q3(plond,plev,pcnst+pnats)    ! kg tracer/kg dry air
C
C--------------------------Local Variables------------------------------
C
      integer i, k                       !  loop counters
C
C-----------------------------------------------------------------------
C
c
c Initialize radon tracer to zero.
c
      if ( trace_test1 .or. trace_test2 .or. trace_test3 ) then
         do k = 1, plev
            do i = 1, plon
               q3(i,k,ixtrct) = 0.0
            end do
         end do
      end if
c
c Initialize conserved unit tracer.
c
      if ( trace_test2 .or. trace_test3 ) then
         do k = 1, plev
            do i = 1, plon
               q3(i,k,ixtrct+1) = 1.0
            end do
         end do
      end if
c
c Initialize strat tracer to 1.e-9 above 100mb
c
      if ( trace_test3 ) then
         do k = 1, plev
            do i = 1, plon
               q3(i,k,ixtrct+2) = 0.0
            end do
         end do
         do k = 1, 5
            do i = 1, plon
               q3(i,k,ixtrct+2) = 1.e-9
            end do
         end do
      end if

      return
      end
C
C#######################################################################
C
      subroutine rndecay( rn, deltat, rnsnk )
C-----------------------------------------------------------------------
C
C Radon decay.
C
C-------------------------Code History----------------------------------
C
C Original version:  B. Eaton, 1995
C Standardized:      T. Acker, Feb 1996
C Reviewed:          
C 
C-----------------------------------------------------------------------
#include <implicit.h>
c-------------------------Parameters------------------------------------ 
#include <pmgrid.h>
c-------------------------Arguments--------------------------------------
C
c input args
      real rn(plond,plev)       ! radon mixing ratio (kg/(kg moist air))
      real deltat               ! time step
C
c output args
      real rnsnk(plond,plev)    ! conversion rate 
c                               !              (kg rn /(s kg moist air))
C
C--------------------------Local Variables------------------------------
C
      integer i                 ! x index
      integer k                 ! z index
C
      real a                    ! lifetime
      parameter( a = 2.10e-6 )
C
C-----------------------------------------------------------------------
C
C   calculate tendencies using Euler Backward
C
      do k = 1,plev
         do i = 1,plon
            rnsnk(i,k) = -rn(i,k)*a / (1. + a*deltat)
         end do
      end do
C
C
      return
C
      end
C
c########################################################################
C
      subroutine rnsfwcrp( rlat, oro, flux )
C-----------------------------------------------------------------------
c
c Set surface fluxes for radon for WCRP95 RN-PB simulation.
c
c  The flux is specified non-zero over land between 60S - 70N, except
c  exclude Greenland.
c
c  Flux strength:
c  60S - 60N:  3.69e-21 kg/m^2/s
c  60N - 70N:  (3.69e-21)/2 kg/m^2/s
c
c  This land source is has been adjusted so that the total radon flux is
c  15 kg/yr for a T42 grid.
c
C-------------------------Code History----------------------------------
C
C Original version:  B. Eaton, 1995
C Standardized:      T. Acker, Feb 1996
C Reviewed:          
C 
C-----------------------------------------------------------------------
#include <implicit.h>
c--------------------------Parameters------------------------------------
#include <pmgrid.h>
c--------------------------Arguments-------------------------------------
c
c Input arguments:
c
      real rlat           ! model latitude in radians
      real oro(plond)     ! ocean=0, land=1, sea ice=2
c
c Output arguments:
c
      real flux(plond)    ! specified radon flux in kg/m^2/s
c
C--------------------------Local Variables------------------------------
c
      integer i      ! loop counter
      integer lon300 ! longitude index at left side of interval 
c                    ! containing 300E.
      real landflx   ! land flux
      real landflxn  ! (land flux)/2
      real rad2deg   ! convert radians to degrees
      real latdeg    ! latitude in degrees
C
C--------------------------Parameters-----------------------------------
C
      parameter( rad2deg = 360. / 6.283185308)
      parameter( lon300 = (300. / 360. * plon) + 1 )
c
c------------------------------------------------------------------------
c
c      landflx = 3.69e-21
      landflx = 3.7796e-21   ! rescaled so total flux is 15 kg/yr (T42)
      landflxn = landflx/2.
c
      latdeg = rlat * rad2deg
c
      do i = 1, plon
         flux(i) = 0.
      end do
c
      if ( latdeg .ge. -60.  .and.  latdeg .le. 60. ) then    ! 60S - 60N
c
         do i = 1, plon
            if( nint(oro(i)) .eq. 1 ) flux(i) = landflx
         end do
c
      else if ( latdeg .gt. 60. .and. latdeg .le. 70 ) then   ! 60N - 70N
c
         do i = 1, lon300                  ! 0 - 300E excludes Greenland
            if( nint(oro(i)) .eq. 1 ) flux(i) = landflxn
         end do
c
      end if
c
      return
c
      end






 
