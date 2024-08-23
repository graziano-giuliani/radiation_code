#include <misc.h>
#include <params.h>
      subroutine cldfrc(pmid    ,rpdeli  ,temp    ,q       ,omga    ,
     $                  cldtop  ,cldbot  ,cloud   ,clc     ,pdel    ,
     $                  cmfmc   ,oro     ,snowh   )
C-----------------------------------------------------------------------
C
C Compute cloud fraction using scheme of J.M.Slingo, 
C as modified by J.J.Hack and J.T.Kiehl
C
C This scheme is based on the operational scheme used in the ECMWF model
C A full description of its development can be found in Slingo (1987),
C which appears in the QJRMS July issue.  A number of modifications have
C been introduced to the original scheme in the following implementation 
C
C---------------------------Code history--------------------------------
C
C Original version: Based on code from J. Slingo 
C Modified:         J. J. Hack, May 1990
C Modified:         J. T. Kiehl, June, 1990
C Modified:         J. J. Hack, January 1991
C Rewritten:        B. P. Briegleb, November, 1991
C Standardized:     J. Rosinski, June 1992
C Last updated:     J. J. Hack, July 1992 
C Reviewed:         J. Kiehl, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: cldfrc.F,v 1.1 1998/04/01 07:21:02 ccm Exp $
c $Author: ccm $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
      real pnot                  ! reference pressure
      parameter (pnot = 1.e5)
C------------------------------Commons----------------------------------
#include <comcon.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pmid(plond,plev)      ! midpoint pressures
      real rpdeli(plond,plev)    ! 1./(pmid(k+1)-pmid(k))
      real temp(plond,plev)      ! temperature
      real q(plond,plev)         ! specific humidity
      real omga(plond,plev)      ! vertical pressure velocity
      real cldtop(plond)         ! top level of convection
      real cldbot(plond)         ! bottom level of convection
      real cmfmc(plond,plev)     ! convective mass flux--m sub c
      real snowh(plond)          ! snow depth (liquid water equivalent)
      real pdel(plond,plev)      ! pressure depth of layer
      real oro(plond)            ! Land/ocean/seaice flag
C
C Output arguments
C
      real cloud(plond,plevp)    ! cloud fraction
      real clc(plond)            ! column convective cloud amount
C
C---------------------------Local workspace-----------------------------
C
      real concld(plond,plev)    ! convective cloud cover
      real cld                   ! intermediate scratch variable (low cld)
      real cld8(plond)           ! low cloud fraction estimate
      real cld9(plond)           ! mid and high cloud fraction estimate
      real cck(plond)            ! convective cloud per level (assuming
C                                  random overlap in convective layer)
      real dthtdp(plond,plev)    ! lapse rate (d theta/dp) below 750 mb
      real dtdpmn(plond)         ! most stable lapse rate below 750 mb
      real dthdp                 ! lapse rate (intermediate variable)
      real es(plond,plev)        ! saturation vapor pressure
      real qs(plond,plev)        ! saturation specific humidity
      real premib                ! bottom pressure bound of middle cloud
      real pretop                ! pressure bounding high cloud
      real rh(plond,plev)        ! relative humidity
      real rhb                   ! intermediate scratch variable
      real rhdif                 ! intermediate scratch variable
      real strat                 ! intermediate scratch variable
      real theta(plond,plev)     ! potential temperature
      real zrth                  ! reciprocal of no. of convective layers
      real bvf                   ! brunt-vaisalla frequency
      real rbvflim               ! bound on inverse of bvf
      real rho                   ! local density (used to calculate bvf)
      real rhlim                 ! local rel. humidity threshold estimate 
      real rhden                 ! intermediate scratch variable
      real rhdif2                ! intermediate scratch variable
      real pdepth                ! intermediate scratch variable
      real stratfac              ! intermediate scratch variable
      real rhminl                ! minimum rh for low stable clouds
      real rhminh                ! minimum rh for high stable clouds
      real mcbar(plond)          ! mean convective scale motion in column
      real dpsum(plond)          ! vertical sum of delta-p (k-1 levels)
      real ccldt(plond)          ! estimate of total convective cloud
      real coef1                 ! coefficient to convert mass flux to mb/d
      real clrsky(plond)         ! temporary used in random overlap calc

      logical lol(plond)         ! region of low level cloud
      logical cldbnd(plond)      ! region below high cloud boundary

      integer i,k                ! longitude, level indices
c
c Statement functions
c
      logical land
      land(i) = nint(oro(i)).eq.1
C
C Set bound for inverse of brunt-vaisalla frequency and minimum relative
C humidity thresholds for stable clouds.  These are the principal 
C "disposable" parameters for the cloud fraction scheme
C
      rbvflim = 1./0.00035
      rhminl = .90
      rhminh = .90
C
C Evaluate potential temperature and relative humidity
C
      call aqsat(temp    ,pmid    ,es      ,qs      ,plond   ,
     $           plon    ,plev    ,1       ,plev    )
      do k=1,plev
        do i=1,plon
          theta(i,k)  = temp(i,k)*(pnot/pmid(i,k))**cappa
          rh(i,k)     = q(i,k)/qs(i,k)
          cloud(i,k)  = 0.
          concld(i,k) = 0.
        end do
      end do
C
C Initialize other temporary variables
C
      do i=1,plon
        cck(i) = 0.0    
        clc(i) = 0.0    
      end do
      coef1 = gravit*864.0    ! conversion to millibars/day
      do i=1,plon
        mcbar(i) = 0.0
        dpsum(i) = 0.0
      end do
c
c Calculate mean convective motion throughout column (in units of mb/day)
c
      do k=1,plev-1
        do i=1,plon
          mcbar(i) = mcbar(i) + max(cmfmc(i,k+1)*coef1,0.)*pdel(i,k)
          dpsum(i) = dpsum(i) + pdel(i,k)
        end do
      end do
c
c Estimate of total convective cloud cover based on mean convective motion
c
      do i=1,plon
        mcbar(i) = max(mcbar(i)/dpsum(i),1.0e-15)
        ccldt(i) = min(0.035*log(1.0+mcbar(i)),0.80)
        if ((cldbot(i) - cldtop(i)).ge.1.0) then
C
C Inverse of depth of convection (depth is expressed in model levels)
C
          zrth = 1.0/(cldbot(i) - cldtop(i))
C
C Compute amount of convective cloud at each level so that
C after random overlap, the total convective cloud cover is ccldt
C
          cck(i) = 1.0 - (1.0 - ccldt(i))**zrth
        end if
      end do
c
c Vertically distribute cloud in convective layer
c
      do k=1,plev-1
        do i=1,plon
          if (k.le.cldbot(i) .and. k.ge.cldtop(i)) then
            concld(i,k) = cck(i)
            rh(i,k) = (rh(i,k) - concld(i,k))/(1.0 - concld(i,k))
          end if
        end do
      end do
C
C Evaluate effective column-integrated convective cloud cover using
C random overlap assumption (for diagnostic purposes only)
C
      do i=1,plon
        clrsky(i) = 1.0
      end do
      do k=plev,1,-1
        do i=1,plon
          clrsky(i) = clrsky(i)*(1. - concld(i,k))
        end do
      end do
      do i=1,plon
        clc(i) = 1. - clrsky(i)
      end do
C
C          ****** Compute layer cloudiness ******
C
C There is effecively no top for high cloud formation (can for all the way
C up to 1mb)
C The bottom of middle level cloud (or the top of low level cloud) is
C arbitrarily define to be 750 mb (premib)
C
      premib = 750.e2
      pretop = 1.0e2                 ! top of cloud layer is at 1 mb
C
C Find most stable level below 750 mb for evaluating stratus regimes
C
      do i=1,plon
        dtdpmn(i) = 0.0     
      end do
      do k=2,plev-2
        do i=1,plon
          if (pmid(i,k).ge.premib) then
            dthdp = 100.0*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
          else
            dthdp = 0.0
          end if
          if (dthdp.lt.dtdpmn(i)) then
            dtdpmn(i) = dthdp
          end if
          dthtdp(i,k) = dthdp
        end do
      end do
      do k=plev-1,plev
        do i=1,plon
          if (0.0.lt.dtdpmn(i)) then
            dtdpmn(i) = 0.0
          end if
          dthtdp(i,k) = 0.0
        end do
      end do
C
C bvf => brunt-vaisalla frequency (approx. 1-sided diff.)
C this stability measure is used to set a local relative humidity 
C threshold when evaluating the fractional area of layered cloud
C
      do 10 k=2,plev-1
        do i=1,plon
          if (dthtdp(i,k).gt.dtdpmn(i)) then
            dthtdp(i,k) = 0.
          end if
          cldbnd(i) = pmid(i,k).ge.pretop
          lol(i) = pmid(i,k).ge.premib
          rho = pmid(i,k)/(rair*temp(i,k))
          bvf = -rho*gravit*gravit*((theta(i,k)-theta(i,k-1))*
     $         rpdeli(i,k-1))/theta(i,k)
          if (cldbnd(i)) then
            rhlim = 0.999 - (1.0-rhminh)*
     $             (1.0-min(1.0,max(0.0,bvf*rbvflim)))
            rhden = 1.0 - rhlim
          else
            rhlim = 0.999
            rhden = 0.001
          end if
          rhdif = (rh(i,k) - rhlim)/rhden
          cld9(i) = min(0.999,(max(rhdif,0.0))**2)
C
C Ignore brunt-vaisalla stability estimate of local relative humidity
C threshold when evaluating low cloud where local vertical motion is 
C less than some prescribed value (see low cloud section below)
C Relative humidity threshold is fixed at rhminl for this case, except
C over snow-free land, where it is reduced by 10%.  This distinction is
C made to account for enhanced cloud drop nucleation ({\it i.e.,} at 
C lower relative humidities) that can occur over CCN rich land areas.
C
          if (lol(i)) then
            if (land(i) .and. (snowh(i).le.0.000001)) then
              rhlim = rhminl - 0.10
            else
              rhlim = rhminl
            endif
            rhdif2 = (rh(i,k) - rhlim)/(1.0-rhlim)
            cld8(i) = min(0.999,(max(rhdif2,0.0))**2)
          else
            cld8(i) = cld9(i)
          end if
        end do
C
C Final evaluation of layered cloud fraction
C
        do i=1,plon
C
C Low cloud: non-zero only if vertical velocity requirements are satisfied
C Current vertical velocity threshold is omega < +50 mb/day with a 50 mb/day 
C linear ramp (other quantities in the class of "disposable" parameters)
C
          if (lol(i)) then
            if (omga(i,k).lt.0.05787) then
              cld = cld8(i)*
     $              min(1.0,max(0.0,(0.05787-omga(i,k))/0.05787))
            else
              cld = 0.0
            end if
C
C Compute cloud associated with low level inversions.
C
            strat = max(0.,min(0.95,-6.67*dthtdp(i,k) - 0.667))
            rhb   = 1.0 - (0.9 - rh(i,k+1))/0.3
            if (rh(i,k+1).lt.0.6) then
              strat = 0.0
            end if
            if (rh(i,k+1).ge.0.6 .and. rh(i,k+1).le.0.9) then
              strat = strat*rhb
            end if
C
C Linear transition from stratus to trade cu as inversion rises.
C Transition starts at 900mb and completes at 750mb (premib)
C
            pdepth = max(pmid(i,k) - 750.e2,0.0)
            stratfac = min(pdepth,150.0e2)/150.e2
            if (dthtdp(i,k).le.-0.125 ) then
              cloud(i,k) = strat*stratfac
            else
              cloud(i,k) = cld
            end if
          else                  ! Middle and high level cloud 
            if( cldbnd(i) ) then
              cloud(i,k) = cld9(i)
            else
              cloud(i,k) = 0.0
            end if
          end if
        end do
   10 continue                  ! k=2,plev-1
C
C Merge convective and layered cloud fraction for total cloud
C
      do k=1,plev
        do i=1,plon
          cloud(i,k) = max(0.0,min(0.999,
     $                 (1.0 - concld(i,k))*cloud(i,k) + concld(i,k)))
        end do
      end do
C
      return
      end
 
