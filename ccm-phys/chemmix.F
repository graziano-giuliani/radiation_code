#include <misc.h>
#include <params.h>
      subroutine chemmix(lat, ps, q, n2o, ch4, cfc11, cfc12)
C-----------------------------------------------------------------------
C
c Specify initial mass mixing ratios of CH4, N2O, CFC11 and CFC12.
C Distributions assume constant mixing ratio in the troposphere
C and a decrease of mixing ratio in the stratosphere. Tropopause
C defined by ptrop. The scale height of the particular trace gas
C depends on latitude. This assumption produces a more realistic
C stratospheric distribution of the various trace gases.
C
C-------------------------Code History----------------------------------
C
C Original version:  J.T.Kiehl Nov 21 1994
C Standardized:      T. Acker, Feb 1996
C Reviewed:          J. Kiehl, April 1996        
C 
C-----------------------------------------------------------------------
C
C $Id: chemmix.F,v 1.1 1998/04/01 07:20:56 ccm Exp $
C
#include <implicit.h>
C-----------------------------Parameters--------------------------------
#include <prgrid.h>
C-----------------------------------------------------------------------
#include <pspect.h>
C------------------------------Commons----------------------------------
#include <comhyb.h>
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comvmr.h>
C-----------------------------Arguments---------------------------------
C
C Input
C
      integer lat                  ! current latitude index
      real ps(plond)               ! surface pressure
C
C Input/Output
C
      real q(plond,plev)           ! water mass mixing ratio
C
C Output
C
      real n2o(plond,plev)         ! nitrous oxide mass mixing ratio
      real ch4(plond,plev)         ! methane mass mixing ratio
      real cfc11(plond,plev)       ! cfc11 mass mixing ratio
      real cfc12(plond,plev)       ! cfc12 mass mixing ratio
C
C--------------------------Local Variables------------------------------
C
      integer i                ! longitude loop index
      integer k                ! level index
C
      real coslat              ! cosine of latitude
      real dlat                ! latitude in degrees
      real pmid(plond,plev)    ! model pressures
      real ptrop               ! pressure level of tropopause
      real pratio              ! pressure divided by ptrop
C
      real xn2o                ! pressure scale height for n2o
      real xch4                ! pressure scale height for ch4
      real xcfc11              ! pressure scale height for cfc11
      real xcfc12              ! pressure scale height for cfc12
C
      real ch40                ! tropospheric mass mixing ratio for ch4
      real n2o0                ! tropospheric mass mixing ratio for n2o
      real cfc110              ! tropospheric mass mixing ratio for cfc11
      real cfc120              ! tropospheric mass mixing ratio for cfc12
C
C-----------------------------------------------------------------------
c$$$For 1990, I get the following:
c$$$
c$$$CFC-11 mixing ratio = 270.1 pptv
c$$$CFC-12 mixing ratio = 465.3 pptv
c$$$
c$$$CFC-11 radiative forcing = 0.0594 W/m2
c$$$CFC-12 radiative forcing = 0.1300 W/m2
c$$$Other CFC radiative forcing = 0.0595 W/m2
c$$$
c$$$Hence, the magic number is "TWO".   
c$$$
c$$$These values compare reasonably well to IPCC.   They have slightly higher
c$$$CFC-11 and 12 mixing ratios in the l990 report (280 and 484), but I think
c$$$mine are better.   They had 0.062 and 0.14 for the corresponding forcings. 
c$$$They also had higher CCl4 and CFC-113 mixing ratios (in one case by quite a
c$$$bit) and had radiative forcing by other CFCs at 0.085.  I think the above
c$$$number can be defended as better.
c$$$
c$$$Just in case you want it, my estimates for l990 CH4 and N2O would be:
c$$$
c$$$CH4 - 1.715 ppmv
c$$$N2O - 310.0 ppbv
c$$$
c$$$but you probably want to use Tom's numbers instead for those.
c$$$
c$$$Thanks,
c$$$Susan
C-----------------------------------------------------------------------
C
C 1990 ch4 vmr from Wigley   : 1.722e-6
C 1990 n2o vmr from Wigley   : 308.4e-9
C 1990 cfc11 vmr from Solomon: 270.1e-12 * 2.
C factor of 2 on cfc11 to account for all other halocarbons
C 1900 cfc12 vmr from Solomon: 465.3e-12
C
      ch40   = rmwch4 * ch4vmr
      n2o0   = rmwn2o * n2ovmr
      cfc110 = rmwf11 * f11vmr
      cfc120 = rmwf12 * f12vmr
C
C Set stratospheric scale height factor for gases
C
      dlat = abs(latdeg(lat))
      if(dlat.le.45.0) then
         xn2o = 0.3478 + 0.00116 * dlat
         xch4 = 0.2353
         xcfc11 = 0.7273 + 0.00606 * dlat
         xcfc12 = 0.4000 + 0.00222 * dlat
      else
         xn2o = 0.4000 + 0.013333  * (dlat - 45)
         xch4 = 0.2353 + 0.0225489 * (dlat - 45)
         xcfc11 = 1.00 + 0.013333  * (dlat - 45)
         xcfc12 = 0.50 + 0.024444  * (dlat - 45)
      end if
C
C Set pressure of tropopause and model layer pressures
C
      coslat = cos(clat(lat))
      ptrop = 250.0e2 - 150.0e2*coslat**2.0

      do k=1,plev
         do i=1,plon
            pmid(i,k) = hyam(k)*ps0 + hybm(k)*ps(i)
         end do
      end do
C
C Determine initial mixing ratios
C
      do k = 1,plev
         do i = 1,plon
            if(pmid(i,k).ge.ptrop) then
              ch4(i,k) = ch40
              n2o(i,k) = n2o0
              cfc11(i,k) = cfc110
              cfc12(i,k) = cfc120
            else
              pratio = pmid(i,k)/ptrop
              ch4(i,k) = ch40 * (pratio)**xch4
              n2o(i,k) = n2o0 * (pratio)**xn2o
              cfc11(i,k) = cfc110 * (pratio)**xcfc11
              cfc12(i,k) = cfc120 * (pratio)**xcfc12
            end if
         end do
      end do
C
C Adjust water using methane mass mxiing ratio
C
      do k=1,plev
         do i=1,plon
            q(i,k) = q(i,k) + 2.*(ch40-ch4(i,k))
         end do
      end do
C
      return
      end
 
