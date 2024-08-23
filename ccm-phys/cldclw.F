#include <misc.h>
#include <params.h>
      subroutine cldclw(zi      ,clwp    ,tpw     ,hl      )
C-----------------------------------------------------------------------
C
C Evaluate cloud liquid water path clwp (g/m**2)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, J. Kiehl, August 1992
C Reviewed:          J. Kiehl April 1996
C
C-----------------------------------------------------------------------
c
c $Id: cldclw.F,v 1.1 1998/04/01 07:20:57 ccm Exp $
c
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real zi(plond,plevp)      ! height at layer interfaces(m)
      real tpw(plond)           ! total precipitable water (mm)
C
C Output arguments
C
      real clwp(plond,plev)     ! cloud liquid water path (g/m**2)
      real hl(plond)            ! liquid water scale height 
      real rhl(plond)           ! 1/hl
C
C---------------------------Local workspace-----------------------------
C
      integer i,k               ! longitude, level indices
      real clwc0                ! reference liquid water concentration (g/m**3)
      real emziohl(plond,plevp) ! exp(-zi/hl)
C
C-----------------------------------------------------------------------
C
C Set reference liquid water concentration
C
      clwc0 = 0.21
C
C Diagnose liquid water scale height from precipitable water
C
      do i=1,plon
        hl(i)  = 700.0*log(max(tpw(i)+1.0,1.0)) 
        rhl(i) = 1.0/hl(i)
      end do
C
C Evaluate cloud liquid water path (vertical integral of exponential fn)
C
      do k=1,plevp
        do i=1,plon
          emziohl(i,k) = exp(-zi(i,k)*rhl(i))
        end do
      end do
      do k=1,plev
        do i=1,plon
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k+1) - emziohl(i,k))
        end do
      end do
C
      return
      end
 
