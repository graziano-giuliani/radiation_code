#include <misc.h>
#include <params.h>
      subroutine cldint(pmid    ,rpdeli  ,t       ,q       ,omga    ,
     $                  cldtop  ,cldbot  ,piln    ,pmln    ,tvm     ,
     $                  zm      ,cld     ,clc     ,clwp    ,emis    ,
     $                  effcld  ,cltot   ,cllow   ,clmed   ,clhgh   ,
     $                  oro     ,rel     ,rei     ,fice    ,pdel    , 
     $                  tpw     ,hl      ,cmfmc   ,ps      ,snowh   )
C-----------------------------------------------------------------------
C
C Interface routine for cloud fraction evaluation
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Ice Included       J.T. Kiehl, January 1993
C Reviewed:          J. Hack, J. Kiehl, August 1992
C Reviewed:          J. Kiehl, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: cldint.F,v 1.1.26.2 1999/01/28 22:37:41 erik Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pmid(plond,plev)     ! midpoint pressures
      real rpdeli(plond,plev)   ! 1./(pmid(k+1)-pmid(k))
      real pdel(plond,plev)     ! pressure depth of layer
      real t(plond,plev)        ! temperature
      real q(plond,plev)        ! specific humidity
      real omga(plond,plev)     ! vertical pressure velocity
      real cldtop(plond)        ! top level of convective activity
      real cldbot(plond)        ! bottom level of convective activity
      real piln(plond,plevp)    ! log of interface pressures
      real pmln(plond,plev)     ! log of midpoint pressures
      real tvm(plond,plev)      ! virtual temperature 
      real zm(plond,plev)       ! height of midpoints (above surface)
      real cmfmc(plond,plev)    ! convective mass flux--m sub c
      real ps(plond)            ! surface pressure
      real snowh(plond)         ! snow depth (liquid water equivalent)
C
C Output arguments
C
      real cld(plond,plevp)     ! cloud fraction
      real clc(plond)           ! column convective cloud fraction
      real clwp(plond,plev)     ! cloud liquid water path
      real emis(plond,plev)     ! cloud emissivity
      real effcld(plond,plevp)  ! effective cloud=cld*emis
      real cltot(plond)         ! diagnosed tot random overlap cld cover
      real cllow(plond)         ! diagnosed low random overlap cld cover
      real clmed(plond)         ! diagnosed med random overlap cld cover
      real clhgh(plond)         ! diagnosed hgh random overlap cld cover
      real oro(plond)           ! Land/ocean/seaice flag
      real rel(plond,plev)      ! effective drop radius (microns)
      real rei(plond,plev)      ! ice effective drop size (microns)
      real fice(plond,plev)     ! fractional ice content within cloud
      real tpw(plond)           ! total precipitable water
      real hl(plond)            ! liquid water scale height 
C
C---------------------------Local workspace-----------------------------
C
      integer i,k               ! longitude,level indices

      real zi(plond,plevp)      ! interface geometric height
      real rgrav                ! inverse gravitational acceleration
C
C-----------------------------------------------------------------------
C
C Compute cloud amount
C
      call cldfrc(pmid    ,rpdeli  ,t       ,q       ,omga    ,
     $            cldtop  ,cldbot  ,cld     ,clc     ,pdel    ,
     $            cmfmc   ,oro     ,snowh   )
C
C Interface heights required for cldclw
C
      call zint(piln    ,pmln    ,rair    ,gravit  ,tvm     ,
     $          zm      ,zi      )
C
C Cloud liquid water path
C
C Begin by diagnosing total preciptable water in column (in mm)
C
      do i=1,plon
        tpw(i) = 0.0
      end do
      rgrav = 1.0/gravit
      do k=1,plev
        do i=1,plon
          tpw(i) = tpw(i) + pdel(i,k)*q(i,k)*rgrav
        end do
      end do
c
      call cldclw(zi      ,clwp    ,tpw     ,hl      )
C
C Cloud particle size and fraction of ice
C
      call cldefr(oro     ,t       ,rel     ,rei     ,fice    ,
     $            ps      ,pmid    )
C
C Cloud emissivity
C
      call cldems(clwp    ,fice    ,rei     ,emis    )
C
C Cloud cover diagnostics (required to get random overlap correct)
C
      call cldsav(cld     ,pmid    ,cltot   ,cllow   ,clmed   ,
     $            clhgh   )
C
C Effective cloud cover
C
      do k=1,plev
        do i=1,plon
          effcld(i,k) = cld(i,k)*emis(i,k)
        end do
      end do
C
C Cloud cover at surface interface always zero (for safety's sake)
C
      do i=1,plon
        effcld(i,plevp) = 0.
        cld(i,plevp)    = 0.
      end do
C
      return
      end
 
