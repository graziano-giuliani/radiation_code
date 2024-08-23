#include <misc.h>
#include <params.h>
      subroutine tphysbc(lat     ,clat    ,clon    ,coslat  ,hbuf    ,
     $                   oro     ,pmidm1  ,pintm1  ,pdelm1  ,omga    ,
     $                   um1     ,vm1     ,tm1     ,qm1     ,psm1    ,
     $                   pmlnm1  ,pilnm1  ,zm      ,asdir   ,asdif   ,
     $                   aldir   ,aldif   ,snowh   ,cnt     ,cnb     ,
     $                   qrs     ,qrl     ,flwds   ,fsns    ,tvm     ,
     $                   thm     ,cmfmc   ,lwup    ,precl   ,precc   ,
     $                   precsl  ,precsc  ,srfrad  ,sols    ,soll    ,
     $                   solsd   ,solld   )
C-----------------------------------------------------------------------
C
C Tendency physics BEFORE coupling to land, sea, and ice models.
C Call physics subroutines and compute the following:
C     o cloud calculations (cloud fraction, emissivity, etc.)
C     o radiation calculations
C Pass surface fields for separate surface flux calculations
C Dump appropriate fields to history file.
C
C----------------------------Code History-------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C                    T. Acker, March 1996
C Reviewed:          J. Kiehl, J. Hack, August 1992
C Include Ice:       J. Kiehl, January 1993
C Reviewed:          J. Kiehl, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: tphysbc.F,v 1.1.26.6 1999/01/28 22:37:46 erik Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
#include <comhst.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat               ! Model latitude index
      real clat                 ! Current latitude(radians)
      real clon(plon)           ! Centered longitude(radians)
      real coslat               ! Cosine latitude
      real hbuf(*)              ! History buffer
      real oro(plond)           ! Sfc type flag (ocean,land,sea ice)
      real pmidm1(plond,plev)   ! Midpoint pressures
      real pintm1(plond,plevp)  ! Interface pressures
      real pdelm1(plond,plev)   ! Layer pressure depth
      real omga(plond,plev)     ! Vertical pressure velocity
      real um1(plond,plev)      ! U-wind (time n-1)
      real vm1(plond,plev)      ! V-wind (time n-1)
      real tm1(plond,plev)      ! Temperature
      real qm1(plond,plev,pcnst+pnats) ! Specific humidity
      real psm1(plond)          ! Surface pressure
      real pmlnm1(plond,plev)   ! ln(pmidm1)
      real pilnm1(plond,plevp)  ! ln(pintm1)
      real zm(plond,plev)       ! Height of midpoints (above surface)
      real asdir(plond)         ! Albedo: shortwave, direct
      real asdif(plond)         ! Albedo: shortwave, diffuse
      real aldir(plond)         ! Albedo: longwave, direct
      real aldif(plond)         ! Albedo: longwave, diffuse
      real cmfmc(plond,plev)    ! Convective mass flux--m sub c
      real lwup(plond)          ! Surface longwave up flux
      real precl(plond)         ! Large scale precip. rate
      real precc(plond)         ! Convective precip. rate
C
C Output arguments
C
      real snowh(plond)         ! Snow depth (liquid water equivalent)
      real cnt(plond)           ! Top level of convective activity
      real cnb(plond)           ! Lowest level of convective activity
      real qrs(plond,plev)      ! Shortwave heating rate
      real qrl(plond,plev)      ! Longwave  heating rate
      real flwds(plond)         ! Surface longwave down flux
      real fsns(plond)          ! Surface solar absorbed flux
      real tvm(plond,plev)      ! Virtual temperature
      real thm(plond,plev)      ! Potential temperature
      real precsl(plond)        ! Large scale snowfall rate
      real precsc(plond)        ! Convective snowfall rate
      real srfrad(plond)        ! Net surface radiative flux (watts/m**2)
      real sols(plond)          ! Direct beam solar rad. onto srf (sw)
      real soll(plond)          ! Direct beam solar rad. onto srf (lw)
      real solsd(plond)         ! Diffuse solar radiation onto srf (sw)
      real solld(plond)         ! Diffuse solar radiation onto srf (lw)
C
C---------------------------Local workspace-----------------------------
C
      integer  i,k              ! Longitude, level indices
      real clwp(plond,plev)     ! Presribed cloud liq. h2o path 
      real rel(plond,plev)      ! Liquid cloud particle effective radius
      real rei(plond,plev)      ! Ice effective drop size (microns)
      real fice(plond,plev)     ! Fractional ice content within cloud
      real cld(plond,plevp)     ! Fraction cloud cover
      real effcld(plond,plevp)  ! Effective cloud=cld*emis
      real emis(plond,plev)     ! Cloud longwave emissivity
      real clc(plond)           ! Total convective cloud (cloud scheme)
      real rpdeli(plond,plev)   ! 1./(pmid(k+1)-pmid(k))
      real rpdel(plond,plev)    ! 1./(pint(k+1)-pint(k))
      real cltot(plond)         ! Diagnostic total cloud cover
      real cllow(plond)         !       "     low  cloud cover
      real clmed(plond)         !       "     mid  cloud cover
      real clhgh(plond)         !       "     hgh  cloud cover
      real coszrs(plond)        ! Cosine solar zenith angle
      real precst(plond)        ! Total large scale and convective snowfall rate
C        
      real actlwp(plond,plev)   ! Actual cloud liquid h2o path 
      real totcwp(plond)        ! Actual column cloud condensed h2o path 
      real totlwp(plond)        ! Actual column cloud liquid h2o path 
      real totiwp(plond)        ! Actual column cloud ice h2o path 
      real clwmr(plond,plev)    ! cloud LIQUID water mixing ratio
      real ciwmr(plond,plev)    ! cloud ICE    water mixing ratio
      real sum1                 ! Sum of actlwp
      real sum2                 ! liquid Sum of actlwp
      real sum3                 ! ice Sum of actlwp
      real cloudup(plond,plev)  ! Cloud fraction during omega < 0.
C                                 
      real tpw(plond)           ! Total precipitable water
      real hl (plond)           ! Liquid water scale height
C
C-----------------------------------------------------------------------
C
C Zero out snowfall before accumulating terms
C
      call resetr(precsl,plond,0.)
      call resetr(precsc,plond,0.)
C
C Cosine solar zenith angle for current time step
C
      call zenith (calday  ,clat    , clon, coszrs  )
C
C Precompute frequently used variables
C
      do k=1,plev
         do i=1,plon
            rpdel(i,k) = 1./pdelm1(i,k)
         end do
      end do
C
C 1./delta-p needed by cloud routines
C
      if (dosw .or. dolw) then
         do k=1,plev-1
            do i=1,plon
               rpdeli(i,k) = 1./(pmidm1(i,k+1) - pmidm1(i,k))
            end do
         end do
C
C Compute cloud properties for input to radiation
C
         call cldint(pmidm1  ,rpdeli  ,tm1     ,qm1     ,omga    ,
     $               cnt     ,cnb     ,pilnm1  ,pmlnm1  ,tvm     ,
     $               zm      ,cld     ,clc     ,clwp    ,emis    ,
     $               effcld  ,cltot   ,cllow   ,clmed   ,clhgh   ,
     $               oro     ,rel     ,rei     ,fice    ,pdelm1  ,
     $               tpw     ,hl      ,cmfmc   ,psm1    ,snowh   )
C
C Dump cloud field information to history tape buffer (diagnostics)
C
         call outfld('CLOUD   ',cld    ,plond,lat,hbuf)
         call outfld('EFFCLD  ',effcld ,plond,lat,hbuf)
         call outfld('CNVCLD  ',clc    ,plond,lat,hbuf)
         call outfld('CLDTOT  ',cltot  ,plond,lat,hbuf)
         call outfld('CLDLOW  ',cllow  ,plond,lat,hbuf)
         call outfld('CLDMED  ',clmed  ,plond,lat,hbuf)
         call outfld('CLDHGH  ',clhgh  ,plond,lat,hbuf)
         call outfld('LWSH    ',hl     ,plond,lat,hbuf)
C
C Special diagnostic cloud water fields:
C
         call outfld('SETLWP  ',clwp   ,plond,lat,hbuf)
C
C Outfield special monthly cloud diagnostics
C      
         if (ninavg(1).eq.'Q') then
            do k=1,plev
               do i=1,plon
                  cloudup(i,k) = 0.
                  if (omga(i,k).lt.0.) then
                     cloudup(i,k) = cld(i,k)
                  end if
               end do
            end do
            call outfld('CLOUDUP ',cloudup ,plond,lat,hbuf)
         end if
C
C Compute actual liquid water quantities
C
         do k=1,plev
           do i=1,plon
             actlwp(i,k) = cld(i,k)*clwp(i,k)
             clwmr(i,k) = (1.0-fice(i,k))*actlwp(i,k)*gravit*rpdel(i,k)
             ciwmr(i,k) = fice(i,k)*actlwp(i,k)*gravit*rpdel(i,k)
           end do
         end do
         call outfld('CLDCWP  ',actlwp ,plond,lat,hbuf)
         call outfld('CLWMR   ',clwmr  ,plond,lat,hbuf)
         call outfld('CIWMR   ',ciwmr  ,plond,lat,hbuf)
         do i=1,plon
           sum1 = 0.
           sum2 = 0.
           sum3 = 0.
           do k=1,plev
             sum1 = sum1 + actlwp(i,k)
             sum2 = sum2 + actlwp(i,k)*(1.0-fice(i,k))
             sum3 = sum3 + actlwp(i,k)*fice(i,k)
           end do
           totcwp(i) = sum1
           totlwp(i) = sum2
           totiwp(i) = sum3
         end do
         call outfld('TOTCWP  ',totcwp ,plond,lat,hbuf)
         call outfld('TOTLWP  ',totlwp ,plond,lat,hbuf)
         call outfld('TOTIWP  ',totiwp ,plond,lat,hbuf)
C
C Complete radiation calculations
C
         call radctl(hbuf    ,clat    ,coslat  ,lat     ,lwup    ,
     $               pmidm1  ,pintm1  ,pmlnm1  ,pilnm1  ,tm1     ,
     $               qm1     ,cld     ,effcld  ,clwp    ,coszrs  ,
     $               asdir   ,asdif   ,aldir   ,aldif   ,fsns    ,
     $               qrs     ,qrl     ,flwds   ,rel     ,rei     ,
     $               fice    ,sols    ,soll    ,solsd   ,solld   ,
     $               oro     )
      end if
C
C Compute net surface radiative flux for use by surface temperature code.
C Note that units have already been converted to mks in RADCTL.  Since
C fsns and flwds are in the buffer, array values will be carried across
C timesteps when the radiation code is not invoked.
C
      do i=1,plon
         srfrad(i) = fsns(i) + flwds(i)
      end do
      call outfld('SRFRAD  ',srfrad,plond,lat,hbuf)
C
C determine whether precipitation, prec, is frozen (snow) or not
C (Use same methodology as in LSM)
C
      do i=1,plon
        if (tm1(i,plev  ).le. (tmelt+2.2) ) then
          precsc(i) = precc(i)
          precsl(i) = precl(i)
        else
          precsc(i) = 0.
          precsl(i) = 0.
        end if
        precst(i) = precsl(i) + precsc(i)
      end do
C
C Save atmospheric fields to force surface models
C
      call srfxfer(lat           ,psm1        ,um1(1,plev) ,
     $                  vm1(1,plev)   ,tm1(1,plev),
     $             qm1(1,plev,1) ,thm(1,plev) ,zm(1,plev)  ,
     $                  pmidm1        ,rhoh2o     ,
     $             rpdel(1,plev) ,fsns        ,precst      )
C
      return
C
      end
 
