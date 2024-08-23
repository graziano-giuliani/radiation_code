#include <misc.h>
#include <params.h>
      subroutine tphysidl(lat     ,hbuf    ,ztodt   ,pmidm1  ,pintm1  ,
     $                    pdelm1  ,pilnm1  ,etamid  ,coslat  ,t2      ,
     $                    fu      ,fv      ,um1     ,vm1     ,tm1     ,
     $                    qm1     ,psomc   ,thm     ,zm      ,pblh    ,
     $                    qpert   ,qrs     ,tpert   ,clat    ,shflx   ,
     $                    taux    ,tauy    ,cflx    ,idlflag )
C-----------------------------------------------------------------------
C
C Tendency physics after coupling to land, sea, and ice models.
C Computes the following: 
C   o Radon surface flux and decay (optional)
C   o Vertical diffusion and planetary boundary layer
C   o Rayleigh friction
C   o Multiple gravity wave drag
C   o Dump appropriate fields to history file.
C
C----------------------------Code History-------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C                    T. Acker, March 1996
C Reviewed:          J. Kiehl, J. Hack, August 1992
C Reviewed:          J. Kiehl April 1996
C
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
      common/combd/t3z,h3z
C
      real t3z(plev,plat),
     $     h3z(plev,plat)
C-----------------------------------------------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <comhyb.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat                 ! Model latitude index
      real hbuf(*)                ! History buffer
      real ztodt                  ! Two times model timestep (2 delta-t)
      real pmidm1(plond,plev)     ! Midpoint pressures
      real pintm1(plond,plevp)    ! Interface pressures
      real pdelm1(plond,plev)     ! Layer pressure depth
      real pilnm1(plond,plevp)    ! ln(pintm1)
      real etamid(plev)           ! midpoint values of eta (a+b)
      real coslat                 ! cosine(latitude)
      real t2(plond,plev)         ! Temperature tendency
      real fu(plond,plev)         ! Accumulated u tendency (from tphys)
      real fv(plond,plev)         ! Accumulated v tendency (from tphys)
      real um1(plond,plev)        ! U horizontal wind component
      real vm1(plond,plev)        ! V horizontal wind component
      real tm1(plond,plev)        ! Temperature
      real qm1(plond,plev,pcnst+pnats)! Specific humidity and constituents
      real psomc(plond,plev)      ! (psm1/pmidm1)**cappa
      real thm(plond,plev)        ! Potential temperature 
      real zm(plond,plev)         ! Height of midpoints (above surface) 
      real pmlnm1(plond,plev)     ! log(pmidm1)
      real zi(plond,plevp)        ! interface geometric height
      real tvm(plond,plev)        ! virtual temperature 
      real qrs(plond,plev)        ! Shortwave heating rate
      real clat                   ! Current latitude(radians)
      integer idlflag           ! Flag to choose which idealized physics
C                               ! to run.
C
C Output arguments
C
      real pblh(plond)            ! Planetary boundary layer height
      real qpert(plond,pcnst)     ! Moisture/constit. perturbation (PBL)
      real tpert(plond)           ! Temperature perturbation (PBL)
      real shflx(plond)           ! Sensible heat flux (w/m^2)
      real taux(plond)            ! X surface stress (zonal)
      real tauy(plond)            ! Y surface stress (meridional)
      real cflx(plond,pcnst)      ! Surface constituent flux (kg/m^2/s)
C
C---------------------------Local workspace-----------------------------
C
      real up1(plond,plev)        ! U-wind after vdiff & rfric
      real vp1(plond,plev)        ! V-wind after vdiff & rfric
      real dtv(plond,plev)        ! T tendency from vertical diffusion
      real duv(plond,plev)        ! U tendency from vertical diffusion
      real dvv(plond,plev)        ! V tendency from vertical diffusion
      real rpdeli(plond,plev)     ! 1./(pmid(k+1)-pmid(k))
      real rpdel(plond,plev)      ! 1./(pint(k+1)-pint(k))
      integer i,k                 ! Longitude, level indices
C
      real sum,coeff(plond)
      real tmp,kf,ka,kaa,ks,kv,kt,trefa,trefc,cossq,cossqsq,sinsq,
     $     onemsig
      real efoldf,efolda,efoldaa,efolds
      real constw,lapsew,p0strat,phi0,dphi0,a0,aeq,apole,pi,acoslat
      real constc,lapsec
      real sigmab
C
C-----------------------------------------------------------------------
C
      if(idlflag .eq. 1) then
C
C-----------------------------------------------------------------------
C
C Held/Suarez IDEALIZED physics algorithm:
C
C   Held, I. M., and M. J. Suarez, 1994: A proposal for the
C   intercomparison of the dynamical cores of atmospheric general
C   circulation models.
C   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
C
C-----------------------------------------------------------------------
C
C Add idealized radiative heating rates to temperature tendency
C
        efoldf =  1.
        efolda = 40.
        efolds =  4.
        sigmab =  0.7
C
        cossq   = coslat*coslat
        cossqsq = cossq*cossq
        sinsq   = sin(clat)*sin(clat)
        trefc   = 315. - 60.*sinsq
        onemsig = 1. - sigmab
C
        ka = 1./(86400.*efolda)
        ks = 1./(86400.*efolds)
C
        do k=1,plev
          if(etamid(k) .gt. sigmab) then
            kt = ka + (ks - ka)*cossqsq*(etamid(k) - sigmab)/onemsig
            tmp   = kt/(1.+ ztodt*kt)
            do i=1,plon
              trefa = (trefc - 10.*cossq*log((pmidm1(i,k)/ps0)))
     $                                       *(pmidm1(i,k)/ps0)**cappa
              trefa    = max(200.,trefa)
              qrs(i,k) = (trefa - tm1(i,k))*tmp
              t2 (i,k) = qrs(i,k)
            end do
          else
            tmp   = ka/(1.+ ztodt*ka)
            do i=1,plon
              trefa = (trefc - 10.*cossq*log((pmidm1(i,k)/ps0)))
     $                                       *(pmidm1(i,k)/ps0)**cappa
              trefa    = max(200.,trefa)
              qrs(i,k) = (trefa - tm1(i,k))*tmp
              t2 (i,k) = qrs(i,k)
            end do
          endif
        end do
C
C Add diffusion near the surface for the wind fields
C
        do k=1,plev
          do i=1,plond
            duv(i,k) = 0.
            dvv(i,k) = 0.
          end do
        end do
        do i=1,plond
          taux(i) = 0.
          tauy(i) = 0.
        end do
C
        kf = 1./(86400.*efoldf)
C
        do k=1,plev
          if(etamid(k) .gt. sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,plon
              duv(i,k) = tmp*um1(i,k)
              dvv(i,k) = tmp*vm1(i,k)
              fu(i,k)  = fu(i,k) + duv(i,k)
              fv(i,k)  = fv(i,k) + dvv(i,k)
            end do
          endif
        end do
      elseif(idlflag .eq. 2) then
C
C-----------------------------------------------------------------------
C
C Modified Held/Suarez IDEALIZED physics algorithm
C (modified with Williamson stratosphere):
C
C   Williamson, D. L., J. G. Olson and B. A. Boville, 1998: A comparison
C   of semi--Lagrangian and Eulerian tropical climate simulations.
C   Mon. Wea. Rev., vol 126, pp. 1001-1012.
C
C-----------------------------------------------------------------------
C
C Add idealized radiative heating rates to temperature tendency
C
        efoldf  =  1.
        efolda  = 40.
        efoldaa = 40.
        efolds  =  4.
        sigmab  =  0.7
C
        cossq   = coslat*coslat
        cossqsq = cossq*cossq
        sinsq   = sin(clat)*sin(clat)
        trefc   = 315. - 60.*sinsq
        onemsig = 1. - sigmab
C
        ka  = 1./(86400.*efolda)
        kaa = 1./(86400.*efoldaa)
        ks  = 1./(86400.*efolds)
C
        pi     = 4.*atan(1.)
        phi0   = 60.*pi/180.
        dphi0  = 15.*pi/180.
        a0     = 2.65/dphi0
        aeq    = 10000.
        apole  = 200.
        lapsew = -3.345e-03
        constw = rair*lapsew/gravit
        lapsec =  2.00e-03
        constc = rair*lapsec/gravit
        acoslat = abs(acos(coslat))
        p0strat = aeq - (aeq - apole)*0.5*
     $            (1. + tanh(a0*(acoslat - phi0)))
        do k=1,plev
          if(etamid(k) .gt. sigmab) then
            kt = ka + (ks - ka)*cossqsq*(etamid(k) - sigmab)/onemsig
            do i=1,plon
              tmp   = kt/(1.+ ztodt*kt)
              trefa = (trefc - 10.*cossq*log((pmidm1(i,k)/ps0)))
     $                                       *(pmidm1(i,k)/ps0)**cappa
              trefa   = max(200.,trefa)
              if(pmidm1(i,k) .lt. 10000.) then
                trefa = 200.*((pmidm1(i,k)/10000.))**constc
                tmp   = kaa/(1.+ ztodt*kaa)
              endif
              if(pmidm1(i,k) .lt. p0strat) then
                trefa = trefa + 
     $                  200.*( ((pmidm1(i,k)/p0strat))**constw - 1. )
                tmp   = kaa/(1.+ ztodt*kaa)
              endif
              qrs(i,k) = (trefa - tm1(i,k))*tmp
              t2 (i,k) = qrs(i,k)
            end do
          else
            do i=1,plon
              tmp   = ka/(1.+ ztodt*ka)
              trefa = (trefc - 10.*cossq*log((pmidm1(i,k)/ps0)))
     $                                       *(pmidm1(i,k)/ps0)**cappa
              trefa   = max(200.,trefa)
              if(pmidm1(i,k) .lt. 10000.) then
                trefa = 200.*((pmidm1(i,k)/10000.))**constc
                tmp   = kaa/(1.+ ztodt*kaa)
              endif
              if(pmidm1(i,k) .lt. p0strat) then
                trefa = trefa + 
     $                  200.*( ((pmidm1(i,k)/p0strat))**constw - 1. )
                tmp   = kaa/(1.+ ztodt*kaa)
              endif
              qrs(i,k) = (trefa - tm1(i,k))*tmp
              t2 (i,k) = qrs(i,k)
            end do
          endif
        end do
C
C Add diffusion near the surface for the wind fields
C
        do k=1,plev
          do i=1,plond
            duv(i,k) = 0.
            dvv(i,k) = 0.
          end do
        end do
        do i=1,plond
          taux(i) = 0.
          tauy(i) = 0.
        end do
C
        kf = 1./(86400.*efoldf)
C
        do k=1,plev
          if(etamid(k) .gt. sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,plon
              duv(i,k) = tmp*um1(i,k)
              dvv(i,k) = tmp*vm1(i,k)
              fu(i,k)  = fu(i,k) + duv(i,k)
              fv(i,k)  = fv(i,k) + dvv(i,k)
            end do
          endif
        end do
      elseif(idlflag .eq. 3) then
C
C---------------------------------------
C Boer/Denis IDEALIZED physics algorithm
C---------------------------------------
        write(6,*) 'TPHYSIDL: should not run model with "idlflag" set ',
     $             'to "3" (experimental code).'
        write(6,*) 'The valid options are 1 or 2.'
        call endrun
C
C Precompute frequently used variables  
C
        do k=1,plev
          do i=1,plon
            rpdel(i,k) = 1./pdelm1(i,k)
          end do
        end do
        do k=1,plev-1
          do i=1,plon
            rpdeli(i,k) = 1./(pmidm1(i,k+1) - pmidm1(i,k))
          end do
        end do
C
        do k=1,plev
          do i=1,plon
            pmlnm1(i,k)  = log(pmidm1(i,k))
          end do
        end do
        call virtem(tm1     ,qm1     ,zvir    ,tvm     )
        call zint  (pilnm1  ,pmlnm1  ,rair    ,gravit  ,tvm     ,
     $              zm      ,zi      )
C
C Add idealized radiative heating rates to temperature tendency
C
        do i = 1,plon
          coeff(i) = ps0/pintm1(i,plevp)
        end do
        do k=1,plev
          sum = 0.
          do i = 1,plon
            sum = sum + tm1(i,k)
          end do
          sum = sum/float(plon)
          do i = 1,plon
            qrs(i,k) = coeff(i)*h3z(k,lat) +
     $                              (1./(75.*86400.))*(t3z(k,lat) - sum)
            t2 (i,k) = qrs(i,k)
          end do
        end do
C
C Add diffusion near the surface for the wind fields
C
        do k=1,plev
          do i=1,plond
            duv(i,k) = 0.
            dvv(i,k) = 0.
          end do
        end do
        do i=1,plond
          taux(i) = 0.
          tauy(i) = 0.
        end do
C
C Vertical diffusion/pbl calculation
C
        call vdintr(hbuf    ,pmidm1  ,pintm1  ,psomc   ,thm     ,
     $              zm      ,rpdel   ,rpdeli  ,um1     ,vm1     ,
     $              tm1     ,taux    ,tauy    ,shflx   ,cflx    ,
     $              ztodt   ,lat     ,qm1     ,fu      ,fv      ,
     $              t2      ,dtv     ,duv     ,dvv     ,up1     ,
     $              vp1     ,pblh    ,tpert   ,qpert   )
      else
        write(6,*) 'TPHYSIDL: flag for choosing desired type of ideali',
     $             'zed physics ("idlflag") is set incorrectly.'
        write(6,*) 'The valid options are 1 or 2.'
        write(6,*) 'idlflag is currently set to: ',idlflag
        call endrun
      endif
C
C Archive idealized temperature tendency
C
      call outfld('QRS     ',qrs     ,plond   ,lat     ,hbuf    )
C
      return
C
      end
 
