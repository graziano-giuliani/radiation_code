#include <misc.h>
#include <params.h>
      subroutine tphysac(lat     ,hbuf    ,ztodt   ,pmidm1  ,pintm1  ,
     $                   pdelm1  ,pilnm1  ,t2      ,fu      ,fv      ,
     $                   um1     ,vm1     ,tm1     ,qm1     ,psomc   ,
     $                   thm     ,zm      ,pblh    ,qpert   ,qrs     ,
     $                   qrl     ,tpert   ,clat    ,shflx   ,taux    ,
     $                   tauy    ,cflx    ,sgh     ,lhflx   ,oro     ,
     $                   ts      ,psm1    ,phis    )
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
c
c $Id: tphysac.F,v 1.1.26.1 1999/01/05 22:39:05 erik Exp $
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
#include <comindx.h>
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
      real t2(plond,plev)         ! Temperature tendency
      real fu(plond,plev)         ! Accumulated u tendency (from tphys)
      real fv(plond,plev)         ! Accumulated v tendency (from tphys)
      real um1(plond,plev)        ! U horizontal wind component
      real vm1(plond,plev)        ! V horizontal wind component
      real tm1(plond,plev)        ! Temperature
      real qm1(plond,plev,pcnst+pnats)! Specific humidity and constituents
      real psm1(plond)            ! Surface pressure
      real psomc(plond,plev)      ! (psm1/pmidm1)**cappa
      real thm(plond,plev)        ! Potential temperature 
      real zm(plond,plev)         ! Height of midpoints (above surface) 
      real qrs(plond,plev)        ! Shortwave heating rate
      real qrl(plond,plev)        ! Longwave heating rate
      real clat                   ! Current latitude(radians)
      real oro(plond)             ! Land/ocean/sea ice flag
      real ts(plond)              ! Land/ocean/sea surface temperature
      real phis(plond)            ! Surface geopotential
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
      real sgh(plond)             ! Std. deviation of orography for gwd
      real lhflx(plond)           ! Latent heat flux (w/m^2)
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
      real srcsnk(plond,plev,pcnst) ! Src/sink for spec. hum. & constit.
      real taugtx(plond)          ! Gravity wave c=0 stress at model top (zonal)
      real taugty(plond)          ! Gravity wave c=0 stress at model top (meridional)
      real taugsx(plond)          ! Gravity wave c=0 stress at surface (zonal)
      real taugsy(plond)          ! Gravity wave c=0 stress at surface (meridional)
      real tor(plond)             ! Torque
      real mom(plond)             ! Angular momementum
      integer i,k                 ! Longitude, level indices
      integer m                   ! Constituent index
      integer npts                ! Number of longitude points
      integer indx(plond)         ! Longitude index array
      integer nptsexc             ! Number of excess points
      integer indxexc(plond)      ! Longitude index array of excess points
C
C-----------------------------------------------------------------------
C
C Initialize temperature tendency: heating rate from TPHYSBC
C
      do k=1,plev
         do i=1,plon
            t2(i,k) = qrs(i,k) + qrl(i,k)
         end do
      end do
C
C Initialize source/sink array to zero for all constituents.
C
      call resetr(srcsnk,plndlv*pcnst,0.)
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
C Check if latent heat flux exceeds the total moisture content of the 
C lowest model layer, thereby creating negative moisture.
C
      npts = plon
      do i=1,plon
        indx(i) = i
      end do

      call qneg4('TPHYSAC '   ,lat           ,indx  ,npts    ,ztodt   ,
     $           qm1(1,plev,1),rpdel(1,plev) ,shflx ,lhflx   ,cflx    ,
     $           indxexc ,nptsexc)
C
C Advected tracers
C
      if ( pcnst.ge.2 ) then
c
C Advected test tracers:
C trace_test1 or trace_test2 or trace_test3:
C   add radon surface flux (from WCRP) and radon decay
C trace_test3 only: 
C   zero strat tracer in bottom layer.
C o ixtrct   = index of radon advected tracer
C o ixtrct+1 = index of conserved unit tracer
C o ixtrct+2 = index of ozone-like tracer
C
         if ( trace_test1 .or. trace_test2 .or. trace_test3 ) then
            call rnsfwcrp( clat, oro, cflx(1,ixtrct) )
            call rndecay( qm1(1,1,ixtrct), ztodt, srcsnk(1,1,ixtrct) )
         endif
         if ( trace_test3) then
            call resetr( qm1(1,plev,ixtrct+2), plon, 0. )
         endif
C
C Advected greenhouse trace gases:
C   o ixtrcg   => index of n2o
C   o ixtrcg+1 => index of ch4
C   o ixtrcg+2 => index of cfc11
C   o ixtrcg+3 => index of cfc12
C
         if (trace_gas) then
            call chem (qm1, srcsnk, cflx, lat, ztodt)
            call outfld ('N2OSRC  ',srcsnk(1,1,ixtrcg  ),plond,lat,hbuf)
            call outfld ('CH4SRC  ',srcsnk(1,1,ixtrcg+1),plond,lat,hbuf)
            call outfld ('CFC11SRC',srcsnk(1,1,ixtrcg+2),plond,lat,hbuf)
            call outfld ('CFC12SRC',srcsnk(1,1,ixtrcg+3),plond,lat,hbuf)
         end if
C
      end if
C
C Add contribution of src/snk tendency to constituent fields
C (including moisture).
C
      do m=1,pcnst
         do k=1,plev
            do i=1,plon
               qm1(i,k,m) = qm1(i,k,m) + ztodt*srcsnk(i,k,m)
            end do
         end do
      end do
C
C Vertical diffusion/pbl calculation
C
      call vdintr(hbuf    ,pmidm1  ,pintm1  ,psomc   ,thm     ,
     $            zm      ,rpdel   ,rpdeli  ,um1     ,vm1     ,
     $            tm1     ,taux    ,tauy    ,shflx   ,cflx    ,
     $            ztodt   ,lat     ,qm1     ,fu      ,fv      ,
     $            t2      ,dtv     ,duv     ,dvv     ,up1     ,
     $            vp1     ,pblh    ,tpert   ,qpert   )
C
C Gravity wave drag
C
      call mgwintr(um1    ,vm1     ,tm1     ,qm1     ,sgh     ,
     $            pmidm1  ,pintm1  ,pdelm1  ,rpdel   ,pilnm1  ,
     $            zm      ,pblh    ,ztodt   ,hbuf    ,lat     ,
     $            clat    ,fu      ,fv      ,t2      ,dtv     ,
     $            taugsx  ,taugsy  ,taugtx  ,taugty   )
C
C Calculate and output total torque
C
      call torque( um1    ,psm1   ,pdelm1  ,phis    ,clat    ,
     $             taux   ,taugsx ,taugtx  ,tor     ,mom     )
      call outfld ('TORQUE  ',tor,plond,lat,hbuf)
      call outfld ('MOMEN   ',mom,plond,lat,hbuf)
C
C Dump diffusive tendencies to history tape buffer
C
      call outfld ('DTV     ',dtv,plond,lat,hbuf)
      call outfld ('DUV     ',duv,plond,lat,hbuf)
      call outfld ('DVV     ',dvv,plond,lat,hbuf)
C
      return
C
      end
 
