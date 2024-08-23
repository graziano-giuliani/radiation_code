#include <misc.h>
#include <params.h>
      subroutine radctl(hbuf    ,xlat    ,coslat  ,lat     ,lwup    ,
     $                  pmid    ,pint    ,pmln    ,piln    ,t       ,
     $                  qm1     ,cld     ,effcld  ,clwp    ,coszrs  ,
     $                  asdir   ,asdif   ,aldir   ,aldif   ,fsns    ,
     $                  qrs     ,qrl     ,flwds   ,rel     ,rei     ,
     $                  fice    ,sols    ,soll    ,solsd   ,solld   ,
     $                  oro     )
C-----------------------------------------------------------------------
C
C Driver for radiation computation.
C
C Radiation uses cgs units, so conversions must be done from
C model fields to radiation fields.
C
C Calling sequence:
C
C     radinp      Converts units of model fields and computes ozone
C                 mixing ratio for solar scheme
C
C     radcsw      Performs solar computation
C       radalb    Computes surface albedos
C       radded    Computes delta-Eddington solution
C       radclr    Computes diagnostic clear sky fluxes
C
C     radclw      Performs longwave computation
C
C       radtpl    Computes path quantities
C       radems    Computes emissivity
C       radabs    Computes absorptivity
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C Modified:          B. Briegleb, March 1995 to add aerosol
C                    to shortwave code
C
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C
C-----------------------------------------------------------------------
c
c $Id: radctl.F,v 1.4.6.6 1999/01/28 22:37:44 erik Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <ptrrgrid.h>
#include <pagrid.h>
#include <pspect.h>
#include <pbtemp.h>
C------------------------------Commons----------------------------------
#include <comengy.h>
C-----------------------------------------------------------------------
#include <crdcon.h>
C-----------------------------------------------------------------------
#include <comhst.h>
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
#include <comsol.h>
C-----------------------------------------------------------------------
#include <combtnam.h>
C-----------------------------------------------------------------------
#include <comindx.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real hbuf(*)              ! history tape buffer,for calls to outfld
      integer lat               ! Latitude row index
      real lwup(plond)          ! Longwave up flux at surface
      real pmid(plond,plev)     ! Model level pressures
      real pint(plond,plevp)    ! Model interface pressures
      real pmln(plond,plev)     ! Natural log of pmid
      real rel(plond,plev)      ! liquid effective drop size (microns)
      real rei(plond,plev)      ! ice effective drop size (microns)
      real fice(plond,plev)     ! fractional ice content within cloud
      real piln(plond,plevp)    ! Natural log of pint
      real t(plond,plev)        ! Model level temperatures
      real qm1(plond,plev,pcnst+pnats) ! Specific humidity and tracers
      real cld(plond,plevp)     ! Fractional cloud cover
      real effcld(plond,plevp)  ! Effective fractional cloud cover
      real clwp(plond,plev)     ! Cloud liquid water path
      real coszrs(plond)        ! Cosine solar zenith angle
      real asdir(plond)         ! albedo shortwave direct
      real asdif(plond)         ! albedo shortwave diffuse
      real aldir(plond)         ! albedo longwave direct
      real aldif(plond)         ! albedo longwave diffuse
      real xlat                 ! current latitude(radians)
      real coslat               ! cosine latitude
      real oro(plond)           ! Land/ocean/sea ice flag
C
C Output solar arguments
C
      real fsns(plond)          ! Surface absorbed solar flux
      real sols(plond)          ! Downward solar rad onto surface (sw direct)
      real soll(plond)          ! Downward solar rad onto surface (lw direct)
      real solsd(plond)         ! Downward solar rad onto surface (sw diffuse)
      real solld(plond)         ! Downward solar rad onto surface (lw diffuse)
      real qrs(plond,plev)      ! Solar heating rate
C
C Output longwave arguments
C
      real qrl(plond,plev)      ! Longwave cooling rate
      real flwds(plond)         ! Surface down longwave flux
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index

      real solin(plond)         ! Solar incident flux
      real fsds(plond)          ! Flux Shortwave Downwelling Surface
      real fsnirt(plond)        ! Near-IR flux absorbed at toa
      real fsnrtc(plond)        ! Clear sky near-IR flux absorbed at toa
      real fsnirtsq(plond)      ! Near-IR flux absorbed at toa >= 0.7 microns
      real fsnt(plond)          ! Net column abs solar flux at model top
      real fsntoa(plond)        ! Net column abs solar flux at top of atmosphere
      real fsntc(plond)         ! Clear sky total column abs solar flux
      real fsnsc(plond)         ! Clear sky surface abs solar flux
      real flnt(plond)          ! Net outgoing lw flux at model top
      real flns(plond)          ! Srf longwave cooling (up-down) flux
      real flntc(plond)         ! Clear sky lw flux at model top
      real flnsc(plond)         ! Clear sky lw flux at srf (up-down)
      real fldsc(plond)         ! Clear sky downwelling lw flux at srf
      real fsdsc(plond)         ! Clear sky solar incident at surface
      real pbr(plond,plevr)     ! Model mid-level pressures (dynes/cm2)
      real pnm(plond,plevrp)    ! Model interface pressures (dynes/cm2)
      real o3vmr(plond,plevr)   ! Ozone volume mixing ratio
      real o3mmr(plond,plevr)   ! Ozone mass mixing ratio
      real plco2(plond,plevrp)  ! Prs weighted CO2 path
      real plh2o(plond,plevrp)  ! Prs weighted H2O path
      real tclrsf(plond,plevrp) ! Total clear sky fraction  level to space
      real eccf                 ! Earth/sun distance factor
      real n2o(plond,plev)      ! nitrous oxide mass mixing ratio
      real ch4(plond,plev)      ! methane mass mixing ratio
      real cfc11(plond,plev)    ! cfc11 mass mixing ratio
      real cfc12(plond,plev)    ! cfc12 mass mixing ratio
      real aermmr(plond,plevr)  ! level aerosol mass mixing ratio
      real rh(plond,plevr)      ! level relative humidity (fraction)
      real lwupcgs(plond)       ! Upward longwave flux in cgs units
      integer irow              ! latitude pair index
      real britemp(plond,pnf_msu)! Infared brightness temperature
      real tb_ir(plond,pnb_hirs)! microwave brightness temperature
      real ts(plond)            ! surface temperature
      real pintmb(plond,plevp)  ! Model interface pressures (hPa)
      integer k                 ! Level subscript index
C
C Declare local arrays to which model input arrays are interpolated here.
C Current default is none since radiation grid = model grid.
C
C--------------------------------------------------------------------------
C
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
C
C Interpolate model input arrays to radiation vertical grid.  Currently this 
C is a do-nothing routine because radiation grid = model grid.
C
      call torgrid(pmid    ,pint    ,pmln    ,piln    ,t       ,
     $             qm1(1,1,1),cld     ,effcld  ,clwp    ,
     $             pmid    ,pint    ,pmln    ,piln    ,t       ,
     $             qm1(1,1,1)  ,cld     ,effcld  ,clwp    )
C
C Interpolate ozone volume mixing ratio to model levels
C
      call radozn(lat     ,pmid    ,o3vmr   )
      call outfld('O3VMR   ',o3vmr ,plond, lat, hbuf)
C
C Set latitude dependent radiation input
C
      call radinp(pmid    ,pint    ,qm1(1,1,1)  ,cld     ,o3vmr   ,
     $            pbr     ,pnm     ,plco2   ,plh2o   ,tclrsf  ,
     $            eccf    ,o3mmr   )
C
C Solar radiation computation
C
      if (dosw) then
C
C Specify aerosol mass mixing ratio
C
        call aermix(pnm     ,aermmr  ,rh      )
        call radcsw(pnm     ,qm1(1,1,1) ,o3mmr,aermmr  ,rh      ,
     $              cld     ,clwp    ,rel     ,rei     ,fice    ,
     $              eccf    ,coszrs  ,asdir   ,asdif   ,aldir   ,
     $              aldif   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $              fsds    ,fsnsc   ,fsntc   ,sols    ,soll    ,
     $              solsd   ,solld   ,fsnirt  ,fsnrtc  ,fsnirtsq,
     $              scon    ,fsdsc   ,fsntoa  )
C
C Convert units of shortwave fields needed by rest of model from CGS to MKS
C
        do i=1,plon
          solin(i)    = solin(i)*1.e-3
          fsds(i)     = fsds(i)*1.e-3
          fsnirt(i)   = fsnirt(i)*1.e-3
          fsnrtc(i)   = fsnrtc(i)*1.e-3
          fsnirtsq(i) = fsnirtsq(i)*1.e-3
          fsnt(i)     = fsnt(i) *1.e-3
          fsntoa(i)   = fsntoa(i) *1.e-3
          fsns(i)     = fsns(i) *1.e-3
          fsntc(i)    = fsntc(i)*1.e-3
          fsnsc(i)    = fsnsc(i)*1.e-3
          fsdsc(i)    = fsdsc(i)*1.e-3
        end do
C
C Dump shortwave radiation information to history tape buffer (diagnostics)
C
        call outfld('SOLIN   ',solin   ,plond,lat,hbuf)
        call outfld('FSDS    ',fsds    ,plond,lat,hbuf)
        call outfld('FSNIRT  ',fsnirt  ,plond,lat,hbuf)
        call outfld('FSNRTC  ',fsnrtc  ,plond,lat,hbuf)
        call outfld('FSNIRTSQ',fsnirtsq,plond,lat,hbuf)
        call outfld('FSNT    ',fsnt    ,plond,lat,hbuf)
        call outfld('FSNTOA  ',fsntoa  ,plond,lat,hbuf)
        call outfld('FSNS    ',fsns    ,plond,lat,hbuf)
        call outfld('FSNTC   ',fsntc   ,plond,lat,hbuf)
        call outfld('FSNSC   ',fsnsc   ,plond,lat,hbuf)
        call outfld('FSDSC   ',fsdsc   ,plond,lat,hbuf)
        call outfld('SOLS    ',sols    ,plond,lat,hbuf)
        call outfld('SOLL    ',soll    ,plond,lat,hbuf) 
        call outfld('SOLSD   ',solsd   ,plond,lat,hbuf) 
        call outfld('SOLLD   ',solld   ,plond,lat,hbuf) 
C
        call flxint(w(irow) ,fsnt   ,fsntil(lat))
        call flxint(w(irow) ,fsns   ,fsnsil(lat))
      end if
C
C Longwave radiation computation
C
      if (dolw) then
c
c Convert upward longwave flux units to CGS
c
        do i=1,plon
          lwupcgs(i) = lwup(i)*1000.
        end do
c
c Do longwave computation. If not implementing greenhouse gas code then
C first specify trace gas mixing ratios. If greenhouse gas code then:
C  o ixtrcg   => indx of advected n2o tracer
C  o ixtrcg+1 => indx of advected ch4 tracer
C  o ixtrcg+2 => indx of advected cfc11 tracer
C  o ixtrcg+3 => indx of advected cfc12 tracer
c
        if (trace_gas) then
           call radclw(lat  ,lwupcgs   ,t         ,qm1(1,1,1)  ,o3vmr,
     $                 pbr  ,pnm       ,pmln      ,piln        ,plco2,
     $                 plh2o,qm1(1,1,ixtrcg), qm1(1,1,ixtrcg+1),
     $                    qm1(1,1,ixtrcg+2), qm1(1,1,ixtrcg+3) ,
     $                 effcld,tclrsf  ,qrl        ,flns        ,flnt ,
     $                 flnsc ,flntc   ,flwds      ,fldsc       )

      else
         call trcmix(pmid, xlat, coslat, n2o, ch4, cfc11, cfc12)

         call radclw(lat     ,lwupcgs ,t       ,qm1(1,1,1)  ,o3vmr   ,
     $               pbr     ,pnm     ,pmln    ,piln        ,plco2   ,
     $               plh2o   ,n2o     ,ch4     ,cfc11       ,cfc12   ,
     $               effcld  ,tclrsf  ,qrl     ,flns        ,flnt    ,
     $               flnsc   ,flntc   ,flwds   ,fldsc   )
      endif
C
C Convert units of longwave fields needed by rest of model from CGS to MKS
C
        do i=1,plon
          flnt(i)    = flnt(i)*1.e-3
          flns(i)    = flns(i)*1.e-3
          flntc(i)   = flntc(i)*1.e-3
          flnsc(i)   = flnsc(i)*1.e-3
          flwds(i)   = flwds(i)*1.e-3
          fldsc(i)   = fldsc(i)*1.e-3
        end do
C
C Dump longwave radiation information to history tape buffer (diagnostics)
C
        call outfld('FLNT    ',flnt   ,plond,lat,hbuf)
        call outfld('FLNTC   ',flntc  ,plond,lat,hbuf)
        call outfld('FLNS    ',flns   ,plond,lat,hbuf)
        call outfld('FLNSC   ',flnsc  ,plond,lat,hbuf)
        call outfld('FLDSC   ',fldsc  ,plond,lat,hbuf)
C
        call flxint(w(irow) ,flnt     ,flntil(lat))
        call flxint(w(irow) ,flns     ,flnsil(lat))
      end if
C
C Brightness temperature calculation as a diagnostic to compare to TOV/MSU satellite data.
C When adding this into the main model development put an if around all
C this so it's only called when output.
C
      if ( dohirs .and. ((nstep .eq. nstepr) 
     $             .or. (mod(nstep-1,ihirsfq) .eq. 0)) )then
        do i=1,plon
          ts(i) = sqrt( sqrt( lwupcgs(i)/stebol ) )
          do k = 1, plev
            pintmb(i,k) = pint(i,k)/100.0   ! Convert from Pa to hPa
          end do
          pintmb(i,plevp) = pint(i,plevp)/100.0   ! Convert from Pa to hPa
        end do
        call hirsrtm( pintmb,     t,    qm1(1,1,1),
     $              o3mmr,    ts,    oro, tb_ir,       britemp  )
        do i = 1, pnb_hirs
          call outfld(hirsname(i),tb_ir(1,i),plon,lat,hbuf)
        end do
        do i = 1, pnf_msu
          call outfld(msuname(i),britemp(1,i),plon,lat,hbuf)
        end do
      end if
C
C Interpolate radiation output arrays to model vertical grid.  Currently this 
C is a do-nothing routine because radiation grid = model grid.
C
      call fmrgrid(qrs     ,qrl     ,
     $             qrs     ,qrl     )
C
      return
      end
 
