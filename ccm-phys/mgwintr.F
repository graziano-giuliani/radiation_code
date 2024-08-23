#include <misc.h>
#include <params.h>
      subroutine mgwintr(u       ,v       ,t        ,q       ,sgh     ,
     $                   pm      ,pi      ,dpm      ,rdpm    ,piln    ,
     $                   zm      ,pblh    ,dt       ,hbuf    ,lat     ,
     $                   rlat    ,ut      ,vt       ,tt      ,dtv     ,
     $                   tau0x   ,tau0y   ,tautx    ,tauty   )
C-----------------------------------------------------------------------
C
C Interface for multiple gravity wave drag parameterization.
C
C-------------------------Code History----------------------------------
C
C Original version:  B. Boville
C Standardized:      T. Acker, Feb 1996
C Reviewed:          B. Boville, April 1996
C 
C-----------------------------------------------------------------------
c
c $Id: mgwintr.F,v 1.1.26.1 1999/01/05 22:39:02 erik Exp $
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
      integer lat                       ! latitude index
      real u(plond,plev)                ! midpoint zonal wind
      real v(plond,plev)                ! midpoint meridional wind
      real t(plond,plev)                ! midpoint temperatures
      real q(plond,plev)                ! midpoint specific humidities
      real sgh(plond)                   ! standard deviation of orography
      real pm(plond,plev)               ! midpoint pressures
      real pi(plond,0:plev)             ! interface pressures
      real piln(plond,0:plev)           ! ln(interface pressures)
      real dpm(plond,plev)              ! midpt. delta p (pi(k)-pi(k-1))
      real rdpm(plond,plev)             ! 1. / (pi(k)-pi(k-1))
      real zm(plond,plev)               ! midpoint heights
      real pblh(plond)                  ! planetary boundary layer height
      real dt                           ! time step
      real dtv(plond,plev)              ! t tendency from vert. diffusion
      real hbuf(*)                      ! history tape buffer
      real rlat                         ! current latitude in radians
C
C Output variables
C
      real ut(plond,plev)               ! zonal wind tendency
      real vt(plond,plev)               ! meridional wind tendency
      real tt(plond,plev)               ! temperature tendency
      real tau0x(plond)                 ! c=0 sfc. stress (zonal)
      real tau0y(plond)                 ! c=0 sfc. stress (meridional)
      real tautx(plond)                 ! c=0 stress at model top (zonal)
      real tauty(plond)                 ! c=0 stress at model top (meridional)
C
C Local workspace
C
      integer i,k                       ! loop indexes
      real utgw(plond,plev)             ! zonal wind tendency
      real vtgw(plond,plev)             ! meridional wind tendency
      real ttgw(plond,plev)             ! temperature tendency
C
C-----------------------------------------------------------------------
C
C Initialize gravity wave drag tendencies to zero
C
      do k=1,plev
         do i=1,plon
            utgw(i,k) = 0.
            vtgw(i,k) = 0.
            ttgw(i,k) = 0.
         end do
      end do
C
C Determine wave sources and properties
C
      call mgwdrag(u       ,v       ,t       ,sgh     ,pm      ,
     $             pi      ,dpm     ,rdpm    ,piln    ,zm      ,
     $             pblh    ,dt      ,rlat    ,utgw    ,vtgw    ,
     $             tau0x   ,tau0y   ,tautx   ,tauty   )
C
      do k = 1, plev
         do i = 1, plon
C
C Compute the temperature tendency from energy conservation.
C
            ttgw(i,k) = -(u(i,k)*utgw(i,k) + v(i,k)*vtgw(i,k))
     $           /(cp*(1.0 + cpvir*q(i,k)))
C
C Add the tendencies to the model arrays.
C
            ut(i,k) = ut(i,k) + utgw(i,k)
            vt(i,k) = vt(i,k) + vtgw(i,k)
            tt(i,k) = tt(i,k) + ttgw(i,k)
            dtv(i,k)= dtv(i,k)+ ttgw(i,k)
         end do
      end do
C
C Write output fields to history file
C
c$$$      call outfld ('TTMGW   ', ttgw, plond, lat, hbuf)
      call outfld ('UTGW    ', utgw, plond, lat, hbuf)
      call outfld ('VTGW    ', vtgw, plond, lat, hbuf)
      call outfld ('TAUGWX  ', tau0x, plond, lat, hbuf)
      call outfld ('TAUGWY  ', tau0y, plond, lat, hbuf)
      call outfld ('TAUTGWX ', tauTx, plond, lat, hbuf)
      call outfld ('TAUTGWY ', tauty, plond, lat, hbuf)

      return
      end


 
