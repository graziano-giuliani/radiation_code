#include <misc.h>
#include <params.h>
      subroutine vdintr(hbuf    ,pmidm1  ,pintm1  ,psomc   ,thm     ,
     $                  zm      ,rpdel   ,rpdeli  ,um1     ,vm1     ,
     $                  tm1     ,taux    ,tauy    ,shflx   ,cflx    ,
     $                  ztodt   ,lat     ,qm1     ,fu      ,fv      ,
     $                  t2      ,dtv     ,duv     ,dvv     ,up1     ,
     $                  vp1     ,pblh    ,tpert   ,qpert   )
C-----------------------------------------------------------------------
C
C interface routine for vertical diffusion and pbl scheme
C
C calling sequence:
C
C    vdinti        initializes vertical diffustion constants
C    pbinti        initializes pbl constants
C     .
C     .
C    vdintr        interface for vertical diffusion and pbl scheme
C      vdiff       performs vert diff and pbl
C        pbldif    boundary layer scheme
C        mvdiff    diffuse momentum
C        qvdiff    diffuse constituents
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: vdintr.F,v 1.1 1998/04/01 07:22:54 ccm Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comtrcnm.h>
C-----------------------------------------------------------------------
#include <comvd.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real hbuf(*)                ! history tape buffer,for calls to outfld
      real pmidm1(plond,plev)     ! midpoint pressures
      real pintm1(plond,plevp)    ! interface pressures
      real psomc(plond,plev)      ! (psm1/pmidm1)**cappa
      real thm(plond,plev)        ! potential temperature midpoints
      real zm(plond,plev)         ! midpoint geopotential height above sfc
      real rpdel(plond,plev)      ! 1./pdel (thickness between interfaces)
      real rpdeli(plond,plev)     ! 1./pdeli (thickness between midpoints)
      real um1(plond,plev)        ! u-wind input
      real vm1(plond,plev)        ! v-wind input
      real tm1(plond,plev)        ! temperature input
      real taux(plond)            ! x surface stress (n)
      real tauy(plond)            ! y surface stress (n)
      real shflx(plond)           ! surface sensible heat flux (w/m2)
      real cflx(plond,pcnst)      ! surface constituent flux (kg/m2/s)
      real ztodt                  ! 2 delta-t

      integer lat                 ! latitude index
C
C Input/output arguments
C
      real qm1(plond,plev,pcnst)  ! initial/final constituent field
      real fu(plond,plev)         ! cumulative u tendency
      real fv(plond,plev)         ! cumulative v tendency
      real t2(plond,plev)         ! cumulative T tendency
C
C Output arguments
C
      real dtv(plond,plev)        ! temperature tendency (heating)
      real duv(plond,plev)        ! u-wind tendency
      real dvv(plond,plev)        ! v-wind tendency
      real up1(plond,plev)        ! u-wind after vertical diffusion
      real vp1(plond,plev)        ! v-wind after vertical diffusion
      real pblh(plond)            ! planetary boundary layer height
      real tpert(plond)           ! convective temperature excess
      real qpert(plond,pcnst)     ! convective humidity and constituent excess
C
C---------------------------Local storage-------------------------------
C
      integer i,k,m               ! longitude,level,constituent indices

      real dqv(plond,plev,pcnst)  ! constituent diffusion tendency
      real ustar(plond)           ! surface friction velocity
      real kvh(plond,plevp)       ! diffusion coefficient for heat
      real kvm(plond,plevp)       ! diffusion coefficient for momentum
      real denom                  ! denominator of expression
      real qp1(plond,plev,pcnst)  ! constituents after vdiff
      real thp(plond,plev)        ! potential temperature after vdiff
      real cgs(plond,plevp)       ! counter-gradient star (cg/flux)
      real rztodt                 ! 1./ztodt
C
C-----------------------------------------------------------------------
C
C Call vertical diffusion code. No setup work is required.
C
      call vdiff(lat     ,um1     ,vm1     ,tm1     ,qm1     ,
     $           pmidm1  ,pintm1  ,rpdel   ,rpdeli  ,ztodt   ,
     $           thm     ,zm      ,taux    ,tauy    ,shflx   ,
     $           cflx    ,up1     ,vp1     ,thp     ,qp1     ,
     $           pblh    ,ustar   ,kvh     ,kvm     ,tpert   ,
     $           qpert   ,cgs     )
C
C Convert the diffused fields back to diffusion tendencies.
C Add the diffusion tendencies to the cummulative physics tendencies,
C except for constituents. The diffused values of the constituents
C replace the input values.
C
      rztodt = 1./ztodt
      do k=ntopfl,plev
        do i=1,plon
          duv(i,k) = (up1(i,k) - um1(i,k))*rztodt
          dvv(i,k) = (vp1(i,k) - vm1(i,k))*rztodt
          fu(i,k)  = fu(i,k) + duv(i,k)
          fv(i,k)  = fv(i,k) + dvv(i,k)
          denom    = cpair*(1. + cpvir*qm1(i,k,1))
          dtv(i,k) = -(um1(i,k)*duv(i,k) + vm1(i,k)*dvv(i,k))/denom +
     $                (thp(i,k)/psomc(i,k) - tm1(i,k))*rztodt
          t2(i,k)  = t2(i,k) + dtv(i,k)
        end do
        do m=1,pcnst
          do i=1,plon
            dqv(i,k,m) = (qp1(i,k,m) - qm1(i,k,m))*rztodt
            qm1(i,k,m) = qp1(i,k,m)
          end do
        end do
      end do
C
C Save the pbl output variables on the history file. The other
C vertical diffusion variables are saved later, after Rayleigh 
C friction tendencies have been added.
C
      call outfld('PBLH    ',pblh ,plond,lat,hbuf)
      call outfld('TPERT   ',tpert,plond,lat,hbuf)
      call outfld('QPERT   ',qpert,plond,lat,hbuf)
      call outfld('USTAR   ',ustar,plond,lat,hbuf)
C
C If you want to uncomment the following outfld calls you must modify
C the argument list to vdiff
C
C      call outfld ('CGH     ',cgh(1,2)  ,plond,nrow,hbuf)
C      call outfld ('CGQ     ',cgq(1,2,1),plond,nrow,hbuf)
      call outfld('KVH     ',kvh(1,2),plond,lat,hbuf)
C      call outfld('KVM     ',kvm(1,2),plond,lat,hbuf)
      call outfld ('CGS     ',cgs(1,2)  ,plond,lat,hbuf)
      call outfld(vdiffnam(1),dqv(1,1,1),plond,lat,hbuf)
c      do m=2,pcnst
c         call outfld(vdiffnam(m),dqv(1,1,m),plond,lat,hbuf)
c      end do
C
      return
      end
 
