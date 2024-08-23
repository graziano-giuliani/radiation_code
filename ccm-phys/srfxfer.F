#include <misc.h>
#include <params.h>
      subroutine srfxfer(lat     ,psm1tmp ,um1     ,vm1     ,tm1     ,
     $                   qm1     ,thm     ,zm      ,pmidm1  ,rhoh2o  ,
     $                   rpdel   ,fsns    ,precst  )
C-----------------------------------------------------------------------
C 
C Transfer atmospheric fields into common block /comsrf/
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath  Dec 1994
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------
c
c $Id: srfxfer.F,v 1.1 1998/04/01 07:22:32 ccm Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comsrf.h>
C-----------------------------------------------------------------------
#include <comtsc.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat               ! Latitude index
C
      real um1(plond)           ! Bottom level u wind
      real vm1(plond)           ! Bottom level v wind
      real tm1(plond)           ! Bottom level temperature
      real qm1(plond)           ! Bottom level specific humidity
      real thm(plond)           ! Bottom level potential temperature
      real zm(plond)            ! Bottom level height above surface
C
      real psm1tmp(plond)       ! Surface pressure
      real pmidm1(plond,plev)   ! Level pressures
      real rpdel(plond)         ! 1./(pint(k+1)-pint(k))
      real fsns(plond)          ! surface absorbed solar flux (W/m2)
      real precst(plond)        ! total snow precipitation (m/s)
      real rhoh2o               ! Density of water
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! Longitude index
C
C-----------------------------------------------------------------------
C
C Stuff global fluxes and state variables into common
C
      do i=1,plon
         tbot(i,lat) = tm1(i)
         thbot(i,lat) = thm(i)
         zbot(i,lat) = zm(i)
         ubot(i,lat) = um1(i)
         vbot(i,lat) = vm1(i)
         qbot(i,lat) = qm1(i)
         pbot(i,lat) = pmidm1(i,plev)
         psm1(i,lat) = psm1tmp(i)
         srfrpdel(i,lat) = rpdel(i)
#if ( defined COUP_SOM )
         absorb(i,lat)   = fsns(i)
         prcsnw(i,lat)   = precst(i)
#endif
#if ( defined COUP_CSM )
         rho(i,lat)   = pmidm1(i,plev)/(rair*tbot(i,lat))
         netsw(i,lat) = srfrad(i,lat) - flwds(i,lat)
#endif
      end do
C    
      return
C
      end






 
