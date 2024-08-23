#include <misc.h>
#include <params.h>
      subroutine qneg4(subnam  ,lat     ,indx    ,npts    ,ztodt   ,
     $                 qbot    ,srfrpdel,shflx   ,lhflx   ,qflx    ,
     $                 indxexc ,nptsexc )
C-----------------------------------------------------------------------
C
C Check if moisture flux into the ground is exceeding the total
C moisture content of the lowest model layer (creating negative moisture
C values).  If so, then subtract the excess from the moisture and
C latent heat fluxes and add it to the sensible heat flux.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson, Aug 1994
C Reviewed:          D. Williamson, March 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: qneg4.F,v 1.1 1998/04/01 07:22:09 ccm Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <comqmin.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      character*8 subnam        ! name of calling routine
C
      integer lat               ! latitude index
      integer indx(plond)       ! longitude index array
      integer npts              ! number of points in indx
C
      real ztodt                ! two times model timestep (2 delta-t)
      real qbot(plond)          ! moisture at lowest model level
      real srfrpdel(plond)      ! 1./(pint(K+1)-pint(K))
C
C Input/Output arguments
C
      real shflx(plond)         ! Surface sensible heat flux (J/m2/s)
      real lhflx(plond)         ! Surface latent   heat flux (J/m2/s)
      real qflx (plond)         ! surface water flux (kg/m^2/s)
C
C Output arguments
C
      integer indxexc(plond)    ! index array of points with excess flux
      integer nptsexc           ! number of points with excess flux
C
C---------------------------Local workspace-----------------------------
C
      integer i,ii              ! longitude indices
      integer iw                ! i index of worst violator
C
      real worst                ! biggest violator
      real excess(plond)        ! Excess downward sfc latent heat flux
C
C-----------------------------------------------------------------------
C
C Compute excess downward (negative) q flux compared to a theoretical
C maximum downward q flux.  The theoretical max is based upon the
C given moisture content of lowest level of the model atmosphere.
C
      nptsexc = 0
CDIR$ IVDEP
      do ii = 1,npts
        i = indx(ii)
        excess(i) = qflx(i) -
     $             (qmin(1) - qbot(i))/(ztodt*gravit*srfrpdel(i))
C
C If there is an excess downward (negative) q flux, then subtract 
C excess from "qflx" and "lhflx" and add to "shflx".
C
        if (excess(i).lt.0.) then
          nptsexc = nptsexc + 1
          indxexc(nptsexc) = i
          qflx (i) = qflx (i) - excess(i)
          lhflx(i) = lhflx(i) - excess(i)*latvap
          shflx(i) = shflx(i) + excess(i)*latvap
        end if
      end do
C
C Write out worst value if excess
C
      if (nptsexc.gt.0) then
        worst = 0.
        do ii=1,nptsexc
          i = indxexc(ii)
          if (excess(i).lt.worst) then
            worst = excess(i)
            iw = i
          end if
        end do
        write(6,9000) subnam,lat,nptsexc,worst,iw
      end if
C
      return
 9000 format(' QNEG4 WARNING from ',a8,', lat = ',i3,';',
     $       ' Max possible LH flx exceeded at ',i4,' points. ',
     $       ' Worst excess = ',1pe12.4,' at i = ',i4)
      end
 
