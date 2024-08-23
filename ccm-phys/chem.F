#include <misc.h>
#include <params.h>
	subroutine chem (qm1, qtn, cflx, lat, dt)
C-----------------------------------------------------------------------
C Apply the interpolated chemical loss rates from the input data to 
C N2O, CH4, CFC11 and CFC12. Set the surface values to a constant.
C-----------------------------------------------------------------------
C
C $Id: chem.F,v 1.1 1998/04/01 07:20:51 ccm Exp $
C
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <comvmr.h>
C-----------------------------------------------------------------------
#include <comchem.h>
C-----------------------------------------------------------------------
#include <comindx.h>
C-----------------------------------------------------------------------
C       
C Arguments:
C       
        integer lat                ! south-north latitude index

        real qm1(plond,plev,pcnst) ! mass mixing ratios
        real qtn(plond,plev,pcnst) ! tendency of mass mixing ratios (output)
	real cflx(plond,pcnst)	   ! Surface constituent flux (kg/m^2/s)
        real dt		           ! time step  
C       
C Local variables
C       
	integer i,k	  	! loop indexes
	real xch4		! new methane mass mixing ratio
        real xn2o		! new nitrous oxide mass mixing ratio
        real xcfc11		! new cfc11 mass mixing ratio
        real xcfc12		! new cfc12 mass mixing ratio
C       
C-----------------------------------------------------------------------
C       
C Apply chemical rate coefficient using time split implicit method. The
C turn the new value back into a tendency. NOTE that water tendency is 
C twice methane tendency. Water is specific humidity which is in mass 
C mixing ratio units. Note that
C  o ixtrcg   => indx of n2o
C  o ixtrcg+1 => indx of ch4
C  o ixtrcg+2 => indx of cfc11
C  o ixtrcg+3 => indx of cfc12
C       
	do k=1,plev-2
	   do i=1,plon
	      xn2o       = qm1(i,k,ixtrcg  ) / (1. + tn2o  (lat,k) * dt)
	      qtn(i,k,ixtrcg  ) =(xn2o   - qm1(i,k,ixtrcg  )) / dt

	      xch4       = qm1(i,k,ixtrcg+1) / (1. + tch4  (lat,k) * dt)
	      qtn(i,k,ixtrcg+1) =(xch4   - qm1(i,k,ixtrcg+1)) / dt

	      qtn(i,k,1) = -2. * rh2och4 * qtn(i,k,ixtrcg+1)

	      xcfc11     = qm1(i,k,ixtrcg+2) / (1. + tcfc11(lat,k) * dt)
	      qtn(i,k,ixtrcg+2) =(xcfc11 - qm1(i,k,ixtrcg+2)) / dt

	      xcfc12     = qm1(i,k,ixtrcg+3) / (1. + tcfc12(lat,k) * dt)
	      qtn(i,k,ixtrcg+3) =(xcfc12 - qm1(i,k,ixtrcg+3)) / dt
	   end do
	end do
C       
C Set the "surface" tendencies (bottom 2 levels) to maintain specified 
C tropospheric concentrations. 
C       
	do k = plev-1, plev
	   do i=1,plon
	      qtn(i,k,ixtrcg  ) =((rmwn2o*n2ovmr) - qm1(i,k,ixtrcg  ))
     $       	   / dt 
	      qtn(i,k,ixtrcg+1) =((rmwch4*ch4vmr) - qm1(i,k,ixtrcg+1))
     $       	   / dt 
	      qtn(i,k,ixtrcg+2) =((rmwf11*f11vmr) - qm1(i,k,ixtrcg+2))
     $	           / dt 
	      qtn(i,k,ixtrcg+3) =((rmwf12*f12vmr) - qm1(i,k,ixtrcg+3))
     $             / dt 
	   end do
	end do
C
C For now set all tracer fluxes to 0
C
	do i=1,plon
	   cflx(i,ixtrcg  ) = 0.
	   cflx(i,ixtrcg+1) = 0.
	   cflx(i,ixtrcg+2) = 0.
	   cflx(i,ixtrcg+3) = 0.
	end do

	return
	end
