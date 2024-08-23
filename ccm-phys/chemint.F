#include <misc.h>
#include <params.h>
	subroutine chemint
C-----------------------------------------------------------------------
C Time interpolate chemical loss rates to current time, reading
C in new monthly data if necessary
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
#include <pagrid.h>
#include <pspect.h>
C-----------------------------------------------------------------------
#if ( defined SPMD )
      integer info
#endif
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comchem.h>
C-----------------------------------------------------------------------
#if ( defined SPMD )
#include <comspmd.h>
#include <mpif.h>
#endif
C-----------------------------------------------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
C
C Local workspace
C
      integer ntmp           ! temporary
      integer j,k            ! indices

      real fact1, fact2      ! time interpolation factors
      real deltat            ! time (days) between interpolating ozone data

C-----------------------------------------------------------------------
C
C SPMD: Master does all the work.  Sends needed info to slaves
C
      if ( masterproc) then
C
C If model time is past current forward timeslice, obtain the next
C timeslice for time interpolation.  Messy logic is for 
C interpolation between December and January (np1.eq.1).
C
         if (calday .gt. cdaytrp .and. .not.
     $        (np1.eq.1 .and. calday.gt.cdaytrm)) then
            np1 = mod(np1,12) + 1
            if (np1.gt.ptrtim) then
               write(6,*)'CHEMINT: Attempt to access bad month'
               call endrun
            end if
            cdaytrm = cdaytrp
            call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
            if (np1.eq.1 .or. calday.le.cdaytrp) then
               ntmp = nm
               nm   = np
               np   = ntmp
	       do j=1,plat
		  do k=1,plev
		     tch4m  (j,k,np) = tch4i  (j,k,np1)
		     tn2om  (j,k,np) = tn2oi  (j,k,np1)
		     tcfc11m(j,k,np) = tcfc11i(j,k,np1)
		     tcfc12m(j,k,np) = tcfc12i(j,k,np1)
		  end do
	       end do
            else
               write(6,*)'CHEMINT: Input data for date',date_tr(np1),
     &              ' sec ',sec_tr(np1), 'does not exceed model date',
     &              mcdate,' sec ',mcsec,' Stopping.'
               call endrun
            end if
         end if
C
C Determine factors for time interpolation.
C
         if (np1.eq.1) then     ! Dec-Jan interpolation
            deltat = cdaytrp + 365. - cdaytrm
            if (calday.gt.cdaytrp) then ! We're in December
               fact1 = (cdaytrp + 365. - calday)/deltat
               fact2 = (calday - cdaytrm)/deltat
            else                ! We're in January
               fact1 = (cdaytrp - calday)/deltat
               fact2 = (calday + 365. - cdaytrm)/deltat
            end if
         else                   ! Non Dec-Jan interpolation
            deltat = cdaytrp - cdaytrm
            fact1 = (cdaytrp - calday)/deltat
            fact2 = (calday - cdaytrm)/deltat
         end if
C
C Check sanity of time interpolation factors to within 32-bit roundoff.
C
         if (abs(fact1+fact2-1.).gt.1.e-6 .or.
     $        fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or.
     $        fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
            write(6,*)'CHEMINT: Bad fact1 and/or fact2=',fact1,fact2
            call endrun
         end if
C
C Do time interpolation
C
	 do j=1,plat
	    do k=1,plev
	       tch4(j,k)   = tch4m  (j,k,nm)*fact1 +
     $                       tch4m  (j,k,np)*fact2
	       tn2o(j,k)   = tn2om  (j,k,nm)*fact1 +
     $		             tn2om  (j,k,np)*fact2
	       tcfc11(j,k) = tcfc11m(j,k,nm)*fact1 +
     $                       tcfc11m(j,k,np)*fact2
	       tcfc12(j,k)  =tcfc12m(j,k,nm)*fact1 +
     $                       tcfc12m(j,k,np)*fact2
	    end do
	 end do
	end if			! end of if-masterproc
#if ( defined SPMD )
	call MPI_BCAST (tch4  , plev*plat,
     $       REALTYPE,
     $       0,MPI_COMM_WORLD,info)
	call MPI_BCAST (tn2o  , plev*plat,
     $       REALTYPE,
     $       0,MPI_COMM_WORLD,info)
         call MPI_BCAST (tcfc11, plev*plat,
     $       REALTYPE,
     $       0,MPI_COMM_WORLD,info)
         call MPI_BCAST (tcfc12, plev*plat,
     $        REALTYPE,
     $        0,MPI_COMM_WORLD,info)
#endif
      
      return
      end
 
