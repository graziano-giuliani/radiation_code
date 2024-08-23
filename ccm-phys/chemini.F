#include <misc.h>
#include <params.h>
      subroutine chemini
C-----------------------------------------------------------------------
C Do initial read of time-variant chemical loss rate frequency dataset, containing
C loss rates as a function of latitude and pressure.  Determine the two
C consecutive months between which the current date lies.  
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
      include 'netcdf.inc'
C-----------------------------------------------------------------------
#include <pmgrid.h>
#include <pagrid.h>
#include <pspect.h>
C-----------------------------------------------------------------------
#if ( defined SPMD )
      integer msgtype, bufid, info
      parameter (msgtype=80000)
#endif
C-----------------------------------------------------------------------
#include <commss.h>
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comhyb.h>
C-----------------------------------------------------------------------
#include <comtim.h> 
C-----------------------------------------------------------------------
#include <comvmr.h>
C-----------------------------------------------------------------------
#include <comchem.h>
C-----------------------------------------------------------------------
#if ( defined SPMD )
#include <comspmd.h>
#endif
C-----------------------------------------------------------------------
C
C Local variables
C
      character*80 trcdata          ! tracer data filename
      character*80 locfn            ! local filename        
      character*(MAXNCNAM) lonnam   ! Name of longitude dimension 
      character*(MAXNCNAM) levnam   ! Name of level dimension 
      character*(MAXNCNAM) latnam   ! Name of latitude dimension 
      character*(MAXNCNAM) timenam  ! Name of time dimension 

      integer ncid_trc              ! netcdf id for input data file
      integer dateid            ! netcdf id for date variable
      integer secid             ! netcdf id for seconds variable
      integer lonid             ! netcdf id for longitude variable
      integer latid             ! netcdf id for latitude variable
      integer levid             ! netcdf id for level variable
      integer timid             ! netcdf id for time variable
      integer tch4id            ! netcdf id for ch4   loss rate
      integer tn2oid            ! netcdf id for n2o   loss rate
      integer tcfc11id          ! netcdf id for cfc11 loss rate
      integer tcfc12id          ! netcdf id for cfc12 loss rate
      integer lonsiz            ! size of longitude dimension on tracer dataset
      integer levsiz            ! size of level dimension on tracer dataset
      integer latsiz            ! size of latitude dimension on tracer dataset
      integer timsiz            ! size of time dimension on tracer dataset
      integer ret               ! return code
      integer j,n,k,nt          ! indices
      integer ki,ko,ji,jo       ! indices 

      real lato(plat)           ! ccm model latitudes (degrees)
      real zo(plev)             ! ccm model heights (m)
      real lati(ptrlat)         ! input data latitudes (degrees)
      real pin(ptrlev)          ! input data pressure values (mbars)
      real zi(ptrlev)           ! input data heights (m)

      real xch4i  (ptrlon,ptrlev,ptrlat,ptrtim) ! input ch4   loss rate coeff
      real xn2oi  (ptrlon,ptrlev,ptrlat,ptrtim) ! input n2o   loss rate coeff
      real xcfc11i(ptrlon,ptrlev,ptrlat,ptrtim) ! input cfc11 loss rate coeff
      real xcfc12i(ptrlon,ptrlev,ptrlat,ptrtim) ! input cfc12 loss rate coeff

      real xch4(ptrlat, ptrlev)   ! input ch4   loss rate coeff indices changed 
      real xn2o(ptrlat, ptrlev)   ! input n2o   loss rate coeff indices changed
      real xcfc11(ptrlat, ptrlev) ! input cfc11 loss rate coeff indices changed
      real xcfc12(ptrlat, ptrlev) ! input cfc12 loss rate coeff indices changed
      
      real xch4lv(ptrlat, plev)   ! input ch4   loss rate coeff interp to ccm levels
      real xn2olv(ptrlat, plev)   ! input n2o   loss rate coeff interp to ccm levels
      real xcfc11lv(ptrlat, plev) ! input cfc11 loss rate coeff interp to ccm levels
      real xcfc12lv(ptrlat, plev) ! input cfc12 loss rate coeff interp to ccm levels
C
C------------------------------Externals--------------------------------
C
      integer lencs               ! Position of right-most non-blank character
C
C-----------------------------------------------------------------------
C
C Initialize
C
      nm = 1
      np = 2
C
C SPMD: Master does all the work.  Sends needed info to slaves
C
      if (masterproc) then
C
C------------------------------------------------------------------------
C Read tracer data
C------------------------------------------------------------------------
C
C Obtain tracer dataset and open file
C
         trcdata = 'noaamisc.r8.nc'
         call getfil(trcdata, datadir, locfn)
         ret = NF_OPEN(locfn, 0, ncid_trc)
         if (ret.eq.0) then
            write(6,*)'CHEMINI: NCOPN returns id ',ncid_trc,
     $           ' for file ',locfn(1:lencs(locfn))
         else
            write(6,*)'CHEMINI: Failed to open file ',
     $           locfn(1:lencs(locfn))
            write(6,*) NF_STRERROR(RET)
            call endrun
         end if
C
C Get dimension info
C
         RET = NF_INQ_DIMID(NCID_TRC, 'lat', LATID)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_DIMID(NCID_TRC, 'lev', LEVID)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_DIMID(NCID_TRC, 'lon', LONID)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_DIMID(NCID_TRC, 'time', TIMID)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
C
C Check dimension info
C
         if (ptrlon.ne.1) then
            write(6,*)'CHEMINI: longitude dependence not implemented'
            call endrun
	 endif

         RET = NF_INQ_DIM(NCID_TRC, lonid, lonnam, lonsiz)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         if (lonsiz .ne. ptrlon) then
            write(6,*)'CHEMINI: lonsiz=',lonsiz,' must = ptrlon=',ptrlon
            call endrun
         end if

         RET = NF_INQ_DIM(NCID_TRC, levid, levnam, levsiz)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         if (levsiz .ne. ptrlev) then
            write(6,*)'CHEMINI: levsiz=',levsiz,' must = ptrlev=',ptrlev
            call endrun
         end if
         
         RET = NF_INQ_DIM(NCID_TRC, latid, latnam, latsiz)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         if (latsiz .ne. ptrlat) then
            write(6,*)'CHEMINI: latsiz=',latsiz,' must = ptrlat=',ptrlat
            call endrun
         end if
         
         RET = NF_INQ_DIM(NCID_TRC, timid, timenam, timsiz)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         if (timsiz .ne. ptrtim) then
            write(6,*)'CHEMINI: timsiz=',timsiz,' must = ptrtim=',ptrtim
            call endrun
         end if
C
C Determine necessary dimension and variable id's
C
         RET = NF_INQ_VARID(NCID_TRC, 'lat', latid)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'lev', levid)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'date', dateid)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'datesec', secid)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'TCH4', tch4id)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'TN2O', tn2oid)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'TCFC11', tcfc11id)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET = NF_INQ_VARID(NCID_TRC, 'TCFC12', tcfc12id)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
C
C Obtain entire date and sec variables. Assume that will always
C cycle over 12 month data.
C
         RET=NF_GET_VAR_INT(ncid_trc, dateid, date_tr)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET=NF_GET_VAR_INT(ncid_trc, secid, sec_tr)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         if (mod(date_tr(1),10000)/100.ne.1) then
            write(6,*)'(CHEMINI): error when cycling data: ',
     $           '1st month must be 1'
            call endrun
         end if
         if (mod(date_tr(ptrtim),10000)/100.ne.12) then
            write(6,*)'(CHEMINI): error when cycling data: ',
     $           'last month must be 12'
            call endrun
         end if
C
C Obtain input data latitude and level arrays.
C
         RET=NF_GET_VAR_DOUBLE(ncid_trc, latid, lati)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET=NF_GET_VAR_DOUBLE(ncid_trc, levid, pin)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
C
C Convert input pressure levels to height (m).
C First convert from millibars to pascals.
C
         do k=1,ptrlev
            pin(k) = pin(k)*100.
	    zi(k) = 7.0e3 * log (1.0e5 / pin(k))
         end do
C
C Convert approximate ccm pressure levels to height (m).
C
	 do k=1,plev
	    zo (k) = 7.0e3 * log (1.0e5 / hypm(k))
	 end do
C
C Convert ccm model latitudes to degrees.
C Input model latitudes already in degrees.
C
         do j=1,plat
            lato(j) = clat(j)*45./atan(1.)
         end do
C
C Obtain all time samples of tracer data.
C
!         strt4(1) = 1
!         strt4(2) = 1
!         strt4(3) = 1
! 	  strt4(4) = 1
!         cnt4(1)  = lonsiz
!         cnt4(2)  = levsiz
!         cnt4(3)  = latsiz
!         cnt4(4)  = timsiz
         RET=NF_GET_VAR_DOUBLE(ncid_trc, tch4id, xch4i)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET=NF_GET_VAR_DOUBLE(ncid_trc, tn2oid, xn2oi)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET=NF_GET_VAR_DOUBLE(ncid_trc, tcfc11id, xcfc11i)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
         RET=NF_GET_VAR_DOUBLE(ncid_trc, tcfc12id, xcfc12i)
         IF (RET .NE. NF_NOERR) CALL HANDLE_ERR('CHEMINI',RET)
C
C Close netcdf file
C
         ret = nf_close(ncid_trc)
C
C------------------------------------------------------------------------
C Interpolate tracer data to model grid
C------------------------------------------------------------------------
C
C Loop over all input times.
C
	 do nt = 1, ptrtim
C
C Remove longitude and time index and switch level and latitude indices
C for the loss coefficients.
C
	    do j=1,ptrlat
	       do k=1,ptrlev
		  xch4  (j,k) = xch4i  (1,k,j,nt)
		  xn2o  (j,k) = xn2oi  (1,k,j,nt)
		  xcfc11(j,k) = xcfc11i(1,k,j,nt)
		  xcfc12(j,k) = xcfc12i(1,k,j,nt)
	       end do
	    end do
C
C Interpolate input data to model levels.
C If the CCM level is outside the range of the input data (this
C can happen only in troposphere) put zero for every latitude.
C Otherwise determine the input data levels bounding the current 
C CCM level and interpolate.
C
            do ko=1,plev
               if (zo(ko) .lt. zi(ptrlev)) then
                  do j=1,ptrlat
                     xch4lv  (j,ko) = 0.0
                     xn2olv  (j,ko) = 0.0
                     xcfc11lv(j,ko) = 0.0
                     xcfc12lv(j,ko) = 0.0
                  end do
                  goto 50
               end if
               do ki=1,ptrlev-1
                  if (zo(ko).lt.zi(ki) .and. zo(ko).ge.zi(ki+1)) then
                     do j=1,ptrlat
                        xch4lv(j,ko) = xch4(j,ki)
     $                       + (xch4(j,ki+1) - xch4(j,ki))
     $                       / (zi(ki+1) - zi(ki))
     $                       * (zo(ko) - zi(ki))
                        xn2olv(j,ko) = xn2o(j,ki)
     $                       + (xn2o(j,ki+1) - xn2o(j,ki))
     $                       / (zi(ki+1) - zi(ki))
     $                       * (zo(ko) - zi(ki))
                        xcfc11lv(j,ko) = xcfc11(j,ki)
     $                       + (xcfc11(j,ki+1) - xcfc11(j,ki))
     $                       / (zi(ki+1) - zi(ki))
     $                       * (zo(ko) - zi(ki))
                        xcfc12lv(j,ko) = xcfc12(j,ki)
     $                       + (xcfc12(j,ki+1) - xcfc12(j,ki))
     $                       / (zi(ki+1) - zi(ki))
     $                       * (zo(ko) - zi(ki))
                     end do
                     goto 50
                  endif
               end do
               write (6,*) '(CHEMINI): Error in vertical interpolation'
               call endrun
   50          continue
            end do
C
C Interpolate input data to model latitudes.
C Determine the input data latitudes bounding the current CCM latitude and 
C interpolate. Use last value from input data if the ccm latitude is
C outside the range of the input data latitudes. 
C
            do jo=1,plat
               if (lato(jo) .le. lati(1)) then
                  do k = 1, plev
                     tch4i(jo,k,nt)   = xch4lv(1,k)
                     tn2oi(jo,k,nt)   = xn2olv(1,k)
                     tcfc11i(jo,k,nt) = xcfc11lv(1,k)
                     tcfc12i(jo,k,nt) = xcfc12lv(1,k)
                  end do
               else if (lato(jo) .ge. lati(ptrlat)) then
                  do k = 1, plev
                     tch4i(jo,k,nt)   = xch4lv(ptrlat,k)
                     tn2oi(jo,k,nt)   = xn2olv(ptrlat,k)
                     tcfc11i(jo,k,nt) = xcfc11lv(ptrlat,k)
                     tcfc12i(jo,k,nt) = xcfc12lv(ptrlat,k)
                  end do
               else
                  do ji=1,ptrlat-1
                     if ( (lato(jo).gt.lati(ji)) .and.
     $                    (lato(jo).le.lati(ji+1))) then
                        do k=1,plev
                           tch4i(jo,k,nt) = xch4lv(ji,k)
     $                          + (xch4lv(ji+1,k) - xch4lv(ji,k))
     $                          / (lati(ji+1)   -  lati(ji))
     $                          * (lato(jo) - lati(ji))
                           tn2oi(jo,k,nt) = xn2olv(ji,k)
     $                          + (xn2olv(ji+1,k) - xn2olv(ji,k))
     $                          / (lati(ji+1)   -  lati(ji))
     $                          * (lato(jo) - lati(ji))
                           tcfc11i(jo,k,nt) = xcfc11lv(ji,k)
     $                          + (xcfc11lv(ji+1,k) - xcfc11lv(ji,k))
     $                          / (lati(ji+1)   -  lati(ji))
     $                          * (lato(jo) - lati(ji))
                           tcfc12i(jo,k,nt) = xcfc12lv(ji,k)
     $                          + (xcfc12lv(ji+1,k) - xcfc12lv(ji,k))
     $                          / (lati(ji+1)   -  lati(ji))
     $                          * (lato(jo) - lati(ji))
                        end do
                        goto 90
                     endif
                  end do
               end if
               write (6,*)'(CHEMINI): Error in horizontal interpolation'
   90          continue
            end do
            
         end do                 ! end loop over time samples
C
C Initial time interpolation between December and January
C
         n = 12
         np1 = 1
         call bnddyi(date_tr(n  ), sec_tr(n  ), cdaytrm)
         call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
         if (calday.le.cdaytrp .or. calday.gt.cdaytrm) then
            do j=1,plat
               do k=1,plev
                  tch4m  (j,k,nm) = tch4i  (j,k,n)
                  tn2om  (j,k,nm) = tn2oi  (j,k,n)
                  tcfc11m(j,k,nm) = tcfc11i(j,k,n)
                  tcfc12m(j,k,nm) = tcfc12i(j,k,n)
                  tch4m  (j,k,np) = tch4i  (j,k,np1)
                  tn2om  (j,k,np) = tn2oi  (j,k,np1)
                  tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                  tcfc12m(j,k,np) = tcfc12i(j,k,np1)
               end do
            end do
            goto 10
         end if
C
C Initial normal interpolation between consecutive time slices.
C
         do n=1,timsiz-1
            np1 = n + 1
            call bnddyi(date_tr(n  ), sec_tr(n  ), cdaytrm)
            call bnddyi(date_tr(np1), sec_tr(np1), cdaytrp)
            if (calday.gt.cdaytrm .and. calday.le.cdaytrp) then
               do j=1,plat
                  do k=1,plev
                     tch4m  (j,k,nm) = tch4i  (j,k,n)
                     tn2om  (j,k,nm) = tn2oi  (j,k,n)
                     tcfc11m(j,k,nm) = tcfc11i(j,k,n)
                     tcfc12m(j,k,nm) = tcfc12i(j,k,n)
                     tch4m  (j,k,np) = tch4i  (j,k,np1)
                     tn2om  (j,k,np) = tn2oi  (j,k,np1)
                     tcfc11m(j,k,np) = tcfc11i(j,k,np1)
                     tcfc12m(j,k,np) = tcfc12i(j,k,np1)
                  end do
               end do
               goto 10
            end if
         end do
         write(6,*)'CHEMINI: Failed to find dates bracketing ',
     $        'mcdate, mcsec=', mcdate, mcsec
         call endrun
C
C Data positioned correctly
C
   10    continue

      endif                     ! end of masterproc

      return
      end
