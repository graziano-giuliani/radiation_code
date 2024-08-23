#include <misc.h>
#include <params.h>
      subroutine mfinti(rair    ,cpair   ,gravit  ,latvap  ,rhowtr  )
C-----------------------------------------------------------------------
C
C Initialize moist convective mass flux procedure common block, cmfmca 
C 
C----------------------------Code History-------------------------------
C
C Original version:  J. Hack
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, G. Taylor, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: mfinti.F,v 1.1 1998/04/01 07:21:54 ccm Exp $
c
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comcmf.h>
C-----------------------------------------------------------------------
#include <comhyb.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real rair              ! gas constant for dry air
      real cpair             ! specific heat of dry air
      real gravit            ! acceleration due to gravity
      real latvap            ! latent heat of vaporization
      real rhowtr            ! density of liquid water (STP)
      integer k              ! vertical level index
C
C-----------------------------------------------------------------------
C
C Initialize physical constants for moist convective mass flux procedure
C
      cp     = cpair         ! specific heat of dry air     
      hlat   = latvap        ! latent heat of vaporization  
      grav   = gravit        ! gravitational constant       
      rgas   = rair          ! gas constant for dry air
      rhoh2o = rhowtr        ! density of liquid water (STP)
C
C Initialize free parameters for moist convective mass flux procedure
C
      c0     = 1.0e-4        ! rain water autoconversion coeff (1/m)
      dzmin  = 0.0           ! minimum cloud depth to precipitate (m)
      betamn = 0.10          ! minimum overshoot parameter
      cmftau = 3600.         ! characteristic adjustment time scale
c
c Limit convection to regions below 40 mb
c
      if (hypi(1) .ge. 4.e3) then
        limcnv = 1
      else
        do k=1,plev
          if (hypi(k).lt.4.e3 .and. hypi(k+1).ge.4.e3) then
            limcnv = k
            goto 10
          end if
        end do
        limcnv = plevp
      end if
   10 if (masterproc) then
        write(6,*)'MFINTI: Convection will be capped at intfc ',limcnv,
     $            ' which is ',hypi(limcnv),' pascals'
      end if
      tpmax  = 1.50          ! maximum acceptable t perturbation (deg C)
      shpmax = 1.50e-3       ! maximum acceptable q perturbation (g/g)
      rlxclm = .true.        ! logical variable to specify that relaxation
c                                time scale should applied to column as 
c                                opposed to triplets individually 
C
C Initialize miscellaneous (frequently used) constants
C
      rhlat  = 1.0/hlat      ! reciprocal latent heat of vaporization  
      rcp    = 1.0/cp        ! reciprocal specific heat of dry air     
      rgrav  = 1.0/grav      ! reciprocal gravitational constant       
C
C Initialize diagnostic location information for moist convection scheme
C
      iloc   = 1             ! longitude point for diagnostic info
      jloc   = 1             ! latitude  point for diagnostic info
      nsloc  = 1             ! nstep value at which to begin diagnostics
C
C Initialize other miscellaneous parameters
C
      tiny   = 1.0e-36       ! arbitrary small number (scalar transport)
      eps    = 1.0e-13       ! convergence criteria (machine dependent) 
C
      return
      end



 
