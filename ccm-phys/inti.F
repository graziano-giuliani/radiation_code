#include <misc.h>
#include <params.h>
      subroutine inti
C-----------------------------------------------------------------------
C
C This routine is supposed to 
C Set constants and call initialization procedures for time independent
C physics/dynamics routines
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
C
C-----------------------------------------------------------------------
c
c $Id: inti.F,v 1.2.10.1 1998/12/18 16:50:33 erik Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <pspect.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <commap.h>
C-----------------------------------------------------------------------
#include <comqfl.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C---------------------------Local workspace-----------------------------
C
      real clrh2o            ! Ratio of latent heat to gas constant (h2o)
C
C-----------------------------------------------------------------------
C
C Initialize constants 
C
      clrh2o = latvap/rh2o   ! Passed to several initialization routines
c
c Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
c
      tmass0 = 98222./gravit
      if(idealdyn) tmass0 = 100000./gravit
      if (masterproc) write(6,800) tmassf,tmass0,qmassf
  800 format('INTI: MASS OF INITIAL DATA BEFORE CORRECTION = '
     $        ,1p,e20.10,/,' DRY MASS WILL BE HELD = ',e20.10,/,
     $        ' MASS OF MOISTURE AFTER REMOVAL OF NEGATIVES = ',e20.10)
C
C Call time independent initialization routines for parameterizations.
C
      call mgwinti(cpair   ,cpwv    ,gravit  ,rair    )
c
      call vdinti (cpwv    ,cpair   ,gravit  ,rair    )
c
      call tsinti (tmelt   ,latvap  ,rair    ,stebol  ,dtime   )
c
      call radini (gravit  ,cpair   ,epsilo  ,stebol  )
c
      call cainti (cappa   ,cpair   ,epsilo  ,gravit  ,latvap  ,
     $             rhoh2o  ,clrh2o  )
c
      call esinti (epsilo  ,latvap  ,latice  ,rh2o    ,cpair   )
c
      call mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o  )
c
      call hdinti (rearth  ,dtime   )
c
      call btempini

      return
      end
 
