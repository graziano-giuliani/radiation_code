#include <misc.h>
#include <params.h>
      subroutine fmrgrid(qrs     ,qrl     ,qrsm    ,qrlm    )
C-----------------------------------------------------------------------
C
C Interpolate model arrays to radiation vertical grid.
C
C-------------------------Code History----------------------------------
C
C Original version:  
C Standardized:      L. Buja, Feb 1996
C Reviewed:          
C 
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
c
c $Id: fmrgrid.F,v 1.1 1998/04/01 07:21:45 ccm Exp $
c
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <ptrrgrid.h>
C-----------------------------Arguments---------------------------------
C
C Input (radiation grid)
C
      real qrs(plond,plevr)      ! Shortwave heating rate
      real qrl(plond,plevr)      ! Longwave heating rate
C
C Output (model grid)
C
      real qrsm(plond,plev)      ! Shortwave heating rate
      real qrlm(plond,plev)      ! Longwave heating rate
C
C-----------------------------------------------------------------------
C
C Code to interpolate goes here.  Do nothing could be coded as a memory
C transfer, but left out for efficiency considerations.
C
      return
      end








 
