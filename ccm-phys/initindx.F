#include <misc.h>
#include <params.h>
      subroutine initindx
C-----------------------------------------------------------------------
C     
C Determines index values for tracers species.
C     
C---------------------------Code history--------------------------------
C     
C Original version:  CCM3
C Standardized:      M. Vertenstein, Aug. 1997
C     
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
#include <comindx.h>
C---------------------------Local variables-----------------------------
C      
      integer padv              ! index pointer to advected tracers
      integer pnad              ! index pointer to non-advected tracers
C
C-----------------------------------------------------------------------
C Determine starting indices
C-----------------------------------------------------------------------
C
C Initialize advected tracer ending index. Note that advected tracer
C index=1 is reserved for Q.
C
      padv = 1
      pnad = pcnst 
C
C Determine starting indices for advected and non-advected user tracers.
C
      if (nusr_adv .gt. 0) then
         ixuadv = padv + 1
         padv = ixuadv + nusr_adv - 1
      endif
      if (nusr_nad .gt. 0) then
         ixunad = pnad + 1
         pnad = ixunad + nusr_nad - 1
      endif
C
C Determine starting index for advected greenhouse gases.
C
      if (trace_gas) then 
         ixtrcg = padv  + 1
         padv = ixtrcg + 3
      endif
C
C Determine starting index for advected test tracers.
C
      if (trace_test1) then
         ixtrct = padv + 1
         padv = ixtrct
      elseif (trace_test2) then
         ixtrct = padv + 1
         padv = ixtrct + 1
      elseif (trace_test3) then
         ixtrct = padv + 1
         padv = ixtrct + 2
      endif
C
C Final consistency check
C
      if (padv .ne. pcnst) then
          write(6,*)'INITINDX: ending advected tracer index ',padv,
     $        ' not equal to pcnst = ',pcnst
         call endrun
      endif
      if (pnad .ne. pcnst+pnats) then
          write(6,*)'INITINDX:ending non-advected tracer index',pnad,
     $        ' not equal to pcnst+pnats = ',pcnst+pnats
         call endrun
      endif
C
      return
      end
