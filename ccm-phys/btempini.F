#include <misc.h>
#include <params.h>
      SUBROUTINE BTEMPINI
C
C Initialize values used for the brightness temperature calculation.
C
C------------------------------Parameters-------------------------------
#include <implicit.h>
#include <pmgrid.h>
C#include <pagrid.h>
#include <pbtemp.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <comvmr.h>
C-----------------------------------------------------------------------
#include <combtemp.h>
C-----------------------------------------------------------------------
#include <combtnam.h>
C------------------------------Local variables--------------------------
      INTEGER I        ! Loop index
C*******************************************************************************************
C     Set constants to values used in the model
C     Leave as they are for now.  Later if this gets committed back to the main 
C     development trunk we should use the values used in the rest of the radiation code.
C
C     GRAV = GRAVIT
C     R = RAIR
C     RCO2 = CO2VMR*RMWCO2
C
C*******************************************************************************************
      msuname(1)  = 'MSU_1   '
      msuname(2)  = 'MSU_2   '
      msuname(3)  = 'MSU_3   '
      msuname(4)  = 'MSU_4   '
      hirsname(1) = 'HIRS_2  '
      hirsname(2) = 'HIRS_4  '
      hirsname(3) = 'HIRS_6  '
      hirsname(4) = 'HIRS_8  '
      hirsname(5) = 'HIRS_10 '
      hirsname(6) = 'HIRS_11 '
      hirsname(7) = 'HIRS_12 '
      PRINT *, 'GRAV, GRAVIT: ', GRAV, GRAVIT
      PRINT *, 'R, RAIR: ', R, RAIR
      PRINT *, 'CO2VMR, RCO2: ', CO2VMR, RCO2/RMWCO2
C
      RETURN
      END

      subroutine hirsrtmfreq
C-----------------------------------------------------------------------
C
C Look through output fields to see if MSU or HIRS brightness
C Temperature fields are asked for.  If they are flag that 
C the HIRSRTM routine should be called and figure out the
C frequency it should be called at.
C
C---------------------------Code history--------------------------------
C
C Original version:  E. Kluzek, December 1998
C
C-----------------------------------------------------------------------
c
c $Id: btempini.F,v 1.1.2.2 1999/01/05 22:38:56 erik Exp $
c $Author: erik $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C-----------------------------------------------------------------------
#include <pbtemp.h>
C------------------------------Commons----------------------------------
#include <combtnam.h>
C-----------------------------------------------------------------------
#include <comhed.h>
C-----------------------------------------------------------------------
#include <comhst.h>
C------------------------------Local Variables--------------------------
      integer itape     ! Tape index
      integer ifld      ! Field index for current tape
      integer ifldtot   ! Field index for all tapes
      integer j         ! HIRS/MSU name index
      integer ihirsmn   ! Minimum frequency
      integer ihirsmx   ! Maximum frequency
      integer itens     ! Flag indicating if field is AVg,instant,max,min, etc.

      ihirsfq = 0
      dohirs  = .FALSE.
C
C Search MSU names
C
      ifldtot = 0
      do itape = 1, mtapes
        do ifld = 1, nflds(itape)
          ifldtot = ifldtot + 1
          do j = 1, pnf_msu
            if ( mcflds(1,ifldtot) .eq. msuname(j) )then
              itens = mflds(1,ifldtot)/10
C
C If Time averaged or max/min field  (anything but instantaneous)
C Need to call every-timestep
C
              if ( itens .ne. 0 ) then
                ihirsfq = 1
              else
C
C Otherwise call at the highest frequency
C (Make sure frequencies are compatable otherwise call each time step)
C
                if ( ihirsfq .eq. 0 )then
                  ihirsfq = nhtfrq(itape)
                else
                  ihirsmn = min(nhtfrq(itape),ihirsfq)
                  ihirsmx = max(nhtfrq(itape),ihirsfq)
                  if ( mod(ihirsmx,ihirsmn) .ne. 0 )then
                    ihirsfq = 1
                  else
                    ihirsfq = ihirsmn
                  end if
                end if
              end if
              dohirs = .TRUE.
            end if
          end do
        end do
      end do
C
C Search HIRS names
C
      ifldtot = 0
      do itape = 1, mtapes
        do ifld = 1, nflds(itape)
          ifldtot = ifldtot + 1
          do j = 1, pnb_hirs
            if ( mcflds(1,ifldtot) .eq. hirsname(j) )then
              itens = mflds(1,ifldtot)/10
C
C If Time averaged or max/min field  (anything but instantaneous)
C Need to call every-timestep
C
              if ( itens .ne. 0 ) then
                ihirsfq = 1
              else
C
C Otherwise call at the highest frequency
C (Make sure frequencies are compatable otherwise call each time step)
C
                if ( ihirsfq .eq. 0 )then
                  ihirsfq = nhtfrq(itape)
                else
                  ihirsmn = min(nhtfrq(itape),ihirsfq)
                  ihirsmx = max(nhtfrq(itape),ihirsfq)
                  if ( mod(ihirsmx,ihirsmn) .ne. 0 )then
                    ihirsfq = 1
                  else
                    ihirsfq = ihirsmn
                  end if
                end if
              end if
              dohirs = .TRUE.
            end if
          end do
        end do
      end do
      if ( dohirs )then
        print *, 'Call HIRSRTM every: ', ihirsfq, ' steps'
      else
        ihirsfq = 1
      end if
C
      return
      end
