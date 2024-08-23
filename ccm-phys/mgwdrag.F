#include <misc.h>
#include <params.h>
      subroutine mgwdrag(u       ,v       ,t       ,sgh     ,pm      ,
     $                   pi      ,dpm     ,rdpm    ,piln    ,zm      ,
     $                   pblh    ,dt      ,rlat    ,ut      ,vt      ,
     $                   tau0x   ,tau0y   ,tautx   ,tauty   )
C-----------------------------------------------------------------------
C
C Driver for multiple gravity wave drag parameterization.
C 
C The parameterization is assumed to operate only where water vapor 
C concentrations are negligible in determining the density.
c
c Note that several parts of the code in this subroutine are not used in
c the standard configuration of ccm3. The standard runs use:
c     PGWV   = 0, only a single (stationary) wave is included;
c     TAULDV = 0, no low level flux divergence is invoked;
c     KBOT= PLEV, a surface source from orography is assumed.
C
C-------------------------Code History----------------------------------
C
C Original version:  B. Boville
C Standardized:      L. Buja,    Feb 1996
C Reviewed:          B. Boville, Feb 1996  
C
C-----------------------------------------------------------------------
c
c $Id: mgwdrag.F,v 1.1.26.1 1999/01/05 22:39:01 erik Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <mgw.h>
C-----------------------------------------------------------------------
C
C Input variables
C
      real u(plond,plev)       ! Midpoint zonal wind
      real v(plond,plev)       ! Midpoint meridional wind
      real t(plond,plev)       ! Midpoint temperatures
      real sgh(plond)          ! Standard deviation of orography
      real pm(plond,plev)      ! Midpoint pressures
      real pi(plond,0:plev)    ! Interface pressures
      real dpm(plond,plev)     ! Midpoint delta p (pi(k)-pi(k-1))
      real rdpm(plond,plev)    ! 1. / (pi(k)-pi(k-1))
      real piln(plond,0:plev)  ! ln(interface pressures)
      real zm(plond,plev)      ! Midpoint heights
      real pblh(plond)         ! Planetary boundary layer height
      real dt                  ! Time step
      real rlat                ! Current latitude in radians
C                              
C Output variables             
C                              
      real ut(plond,plev)      ! Zonal wind tendency
      real vt(plond,plev)      ! Meridional wind tendency
      real tau0x(plond)        ! c=0 sfc. stress (zonal)
      real tau0y(plond)        ! c=0 sfc. stress (meridional)
      real tautx(plond)        ! c=0 stress at model top (zonal)
      real tauty(plond)        ! c=0 stress at model top (meridional)
C
C--------------------------Local Variables------------------------------
C
      integer i,k,l            ! Loop indexes
      integer kldv(plond)      ! Top intfc index of low lvl stress div
      integer kldvmn           ! Min value of kldv
      integer ksrc(plond)      ! Index of top interface of source region
      integer ksrcmn           ! Min value of ksrc

      real dtdp 
      real n2                           ! Brunt-Vaisalla frequency squared
      real rhoi(plond,0:plev)           ! Interface density
      real ni(plond,0:plev)             ! Interface Brunt-Vaisalla frequency
      real ti(plond,0:plev)             ! Interface temperature
      real nm(plond,plev)               ! Midpoint Brunt-Vaisalla frequency
      real sghmax                       ! Max orographic sdv to use
      real tau(plond,-pgwv:pgwv,0:plev) ! Wave Reynolds stress
      real ubi(plond,0:plev)            ! Projection of wind at interfaces
      real ubm(plond,plev)              ! Projection of wind at midpoints
      real ubt(plond,plev)              ! ubar tendency
      real ubtl                         ! ubar tendency from wave l
      real xv(plond)                    ! Unit vectors of source wind (x)
      real yv(plond)                    ! Unit vectors of source wind (y)
      real ubtlsat                      ! Saturation tendency
      real d(plond)                     ! "total" diffusivity 
      real dsat(plond,-pgwv:pgwv)       ! Saturation diffusivity
      real dscal                        ! Fraction of dsat to use
      real mi                           ! Imaginary part of vertical wavenumber
      real tauback                      ! Background stress at c=0
      real taudmp                       ! Stress after damping
      real taumax(plond)                ! Max(tau) for any l
      real tauoro(plond)                ! c=0 stress from orography
      real tausat(plond,-pgwv:pgwv)     ! Saturation stress
      real ubmc(plond,-pgwv:pgwv)       ! (ub-c)
      real ubmc2                        ! (ub-c)**2
      real disp(plond)                  ! twice standard deviation of orography
C
C$$$      real zldv(plond)   ! Top of the low level stress divergence region
      real rdpldv(plond)  ! 1/dp across low level divergence region
      real nsrc(plond)    ! b-f frequency averaged over source region
      real psrc(plond)    ! Interface pressure at top of source region
      real rsrc(plond)    ! Density averaged over source region
      real usrc(plond)    ! u wind averaged over source region
      real vsrc(plond)    ! v wind averaged over source region
C
      real tauldv         ! Fraction of sfc stress which diverges
C                         !  between Sfc and zldv.
#if ( defined HP )
      real tmp	
#endif
C
C-----------------------------------------------------------------------
C
c move to mgwinti later
C
      tauldv = TAULDV
C
C Determine twice standard deviation of orography
C
      do i=1,plon
         disp(i) = 2.0*sgh(i)
      end do
C
C-----------------------------------------------------------------------
C Determine the interface densities and Brunt-Vaisala frequencies. 
C-----------------------------------------------------------------------
C
C The top interface values are calculated assuming an isothermal atmosphere 
C above the top level.
C
      k = 0
      do i = 1, plon
         ti(i,k) = t(i,k+1)
         rhoi(i,k) = pi(i,k) / (r*ti(i,k))
         ni(i,k) = sqrt (g*g / (cp*ti(i,k)))
      end do
      do k = 1, kbot-1
         do i = 1, plon
            ti(i,k) = 0.5 * (t(i,k) + t(i,k+1))
            rhoi(i,k) = pi(i,k) / (r*ti(i,k))
            dtdp = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
            n2 = g*g/ti(i,k) * (1./cp - rhoi(i,k)*dtdp)
#if ( defined HP )
            tmp = max(n2min,n2)
            ni(i,k) = sqrt (tmp)
#else
            ni(i,k) = sqrt (max (n2min, n2))
#endif
         end do
      end do
      if (kbot .lt. plev) then
         k = kbot
         do i = 1, plon
            ti(i,k) = 0.5 * (t(i,k) + t(i,k+1))
            rhoi(i,k) = pi(i,k) / (r*ti(i,k))
            dtdp = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
            n2 = g*g/ti(i,k) * (1./cp - rhoi(i,k)*dtdp)
#if ( defined HP )
            tmp = max(n2min,n2)
            ni(i,k) = sqrt (tmp)
#else
            ni(i,k) = sqrt (max (n2min, n2))
#endif
         end do
      else
         k = kbot
         do i = 1, plon
            ti(i,k) = t(i,k)
            rhoi(i,k) = pi(i,k) / (r*ti(i,k))
            ni(i,k) = ni(i,k-1)
         end do
      end if         
C
C-----------------------------------------------------------------------
C Determine the midpoint Brunt-Vaisala frequencies.
C-----------------------------------------------------------------------
C
      do k=1,kbot
         do i=1,plon
            nm(i,k) = 0.5 * (ni(i,k-1) + ni(i,k))
         end do
      end do
C
C-----------------------------------------------------------------------
C Determine the source layer wind and unit vectors, then project winds.
C-----------------------------------------------------------------------
C
C If the orographic term is not included, just use the bottom interface
C wind speed and direction (unit vector).
C
      if (kbot .lt. plev) then
C
C The following "k=kbot" assignment line is to fool subscript checking
C on the SUN compiler
C
         k = kbot
         do i = 1, plon
            ksrc(i) = k
            kldv(i) = k
            usrc(i) = 0.5*(u(i,k+1)+u(i,k))
            vsrc(i) = 0.5*(v(i,k+1)+v(i,k))
            ubi(i,kbot) = sqrt (usrc(i)**2 + vsrc(i)**2)
            xv(i) = usrc(i) / ubi(i,k)
            yv(i) = vsrc(i) / ubi(i,k)
         end do
C
C If the orographic term is included, average over the depth of the
C orographic standard deviation. Reduces to the bottom midpoint values
C when disp=0, such as over ocean.
C Also determine the depth of the low level stress divergence region, as
C the max of the boundary layer depth and the source region depth. This
C can be done here if the stress magnitude does not determine the depth,
C otherwise it must be done below.
C
      else
         k = plev
         do i = 1, plon
            ksrc(i) = k-1
            kldv(i) = k-1
            psrc(i) = pi(i,k-1)
            rsrc(i) = pm(i,k)/(r*t(i,k)) * dpm(i,k)
            usrc(i) = u(i,k) * dpm(i,k)
            vsrc(i) = v(i,k) * dpm(i,k)
            nsrc(i) = nm(i,k)* dpm(i,k)
         end do
         do k = plev-1, plev/2, -1
            do i = 1, plon
               if (disp(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
                  ksrc(i) = k-1
                  kldv(i) = k-1
                  psrc(i) = pi(i,k-1)
                  rsrc(i) = rsrc(i)
     $                 + pm(i,k) / (r*t(i,k))* dpm(i,k)
                  usrc(i) = usrc(i) + u(i,k) * dpm(i,k)
                  vsrc(i) = vsrc(i) + v(i,k) * dpm(i,k)
                  nsrc(i) = nsrc(i) + nm(i,k)* dpm(i,k)
               elseif (pblh(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
                  kldv(i) = k-1
               end if
            end do
         end do
         do i = 1, plon
            rsrc(i) = rsrc(i) / (pi(i,plev) - psrc(i))
            usrc(i) = usrc(i) / (pi(i,plev) - psrc(i))
            vsrc(i) = vsrc(i) / (pi(i,plev) - psrc(i))
            nsrc(i) = nsrc(i) / (pi(i,plev) - psrc(i))

            ubi(i,plev) = sqrt (usrc(i)**2 + vsrc(i)**2)
            xv(i) = usrc(i) / ubi(i,plev)
            yv(i) = vsrc(i) / ubi(i,plev)
         end do
      end if
C
C Project the local wind at midpoints onto the source wind.
C
      do k = 1, kbot
         do i = 1, plon
            ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
         end do
      end do
C
C Compute the interface wind projection by averaging the midpoint winds.
C Use the top level wind at the top interface.
C
      do i = 1, plon
         ubi(i,0) = ubm(i,1)
      end do
      do k = 1, kbot-1
         do i = 1, plon
            ubi(i,k) = 0.5 * (ubm(i,k) + ubm(i,k+1))
         end do
      end do
C
C-----------------------------------------------------------------------
C Gravity wave sources
C-----------------------------------------------------------------------
C
C
C Determine the orographic c=0 source term following McFarlane (1987).
C Set the source top interface index to plev, if the orographic term is
C zero.
C
      do i = 1, plon
         tauoro(i) = 0.
      end do
      if (kbot .eq. plev) then
         do i = 1, plon
            if ((ubi(i,plev) .gt. orovmin) .and.
     $           (disp(i) .gt. orohmin)) then
               sghmax = orofc2 * (ubi(i,plev) / nsrc(i))**2
               tauoro(i) = oroeko2 * min(disp(i)**2, sghmax) *
     $              rsrc(i) * nsrc(i) * ubi(i,plev)
            else
               ksrc(i) = plev
               kldv(i) = plev
            end if
         end do
      end if
C
C Determine the background stress at c=0
C
#if ( defined HP )
      tmp=2.*rlat
      if (rlat .ge. 0.) then
         tauback = taunor * sin(tmp)
      else
         tauback = -tausou * sin(tmp)
      end if
      tmp=-(rlat/taueqw)**2
      tauback = tauback + taueq * exp(tmp)
#else
      if (rlat .ge. 0.) then
         tauback = taunor * sin(2.*rlat)
      else
         tauback = tausou * abs(sin(2.*rlat))
      end if
      tauback = tauback + taueq * exp(-(rlat/taueqw)**2)
#endif
      tauback = tauback * tauscal
C
C Set the phase speeds and wave numbers in the direction of the source
C wind.  Set the source stress magnitude (positive only, note that the
C sign of the stress is the same as (c-u).
C
      do l = 1, pgwv
         do i = 1, plon
            tau(i, l,kbot) = tauback * exp(-(c(l)/30.)**2)
            tau(i,-l,kbot) = tau(i, l,kbot)
         end do
      end do
      do i = 1, plon
         tau(i,0,kbot) = max (tauback, tauoro(i))
      end do
C
C Find the top interface of the low level stress divergence region
C according to the maximum depth of three criterion.
C 1. source region depth
C 2. planetary boundary layer depth
C 3. 10 * (u_*) / N where u_* is defined from the gravity wave stresss
C                             = sqrt(tau/rho) using source region values
C
C$$$      if (kbot .lt. plev) then
C$$$         do i = 1, plon
C$$$            kldv(i) = kbot
C$$$         end do
C$$$      else
C$$$         do i = 1, plon
C$$$            zldv(i) = max (pblh(i), disp(i)
C$$$            zldv(i) = max (zdlv(i),
C$$$     $           zldvcon * sqrt(tau(i,0,k)/rsrc(i)) / nsrc(i))
C$$$         end do
C$$$         kldv(i) = plev-1
C$$$         do k = plev-1, plev/2, -1
C$$$            do i = 1, plon
C$$$               if (zldv(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
C$$$                  kldv(i)  = k-1
C$$$            end do
C$$$         end do
C$$$      end if
C
C Determine the min value of kldv and ksrc for limiting later loops and
C the pressure at the top interface of the low level stress divergence
C region.
C
      ksrcmn = plev
      kldvmn = plev
      if (kbot .eq. plev) then
         do i = 1, plon
            ksrcmn = min(ksrcmn, ksrc(i))
            kldvmn = min(kldvmn, kldv(i))
            if (kldv(i) .ne. plev) then
               rdpldv(i) = 1. / (pi(i,kldv(i)) - pi(i,plev))
            end if
         end do
      end if
      if (tauldv .le. 0.) kldvmn = plev
C
C-----------------------------------------------------------------------
C Compute the stress profiles and diffusivities
C-----------------------------------------------------------------------
C
C Loop from bottom to top to get stress profiles      
C
      do k = kbot-1, 0, -1
C
C Determine the absolute value of the saturation stress and the
C   diffusivity for each wave.
C Define critical levels where the sign of (u-c) changes between
C   interfaces.
C
         do i = 1, plon
            d(i) = dback
         end do
         do l = -pgwv, pgwv
            do i = 1, plon
               ubmc(i,l) = ubi(i,k) - c(l)
#if ( defined HP )
	       tmp=ubmc(i,l)
               tausat(i,l) = abs (efkw * rhoi(i,k) * tmp*tmp*tmp
     $              / (2.*ni(i,k)) )
#else
               tausat(i,l) = abs (efkw * rhoi(i,k) * ubmc(i,l)**3
     $              / (2.*ni(i,k)) )
#endif
               if (tausat(i,l) .le. taumin)
     $              tausat(i,l) = 0.0
               if (ubmc(i,l) / (ubi(i,k+1) - c(l)) .le. 0.0)
     $              tausat(i,l) = 0.0
               dsat(i,l) = (ubmc(i,l) / ni(i,k))**2 *
     $              (efkw * ubmc(i,l)**2 / (rog * ti(i,k) * ni(i,k))
     $              - alpha(k))
               dscal = min (1.0, tau(i,l,k+1) / (tausat(i,l)+taumin))
               d(i) = max( d(i), dscal * dsat(i,l))
            end do
         end do
C
C Compute stress for each wave. The stress at this level is the min of
C the saturation stress and the stress at the level below reduced by
C damping.  The sign of the stress must be the same as at the level
C below.
C
         do l = -pgwv, pgwv
            do i = 1, plon
               ubmc2 = max(ubmc(i,l)**2, ubmc2mn)
               mi = ni(i,k) / (2. * kwv * ubmc2) *
     $              (alpha(k) + ni(i,k)**2/ubmc2 * d(i))
#if ( defined HP )
      tmp  = -2.*mi*rog*t(i,k+1)*(piln(i,k+1)-piln(i,k))
               taudmp = tau(i,l,k+1)
     $              * exp(tmp)
#else
               taudmp = tau(i,l,k+1)
     $              * exp(-2.*mi*rog*t(i,k+1)*(piln(i,k+1)-piln(i,k)))
#endif
               if (taudmp .le. taumin) taudmp = 0.
               tau(i,l,k) = min (taudmp, tausat(i,l))
            end do
         end do
C
C The orographic stress term does not change across the source region
C Note that k ge ksrcmn cannot occur without an orographic source term
C
         if (k .ge. ksrcmn) then
            do i = 1, plon
               if (k .ge. ksrc(i)) then
                  tau(i,0,k) = tau(i,0,plev) 
               end if
            end do
         end if
C
C Require that the orographic term decrease linearly (with pressure)
C within the low level stress divergence region. This supersedes the
C above requirment of constant stress within the source region.
C
         if (k .ge. kldvmn) then
            do i = 1, plon
               if (k .ge. kldv(i)) then
                  tau(i,0,k) = min (tau(i,0,k), tau(i,0,plev)  *
     $                 (1. - tauldv * (pi(i,k)-pi(i,plev)) * rdpldv(i)))
               end if
            end do
         end if
C
C Apply lower bounds to the stress if pgwv > 0.
C
         if (pgwv .ge. 1) then
C
C Determine the max value of tau for any l
C
            do i = 1, plon
               taumax(i) = tau(i,-pgwv,k)
            end do
            do l = -pgwv+1, pgwv
               do i = 1, plon
                  taumax(i) = max(taumax(i), tau(i,l,k))
               end do
            end do
            do i = 1, plon
               taumax(i) = mxrange * taumax(i)
            end do
C
C Set the min value of tau for each wave to the max of mxrange*taumax or
C mxasym*tau(-c)
C
            do l = 1, pgwv
               do i = 1, plon
                  tau(i, l,k) = max(tau(i, l,k), taumax(i))
                  tau(i, l,k) = max(tau(i, l,k), mxasym*tau(i,-l,k))
                  tau(i,-l,k) = max(tau(i,-l,k), taumax(i))
                  tau(i,-l,k) = max(tau(i,-l,k), mxasym*tau(i, l,k))
               end do
            end do
C
C The following "l=0"assignment line is to fool subscript checking on 
C the SUN compiler
C
            l = 0
            do i = 1, plon
               tau(i,l,k) = max(tau(i,l,k),
     $              mxasym * 0.5 * (tau(i,l-1,k) + tau(i,l+1,k)))
            end do

         end if

      end do
C
C Set the top interface stress to zero. This is a poor approximation to
C the effect of level above, but may be better than letting the stress
C pass out of the model top and never act at all.
C
C$$$      do l = -pgwv, pgwv
C$$$         do i = 1, plon
C$$$            tau(i,l,0) = 0.0
C$$$         end do
C$$$      end do
C
C-----------------------------------------------------------------------
C Compute the tendencies from the stress divergence.
C-----------------------------------------------------------------------
C
C Loop over levels from top to bottom
C
      do k = 1, kbot
C
C Accumulate the mean wind tendency over wavenumber.
C
         do i = 1, plon
            ubt (i,k) = 0.0
         end do
         do l = -pgwv, pgwv
            do i = 1, plon
C
C Determine the wind tendency including excess stress carried down from
C above.
C
               ubtl = g * (tau(i,l,k)-tau(i,l,k-1)) * rdpm(i,k)
C
C Require that the tendency be no larger than the analytic solution for
C a saturated region [proportional to (u-c)^3].
C
#if ( defined HP )
	       tmp=(c(l)-ubm(i,k))
               ubtlsat = efkw * abs(tmp*tmp*tmp)
#else
               ubtlsat = efkw * abs((c(l)-ubm(i,k))**3)
#endif
     $              /(2.*rog*t(i,k)*nm(i,k))
C$$$     $              /(2.*7000.     *nm(i,k))
C$$$     $              /(2.*8781.42   *nm(i,k))
               if (ubtl .gt. ubtlsat) ubtl = ubtlsat
C
C Apply tendency limits to maintain numerical stability.
C 1. du/dt < |c-u|/dt so u-c can't change sign: u^n+1 = u^n + du/dt * dt
C 2. du/dt < tndmax   so ridicuously large tendency are not permitted
C
               ubtl = min(ubtl, 0.5 * abs(c(l)-ubm(i,k)) / dt)
               ubtl = min(ubtl, tndmax)
C
C Accumulate the mean wind tendency over wavenumber.
C
               ubt (i,k) = ubt (i,k) + sign(ubtl, c(l)-ubm(i,k))
C
C Redetermine the effective stress on the interface below from the wind 
C tendency. If the wind tendency was limited above, then the new stress
C will be smaller than the old stress and will cause stress divergence in
C the next layer down. This has the effect of smoothing large stress 
C divergences downward while conserving total stress.
C
               tau(i,l,k) = tau(i,l,k-1) + ubtl * dpm(i,k) / g
            end do
         end do
C
C Project the mean wind tendency onto the components.
C
         do i = 1, plon
            ut(i,k) = ubt(i,k) * xv(i)
            vt(i,k) = ubt(i,k) * yv(i)
         end do
C
C End of level loop
C
      end do
C
C-----------------------------------------------------------------------
C Project the c=0 stress in the direction of the source wind for
C recording on the output file.
C-----------------------------------------------------------------------
C
      do i = 1, plon
         tau0x(i) = tau(i,0,kbot) * xv(i)
         tau0y(i) = tau(i,0,kbot) * yv(i)
      end do
      do i = 1, plon
         tautx(i) = tau(i,0,0) * xv(i)
         tauty(i) = tau(i,0,0) * yv(i)
      end do
C
      return
      end

 
