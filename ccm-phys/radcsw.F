#include <misc.h>
#include <params.h>
      subroutine radcsw(pint    ,h2ommr  ,o3mmr   ,aermmr  ,rh      ,
     $                  cld     ,clwp    ,rel     ,rei     ,fice    ,
     $                  eccf    ,coszrs  ,asdir   ,asdif   ,aldir   ,
     $                  aldif   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $                  fsds    ,fsnsc   ,fsntc   ,sols    ,soll    ,
     $                  solsd   ,solld   ,fsnirt  ,fsnrtc  ,fsnirtsq,
     $                  scon    ,fsdsc   ,fsntoa  )
C-----------------------------------------------------------------------
C
C Solar radiation code
C
C Basic method is Delta-Eddington as described in:
C
C    Briegleb, Bruce P., 1992: Delta-Eddington
C    Approximation for Solar Radiation in the NCAR Community Climate Model,
C    Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C Two changes to the basic method described above are: (1) the distinction
C between liquid and ice particle clouds, and (2) the addition of an
C aerosol with sulfate radiative properties.
C
C Divides solar spectrum into 18 intervals from 0.2-5.0 micro-meters.
C solar flux fractions specified for each interval. allows for
C seasonally and diurnally varying solar input.  Includes molecular,
C cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud, 
C and surface absorption. Computes delta-eddington reflections and
C transmissions assuming homogeneously mixed layers. Adds the layers 
C assuming scattering between layers to be isotropic, and distinguishes 
C direct solar beam from scattered radiation.
C
C Longitude loops are broken into 1 or 2 sections, so that only daylight
C (i.e. coszrs > 0) computations are done.
C
C Note that an extra layer above the model top layer is added.
C
C cgs units are used.
C
C Special diagnostic calculation of the clear sky surface and total column
C absorbed flux is also done for cloud forcing diagnostics.
C
C
C---------------------------Code history--------------------------------
C
C Modified March 1995 to add aerosols
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C Reviewed:          J. Kiehl, April 1996
C Reviewed:          B. Briegleb, May 1996
C 19th Band Added:   W. Collins March 1997
C
C-----------------------------------------------------------------------
c
c $Id: radcsw.F,v 1.1.28.4 1999/01/28 22:37:44 erik Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <prgrid.h>
C-----------------------------------------------------------------------
      integer nspint           ! Num of spctrl intervals across solar spectrum

      parameter ( nspint = 19 )
C-----------------------Constants for new bands--------------------------
      real*8
     $     V_RAYTAU_35,
     $     V_RAYTAU_64,
     $     V_ABO3_35,
     $     V_ABO3_64,
     $     V_KSA_35,
     $     V_KSA_64,
     $     V_GSA_35,
     $     V_GSA_64
      parameter(
     $     V_RAYTAU_35 = 0.155208, 
     $     V_RAYTAU_64 = 0.0392,   
     $     V_ABO3_35 = 2.4058030e+01,  
     $     V_ABO3_64 = 2.210e+01,  
     $     V_KSA_35 = 5.64884,       
     $     V_KSA_64 = 3.6771,       
     $     V_GSA_35 = .699326,      
     $     V_GSA_64 = .663642
     $     )
C------------------------------Commons----------------------------------
#include <crdcon.h>
C
C Input arguments
C
      real pint(plond,plevp)   ! Interface pressure
      real h2ommr(plond,plev)  ! Specific humidity (h2o mass mix ratio)
      real o3mmr(plond,plev)   ! Ozone mass mixing ratio
      real aermmr(plond,plev)  ! Aerosol mass mixing ratio
      real rh(plond,plev)      ! Relative humidity (fraction)
C
      real cld(plond,plevp)    ! Fractional cloud cover
      real clwp(plond,plev)    ! Layer liquid water path
      real rel(plond,plev)     ! Liquid effective drop size (microns)
      real rei(plond,plev)     ! Ice effective drop size (microns)
      real fice(plond,plev)    ! Fractional ice content within cloud
C
      real eccf                ! Eccentricity factor (1./earth-sun dist ** 2)
      real coszrs(plond)       ! Cosine solar zenith angle
      real asdir(plond)        ! 0.2-0.7 micro-meter srfc alb to direct rad
      real aldir(plond)        ! 0.7-5.0 micro-meter srfc alb to direct rad
      real asdif(plond)        ! 0.2-0.7 micro-meter srfc alb to diffuse rad
      real aldif(plond)        ! 0.7-5.0 micro-meter srfc alb to diffuse rad

      real scon                ! solar constant 
C
C Output arguments
C
      real solin(plond)        ! Incident solar flux
      real qrs(plond,plev)     ! Solar heating rate
      real fsns(plond)         ! Surface absorbed solar flux
      real fsnt(plond)         ! Total column absorbed solar flux
      real fsntoa(plond)       ! Total column absorbed solar flux at top of atmosphere
      real fsntoac(plond)      ! Clear sky total column absorbed solar flux at top of atmosphere
      real fsds(plond)         ! Flux Shortwave Downwelling Surface
      real fsdsc(plond)        ! Clear sky flux Shortwave incident at surface
C
      real fsnsc(plond)        ! Clear sky surface absorbed solar flux
      real fsntc(plond)        ! Clear sky total column absorbed solar flx
      real sols(plond)         ! Direct solar rad incident on surface (< 0.7)
      real soll(plond)         ! Direct solar rad incident on surface (>= 0.7)
      real solsd(plond)        ! Diffuse solar rad incident on surface (< 0.7)
      real solld(plond)        ! Diffuse solar rad incident on surface (>= 0.7)
      real fsnirt(plond)       ! Near-IR flux absorbed at toa
      real fsnrtc(plond)       ! Clear sky near-IR flux absorbed at toa
      real fsnirtsq(plond)     ! Near-IR flux absorbed at toa >= 0.7 microns
C
C------------------------------Externals--------------------------------
C
      integer   isrchfgt       ! Search for first array element > 0
      integer   isrchfle       ! Search for first array element < 0
C
C---------------------------Local variables-----------------------------
C
      integer ns               ! Spectral loop index
      integer i                ! Longitude loop index
      integer k                ! Level loop index
      integer n                ! Loop index for daylight
      integer nloop            ! Number of daylight loops
      integer is(2)            ! Daytime start indices
      integer ie(2)            ! Daytime end indices
      integer indxsl           ! Index for cloud particle properties
C
C A. Slingo's data for cloud particle radiative properties (from 'A GCM
C Parameterization for the Shortwave Properties of Water Clouds' JAS
C vol. 46 may 1989 pp 1419-1427)
C
      real abarl(4)            ! A coefficient for extinction optical depth
      real bbarl(4)            ! B coefficient for extinction optical depth
      real cbarl(4)            ! C coefficient for single particle scat albedo
      real dbarl(4)            ! D coefficient for single particle scat albedo
      real ebarl(4)            ! E coefficient for asymmetry parameter
      real fbarl(4)            ! F coefficient for asymmetry parameter

      save abarl, bbarl, cbarl, dbarl, ebarl, fbarl

      data abarl/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
      data bbarl/ 1.305    , 1.346    ,1.454    ,1.641    /
      data cbarl/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
      data dbarl/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
      data ebarl/ 0.829    , 0.794    ,0.754    ,0.826    /
      data fbarl/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/

      real abarli              ! A coefficient for current spectral interval
      real bbarli              ! B coefficient for current spectral interval
      real cbarli              ! C coefficient for current spectral interval
      real dbarli              ! D coefficient for current spectral interval
      real ebarli              ! E coefficient for current spectral interval
      real fbarli              ! F coefficient for current spectral interval
C
C Caution... A. Slingo recommends no less than 4.0 micro-meters nor
C greater than 20 micro-meters
C
c ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
c
      real abari(4)            ! a coefficient for extinction optical depth
      real bbari(4)            ! b coefficient for extinction optical depth
      real cbari(4)            ! c coefficient for single particle scat albedo
      real dbari(4)            ! d coefficient for single particle scat albedo
      real ebari(4)            ! e coefficient for asymmetry parameter
      real fbari(4)            ! f coefficient for asymmetry parameter

      save abari, bbari, cbari, dbari, ebari, fbari

      data abari/ 3.448e-03, 3.448e-03,3.448e-03,3.448e-03/
      data bbari/ 2.431    , 2.431    ,2.431    ,2.431    /
      data cbari/ 1.00e-05 , 1.10e-04 ,1.861e-02,.46658   /
      data dbari/ 0.0      , 1.405e-05,8.328e-04,2.05e-05 /
      data ebari/ 0.7661   , 0.7730   ,0.794    ,0.9595   /
      data fbari/ 5.851e-04, 5.665e-04,7.267e-04,1.076e-04/

      real abarii              ! A coefficient for current spectral interval
      real bbarii              ! B coefficient for current spectral interval
      real cbarii              ! C coefficient for current spectral interval
      real dbarii              ! D coefficient for current spectral interval
      real ebarii              ! E coefficient for current spectral interval
      real fbarii              ! F coefficient for current spectral interval
C
      real delta               ! Pressure (atmospheres) for stratos. h2o limit
      real o2mmr               ! O2 mass mixing ratio:

      save delta, o2mmr

      data delta /  1.70e-3 /
      data o2mmr / .23143 /

      real albdir(plond)   ! Current spc intrvl srf alb to direct rad
      real albdif(plond)   ! Current spc intrvl srf alb to diffuse rad
C
C Next series depends on spectral interval
C
      real frcsol(nspint)  ! Fraction of solar flux in each spectral interval
      real wavmin(nspint)  ! Min wavelength (micro-meters) of interval
      real wavmax(nspint)  ! Max wavelength (micro-meters) of interval
      real raytau(nspint)  ! Rayleigh scattering optical depth
      real abh2o(nspint)   ! Absorption coefficiant for h2o (cm2/g)
      real abo3 (nspint)   ! Absorption coefficiant for o3  (cm2/g)
      real abco2(nspint)   ! Absorption coefficiant for co2 (cm2/g)
      real abo2 (nspint)   ! Absorption coefficiant for o2  (cm2/g)
      real ph2o(nspint)    ! Weight of h2o in spectral interval
      real pco2(nspint)    ! Weight of co2 in spectral interval
      real po2 (nspint)    ! Weight of o2  in spectral interval
      real nirwgt(nspint)  ! Weights for intervals to simulate satellite filter
      real wgtint          ! Weight for specific spectral interval

      save frcsol ,wavmin ,wavmax ,raytau ,abh2o ,abo3 ,
     $     abco2  ,abo2   ,ph2o   ,pco2   ,po2   ,nirwgt

      data frcsol / .001488, .001389, .001290, .001686, .002877,
     $              .003869, .026336, .360739, .065392, .526861, 
     $              .526861, .526861, .526861, .526861, .526861, 
     $              .526861, .006239, .001834, .001834/
C
C weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans
C
      data nirwgt /  0.0,   0.0,   0.0,      0.0,   0.0,
     $               0.0,   0.0,   0.0, 0.320518,   1.0,  1.0,
     $               1.0,   1.0,   1.0,      1.0,   1.0,
     $               1.0,   1.0,   1.0 /

      data wavmin / .200,  .245,  .265,  .275,  .285,
     $              .295,  .305,  .350,  .640,  .700,  .701,
     $              .701,  .701,  .701,  .702,  .702,
     $             2.630, 4.160, 4.160/

      data wavmax / .245,  .265,  .275,  .285,  .295,
     $              .305,  .350,  .640,  .700, 5.000, 5.000,
     $             5.000, 5.000, 5.000, 5.000, 5.000,
     $             2.860, 4.550, 4.550/

      data raytau / 4.020, 2.180, 1.700, 1.450, 1.250,
     $              1.085, 0.730, V_RAYTAU_35, V_RAYTAU_64, 0.020, 
     $              .0001, .0001, .0001, .0001, .0001, .0001,
     $              .0001, .0001, .0001/
C
C Absorption coefficients
C
      data abh2o /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .002,    
     $                .035,     .377,    1.950,   9.400,  44.600, 
     $             190.000,     .000,     .000,    .000/

      data abo3  /
     $ 5.370e+04, 13.080e+04,  9.292e+04, 4.530e+04, 1.616e+04,
     $ 4.441e+03,  1.775e+02, V_ABO3_35, V_ABO3_64,      .000,      
     $      .000,   .000    ,   .000   ,   .000   ,      .000,
     $      .000,   .000    ,   .000   ,   .000    /

      data abco2  /    .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,    .000,    .000,    .000,    
     $                 .000,     .000,     .000,    .000,    .000, 
     $                 .000,     .094,     .196,   1.963/

      data abo2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,1.11e-05,6.69e-05,    
     $                .000,     .000,     .000,    .000,    .000,    
     $                .000,     .000,    .000,    .000/
C
C Spectral interval weights
C
      data ph2o  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .505,    
     $                .210,     .120,     .070,    .048,    .029,    
     $                .018,     .000,     .000,    .000/

      data pco2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .000,    
     $                .000,     .000,     .000,    .000,    .000,    
     $                .000,    1.000,     .640,    .360/

      data po2   /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,     .000,   1.000,   1.000,    
     $                .000,     .000,     .000,    .000,    .000,    
     $                .000,     .000,     .000,    .000/
C
C Diagnostic and accumulation arrays; note that sfltot, fswup, and
C fswdn are not used in the computation,but are retained for future use.
C
      real solflx(plond)         ! Solar flux in current interval
      real sfltot(plond)         ! Spectrally summed total solar flux
      real totfld(plond,0:plev)  ! Spectrally summed flux divergence
      real fswup(plond,0:plevp)  ! Spectrally summed up flux
      real fswdn(plond,0:plevp)  ! Spectrally summed down flux
C
C Cloud radiative property arrays
C
      real tauxcl(plond,0:plev)  ! water cloud extinction optical depth
      real tauxci(plond,0:plev)  ! ice cloud extinction optical depth
      real wcl(plond,0:plev)     ! liquid cloud single scattering albedo
      real gcl(plond,0:plev)     ! liquid cloud asymmetry parameter
      real fcl(plond,0:plev)     ! liquid cloud forward scattered fraction
      real wci(plond,0:plev)     ! ice cloud single scattering albedo
      real gci(plond,0:plev)     ! ice cloud asymmetry parameter
      real fci(plond,0:plev)     ! ice cloud forward scattered fraction
C
C Aerosol radiative property arrays
C
      real tauxar(plond,0:plev)  ! aerosol extinction optical depth
      real wa(plond,0:plev)      ! aerosol single scattering albedo
      real ga(plond,0:plev)      ! aerosol assymetry parameter
      real fa(plond,0:plev)      ! aerosol forward scattered fraction
      real tauaer(plond)         ! total column aerosol extinction
      real waer(plond)           ! aerosol single scattering albedo
      real gaer(plond)           ! aerosol asymmetry parameter
      real faer(plond)           ! aerosol forward scattering fraction
C
C Sulphate aerosol properties taken from:
C
C    Kiehl, J.T., B.P.Briegleb, 1993. The Relative Roles of Sulfate Aerosols
C    and Greenhouse Gases in Climate Forcing. Science, Vol. 260, pp. 311-314.
C
      real ksa(nspint)           ! aerosol spectral mass absorption coeff(m2/g)
      real wsa(nspint)           ! aerosol spectral single scattering albedo
      real gsa(nspint)           ! aerosol spectral asymmetry parameter
C
      data ksa /11.1163, 10.5472, 10.2468, 10.0392,  9.8292,
     $           9.6199,  9.0407,   V_KSA_35, V_KSA_64,  1.9169,  
     $           0.3780,  0.3780,  0.3780,  0.3780,  0.5704,  
     $           0.5704,  0.5704,  0.5704,  0.5704 /

      data wsa / .999999, .999999, .999999, .999999, .999999,
     $           .999999, .999999, .999999, .999999, .999991, 
     $           .989772, .989772, .989772, .989772, .847061, 
     $           .847061, .847061, .847061, .847061 /

      data gsa / .719161, .719012, .718453, .717820, .716997,
     $           .715974, .712743, V_GSA_35, V_GSA_64, .618115, 
     $           .485286, .485286, .485286, .485286, .295557, 
     $           .295557, .295557, .295557, .295557 /

C
C Other variables and arrays needed for aerosol:
C
      real rhfac                 ! multiplication factor for kaer
      real rhpc                  ! level relative humidity in %

      real a0                    ! constant in rh mult factor
      real a1                    ! constant in rh mult factor
      real a2                    ! constant in rh mult factor
      real a3                    ! constant in rh mult factor

      save a0,a1,a2,a3

      data a0 / -9.2906106183    /
      data a1 /  0.52570211505   /
      data a2 / -0.0089285760691 /
      data a3 /  5.0877212432e-05/
C
C Various arrays and other constants:
C
      real pflx(plond,0:plevp)   ! Interface press, including extra layer
      real zenfac(plond)         ! Square root of cos solar zenith angle
      real sqrco2                ! Square root of the co2 mass mixg ratio
      real tmp1                  ! Temporary constant array
      real tmp2                  ! Temporary constant array
      real pdel                  ! Pressure difference across layer
      real path                  ! Mass path of layer
      real ptop                  ! Lower interface pressure of extra layer
      real ptho2                 ! Used to compute mass path of o2
      real ptho3                 ! Used to compute mass path of o3
      real pthco2                ! Used to compute mass path of co2
      real pthh2o                ! Used to compute mass path of h2o
      real h2ostr                ! Inverse square root h2o mass mixing ratio
      real wavmid                ! Spectral interval middle wavelength
      real trayoslp              ! Rayleigh optical depth/standard pressure
      real tmp1l                 ! Temporary constant array
      real tmp2l                 ! Temporary constant array
      real tmp3l                 ! Temporary constant array
      real tmp1i                 ! Temporary constant array
      real tmp2i                 ! Temporary constant array
      real tmp3i                 ! Temporary constant array
      real rdenom                ! Multiple scattering term
      real psf                   ! Frac of solar flux in spect interval
      real gocp                  ! Gravity/cp
C
C Layer absorber amounts; note that 0 refers to the extra layer added
C above the top model layer
C
      real uh2o(plond,0:plev)    ! Layer absorber amount of h2o
      real uo3(plond,0:plev)     ! Layer absorber amount of  o3
      real uco2(plond,0:plev)    ! Layer absorber amount of co2
      real uo2(plond,0:plev)     ! Layer absorber amount of  o2
      real uaer(plond,0:plev)    ! Layer aerosol amount 
C
C Total column absorber amounts:
C
      real uth2o(plond)          ! Total column  absorber amount of  h2o
      real uto3(plond)           ! Total column  absorber amount of  o3
      real utco2(plond)          ! Total column  absorber amount of  co2
      real uto2(plond)           ! Total column  absorber amount of  o2
      real utaer(plond)          ! Total column  aerosol
C
C These arrays are defined for plev model layers; 0 refers to the extra
C layer on top:
C
      real rdir(plond,0:plev)    ! Layer reflectivity to direct rad
      real rdif(plond,0:plev)    ! Layer reflectivity to diffuse rad
      real tdir(plond,0:plev)    ! Layer transmission to direct rad
      real tdif(plond,0:plev)    ! Layer transmission to diffuse rad
      real explay(plond,0:plev)  ! Solar beam exp transmission for layer
      real flxdiv(plond,0:plev)  ! Flux divergence for layer
C
C These arrays are defined at model interfaces; 0 is the top of the
C extra layer above the model top; plevp is the earth surface:
C
      real rupdir(plond,0:plevp) ! Ref to dir rad for layers below
      real rupdif(plond,0:plevp) ! Ref to dif rad for layers below
      real rdndif(plond,0:plevp) ! Ref to dif rad for layers above
      real exptdn(plond,0:plevp) ! Solar beam exp down transm from top
      real tottrn(plond,0:plevp) ! Total transmission for layers above
      real fluxup(plond,0:plevp) ! Up   flux at model interface
      real fluxdn(plond,0:plevp) ! Down flux at model interface
C
C-----------------------------------------------------------------------
C
C Initialize output fields:
C
      do i=1, plon
        fsds(i)     = 0.0
        fsnirt(i)   = 0.0
        fsnrtc(i)   = 0.0
        fsnirtsq(i) = 0.0
        fsnt(i)     = 0.0
        fsntoa(i)   = 0.0
        fsntoac(i)  = 0.0
        fsns(i)     = 0.0
        solin(i)    = 0.0
        fsnsc(i)    = 0.0
        fsdsc(i)    = 0.0
        fsntc(i)    = 0.0
        sols(i)     = 0.0
        soll(i)     = 0.0
        solsd(i)    = 0.0
        solld(i)    = 0.0
      end do
      do k=1, plev
        do i=1, plon
          qrs(i,k) = 0.0
        end do
      end do
C
C Compute starting, ending daytime loop indices:
C
      nloop = 0
      is(1) = isrchfgt(plon,coszrs,1,0.0)
C
C If night everywhere, return:
C
      if(is(1).gt.plon) return
      ie(1) = isrchfle(plon-is(1),coszrs(is(1)+1),1,0.0) + is(1) - 1
      nloop = 1
C
C Possibly 2 daytime loops needed:
C
      if (ie(1).ne.plon) then
        is(2) = isrchfgt(plon-ie(1),coszrs(ie(1)+1),1,0.0) + ie(1)
        if(is(2).le.plon) then
          nloop = 2
          ie(2) = plon
        end if
      end if
C
C Define solar incident radiation and interface pressures:
C
      do n=1,nloop
        do i=is(n),ie(n)
          solin(i) = scon*eccf*coszrs(i)
          pflx(i,0) = 0.
        end do
      end do
      do k=1,plevp
        do n=1,nloop
          do i=is(n),ie(n)
            pflx(i,k) = pint(i,k)
          end do
        end do
      end do
C
C Compute optical paths:
C
      tmp1   = 0.5/(gravit*sslp)
      sqrco2 = sqrt(co2mmr)
      do n=1,nloop
        do i=is(n),ie(n)
          ptop      = pflx(i,1)
          ptho2     = o2mmr * ptop / gravit
          ptho3     = o3mmr(i,1) * ptop / gravit
          pthco2    = sqrco2 * (ptop / gravit)
          h2ostr    = sqrt( 1. / h2ommr(i,1) )
          zenfac(i) = sqrt(coszrs(i))
          pthh2o    = ptop**2*tmp1 + (ptop*rga)*(h2ostr*zenfac(i)*delta)
          uh2o(i,0) = h2ommr(i,1)*pthh2o
          uco2(i,0) = zenfac(i)*pthco2
          uo2 (i,0) = zenfac(i)*ptho2
          uo3 (i,0) = ptho3
          uaer(i,0) = 0.0
        end do
      end do
C
      tmp2 = delta/gravit
      do k=1,plev
        do n=1,nloop
          do i=is(n),ie(n)
            pdel   = pflx(i,k+1) - pflx(i,k)
            path   = pdel / gravit
            ptho2  = o2mmr * path
            ptho3  = o3mmr(i,k) * path
            pthco2 = sqrco2 * path
            h2ostr = sqrt(1.0/h2ommr(i,k))
            pthh2o = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 +
     $              pdel*h2ostr*zenfac(i)*tmp2
            uh2o(i,k) = h2ommr(i,k)*pthh2o
            uco2(i,k) = zenfac(i)*pthco2
            uo2 (i,k) = zenfac(i)*ptho2
            uo3 (i,k) = ptho3
C
C Adjust aerosol amount by relative humidity factor:
C
            if( rh(i,k) .gt. .90 ) then
              rhfac = 2.8
            else if (rh(i,k) .lt. .60 ) then
              rhfac = 1.0
            else
              rhpc  = 100. * rh(i,k)
              rhfac = (a0 + a1*rhpc + a2*rhpc**2 + a3*rhpc**3)
            endif
            uaer(i,k) = aermmr(i,k)*rhfac*path
          end do
        end do
      end do
C
C Compute column absorber amounts for the clear sky computation:
C
      do n=1,nloop
        do i=is(n),ie(n)
          uth2o(i) = 0.0
          uto3(i)  = 0.0
          utco2(i) = 0.0
          uto2(i)  = 0.0
          utaer(i) = 0.0
        end do
      end do
      do k=1,plev
        do n=1,nloop
          do i=is(n),ie(n)
            uth2o(i) = uth2o(i) + uh2o(i,k)
            uto3(i)  = uto3(i)  + uo3(i,k)
            utco2(i) = utco2(i) + uco2(i,k)
            uto2(i)  = uto2(i)  + uo2(i,k)
            utaer(i) = utaer(i) + uaer(i,k)
          end do
        end do
      end do
C
C Initialize spectrally integrated totals:
C
      do k=0,plev
        do i=1,plon
          totfld(i,k) = 0.0
          fswup (i,k) = 0.0
          fswdn (i,k) = 0.0
        end do
      end do
      do i=1,plon
        sfltot(i)       = 0.0
        fswup (i,plevp) = 0.0
        fswdn (i,plevp) = 0.0
      end do
C
C Set cloud properties for top (0) layer; so long as tauxcl is zero,
C there is no cloud above top of model; the other cloud properties
C are arbitrary:
C
      do n=1,nloop
        do i=is(n),ie(n)
          tauxcl(i,0) = 0.
          wcl(i,0)     = 0.999999
          gcl(i,0)     = 0.85
          fcl(i,0)     = 0.725
          tauxci(i,0) = 0.
          wci(i,0)     = 0.999999
          gci(i,0)     = 0.85
          fci(i,0)     = 0.725
C
C Aerosol 
C
          tauxar(i,0) = 0.
          wa(i,0)      = 0.925
          ga(i,0)      = 0.850
          fa(i,0)      = 0.7225
        end do
      end do
C
C Begin spectral loop
C
      do 100 ns=1,nspint
         wgtint = nirwgt(ns)
C
C Set index for cloud particle properties based on the wavelength,
C according to A. Slingo (1989) equations 1-3:
C Use index 1 (0.25 to 0.69 micrometers) for visible
C Use index 2 (0.69 - 1.19 micrometers) for near-infrared
C Use index 3 (1.19 to 2.38 micrometers) for near-infrared
C Use index 4 (2.38 to 4.00 micrometers) for near-infrared
C
C Note that the minimum wavelength is encoded (with .001, .002, .003)
C in order to specify the index appropriate for the near-infrared
C cloud absorption properties
C
        if(wavmax(ns) .le. 0.7) then
          indxsl = 1
        else if(wavmin(ns) .eq. 0.700) then
          indxsl = 2
        else if(wavmin(ns) .eq. 0.701) then
          indxsl = 3
        else if(wavmin(ns) .eq. 0.702 .or. wavmin(ns) .gt. 2.38) then
          indxsl = 4
        end if
C
C Set cloud extinction optical depth, single scatter albedo,
C asymmetry parameter, and forward scattered fraction:
C
        abarli = abarl(indxsl)
        bbarli = bbarl(indxsl)
        cbarli = cbarl(indxsl)
        dbarli = dbarl(indxsl)
        ebarli = ebarl(indxsl)
        fbarli = fbarl(indxsl)
c
        abarii = abari(indxsl)
        bbarii = bbari(indxsl)
        cbarii = cbari(indxsl)
        dbarii = dbari(indxsl)
        ebarii = ebari(indxsl)
        fbarii = fbari(indxsl)
        do k=1,plev
          do n=1,nloop
            do i=is(n),ie(n)
c
c liquid
c
              tmp1l = abarli + bbarli/rel(i,k)
              tmp2l = 1. - cbarli - dbarli*rel(i,k)
              tmp3l = fbarli*rel(i,k)
c
c ice
c
              tmp1i = abarii + bbarii/rei(i,k)
              tmp2i = 1. - cbarii - dbarii*rei(i,k)
              tmp3i = fbarii*rei(i,k)
C
C Cloud fraction incorporated into cloud extinction optical depth
C
              tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
     $                     *cld(i,k)*sqrt(cld(i,k))
              tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
     $                     *cld(i,k)*sqrt(cld(i,k))
C
C Do not let single scatter albedo be 1; delta-eddington solution
C for non-conservative case:
C
              wcl(i,k) = min(tmp2l,.999999)
              gcl(i,k) = ebarli + tmp3l
              fcl(i,k) = gcl(i,k)*gcl(i,k)
C
              wci(i,k) = min(tmp2i,.999999)
              gci(i,k) = ebarii + tmp3i
              fci(i,k) = gci(i,k)*gci(i,k)
C
C Set aerosol properties
C Conversion factor to adjust aerosol extinction (m2/g)
C
              tauxar(i,k) = 1.e4 * ksa(ns) * uaer(i,k)
C
              wa(i,k)     = wsa(ns)
              ga(i,k)     = gsa(ns)
              fa(i,k)     = gsa(ns)*gsa(ns)
C
              waer(i)     = wa(i,k)
              gaer(i)     = ga(i,k)
              faer(i)     = fa(i,k)
            end do
          end do
        end do
C
C Set reflectivities for surface based on mid-point wavelength
C
        wavmid = 0.5*(wavmin(ns) + wavmax(ns))
C
C Wavelength less  than 0.7 micro-meter
C
        if (wavmid .lt. 0.7 ) then
          do n=1,nloop
            do i=is(n),ie(n)
              albdir(i) = asdir(i)
              albdif(i) = asdif(i)
            end do
          end do
C
C Wavelength greater than 0.7 micro-meter
C
        else
          do n=1,nloop
            do i=is(n),ie(n)
              albdir(i) = aldir(i)
              albdif(i) = aldif(i)
            end do
          end do
        end if
        trayoslp = raytau(ns)/sslp
C
C Layer input properties now completely specified; compute the
C delta-Eddington solution reflectivities and transmissivities
C for each layer, starting from the top and working downwards:
C
        call radded(coszrs   ,trayoslp,pflx   ,abh2o(ns),abo3(ns),
     $              abco2(ns),abo2(ns),uh2o   ,uo3      ,uco2    ,
     $              uo2      ,tauxcl  ,wcl    ,gcl      ,fcl     ,
     $              tauxci   ,wci     ,gci    ,fci      ,tauxar  ,
     $              wa       ,ga      ,fa     ,nloop    ,is      ,
     $              ie       ,rdir    ,rdif   ,tdir     ,tdif    ,
     $              explay   ,exptdn  ,rdndif ,tottrn   )
C
C Compute reflectivity to direct and diffuse radiation for layers below
C by adding succesive layers starting from the surface and working
C upwards:
C
        do n=1,nloop
          do i=is(n),ie(n)
            rupdir(i,plevp) = albdir(i)
            rupdif(i,plevp) = albdif(i)
          end do
        end do
        do k=plev,0,-1
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) +
     $                      rupdif(i,k+1)*tdif(i,k)**2*rdenom
            end do
          end do
        end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
        do k=0,plevp
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
              fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $               exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
C
C Compute flux divergence in each layer using the interface up and down
C fluxes:
C
        do k=0,plev
          do n=1,nloop
            do i=is(n),ie(n)
              flxdiv(i,k) = (fluxdn(i,k  ) - fluxdn(i,k+1)) +
     $                      (fluxup(i,k+1) - fluxup(i,k  ))
            end do
          end do
        end do
C
C Monochromatic computation completed; accumulate in totals; adjust
C fraction within spectral interval to allow for the possibility of
C sub-divisions within a particular interval:
C
        psf = 1.0
        if(ph2o(ns).ne.0.) psf = psf*ph2o(ns)
        if(pco2(ns).ne.0.) psf = psf*pco2(ns)
        if(po2 (ns).ne.0.) psf = psf*po2 (ns)
        do n=1,nloop
          do i=is(n),ie(n)
            solflx(i)  = solin(i)*frcsol(ns)*psf
            fsnt(i) = fsnt(i) + solflx(i)*(fluxdn(i,1) - fluxup(i,1))
            fsntoa(i) = fsntoa(i) + solflx(i)*(fluxdn(i,0) 
     $                  - fluxup(i,0))
            fsns(i) = fsns(i) + solflx(i)*
     $               (fluxdn(i,plevp) - fluxup(i,plevp))
            sfltot(i)  = sfltot(i) + solflx(i)
            fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(i,0)
            fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(i,0)
C
C Down spectral fluxes need to be in mks; thus the .001 conversion factors
C
            if (wavmid .lt. 0.7) then
              sols(i) = sols(i) + exptdn(i,plevp)*solflx(i)*0.001
              solsd(i) = solsd(i) + (fluxdn(i,plevp) -
     $                  exptdn(i,plevp)) * solflx(i)*0.001 
            else
              soll(i) = soll(i) + exptdn(i,plevp)*solflx(i)*0.001
              solld(i) = solld(i) + (fluxdn(i,plevp) -
     $                  exptdn(i,plevp)) * solflx(i)*0.001 
              fsnirtsq(i) = fsnirtsq(i) + 
     $                      solflx(i)*(fluxdn(i,0) - fluxup(i,0))
            end if
              fsnirt(i) = fsnirt(i) + 
     $                    wgtint * solflx(i)*
     $                    (fluxdn(i,0) - fluxup(i,0))

C
          end do
        end do
        do k=0,plev
          do n=1,nloop
            do i=is(n),ie(n)
              totfld(i,k)  = totfld(i,k)  + solflx(i)*flxdiv(i,k)
              fswup(i,k+1) = fswup(i,k+1) + solflx(i)*fluxup(i,k+1)
              fswdn(i,k+1) = fswdn(i,k+1) + solflx(i)*fluxdn(i,k+1)
            end do
          end do
        end do
C
C
C Following code is the diagnostic clear sky computation:
C
C Compute delta-Eddington solution reflectivities and transmissivities
C for the entire column; note, for convenience, we use the same
C reflectivity and transmissivity arrays as for the full calculation
C above, where 0 for layer quantities refers to the entire atmospheric
C column, and where 0 for interface quantities refers to top of atmos-
C phere, while 1 refers to the surface:
C
C
C Compute total column aerosol optical depth:
C
        do n=1,nloop
          do i=is(n),ie(n)
C
C Conversion factor to adjust aerosol extinction (m2/g)
C
            tauaer(i) = 1.e4 * ksa(ns) * utaer(i)
          end do
        end do
        call radclr(coszrs   ,trayoslp,pflx    ,abh2o(ns),abo3(ns) ,
     $              abco2(ns),abo2(ns),uth2o   ,uto3     ,utco2    ,
     $              uto2     ,tauaer  ,waer    ,gaer     ,faer     ,
     $              nloop    ,is      ,ie      ,rdir     ,rdif     ,
     $              tdir     ,tdif    ,explay  ,exptdn   ,rdndif   ,
     $              tottrn   )
C
C Compute reflectivity to direct and diffuse radiation for entire
C column; 0,1 on layer quantities refers to two effective layers
C overlying surface; 0 on interface quantities refers to top of column;
C 2 on interface quantities refers to the surface:
C
        do n=1,nloop
          do i=is(n),ie(n)
            rupdir(i,2) = albdir(i)
            rupdif(i,2) = albdif(i)
          end do
        end do
C
        do k=1,0,-1
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) +
     $                        rupdif(i,k+1)*tdif(i,k)**2*rdenom
            end do
          end do
        end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
        do k=0,2
          do n=1,nloop
            do i=is(n),ie(n)
              rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
              fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
C
        do n=1,nloop
          do i=is(n),ie(n)
            fsntc(i) = fsntc(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
            fsnsc(i) = fsnsc(i) + solflx(i)*(fluxdn(i,2)-fluxup(i,2))
            fsdsc(i) = fsdsc(i) + solflx(i)*fluxdn(i,2)
            fsnrtc(i) = fsnrtc(i) + 
     $                  wgtint * solflx(i) *
     $                  (fluxdn(i,0) - fluxup(i,0))
          end do
        end do
C
C End of clear sky calculation
C
  100 continue                  ! End of spectral interval loop
C
C Compute solar heating rate (k/s)
C
      gocp = gravit/cpair
      do k=1,plev
        do n=1,nloop
          do i=is(n),ie(n)
            qrs(i,k) = -gocp*totfld(i,k)/(pint(i,k) - pint(i,k+1))
          end do
        end do
      end do
c
c Set the downwelling flux at the surface 
c
      do i=1,plon
        fsds(i) = fswdn(i,plevp)
      end do
C
      return
      end
