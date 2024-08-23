#include <misc.h>
#include <params.h>
      SUBROUTINE HIRSRTM(PP, TT, RMIX, O3MIX, TS, ORO,
     $                   TB_IR,BRITEMP)
C -- ------------------------------------------------------------------2
c
c                  ***       VERSION 2.0        ***
c
c This subroutine calculates brightness temperatures at the top of the
c atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
c channels (1,2,3,4).
c
c On Input:
c    pp      -  level pressure (hPa) at which fluxes are evaluated
c               from top of the atmosphere to surface
c    tt      -  layer temperatures (K)
c    rmix    -  layer H2O mass mixing ratio in kg/kg
c    o3mix   -  layer ozone mass mixing ratio kg/kg 
c    ts      -  surface temperature (K)
c    oro     - land-sea flag (sea=0, land=1)
c
c On Ouput:
c    tb_ir   -  infrared brightness temperatures
c    britemp -  microwave brightness temperatures
c
c
c  A flag to include the 4 MSU channels can be switched on (msu_flag=1)
c  and off (msu_flag=0). To decrease the amount of computer time, the
c  microwave routine is changed to a lookup table with almost the same 
c  accuracy as the original routine.
c
c **  last revised 3/31/97 by Richard Engelen **
c
c   This version differs from original version:
c
c     1.  New NOAA10 coefficients
c     2.  Continuum added
c     3.  Any level exceeding 100% RH is changed to 100% RH
c     4.  New channels (2,4,6,8,10,11,12)
c
      IMPLICIT NONE
C -- ------------------------------------------------------------------2

C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <pbtemp.h>
C------------------------------Commons----------------------------------
#include <combtemp.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      REAL  PP(plond,plevp)   ! Pressure for each column and layer (hPa)
      REAL  TT(plond,plev)    ! Temperature of each column and layer (K)
      REAL  TS(plond)         ! Surface temperature (K)
      REAL  RMIX(plond,plev)  ! Water vapor mass mixing ratio (kg/kg)
      REAL  O3MIX(plond,plev) ! Ozone mass mixing ratio (kg/kg)
      REAL  ORO(plond)        ! Land surface flag, sea=0, land=1
C
C Output arguments
C
      REAL  BRITEMP(plond,pnf_msu)   ! HIRS brightness temperatures
      REAL  TB_IR(plond,pnb_hirs)    ! MSU brightness temperatures
C------------------------------Local variables--------------------------
      REAL  UH2O(plev)
      REAL  UO3(plev)
      REAL  A_HIRS(pnb_hirs)
      REAL  B_HIRS(pnb_hirs)
      REAL  XH2O(pnb_hirs)
      REAL  YH2O(pnb_hirs)
      REAL  XO3(pnb_hirs)
      REAL  YO3(pnb_hirs)
      REAL  XCO2(pnb_hirs)
      REAL  YCO2(pnb_hirs)
      REAL  B_IR(pnb_hirs)
      REAL  OTRANS(pnb_hirs)
      REAL  TBAND(pnb_hirs)
      REAL  SCOEF10(pnb_hirs)
      REAL  FCOEF10(pnb_hirs)
      REAL  CWN(pnb_hirs)
      REAL  DTRANS(pnb_hirs)
      REAL  RAD_LAY(pnb_hirs)
      REAL  RADIR(pnb_hirs)
      REAL  RAD(pnb_hirs)
      REAL  RAD2(pnb_hirs)
      REAL  RAD3(pnb_hirs)
      REAL  REFL(pnb_hirs)
      REAL  OTR_MW(pnf_msu)
      REAL  TAU(pnf_msu)
      REAL  TRANS(pnf_msu)
      REAL  OTR_MW2(pnf_msu)
      REAL  TAU2(pnf_msu)
      REAL  TRANS2(pnf_msu)
C     REAL  BRITEMP1(plond,pnf_msu)
C     REAL  BRITEMP2(plond,pnf_msu)
C     REAL  BRITEMP3(plond,pnf_msu)
      REAL  FREQ(pnf_msu)
      REAL  UPATH_H2O
      REAL  UPATH_CO2
      REAL  UPATH_O3
      REAL  UCONT1
      REAL  UCONT2
      REAL  PPATH_H2O
      REAL  PPATH_CO2
      REAL  PPATH_O3
      REAL  SFCTEMP                     ! Surface temperature
      REAL  DP                          ! Pressure depth of layer
      REAL  TLAY                        ! Layer temperature
      REAL  PLAY                        ! Pressure
      REAL  RLAY                        ! Water vapor mixing ratio
      REAL  O3LAY                       ! Ozone mixing ratio.
      REAL  RHOAIR                      ! Layer density
      REAL  E                           ! Saturation pressure
      REAL  DELZ                        ! Height of layer
      REAL  DP2                         ! Pressure depth of layer below
      REAL  TLAY2
      REAL  PLAY2
      REAL  RLAY2
      REAL  RHOAIR2
      REAL  E2
      REAL  DELZ2
      REAL  B_MW
      REAL  B_MW2
      REAL  DTR_MW
      REAL  UCO2
      REAL  PW
      REAL  PSC_H2O                   ! Partial pressure of H2O
      REAL  PSC_CO2                   ! Partial pressure of CO2
      REAL  PSC_O3                    ! Partial pressure of O3
      REAL  T_CONT
      REAL  T_H2O
      REAL  T_CO2
      REAL  T_O3
      REAL  ABS(pnf_msu)              ! Absorption
      REAL  ABS2(pnf_msu)
      REAL  DTR_MW2
C
      REAL T_MALK, PLNCK, BTEMP
      EXTERNAL T_MALK, PLNCK, BTEMP
      INTEGER I, IB, IF, ICOL    ! Loop control

C------------------------------Data Statments---------------------------
      DATA SCOEF10/3.09,2.64,2.28,1.004,0.429,1.945,11.95/
      DATA FCOEF10/.00434,.00301,.0018,4.33e-5,0.000393,0.0738,1.110/
      DATA A_HIRS/.0183,-.00203,.0653,0.21797,0.29846,0.04612,0.06453/
      DATA B_HIRS/.99992,.99994,.9998,.99957,.9996,.99963,1.0006/
      DATA CWN/680.23,704.33,733.13,899.5,1224.07,1363.32,1489.42/
      DATA FREQ/50.31,53.73,54.96,57.95/


      DATA XH2O/0.41,0.52,0.17,0.05,0.68,5.02,17.24/
      DATA YH2O/0.70,1.41,0.14,0.02,1.09,125.9,1194.5/

      DATA XCO2/64.32,13.29,5.04,0.05,0.03,0.25,0.01/
      DATA YCO2/1325.2,121.01,14.36,0.01,0.03,0.02,0.01/

      DATA XO3/32.83,28.9,27.44,0.01,0.68,0.37,0.01/
      DATA YO3/77.44,66.7,67.44,1.6,4.17,0.07,0.01/

      REAL EPS
      DATA EPS/0.5/

      DO 60 icol=1,plon

        SFCTEMP=TS(icol)

        DO IB=1,pnb_hirs
          RADIR(IB)=0.
        ENDDO

        UPATH_H2O=0.0
        UPATH_CO2=0.0
        UPATH_O3=0.0
        UCONT1=0.
        UCONT2=0.
        PPATH_H2O=0.0
        PPATH_CO2=0.0
        PPATH_O3=0.
C       IF(MSU_FLAG.EQ.1) THEN
          DO 11 IF=1,pnf_msu
             TAU(IF)=0.
             TRANS(IF)=0.
             RAD(IF)=0.
             OTR_MW(IF)=1.
             TAU2(IF)=0.
             TRANS2(IF)=0.
             RAD2(IF)=0.
             OTR_MW2(IF)=1.
11        CONTINUE
C       ENDIF
        DO 12 IB=1,pnb_hirs
           TBAND(IB)=1.0
           OTRANS(IB)=1.
12      CONTINUE


        DO 100 I=1,plev
      
          DP=PP(icol,I+1)-PP(icol,I)
          TLAY=TT(icol,I)
          PLAY=SQRT(PP(icol,I)*PP(icol,I+1))
          RLAY=RMIX(icol,I)
          O3LAY=O3MIX(icol,I)
          
          RHOAIR=PLAY*100./(R*TT(icol,I))
          E = PLAY*RLAY/(RLAY+0.6220)
          DELZ=DP*0.1/(RHOAIR*GRAV)
          
C         IF(MSU_FLAG.EQ.1) THEN
            DP2=PP(icol,plev+2-I)-PP(icol,plev+1-I)
            TLAY2=TT(icol,plev+1-I)
            PLAY2=SQRT(PP(icol,plev+1-I)*PP(icol,plev+2-I))
            RLAY2=RMIX(icol,plev+1-I)

            RHOAIR2=PLAY2*100./(R*TT(icol,plev+1-I))
            E2 = PLAY2*RLAY2/(RLAY2+0.6220)
            DELZ2=DP2*0.1/(RHOAIR2*GRAV)
          
C
C         MICROWAVE TRANSFER
C
            CALL LOOKUP ( PLAY-E, TLAY, ABS)
            CALL LOOKUP (PLAY2-E2, TLAY2, ABS2)
            DO 10 IF=1,pnf_msu
              TAU(IF)=TAU(IF)+ABS(IF)*DELZ
              TRANS(IF) = EXP(-TAU(IF))
              B_MW = 1.47445E-23*FREQ(IF)**3/(EXP(0.047981*FREQ(IF)
     $               /TLAY)-1)
              DTR_MW = OTR_MW(IF)-TRANS(IF)
              RAD(IF)=RAD(IF) + B_MW*DTR_MW
              OTR_MW(IF)=TRANS(IF)
           
              TAU2(IF)=TAU2(IF)+ABS2(IF)*DELZ2
              TRANS2(IF) = EXP(-TAU2(IF))
              B_MW2 = 1.47445E-23*FREQ(IF)**3/(EXP(0.047981*FREQ(IF)
     $             /TLAY2)-1)
              DTR_MW2 = OTR_MW2(IF)-TRANS2(IF)
              RAD2(IF)=RAD2(IF) + B_MW2*DTR_MW2
              OTR_MW2(IF)=TRANS2(IF)
10          CONTINUE
C         ENDIF

C
C                IR TRANSFER
C
          UH2O(I)=RLAY*100.*DP/GRAV
          UO3(I)=O3LAY*100.*DP/GRAV
          UCO2=RCO2*100.*DP/GRAV
          PW=PLAY*RLAY*28.97/18.
          UCONT1=UCONT1+UH2O(I)*((PLAY-PW)/1013.)*(296./TLAY)
          UCONT2=UCONT2+UH2O(I)*(PW/1013.)*(296./TLAY)

          UPATH_H2O=UPATH_H2O+UH2O(I)
          PPATH_H2O=PPATH_H2O+UH2O(I)*PLAY
          PSC_H2O=PPATH_H2O/UPATH_H2O

          UPATH_CO2=UPATH_CO2+UCO2
          PPATH_CO2=PPATH_CO2+UCO2*PLAY
          PSC_CO2=PPATH_CO2/UPATH_CO2

          UPATH_O3=UPATH_O3+UO3(I)
          PPATH_O3=PPATH_O3+UO3(I)*PLAY
          PSC_O3=PPATH_O3/UPATH_O3

          DO 15 IB=1,pnb_hirs
            T_CONT=EXP(-SCOEF10(IB)*UCONT2-FCOEF10(IB)*UCONT1)

            T_H2O=t_malk(UPATH_H2O,PSC_H2O,XH2O(IB),YH2O(IB))
            T_CO2=t_malk(UPATH_CO2,PSC_CO2,XCO2(IB),YCO2(IB))
            T_O3=t_malk(UPATH_O3,PSC_O3,XO3(IB),YO3(IB))
            TBAND(IB)=T_CO2*T_O3*T_H2O*T_CONT


            B_IR(IB)=PLNCK(CWN(IB),A_HIRS(IB),B_HIRS(IB),TLAY)
            DTRANS(IB)=OTRANS(IB)-TBAND(IB)
            RAD_LAY(IB)=  B_IR(IB)*DTRANS(IB)
            RADIR(IB)= RADIR(IB)+ RAD_LAY(IB)
            OTRANS(IB) = TBAND(IB)

  15      CONTINUE
 100    CONTINUE   ! End of loop over vertical levels


C
C    Add in the surface contribution to the radiance, including
C    reflection
C
        
C       IF(MSU_FLAG.EQ.1) THEN
          IF(ORO(ICOL).eq.0)THEN
            EPS=0.65      ! Ocean
          ELSE
            EPS=0.9    ! Land or sea-ice
          END IF
      
          DO 16 IF=1,pnf_msu
            B_MW = 1.47445E-23*FREQ(IF)**3/(EXP(0.047981*FREQ(IF)/
     $         SFCTEMP)-1)
            RAD(IF) = RAD(IF) + EPS*TRANS(IF)*B_MW
            REFL(IF)=(1-EPS)*RAD2(IF)*TRANS(IF)
C           BRITEMP1(ICOL,IF) = 0.047981*FREQ(IF)/LOG(1.0 
C    $                + 1.47445E-23*FREQ(IF)**3/RAD(IF))
C           BRITEMP2(ICOL,IF) = 0.047981*FREQ(IF)/LOG(1.0 
C    $                + 1.47445E-23*FREQ(IF)**3/RAD2(IF))
C           BRITEMP3(ICOL,IF) = 0.047981*FREQ(IF)/LOG(1.0 
C    $                + 1.47445E-23*FREQ(IF)**3/REFL(IF))
     
            RAD3(IF) = RAD(IF) + REFL(IF)
            BRITEMP(ICOL,IF) = 0.047981*FREQ(IF)/LOG(1.0 
     $                + 1.47445E-23*FREQ(IF)**3/RAD3(IF))
16        CONTINUE
C       ENDIF
      
        DO 17 IB=1,pnb_hirs
          B_IR(IB)=PLNCK(CWN(IB),A_HIRS(IB),B_HIRS(IB),SFCTEMP)
          RADIR(IB)=RADIR(IB)+B_IR(IB)*TBAND(IB)
          TB_IR(icol,IB)=BTEMP(CWN(IB),A_HIRS(IB),B_HIRS(IB),RADIR(IB))
17      CONTINUE


60    CONTINUE   ! End of loop over columns


      END
C
C==================================================================
C
      REAL FUNCTION BTEMP(WVN,A,B,CHNRAD)
C
C*  Calculates the brightness temperature given the channel radiance.
c   Uses Planck function tuned to TOVS frequencies
c
      IMPLICIT NONE
C------------------------------Arguments--------------------------------
      REAL WVN
      REAL CHNRAD
      REAL A
      REAL B
C------------------------------Local variables--------------------------
      REAL C1
      REAL C2
C
C   Planck function parameters
C   NOTE: Parameters A and B are temperature correction factors
c   which are dependent on the channel and satellite.  These
c   parameters were extracted from the RTTOVS model.
c
      PARAMETER( C1=1.191066E-08 )
      PARAMETER( C2=1.438833 )

      BTEMP=(C2*WVN/ALOG(C1*WVN**3/CHNRAD+1.)-A)/B

      RETURN
      END
C
C==================================================================
C
      REAL FUNCTION PLNCK(WVN,A,B,T)

C
C Planck function
C
      IMPLICIT NONE
      REAL WVN,A,B,T,C1,C2

      PARAMETER( C1=1.191066E-08 )
      PARAMETER( C2=1.438833 )

      PLNCK=C1*(WVN)**3/(EXP(C2*WVN/(A+B*T))-1.)

      RETURN
      END
C
C===================================================================
C
      REAL FUNCTION T_MALK(U,P,X,Y)

      IMPLICIT NONE
      REAL U, P, X, Y
      REAL P0, DNU, PI, B, BP

      PARAMETER( P0=1013. )
      PARAMETER( DNU=10. )

      PI=ACOS(-1.)
      B=4.*X**2/(PI*Y*DNU)
      BP=PI*B*P/P0
      T_MALK = EXP(-0.5*BP*(SQRT(1.+4.*Y*U/(DNU*BP))-1.))

      RETURN
      END
C
C===================================================================
C
      SUBROUTINE LOOKUP(P,T,ABS)

      IMPLICIT NONE 
#include <pbtemp.h>
      REAL P, T, ABS(pnf_msu)
C
      INTEGER IF, I, J, N, M
      parameter( N=17 )
      parameter( M=16 )
      REAL XX(N), YY(M)
      REAL ZZ1(N,M), ZZ2(N,M), ZZ3(N,M), ZZ4(N,M)
      REAL T1, T2
      DATA YY/5.0,7.5,10.,25.,50.,100.,200.,300.,
     $        400.,500.,600.,700.,800.,900.,1000,1050./
      data XX/160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,
     $        260.,270.,280.,290.,300.,310.,320./
      data ZZ1/0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     $         0.0000,0.0000,0.0002,0.0002,0.0001,0.0001,0.0001,
     $         0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,
     $         0.0000,0.0000,0.0000,0.0000,0.0000,0.0008,0.0007,
     $         0.0006,0.0005,0.0004,0.0004,0.0003,0.0003,0.0003,
     $         0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001,
     $         0.0001,0.0032,0.0027,0.0023,0.0020,0.0017,0.0015,
     $         0.0013,0.0011,0.0010,0.0009,0.0008,0.0008,0.0007,
     $         0.0006,0.0006,0.0006,0.0005,0.0130,0.0108,0.0091,
     $         0.0078,0.0068,0.0059,0.0052,0.0046,0.0041,0.0037,
     $         0.0033,0.0030,0.0028,0.0025,0.0023,0.0022,0.0020,
     $         0.0292,0.0244,0.0206,0.0176,0.0152,0.0133,0.0117,
     $         0.0103,0.0092,0.0083,0.0074,0.0067,0.0061,0.0056,
     $         0.0052,0.0048,0.0045,0.0521,0.0434,0.0367,0.0314,
     $         0.0271,0.0237,0.0208,0.0184,0.0164,0.0147,0.0132,
     $         0.0120,0.0109,0.0100,0.0092,0.0085,0.0079,0.0818,
     $         0.0681,0.0576,0.0493,0.0426,0.0371,0.0326,0.0289,
     $         0.0257,0.0230,0.0207,0.0188,0.0170,0.0156,0.0143,
     $         0.0132,0.0122,0.1182,0.0985,0.0833,0.0712,0.0616,
     $         0.0537,0.0472,0.0418,0.0372,0.0333,0.0300,0.0271,
     $         0.0246,0.0225,0.0206,0.0189,0.0175,0.1617,0.1348,
     $         0.1139,0.0975,0.0842,0.0734,0.0645,0.0571,0.0508,
     $         0.0455,0.0409,0.0370,0.0336,0.0306,0.0281,0.0258,
     $         0.0238,0.2123,0.1770,0.1496,0.1280,0.1106,0.0965,
     $         0.0848,0.0750,0.0667,0.0597,0.0537,0.0486,0.0441,
     $         0.0402,0.0368,0.0338,0.0312,0.2701,0.2253,0.1905,
     $         0.1630,0.1408,0.1228,0.1079,0.0954,0.0849,0.0760,
     $         0.0684,0.0618,0.0560,0.0511,0.0467,0.0429,0.0396,
     $         0.3353,0.2798,0.2366,0.2024,0.1749,0.1525,0.1340,
     $         0.1185,0.1055,0.0944,0.0849,0.0767,0.0695,0.0634,
     $         0.0579,0.0532,0.0490,0.3707,0.3094,0.2617,0.2239,
     $         0.1935,0.1687,0.1482,0.1311,0.1166,0.1043,0.0938,
     $         0.0848,0.0769,0.0700,0.0640,0.0588,0.0541/
      data ZZ2/0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,
     $         0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,
     $         0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,
     $         0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,
     $         0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0002,
     $         0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,
     $         0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,
     $         0.0002,0.0002,0.0010,0.0010,0.0010,0.0011,0.0011,
     $         0.0012,0.0012,0.0013,0.0013,0.0014,0.0014,0.0014,
     $         0.0015,0.0015,0.0015,0.0015,0.0015,0.0040,0.0038,
     $         0.0038,0.0039,0.0041,0.0042,0.0044,0.0046,0.0048,
     $         0.0050,0.0051,0.0053,0.0054,0.0055,0.0056,0.0056,
     $         0.0057,0.0144,0.0135,0.0131,0.0130,0.0132,0.0135,
     $         0.0140,0.0145,0.0150,0.0155,0.0161,0.0166,0.0170,
     $         0.0174,0.0177,0.0180,0.0182,0.0525,0.0473,0.0439,
     $         0.0417,0.0405,0.0399,0.0399,0.0402,0.0407,0.0414,
     $         0.0421,0.0429,0.0437,0.0445,0.0452,0.0458,0.0464,
     $         0.1135,0.1005,0.0913,0.0848,0.0802,0.0772,0.0752,
     $         0.0740,0.0734,0.0732,0.0733,0.0736,0.0739,0.0744,
     $         0.0749,0.0753,0.0757,0.1971,0.1730,0.1553,0.1422,
     $         0.1326,0.1255,0.1202,0.1164,0.1137,0.1118,0.1104,
     $         0.1095,0.1088,0.1084,0.1081,0.1079,0.1077,0.3029,
     $         0.2645,0.2358,0.2140,0.1975,0.1848,0.1750,0.1675,
     $         0.1617,0.1572,0.1537,0.1509,0.1487,0.1468,0.1453,
     $         0.1440,0.1429,0.4307,0.3749,0.3325,0.3000,0.2747,
     $         0.2550,0.2395,0.2272,0.2174,0.2095,0.2030,0.1978,
     $         0.1934,0.1897,0.1865,0.1837,0.1812,0.5801,0.5037,
     $         0.4452,0.3998,0.3642,0.3360,0.3134,0.2953,0.2805,
     $         0.2684,0.2584,0.2500,0.2429,0.2368,0.2316,0.2269,
     $         0.2228,0.7503,0.6505,0.5734,0.5131,0.4655,0.4273,
     $         0.3966,0.3715,0.3509,0.3339,0.3196,0.3075,0.2972,
     $         0.2882,0.2805,0.2736,0.2674,0.9408,0.8148,0.7168,
     $         0.6396,0.5782,0.5288,0.4887,0.4557,0.4284,0.4056,
     $         0.3864,0.3700,0.3560,0.3437,0.3331,0.3236,0.3151,
     $         1.1504,0.9956,0.8746,0.7787,0.7021,0.6400,0.5893,
     $         0.5475,0.5126,0.4834,0.4586,0.4374,0.4191,0.4032,
     $         0.3892,0.3768,0.3657,1.2620,1.0919,0.9586,0.8528,
     $         0.7679,0.6991,0.6427,0.5961,0.5572,0.5245,0.4967,
     $         0.4728,0.4523,0.4343,0.4186,0.4046,0.3921/
      data ZZ3/0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001,
     $         0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,
     $         0.0001,0.0001,0.0001,0.0004,0.0004,0.0003,0.0003,
     $         0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,
     $         0.0003,0.0003,0.0003,0.0003,0.0002,0.0002,0.0006,
     $         0.0006,0.0006,0.0006,0.0006,0.0006,0.0006,0.0006,
     $         0.0006,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,
     $         0.0004,0.0004,0.0040,0.0039,0.0038,0.0038,0.0037,
     $         0.0037,0.0036,0.0035,0.0034,0.0033,0.0032,0.0031,
     $         0.0030,0.0029,0.0028,0.0027,0.0026,0.0154,0.0150,
     $         0.0147,0.0145,0.0143,0.0141,0.0139,0.0136,0.0133,
     $         0.0130,0.0126,0.0122,0.0119,0.0115,0.0111,0.0107,
     $         0.0103,0.0543,0.0528,0.0518,0.0511,0.0504,0.0498,
     $         0.0491,0.0484,0.0475,0.0465,0.0454,0.0443,0.0431,
     $         0.0419,0.0406,0.0393,0.0380,0.1702,0.1617,0.1560,
     $         0.1521,0.1492,0.1470,0.1449,0.1430,0.1409,0.1387,
     $         0.1363,0.1338,0.1310,0.1281,0.1251,0.1219,0.1187,
     $         0.3229,0.3005,0.2846,0.2729,0.2642,0.2573,0.2516,
     $         0.2467,0.2421,0.2377,0.2334,0.2290,0.2246,0.2200,
     $         0.2153,0.2104,0.2055,0.5089,0.4673,0.4363,0.4128,
     $         0.3946,0.3801,0.3681,0.3579,0.3489,0.3407,0.3330,
     $         0.3257,0.3185,0.3115,0.3045,0.2975,0.2906,0.7263,
     $         0.6605,0.6104,0.5716,0.5410,0.5162,0.4957,0.4783,
     $         0.4631,0.4496,0.4373,0.4258,0.4149,0.4046,0.3945,
     $         0.3848,0.3752,0.9737,0.8791,0.8060,0.7486,0.7028,
     $         0.6654,0.6344,0.6080,0.5852,0.5651,0.5469,0.5303,
     $         0.5149,0.5005,0.4867,0.4735,0.4609,1.2496,1.1217,
     $         1.0219,0.9429,0.8792,0.8271,0.7836,0.7467,0.7149,
     $         0.6870,0.6621,0.6395,0.6188,0.5995,0.5815,0.5645,
     $         0.5482,1.5516,1.3866,1.2568,1.1532,1.0693,1.0004,
     $         0.9428,0.8939,0.8518,0.8150,0.7824,0.7531,0.7263,
     $         0.7018,0.6789,0.6575,0.6373,1.8769,1.6714,1.5088,
     $         1.3782,1.2720,1.1844,1.1111,1.0489,0.9954,0.9488,
     $         0.9076,0.8707,0.8374,0.8069,0.7788,0.7527,0.7282,
     $         2.2222,1.9737,1.7758,1.6162,1.4859,1.3780,1.2876,
     $         1.2109,1.1450,1.0876,1.0371,0.9921,0.9516,0.9147,
     $         0.8809,0.8497,0.8206,2.4013,2.1306,1.9144,1.7396,
     $         1.5966,1.4781,1.3787,1.2944,1.2219,1.1588,1.1034,
     $         1.0541,1.0098,0.9696,0.9328,0.8989,0.8673/
      data ZZ4/0.0027,0.0023,0.0019,0.0016,0.0014,0.0012,0.0011,
     $         0.0009,0.0008,0.0007,0.0006,0.0006,0.0005,0.0005,
     $         0.0004,0.0004,0.0003,0.0059,0.0050,0.0043,0.0037,
     $         0.0032,0.0027,0.0024,0.0021,0.0018,0.0016,0.0014,
     $         0.0013,0.0011,0.0010,0.0009,0.0008,0.0007,0.0105,
     $         0.0089,0.0076,0.0065,0.0056,0.0049,0.0042,0.0037,
     $         0.0033,0.0029,0.0025,0.0023,0.0020,0.0018,0.0016,
     $         0.0015,0.0013,0.0649,0.0550,0.0468,0.0402,0.0347,
     $         0.0301,0.0262,0.0230,0.0202,0.0178,0.0158,0.0141,
     $         0.0126,0.0113,0.0101,0.0091,0.0082,0.2465,0.2097,
     $         0.1794,0.1544,0.1336,0.1162,0.1015,0.0891,0.0785,
     $         0.0695,0.0617,0.0550,0.0491,0.0441,0.0396,0.0358,
     $         0.0323,0.8274,0.7127,0.6168,0.5364,0.4684,0.4108,
     $         0.3616,0.3196,0.2834,0.2522,0.2251,0.2015,0.1810,
     $         0.1630,0.1471,0.1332,0.1208,2.1322,1.8745,1.6542,
     $         1.4649,1.3016,1.1601,1.0371,0.9298,0.8358,0.7532,
     $         0.6805,0.6163,0.5593,0.5088,0.4638,0.4236,0.3876,
     $         3.2645,2.8943,2.5762,2.3012,2.0625,1.8542,1.6717,
     $         1.5114,1.3698,1.2446,1.1334,1.0344,0.9460,0.8668,
     $         0.7958,0.7319,0.6743,4.2737,3.8017,3.3956,3.0442,
     $         2.7385,2.4714,2.2370,2.0305,1.8478,1.6858,1.5415,
     $         1.4126,1.2972,1.1936,1.1004,1.0162,0.9400,5.2092,
     $         4.6423,4.1542,3.7314,3.3633,3.0414,2.7585,2.5090,
     $         2.2882,2.0920,1.9171,1.7607,1.6205,1.4945,1.3809,
     $         1.2782,1.1851,6.0871,5.4326,4.8683,4.3790,3.9525,
     $         3.5791,3.2507,2.9608,2.7039,2.4754,2.2716,2.0892,
     $         1.9256,1.7782,1.6453,1.5251,1.4161,6.9131,6.1782,
     $         5.5437,4.9928,4.5120,4.0906,3.7196,3.3917,3.1008,
     $         2.8419,2.6106,2.4035,2.2175,2.0500,1.8987,1.7617,
     $         1.6374,7.6903,6.8822,6.1833,5.5755,5.0445,4.5784,
     $         4.1676,3.8041,3.4813,3.1937,2.9366,2.7061,2.4989,
     $         2.3121,2.1432,1.9903,1.8514,8.4213,7.5468,6.7890,
     $         6.1290,5.5515,5.0440,4.5961,4.1994,3.8467,3.5321,
     $         3.2507,2.9981,2.7708,2.5657,2.3801,2.2119,2.0590,
     $         9.1086,8.1740,7.3627,6.6548,6.0345,5.4886,5.0062,
     $         4.5785,4.1978,3.8579,3.5535,3.2801,3.0338,2.8114,
     $         2.6100,2.4272,2.2610,9.4365,8.4743,7.6380,6.9077,
     $         6.2673,5.7033,5.2047,4.7622,4.3683,4.0163,3.7009,
     $         3.4175,3.1621,2.9314,2.7224,2.5326,2.3600/
                    
      IF(P.le.5) THEN
        DO IF = 1, pnf_msu
          ABS(IF)=0.0
        END DO
        RETURN
      ENDIF
      
      CALL LOCATE(XX,N,T,I)
      CALL LOCATE(YY,M,P,J)
            
      T1=(T-XX(I))/(XX(I+1)-XX(I))
      T2=(P-YY(J))/(YY(J+1)-YY(J))
      ABS(1)=(1-T1)*(1-T2)*ZZ1(I,J)+T1*(1-T2)*ZZ1(I+1,J)+
     $                T1*T2*ZZ1(I+1,J+1)+(1-T1)*T2*ZZ1(I,J+1)
      ABS(2)=(1-T1)*(1-T2)*ZZ2(I,J)+T1*(1-T2)*ZZ2(I+1,J)+
     $                T1*T2*ZZ2(I+1,J+1)+(1-T1)*T2*ZZ2(I,J+1)
      ABS(3)=(1-T1)*(1-t2)*ZZ3(I,J)+T1*(1-T2)*ZZ3(I+1,J)+
     $                T1*T2*ZZ3(I+1,J+1)+(1-T1)*T2*ZZ3(I,J+1)
      ABS(4)=(1-T1)*(1-T2)*ZZ4(I,J)+T1*(1-T2)*ZZ4(I+1,J)+
     $                T1*T2*ZZ4(I+1,J+1)+(1-T1)*T2*ZZ4(I,J+1)
           
      RETURN
      END
C
C===================================================================
C
      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT NONE
      INTEGER N, J
      REAL XX(N), X
C
      INTEGER JL, JU, JM
C
      JL=0
      JU=N+1
      DO WHILE (JU-JL.GT.1)
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      END DO
      J=JL
      RETURN
      END


