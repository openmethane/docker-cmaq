
C***********************************************************************
C   Portions of Models-3/CMAQ software were developed or based on      *
C   information from various groups: Federal Government employees,     *
C   contractors working on a United States Government contract, and    *
C   non-Federal sources (including research institutions).  These      *
C   research institutions have given the Government permission to      *
C   use, prepare derivative works, and distribute copies of their      *
C   work in Models-3/CMAQ to the public and to permit others to do     *
C   so.  EPA therefore grants similar permissions for use of the       *
C   Models-3/CMAQ software, but users are requested to provide copies  *
C   of derivative works to the Government without restrictions as to   *
C   use by others.  Users are responsible for acquiring their own      *
C   copies of commercial software associated with Models-3/CMAQ and    *
C   for complying with vendor requirements.  Software copyrights by    *
C   the MCNC Environmental Modeling Center are used with their         *
C   permissions subject to the above restrictions.                     *
C***********************************************************************

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/cloud/cloud_acm_ae5_txhgsim/hg_aqchem_data.F,v 1.1.1.1 2010/06/14 16:03:03 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      MODULE HG_AQCHEM_DATA
C-----------------------------------------------------------------------
C Function: Sets up Data for mercury cloud chemistry and contain subroutine to
C           calculate rates and minimum time steps of its reactions.

C Revision History:
C     27 Aug 2008 B.Hutzell: initial implementation
C     06 Jul 2009 J.Bash:    Replaced the Hg(II) reduction by HO2 with the 
C                            reduction mechanism by dicaboxylic acids of 
C                            Si and Ariya 2008 ES&T 
C-----------------------------------------------------------------------


      IMPLICIT NONE

      INTEGER, PARAMETER  :: NPR    = 10  ! number of Hg reactions exclude adsorb/absorb

      INTEGER, PARAMETER  :: NUMOX_v46  = 5  ! number of oxidizing reactions

      INTEGER, PARAMETER  :: NHgRXN = 24 ! number of Hg reactions

      INTEGER, PARAMETER  :: NRXN   = NUMOX_v46 + NHgRXN     ! number of all reactions


C INDICES FOR MERCURIC REDOX REACTIONS

      INTEGER          :: IHg        ! counter index over Hg rxns
      INTEGER, SAVE    :: IHg_ox = 1 ! Hg(aq) + O3(aq) => HgO(aq)    '
      INTEGER, SAVE    :: IHgSO3 = 2 ! HgSO3 => Hg(aq) + S(IV) van Loon et al.      '
      INTEGER, SAVE    :: IHGHY  = 3 ! Hg(OH)2 => Hg + products      '
      INTEGER, SAVE    :: IOHRAD = 4 ! Hg(aq) + OHRAD(aq) => Hg(II)  '
      INTEGER, SAVE    :: ICLI   = 5 ! oxidation by CL2 dissocations products (HOCL and OCL)
      INTEGER, SAVE    :: IORGC  = 6 ! Reduction of all Hg(II) species by oxalate RXN 6
C                                      Hg(II) + R(CO2)2 =>  Hg+ + products
C                                      Hg+    + HORCO2  =>  Hg(aq) + products
      INTEGER, SAVE    :: IHgdisulf = 7  ! Reduction of Hg(SO3)2 by oxalate RXN 7
      INTEGER, SAVE    :: IHgOHp    = 8  ! Reduction of HgOHp    by oxalate RXN 8
      INTEGER, SAVE    :: IHgOHCL   = 9  ! Reduction of HgOHCL   by oxalate RXN 9
      INTEGER, SAVE    :: IHgCl2    = 10  ! Reduction of HgCL2   by oxalate RXN 9
      INTEGER, SAVE    :: ISHgCl2    = NPR + 1  ! index for HgCl2    sorption
      INTEGER, SAVE    :: ISHgSO3    = NPR + 3  ! index for HgSO3    sorption
      INTEGER, SAVE    :: ISHgHY     = NPR + 5  ! index for HgHY     sorption
      INTEGER, SAVE    :: ISHgdisulf = NPR + 7  ! index for Hg(SO3)2-- sorption
      INTEGER, SAVE    :: ISHgOHp    = NPR + 9  ! index for HgOH     sorption
      INTEGER, SAVE    :: ISHgOHCl   = NPR + 11 ! index for HgOHCl   sorption
      INTEGER, SAVE    :: ISHgII     = NPR + 13 ! index for Hg(II)   sorption
      INTEGER, SAVE    :: IDHgCl2    = NPR + 2  ! index for HgCl2    desorption
      INTEGER, SAVE    :: IDHgSO3    = NPR + 4  ! index for HgSO3    desorption
      INTEGER, SAVE    :: IDHgHY     = NPR + 6  ! index for HgHY     desorption
      INTEGER, SAVE    :: IDHgdisulf = NPR + 8  ! index for Hg(SO3)2-- desorption
      INTEGER, SAVE    :: IDHgOHp    = NPR + 10 ! index for HgOH     desorption
      INTEGER, SAVE    :: IDHgOHCl   = NPR + 12 ! index for HgOHCl   desorption
      INTEGER, SAVE    :: IDHgII     = NPR + 14 ! index for Hg(II)   desorption

      REAL         TMPVAR          ! temporary variable for calculations

      REAL         COSINE_ZENITH   ! solar zenith anagle

      REAL         ClI             ! Cl(I) conc in cloudwater (mol/liter), from Cl2(aq)

      REAL         DTWHG_MIN       ! safe timestep for mercury chemistry
      REAL         DTWHG( NRXN )   ! timesteps for mercury chemistry

      REAL         HgII            ! Hg(II) conc in cloudwater (mol/liter)

      REAL         HgH             ! Henry's Law Constant for Hg
      REAL         HgL             ! Hg conc in cloudwater (mol/liter)

      REAL         HgOHp           ! HgOH+ conc in cloudwater (mol/liter)
      REAL         HgOHp1          ! First dissociation constant for HgOH+
      REAL         HgOHp1I         ! Inverse HgOHp1

      REAL         HgOHCl          ! HgOHCl conc in cloudwater (mol/liter)
      REAL         HgOHCl1         ! First dissociation constant for HgOHCl
      REAL         HgOHCl1I        ! Inverse HgOHCl1

      REAL         HgCl21          ! First dissociation constant for HgCl2
      REAL         HgCl21I         ! Inverse HgCl21
      REAL         HgCl2H          ! Henry's Law Constant for HgCl2
      REAL         HeffHgCl2       ! Effective Henry's Law Constant for HgCl2
      REAL         HgCl21H         ! HgCl21 * HgCl2H
      REAL         HgCl2L          ! HgCl2(aq) conc in cloudwater (mol/liter)

      REAL         HgHY            ! Hg(OH)2 conc in cloudwater (mol/liter)
      REAL         HgHY1           ! Equilibrium constant for Hg(OH)2
      REAL         HgHY1I          ! Inverse HgHY1

      REAL         HgSO3           ! HgSO3 conc in cloudwater (mol/liter)
      REAL         HgSO31          ! Equilibrium constant for HgSO3
      REAL         HgSO31I         ! Inverse HgSO31

      REAL         Hgdisulf        ! Hg(SO3)2-- conc in cloudwater (mol/liter)
      REAL         Hgdisulf1       ! Equilibrium constant for Hg(SO3)2--
      REAL         Hgdisulf1I      ! Inverse Hgdisulf1

      REAL         K1K2I           ! HgSO31I * Hgdisulf1I
      REAL         K4K5I           ! HgOHp1I * HgHY1I
      REAL         K4K6I           ! HgOHp1I * HgOHCl1I
      REAL         K6a             ! Rate constant used in ox of Hg by chlorine
      REAL         K6b             ! Rate constant used in ox of Hg by chlorine
      REAL         KORGC            ! Rate constant used in redux of all HgII by R(CO2)2

      REAL         PHGAKNA         ! aitken mercury aerosol in water (moles/L)
      REAL         PHGACCA         ! accum  mercury aerosol in water (moles/L)

      REAL         PHgCl20         ! total HgCl2 partial pressure (atm)
      REAL         PHgCl2F         ! gas only HgCl2 partial pressure (atm)

      REAL         PHG_AEROSOL_BAK ! previous sorbed Hg(II) in liquid phase
      REAL         PHG_AEROSOL_NOW ! current  sorbed Hg(II) in liquid phase
      REAL         PHG_AEROSOL_DEL ! PHG_AEROSOL_NOW - PHG_AEROSOL_BAK

      REAL         RHg6            !     ''        Hg with chlorine
      REAL         RHg_ox          !     ''        Hg ox by O3
      REAL         RHgSO3          !     ''
      REAL         RHgSO3I         !     ''        , inverse
      REAL         RHgHY           !     ''
      REAL         RHgHYI          !     ''        , inverse
      REAL         ROHRAD          !     ''

      REAL         SOVERL          ! Ratio of sorbed to dissolved Hg(II)
      REAL         SHgCl2          ! Sorbed HgCl2 conc in cloudwater (mol/liter)
      REAL         SHgSO3          ! Sorbed HgSO3 conc in cloudwater (mol/liter)
      REAL         SHgdisulf       ! Sorbed Hg(SO3)2-- conc in cloudwater (mol/liter)
      REAL         SHgOHp          ! Sorbed HgOHp conc in cloudwater (mol/liter)
      REAL         SHgHY           ! Sorbed Hg(OH)2 conc in cloudwater (mol/liter)
      REAL         SHgOHCl         ! Sorbed HgOHCl conc in cloudwater (mol/liter)
      REAL         SHgII           ! Sorbed HgII conc in cloudwater (mol/liter)
      REAL         SORBED_HG_INIT         ! Initial total sorbed Hg in cloudwater (mol/liter)

      REAL         AHgCl2    ! additional mass from RGMS to add to SHgCl2
      REAL         AHgSO3    ! additional mass from RGMS to add to SHgSO3
      REAL         AHgdisulf ! additional mass from RGMS to add to SHgdisulf
      REAL         AHgOHp    ! additional mass from RGMS to add to SHgOHp
      REAL         AHgHY     ! additional mass from RGMS to add to SHgHY
      REAL         AHgOHCl   ! additional mass from RGMS to add to SHgOHCl
      REAL         AHgII     ! additional mass from RGMS to add to SHgII


      REAL         XLHg            !
      REAL         XLHgCl2         !

      REAL         DHgDT( NHgRXN ) ! rate of Hg spcs prod in cld (mol/liter/sec)
      REAL         DHg  ( NHgRXN ) ! Hg species product produced over tstep DTWHG(0)
c      CHARACTER*30 RXNSTRING ( NHgRXN ) ! Hg reaction description
C...........INSORB is a logical variable used to initialize the fraction
C           of aqueous Hg(II) sorbed to suspended carbon

      LOGICAL      INSORB

c...this section of code based on phot.F method of retrieving lat and
c   lon from the grid coordinate file

C      REAL, ALLOCATABLE, SAVE :: LAT( :,: ) ! north lat in deg (cross pt.)
C      REAL, ALLOCATABLE, SAVE :: LON( :,: ) ! west long in deg (cross pt.)
C      INTEGER ALLOCSTAT
      CHARACTER( 16 ) :: VARNM

      INTEGER      GXOFF, GYOFF            ! global origin offset from file
C for INTERPX
      INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2

c these vars based on jproc method of zenith calculation

      real    csz  ! FRACTION OF THE SOLAR CONSTANT

C     XLON   - LONGITUDE (west longitudes > 0)
C     XLAT   - LATITUDE  (south latitudes 0 north latitudes)
      CONTAINS

      SUBROUTINE INIT_AQCHEM_HG(TEMP, WCAVG, DTRMV)


      IMPLICIT NONE

C      INCLUDE 'CONST.EXT'          ! constants
      INCLUDE SUBST_CONST        ! commonly used constants

      REAL TEMP     ! AIR TEMP (K)
      REAL WCAVG    ! liquid water content (kg/m3)
      REAL DTRMV    ! Default Timestep (s)


      REAL XL       ! conversion factor (liter-atm/mol)


C...........EXTERNAL FUNCTIONS and their descriptions:

      REAL          HLCONST
      EXTERNAL      HLCONST


      real    xlat    ! latitude  of single cell
      real    xlon    ! longitude of single cell
      INTEGER      DAY                ! julian day of year
      INTEGER      YEAR               ! year (yyyy)
      INTEGER      MONTH              ! month counter
      INTEGER ihr,imin,isec ! time in hrs, min, seconds
      INTEGER idate         ! yyyymmdd
      real    gmt  ! utc, military time in decimal hours
      real    gw(2),aw(12),fac,wk,ga(2),aa(12)

      LOGICAL, SAVE   :: SCNDTIME   = .FALSE. ! Flag for 2ND call to aqchem for testing
      REAL, PARAMETER :: H2ODENS    = 1000.0  ! density of water (kg/m3) at 20 C and 1 ATM



        scndtime = .TRUE.

      HgH    = HLCONST( 'HG              ', TEMP, .FALSE., 0.0 )
      HgCl2H = HLCONST( 'HGIIGAS         ', TEMP, .FALSE., 0.0 )

      XL = MOLVOL*(WCAVG/H2ODENS)*(TEMP/STDTEMP)   ! conversion factor (l-atm/mol)

      XLHg    = HgH   * XL
      XLHgCl2 = HgCl2H* XL

C...dissociation constant for dissolved mercury species

c     K1:
      HgSO31 = 2.e-13      ! (M)
      HgSO31I= 1./HgSO31   !  (M**-1)
c     K2:
      Hgdisulf1 = 4.e-12   ! (M)
      Hgdisulf1I= 1./Hgdisulf1
c     K3:
      HgCl21  = 1.E-14     ! (M**2)
      HgCl21I = 1./HgCl21  ! (M**-2)
      HgCl21H = HgCl21 * HgCl2H   !RB: HgCl21H is never used
c     K4:
      HgOHp1   =  2.51e-11 ! (M)
      HgOHp1I   =  1./HgOHp1
c     K5:
      HgHY1  = 6.31e-12    ! (M)
      HgHY1I = 1./HgHY1
c     K6:
      HgOHCl1 = 3.72e-8    ! (M)
      HgOHCl1I = 1./3.72e-8

      K1K2I = HgSO31I * Hgdisulf1I
      K4K5I = HgOHp1I * HgHY1I
      K4K6I = HgOHp1I * HgOHCl1I


C...Hg reaction rates  RXN
C  Hg(aq) + O3(aq) => HgO(aq)     RXN 1  4.7E7
      RHg_ox  = 4.7E+7

C  HgSO3 => Hg(aq) + S(IV)        RXN 2 from Van Loon et al.
      RHgSO3  =  temp * EXP( (31.971*temp - 12595.)/temp )
      RHgSO3I = 1./RHgSO3
      IHgSO3  = 2

C  Hg(OH)2 => Hg + products        RXN 3   3.0E-7
      if ( COSINE_ZENITH .gt. 0. ) then
        RHgHY  = 6.0E-7 * COSINE_ZENITH   ! csz = 1 with sun at zenith
        RHgHYI = 1./RHgHY
      else
        RHgHY  = 0.
        RHgHYI = 0.               ! not used if rhghy = 0
      endif

C  Hg(aq) + OHRAD(aq) => Hg(II)     RXN 4  2.0E9
      ROHRAD = 2.0E9

C  OXIDATION OF Hg(aq) BY CHLORINE (HOCl and OCl-)  RXN 5

C  Hg(aq) + HOCl(aq) => Hg(II) + products
C  Hg(aq) + OCl-   --(H+)-->  Hg(II) + products
C  HOCl <=> H+ + OCl-                   K = 10**-7.5

      K6a     = 2.09e6
      K6b     = 1.99e6

C Reduction of all Hg(II) species by R(CO2)2      RXN 6
C  Hg(II) + R(CO2)2 =>  Hg+ + products
C  Hg+    + HORCO2 =>  Hg(aq) + products
C  Overall: Hg(II) + R(CO2)2 => Hg(aq) + products
      if ( COSINE_ZENITH .gt. 0.0 ) then
         KORGC   = 1.2e4 * COSINE_ZENITH        ! 1/(M s)
      else 
         KORGC   = 0.0
      endif

C ADSORPTION AND DESORPTION:  RXN 7-20
C
C All Hg(II) species sorb/desorb at the same time rate
C
C HgCl2L => SHgCl2                 RXN 7  adsorb = RHgad
C SHgCl2 => HgCl2L                 RXN 8  desorb = RHgde
C HgSO3  => SHgSO3                 RXN 9  adsorb = RHgad
C SHgSO3 => HgSO3                  RXN 10 desorb = RHgde
C HgHY   => SHgHY                  RXN 11 adsorb = RHgad
C SHgHY  => HgHY                   RXN 12 desorb = RHgde
C Hgdisulf  => SHgdisulf           RXN 13 adsorb = RHgad
C SHgdisulf => Hgdisulf            RXN 14 desorb = RHgde
C HgOHp     => SHgOHp              RXN 15 adsorb = RHgad
C SHgOHp    => HgOHp               RXN 16 desorb = RHgad
C HgOHCl    => SHgOHCl             RXN 17 adsorb = RHgad
C SHgOHCl   => HgOHCl              RXN 18 desorb = RHgad
C HgII      => SHgII               RXN 19 adsorb = RHgad
C SHgII     => HgII                RXN 20 desorb = RHgad

       DTWHG(:)  = DTRMV
       DTWHG_MIN = DTRMV

       DO IHg = 1, NHgRXN

         DHgDT( IHg ) = 0.0
         DHg  ( IHg ) = 0.0

       END DO

       SORBED_HG_INIT  = 0.0
c...set history rgms to 0 hg before the time loop begins.
       PHG_AEROSOL_BAK  = 0.0
 

       RETURN
       END SUBROUTINE 
C
       SUBROUTINE MERCURY_RATES(WCAVG, DTRMV, EC, O3L, HPLUS, OHRAD, ORGC,
     &                          HOCL, OCL)

       IMPLICIT NONE

C       INCLUDE 'CONST.EXT'          ! constants
      INCLUDE SUBST_CONST        ! commonly used constants


       REAL WCAVG  ! Liquid water content   (kg/m3)
       REAL EC     ! elemental carbon acc+akn aerosol in cloudwater (mol/liter)
C       REAL DTRMV  ! Default Timestep (s)
       REAL DTRMV   ! Minimum time step required
       REAL O3L    ! ozone dissolved in cloud water  (mol/liter)
       REAL HPLUS  ! hydrogen ion concentration (mol/liter)
       REAL OHRAD  ! OH ion concentration (mol/liter)
       REAL ORGC   ! Assumed Oxalic acid concentration (mol/liter)
       REAL HOCL   ! HOCL in cloud water  (mol/liter)
       REAL OCL    ! OCL in cloud water  (mol/liter)


       REAL SovD   ! sorbed Hg(II) / desorbed Hg(II)
       REAL SF     ! fraction of aqueous Hg(II) sorbed to E.C.
       REAL ECw    ! elem. carbon suspended in water (g / liter)
       REAL EC5    ! EC * 5

       REAL RHgad    ! HG overall adsorption
       REAL RHgadI   !     ''        , inverse
       REAL RHgde    ! HG overall desorption
       REAL RHgdeI   !     ''        , inverse
       REAL KeqHOCl  ! equilib constant used in ox of Hg by chlorine
       REAL mincarb  ! 1.0E-18/WCAVG in g/Liter
       
       REAL DHGPROLOSS ! sum of mercury species production and loss terms over dt
       REAL DHg_SCALE  ! Rate scale to preserve mass under S controled time step

       REAL, PARAMETER :: Kp    = 900.0  ! sorption coeff. from Seigneur et al.(1998)
       REAL, PARAMETER :: TFOLD = 3600.0 ! e-folding time for adsorption equilibrium (sec)

C        DTWHG(:) = DTRMV
C
C...first Hg reaction:  Hg(aq) + O3(aq) => HgO(aq)

        DHgDT(IHg_ox) = RHg_ox * O3L * HgL
c        IF ((DHgDT(IHg_ox) .EQ. 0.0).OR.( HgL .LE. 1.0E-30 ) ) THEN
c          DTWHG(IHg_ox + NUMOX_v46) = DTRMV
c        ELSE
cC         at completion, dtw*dhgdt=min(O3L, HgL) (limiting reagent)
c          DTWHG(IHg_ox + NUMOX_v46) = 0.1 * (MIN(O3L, HgL))/DHgDT(IHg_ox )
c        END IF

C...second Hg reaction:  HgSO3 => Hg(aq) + S(IV)
!WTH:(08/14/09): includes ORGC reduction reaction
        DHgDT(IHgSO3) = (RHgSO3 + KORGC * ORGC) * HgSO3
        IF ( ( DHgDT(IHgSO3) .EQ. 0.0 ).OR.(HgSO3 .LE. 1.0E-30) ) THEN
          DTWHG(IHgSO3 + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IHgSO3 + NUMOX_v46 ) = .1* RHgSO3I
       END IF


C...third Hg reaction:  Hg(OH)2 => Hg + products
!WTH:(08/14/09): includes ORGC reduction reaction
        DHgDT(IHgHY) = (RHgHY + KORGC * ORGC) * HgHY
        IF ( ( DHgDT(IHgHY) .EQ. 0.0 ) .OR. ( HgHY .LE. 1.0E-30 ) ) THEN
          DTWHG(IHgHY + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IHgHY + NUMOX_v46 ) = .1* RHgHYI
        END IF

C...fourth Hg reaction:  Hg(aq) + OHRAD(aq) => Hg(II)
        DHgDT(IOHRAD) = ROHRAD * OHRAD * HgL

c        IF ( ( DHgDT(IOHRAD) .EQ. 0.0 ) .OR. ( HgL .LE. 1.0E-30 ) .OR.
c     &      ( OHRAD .LE. 1.0E-30 ) ) THEN
c          DTWHG(IOHRAD + NUMOX_v46 ) = DTRMV
c        ELSE
cC         at completion, dtw*dhgdt=min(OHRAD, HgL) (limiting reagent)
c          DTWHG(IOHRAD + NUMOX_v46 ) = 0.1 * (MIN(OHRAD, HgL))/DHgDT(IOHRAD)
c        END IF

C...fifth Hg reaction:  OXIDATION OF Hg(aq) BY CHLORINE (HOCl and OCl-)
C
C  Hg(aq) + HOCl(aq) => Hg(II) + products
C  Hg(aq) + OCl-   --(H+)-->  Hg(II) + products
C  HOCl <=> H+ + OCl-                   K = 10**-7.5
C
        KeqHOCl = 3.1623e-8    ! equals 10**(-7.5)

        RHg6 =  K6a*HOCL + K6b*OCL
        DHgDT(IClI) = RHg6 *  HgL 
c        IF ( ( DHgDT(IClI) .EQ. 0.0 ) .OR. ( HgL .LE. 1.0E-30 ) ) THEN
c          DTWHG(IClI + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IClI+ NUMOX_v46 ) = 0.1 * (MIN(HgL, ClI))/DHgDT(IClI)
c        END IF
c        if(DTWHG(IClI+ NUMOX_v46 ) .lt. 0.1 ) then  !temporary check
c           DTWHG(IClI+ NUMOX_v46 ) = DTRMV
c           print*,"Warning: DTW too short for Chlorine ox. of Hg(0)"
c        endif


C...sixth Hg reaction:  REDUCTION OF Hg(II) by R(CO2)2
C
C  Hg(II) + R(CO2)2 =>  Hg+ + products
C  Hg+    + HORCO2 =>  Hg(aq) + products
C  Overall: Hg(II) + R(CO2)2 => Hg(aq) + products
C
C        TMPVAR =HgII+HgHY+HgSO3+Hgdisulf+HgOHp+HgOHCl+HgCl2L

        DHgDT(IORGC) = KORGC * ORGC *  HgII 
        IF ( ( DHgDT(IORGC) .EQ. 0.0 ) .OR.( HgII .LE. 1.0E-30 ) ) THEN
          DTWHG(IORGC + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IORGC+ NUMOX_v46 ) = 0.1 * (MIN(ORGC, HgII ))/DHgDT(IORGC)
        END IF
	
        DHgDT(IHgdisulf) = KORGC * ORGC * Hgdisulf 
        IF ( ( DHgDT(IHgdisulf) .EQ. 0.0 ) .OR.( Hgdisulf .LE. 1.0E-30 ) ) THEN
          DTWHG(IHgdisulf + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IHgdisulf+ NUMOX_v46 ) = 0.1 * (MIN(ORGC, Hgdisulf ))/DHgDT(IHgdisulf)
        END IF

        DHgDT(IHgOHp) = KORGC * ORGC * HgOHp
        IF ( ( DHgDT(IHgOHp) .EQ. 0.0 ) .OR.( HgOHp .LE. 1.0E-30 ) ) THEN
          DTWHG(IHgOHp + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IHgOHp+ NUMOX_v46 ) = 0.1 * (MIN(ORGC, HgOHp ))/DHgDT(IHgOHp)
        END IF

        DHgDT(IHgOHCL) = KORGC * ORGC * HgOHCL
        IF ( ( DHgDT(IHgOHCL) .EQ. 0.0 ) .OR.( HgOHCL .LE. 1.0E-30 ) ) THEN
          DTWHG(IHgOHCL + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IHgOHCL+ NUMOX_v46 ) = 0.1 * (MIN(ORGC, HgOHCL ))/DHgDT(IHgOHCL)
        END IF

	DHgDT(IHgCl2) = KORGC * ORGC * HgCl2L
        IF ( ( DHgDT(IHgCl2) .EQ. 0.0 ) .OR.( HgCl2L .LE. 1.0E-30 ) ) THEN
          DTWHG(IHgCl2 + NUMOX_v46 ) = DTRMV
        ELSE
          DTWHG(IHgCl2+ NUMOX_v46 ) = 0.1 * (MIN(ORGC, HgCl2L ))/DHgDT(IHgCl2)
        END IF

C...new sorption/desorption code to allow adjustable elemental carbon
C   aerosol air concentrations (Russ Bullock, 09/19/2000)
c   using carbon aerosol within water

c EC here in moles/liter of water
c
c times 5 to estimate amount of elem carbon avail if 5% of pm2.5 rather
c than 1% had been set to elem carbon originally in aero_driver.F of

c      EC5 = max(0.0, (EC * 5.0)) ! moles/L
      EC5 = max(0.0, EC) ! moles/L

c set min carb as equivilant to 1e-12 microg/m**3, and convert this
c to grams carbon/Liter of water.  1e-12 microg/m**3 is
c 1e-18g/m**3, which is (1e-18g/m**3)/WCAVG grams/Liter,
c where WCAVG is liquid water content (kg/m3)
c (1e-18g carb)/m**3 / [(?Kg H20)/m**3] = (1e-18g carb) /(?Kg H20)
c = (1e-18g carb) /(?Liters H20)
c (density of H20 is assumed to be close to 1kg/Liter)


      ECw     = EC5 * 12.01      ! moles/liter * 12g per mole = g/liter
      mincarb = 1.0E-18/WCAVG    ! in g/L
      ECw     = max(ECw,mincarb) ! forced min of 1e-12ug/m**3 as g/liter

C      Kp = 900.               ! adapted from Seigneur et al. (1998)
C      TFOLD = 3600.           ! e-folding time for adsorption equilibrium

      SovD = Kp*ECw           ! sorbed Hg(II) / desorbed Hg(II)
      SF   = SovD/(1.0+SovD)  ! fraction sorbed at equilibrium

      RHgad  = SF/TFOLD        ! adsorb rate for e-folding time
      RHgadI = 1./RHgad
      RHgde  = RHgad/SovD      ! desorb rate to yield SF at equilibrium
      RHgdeI = 1./RHgde



C Adsorption of Hg Species HgSO3 -> SHgSO3
        DHgDT(ISHgSO3) = RHgad * HgSO3
c        IF ((DHgDT(ISHgSO3) .EQ. 0.0 ).OR.(HgSO3  .LE. 1.0E-30 ))THEN
c          DTWHG(ISHgSO3 + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgSO3 + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgSO3 -> HgSO3
        DHgDT(IDHgSO3) = RHgde * SHgSO3
c        IF ( ( DHgDT(IDHgSO3).EQ.0.0).OR.(SHgSO3.LE.1.0E-30 ) ) THEN
c          DTWHG(IDHgSO3 + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgSO3 + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C Adsorption of Hg Species HgCl2L -> SHgCl2
        DHgDT(ISHgCl2) = RHgad * HgCl2L
c        IF ((DHgDT(ISHgCl2) .EQ. 0.0) .OR. (HgCl2L .LE. 1.0E-30)) THEN
c          DTWHG(ISHgCl2 + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgCl2 + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgCl2 -> HgCl2L
        DHgDT(IDHgCl2) = RHgde * SHgCl2
c        IF ((DHgDT(IDHgCl2) .EQ. 0.0 ).OR.( SHgCl2 .LE. 1.0E-30 ))THEN
c          DTWHG(IDHgCl2 + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgCl2 + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C Adsorption and Loss of Hg Species Hgdisulf -> SHgdisulf
        DHgDT(ISHgdisulf) = RHgad * Hgdisulf
c        IF ((DHgDT(ISHgdisulf).EQ.0.0 ).OR.(Hgdisulf.LE.1.0E-30))THEN
c          DTWHG(ISHgdisulf + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgdisulf + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgdisulf -> Hgdisulf
        DHgDT(IDHgdisulf) = RHgde * SHgdisulf
c        IF ((DHgDT(IDHgdisulf).EQ.0.0).OR.(SHgdisulf.LE.1.0E-30))THEN
c          DTWHG(IDHgdisulf + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgdisulf + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C Adsorption and Loss of Hg Species HgOHp    -> SHgOHp
        DHgDT(ISHgOHp) = RHgad * HgOHp
c        IF ( (DHgDT(ISHgOHp) .EQ. 0.0 ).OR.( HgOHp .LE. 1.0E-30))THEN
c          DTWHG(ISHgOHp + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgOHp + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgOHp    -> HgOHp
        DHgDT(IDHgOHp) = RHgde * SHgOHp
c        IF (( DHgDT(IDHgOHp) .EQ. 0.0 ).OR.( SHgOHp .LE.1.0E-30) )THEN
c          DTWHG(IDHgOHp + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgOHp + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C Adsorption and Loss of Hg Species HgOHCl   -> SHgOHCl
        DHgDT(ISHgOHCl) = RHgad * HgOHCl
c        IF ((DHgDT(ISHgOHCl).EQ.0.0 ).OR.( HgOHCl .LE. 1.0E-30 ) )THEN
c          DTWHG(ISHgOHCl + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgOHCl + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgOHCl   -> HgOHCl
        DHgDT(IDHgOHCl) = RHgde * SHgOHCl
c        IF (( DHgDT(IDHgOHCl).EQ.0.0).OR.(SHgOHCl.LE. 1.0E-30 ) )THEN
c          DTWHG(IDHgOHCl + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgOHCl + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C Adsorption of Hg Species HgHY -> SHgHY
        DHgDT(ISHgHY) = RHgad * HgHY
c        IF ((DHgDT(ISHgHY).EQ.0.0 ) .OR. ( HgHY .LE. 1.0E-30 ) ) THEN
c          DTWHG(ISHgHY + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgHY + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgHY -> HgHY
        DHgDT(IDHgHY) = RHgde * SHgHY
c        IF ((DHgDT(IDHgHY) .EQ. 0.0 ).OR.( SHgHY .LE. 1.0E-30 ) ) THEN
c          DTWHG(IDHgHY + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgHY + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C Adsorption of Hg Species HgII -> SHgII
        DHgDT(ISHgII) = RHgad * HgII
c        IF (( DHgDT(ISHgII).EQ.0.0 ) .OR. ( HgII .LE. 1.0E-30 ) ) THEN
c          DTWHG(ISHgII + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(ISHgII + NUMOX_v46 ) = .1* RHgadI
c        END IF

C Desorption of Hg Species SHgII -> HgII
        DHgDT(IDHgII) = RHgde * SHgII
c        IF ((DHgDT(IDHgII) .EQ. 0.0 ).OR.( SHgII .LE. 1.0E-30 ) ) THEN
c          DTWHG(IDHgII + NUMOX_v46 ) = DTRMV
c        ELSE
c          DTWHG(IDHgII + NUMOX_v46 ) = .1* RHgdeI
c        END IF

C      DTWHG_MIN = MINVAL(DTWHG)
C Don't allow any aqeuous mercury species concentrations to become negative
      ! If the losses of HgII are greater than the produciton terms scale
      ! the losses to zero out the HgII concentration in the time step
      
      DTWHG_MIN = MINVAL(DTWHG)
      
      END SUBROUTINE MERCURY_RATES

      END MODULE HG_AQCHEM_DATA
