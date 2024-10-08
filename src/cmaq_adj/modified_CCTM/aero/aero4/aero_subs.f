
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

C RCS file, release, date & time of last delta, author, state, [locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero4/aero_subs.f,v 1.1.1.1 2010/06/14 16:02:58 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

c  This collection of subroutines is version 4 of the CMAQ aerosol 
c  component (AERO4).  It differs from the RPM version in that the 
c  geometric standard deviations of the Aitken and accumulation modes
c  are now variable.  This is made possible by adding the second moment
c  (modal surface area divided by PI) as a predicted quantity.
c
c                                 Coded by Dr. Francis S. Binkowski

c //////////////////////////////////////////////////////////////////
c  SUBROUTINE AEROPROC advances the number, second moment, and mass
c   concentrations for each mode over the time interval DT.  Other 
c   than CBLK, ORGPROD, and VAPORS, all arguments are scalar.
c
c  KEY SUBROUTINES CALLED:
c     GETPAR, ORGAER3, HETCHEM, EQL3, HCOND3, NEWPART3, GETCOAGS,
c     INTERCOAG_GH, INTRACOAG_GH
c
c  KEY FUNCTIONS CALLED:  ERF, ERFC, GETAF
c
c  REVISION HISTORY:
c     Coded in December 1999 by Dr. Francis S. Binkowski
c      Modified from older versions used in CMAQ, eliminated separate
c      routines for coagulation and growth by including these processes
c      in line.  This code is now a scalar code and has no internal
c      spatial arrays.
c
c FSB 05/17/00  new version of RPMARES included
c
c FSB 05/30/00  Fixed minor bug in awater.f
c
c FSB Correction to extinction coefficient in getbext.
c
c FSB 07/21/00 corrected units on ORGRATES and ORGBRATE_IN from 
c      ppm/min to ppm/sec in AEROPROC and AEROSTEP.
c
c FSB 11/30/00 following changes from 07/28/00
c     Fixed problem for new emissions file version
c     Combined emissions for M3
c     Used a fixed value of Dpmin for plotting
c     Added variables OMEGA_AT & OMEGA_AC for partitioning
c     Eliminated the restriction on relative humidity for nucleation.
c     Added a branch in EQL for very low relative humidity (<1%).
c
c FSB Following changes to RPMARES:
c     Number of iterations reduced from 150 to 50.
c     Iterations are used only if 0.5 < = RATIO.
c     In calculating the molality of the bisulfate ion, a
c     MAX(1.0e-10, MSO4 ) is used.
c
c FSB 08/08/01 Changes to NEWPART to correct mass rate and to AEROSTEP
c     to correct sulfate, include emissions, and trap problem in
c     accumulation mode number calculation.
c
c FSB 09/19/01 Version AE3, major changes
c     Organics are done with Dr. Benedikt Schell's approach. (ORGAER3)
c     Particle production uses Kulmala approach (NEWPART3)
c     Condensational factors for sulfate and organics are calculated
c      separately.
c     Emissions are assumed to be input in vertical diffusion.
c     Include files replaced by Fortran 90 Modules.
c     The modules also contain subroutines.
c
c FSB 10/24/01 Added the treatment of gaseous N2O5 -> aerosol HNO3
c
c FSB 10/25/01 Changed the mass transfer calculation for Aitken to 
c     accumulation mode by coagulation as recommended by Dr. Benedikt
c     Schell.
c
c SJR 04/24/02 Replaced thermodynamic code RPMARES with ISORROPIA
c     --Shawn J. Roselle
c
c GLG 04/04/03 Modifications to allow for evaporation of semi-volatile
c     organics from aerosol phase.  --Gerald L. Gipson
c
c FSB 11/18/03 Corrections to sulfate condensation.  Previously, 
c     SCONDRATE was undefined when SO4RATE=0 and SCONDRATE was 
c     negative when SO4RATE < DMDT_so4.
c
c PVB 01/08/04 Several changes in preparation for simulation of 2001
c      calendar year --Dr. Prakash V. Bhave
c    -Added interface to new subroutine, GETCOAGS, for calculating 
c      coagulation rates.  GETCOAGS is used instead of Gauss-Hermite
c      quadrature for computational efficiency, by setting
c      FASTCOAG_FLAG = .TRUE.  --PVB
c    -Removed SOA from the definition of "DRY" aerosol.  Aerosol surface
c      area is now transported without SOA.  See notes in GETPAR, AERO,
c      and AERO_DEPV subroutines.  --SJR
c    -Moved EQL3 call from AEROPROC to AEROSTEP, immediately following
c      the ORGAER3 call.  This is a side effect of transporting aerosol 
c      surface area without SOA.  --SJR
c    -New subroutine, HCOND3, to calculate condensation rates for
c      2nd and 3rd moments.  Results are unchanged.  --FSB
c    -Revised method of calculating SOA.  Partition SOA to the modes 
c      in proportion to the amounts of total organic mass (SOA plus 
c      primary) in each mode.  Modal geometric standard deviations are
c      now preserved during SOA condensation and evaporation. --FSB
c    -Combined the former subroutines AEROPROC and AEROSTEP into this
c      subroutine; retained the name AEROPROC.  --PVB
c
c PVB 09/21/04 added in-line documentation with input from FSB.
c     Changed MWH2SO4 from 98.07354 to 98.0 g/mol for consistency with
c     the mechanism files.
c
c PVB 09/27/04 removed the IF(XM3.GT.0.0) mode-merging precondition 
c     because it caused significant erroneous mode crossover.  Fix
c     suggested by Dr. Chris Nolte.
c
c PVB 01/19/05 Added SO4RATE to the EQL3 call vector.  This is necessary
c     to ensure that the gas and inorganic fine PM (i+j) concentrations
c     are in thermodynamic equilibrium at the end of each time step.
c
c PVB 05/02/05 Modified ERF statement function for negative arguments such 
c     that erf(-x) = -erf(x).  Previous version had erf(-x) = erf(x).
c
c PVB 04/06/06 Added GAMMA_N2O5 to the AEROPROC and EQL3 call vectors, 
c     so it can be written to the aerosol diagnostic file.
c
c PVB 11/02/07 Moved heterogenous N2O5 chemistry from EQL3 to a new
c     subroutine, HETCHEM.
c
c  REFERENCES:
c   1. Binkowski, F.S. and U. Shankar, The regional particulate matter
c      model 1. Model description and preliminary results, J. Geophys.
c      Res., Vol 100, No D12, 26101-26209, 1995.
c
c   2. Binkowski, F.S. Aerosols in Models-3 CMAQ, Chapter 10 of Science
c      Algorithms of the EPA Models-3 Community Multiscale Air Quality
c      (CMAQ) Modeling System, EPA/R-99/030, March 1999.
c      Available at: http://www.epa.gov/asmdnerl/models3
c
c   3. Binkowski, F.S. and S.J. Roselle, Models-3 Community 
c      Multiscale Air Quality (CMAQ) model aerosol component 1:
c      Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
c      doi:10.1029/2001JD001409, 2003.
c             
c   4. Reid, Prausnitz, and Poling The Properties of Gases and
c      Liquids, 4th edition, McGraw-Hill, 1987, pp 587-588
c
c   5. Jiang, W. and H. Roth, A detailed review and analysis of 
c      science, algorithms, and code in the aerosol components of
c      Models-3/CMAQ.  1. Kinetic and thermodynamic processes in the
c      AERO2 module.  Report Number PET-1534-03S, 2003.
c
c   6. Bhave, P.V., S.J. Roselle, F.S. Binkowski, C.G. Nolte, S. Yu,
c      G.L. Gipson, and K.L. Schere, CMAQ aerosol module development:
c      recent enhancements and future plans, Paper No. 6.8, CMAS Annual
c      Conference, Chapel Hill, NC, 2004.

      SUBROUTINE AEROPROC( NSPCSDA,
     &                     CBLK, DT, LAYER,
     &                     AIRTEMP,AIRPRS,AIRDENS,AIRRH,
     &                     SO4RATE,
     &                     ORGPROD,NPSPCS, VAPORS, NCVAP,
     &                     XLM, AMU,
     &                     DGATK, DGACC, DGCOR,
     &                     XXLSGAT, XXLSGAC,
     &                     PMASSAT, PMASSAC,PMASSCO,
     &                     PDENSAT,PDENSAC,PDENSCO,
     &                     GAMMA_N2O5,
     &                     LOGDEV )

      USE AERO_INFO_AE4    ! Module containing coagulation subroutines,
                           ! GETPAR, GETAF, CBLK indices, and several 
                           ! constants

      IMPLICIT NONE

C *** ARGUMENTS

      REAL DT              ! time step [sec]
      INTEGER LAYER        ! model layer index

      INTEGER NSPCSDA      ! number of variables in CBLK
      REAL CBLK( NSPCSDA ) ! main array of variables

      REAL AIRTEMP         ! Air temperature [ K ]
      REAL AIRPRS          ! Air pressure in [ Pa ]
      REAL AIRDENS         ! Air density [ kg / m**3 ]
      REAL AIRRH           ! Fractional relative humidity

      REAL SO4RATE         ! Gas-phase sulfate production rate 
                           !  [ ug / m**3 s ]

      INTEGER NPSPCS       ! Number of SOA-producing ROGs
      REAL ORGPROD(NPSPCS) ! change in ROG mixing ratio (ppm)

      INTEGER NCVAP        ! Number of partitioning SVOCs
      REAL VAPORS( NCVAP ) ! partitioning SVOC species concentrations

      REAL XLM             ! atmospheric mean free path [ m ]
      REAL AMU             ! atmospheric dynamic viscosity [ kg/m s ]

      REAL DGATK           ! Aitken mode geometric mean diameter [m]
      REAL DGACC           ! accumulation mode geometric mean diam [m]
      REAL DGCOR           ! coarse mode geometric mean diameter [m]

      REAL XXLSGAT         ! natural log of geometric standard
      REAL XXLSGAC         !  deviations

      REAL PMASSAT         ! Aitken mode mass concentration [ ug/m**3 ]
      REAL PMASSAC         ! accumulation mode mass conc [ ug/m**3 ]
      REAL PMASSCO         ! coarse mode mass concentration [ ug/m**3 ]

      REAL PDENSAT         ! Aitken mode particle density [ kg / m**3 ]
      REAL PDENSAC         ! accumulation mode density [ kg / m**3 ]
      REAL PDENSCO         ! coarse mode particle density [ kg / m**3 ]

      REAL GAMMA_N2O5      ! N2O5 heterogeneous reaction probability [ ]

      INTEGER LOGDEV       ! unit number of log file


C *** LOCAL VARIABLES

      CHARACTER* 16  PNAME
      PARAMETER( PNAME = 'AEROPROC        ')
      INTEGER    SPC        ! Loop index

      REAL      SQRT_TEMP   ! square root of ambient temperature

C *** variables for coagulation processes (all double precision)

      REAL*8    LAMDA       ! mean free path [ m ]

C *** Intramodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
      REAL*8    BATAT( 2 )  ! Aitken mode
      REAL*8    BACAC( 2 )  ! accumulation mode

C *** Intermodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
      REAL*8    BATAC( 2 )  ! Aitken to accumulation
      REAL*8    BACAT( 2 )  ! accumulation from Aitken

c *** Intermodal coagulation rate [ m**3/s ] ( 3rd moment )
      REAL*8    c3ij        ! Aitken to accumulation

C *** 3rd moment intermodal transfer rate by coagulation
      REAL*8    C30ATAC     ! Aitken to accumulation

c *** Near Continnuum regime (independent of mode)
      REAL*8    KNC         ! KNC = TWO3 * BOLTZ *  AIRTEMP / AMU

c *** Free Molecular regime (depends upon modal density)
      REAL*8    KFMAT       ! KFMAT = SQRT(3.0*BOLTZ*AIRTEMP/PDENSAT)
      REAL*8    KFMAC       ! KFMAC = SQRT(3.0*BOLTZ*AIRTEMP/PDENSAC)
      REAL*8    KFMATAC     ! KFMATAC = SQRT( 6.0 * BOLTZ * AIRTEMP /
                            !                ( PDENSAT + PDENSAC ) )

c *** Parameters of the aerosol size distributions
      REAL*8 DGATK_D        ! double precision diameters
      REAL*8 DGACC_D
      REAL*8 SGATK_D        ! double precision standard deviations
      REAL*8 SGACC_D
      REAL*8 XXLSGAT_D      ! double precision natural logs of 
      REAL*8 XXLSGAC_D      !  standard deviations


C *** variables for advancing concentrations one time step

      REAL*8 A, B, C
      REAL*8 M1, M2, Y0, Y
      REAL*8 DHAT, P, PEXPDT, EXPDT
      REAL*8 LOSS, PROD, POL, LOSSINV
      REAL SCONDRATE ! SO4 condensation rate [ ug/m**3 s ]
      REAL TSO4, OLDSO4, NEWSO4, TMASS
      REAL FACTRANS ! special factor to compute mass transfer
      REAL M20      ! for initial condidtions in time stepping


C *** Variables for new particle formation:

      REAL XH2SO4     ! steady state H2SO4 concentration
      REAL*8 DMDT_so4 ! particle mass production rate [ ug/m**3 s ]
      REAL*8 DNDT     ! particle number production rate [ # / m**3 s ]
      REAL*8 DM2DT    ! second moment production rate [ m**2 / m**3 s]


C *** variables for secondary organic aerosol production

      REAL ORGRATE    ! [ ug / m**3 s ]
      REAL ORGBRATE   ! [ ug / m**3 s ]

      REAL OLDSOA_A   ! anthropogenic SOA at beginning of time step
      REAL OLDSOA_B   ! biogenic SOA at beginning of time step
      REAL OLDPOA     ! primary organics (does not change)
      REAL SUMPROD    ! total conc of SVOC produced during time 
                      !  step [ ug / m**3 ]  (not used)

      REAL SOA_A      ! anthropogenic SOA at end of time step [ug/m^3]
      REAL SOA_B      ! biogenic SOA at end of time step [ug/m^3]

C     In following definitions, *_I = Aitken mode and *_J = accumulation

      REAL OLD_M3_I, OLD_M3_J  ! input values of third moments
                               ! [ mom / m**3 ] 
      REAL OLD_M2_I, OLD_M2_J  ! input values of second moments
                               ! [ mom / m**3 ]
      REAL OMASS_I, OMASS_J    ! input values of total organic mass
                               ! [ ug / m**3 ]
      REAL FRACI, FRACJ        ! fractions of organic mass in each mode
      REAL NEW_M3_I, NEW_M3_J  ! updated third moments
                               ! [ mom / m**3 ]
      REAL NEW_M2_I, NEW_M2_J  ! updated second moments
                               ! [ mom / m**3 ]
                               

C *** variables for calculating condensational growth

c ***  size-dependent term in the condensational-growth expressions
C      defined in Equations A13-A14 of [Binkowski & Shankar,1995]

      REAL*8 FCONCAT_SO4( 2 )  ! Aitken mode 2nd and 3rd moments
      REAL*8 FCONCAC_SO4( 2 )  ! Accumulation mode 2nd and 3rd moments

      REAL*8 FCONCAT_ORG( 2 )  ! Aitken mode 2nd and 3rd moments
      REAL*8 FCONCAC_ORG( 2 )  ! Accumulation mode 2nd and 3rd moments

C *** modal partition factors [ dimensionless ]
C      defined in Equations A17-A18 of [Binkowski & Shankar,1995]

      REAL*8 OMEGA_AT_SO4  ! Aitken mode 2nd and 3rd moments
      REAL*8 OMEGA_AC_SO4  ! Accumulation mode 2nd and 3rd moments

      REAL*8 OMEGA_AT_ORG  ! Aitken mode 2nd and 3rd moments
      REAL*8 OMEGA_AC_ORG  ! Accumulation mode 2nd and 3rd moments

      REAL FCONCM1_SO4   ! reciprocals of total condensational rates
      REAL FCONCM1_ORG   ! 

      REAL CHEMRATE_SO4  ! sulfate production rate for third moment
      REAL CHEMRATE_ORG  ! organic production rate for third moment

c *** growth rates for 2nd and 3rd moments due to condensation
c     of precursor vapor on existing particles [ mom/m**3 s ]
      REAL*8 CGRAT( 2 )    ! Aitken mode
      REAL*8 CGRAC( 2 )    ! Accumulation

      REAL DIFFSULF ! molecular diffusiviity for sulfuric acid
                    ! [ m**2 /sec ]; calculated using Eqn 11-4.4 of 
                    ! Reid et al. (1987) at 273.16 K and 101325 Pa
       PARAMETER( DIFFSULF = 9.36E-06 )
      REAL DIFFORG
       PARAMETER( DIFFORG  = 5.66E-06  )
      REAL DIFFCORR ! Correction to DIFFSULF & DIFFORG for pressure

      REAL MWH2SO4  ! molecular weight of H2SO4 [ kg/mole ]
       PARAMETER( MWH2SO4 = 98.0E-3  )
      REAL MWORG
       PARAMETER( MWORG = 167.0E-3 )

      REAL    COFCBAR_SO4 ! Temperature-independent coefficients
      REAL    COFCBAR_ORG !  for calculating molecular velocities
                          !  = sqrt((8*Rgas)/(pi*MW))
       SAVE    COFCBAR_SO4, COFCBAR_ORG

C *** inputs to subroutine HCOND3
      REAL*8  AM0AT, AM0AC ! zeroeth moments
      REAL*8  AM1AT, AM1AC ! first moments
      REAL*8  AM2AT, AM2AC ! second moments
      REAL    DV_SO4 ! molecular diffusivity of H2SO4 vapor 
                     ! after correction for ambient conditions
      REAL    DV_ORG ! molecular diffusivity of ORG vapors
                     ! after correction for ambient conditions
      REAL    CBAR_SO4 ! molecular velocity of H2SO4                      
      REAL    CBAR_ORG ! molecular velocity of secondary organics
      REAL ALPHSULF ! Accommodation coefficient for sulfuric acid
       PARAMETER ( ALPHSULF = 1.0 )
      REAL ALPHORG  ! Accommodation coefficient for secondary organics
       PARAMETER ( ALPHORG = 1.0 )

C *** Variables for mode merging
      REAL        AAA, XNUM, XM2, XM3,XXM2, XXM3
      REAL        FNUM, FM2, FM3, PHNUM, PHM2, PHM3
      REAL        ERF, ERFC  ! Error and complementary error function
      REAL        XX         ! dummy argument for ERF and ERFC

C *** Define the minimum atmospheric concentrations

      REAL CONMIN         ! minimum concentration [ppm or ug/m3]
       PARAMETER( CONMIN = 1.0E-30 )

      REAL AEROCONCMIN_AC ! minimum sulfate concentration in
                          !  acccumulation mode [ ug/m**3 ]
                          !  changed to 1 pg/m3 by FSB 12/13/99
       PARAMETER ( AEROCONCMIN_AC = 1.0E-6 )

      REAL AEROCONCMIN_AT ! minimum sulfate concentration in
                          !  Aitken mode [ ug/m**3 ]; smaller than
                          !  any reported tropospheric concentration
       PARAMETER ( AEROCONCMIN_AT = 1.0E-6 * AEROCONCMIN_AC )

      REAL AEROCONCMIN_CO ! minimum coarse mode mass concentration
                          !  [ ug/m**3 ]; set such that minimum number
                          !  concentration = 1.0 #/m**3
       PARAMETER ( AEROCONCMIN_CO = 1.889544E-05 )

      REAL*8 NUMMIN_AT ! minimum number conc in Aitken mode [#/m**3]
      REAL*8 NUMMIN_AC ! minimum number conc in accumulation mode
      REAL*8 NUMMIN_C  ! minimum number conc in coarse mode [#/m**3]

      REAL*8 M2MIN_AT  ! minimum 2nd moment conc in Aitken mode [m2/m3]
      REAL*8 M2MIN_AC  ! minimum 2nd moment conc in accumulation mode

C *** Miscellaneous constants

      REAL*8    ONE
       PARAMETER ( ONE = 1.0D0 )
      REAL*8    ONE3D
       PARAMETER( ONE3D = 1.0D0 / 3.0D0 )
      REAL*8    TWO3D
       PARAMETER( TWO3D = 2.0D0 * ONE3D )

C *** Logical flags (default values set within this subroutine)

      LOGICAL M3_WET_FLAG   ! flag to include water in the 3rd
                            !  moment calculation

      LOGICAL FASTCOAG_FLAG ! flag to obtain coagulation coefficients
                            !  by analytical approximation (TRUE) or 
                            !  by Gauss-Hermite quadrature (FALSE)

      LOGICAL FIRSTIME
       DATA   FIRSTIME / .TRUE./
       SAVE FIRSTIME,
     &                NUMMIN_AT, NUMMIN_AC,NUMMIN_C,
     &                M2MIN_AT, M2MIN_AC

C     :::::::::::::::::::::::::::::::::::::
c *** Statement function given for error function. Source is
c     Meng, Z., and J.H.Seinfeld (1994) On the source of the
c      submicrometer droplet mode of urban and regional aerosols.
c      Aerosol Sci. and Technology, 20:253-265.
c      They cite Reasearch & Education Asociation (REA),
c      (1991) Handbook of Mathematical, Scientific, and Engineering
c      Formulas, Tables, Functions, Graphs, Transforms: REA,
c      Piscataway, NJ. p. 493.
c     Inserted the SIGN(1.0,XX) multiplier so that erf(-x) = -erf(x)

       ERF(XX) = SIGN(1.0,XX) * SQRT(1.0 - EXP( -4.0 * XX * XX / PI ) )
       ERFC(XX) = 1.0 - ERF(XX)
C     ::::::::::::::::::::::::::::::::::::::::

C ---------------------------- Begin code -------------------------------

      IF( FIRSTIME ) THEN
      
C *** Compute these once and they will all be saved in COMMON

c *** coarse mode has a fixed standard deviation

         XXLSGCO = LOG( SGINICO)

         EC1    = EXP( 0.125 * XXLSGCO ** 2 )
         ESC04  = EC1 ** 4
         ESC08  = ESC04 * ESC04
         ESC12  = ESC04 * ESC04 * ESC04
         ESC16  = ESC08 * ESC08
         ESC20  = ESC16 * ESC04
         ESC24  = ESC12 * ESC12
         ESC28  = ESC20 * ESC08
         ESC32  = ESC16 * ESC16
         ESC36  = ESC16 * ESC20
         ESC64  = ESC32 * ESC32
         ESCM20 = 1.0 / ESC20
         ESCM32 = 1.0 / ESC32

C *** calculate minimum values for number and 2nd moment.

         NUMMIN_AT = SO4FAC * AEROCONCMIN_AT /
     &            ( DGINIAT ** 3 * EXP(4.5 * LOG(SGINIAT)**2))

         NUMMIN_AC = SO4FAC * AEROCONCMIN_AC /
     &            ( DGINIAC ** 3 * EXP(4.5 * LOG(SGINIAC)**2))

         M2MIN_AT = NUMMIN_AT * DGINIAT ** 2 *
     &              EXP(2.0 * LOG(SGINIAT)**2)

         M2MIN_AC = NUMMIN_AC * DGINIAC ** 2 *
     &              EXP(2.0 * LOG(SGINIAC)**2)

         NUMMIN_C = ANTHFAC * AEROCONCMIN_CO / ( DGINICO**3 * ESC36)

C *** set constant part of molecular velocities

         COFCBAR_SO4 = SQRT( 8.0 * RGASUNIV/( PI * MWH2SO4 ))
         COFCBAR_ORG = SQRT( 8.0 * RGASUNIV/( PI * MWORG ))

         FIRSTIME = .FALSE.

      END IF  ! firstime

C *** calculate square root of the ambient temperature for later use

      SQRT_TEMP = SQRT( AIRTEMP)

C *** Calculate mean free path [ m ]:
C     6.6328E-8 is the sea level value given in Table I.2.8
C     on page 10 of U.S. Standard Atmosphere 1962

      XLM = 6.6328E-8 * P0 * AIRTEMP  / ( T0 * AIRPRS )

C *** Calculate dynamic viscosity [ kg m**-1 s**-1 ]:
C     U.S. Standard Atmosphere 1962 page 14 expression
C     for dynamic viscosity is:
c     dynamic viscosity =  beta * T * sqrt(T) / ( T + S)
c     where beta = 1.458e-6 [ kg sec^-1 K**-0.5 ], s = 110.4 [ K ].

      AMU = 1.458E-6 * AIRTEMP * SQRT_TEMP / ( AIRTEMP + 110.4 )

C...set minimums for coarse mode

      CBLK( VANTHA ) = MAX( AEROCONCMIN_CO, CBLK( VANTHA ) )
      CBLK( VCOR0 )  = MAX( REAL(NUMMIN_C), CBLK( VCOR0 ) )

c *** Outside AEROPROC, the surface area concentration [ m**2 / m**3 ]
c     of the Aitken and accumulation modes is tracked.  Internally,
c     the 2nd moment concentration is the variable of interest.
c       Surface area = pi * 2nd moment.

      CBLK(VAT2) = CBLK(VSURFAT) / PI
      CBLK(VAC2) = CBLK(VSURFAC) / PI

C *** Update the third moments, geometric mean diameters, geometric 
c     standard deviations, modal mass totals, and modal particle 
c     densities, based on the concentrations of M2, M0, and speciated 
c     masses resulting from transport, cloud processing, and gas-phase 
c     chemistry.  Ignore H2O and SOA (i.e., M3_WET_FLAG = .FALSE.) 
c     because those species are not transported with the 2nd moment and
c     their concentrations will be updated subsequently to establish
c     gas/particle equilibrium.

      M3_WET_FLAG = .FALSE.

      CALL GETPAR( NSPCSDA, CBLK,
     &             PMASSAT, PMASSAC, PMASSCO,
     &             PDENSAT, PDENSAC, PDENSCO,
     &             DGATK, DGACC, DGCOR,
     &             XXLSGAT, XXLSGAC,
     &             M3_WET_FLAG    )

c *** SECONDARY ORGANICS
c      Update the secondary organic aerosol (SOA) mass concentrations
c     by equilibration of the aerosol with the appropriate SVOCs.
c     Partitioning of SOA to each mode is in proportion to the total
c     organic aerosol mass (primary + secondary) in each mode and 
c     is assumed to have no affect on the geometric standard deviations
c     of the modes.

c *** Calculate the fractions of organic mass in each mode

      OMASS_I = CBLK( VORGAI ) + CBLK( VORGPAI ) + CBLK( VORGBAI )
      OMASS_J = CBLK( VORGAJ ) + CBLK( VORGPAJ ) + CBLK( VORGBAJ )
      FRACI   = OMASS_I / ( OMASS_I + OMASS_J)
      FRACJ   = MAX( 0.0, 1.0 - FRACI )

C *** Set inputs for ORGAER3

      OLDSOA_A = CBLK( VORGAI )  + CBLK( VORGAJ )
      OLDSOA_B = CBLK( VORGBAI ) + CBLK( VORGBAJ )
      OLDPOA   = CBLK( VORGPAI ) + CBLK( VORGPAJ )

C *** Equilibrate SVOCs with the aerosol phase

      call ORGAER3(DT, LAYER, AIRTEMP, AIRPRS,
     &             ORGPROD, NPSPCS, VAPORS, NCVAP,
     &             OLDSOA_A, OLDSOA_B, OLDPOA,
     &             SOA_A, SOA_B, SUMPROD, LOGDEV )

c *** Calculate average production rates for this time step from the
c     values of SOA_A and SOA_B.  These rates are non-negative.

      ORGRATE  = SOA_A / DT
      ORGBRATE = SOA_B / DT

c *** Repartition SOA_A, and SOA_B to the modes in the
c     same proportion as the original modal concentrations

      CBLK( VORGAI )  = MAX( FRACI * SOA_A, CONMIN )
      CBLK( VORGBAI ) = MAX( FRACI * SOA_B, CONMIN )
      
      CBLK( VORGAJ )  = MAX( FRACJ * SOA_A, CONMIN )
      CBLK( VORGBAJ ) = MAX( FRACJ * SOA_B, CONMIN )

C *** Update 3rd moments

      OLD_M2_I = CBLK( VAT2 )
      OLD_M2_J = CBLK( VAC2 )
      OLD_M3_I = CBLK( VAT3 )
      OLD_M3_J = CBLK( VAC3 )

      NEW_M3_I = OLD_M3_I +
     &           ORGFAC * CBLK( VORGAI  ) +
     &           ORGFAC * CBLK( VORGBAI )
      CBLK( VAT3 ) = NEW_M3_I

      NEW_M3_J = OLD_M3_J +
     &           ORGFAC * CBLK( VORGAJ  ) +
     &           ORGFAC * CBLK( VORGBAJ )
      CBLK( VAC3 ) = NEW_M3_J

C *** Update 2nd moments and geometric mean diameters assuming that 
c     the SOA condensation/evaporation does not affect the geometric
c     standard deviations

      NEW_M2_I = OLD_M2_I * ( NEW_M3_I / OLD_M3_I ) ** TWO3
      CBLK( VAT2 ) = NEW_M2_I

      NEW_M2_J = OLD_M2_J * ( NEW_M3_J / OLD_M3_J ) ** TWO3
      CBLK( VAC2 ) = NEW_M2_J

      CBLK( VDGAT ) = CBLK( VDGAT ) * SQRT( NEW_M2_I / OLD_M2_I )
      CBLK( VDGAC ) = CBLK( VDGAC ) * SQRT( NEW_M2_J / OLD_M2_J )      
      

c *** HETEROGENEOUS CHEMISTRY
c      Update the N2O5 and HNO3 concentrations to account for 
c     heterogenous nitrate formation.

      CALL HETCHEM ( NSPCSDA, CBLK, GAMMA_N2O5, AIRTEMP, AIRRH, 
     &               AIRPRS, DT )
      

c *** WATER, AMMONIUM, and NITRATE
c      Update the H2O, NH4, and NO3 concentrations to account for 
c     heterogenous NO3 production and thermodynamic equilibrium with
c     gas-phase species.  See notes in SUBROUTINE EQL3.

      CALL EQL3( NSPCSDA, CBLK, SO4RATE, AIRTEMP, AIRRH, DT )


C *** Update the third moments, geometric mean diameters, geometric 
c     standard deviations, modal mass totals, and modal particle 
c     densities, based on the new 2nd moment and speciated mass 
c     concentrations calculated in SUBROUTINE EQL3.  Include H2O
c     and SOA in these updates (i.e., M3_WET_FLAG = .TRUE.) because
c     the wet distribution is needed when solving the modal dynamics
c     equations.

      M3_WET_FLAG = .true.

      CALL GETPAR( NSPCSDA, CBLK,
     &             PMASSAT, PMASSAC, PMASSCO,
     &             PDENSAT, PDENSAC, PDENSCO,
     &             DGATK, DGACC, DGCOR,
     &             XXLSGAT, XXLSGAC,
     &             M3_WET_FLAG     )


C *** CONDENSATIONAL GROWTH (part 1)
c      Calculate intermediate variables needed to determine the 2nd and
c     3rd moment condensational-growth rates.  3rd moment terms are 
c     needed for the calculation of new particle production.  See 
c     Section 3.3 of Jiang & Roth (2003) for a detailed discussion.

C *** Correct diffusivities for temperature and pressure
      DIFFCORR = ( P0 / AIRPRS ) * ( AIRTEMP / 273.16 ) ** 1.75
      DV_SO4 = DIFFSULF * DIFFCORR
      DV_ORG = DIFFORG  * DIFFCORR

C *** Calculate molecular velocities (temperature dependant)
      CBAR_SO4 = COFCBAR_SO4 * SQRT_TEMP
      CBAR_ORG = COFCBAR_ORG * SQRT_TEMP

C *** Set zeroth moments
      AM0AT = CBLK(VAT0)
      AM0AC = CBLK(VAC0)

C *** Calculate first moments using Equation 4 of Binkowski & Shankar
c     (1995) or Equation 3 of Binkowski and Roselle (2003).

      AM1AT = CBLK( VAT0 ) * DGATK * EXP( 0.5 * XXLSGAT * XXLSGAT )
      AM1AC = CBLK( VAC0 ) * DGACC * EXP( 0.5 * XXLSGAC * XXLSGAC )

C *** Set second moments
      AM2AT = CBLK(VAT2)
      AM2AC = CBLK(VAC2)
      
C *** Calculate the size-dependent terms in the condensational-
c     growth factor expressions for sulfate and organics using 
c     Equations A13-A14 of Binkowski & Shankar (1995).  FCONCAT_* and
c     FCONCAC_* are the "I_ki" and "I_kj" terms in Equations 7a and 7b
c     of Binkowski & Shankar (1995).

      CALL HCOND3(AM0AT,AM1AT,AM2AT,DV_SO4,ALPHSULF, CBAR_SO4, 
     &            FCONCAT_SO4 )   ! Equation A13
      CALL HCOND3(AM0AC,AM1AC,AM2AC,DV_SO4,ALPHSULF, CBAR_SO4, 
     &            FCONCAC_SO4 )   ! Equation A14
      CALL HCOND3(AM0AT,AM1AT,AM2AT,DV_ORG,ALPHORG, CBAR_ORG,
     &            FCONCAT_ORG )   ! Equation A13
      CALL HCOND3(AM0AC,AM1AC,AM2AC,DV_ORG,ALPHORG, CBAR_ORG,
     &            FCONCAC_ORG )   ! Equation A14

c *** Calculate the fraction of condensing material injected into each
c     mode.  These are the "omega" factors in Equations 7a and 7b of
c     Binkowski & Shankar (1995).  The i-mode factors are calculated
c     using Equation A17 of Binkowski & Shankar (1995).  The j-mode 
c     factors are calculated by difference, to avoid mass conservation
c     violations arising from numerical error.
 
      FCONCM1_SO4  = 1.0 / ( FCONCAT_SO4(2) + FCONCAC_SO4(2) )
      OMEGA_AT_SO4 = FCONCM1_SO4 * FCONCAT_SO4(2)
      OMEGA_AC_SO4 = 1.0 - OMEGA_AT_SO4
            
      FCONCM1_ORG  = 1.0 / ( FCONCAT_ORG(2) + FCONCAC_ORG(2) )
      OMEGA_AT_ORG = FCONCM1_ORG * FCONCAT_ORG(2)
      OMEGA_AC_ORG = 1.0 - OMEGA_AT_ORG


C *** NEW PARTICLE PRODUCTION
c      Calculate the new particle production rate due to binary
c     nucleation of H2O and H2SO4.  These calculations are performed 
c     only when the gas-phase production rate of H2SO4 (i.e., SO4RATE) 
c     is non-zero.  The condensation rate of H2SO4 is calculated as the
c     gas-phase production rate minus the new particle production rate.

C *** initialize variables
      DMDT_so4  = 0.0D0
      DNDT      = 0.0D0
      DM2DT     = 0.0D0
      SCONDRATE = 0.0

c *** produce new particles only during time steps when the gas-phase 
c     production rate of H2SO4 is non-zero

      IF( SO4RATE .NE. 0.0 ) THEN

c *** adjust sulfuric acid vapor concentration to a value in
c     equilibrium with the production of new particles and the
c     condensation of sulfuric acid vapor on existing particles, based 
c     on Equations A21 and A23 of Binkowski & Shankar (1995).

       XH2SO4 = SO4RATE / ( FCONCAT_SO4( 2 ) + FCONCAC_SO4( 2 ) )
       XH2SO4 = MAX(XH2SO4, 1.0E-30)
       CBLK( VSULF ) = XH2SO4

c *** calculate new particle production rate for 0th, 2nd,
c     & 3rd moments

       CALL NEWPART3 ( AIRRH, AIRTEMP, XH2SO4, SO4RATE,
     &                 DNDT, DMDT_so4, DM2DT )

C *** calculate sulfate condensation rate as the gas-phase production 
c     rate minus the new particle production rate, following Equation
c     3.23 of Jiang & Roth (2003).

       SCONDRATE = MAX( SO4RATE - DMDT_so4, 0.0D0 )

      END IF ! check on SO4RATE


C *** CONDENSATIONAL GROWTH (part 2)
c      Calculate the 2nd and 3rd moment condensational-growth rates.

C *** calculate third moment production rates [ m**3 m**-3 s-1], based 
c     on Equation 8 of Binkowski & Shankar (1995)

      CHEMRATE_SO4 = SO4FAC * SCONDRATE
      CHEMRATE_ORG = ORGFAC * ( ORGRATE + ORGBRATE )

c *** calculate 2nd moment condensational-growth rates due to SO4
c      Note: 2nd moment is already adjusted for SOA, NH4, NO3, & H2O
c
c     The following code exactly implements equations (7a) & (7b)
c     of Binkowski & Shankar (1995) but is inefficient in requiring
c     two extra divisions
c       CGRAT(1) = TWO3 * CHEMRATE_SO4
c     &                 * OMEGA_AT * FCONCAT_SO4(1) / FCONCAT_SO4(2)
c       CGRAC(1) = TWO3 * CHEMRATE_SO4
c     &                 * OMEGA_AC * FCONCAC_SO4(1) / FCONCAC_SO4(2)
c
c     After some simple algebra, the following code is equivalent

      CGRAT(1) = TWO3 * CHEMRATE_SO4 * FCONCAT_SO4(1) * FCONCM1_SO4
      CGRAC(1) = TWO3 * CHEMRATE_SO4 * FCONCAC_SO4(1) * FCONCM1_SO4

C *** calculate 3rd moment condensational-growth rates using Equations
c     7a and 7b of Binkowski & Shankar (1995).  The 3rd moment growth
c     rates are used only to check whether mode merging is necessary.
c
c      Note: 3rd moment includes organics for determination of possible
c      mode merging, whereas 2nd moment is already adjusted for organics
c      3rd moment growth rates due to NH4, NO3, & H2O are not currently
c      considered in the determination of possible mode merging

      CGRAT(2) = CHEMRATE_SO4 * OMEGA_AT_SO4
     &         + CHEMRATE_ORG * OMEGA_AT_ORG
      CGRAC(2) = CHEMRATE_SO4 * OMEGA_AC_SO4
     &         + CHEMRATE_ORG * OMEGA_AC_ORG


C *** COAGULATION
c      Calculate coagulation coefficients using a method dictated by
c     the value of FASTCOAG_FLAG.  If TRUE, the computationally-
c     efficient GETCOAGS routine is used.  If FALSE, the more intensive
c     Gauss-Hermite numerical quadrature method is used.  See Section 
c     2.1 of Bhave et al. (2004) for further discussion.

      FASTCOAG_FLAG = .TRUE.

C *** set atmospheric mean free path in double precision
      LAMDA    = XLM

C *** calculate term used in Equation A6 of Binkowski & Shankar (1995)
      KNC      = TWO3 * BOLTZ *  AIRTEMP / AMU

C *** calculate terms used in Equation A5 of Binkowski & Shankar (1995)
      KFMAT    = SQRT( 3.0 * BOLTZ * AIRTEMP / PDENSAT )
      KFMAC    = SQRT( 3.0 * BOLTZ * AIRTEMP / PDENSAC )
      KFMATAC  = SQRT( 6.0 * BOLTZ * AIRTEMP / ( PDENSAT + PDENSAC ) )
      
C *** transfer of number to accumulation mode from Aitken mode is zero
      BACAT(1) = 0.0

      IF ( FASTCOAG_FLAG ) THEN ! Solve coagulation analytically

C *** set geometric mean diameters, geometric standard deviations, and
c     ln(GSD) in double precision
        DGATK_D = DGATK
        DGACC_D = DGACC
        XXLSGAT_D = XXLSGAT
        XXLSGAC_D = XXLSGAC
        SGATK_D = EXP( XXLSGAT )
        SGACC_D = EXP( XXLSGAC )

c *** calculate intermodal and intramodal coagulation coefficients
c     for zeroth and second moments, and intermodal coagulation 
c     coefficient for third moment
        CALL GETCOAGS( LAMDA, KFMATAC, KFMAT, KFMAC, KNC,
     &                 DGATK_D, DGACC_D, SGATK_D, SGACC_D, 
     &                 XXLSGAT_D,XXLSGAC_D, 
     &                 BATAT(2), BATAT(1), BACAC(2), BACAC(1),
     &                 BATAC(2), BACAT(2), BATAC(1), C3IJ )
      
      ELSE                 ! Use Gauss-Hermite numerical quadrature

c *** calculate Aitken-mode intramodal coagulation coefficients
c     for zeroth and second moments
        CALL INTRACOAG_GH(LAMDA, KFMAT, KNC, DGATK,
     &                    XXLSGAT, BATAT(2), BATAT(1) )

C *** calculate accumulation-mode intramodal coagulation coefficients
c     for zeroth and second moments
        CALL INTRACOAG_GH(LAMDA, KFMAC, KNC, DGACC,
     &                    XXLSGAC, BACAC(2), BACAC(1) )

C *** calculate intermodal coagulation coefficients for zeroth, second, 
c     and third moments
        CALL INTERCOAG_GH(LAMDA, KFMATAC, KNC, DGATK, DGACC,
     &                    XXLSGAT, XXLSGAC,
     &                    BATAC(2), BACAT(2), BATAC(1), C3IJ )
      
      END IF

C *** calculate 3rd moment intermodal transfer rate by coagulation
      C30ATAC = C3IJ * CBLK( VAC0 ) * CBLK( VAT0 )


C *** TAKE ONE FORWARD TIME STEP - Solve Modal Dynamics Equations
c     This code implements Section 1.4 of Binkowski and Roselle (2003)
c     with one notable exception.  Because emissions are treated in 
c     CMAQ's vertical diffusion routine, they do not appear in the 
c     following equations.
c
c     M2 is updated before M0 because the intermodal transfer rate of 
c     M2 is a function of the number concentrations.  In contrast, 
c     production and loss rates of M0 are independent of M2.  Advancing
c     M2 before M0 avoids operator splitting within the modal-dynamic-
c     equation solution.  A similar rearrangement would be necessary
c     for the M3 update, but the dependence of M3 on number 
c     concentrations already is accounted for in the C30ATAC term.

C *** UPDATE SECOND MOMENT
C     For each lognormal mode, solve equations of form:
C       dM2/dt = P2 - L2*M2   ! if L2 > 0
C     with solution
c       M2(t) = P2/L2 + ( M2(t0) - P2/L2 ) * exp( -L2*dt )
c     OR
c       dM2/dt = P2           ! if L2 = 0
c     with solution
c       M2(t) = M2(t0) + P2*dt

C *** Aitken mode: initial value of M2

      M20 = CBLK( VAT2 )

c *** Production of 2nd moment in Aitken mode is due to new 
c     particle production and condensational growth

      PROD = CGRAT(1) + DM2DT

C *** Loss of 2nd moment from Aitken mode is due to intermodal
c     coagulation with accumulation mode and intramodal coagulation

      LOSS = (
     &       (
     &       BATAT(2) * CBLK( VAT0 ) +
     &       BATAC(2) * CBLK( VAC0 )
     &       ) * CBLK( VAT0 )
     &       ) / M20

C *** Solve for M2_Aitken based on PROD and LOSS during this time step
c     Note: LOSS is assumed to be non-negative.

      IF ( LOSS .GT. 0.0) THEN
        POL = PROD / LOSS
        Y = POL + ( M20 - POL ) * EXP( -LOSS * DT )
      ELSE
        Y = M20 + PROD * DT
      END IF ! test on loss

C *** Transfer new value of M2_Aitken to the CBLK array

      CBLK( VAT2 ) = MAX( M2MIN_AT, Y )


C *** Accumulation mode: initial value of M2

      M20 = CBLK( VAC2 )

c *** Production of 2nd moment in accumulation mode is due to 
c     condensational growth and intermodal coagulation with 
c     Aitken mode

      PROD = CGRAC(1) +
     &       BACAT(2) * CBLK( VAC0 ) * CBLK( VAT0 )

C *** Loss of 2nd moment from accumulation mode is due only to 
c     intramodal coagulation

      LOSS =  ( BACAC(2) * CBLK(VAC0) * CBLK(VAC0) ) / M20

C *** Solve for M2_accum based on PROD and LOSS during this time step
c     Note: LOSS is assumed to be non-negative.

      IF ( LOSS .GT. 0.0) THEN
        POL = PROD / LOSS
        Y = POL + ( M20 - POL ) * EXP( -LOSS * DT )
      ELSE
        Y = M20 + PROD * DT
      END IF ! test on loss

C *** Transfer new value of M2_accum to CBLK array

      CBLK( VAC2 ) = MAX( M2MIN_AC, Y )

c *** end of update for second moment


C *** UPDATE ZEROTH MOMENT (i.e. number concentration)

C *** Aitken mode: initial value of M0

      Y0 = CBLK( VAT0 )

c     The rate of change for M0_Aitken is described in Equation 8a of 
c     Binkowski & Roselle (2003).  This is a Riccati-type equation 
c     of the form:
c       dY/dt = C - A * Y**2 - B * Y
c     where

      A = BATAT(1)                 ! intramodal coagulation
      B = BATAC(1) * CBLK( VAC0 )  ! intermodal coagulation
      C = DNDT                     ! new particle production

C *** Solve for M0_Aitken using analytical solutions to the Riccati 
c     equation.  As described in Paragraph 19 of Binkowski & Roselle
c     (2003), the Riccati equation has different analytical solutions 
c     depending on whether the coefficient C is zero or non-zero.
c     Note: C is assumed to be non-negative in the following code.

      IF( C .GT. 0.0D0 ) THEN
c       Note: The Riccati-equation solution for the case where C is 
c             non-zero is printed incorrectly in Binkowski & Roselle
c             (2003).  The exp(delta*t) should read exp(-delta*t).
c             The r2 term should be defined as -(b+delta)/2 instead
c             of (b+delta)/2.  The following code implements the
c             solution correctly.  The comments in the right-hand
c             margin correspond to names assigned to these terms
c             in Binkowski & Roselle (2003).

        DHAT = SQRT( B * B + 4.0D0 * A * C )      ! delta
        M1 = 2.0D0 * A * C / ( B + DHAT )         ! r1
        M2 = - 0.5D0 * ( B + DHAT )               ! r2
        P =  - ( M1 - A  * Y0 ) / (M2 - A * Y0 )  ! gamma
        PEXPDT = P * EXP( -DHAT * DT )
        Y = (  M1 + M2 * PEXPDT ) /
     &        ( A * (1.0D0 + PEXPDT ) )           ! N(t) solution

      ELSE ! for c = 0.0d0 solution simplifies to the following:

        EXPDT = EXP( - B * DT )
        IF ( EXPDT .LT. 1.0D0 ) THEN
          Y = B * Y0 * EXPDT / ( B + A * Y0 * (1.0D0 - EXPDT ) )
        ELSE
          Y = Y0    ! solution in the limit that B approaches zero
        END IF

      END IF

C *** Transfer new value of M0_Aitken to the CBLK array

      CBLK( VAT0 )  = MAX( NUMMIN_AT, Y )

C *** Accumulation mode: initial value of M0

      Y0 = CBLK( VAC0 )

c     The rate of change for M0_accum is described in Equation 8b of 
c     Binkowski & Roselle (2003), except the coefficient C is zero
c     because emissions are treated outside the CMAQ aerosol module.
c     The equation reduces to the form: dY/dt = -A * Y**2 , where

      A = BACAC(1)                 ! intramodal coagulation

C *** Solve for M0_accum using Smoluchowski's solution

      Y = Y0 / ( ONE + A * Y0 * DT)

C *** Transfer new value of M0_accum to the CBLK array

      CBLK( VAC0 ) = MAX( NUMMIN_AC, Y )

c *** end of update for zeroth moment


C *** UPDATE MASS CONCENTRATIONS (for each species)
c     The following procedure is described in Paragraphs 21-23 
c     of Binkowski & Roselle (2003), except the Ei,n and Ej,n terms 
c     are excluded here because emissions are treated outside the 
c     CMAQ aerosol module.

c     Aitken mode mass concentration rates of change are of the form:
c       dc/dt = P - L*c    ! Equation 9a of Binkowski & Roselle (2003)
c     with solution
c       c(t0 + dt) = P/L + ( c(t0) - P/L ) * exp(-L*dt)

c     For all species, loss of Aitken mode mass is due to intermodal 
c     coagulation.
c       LOSSn = PI/6 * RHOn * C30ATAC / MASSn
c       RHOn  = MASSn / (M3 * PI/6)
c     When above equations are combined, the PI/6 terms cancel yielding
c       LOSSn = C30ATAC / M3
c     where LOSSn is the loss rate of species n, RHOn is the mass of 
c     species n per unit of particle volume, C30ATAC is the 3rd moment
c     loss rate due to intermodal coagulation, MASSn is the mass 
c     concentration of species n, and M3 is the 3rd moment 
c     concentration.

      LOSS =  C30ATAC / CBLK( VAT3 )

c     Setup extra variables to solve for Aitken mode mass concentrations

      FACTRANS = LOSS * DT
      EXPDT = EXP(- FACTRANS)
      LOSSINV = 1.0 / LOSS

C *** SULFATE

      OLDSO4 = CBLK( VSO4AI ) + CBLK( VSO4AJ )
      NEWSO4 =( SCONDRATE + DMDT_so4 ) * DT
      TSO4 =  OLDSO4 + NEWSO4

c     Production of SO4_Aitken is from nucleation and condensation

      PROD = SCONDRATE * OMEGA_AT_SO4 + DMDT_so4

c     Solve for SO4_Aitken.  Calculate SO4_accum by difference.

      POL = PROD * LOSSINV
      CBLK( VSO4AI ) = POL + ( CBLK( VSO4AI ) - POL ) * EXPDT
      CBLK( VSO4AI ) = MAX( AEROCONCMIN_AT, CBLK( VSO4AI ) )
      CBLK( VSO4AJ ) = TSO4 - CBLK ( VSO4AI )
      CBLK( VSO4AJ ) = MAX(AEROCONCMIN_AC, CBLK( VSO4AJ ) )

C *** SECONDARY INORGANICS: AMMONIUM, NITRATE, and WATER
c     These species are produced by condensation, but their mass
c     concentrations have been updated in subroutine EQL3.  Therefore,
c     only the mass transfer due to intermodal coagulation is 
c     treated here.

      TMASS = CBLK( VNH4AI ) + CBLK( VNH4AJ )
      CBLK( VNH4AI ) = CBLK( VNH4AI ) * EXPDT
      CBLK( VNH4AJ ) = MAX( CONMIN, TMASS - CBLK( VNH4AI ) )

      TMASS = CBLK( VNO3AI ) + CBLK( VNO3AJ )
      CBLK( VNO3AI ) = CBLK( VNO3AI )  * EXPDT
      CBLK( VNO3AJ ) = MAX( CONMIN, TMASS - CBLK( VNO3AI ) )

      TMASS = CBLK( VH2OAI ) + CBLK( VH2OAJ )
      CBLK( VH2OAI ) = CBLK( VH2OAI ) * EXPDT
      CBLK( VH2OAJ ) = MAX( CONMIN, TMASS - CBLK( VH2OAI) )

C *** SECONDARY ORGANICS: ANTHROPOGENIC & BIOGENIC
c     These species are produced by condensation, but their mass
c     concentrations were updated immediately after the call to
c     subroutine ORGAER3.  Therefore, only the mass transfer due to 
c     intermodal coagulation is treated here.

      TMASS = SOA_A
      CBLK( VORGAI) = CBLK( VORGAI) * EXPDT
      CBLK( VORGAI) = MAX( CONMIN, CBLK( VORGAI ) )
      CBLK( VORGAJ) = MAX( CONMIN, TMASS - CBLK( VORGAI ) )

      TMASS = SOA_B
      CBLK( VORGBAI) = CBLK( VORGBAI) * EXPDT 
      CBLK( VORGBAI) = MAX( CONMIN, CBLK( VORGBAI ) )
      CBLK( VORGBAJ) = MAX( CONMIN, TMASS - CBLK( VORGBAI ) )

C *** PRIMARY SPECIES: ANTHROPOGENIC ORGANICS, EC, and OTHER PM2.5
c     Because emissions are treated in a different module, there is 
c     no production term for these species.  Only the mass transfer 
c     due to intermodal coagulation is treated here.

      TMASS = CBLK( VORGPAI ) + CBLK( VORGPAJ)
      CBLK( VORGPAI ) = CBLK( VORGPAI ) * EXPDT
      CBLK( VORGPAI ) = MAX(CONMIN, CBLK( VORGPAI ) )
      CBLK( VORGPAJ ) = MAX( CONMIN, TMASS - CBLK( VORGPAI ) )

      TMASS = CBLK( VECI) + CBLK( VECJ )
      CBLK( VECI ) = CBLK( VECI ) * EXPDT
      CBLK( VECI ) = MAX(CONMIN, CBLK( VECI ) )
      CBLK( VECJ ) = MAX( CONMIN, TMASS - CBLK( VECI ) )

      TMASS = CBLK( VP25AI ) + CBLK( VP25AJ)
      CBLK( VP25AI ) = CBLK( VP25AI) * EXPDT
      CBLK( VP25AI ) = MAX( CONMIN, CBLK( VP25AI) )
      CBLK( VP25AJ ) = MAX( CONMIN, TMASS - CBLK( VP25AI ) )

c *** end of update for species mass concentrations


C *** MODE MERGING
c     This code implements Section 1.5 of Binkowski and Roselle (2003).
c     If the Aitken mode mass is growing faster than accumulation mode
c     mass and the Aitken mode number concentration exceeds the 
c     accumulation mode number concentration, then modes are merged by
c     renaming.

      IF( CGRAT(2) .GT. CGRAC(2) .AND.
     &    CBLK( VAT0) .GT. CBLK( VAC0) ) THEN

C *** Before mode merging, update the third moments, geometric mean
c     diameters, geometric standard deviations, modal mass totals, and
c     particle densities, based on the new concentrations of M2, M0, and
c     speciated masses calculated above.

           CALL GETPAR( NSPCSDA, CBLK,
     &                  PMASSAT, PMASSAC, PMASSCO,
     &                  PDENSAT, PDENSAC, PDENSCO,
     &                  DGATK, DGACC, DGCOR,
     &                  XXLSGAT, XXLSGAC,
     &                  M3_WET_FLAG    )

C *** Calculate AAA = ln( Dij / DGATK ) / ( SQRT2 * XXLSGAT ), where Dij
c     is the diameter at which the Aitken-mode and accumulation-mode 
c     number distributions intersect (i.e., overlap).  AAA is equivalent
c     to the "Xnum" term described below Equation 10a by Binkowski and 
c     Roselle (2003).

           AAA = getaf(
     &            CBLK(  VAT0),
     &            CBLK(  VAC0),
     &            DGATK ,
     &            DGACC ,
     &            XXLSGAT ,
     &            XXLSGAC ,
     &            SQRT2       )

C *** Ensure that Xnum is large enough so that no more than half of
c     the Aitken mode mass is merged into the accumulation mode during 
c     any given time step.  This criterion is described in Paragraph 26
c     of Binkowski and Roselle (2003).

           XXM3 = 3.0 * XXLSGAT  / SQRT2
           XNUM = MAX( AAA, XXM3 )

C *** Factors used in error function calls for M2 and M3 mode merging

           XXM2 = TWO3 * XXM3
           XM2  = XNUM - XXM2 ! set up for 2nd moment transfer
           XM3  = XNUM - XXM3 ! set up for 3rd moment and mass transfers

C *** Calculate the fractions of the number, 2nd, and 3rd moment 
c     distributions with diameter greater than the intersection diameter
           FNUM  = 0.5 * ERFC(XNUM)            ! Eq 10a of B&R 2003
           FM2   = 0.5 * ERFC(XM2)             ! Eq 10b of B&R 2003
           FM3   = 0.5 * ERFC(XM3)             ! Eq 10b of B&R 2003

C *** Calculate the fractions of the number, 2nd, and 3rd moment
c     distributions with diameters less than the intersection diameter.
           PHNUM = 0.5 * ( 1.0 + ERF( XNUM) )  ! Eq 10c of B&R 2003
           PHM2  = 0.5 * ( 1.0 + ERF( XM2 ) )  ! Eq 10d of B&R 2003
           PHM3  = 0.5 * ( 1.0 + ERF( XM3 ) )  ! Eq 10d of B&R 2003

C *** Update accumulation-mode moment concentrations using
c     Equations 11a - 11c of Binkowski and Roselle (2003).

           CBLK( VAC0 ) = CBLK( VAC0 ) + CBLK( VAT0 ) * FNUM
           CBLK( VAC2 ) = CBLK( VAC2 ) + CBLK( VAT2 ) * FM2
           CBLK( VAC3 ) = CBLK( VAC3 ) + CBLK( VAT3 ) * FM3

C *** Update Aitken-mode moment concentrations using
c     Equations 11d - 11f of Binkowski and Roselle (2003).

           CBLK( VAT0 ) = CBLK( VAT0 ) * PHNUM
           CBLK( VAT2 ) = CBLK( VAT2 ) * PHM2
           CBLK( VAT3 ) = CBLK( VAT3 ) * PHM3

C *** Rename masses of each species from Aitken mode to acumulation mode
c     using Equation 11b of Binkowski and Roselle (2003).

           CBLK( VSO4AJ ) = CBLK( VSO4AJ)   + CBLK( VSO4AI )   * FM3
           CBLK( VNH4AJ ) = CBLK( VNH4AJ )  + CBLK( VNH4AI )   * FM3
           CBLK( VNO3AJ ) = CBLK( VNO3AJ )  + CBLK( VNO3AI )   * FM3
           CBLK( VH2OAJ ) = CBLK( VH2OAJ )  + CBLK( VH2OAI )   * FM3
           CBLK( VORGAJ ) = CBLK( VORGAJ )  + CBLK( VORGAI )   * FM3
           CBLK( VORGPAJ )= CBLK( VORGPAJ ) + CBLK( VORGPAI )  * FM3
           CBLK( VORGBAJ )= CBLK( VORGBAJ ) + CBLK( VORGBAI )  * FM3
           CBLK( VP25AJ ) = CBLK( VP25AJ )  + CBLK( VP25AI )   * FM3
           CBLK( VECJ )   = CBLK( VECJ )    + CBLK( VECI )     * FM3

C *** Update Aitken-mode species masses for loss to accumulation mode 
c     using Equation 11e of Binkowski and Roselle (2003).

           CBLK( VSO4AI )   = CBLK( VSO4AI )   * PHM3
           CBLK( VNH4AI )   = CBLK( VNH4AI )   * PHM3
           CBLK( VNO3AI )   = CBLK( VNO3AI )   * PHM3
           CBLK( VH2OAI )   = CBLK( VH2OAI )   * PHM3
           CBLK( VORGAI )   = CBLK( VORGAI )   * PHM3
           CBLK( VORGPAI )  = CBLK( VORGPAI )  * PHM3
           CBLK( VORGBAI )  = CBLK( VORGBAI )  * PHM3
           CBLK( VP25AI )   = CBLK( VP25AI )   * PHM3
           CBLK( VECI )     = CBLK( VECI)      * PHM3

      END IF ! end check on necessity for merging

c *** end of update for mode merging


C *** Update the third moments, geometric mean diameters, geometric 
c     standard deviations, modal mass totals, and particle densities,
c     based on the final concentrations of M2, M0, and speciated masses
c     after mode merging is complete.

      CALL GETPAR( NSPCSDA, CBLK,
     &             PMASSAT, PMASSAC, PMASSCO,
     &             PDENSAT, PDENSAC, PDENSCO,
     &             DGATK, DGACC, DGCOR,
     &             XXLSGAT, XXLSGAC,
     &             M3_WET_FLAG    )


C *** Set minimum value for all concentrations in the CBLK array

      DO SPC = 1, NSPCSDA
            CBLK( SPC )  = MAX( CBLK( SPC ), CONMIN )
      END DO


C *** Update surface area concentration by multiplying 2nd moment by PI

      CBLK(VSURFAT) =  PI * CBLK(VAT2)
      CBLK(VSURFAC) =  PI * CBLK(VAC2)
      
      
      RETURN

      END  SUBROUTINE AEROPROC


c /////////////////////////////////////////////////////////////////////
c  SUBROUTINE EQL3 calculates the heterogeneous production of HNO3 from
c    N2O5 and then distributes the ammonia/ammonium, nitric 
c    acid/nitrate, and water, between the gas and aerosol phases as a 
c    function of total sulfate, total ammonia, total nitrate, relative
c    humidity, and temperature.
c
c  KEY SUBROUTINES CALLED: ISOROPIA
c
c  KEY FUNCTIONS CALLED: none
c
c  REVISION HISTORY:
c     Prototype coded in January 1995 by Uma Shankar and Carlie Coates
c
c US  08/1995  Revised to calculate air density in statement function
c              and collect met variable stmt funcs in one include file
c
c FSB 07/26/96 Revised to use the block concept
c
c FSB 12/18/96 Revised to do do i-mode calculation
c
c FSB 04/19/99 Revised to use new RPMARES algorithm
c
c DW  01/21/99 Revised subroutine ACTCOF --David Wong
c
c FSB 07/06/99 Revised to include calculation of third moments and 
c              adjustments of second moments
c
c FSB 07/08/99 Revised to apportion new water, ammonium, and nitrate
c              in terms of sulfate + ammonium + nitrate in each mode
c              before equilibration.
c
c FSB 10/23/01 Revised to include the heterogeneous reaction of N2O5 
c              with aerosol to produce HNO3, based on Pleim et al 
c              (1995).  AIRPRS & DT added to call vector. These 
c              modifications assume that GETPAR has been called prior 
c              to calling EQL3.  It is also assumed that AH2OI and AH2OJ 
c              have been added to the transport routines.
c
c SJR 04/24/02 Revised to use ISORROPIA in place of RPMARES for
c              thermodynamic equilibrium calculations
c 
c GLG 08/15/02 Revised to use radius instead of diameter in calculation 
c              of N2O5->HNO3 rate constant
c
c GLG 03/10/03 Revised to use composition-dependent gamma from Riemer
c              et al. (2003)
c
c SJR 03/14/03 Revised to use the effective diameter in the calculation 
c              of N2O5->HNO3 rate constant instead of the geometric
c              mean diameter
c
c SJR 04/15/03 Corrected units in the HNO3 yield from the heterogeneous 
c              N2O5 Rxn of Riemer et al. (2003)
c
c PVB 09/21/04 Adjusted MWNO3, MWHNO3, MWSO4, MWNH3, MWNH4, and MWN2O5
c              to be consist with mechanism files.  Added in-line
c              documentation.
c
c PVB 01/19/05 Added SO4RATE to the call vector.  Adjusted TSO4 to 
c              include aerosol SO4 plus newly-formed SO4 from gas-phase
c              production (i.e., SO4RATE*DT).  This is necessary to  
c              ensure that the gas and inorganic fine PM (i+j) 
c              concentrations are in thermodynamic equilibrium at the 
c              end of each time step.
c
c PVB 06/13/05 Added TNA and TCL to ISOROPIA call; put all chloride
c              in accumulation mode
c
c PVB 06/22/05 Added logic to avoid mass conservation problems in
c              ISORROPIA-ISRP3F with low NH4 and Cl concentrations
c
c GS  04/04/06 Revised to use T,RH-dependent gamma from Evans & Jacob 
c              (2005).  Retained ratio of GAMMA1/GAMMA2 from Riemer et
c              al. (2003).
c
c PVB 04/06/06 Added GAMMA to the call vector, so it can be written
c              to the aerosol diagnostic file.
c
c PVB 05/25/06 Set upper limit for the RH input to ISORROPIA (i.e., RHI) 
c              at 95%, based on guidance from Dr. Thanos Nenes of Georgia
c              Tech and Uma Shankar of Univ. of North Carolina.  Higher
c              RH levels yield extremely large AH2O outputs from ISORROPIA, 
c              indicative of cloud droplet or fog formation.
c
c JOY 04/04/07 Optimized GAMMA calculation; initialized GAMMA in case
c              RH < 1%; Note: some compilers recognize GAMMA as as
c              intrinsic GAMMA function.
c
c PVB 11/02/07 Moved calculation of GAMMA to a new subroutine, HETCHEM.
c
c PVB 11/03/07 Modified code for transferring ISOROPIA output values to
c              local variables, based on the recommendation of Fang-Yi
c              Cheng and Andrey Martynenko of Univ.of Houston.  The 
c              revised code ensures mass conservation of TNH4, TNO3, and
c              TCL, and allows users to switch between stable and 
c              metastable calculations (i.e., the CNTRL(2) flag) without
c              modifying any code that follows the ISOROPIA call.
c
c  REFERENCES:
c   1. Nenes, A., ISORROPIA v1.5 Reference Manual, 2004.
c      Available at http://nenes.eas.gatech.edu/ISORROPIA/
c
c   2. Riemer, N., H. Vogel, B. Vogel, B. Schell, I. Ackermann, C.
c      Kessler, and H. Hass, Impact of the heterogeneous hydrolysis
c      of N2O5 on chemistry of nitrate aerosol formation in the lower
c      troposphere under photosmog conditions.  J. Geophys. Res., Vol 
c      108, No D4, 4144, doi:10.1029/2002JD002436, 2003.
c
c   3. Pleim, J.E., F.S. Binkowski, J.K.S. Ching, R.L. Dennis, and N.V.
c      Gallani, 1995, An improved representation of the reaction of 
c      N2O5 on aerosols for mesoscale air quality models.  In "Regional 
c      Photochemical Measurement and Modeling Studies, Vol 2 - Results
c      and Status of Modeling," Eds A.J. Ranzieri and P.A. Solomon, pp 
c      904-913.
c
c   4. Evans, M.J. and D.J. Jacob, Impact of new laboratory studies of 
c      N2O5 hydrolysis on global model budgets of troposphreic nitrogen
c      oxides, ozone, and OH.  Geophys. Res. Lett., 32, L09813, 
c      doi:10.1029/2005GL022469

      SUBROUTINE EQL3( NSPCSDA, CBLK, SO4RATE, AIRTEMP, AIRRH, DT )

      USE AERO_INFO_AE4    ! Module containing CBLK indices and 
                           ! several constants
      IMPLICIT NONE

C *** ARGUMENTS

      INTEGER NSPCSDA        ! number of species in CBLK
      REAL CBLK( NSPCSDA  )  ! main array of variables
      REAL SO4RATE           ! Gas-phase SO4 prod rate [ug/m3/s]
      REAL AIRTEMP           ! Air temperature [ K ]
      REAL AIRRH             ! Fractional relative humidity
      REAL DT                ! Synchronization time step

C *** PARAMETERS

C *** molecular weights, hardcoded to match values in mechanism files

      REAL        MWNO3            ! molecular weight for NO3
      PARAMETER ( MWNO3  = 62.0 )

      REAL        MWHNO3           ! molecular weight for HNO3
      PARAMETER ( MWHNO3 = 63.0 )

      REAL        MWSO4            ! molecular weight for SO4
      PARAMETER ( MWSO4  = 96.0 )

      REAL        MWNH3            ! molecular weight for NH3
      PARAMETER ( MWNH3  = 17.0 )

      REAL        MWNH4            ! molecular weight for NH4
      PARAMETER ( MWNH4  = 18.0 )

      REAL        MWNA             ! molecular weight for Na
      PARAMETER ( MWNA   =  23.0 )
    
      REAL        MWCL             ! molecular weight for Cl  
      PARAMETER ( MWCL   =  35.0 )

      REAL        MWHCL            ! molecular weight for HCl 
      PARAMETER ( MWHCl  =  36.0 )

C *** conversion factors for thermodynamics vars

      REAL        FAERNH4          ! for ug  -> mole
      PARAMETER ( FAERNH4  = 1.0E-6 / MWNH4 )

      REAL        FAERNO3          ! for ug  -> mole
      PARAMETER ( FAERNO3  = 1.0E-6 / MWNO3 )

      REAL        FAERNH3          ! for ug  -> mole
      PARAMETER ( FAERNH3  = 1.0E-6 / MWNH3 )

      REAL        FAERHNO3         ! for ug  -> mole
      PARAMETER ( FAERHNO3 = 1.0E-6 / MWHNO3 )

      REAL        FAERSO4          ! for ug  -> mole
      PARAMETER ( FAERSO4  = 1.0E-6 / MWSO4 )

      REAL        FAERH2O          ! for ug  -> mole
      PARAMETER ( FAERH2O  = 1.0E-6 / MWWAT )

      REAL        FAERNA           ! for ug  -> mole
      PARAMETER ( FAERNA   = 1.0E-6 / MWNA  )

      REAL        FAERCL           ! for ug  -> mole
      PARAMETER ( FAERCL   = 1.0E-6 / MWCL  )

      REAL        FAERHCL          ! for ug  -> mole
      PARAMETER ( FAERHCL  = 1.0E-6 / MWHCL )

      REAL        DFNH3           ! for mole  -> ug
      PARAMETER ( DFNH3    = 1.0 / FAERNH3 )

      REAL        DFHNO3           ! for mole  -> ug
      PARAMETER ( DFHNO3   = 1.0 / FAERHNO3 )

      REAL        DFHCL           ! for mole  -> ug
      PARAMETER ( DFHCL    = 1.0 / FAERHCL )

      REAL        DFH2O            ! for mole  -> ug
      PARAMETER ( DFH2O    = 1.0 / FAERH2O )

      REAL        DFNH4            ! for mole  -> ug
      PARAMETER ( DFNH4    = 1.0 / FAERNH4 )

      REAL        DFNO3            ! for mole  -> ug
      PARAMETER ( DFNO3    = 1.0 / FAERNO3 )

      REAL        DFCL            ! for mole  -> ug
      PARAMETER ( DFCL     = 1.0 / FAERCL  )

C *** floor value for concentrations

      REAL      CONMIN ! concentration lower limit [ ug/m**3 ]
      PARAMETER ( CONMIN = 1.0E-30 )


C *** LOCAL VARIABLES

C *** chemical species concentrations [ug/m3]

      REAL      ASO4     ! i+j mode sulfate
      REAL      ANO3     ! i+j mode nitrate
      REAL      ANH4     ! i+j mode ammonium
      REAL      AH2O     ! i+j mode water
      REAL      GNH3     ! gas-phase ammonium
      REAL      GNO3     ! gas-phase nitric acid
      REAL      ANA      ! i+j mode sodium
      REAL      ACL      ! i+j mode chloride
      REAL      GCL      ! gas-phase hydrochloric acid

C *** variables to calculate mass fractions of SO4 in each mode

      REAL      SMASS_I, SMASS_J
      REAL      FRACI, FRACJ

C *** 2nd and 3rd moments before equilibration (without H2O)

      REAL      OLD_M3_I, OLD_M3_J
      REAL      OLD_M2_I, OLD_M2_J

C *** chemical concentrations for ISORROPIA [mol/m3]

      REAL      TSO4     ! i+j mode SO4 + newly-formed gas-phase SO4
      REAL      TNH4     ! i+j mode NH4 + gaseous NH3
      REAL      TNO3     ! i+j mode NO3 + gaseous HNO3
      REAL      TNA      ! i+j mode Na
      REAL      TCL      ! i+j mode Cl  + gaseous HCl

C *** ISORROPIA input variables

      REAL*8 WI(5)                   ! species array
      REAL*8 RHI                     ! relative humidity
      REAL*8 TEMPI                   ! temperature
      REAL*8 CNTRL(2)                ! control parameters 
      
C *** ISORROPIA output variables
      
      REAL*8 WT(5)                   ! species output array
      REAL*8 GAS(3)                  ! gas-phase   "     " 
      REAL*8 AERLIQ(12)              ! liq aerosol "     " 
      REAL*8 AERSLD(9)               ! solid "     "     " 
      CHARACTER*15 SCASI             ! subcase number output
      REAL*8 OTHER(6)                ! supplmentary output array

C *** Variables to account for mass conservation violations in ISRP3F

      LOGICAL   TRUSTNH4             ! false if ISOROPIA's partitioning
                                     !  of NH4/NH3 is to be ignored
      LOGICAL   TRUSTCL              ! false if ISOROPIA's partitioning
                                     !  of Cl/HCl is to be ignored

C *** 2nd and 3rd moments after equilibration (with H2O)

      REAL      NEW_M3_I, NEW_M3_J
      REAL      NEW_M2_I, NEW_M2_J

C.......................................................................


C *** For extremely low relative humidity ( less than 1% ) set the
c     water content to a minimum value and skip the calculation

      IF ( AIRRH .LT. 0.01 ) THEN
         CBLK( VH2OAI ) = CONMIN
         CBLK( VH2OAJ ) = CONMIN
         RETURN
      END IF

C *** calculate total aerosol ammonium, nitrate, and sulfate [ug/m3]

      ANA  = CBLK( VNAJ )   + CBLK( VNAI )
      ANH4 = CBLK( VNH4AJ ) + CBLK( VNH4AI )
      ANO3 = CBLK( VNO3AJ ) + CBLK( VNO3AI )
      ASO4 = CBLK( VSO4AJ ) + CBLK( VSO4AI )
      ACL  = CBLK( VCLJ )   + CBLK( VCLI )

c *** fetch vapor-phase concentrations [ug/m3]

      GNO3  = CBLK( VHNO3 )
      GNH3  = CBLK( VNH3 )
      GCL   = CBLK( VHCL )

c *** calculate the fractions of aerosol SO4 in the i-mode and j-mode

      SMASS_I = CBLK( VSO4AI )
      SMASS_J = CBLK( VSO4AJ )
      FRACI = SMASS_I / ( SMASS_I + SMASS_J )
      FRACJ = MAX( 0.0, 1.0 - FRACI )

C *** capture values of "dry" 2nd and 3rd moments before equilibration
c     the folowing code assumes that GETPAR has been called with
c     M3_WET_FLAG set to .FALSE. and that the 2nd and 3rd moments have
c     been adjusted for the new SOA.

      OLD_M3_I = CBLK( VAT3 )
      OLD_M3_J = CBLK( VAC3 )
      OLD_M2_I = CBLK( VAT2 )
      OLD_M2_J = CBLK( VAC2 )


C *** GAS/PARTICLE THERMODYNAMIC EQUILIBRIUM
c      This section of code calculates the distribution of ammonia/
c    ammonium, nitric acid/nitrate, and water between the gas and 
c    aerosol phases as a function of total sulfate, total ammonia, 
c    total nitrate, relative humidity, and temperature.  It is assumed
c    that the aerosol is entirely aqueous (i.e., metastable assumption),
c    irrespective of ambient relative humidity.

C...gas+aerosol species concentrations [mol/m3]

      TNA  = ANA  * FAERNA 
      TSO4 = (ASO4 + SO4RATE * DT) * FAERSO4
      TNH4 = GNH3 * FAERNH3 + ANH4 * FAERNH4 
      TNO3 = ANO3 * FAERNO3 + GNO3 * FAERHNO3 
      TCL  = ACL  * FAERCL  + GCL  * FAERHCL  

C...double precision vars for ISORROPIA

      WI( 1 ) = TNA 
      WI( 2 ) = TSO4
      WI( 3 ) = TNH4
      WI( 4 ) = TNO3
      WI( 5 ) = TCL 

      RHI   = MIN( 0.95, AIRRH )
      TEMPI = AIRTEMP

      CNTRL( 1 ) = 0.D0               ! forward problem
      CNTRL( 2 ) = 1.D0               ! aerosol in metastable state

C...set flags to account for mass conservation violations in ISRP3F

      TRUSTCL  = .TRUE.
      TRUSTNH4 = .TRUE.
      IF ( TNA+TCL .LT. 1.E-20 ) THEN
         TRUSTCL = .FALSE.
      ELSE
         IF ( TNH4 .LT. 1.E-10 ) TRUSTNH4 = .FALSE.
         IF ( TCL  .LT. 1.E-10 ) TRUSTCL  = .FALSE.
      ENDIF

C...calculate gas/particle equilibrium

      CALL ISOROPIA ( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,
     &                AERSLD, SCASI, OTHER )

C...convert output back to single precision and [ug/m3] units
c    Note: all ISOROPIA outputs are in [mol/m3] units
c     AERLIQ(03) - NH4+(aq)
c     AERLIQ(04) - Cl-(aq)
c     AERLIQ(07) - NO3-(aq)
c     AERLIQ(08) - H2O
c     AERLIQ(09) - NH3(aq) (undissociated)
c     AERLIQ(10) - HCl(aq) (undissociated)
c     AERLIQ(11) - HNO3(aq) (undissociated)
c     GAS(1) - NH3
c     GAS(2) - HNO3
c     GAS(3) - HCl
c     AERSLD = 0.D0 because of metastable assumption

      AH2O = AERLIQ( 8 ) * DFH2O
      ANO3 = ( WT( 4 ) - GAS( 2 ) ) * DFNO3
      GNO3 = GAS( 2 ) * DFHNO3
      IF ( TRUSTNH4 ) THEN
         ANH4 = ( WT( 3 ) - GAS( 1 ) ) * DFNH4
         GNH3 = GAS( 1 ) * DFNH3
      ENDIF
      IF ( TRUSTCL ) THEN
         ACL  = ( WT( 5 ) - GAS( 3 ) ) * DFCL
         GCL  = GAS( 3 ) * DFHCL
      ENDIF


C *** UPDATE CBLK ARRAY
c      Assume that NH4, NO3, and H2O, exhibit the same size 
c     distribution as SO4.  In other words, these species partition 
c     to the Aitken and accumulation modes in proportion to the
c     pre-existing SO4 mass.  Ensure that all species remain above
c     the minimum concentration.  Update 3rd moments to account for
c     changes in the NH4, NO3, and H2O concentrations.  Finally, 
c     adjust 2nd moments under the assumption that the condensation 
c     and evaporation of NH4, NO3, and H2O, do not affect the 
c     geometric standard deviation of either mode.

c *** update Aitken mode

      CBLK( VH2OAI ) = MAX( FRACI * AH2O, CONMIN )
      CBLK( VNH4AI ) = MAX( FRACI * ANH4, CONMIN )
      CBLK( VNO3AI ) = MAX( FRACI * ANO3, CONMIN )
      CBLK( VCLI )   = CONMIN

C *** update accumulation mode:

      CBLK( VH2OAJ ) = MAX( FRACJ * AH2O, CONMIN )
      CBLK( VNH4AJ ) = MAX( FRACJ * ANH4, CONMIN )
      CBLK( VNO3AJ ) = MAX( FRACJ * ANO3, CONMIN )
      CBLK( VCLJ )   = MAX( ACL - CONMIN, CONMIN )

C *** update gas / vapor phase

      CBLK( VNH3  ) = MAX( GNH3,  CONMIN )
      CBLK( VHNO3 ) = MAX( GNO3,  CONMIN )
      CBLK( VHCL  ) = MAX( GCL,   CONMIN )

c *** update third moments
c      Note: new values of 3rd moment are not transferred to CBLK 
c      because GETPAR is called immediately after EQL3 and the
c      CBLK update occurs there.

      NEW_M3_I = MAX( CONMIN, ( SO4FAC  * CBLK( VSO4AI  ) +
     &                          NH4FAC  * CBLK( VNH4AI  ) +
     &                          H2OFAC  * CBLK( VH2OAI  ) +
     &                          NO3FAC  * CBLK( VNO3AI  ) +
     &                          ORGFAC  * CBLK( VORGAI  ) +
     &                          ORGFAC  * CBLK( VORGPAI ) +
     &                          ORGFAC  * CBLK( VORGBAI ) +
     &                          ANTHFAC * CBLK( VP25AI  ) +
     &                          ANTHFAC * CBLK( VECI    ) +
     &                          SEASFAC * CBLK( VNAI    ) +
     &                          SEASFAC * CBLK( VCLI    ) ) )

      NEW_M3_J = MAX( CONMIN, ( SO4FAC  * CBLK( VSO4AJ  ) +
     &                          NH4FAC  * CBLK( VNH4AJ  ) +
     &                          H2OFAC  * CBLK( VH2OAJ  ) +
     &                          NO3FAC  * CBLK( VNO3AJ  ) +
     &                          ORGFAC  * CBLK( VORGAJ  ) +
     &                          ORGFAC  * CBLK( VORGPAJ ) +
     &                          ORGFAC  * CBLK( VORGBAJ ) +
     &                          ANTHFAC * CBLK( VP25AJ  ) +
     &                          ANTHFAC * CBLK( VECJ    ) +
     &                          SEASFAC * CBLK( VNAJ    ) +
     &                          SEASFAC * CBLK( VCLJ    ) ) )

C *** adjust second moments assuming that the geometric standard
c     deviations are unaffected by volatile inorganic species

      NEW_M2_I = OLD_M2_I * ( NEW_M3_I / OLD_M3_I ) ** TWO3
      CBLK(   VAT2 ) = NEW_M2_I

      NEW_M2_J = OLD_M2_J * ( NEW_M3_J / OLD_M3_J ) ** TWO3
      CBLK(   VAC2 ) = NEW_M2_J
      

      RETURN

      END SUBROUTINE EQL3


c /////////////////////////////////////////////////////////////////////
c  SUBROUTINE NEWPART3 calculates the new particle production rate
c    due to binary nucleation of H2SO4 and H2O vapor.  The nucleation
c    rate is a parameterized function of temperature, relative humidity,
c    and the vapor-phase H2SO4 concentration, following the work of
c    Kulmala et al (1998).  This rate is then used to determine the 
c    production rates of aerosol number, 2nd moment, and aerosol
c    sulfate, following the description in Section 1.2 of Binkowski & 
c    Roselle (2003), except the new particles are assumed to be of 2.0
c    nm diameter instead of 3.5 nm.
c
c  KEY SUBROUTINES/FUNCTIONS CALLED: none
c
c  REVISION HISTORY:
c
c FSB 11/29/99 Extracted code from PARTPROD_VA subroutine of RPM
c
c SJR ??/??/?? New call vector and limited the mass production rate
c
c FSB 04/11/02 Decreased the diameter of new particles from 3.5 nm to
c     2 nm.  New particles are monodispersed.
c
c PVB 09/21/04 Changed MWH2SO4 from 98.07354 to 98.0 g/mol for 
c     consistency with the mechanism files.  Added in-line documentation
c     with input from FSB.
c
c  REFERENCES:
c   1. Kulmala, M., A. Laaksonen, and L. Pirjola, Parameterizations for
c      sulfuric acid/water nucleation rates.  J. Geophys. Res., Vol 103, 
c      No D7, 8301-8307, 1998.
c
c   2. Binkowski, F.S. and S.J. Roselle, Models-3 Community 
c      Multiscale Air Quality (CMAQ) model aerosol component 1:
c      Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
c      doi:10.1029/2001JD001409, 2003.

      SUBROUTINE NEWPART3 ( RH,Temp,XH2SO4,SO4RATE,
     &                      DNDT,DMDT_so4, DM2DT )

      USE AERO_INFO_AE4    ! Module containing useful constants
                           ! including RGASUNIV, AVO, and RHOSO4
      implicit none

C *** INPUT ARGUMENTS

      real RH         ! fractional relative humidity
      real Temp       ! ambient temperature [ K ]
      real XH2SO4     ! sulfuric acid concentration [ ug / m**3 ]
      real SO4RATE    ! gas-phase H2SO4 production rate [ ug / m**3 s ]

c *** OUTPUT ARGUMENTS

      real*8 DNDT     ! particle number production rate [ m^-3 s^-1 ]
      real*8 DMDT_so4 ! SO4 mass production rate [ ug / m**3 s ]
      real*8 DM2DT    ! second moment production rate [ m**2 / m**3 s ]

c *** PARAMETERS

      character*16 PNAME
       parameter( PNAME = 'NEWPART         ')

c *** conversion factors:

      real MWh2so4   ! molar mass of h2so4 [ g/mole ]
       parameter( MWh2so4 = 98.0 )

      real ugm3_ncm3 ! micrograms/m**3 to number/cm**3
       parameter( ugm3_ncm3 = (AVO / MWh2so4) * 1.0e-12 )

c *** particle size parameters:

      real d20   ! diameter of a new particle [cm]
       parameter( d20 = 2.0E-07)

      real d20sq ! new-particle diameter squared [cm**2]
       parameter( d20sq = d20 * d20 )

      real m2_20 ! new-particle diameter squared [m**2]
       parameter( m2_20 = 1.0e-4 * d20sq  )

      real v20   ! volume of a new particle [cm**3]
       parameter( v20 = PI6 * d20 * d20sq  )

      real*8 sulfmass  ! mass of a new particle [ug]
       parameter(sulfmass = 1.0e3 * RHOSO4 * v20 )

      real*8 sulfmass1 ! inverse of sulfmass [ug**-1]
       parameter(sulfmass1 = 1.0D0 / sulfmass )


c *** LOCAL VARIABLES used to determine nucleation rate

      real Nwv     ! water vapor concentration [ 1/cm**3 ]
      real Nav0    ! saturation vapor conc of H2SO4 [ 1/cm**3 ]
      real Nav     ! H2SO4 vapor concentration [ 1/cm**3 ]
      real Nac     ! critical H2SO4 vapor concentration [ 1/cm**3 ]
      real RA      ! fractional relative acidity
      real delta   ! temperature-correction term
      real Xal     ! H2SO4 mole fraction in the critical nucleus
      real Nsulf   ! see usage
      real*8 chi   ! factor to calculate Jnuc
      real*8 JnucK ! nucleation rate [ cm ** -3  s ** -1 ]
                   ! (Kulmala et al.)
      real tt      ! dummy variable for statement functions

c *** STATEMENT FUNCTIONS

      real ph2so4, ph2o ! arithmetic statement functions for saturation
                        ! vapor pressures of h2so4 and h2o [Pa] taken 
                        ! from Appendix of Kulmala et al. (1998), p8306

      ph2o(tt) = exp(77.34491296 -7235.424651 / tt
     &               - 8.2 * log (tt) + tt * 5.7113e-03 )

      ph2so4(tt) = exp(27.78492066 - 10156.0/tt)

C ---------------------------- Begin code ------------------------------

C *** initialize variables
      DNDT     = 0.0D0
      DMDT_so4 = 0.0D0
      DM2DT    = 0.0D0

C *** NUCLEATION RATE
c      Calculate the sulfuric acid/water nucleation rate.  This code 
c     implements Section 3 of Kulmala et al. (1998).  Note that all
c     variables in the Kulmala parameterization are in cgs units.

c *** calculate water vapor concentration [1/cm**3] using ambient RH
c      and the formula in Appendix of Kulmala et al. (1998), p8306.

      Nwv = RH * ph2o(Temp) / ( RGASUNIV * Temp ) * AVO * 1.0 e-6

c *** calculate saturation vapor concentration of H2SO4 [1/cm**3]
c       using formula in the Appendix of Kulmala et al. (1998), p8306.

      Nav0 = ph2so4(Temp) / ( RGASUNIV * Temp ) * AVO * 1.0 e-6

c *** convert ambient H2SO4 vapor concentration into [1/cm**3] units

      Nav = ugm3_ncm3 * XH2SO4

c *** calculate critical concentration of H2SO4 vapor needed to produce
c      1 particle/(cm**3 s) using Equation 18 of Kulmala et al (1998)

      Nac = exp( - 14.5125 + 0.1335 * Temp
     &           - 10.5462 * RH + 1958.4 * RH  / Temp )

c *** calculate relative acidity, defined as the ambient concentration
c      divided by the saturation concentration of H2SO4 vapor

      RA = Nav / Nav0

c *** calculate temperature correction factor using Equation 22 of
c       Kulmala et al (1998)

      delta = 1.0 + ( Temp - 273.15) / 273.15

c *** calculate mole fraction of H2SO4 in the critical nucleus using
c      Equation 17 of Kulmala et al (1998)

      Xal = 1.2233 - 0.0154 * RA / ( RA + RH ) + 0.0102 * log (Nav)
     &    - 0.0415 * log( Nwv) + 0.0016 * Temp

C *** calculate Nsulf as defined in Equation 21 of Kulmala et al (1998)

      Nsulf = log( Nav / Nac)

c *** calculate natural log of nucleation rate using Equation 20 of
c      Kulmala et al (1998)

      chi = 25.1289 * Nsulf - 4890.8 * Nsulf / Temp
     &    - 1743.3 / Temp - 2.2479 * delta * Nsulf * RH
     &    + 7643.4 * Xal / Temp
     &    - 1.9712 * Xal * delta / RH

c *** calculate nucleation rate using Eq 19 of Kulmala et al (1998)

      JnucK = exp(chi) ! [ # / cm**3 s ]


C *** MOMENT PRODUCTION RATES
c      Calculate the production rate of number, sulfate mass, and 2nd
c     moment, due to the H2SO4/H2O nucleation assuming that each critical
c     nucleus grows instantaneously to 2.0 nm.  This code follows Section
c     1.2 of Binkowski & Roselle (2003), except the assumed particle
c     diameter has been changed from 3.5 to 2.0 nm.

C *** convert production rate of number to [ # / m**3 s]

      DNDT = (1.0e06) * JnucK

C *** calculate mass production rate [ ug / (m**3 s) ] analogous to
c     Equation 6a of Binkowski & Roselle (2003).  Set the upper limit 
c     of the mass production rate as the gas-phase production rate of 
c     H2SO4, and adjust the number production rate accordingly.

      DMDT_so4 = sulfmass * DNDT
      IF ( DMDT_so4 .GT. SO4RATE ) THEN
        DMDT_so4 = SO4RATE
        DNDT = DMDT_SO4 * sulfmass1
      END IF

C *** calculate the production rate of 2nd moment [ m**2 / (m**3 s) ]
c     This is similar to Equation 6c of Binkowski & Roselle (2003),
c     except the factor of PI is removed and the assumed particle 
c     diameter is different.

      DM2DT =  DNDT * m2_20

      return
      
      end subroutine NEWPART3


c //////////////////////////////////////////////////////////////////
c  SUBROUTINE ORGAER3 implements Dr. Benedikt Schell's method of
c   calculating secondary organic aerosol, which is based on the 
c   absorptive partitioning model of Pankow that was extended by 
c   Odum,1996.
c
c   Input includes the concentrations of ROG's that reacted
c   (ORGPROD) during the time step, the total gas + aerosol
c   concentration of each semi-volatile organic compound (SVOCS),
c   and the concentrations of primary and secondary organic
c   aerosols at the beginning of the time step (OLDPOA, OLDSOA_A,
c   and OLDSOA_B).  The subroutine partitions the semi-volatile
c   species into the gas and aerosol phases and returns updated
c   values of the total semi-volatile organic concentrations
c   (SVOCS), and the anthropogenic and biogenic SOA (SOA_A, SOA_B).
c
c   Saturation vapor pressures (cstar) and mass-based stoichiometric
c   yield coefficients (alpha) are obtained from either smog chamber
c   experiments or from published estimates in cases where smog
c   chamber data are unavailable (e.g. alkanes, cresol).  Saturation
c   vapor pressures are modified as a function of temperature using
c   eqn 6 of Sheehan & Bowman.
c  
c   If the pre-existing organic aerosol concentration is zero,
c   gas/aerosol equilibrium is established only after the organic gas 
c   concentration reaches the threshold value defined in eqn 9 of 
c   Schell et al.  Until this threshold value is reached, organic
c   vapors do not partition to the aerosol phase.
c
c   Once the organic gas/aerosol equilibrium has been established,
c   gas and aerosol-phase concentrations of each condensible species
c   are calculated iteratively using a globally convergent variation
c   of Newton's method (SUBROUTINE NEWT), as described in eqn 8 of
c   Schell et al.
c
c  KEY SUBROUTINES/FUNCTIONS CALLED:  NEWT
c
c  REVISION HISTORY:
c     coded August 1, 2001 by Dr. Francis S. Binkowski
c
c     Revised April 4, 2003 by Gerald Gipson to allow for evaporation
c     of organics from aerosols. Now total vapor + aerosol phase is
c     repartitioned at each time step and totorgsw ( Mo ) does not
c     include oldsoa.
c
c     Revised July 14, 2003 by Dr. Prakash V. Bhave
c     - changed cstar(2,3) from 10.103 & 90.925 to 111.11 & 1000.0
c       because smog chamber data of Kalberer et al. were collected
c       at 298 K (not 310 K, as was previously assumed)
c     - changed mw_vap(9,10) from 184 g/mol to 177 g/mol, to be
c       consistent with mwsoa_b
c     - modified threshold criteria for establishing gas/particle 
c       equilibrium by removing the loose criterion involving "mtot"
c     - changed variable names to reflect that the combined vapor +
c       aerosol concentrations are now being repartitioned during
c       each time step (not just the newly formed SVOC's)
c     - added documentation and removed extraneous lines of code
c
c     Revised December 4, 2003 by Dr. Francis S. Binkowski
c     - output variables ORGRATE and ORGBRATE removed and replaced
c       by SOA_A and SOA_B, the newly equilibrated values of
c       Anthropogenic and Biogenic SOA, respectively.  These are non-
c       negative values.
c     - variable jj also removed
c
c     Revised January 8, 2004 by Dr. Prakash V. Bhave
c     - removed the output variable YIELD.  It has no physical meaning
c       after the 12/04/2003 revisions.
c
c     Revised January 12, 2004 by Dr. Chris G. Nolte
c     - for computational efficiency, modified the initial caer guess
c       used as input to NEWT.  If NEWT returns check .eq. true, then
c       NEWT is called again with a guess of caer = 0.5*ctot
c     - removed ITS parameter from NEWT call vector
c     - fixed bug where concentrations less than TOLMIN (i.e., 1.0E-12)
c       were reset to 1.0e-30
c     - removed extraneous code related to "Pandis method" of SVOC
c       partitioning when threshold criterion is not met (i.e., 
c       insufficient organic matter to establish gas/particle 
c       equilibrium)  ** results unaffected by this change

c  REFERENCES:
c   1. Schell, B., I. J. Ackermann, H. Hass, F. S. Binkowski, and
c      A. Abel, Modeling the formation of secondary organic aerosol
c      within a comprehensive air quality modeling system, J. Geophys.
c      Res., Vol 106, No D22, 28275-28293, 2001.
c
c   2. Pankow, J. F., An absorption model of gas/particle partitioning
c      of organic compounds in the atmosphere, Atmos. Environ., Vol 28, 
c      No 2, 185-188, 1994.
c
c   3. Odum, J. R., T. Hoffmann, F. Bowman, D. Collins, R. C. Flagan,
c      and J. H. Seinfeld, Gas/particle partitioning and secondary
c      organic aerosol yields, Environ. Sci. Technol., Vol 30, No 8, 
c      2580-2585, 1996.
c
c   4. Sheehan, P. E. and F. M. Bowman, Estimated effects of temperature
c      on secondary organic aerosol concentrations, Environ. Sci.
c      Technol., Vol 35, No 11, 2129-2135, 2001.
c
c   5. Strader, R., F. Lurmann, and S. N. Pandis, Evaluation of  
c      secondary organic aerosol formation in winter, Atmos. Environ.,
c      Vol 33, 4849-4863, 1999.
c
c   6. Kalberer, M., J. Yu, D. R. Cocker, R. C. Flagan, and J. H.
c      Seinfeld, Aerosol formation in the cyclohexene-ozone system,
c      Environ. Sci. Technol., Vol 34, No 23, 4894-4901, 2000.
c
c   7. Odum, J. R., T. P. W. Jungkamp, R. J. Griffin, R. C. Flagan, 
c      and J. H. Seinfeld, The atmospheric aerosol-forming potential
c      of whole gasoline vapor, Science, Vol 276, 96-99, April 4, 1997.
c
c   8. Griffin, R. J., D. R. Cocker III, R. C. Flagan, and J. H. 
c      Seinfeld, Organic aerosol formation from the oxidation of 
c      biogenic hydrocarbons, J. Geophys. Res., Vol 104, No D3, 
c      3555-3567, 1999.
c
c   9. Bian, F. and F. M. Bowman, Theoretical method for lumping 
c      multicomponent secondary organic aerosol mixtures, Environ.
c      Sci. Technol., Vol 36, No 11, 2491-2497, 2002.
c
c  10. Tao, Y. and P. H. McMurry, Vapor pressures and surface free
c      energies of C14-C18 monocarboxylic acids and C5 and C6 
c      dicarboxylic acids, Environ. Sci. Technol., Vol 23, 1519-1523,
c      1989.

       subroutine ORGAER3(DT, LAYER, AIRTEMP, AIRPRS,
     &                    ORGPROD, NPSPCS, SVOCS, NCVAP,
     &                    OLDSOA_A, OLDSOA_B, OLDPOA,
     &                    SOA_A, SOA_B, SUMPROD, LOGDEV )

       use soa_NEWT      ! module containing all information and 
                         ! subroutines used in SUBROUTINE NEWT
       implicit none


c     Input variables

       real DT           ! timestep [ s ]
       integer  LAYER    ! model layer number
       real AIRTEMP      ! air temperature [ K ]
       real AIRPRS       ! air pressure [ Pa ]

       integer NPSPCS    ! number of aerosol precursor species
       real ORGPROD(NPSPCS) ! change in precursor mixing ratio [ ppm ]

       integer NCVAP     ! number of SVOC species that partition
       real SVOCS(NCVAP) ! semi-volatile conc [ ppm ] (also output)

       real OLDSOA_A     ! input SOA of anthropogenic origin
                         !   [ ug / m**3 ]
       real OLDSOA_B     ! input SOA of biogenic origin
                         !   [ ug / m**3 ]
       real OLDPOA       ! input primary organic aerosol
                         !   [ ug / m**3 ]

c     Output variables

       real SOA_A        ! new value of anthropogenic
                         !   SOA [ ug / ( m**3 ) ]
       real SOA_B        ! new value of biogenic
                         !   SOA [ ug / ( m**3 ) ]
       real SUMPROD      ! total conc of new SOA produced during
                         !   time step [ ug / m**3 ]

c     Internal variables
c *** NOTE! dimensions of internal arrays are set to constants,
c     nprec and NP, which are defined in the soa_NEWT module

       integer LOGDEV    ! unit number for log file

       logical first_time
        data first_time / .true. /
        save first_time

       integer ii    ! loop index

c  Species order within the ORGPROD and mw_rog arrays:
c        (1)  "long" alkanes ( alk5 in saprc99)
c        (2)  internal alkenes ( cyclohexene. ole2 in saprc99 )
c        (3)  aromatics like xylene  (aro2 in saprc99)
c        (4)  aromatics like cresol (cres in saprc99)
c        (5)  aromatics like toluene (aro1 in saprc99)
c        (6)  monoterpenes (trp1 in saprc99)
c
c  Species order within SVOCS, alpha, cstar, and mw_vap arrays:
c        (1) alkane
c        (2) alkene_1
c        (3) alkene_2
c        (4) xylene_1
c        (5) xylene_2
c        (6) cresol
c        (7) toluene_1
c        (8) toluene_2
c        (9) monoterpene_1
c       (10) monoterpene_2

       integer nbio1, nbio2      ! indices for biogenic SOA
        parameter( nbio1 = 9, nbio2 = 10 )

       integer nole1, nole2      ! indices for olefin-derived SOA
        parameter( nole1 = 2, nole2 = 3 )

c  Gas/aerosol partitioning parameters
c
c   The effective saturation concentrations (cstar) of the following
c   species are approximated as the reciprocal of the partition
c   coefficients (Kom) determined from smog chamber experiments (with
c   representative smog chamber temperatures given in parentheses)
c      alkene    cyclohexene experiments of Kalberer et al. (298 K)
c      xylene    low-yield aromatic curve of Odum,1997 (310 K)
c      toluene   high-yield aromatic curve of Odum,1997 (310 K)
c   then adjusted to 298 K using the method of Sheehan and Bowman.
c   Mass-based stoichiometry coefficients (alpha) for the above species
c   are taken directly from the referenced literature.
c
c   The cstar values for alk5 and cresol are obtained from estimates
c   in Strader et al. for the AAR4 and CRES species, respectively.  
c   Strader's values are tabulated at 281.5 K, so they too are adjusted
c   to 298 K using the Sheehan and Bowman method.  The alpha values for
c   alkane and cresol are taken from Strader et al. and converted to
c   mass-based stoichiometry coefficients by assuming molecular weights of
c   114 and 108 g/mol for higher alkanes and cresol, respectively.
c
c   The monoterpene partitioning parameters are obtained from a weighted
c   average of 5 sets of smog chamber experiments reported by Griffin
c   et al.  The emissions-based weighting factors and Griffin's smog
c   chamber-based partitioning parameters are:
c      Compound     wght  alpha1   Kom1   alpha2   Kom2
c     ----------    ----  ------  ------  ------  ------
c      a-pinene     0.4    .038    .171    .326   .0040
c      b-pinene     0.25   .13     .044    .406   .0049
c      d3-carene    0.15   .054    .043    .517   .0042
c      sabinene     0.1    .067    .258    .399   .0038
c      limonene     0.1    .239    .055    .363   .0053
c   A smog chamber temperature of 313 K is assumed for all monoterpene
c   experiments and the Kom values are adjusted to 298 K accordingly.
c   The monoterpene parameter averaging follows eqn's 8 and 10 of Bian
c   and Bowman, except the alpha values in those equations are first
c   multiplied by the appropriate weighting factor listed above, to
c   account for the relative abundances of each monoterpene.
c
c   The enthalpy of vaporization (used for the temperature adjustments
c   to cstar) is assumed to be 156 kJ/mol for all SVOC's, based on the
c   data of Tao and McMurry.

       real alpha(NP)  ! mass-based stoichiometric coefficients
                       ! [ ug / m **3 ] / [ ug / m**3 ]
        data alpha / 0.0718, 0.36,  0.32,  0.038,  0.167,
     &               0.05,   0.071, 0.138, 0.0864, 0.3857 /

       real  cstar(NP) ! effective saturation concentrations
                       ! of SVOC's [ ug / m**3 ] at 298 K
        ! data cstar / 0.3103, 10.103, 90.925, 2.165, 64.946,   ! 2002
        !              0.2611,  1.716, 47.855, 0.865, 11.804 /  ! release
        data cstar / 0.3103, 111.11, 1000.0, 2.165, 64.946,
     &               0.2611,  1.716, 47.855, 0.865, 11.804 /

       real h_vap(NP) ! enthalpy of vaporization [ J / mol ]
        data h_vap / 10 * 156.0e3 /

       real rgas1   ! reciprocal of universal gas constant
        parameter( rgas1 = 1.0 /  8.314510 )

       real hfac  ! fixed value of h_vap * rgas1
        parameter( hfac = 156.0e3 * rgas1 )

c  Molecular weights
c   ROG molecular weights
c    (1) ALK5 - MW of HC8AER species in RADM2 (114 g/mol)
c    (2) OLE2 - MW of cyclohexene (82 g/mol)
c    (3) ARO2 - MW of xylene (106 g/mol)
c    (4) CRES - MW of cresol (108 g/mol)
c    (5) ARO1 - MW of toluene (92 g/mol)
c    (6) TRP1 - MW of alpha-pinene (136 g/mol)
c   SVOC & SOA molecular weights
c     150 g/mol for anthropogenic reaction products is somewhat arbitrary
c     177 g/mol for all biogenic products is based on average of
c     pinionaldehyde (168 g/mol) and pinonic acid (186 g/mol)
c   POA molecular weight is set to 220 g/mol, similar to a C15 compound

       real mw_rog(nprec) ! Molecular weight of reactive organic gases
                          ! that are SOA precursors [ g / mol ]
        data mw_rog / 114.0, 82.0, 106.0, 108.0, 92.0, 136.0 /

       real mw_vap(NP)    ! Molecular weights of SVOCs [ g / mol ]
        !data mw_vap / 8 * 150.0, 2 * 184.0  /     ! 2002 release
        data mw_vap / 8 * 150.0, 2 * 177.0  /

       real mw_vap_m1(NP) ! Inverse SVOC MW's [ mol / g ]
        save mw_vap_m1

       real mwpoa         ! Molecular weight of pre-exisiting POA [ g / mol ]
        parameter( mwpoa = 220.0 )

       real mwsoa_a       ! Molecular weight of pre-existing anthropogenic
                          ! SOA [ g / mol ]
        parameter( mwsoa_a = 150.0 )

       real mwsoa_b       ! Molecular weight of pre-existing biogenic
                          ! SOA [ g / mol ]
        parameter( mwsoa_b = 177.0 )

c     Reference temperature and pressure

       real tref          ! reference temperature
        parameter( tref = 298.0 ) ! [ K ]
       real trefm1        ! inverse of reference temperature
        parameter( trefm1 = 1.0 / tref )

       real pref          ! reference pressure [ Pa ]
        parameter( pref = 101325.0 )
       real prefm1        ! inverse of reference pressure
        parameter( prefm1 = 1.0 / pref )

c     Unit conversion factors at reference temperature and pressure

       real vap_ppm2ug(NP) ! [ ppm per ug/m3 ] for SVOC's
       real vap_ug2ppm(NP) ! [ ug/m3 per ppm ] for SVOC's
        save vap_ppm2ug, vap_ug2ppm

       real drog_ppm2ug(nprec) ! [ ppm per ug/m3 ] for ROG's
        save drog_ppm2ug

c     Variables used to adjust cstar as a function of temperature

       real convfac_298    ! P/ RT at 1 atm and 298 K [ mole / m ** 3 ]
        parameter( convfac_298 = pref * rgas1 * trefm1 )
       real convfac, convfacm1

       real tt1,tt2       ! temperature-related factors
       real tempcorr      ! temperature correction factor for cstar

c     Variables used in organic equilibrium calculations

       real drog(nprec)   ! change in precursor conc [ ug / m**3 ]
       real totrog(NP)    ! drog conc mapped to each SVOC [ ug / m**3 ]
       real c0(NP)        ! cstar at AIRTEMP [ ug / m**3 ]
       real ctoti(NP)     ! SVOC conc before the current time step
                          !  [ ug / m**3 ]
       real prod(NP)      ! SVOC generated during current time step
                          !  [ ug / m**3 ]
       real ctotf(NP)     ! SVOC conc after the current time step
                          !  [ ug / m**3 ]
       real caer(NP)      ! SVOC conc in aerosol phase [ ug / m**3 ]
       real totsoa        ! total SOA conc after time step [ ug / m**3 ]
       real totsoa_a      ! anthropogenic SOA conc after time step
       real totsoa_b      ! biogenic SOA conc after time step
       real totorgsw      ! molar concentration of POA [ u-mole / m**3 ]
       real totorg        ! SOA + POA before time step [ u-mole / m**3 ]
       real threshold     ! criterion for establishing Gas/Part equil.
       real threshmin     ! small positive number
        parameter( threshmin = 1.0 e-19 )
       real conmin        ! concentration lower limit
        parameter ( conmin = 1.0e-30 )

      real CTOLMIN
      parameter ( CTOLMIN = 1.E-06 )

c     Variable internal to NEWT subroutine
       logical check      ! flag to indicate if NEWT subroutine
                          !  converged to a spurious root

c ---------------------------------------------------------------------

       if ( first_time )  then

C *** Set unit conversion and inverse mw constants

         do ii = 1, NPSPCS
            drog_ppm2ug(ii) =  mw_rog(ii) * convfac_298
         end do

         do ii = 1, NCVAP
            vap_ppm2ug(ii) =  mw_vap(ii) * convfac_298
            vap_ug2ppm(ii) = 1.0 / vap_ppm2ug(ii)
            mw_vap_m1(ii) = 1.0 / mw_vap(ii)
         end do

         first_time = .false.

       end if ! first_time

c ---------------------------------------------------------------------

c *** Begin solution code

C *** set temperature factors

       tt1 = tref / AIRTEMP
       tt2 = trefm1 - 1.0 / AIRTEMP
       convfac = tt1 * AIRPRS * prefm1

c *** initialize drog from ORGPROD and change units to [ ug / m**3 ]
c     calculate SUMPROD

       SUMPROD = 0.0
       do ii = 1, NPSPCS
         drog(ii) = ORGPROD(ii) * drog_ppm2ug(ii) * convfac
         SUMPROD  = SUMPROD + drog(ii)
       end do

C *** explicit assignment of drog to totrog

       totrog( 1) = drog(1) ! alkane
       totrog( 2) = drog(2) ! alkene_1
       totrog( 3) = drog(2) ! alkene_2
       totrog( 4) = drog(3) ! xylene_1
       totrog( 5) = drog(3) ! xylene_2
       totrog( 6) = drog(4) ! cresol
       totrog( 7) = drog(5) ! toluene_1
       totrog( 8) = drog(5) ! toluene_2
       totrog( 9) = drog(6) ! monoterpene_1
       totrog(10) = drog(6) ! monoterpene_2

c *** to save operations because h_vap is currently constant
c     set a fixed value for tempcorr outside the loop

       tempcorr = tt1 * exp( hfac * tt2 )

c *** set SVOC concentrations  [ ppm ] -> [ ug / m**3 ]
c     and initialize ctoti from SVOCS

       threshold = 0.0

c *** initialize totorgsw and totorg
       totorgsw = OLDPOA   / mwpoa
       totorg   = totorgsw + OLDSOA_A / mwsoa_a + OLDSOA_B / mwsoa_b

c *** Initial guess of caer is computed as follows:
c     From eqn (8) of Schell et al. (2001)
c       caer = ctotf - c0 * (caer/MW) / totorg
c     Assuming totorg doesn't change during the timestep,
c       caer * (1 + c0/MW / totorg) = ctotf
c       caer = ctotf / ( 1 + c0/MW / totorg)

       do ii = 1, NCVAP
         ! tempcorr = tt1 * exp( h_vap(ii) * rgas1 * tt2)
         c0(ii)     = cstar(ii) * tempcorr
         ctoti(ii)  = SVOCS(ii) * vap_ppm2ug(ii) * convfac
         prod(ii)   = alpha(ii) * totrog(ii)
         ctotf(ii)  = ctoti(ii) + prod(ii)
         threshold  = threshold +  ctotf(ii) / c0(ii)
         caer(ii) = ctotf(ii) * totorg /            ! Initial
     &             (totorg + c0(ii)*mw_vap_m1(ii))  !  guess
       end do

c *** check if gas/particle equilibrium can be established
       if ( (threshold .gt. 1.0) .or. (totorgsw .gt. threshmin) ) then

c *** calculate new SOA by partitioning. This method uses a globally
c     convergent Newton-Raphson method coded by Dr Benedikt Schell
c     to solve the nonlinear quadratic system shown in eqn 8 of 
c     Schell et al:
c        A(i)  * caer(i) ** 2 + B * caer(i) + C(i) = 0.0,
c        where B(i) contains the sum of all caer(j),
c        for j not equal to i.

         call NEWT( LAYER, caer, NCVAP, check,
     &              ctotf, c0, mw_vap_m1, totorgsw )
         if( check ) then
c  Try again with initial guess of 50/50 gas/aerosol split.
           do ii = 1, NCVAP
             caer(ii) = 0.5 * ctotf(ii)
           end do
           call NEWT( LAYER, caer, NCVAP, check,
     &                 ctotf, c0, mw_vap_m1, totorgsw )
           if (check) then                           
             write(LOGDEV,*) ' problem in NEWT at LAYER = ', LAYER
           end if            
         end if

c *** Constrain caer to values between CONMIN and ctotf
         do ii = 1, NCVAP
           if (caer(ii) .lt. CONMIN ) then
            If(caer(ii) .lt. 0.0 ) then
              write(logdev,*) ' caer negative '
            end if
            caer(ii) = CONMIN
           end if

!          if( ( caer(ii) - ctotf(ii) ) .gt. TOLMIN ) then
           if( ( caer(ii) - ctotf(ii) ) .gt. CTOLMIN ) then
            write(LOGDEV,*) ' ii = ', ii
            write(LOGDEV,*) '  caer(ii) = ', caer(ii)
            write(LOGDEV,*) ' ctotf(ii) = ', ctotf(ii)
            caer(ii) = ctotf(ii)
           end if
         end do

       else

C *** threshold not exceeded; no material transferred to aerosol phase
         do ii = 1, NCVAP
           caer(ii) = 0.0
         end do

       end if    ! check on equilibrium threshold

c *** Calculate total SOA concentrations at end of time step

       totsoa = 0.0
       do ii= 1, NCVAP
        totsoa = totsoa + caer(ii)
       end do

       totsoa_b  = caer(nbio1) + caer(nbio2) 
       totsoa_a  = max(0.0, totsoa - totsoa_b )
       
C FSB December 04, 2003, ORGRATE and ORGBRATE eliminated in favor
c     of returning the newly equilibrated concentrations.

       SOA_A  = totsoa_a ! total anthropogenic SOA
       SOA_B  = totsoa_b ! total biogenic SOA


C *** convert SVOCS to [ ppm ] for output
c..04/03 update - save total conc (gas + aerosol) for
c..partitioning at next time step

       convfacm1 = 1.0 / convfac
       do ii = 1, NCVAP
         SVOCS(ii) = ctotf(ii) * vap_ug2ppm(ii) * convfacm1
       end do

       return
       end   subroutine orgaer3


c /////////////////////////////////////////////////////////////////////
c  SUBROUTINE HCOND3 calculates the size-dependent term in the 
c   condensational-growth rate expression for the 2nd and 3rd moments of
c   a lognormal aerosol mode using the harmonic mean method.  This code 
c   follows Section A2 of Binkowski & Shankar (1995).
c
c  KEY SUBROUTINES/FUNCTIONS CALLED:  none
c
c  REVISION HISTORY:
c     coded November 7, 2003 by Dr. Francis S. Binkowski
c
c     Revised November 20, 2003 by F. Binkowski to have am1 and 
c     am2 as inputs

c  REFERENCE:
c   1. Binkowski, F.S. and U. Shankar, The regional particulate matter
c      model 1. Model description and preliminary results, J. Geophys.
c      Res., Vol 100, No D12, 26101-26209, 1995.

      subroutine HCOND3( am0, am1, am2, Dv, alpha, cbar, F )

      implicit none
      
c *** input:

      real*8 am0  ! zeroth moment of mode  [ # / m**3 ]
      real*8 am1  ! first moment of mode   [ m / m**3 ]
      real*8 am2  ! second moment of mode  [ m**2 / m**3 ]
      real Dv     ! molecular diffusivity of the
                  ! condensing vapor  [ m**2 / s ]
      real alpha  ! accommodation coefficient           
      real cbar   ! kinetic velocity of condensing vapor [ m / s ]   
              
c *** output: 
                
      real*8 F(2) ! size-dependent term in condensational-growth rate
                  ! F(1) = 2nd moment [m^2 / m^3 s]
                  ! F(2) = 3rd moment [m^3 / m^3 s]

c *** internal

      real*8 gnc2 ! integrals used to calculate F(1) [m^2 / m^3 s]
      real*8 gfm2 ! 

      real*8 gnc3 ! integrals used to calculate F(2) [m^3 / m^3 s]
      real*8 gfm3 ! 
                  
      real*8 pi
       parameter( pi = 3.14159265358979324d0 )
      real*8 twopi
       parameter( twopi = 2.0d0 * pi )
      real*8 pi4
       parameter( pi4 = 0.25d0 * pi )
                   
c *** start code

C *** Implement equation A15 of Binkowski & Shankar (1995) for the
c     2nd and 3rd moments of a lognormal mode of arbitrary size.  

      gnc2 = twopi * Dv * am0          ! 2nd moment, near-continuum
      gnc3 = twopi * Dv * am1          ! 3rd moment, near-continuum
      gfm2 = pi4 * alpha * cbar * am1  ! 2nd moment, free-molecular
      gfm3 = pi4 * alpha * cbar * am2  ! 3rd moment, free-molecular

C *** Implement equation A13 of Binkowski & Shankar (1995) for a
c     lognormal mode of arbitrary size.  These are the size-dependent
c     terms in the condensational-growth rate expression, given in 
c     equation 7a of B&S (1995).

      F(1) = gnc2 * gfm2 / ( gnc2 + gfm2 )  ! 2nd moment
      F(2) = gnc3 * gfm3 / ( gnc3 + gfm3 )  ! 3rd moment
            
      return 
      
      end subroutine HCOND3


c /////////////////////////////////////////////////////////////////////
C *** this routine calculates the  sedimentation
c     velocities for the coarse mode.
C     coded 5/15/97 by Dr. Francis S. Binkowski.
c     modified to eliminate INCLUDE files  3/3/99 - FSB

       SUBROUTINE GETVSED(   NASPCSSED,
     &                    AIRTEMP, AIRDENS,
     &                    XLM, AMU,
     &                    DGCOR,
     &                    PDENSCO,
     &                    VSED            )

C *** calculate size-averaged particle dry sedimentation velocities.

c  REFERENCE:
c      Binkowski, F.S. and U. Shankar, The regional particulate matter
c      model 1. Model description and preliminary results, J. Geophys.
c      Res., Vol 100, No D12, 26101-26209, 1995.
c

      IMPLICIT NONE

C *** Input

      INTEGER NASPCSSED   ! number of sedimentation velocities

C Meteorological information in blocked arays:

      REAL AIRTEMP     ! Air temperature [ K ]
      REAL AIRDENS   ! Air density  [ kg m^-3 ]

C *** atmospheric properties

      REAL XLM   ! atmospheric mean free path [ m ]
      REAL AMU   ! atmospheric dynamic viscosity [ kg/m s ]

C *** aerosol properties:

      REAL DGCOR  ! coarse mode geometric mean diameter [ m ]

c *** average modal particle densities  [ kg/m**3 ]

      REAL PDENSCO  ! average particle density in coarse mode

C ***  sedimentation velocities

      REAL VSED( NASPCSSED)  ! sedimentation  velocity [ m/s ]



        INTEGER VSNCOR  ! index for coarse mode number
         PARAMETER( VSNCOR = 1)

        INTEGER VSMCOR  ! index for coarse mass
         PARAMETER( VSMCOR = 2)

      REAL DCONST2, DCONST3C
      REAL KNCOR                    ! coarse mode Knudsen Number

      REAL BHAT
      PARAMETER( BHAT =  1.246 )
                  ! Constant from Cunningham slip correction.

      REAL        GRAV ! mean gravitational acceleration [ m/sec**2 ]
                     ! NOTE: Value is now mean of polar and equatorial
                     ! values. Source: CRC Handbook (76th Ed) page 14-6.
      PARAMETER ( GRAV = 9.80622  )

C Scalar variables for fixed standard deviations.

      REAL    XXLSGCO          ! log(SGINICO )

      REAL    L2SGINICO        ! log(SGINICO ) ** 2

      REAL    EC1             ! coarse mode exp( log^2( sigmag )/8 )

      REAL    ESC04           ! coarse       "
      SAVE    ESC04

      REAL    ESC08           ! coarse       "
      SAVE    ESC08

      REAL    ESC16           ! coarse       "
      SAVE    ESC16

      REAL    ESC20           ! coarse       "
      SAVE    ESC20

      REAL    ESC28           ! coarse       "
      SAVE    ESC28

      REAL    ESC32           ! coarse       "
      SAVE    ESC32

      REAL    ESC36           ! coarse       "
      SAVE    ESC36

      REAL    ESC64           ! coarse       "
      SAVE    ESC64

      REAL    ESCM20          ! coarse       "
      SAVE    ESCM20

      REAL    ESCM32          ! coarse       "
      SAVE    ESCM32

      REAL        SGINICO        ! fixed l sigma-G for coarse mode
      PARAMETER ( SGINICO = 2.2)

      LOGICAL FIRSTIME
      SAVE FIRSTIME
      DATA FIRSTIME / .TRUE. /

C

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         XXLSGCO = LOG( SGINICO )
         L2SGINICO = XXLSGCO ** 2

         EC1   = EXP( 0.125 * L2SGINICO )

         ESC04  = EC1 ** 4

         ESC08  = ESC04 * ESC04

         ESC16  = ESC08 * ESC08

         ESC20  = ESC16 * ESC04

         ESC28  = ESC20 * ESC08

         ESC32  = ESC16 * ESC16

         ESC36  = ESC16 * ESC20

         ESC64  = ESC32 * ESC32

         ESCM20 = 1.0 / ESC20

         ESCM32 = 1.0 / ESC32

      END IF

C *** begin code ------------------------------------------------------*

C ***  calculate  sedimentation velocities
c FSB See Section A4 Dry Deposition of Binkowski & Shankar (1995)
            KNCOR = 2.0 * XLM / DGCOR

            DCONST2 = GRAV / ( 18.0 * AMU )

            DCONST3C = DCONST2 * PDENSCO * DGCOR * DGCOR

C *** coarse mode sedimentation velocity for number

           VSED(  VSNCOR) = DCONST3C
     &        * ( ESC16 + BHAT * KNCOR * ESC04 )


C *** coarse mode sedimentation velocity for mass

            VSED( VSMCOR) = DCONST3C
     &        * ( ESC64 + BHAT * KNCOR * ESC28 )

      END SUBROUTINE GETVSED
      
      
c ///////////////////////////////////////////////
      SUBROUTINE VSNdist5( Ni, dgni, xlsgi,
     &                     Nj, dgnj, xlsgj,
     &                     Nc, dgnc, xlsgc,
     &                     dia, vol, num, sfc,
     &                     NPOINTS )

c *** subroutine to generate log-normal number, surface area,
c     and volume  distributions from the
c     total number, geometric mean diameter, and log of
c     the geometric standard deviation for each of the three modes.
c     Aitken mode:       Ni, dgni, xlsgi
c     accumulation mode: Nj, dgnj, xlsgj
c     coarse mode:       Nc, dgnc, xlsgc
c *** the results are returned in arrays
c     dia - diameter
c     vol - volume
c     num - number
c     sfc - surface area
c     which are of dimension NPOINTS.
C     NOTE:  NPOINTS should be evenly divisbile by 4.

c     this version coded by Dr. Francis S. Binkowski
c     July 18, 1999.

c *** NOTE: This version has xlsgi, xlsgj & xlssgc as inputs

c *** This version is developed from VSNdist3.f and now includes
c     the coarse mode, i.e. it is now trimodal.

c *** FURTHER NOTE: this version now uses base 10 logs for
c     coni, conj, and conc. This rescales the size plots to the more
c     prevalent usage in the literature.

C *** This code  plots the size distributions by interpolating
c     particle diameters, Dp, over the four domains of diameter
c     defined as:

c *** Domain one   Dpmin  -> dgni
c     Domain two   dgni   -> dgnj
c     Domain three dgnj   -> dgnc
c     Domain four  dgnc   -> Dpmax


      implicit  none

C *** input:

c *** input units

      real  dgni,dgnj, dgnc    ! modal geometric
                               ! mean diameters [ um ]

      real  Ni,Nj,Nc           ! modal number
                               ! concentrations [ # / cm **3 ]

      real  xlsgi,xlsgj, xlsgc ! natural log of geometric
                               ! standard deviations

      integer NPOINTS          ! the number of points in the
                               ! output arrays

c *** output:

c *** distributions:

      real dia(NPOINTS),
     &     vol(NPOINTS),
     &     num(NPOINTS),
     &     sfc(NPOINTS)  ! standard cgs units

C *** local:
      integer ii, jj, kk

      integer MPOINTS ! The number of diameters assigned to each domain
c *** MPOINTS is set to NPOINTS / below

       real ai, aj, ac, ak,
     &      coni, conj, conc, p, q, r,
     &      Dpmin, Dpmax, xl(5),
     &      n, S, V, Dp, xdp,dxdp

       real pfac


      real pi
      parameter(  pi = 3.1415926536 )

      real const1
      parameter( const1 = 0.918599  )
c      const1 = log(10.0) / sqrt( 2.0 * pi )
c      this converts xlsgi, xlsgj, & xlsgc
c      to common logs.

c *** start code

      MPOINTS = NPOINTS / 4

      pfac = 1.0 / float(MPOINTS)

c *** Set up general constants for each mode:


      coni = const1 * Ni / xlsgi
      conj = const1 * Nj / xlsgj
      conc = const1 * Nc / xlsgc


C *** calculate the minimum and maximum diameters
c *** the minimum diameter is three times xlsgi smaller than dgni
c
c      Dpmin = dgni * exp( -3.0 * xlsgi )

C *** changed 10/24/2000 FSB set Dpmin to a constant

      Dpmin = 0.01

c *** Dpmax set to a constant of 10.0  [ um ]
      Dpmax = 10.0

c *** set boundaries of the four domains


      xl(1) = Dpmin
      xl(2) = dgni
      xl(3) = dgnj
      xl(4) = dgnc
      xl(5) = Dpmax

c *** initialize index for storage of results

       ii = 1

c *** loop over the four domains

      do kk = 1,4


C *** set initial diameter

      xdp = xl(kk)

c *** set difference in diameters

      dxdp = xl(kk+1) - xl(kk)

      do jj = 1, MPOINTS

c *** interpolate the diameters to Dp

      ak = pfac * float(jj)
      Dp = xdp + ak * dxdp


C *** calculate arguments of the exponentials at Dp:
      p = log( Dp / dgni ) / xlsgi
      q = log( Dp / dgnj ) / xlsgj
      r = log( Dp / dgnc ) / xlsgc


C *** calculate the exponentials:
      ai = exp( -0.5 * p * p )
      aj = exp( -0.5 * q * q)
      ac = exp( -0.5 * r * r )


C *** calculate the individual values of the distributions at Dp:

C *** number at Dp

      n = ( coni * ai + conj * aj + conc * ac)


C *** surface area at Dp
      S = pi * Dp * Dp * n


C *** Volume at Dp
      V = Dp * S / 6.0


C *** store these values in appropriate arrays

      num(ii) = n
      sfc(ii) = S
      vol(ii) = V
      dia(ii) = Dp

c *** increment storage index

      ii = ii + 1

      end do ! loop over points

      end do ! loop over four domains

      return
      end     SUBROUTINE VSNdist5

c ///////////////////////////////////////////////////

      SUBROUTINE GETVISBY( NSPCSDA,
     &               CBLK,AIRRH,
     &               DCV1, EXT1, DCV2, EXT2,
     &               DGATK, DGACC, DGCOR,
     &               XXLSGAT, XXLSGAC,
     &               PMASSAT, PMASSAC,PMASSCO )

c
c *** this routine calculates the Pitchford & Malm visual index,
c     deciview, using the Evans & Fournier extinction with 20 point
c     Gauss_Hermite numerical quadrature to integrate the extincetion
c     over the log normal size distribution.
C     This new method reproduces the results of Willeke and Brockmann
c     very very closely.
c     The Heintzenberg and Baker (h&b) method used previously has been
c     replaced.
c *** This also routine calculates the Pitchford & Malm visual index,
c     deciview, using reconstructed extinction values from the IMPROVE
c      monitoring network

c *** references:

c      Evans, B.T.N.  and G.R. Fournier, Simple approximations to
c       extinction efficiency valid over all size parameters,
C       Applied Optics, 29, 4666 - 4670.

c     Heintzenberg, j. and M. Baker, Applied Optics, vol. 15, no. 5,
c     pp 1178-1181, 1976.
c      correction in Applied Optics August 1976 does not affect this code.

c      Sisler, J. Fax dated March 18, 1998 with IMPROVE information.

c      Pitchford, M. and W. Malm, Atmos. Environ.,vol 28,no.5,
c       pp 1049-1054, 1994.

c      Willeke,K. & J. E. Brockmann, Atmos. Environ., vol.11,
c       pp 995 - 999, 1977.
c
coding history:
c    who           when     what
c ------------------------------------------------------------------
c   F.S. Binkowski 5/9/95
c             coded this version to use the h&b approximation
c             and alfa and c obtained from fits to adt results
c             b was fit from comparisons to willeke & brockmann.
c   F.S. Binkowski 9/12/95  modified code to push j-loop inside.
c   F.S. Binkowski 5/13/97  made models-3 version
c   F.S. Binkowski 4/20/98  changed to Evans and Fournier approach
c                           for extinction
c   F.S. Binkowski 4/20/98  began code for reconstructed method
c   F.S. Binkowski 4/24/98  merged codes for both methods
c   F.S. Binkowski 5/20/98  corrected CONSTL
C   F.S. Binkowski 3/9/99   Modified for variable XXLSGAT and XXLSGAC
c ////////
C NOTE: This subroutine is dependent upon the variables defined in the
c       IF( FIRSTIME ) section of aeroproc.f for implementation of
c                      coarse mode contribution.
c ///////

      USE       AERO_INFO_AE4

      IMPLICIT NONE


      INTEGER NSPCSDA           ! number of species in CBLK
      REAL CBLK( NSPCSDA  )  ! main array of variables
                                !( INPUT and OUTPUT )
      REAL AIRRH          !  fractional relative humidity

      REAL DCV1   ! block deciview (Mie)
      REAL EXT1   ! block extinction [ km**-1 ] (Mie)

      REAL DCV2   ! block deciview (Reconstructed)
      REAL EXT2   ! block extinction [ km**-1 ]
                               ! (Reconstructed)

C *** modal geometric mean diameters: [ m ]

      REAL DGATK          ! nuclei mode
      REAL DGACC          ! accumulation mode
      REAL DGCOR          ! coarse mode

c *** log of modal geometric standard deviation

      REAL XXLSGAT         ! Aitken mode
      REAL XXLSGAC         ! accumulation mode

C *** aerosol properties:

c *** Modal mass concentrations [ ug m**3 ]

      REAL PMASSAT         ! Aitken mode
      REAL PMASSAC         ! accumulation mode
      REAL PMASSCO         ! coarse mode

C *** internal variables:

      INTEGER IRH ! percent realtive humidity as an
                  ! integer used for index

      REAL        CONMIN
       PARAMETER ( CONMIN = 1.0E-30 ) ! concentration lower limit

      REAL LAM ! wavelenght of visible light [ m ]
       PARAMETER( LAM = 0.55E-6 )


      REAL CONSTL
c       PARAMETER ( CONSTL = 1.5E3 * PI / LAM  ) !
C *** constant CONSTL corrected 5/20/98

        PARAMETER( CONSTL = 1.0E3 * PI6 / LAM ) ! 1.0e3  to get km**-1

      REAL CONST3
       PARAMETER( CONST3 = PI / LAM ) ! Changed 3/9/99 FSB

      REAL WFRAC  ! water mass fraction
      REAL NR, NI ! real and imaginary parts of the refracive index
      REAL ALFVN, ALFVA, ALFVC ! Mie parameters for modal mass median
      !                           diameters
      REAL BBEXT  ! dimensionless extinction coefficient
      REAL    BEXTN, BEXTA, BEXTC
                    ! Modal extinction coefficients [ 1/km ]

      REAL SCALE  ! factor to rescale units from [ 1/Mm ] to [ 1/km ]
      PARAMETER ( SCALE = 1.0E-03 )
      REAL FRH ! RH correction factor for sulfate and nitrate aerosols

      REAL RAY, RAY1 ! standard value for Rayleigh
                     ! extinction [ 1/km ] and its reciprocal.
      PARAMETER ( RAY = 0.01 )
      PARAMETER ( RAY1 = 1.0 / RAY )

      REAL humfac(99) ! humidity scaling factors at 1% RH values
         DATA humfac/
     &   1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,
     &   1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0000E+00,
     &   1.0000E+00,  1.0000E+00,  1.0000E+00,  1.0001E+00,  1.0001E+00,
     &   1.0004E+00,  1.0006E+00,  1.0024E+00,  1.0056E+00,  1.0089E+00,
     &   1.0097E+00,  1.0105E+00,  1.0111E+00,  1.0115E+00,  1.0118E+00,
     &   1.0122E+00,  1.0126E+00,  1.0130E+00,  1.0135E+00,  1.0139E+00,
     &   1.0173E+00,  1.0206E+00,  1.0254E+00,  1.0315E+00,  1.0377E+00,
     &   1.0486E+00,  1.0596E+00,  1.0751E+00,  1.0951E+00,  1.1151E+00,
     &   1.1247E+00,  1.1343E+00,  1.1436E+00,  1.1525E+00,  1.1615E+00,
     &   1.1724E+00,  1.1833E+00,  1.1955E+00,  1.2090E+00,  1.2224E+00,
     &   1.2368E+00,  1.2512E+00,  1.2671E+00,  1.2844E+00,  1.3018E+00,
     &   1.3234E+00,  1.3450E+00,  1.3695E+00,  1.3969E+00,  1.4246E+00,
     &   1.4628E+00,  1.5014E+00,  1.5468E+00,  1.5992E+00,  1.6516E+00,
     &   1.6991E+00,  1.7466E+00,  1.7985E+00,  1.8549E+00,  1.9113E+00,
     &   1.9596E+00,  2.0080E+00,  2.0596E+00,  2.1146E+00,  2.1695E+00,
     &   2.2630E+00,  2.3565E+00,  2.4692E+00,  2.6011E+00,  2.7330E+00,
     &   2.8461E+00,  2.9592E+00,  3.0853E+00,  3.2245E+00,  3.3637E+00,
     &   3.5743E+00,  3.7849E+00,  4.0466E+00,  4.3594E+00,  4.6721E+00,
     &   5.3067E+00,  5.9412E+00,  6.9627E+00,  8.3710E+00,  9.7793E+00,
     &   1.2429E+01,  1.5078E+01,  1.8059E+01,  2.1371E+01/

         SAVE humfac

C NOTE: In the following calculations, the contribution from the
c        coarse mode is ignored.

c *** calculate the  mass fraction of aerosol water

        WFRAC = MIN( ( CBLK( VH2OAJ) + CBLK( VH2OAI) )
     &              / ( PMASSAT  + PMASSAC  ), 1.0)

c *** interpolate between "dry" state with m = 1.5 - 0.01i
c      and pure water particle with  m = 1.33 - 0.0i as a function of
c       wfrac

        NR = 1.5 - 0.17 * WFRAC   ! real part of refracive index
        NI = 0.01 * (1.0 - WFRAC) ! imaginary part of refractive index

c *** set up Mie parameters for Volume ( mass median diameter)

        ALFVN = CONST3 * DGATK  *
     &        EXP( 3.0 * XXLSGAT  * XXLSGAT  )

        ALFVA = CONST3 * DGACC  *
     &        EXP( 3.0 * XXLSGAC  * XXLSGAC  )

        ALFVC = CONST3 * DGCOR  * ESC24

C *** Call extinction routines

        CALL getbext(NR, NI, ALFVN, XXLSGAT , BBEXT)
        BEXTN = CONSTL * CBLK( VAT3) * BBEXT

        CALL getbext(NR, NI, ALFVA, XXLSGAC , BBEXT)
        BEXTA = CONSTL * CBLK( VAC3) * BBEXT

C        CALL getbext(NR, NI, ALFVC, XXLSGCO, BBEXT)
C        BEXTC = CONSTL * CBLK( VCOR3) * BBEXT
         BEXTC = 0.0

      EXT1  = BEXTN + BEXTA + RAY

      DCV1  = 10.0 * LOG ( EXT1 * RAY1 )

c     note if EXT1 < 0.01 then DCV1 is negative.
c     this implies that visual range is greater than the Rayleigh limit.
c     The definition of deciviews is based upon the Rayleigh limit
c     being the maximum visual range
c     thus, set a floor of 0.0 on DCV1.

      DCV1  = MAX(0.0,DCV1 )


c *** begin  IMPROVE reconstructed method


      IRH = INT( 100.0 * AIRRH  ) ! truncate relative humidity to
                                    ! nearest integer
      IRH = MIN(99, IRH) ! set maximum value on IRH
      IRH = MAX(1,IRH) ! set minimum value on IRH

      FRH = humfac(IRH) ! set humidity correction

C *** NOTE in the following the fine primary mass "other" is
c          treated as though it were fine mass soil.

      EXT2  = SCALE * (
     &            3.0 * FRH * ( CBLK(  VSO4AI) + CBLK(  VSO4AJ)
     &                      +   CBLK(  VNO3AI) + CBLK( VNO3AJ)
     &                      +   CBLK(  VNH4AI) + CBLK(  VNH4AJ) )
     &         + 4.0 * (        CBLK(  VORGAI) + CBLK(  VORGAJ)
     &                      +   CBLK(  VORGPAI) + CBLK( VORGPAJ)
     &                      +   CBLK(  VORGBAI) + CBLK( VORGBAJ) )
     &         + 10.0 * (       CBLK( VECI)  +  CBLK( VECJ) )
     &         + 1.0 * (        CBLK( VP25AI) + CBLK(  VP25AJ) )
c     &         + 0.6 * (        CBLK( VSOILA) + CBLK(  VANTHA) )
     &                                                             )
     &        + RAY

      DCV2  = 10.0 * LOG ( EXT2  * RAY1 )

c     note if EXT2 < RAY then DCV2 is negative.
c     this implies that visual range is greater than the Rayleigh limit.
c     The definition of deciviews is based upon the Rayleigh limit
c     being the maximum visual range
c     thus, set a floor of 0.0 on BLKDCV.

      DCV2  = MAX(0.0,DCV2 )

      RETURN
      END SUBROUTINE GETVISBY


c ************************************************
       subroutine getbext(nr, ni, alfv, xlnsig, bext)

c      calculates the extinction coefficient normalized by wavelength
c      and total particle volume concentration for a log normal
c      particle distribution with the logarithm of the
c     geometric standard deviation given by xlnsig.

c *** does gauss-hermite quadrature of Qext / alfa
c  over log normal distribution
c     using 20 symmetric points
c
      implicit none

      real nr, ni  ! indices of refraction
      real alfv    ! Mie parameter for dgv
      real xlnsig  ! log of geometric standard deviation
      real bext    ! normalized extinction coefficient
      real aa, aa1 ! see below for definition

      real alfaip, alfaim   ! Mie parameters at abscissas

      real  xxqalf  ! function to calculate the extinction ceofficient
      real qalfip, qalfim   ! extinction efficiencies at abscissas

      real pi
      parameter( pi = 3.14159265)

      real sqrtpi
      parameter( sqrtpi = 1.772454 )

      real sqrtpi1
      parameter( sqrtpi1 = 1.0 / sqrtpi )

      real sqrt2
      parameter( sqrt2 = 1.414214 )

      real three_pi_two
      parameter( three_pi_two = 3.0 * pi / 2.0 )

      real const
      parameter( const = three_pi_two * sqrtpi1 )


      integer i
      real sum, xi,wxi,xf

      integer n    ! one-half the number of abscissas
      parameter ( n = 10 )
      real ghxi(10) ! Gauss-Hermite abscissas
      real ghwi(10) ! Gauss-Hermite weights

c *** the following weights and abscissas are from abramowitz and
c     stegun, page 924
c *** tests show that 20 point is adquate.

      data ghxi/0.245341,0.737474,1.234076,1.738538,2.254974,
     & 2.788806,3.347855,3.944764,4.603682,5.387481/

      data ghwi/4.622437e-1,2.866755e-1,1.090172e-1,2.481052e-2,
     & 3.243773e-3,2.283386e-4,7.802556e-6,1.086069e-7,
     4 4.399341e-10,2.229394e-13/

      sum=0.0

      aa = 1.0 / ( sqrt2 * xlnsig )
      aa1 = sqrt2 * xlnsig ! multiplication cheaper
                           ! than another division

      do  i = 1, n

       xi      = ghxi(i)
       wxi     = ghwi(i)
       xf      = exp( xi * aa1 )
       alfaip  = alfv * xf
       alfaim  = alfv / xf
       qalfip  = xxqalf(alfaip, nr, ni)
       qalfim  = xxqalf(alfaim, nr, ni)

       sum = sum + wxi * ( qalfip + qalfim )

      end do ! i

c fsb      bext = const * aa * sum
      bext = const * sum ! corrected 07/21/2000 FSB
                         ! found by Rokjin Park

      return
      end subroutine getbext
c ________________________
       real function xxqalf(alfa,nr,ni)
c ***  compute the extinction efficiency divided by the Mie parameter
c      reference:
c      Evans, B.T.N.  and G.R. Fournier, Simple approximations to
c      extinction efficiency valid over all size parameters,
C      Applied Optics, 29, 4666 - 4670.

       implicit none

       real nr, ni  ! real and imaginary parts of index of refraction
       real qextray
       real qextadt
       real qextef



       real tt        !  edge effect factor
       real mu, mum1  !  exponents in formula
       real aa        !  first coefficient in mu
       real gg        !  second coefficient in mu
       real nrm1, sqrtni
       real alfa
       real alfm23   ! functions of alfa (mie parameter)
       real three5, three4, two3
        parameter( three5    = 3.0 / 5.0   ,
     &             three4    = 3.0 / 4.0   ,
     &             two3      = 2.0 / 3.0   )


       nrm1   = nr - 1.0
       sqrtni = sqrt (ni )

       call adtqext(alfa,nr,ni, qextadt)
       call pendrfx(alfa,nr,ni,qextray)

       if( alfa .gt. 0.5 ) then

         alfm23  = 1.0 / alfa ** two3

         tt      = 2.0 - exp( -alfm23 )

         aa      = 0.5 + ( nrm1 - two3 * sqrtni - 0.5 * ni ) +
     &           ( nrm1 + two3 * ( sqrtni - 5.0  * ni )  ) ** 2

         gg      = three5 - three4 * sqrt ( nrm1) + 3.0 * nrm1 ** 4 +
     &            25.0 * ni / ( 6.0 * ni + 5.0 * nrm1 )

         mu      = aa + gg / alfa
         mum1    = - 1.0 / mu

         qextef  = qextray *
     &           ( 1.0 + (qextray /( qextadt * tt)) ** mu ) ** mum1

        else

         qextef  = qextray ! Use Rayleigh extinction for
                           ! really small alfa's

        end if ! check on alfa

        xxqalf  =  qextef / alfa

        return
       end function xxqalf
c ______________________________________
      subroutine adtqext(alfa, nr, ni, QEXT)

c *** van de Hulst approximation for QEXT.
c ***  This approximation is known as Anomalous Diffraction Theory (ADT)

c *** originally coded by Dr Francis  S. Binkowski,
c        AMAB/MD/ESRL RTP,N.C. 27711 28 July 1977.
c       corrected 7/19/90 by fsb.
c       revised 1/8/98 by FSB


C *** reference:
c         van de Hulst- Light Scattering by Small Particles,
c         Dover,1981 page 179. Original edition was Wiley, 1957.

        implicit none
        real alfa   ! Mie parameter
        real nr, ni ! real and imaginary parts of the index of refraction
        real QEXT   ! extinction efficiency for a sphere
        real z, tanb, b, cos2b, v1, v2, x, expmx, expm2x, cs1, cs2
        real twob, cosb
        real nr1

        nr1    = nr - 1.0
        z     = 2.0 * alfa * nr1
        tanb  = ni/nr1
        b     = atan(tanb)
        cosb  = cos(b)
        twob  = 2.0 * b
        cos2b = cos(twob)
        v1     = 5.0 * nr1
        v2     = 4.08 / (1.0 + 3.0 * tanb)
        x      = z * tanb
        expmx  = exp( -x )
        expm2x = exp( -2.0 * x )
        cs1    = cosb / z
        cs2    = cs1 * cs1

	    QEXT   = 2.0 - 4.0 * cs1 * expmx * sin( z - b) +
     &           4.0 * cs2 * ( cos2b - expmx * cos(z - twob))
        return
        end subroutine adtqext
c __________________________________________________
      subroutine pendrfx(alfa, nr, ni, QEXT)
c *** calculates the Mie efficiencies for extinction, scattering and
c     absorption using Penndorf's approximations for small
c     values of alfa.
c
c  input:
c       nr        the real part of the refractive index.
c       ni        the imaginary part of the refractive index
c       alfa      mie parameter

c
c  output:
c       QEXT      extinction efficiency

c
c *** reference:
c       Penndorf, r., Scattering and extinction coefficients for small
c       aerosols, J. Atmos. Sci., 19, p 193, 1962.
c
c *** coded by Dr Francis  S.  Binkowski,
c       AMAB/MD/ESRL RTP,N.C. 27711 28 July 1977.
c       corrected 7/19/90 by FSB
c       modified 30 september 1992 by FSB
c       modified 1/6/98 by FSB

        implicit none

        real alfa, nr, ni
        real QEXT
        real alf2,alf3,alf4

        real a1,a2,a3
        real xnr,xni,xnr2,xni2,xnri,xnri2,xnrmi
        real xri,xri2,xri36,xnx,xnx2
	    real z1,z12,z2,xc1

        xnr   = nr
        xni   = ni
        xnr2  = xnr   * xnr
        xni2  = xni   * xni
        xnri  = xnr2  + xni2
        xnri2 = xnri  * xnri
        xnrmi = xnr2  - xni2
        xri   = xnr   * xni
        xri2  = xri   * xri
        xri36 = 36.0  * xri2
        xnx   = xnri2 + xnrmi - 2.0
        xnx2  = xnx   * xnx

        z1    = xnri2 + 4.0 * xnrmi + 4.0
        z12   = z1    * z1
        z2    = 4.0   * xnri2 + 12.0 * xnrmi + 9.0
        xc1   = 8.0   / (3.0*z12)

        alf2  = alfa  * alfa
        alf3  = alfa  * alf2
        alf4  = alfa  * alf3

        a1    = 24.0  * xri / z1

        a2    = 4.0   * xri / 15.0 + 20.0 * xri/(3.0 * z2) +
     &              4.8   * xri * (7.0 * xnri2 +
     &              4.0   * ( xnrmi - 5.0) ) / z12

        a3    = xc1   * (xnx2 - xri36)

	    QEXT  = a1    * alfa + a2 * alf3 + a3 * alf4

        return
        end subroutine pendrfx
         ! end of visibility codes
c //////////////////////////////////////

       SUBROUTINE GETVDEP( N_AE_DEP_SPC,
     &                    AIRTEMP, AIRDENS,
     &                    XLM, AMU,
     &                    WSTAR, USTAR, RA,
     &                    DGATK, DGACC, DGCOR,
     &                    XXLSGAT, XXLSGAC,
     &                    PDENSAT, PDENSAC, PDENSCO,
     &                    VDEP            )

c *** calculate deposition velocity for Aitken, accumulation, and
c     coarse modes.
c     Reference:
c     Binkowski F. S., and U. Shankar, The regional particulate
c     model 1. Model description and preliminary results.
c     J. Geophys. Res., 100, D12, 26191-26209, 1995.

      IMPLICIT NONE

      INTEGER N_AE_DEP_SPC

C *** arguments:

c *** input

C Meteorological information in blocked arays:

      REAL AIRTEMP     ! Air temperature [ K ]
      REAL AIRDENS   ! Air density  [ kg m^-3 ]

C *** atmospheric properties

      REAL XLM    ! atmospheric mean free path [ m ]
      REAL AMU    ! atmospheric dynamic viscosity [ kg/m s ]

c *** Planetary boundary laryer (PBL) variables:

      REAL WSTAR  ! convective velocity scale [ m s**-1 ]
      REAL USTAR  ! friction velocity [ m s**-1 ]
      REAL RA     ! aerodynamic resistance [ s m**-1 ]

C *** aerosol properties:

C *** modal geometric mean diameters: [ m ]

      REAL DGATK     ! Atken mode
      REAL DGACC     ! accumulation mode
      REAL DGCOR     ! coarse mode

c *** log of modal geometric standard deviations

      REAL XXLSGAT         ! Aitken mode
      REAL XXLSGAC         ! accumulation mode


c *** average modal particle densities  [ kg/m**3 ]

      REAL PDENSAT         ! Aitken mode
      REAL PDENSAC         ! accumulation mode
      REAL PDENSCO         ! coarse mode

C *** Output

      REAL VDEP( N_AE_DEP_SPC ) ! deposition  velocity [ m/s ]

C *** internal:

      INTEGER VDNATK            ! Aitken mode number
      PARAMETER ( VDNATK = 1 )

      INTEGER VDNACC            ! accumulation mode number
      PARAMETER ( VDNACC = 2 )

      INTEGER VDNCOR            ! coarse mode number
      PARAMETER ( VDNCOR = 3 )

      INTEGER VDMATK            ! Aitken mode mass
      PARAMETER ( VDMATK = 4 )

      INTEGER VDMACC           !  accumulation mode mass
      PARAMETER ( VDMACC = 5 )

      INTEGER VDMCOR            ! coarse mode mass
      PARAMETER ( VDMCOR = 6 )

      INTEGER VDSATK           ! Aitken mode surface area
      PARAMETER ( VDSATK = 7 )

      INTEGER VDSACC           ! accumulation mode surface area
      PARAMETER ( VDSACC = 8 )

C model Knudsen numbers

      REAL KNATK   ! Aitken mode Knudsen number
      REAL KNACC   ! accumulation "
      REAL KNCOR   ! coarse mode

C modal particle diffusivities for number and 3rd moment, or mass:

      REAL DCHAT0AT , DCHAT0AC , DCHAT0C
      REAL DCHAT2AT , DCHAT2AC
      REAL DCHAT3AT , DCHAT3AC , DCHAT3C

C     modal sedimentation velocities for number, surface area,
C         and 3rd moment, or mass:

      REAL VGHAT0AT , VGHAT0AC , VGHAT0C
      REAL VGHAT2AT , VGHAT2AC
      REAL VGHAT3AT , VGHAT3AC , VGHAT3C

      REAL DCONST1, DCONST1AT, DCONST1AC, DCONST1C
      REAL DCONST2, DCONST3AT, DCONST3AC, DCONST3C
      REAL SC0AT, SC0AC, SC0C  ! Schmidt numbers for number
      REAL SC2AT, SC2AC        ! Schmidt numbers for 2ATD MOMENT
      REAL SC3AT, SC3AC, SC3C  ! Schmidt numbers for 3rd moment
      REAL ST0AT, ST0AC        ! Stokes numbers for number
      REAL ST2AT, ST2AC        ! Stokes numbers for 2nd moment
      REAL ST3AT, ST3AC        ! Stokes numbers for 3rd moment
      REAL RD0AT, RD0AC, RD0C  ! canopy resistance for number
      REAL RD2AT, RD2AC        ! canopy resistance for 2nd moment
      REAL RD3AT, RD3AC, RD3C  ! canopy resisteance for 3rd moment
      REAL UTSCALE           ! scratch function of USTAR and WSTAR
      REAL NU                ! kinematic viscosity [ m**2 s**-1 ]
      REAL USTFAC            ! scratch function of USTAR, NU, and GRAV
      REAL       BHAT        ! Constant from Cunningham slip correction
      PARAMETER( BHAT = 1.246 )

      REAL      PI ! PI (single precision 3.141593)
      PARAMETER ( PI = 3.141593 )

      REAL        PI6           ! PI/6
      PARAMETER ( PI6 = PI / 6.0  )

      REAL        THREEPI       !  3*PI
      PARAMETER ( THREEPI  = 3.0 * PI )

      REAL        ONE3          !  1/3
      PARAMETER ( ONE3     = 1.0 / 3.0 )

      REAL        TWO3          !  2/3
      PARAMETER ( TWO3     = 2.0 / 3.0 )

      REAL        AVO ! Avogadro's Constant [ 1/mol ]
      PARAMETER ( AVO = 6.0221367 E23  )

      REAL        RGASUNIV ! universal gas constant [ J/mol-K ]
      PARAMETER ( RGASUNIV = 8.314510  )

      REAL        BOLTZ         ! Boltzmann's Constant [ J/K]
      PARAMETER ( BOLTZ = RGASUNIV / AVO )

      REAL        GRAV ! mean gravitational acceleration [ m/sec**2 ]
C FSB               NOTE: Value is now mean of polar and equatorial
                   ! values. Source: CRC Handbook (76th Ed) page 14-6.
      PARAMETER ( GRAV = 9.80622  )

C Scalar variables for fixed standard deviations.

      REAL    XXLSGCO          ! log(SGINICO )

      REAL    L2SGINICO        ! log(SGINICO ) ** 2

      REAL    EC1             ! coarse mode exp( log^2( sigmag )/8 )

      REAL    ESC04           ! coarse       "
      SAVE    ESC04

      REAL    ESC08           ! coarse       "
      SAVE    ESC08

      REAL    ESC16           ! coarse       "
      SAVE    ESC16

      REAL    ESC20           ! coarse       "
      SAVE    ESC20

      REAL    ESC28           ! coarse       "
      SAVE    ESC28

      REAL    ESC32           ! coarse       "
      SAVE    ESC32

      REAL    ESC36           ! coarse       "
      SAVE    ESC36

      REAL    ESC64           ! coarse       "
      SAVE    ESC64

      REAL    ESCM20          ! coarse       "
      SAVE    ESCM20

      REAL    ESCM32          ! coarse       "
      SAVE    ESCM32

      REAL        SGINICO    ! fixed l sigma-G for coarse mode
      PARAMETER ( SGINICO = 2.2)

      REAL        DGINIAT    ! initial mean diam. for Aitken mode [ m ]
      PARAMETER ( DGINIAT = 0.01E-6  )


C Scalar variables for  VARIABLE standard deviations.

      REAL    L2SGAT, L2SGAC  ! see usage

      REAL    EAT1         ! Aitken mode exp( log^2( sigmag )/8 )
      REAL    EAC1         ! accumulation mode exp( log^2( sigmag )/8 )

      REAL    ESAT04
      REAL    ESAC04

      REAL    ESAT08
      REAL    ESAC08

      REAL    ESAT12
      REAL    ESAC12

      REAL    ESAT16
      REAL    ESAC16

      REAL    ESAT20
      REAL    ESAC20

      REAL    ESAT28
      REAL    ESAC28

      REAL    ESAT32
      REAL    ESAC32

      REAL    ESAT36
      REAL    ESAC36

      REAL    ESAT48
      REAL    ESAC48

      REAL    ESAT64
      REAL    ESAC64

      REAL    ESATM12
      REAL    ESACM12

      REAL    ESATM16
      REAL    ESACM16

      REAL    ESATM20
      REAL    ESACM20

      REAL    ESATM32
      REAL    ESACM32


      LOGICAL FIRSTIME
      SAVE FIRSTIME
      DATA FIRSTIME / .TRUE. /

C

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         XXLSGCO = LOG( SGINICO )
         L2SGINICO = XXLSGCO ** 2
         EC1   = EXP( 0.125 * L2SGINICO )
         ESC04  = EC1 ** 4
         ESC08  = ESC04 * ESC04
         ESC16  = ESC08 * ESC08
         ESC20  = ESC16 * ESC04
         ESC28  = ESC20 * ESC08
         ESC32  = ESC16 * ESC16
         ESC36  = ESC16 * ESC20
         ESC64  = ESC32 * ESC32
         ESCM20 = 1.0 / ESC20
         ESCM32 = 1.0 / ESC32

      END IF

         KNATK = 2.0 * XLM / DGATK
         KNACC = 2.0 * XLM / DGACC
         KNCOR = 2.0 * XLM / DGCOR

C *** Calculate functions of variable standard deviation.

         L2SGAT = XXLSGAT * XXLSGAT
         L2SGAC = XXLSGAC * XXLSGAC

         EAT1   = EXP( 0.125 * L2SGAT )
         EAC1   = EXP( 0.125 * L2SGAC )

         ESAT04  = EAT1 ** 4
         ESAC04  = EAC1 ** 4

         ESAT08  = ESAT04 * ESAT04
         ESAC08  = ESAC04 * ESAC04

         ESAT12  = ESAT04 * ESAT08
         ESAC12  = ESAC04 * ESAC08

         ESAT16  = ESAT08 * ESAT08
         ESAC16  = ESAC08 * ESAC08

         ESAT20  = ESAT16 * ESAT04
         ESAC20  = ESAC16 * ESAC04

         ESAT28  = ESAT20 * ESAT08
         ESAC28  = ESAC20 * ESAC08

         ESAT32  = ESAT16 * ESAT16
         ESAC32  = ESAC16 * ESAC16

         ESAT36  = ESAT16 * ESAT20
         ESAC36  = ESAC16 * ESAC20

         ESAT48  = ESAT36 * ESAT12
         ESAC48  = ESAC36 * ESAC12

         ESAT64  = ESAT32 * ESAT32
         ESAC64  = ESAC32 * ESAC32

C *** calculate inverses:

         ESATM12 = 1.0 / ESAT12
         ESACM12 = 1.0 / ESAC12

         ESATM16 = 1.0 / ESAT16
         ESACM16 = 1.0 / ESAC16

         ESATM20 = 1.0 / ESAT20
         ESACM20 = 1.0 / ESAC20

         ESATM32 = 1.0 / ESAT32
         ESACM32 = 1.0 / ESAC32

         DCONST1 = BOLTZ * AIRTEMP /
     &              ( THREEPI * AMU )
         DCONST1AT = DCONST1 / DGATK
         DCONST1AC = DCONST1 / DGACC
         DCONST1C = DCONST1 / DGCOR
         DCONST2  = GRAV / ( 18.0 * AMU )
         DCONST3AT = DCONST2 * PDENSAT * DGATK * DGATK
         DCONST3AC = DCONST2 * PDENSAC * DGACC * DGACC
         DCONST3C  = DCONST2 * PDENSCO * DGCOR * DGCOR

C      Aitken mode

         DCHAT0AT  = DCONST1AT
     &                    * ( ESAT04  + BHAT * KNATK * ESAT16 )
         DCHAT2AT  = DCONST1AT
     &                    * ( ESATM12  + BHAT * KNATK * ESATM16 )
         DCHAT3AT  = DCONST1AT
     &                    * ( ESATM20 + BHAT * KNATK * ESATM32 )
         VGHAT0AT  = DCONST3AT
     &                    * ( ESAT16  + BHAT * KNATK * ESAT04 )
         VGHAT2AT  = DCONST3AT
     &                    * ( ESAT48  + BHAT * KNATK * ESAT20 )
         VGHAT3AT  = DCONST3AT
     &                    * ( ESAT64  + BHAT * KNATK * ESAT28 )

C     accumulation mode

         DCHAT0AC  = DCONST1AC
     &                    * ( ESAC04  + BHAT * KNACC * ESAC16 )
         DCHAT2AC  = DCONST1AC
     &                    * ( ESACM12  + BHAT * KNACC * ESACM16 )
         DCHAT3AC  = DCONST1AC
     &                    * ( ESACM20 + BHAT * KNACC * ESACM32 )
         VGHAT0AC  = DCONST3AC
     &                    * ( ESAC16  + BHAT * KNACC * ESAC04 )
         VGHAT2AC  = DCONST3AC
     &                    * ( ESAC48  + BHAT * KNACC * ESAC20 )
         VGHAT3AC  = DCONST3AC
     &                    * ( ESAC64  + BHAT * KNACC * ESAC28 )

C     coarse mode

         DCHAT0C  = DCONST1C
     &                    * ( ESC04  + BHAT * KNCOR * ESC16 )
         DCHAT3C  = DCONST1C
     &                    * ( ESCM20 + BHAT * KNCOR * ESCM32 )
         VGHAT0C  = DCONST3C
     &                    * ( ESC16  + BHAT * KNCOR * ESC04 )
         VGHAT3C  = DCONST3C
     &                    * ( ESC64  + BHAT * KNCOR * ESC28 )

c now calculate the deposition  velocities

         NU = AMU / AIRDENS
         USTFAC = USTAR * USTAR / ( GRAV * NU )
         UTSCALE = USTAR
     &           + 0.24 * WSTAR * WSTAR
     &           /        USTAR

C first do 0th moment for the deposition of number number

C  Aitken mode

         SC0AT = NU / DCHAT0AT
         ST0AT = MAX( VGHAT0AT * USTFAC, 0.01 )
         RD0AT = 1.0 / ( UTSCALE *
     &            ( SC0AT ** ( -TWO3 ) + 10.0 ** ( -3.0 / ST0AT ) ) )

         VDEP( VDNATK ) = VGHAT0AT
     &                        + 1.0 / (
     &     RA + RD0AT + RD0AT * RA * VGHAT0AT
     &                                )

C accumulation mode

         SC0AC = NU / DCHAT0AC
         ST0AC = MAX ( VGHAT0AC * USTFAC, 0.01 )
         RD0AC = 1.0 / ( UTSCALE *
     &             ( SC0AC ** ( -TWO3 ) + 10.0 ** ( -3.0 / ST0AC ) ) )

         VDEP( VDNACC ) = VGHAT0AC
     &                        + 1.0 / (
     &     RA + RD0AC + RD0AC * RA * VGHAT0AC
     &                                )

c coarse mode

         SC0C = NU / DCHAT0C

         RD0C = 1.0 / ( UTSCALE *
     &                ( SC0C ** ( -TWO3 ) ) )
                         ! eliminate impaction term

         VDEP( VDNCOR ) = VGHAT0C
     &                         + 1.0 / (
     &     RA + RD0C + RD0C * RA * VGHAT0C
     &                                 )

C now do 2nd moment for the deposition of surface area

C  Aitken mode

         SC2AT = NU / DCHAT2AT
         ST2AT = MAX( VGHAT2AT * USTFAC, 0.01 )
         RD2AT = 1.0 / ( UTSCALE *
     &            ( SC2AT ** ( -TWO3 ) + 10.0 ** ( -3.0 / ST2AT ) ) )

         VDEP( VDSATK ) = VGHAT2AT + 1.0
     &                       / (
     &     RA + RD2AT + RD2AT * RA * VGHAT2AT
     &                         )

C accumulation mode

         SC2AC = NU / DCHAT2AC
         ST2AC = MAX( VGHAT2AC * USTFAC , 0.01 )
         RD2AC = 1.0 / ( UTSCALE *
     &           ( SC2AC ** ( -TWO3 ) + 10.0 ** ( -3.0 / ST2AC ) ) )


         VDEP( VDSACC ) =  VGHAT2AC + 1.0
     &        / ( RA + RD2AC +
     &          RD2AC * RA * VGHAT2AC
     &                                    )

C now do 3rd moment for the deposition of mass

C  Aitken mode

         SC3AT = NU / DCHAT3AT
         ST3AT = MAX( VGHAT3AT * USTFAC, 0.01 )
         RD3AT = 1.0 / ( UTSCALE *
     &             ( SC3AT ** ( -TWO3 ) + 10.0 ** ( -3.0 / ST3AT ) ) )

         VDEP( VDMATK ) = VGHAT3AT + 1.0
     &                       / (
     &     RA + RD3AT + RD3AT * RA * VGHAT3AT
     &                         )

C accumulation mode

         SC3AC = NU / DCHAT3AC
         ST3AC = MAX( VGHAT3AC * USTFAC , 0.01 )
         RD3AC = 1.0 / ( UTSCALE *
     &            ( SC3AC ** ( -TWO3 ) + 10.0 ** ( -3.0 / ST3AC ) ) )


         VDEP( VDMACC ) =  VGHAT3AC + 1.0
     &        / ( RA + RD3AC +
     &          RD3AC * RA * VGHAT3AC
     &                                           )

c coarse mode

         SC3C = NU / DCHAT3C

         RD3C = 1.0 / ( UTSCALE *
     &                ( SC3C ** ( -TWO3 ) ) )
                        ! eliminate impaction term

         VDEP( VDMCOR ) = VGHAT3C + 1.0
     &                        / (
     &     RA + RD3C + RD3C * RA * VGHAT3C
     &                          )

      RETURN

      END  SUBROUTINE GETVDEP
c /////////////////////////////////////////////////////////////////////
