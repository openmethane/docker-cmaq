
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/cloud/cloud_acm_ae5_txhgsim/scavwdep.F,v 1.1.1.1 2010/06/14 16:03:04 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE SCAVWDEP ( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                      CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                      REMOV, REMOVAC, ALFA0, ALFA2, ALFA3 )
C-----------------------------------------------------------------------
C
C  DESCRIPTION: Compute simplistic incloud scavenging and wet removal
C
C  Revision History:
C      No   Date   Who	What
C      -- -------- ---  -----------------------------------------
C       0 01/15/98 sjr  created program
C       1 07/27/99 sjr  added scavenging coefficient scaling factors
C       3 Dec 00   Jeff move CGRID_MAP into f90 module
C       4 Dec 02   sjr  revised calls to HLCONST
C
C  Called by:  RADMCLD and RESCLD
C
C  Calls the following subroutines:  GETALPHA
C
C  Calls the following functions:  HLCONST
C
C  ARGUMENTS    TYPE      I/O       DESCRIPTION
C  ---------   -------  ------  --------------------------------
C    JDATE     integer   input  current model julian date (yyyyddd)
C    JTIME     integer   input  current model time (hhmmss)
C    WTBAR      real     input  avg total water content (kg/m3)
C    WCBAR      real     input  avg liquid water content (kg/m3)
C    TBARC      real     input  avg cloud temperature (K)
C    PBARC      real     input  avg cloud pressure (Pa)
C    CTHK1      real     input  cloud thickness (m)
C    AIRM       real     input  total air mass (moles/m2) in cloudy air
C    PRATE1     real     input  precipitation rate (mm/hr)
C    TAUCLD     real     input  cloud lifetime (s)
C    POLC       real     input  ave vert conc incloud (moles sp/ mole air)
C    CEND       real    output  ending incloud conc (moles/mole)
C    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
C    REMOVAC    real    output  variable storing H+ deposition
C    ALFA0      real    output  scav coef for aitken aerosol number
C    ALFA3      real    output  scav coef for aitken aerosol mass
C
C-----------------------------------------------------------------------

      USE CGRID_SPCS   ! CGRID species number and offsets

      IMPLICIT NONE

C...........INCLUDES and their descriptions

      INCLUDE SUBST_GC_SPC              ! gas chemistry species table
      INCLUDE SUBST_AE_SPC              ! aerosol species table
      INCLUDE SUBST_NR_SPC              ! non-reactive species table
      INCLUDE SUBST_TR_SPC              ! tracer species table
      INCLUDE SUBST_GC_SCAV             ! gas scavenging maping table
      INCLUDE SUBST_AE_SCAV             ! aerosol scavenging maping table
      INCLUDE SUBST_NR_SCAV             ! non-reactive scavenging maping table
      INCLUDE SUBST_TR_SCAV             ! tracer scavenging maping table

      INCLUDE SUBST_CONST               ! constants

      CHARACTER*120 XMSG                ! Exit status message
      DATA          XMSG / ' ' /

C...........PARAMETERS and their descriptions:

      REAL          H2ODENS             ! density of water at 20 C and 1 ATM
      PARAMETER   ( H2ODENS = 1000.0 )  ! (kg/m3)

      INTEGER       MXSPCS              ! Number of species in CGRID
      PARAMETER   ( MXSPCS = N_GC_SPCD
     &                     + N_AE_SPC
     &                     + N_NR_SPC
     &                     + N_TR_SPC )

      INTEGER       N_CGRID_SCAVD       ! # of species scavenged
      PARAMETER   ( N_CGRID_SCAVD = N_GC_SCAVD
     &                            + N_AE_SCAV
     &                            + N_NR_SCAV
     &                            + N_TR_SCAV )

      INTEGER       N_CGRID_SCAV        ! # of species scavenged
      PARAMETER   ( N_CGRID_SCAV = N_GC_SCAV
     &                           + N_AE_SCAV
     &                           + N_NR_SCAV
     &                           + N_TR_SCAV )

      REAL         TWOTHIRDS            ! 2/3
      PARAMETER  ( TWOTHIRDS = 2.0 / 3.0 )

      REAL         KGPG             ! kilograms per gram
      PARAMETER  ( KGPG = 1.0E-03 )

C...........ARGUMENTS and their descriptions

      INTEGER       JDATE               ! current model date, coded YYYYDDD
      INTEGER       JTIME               ! current model time, coded HHMMSS
      REAL          WTBAR               ! total wat cont (kg/m2) int. thru cloud depth
      REAL          WCBAR               ! liq water content of cloud (kg/m3)
      REAL          TBARC               ! mean cloud temp (K)
      REAL          PBARC               ! mean cloud pressure (Pa)
      REAL          CTHK1               ! aq chem calc cloud thickness
      REAL          AIRM                ! total air mass (moles/m2) in cloudy air
      REAL          PRATE1              ! storm rainfall rate (mm/hr)
      REAL          TAUCLD              ! cloud lifetime
      REAL          ALFA0               ! scav coef for aitken aerosol number
      REAL          ALFA2               ! scav coef for aitken aerosol sfc area
      REAL          ALFA3               ! scav coef for aitken aerosol mass
      REAL          REMOVAC             ! variable storing H+ deposition
      REAL          POLC   ( MXSPCS )   ! ave conc incloud (moles/mole)
      REAL          CEND   ( MXSPCS )   ! ending conc (moles/mole)
      REAL          REMOV  ( MXSPCS )   ! moles/m2 or mm*mol/lit scavenged

C...........SCRATCH LOCAL VARIABLES and their descriptions:

      LOGICAL       FIRSTIME            ! flag for first pass thru
      DATA          FIRSTIME / .TRUE. /
      SAVE          FIRSTIME

      CHARACTER*16  PNAME               ! program name
      DATA          PNAME / 'SCAVWDEP' /
      SAVE          PNAME
      CHARACTER*16  CGRID_SCAV( N_CGRID_SCAVD )  ! CGRID species scavenged
      SAVE          CGRID_SCAV

      INTEGER       I
      INTEGER       N_NUMAKN            ! # aitken aerosol number species
      SAVE          N_NUMAKN
      INTEGER       N_MASAKN            ! # aitken aerosol mass species
      SAVE          N_MASAKN
      INTEGER       N_SRFAKN            ! # aitken aerosol sfc area species
      SAVE          N_SRFAKN
      INTEGER       PNTR                ! relative pointer variable
      INTEGER       SPC                 ! liquid species loop counter
      INTEGER       VAR                 ! variable loop counter
      INTEGER       CGRID_SCAV_MAP( N_CGRID_SCAVD ) ! CGRID map to scavenged spc
      SAVE          CGRID_SCAV_MAP
      INTEGER       L_NUMAKN( MXSPCS )  ! pointers to aitken aerosol #
      SAVE          L_NUMAKN
      INTEGER       L_MASAKN( MXSPCS )  ! pointers to aitken aerosols
      SAVE          L_MASAKN
      INTEGER       L_SRFAKN( MXSPCS )  ! pntrs to aitken aerosol surface area
      SAVE          L_SRFAKN

      REAL          ALFA                ! scavenging coefficient (1/s)
      REAL          KH                  ! Henry's law constant (mol/l/atm)
      REAL          NUMAKN              ! Aitken mode aerosol # (#/m3)
      REAL          MASAKN              ! Total Aitken mode mass (ug/m3)
      REAL          SRFAKN
      REAL          ONE_OVER_TWASH      ! 1 / TWASH
      REAL          HPLUS               ! hydrogen ion concentration (mol/l)
      DATA          HPLUS / 1.0E-4 /    ! typical value
      SAVE          HPLUS
      REAL          RHOAIR              ! air density in kg/m3
      REAL          RTCH                ! chemical gas const times temp
      REAL          TWASH               ! washout time for clouds (sec) with low liq wat content
      REAL          TWF                 ! washout scaling factor (mol/l/atm)
      REAL          CGRID_SCAV_FAC( N_CGRID_SCAVD )  ! CGRID scav coef factors
      SAVE          CGRID_SCAV_FAC

C...........EXTERNAL FUNCTIONS and their descriptions:

      INTEGER       INDEXN
      REAL          HLCONST
      EXTERNAL      INDEXN, HLCONST

C-----------------------------------------------------------------------
C  begin body of subroutine SCAVWDEP

C...INITIALIZATION SCAVWDEP module:
C...  event-statistics variables.  Open output files.

      IF ( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

C...first check to make sure that some species in CGRID were specified
C...  for scavenging, otherwise notify the user and return

        IF ( N_CGRID_SCAV .LE. 0 ) THEN
          XMSG = 'No species were specified for scavenging by cloud ' //
     &           'or rain water...SCAVENGING WILL NOT BE PERFORMED!'
          CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
          RETURN
        END IF

        CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C...prepare indices for scavenged species

C...  load the CGRID to scavenged species pointers for the gases

        SPC = 0
        DO VAR=1, N_GC_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = GC_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = GC_SCAV_MAP( VAR ) + GC_STRT - 1
          CGRID_SCAV_FAC( SPC ) = GC_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the aerosols

        DO VAR=1, N_AE_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = AE_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = AE_SCAV_MAP( VAR ) + AE_STRT - 1
          CGRID_SCAV_FAC( SPC ) = AE_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the non-reactives

        DO VAR=1, N_NR_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = NR_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = NR_SCAV_MAP( VAR ) + NR_STRT - 1
          CGRID_SCAV_FAC( SPC ) = NR_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the tracers

        DO VAR=1, N_TR_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = TR_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = TR_SCAV_MAP( VAR ) + TR_STRT - 1
          CGRID_SCAV_FAC( SPC ) = TR_SCAV_FAC( VAR )
        END DO

C...create the pointers from CGRID to the species needed by AQCHEM

        N_NUMAKN = INDEXN( 'NUM_AITKEN      ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_NUMAKN )
        N_MASAKN = INDEXN( 'AITKEN', N_CGRID_SCAV, CGRID_SCAV, L_MASAKN )
        N_SRFAKN = INDEXN( 'SRF_AITKEN      ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_SRFAKN )

      END IF

C...for subsequent calls, check to make sure some species were
C...  specified, otherwise there is no need to perform scavenging

      IF ( N_CGRID_SCAV .LE. 0 ) THEN
        RETURN
      END IF

      RTCH = ( MOLVOL / STDTEMP ) * TBARC
      TWASH = WTBAR * 1000.0 * CTHK1 * 3600.0
     &      / ( H2ODENS * AMAX1( 1.0E-20, PRATE1 ) )
ccc          TWASH = AMAX1( TWASH, TAUCLD / 60.0 )   ! <different units?sec&min
      TWASH = AMAX1( TWASH, TAUCLD )
      ONE_OVER_TWASH = 1.0 / TWASH
      TWF = H2ODENS / ( WTBAR * RTCH )
      REMOVAC = 0.0

      RHOAIR = ( AIRM / CTHK1 ) * MWAIR * KGPG

C...compute total Aitken mode number (#/m3)

      NUMAKN = 0.0

      DO I = 1, N_NUMAKN
        PNTR = CGRID_SCAV_MAP( L_NUMAKN( I ) )
        NUMAKN = NUMAKN + ( POLC( PNTR ) * AIRM / CTHK1 )
      END DO

C...compute total Aitken mode mass (ug/m3)

      MASAKN = 0.0

      DO I = 1, N_MASAKN
        PNTR = CGRID_SCAV_MAP( L_MASAKN( I ) )
        IF (( INDEX( CGRID_SCAV( L_MASAKN( I ) ), 'NUM' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASAKN( I ) ), 'SRF' ) .EQ. 0 ) .AND.
     &      ( INDEX( CGRID_SCAV( L_MASAKN( I ) ), 'H2O' ) .EQ. 0 ) ) THEN
          MASAKN = MASAKN + ( POLC( PNTR ) * AIRM / CTHK1
     &           * AE_MOLWT( PNTR - AE_STRT + 1 ) * 1.0E6 )
        END IF
      END DO

C...compute total Aitken mode surface area (m2/m3)

      SRFAKN = 0.0

      DO I = 1, N_SRFAKN
        PNTR = CGRID_SCAV_MAP( L_SRFAKN( I ) )
        SRFAKN = SRFAKN + ( POLC( PNTR ) * AIRM / CTHK1 )
      END DO

C...compute the scavenging coefficients for aitken mode aerosol mass
C...  and number
C...  NOTE:  for now, scavenging coefficients are computed for only
C...         the liquid water content, not on the total water content
C...         therefore, no ice phase scavenging is considered at this
C...         time, but it should be added in the future!

      CALL GETALPHA ( NUMAKN, MASAKN, SRFAKN, WCBAR, TBARC, PBARC,
     &                RHOAIR, ALFA0, ALFA2, ALFA3 )

C...gas scavenging and wet deposition

      SPC = 0

      DO VAR = 1, N_GC_SCAV
        SPC = SPC + 1
        PNTR = CGRID_SCAV_MAP( SPC )
        KH = HLCONST( CGRID_SCAV( SPC ), TBARC, .TRUE., HPLUS )
        IF ( KH .GT. 0.0 ) THEN
          ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH / ( 1.0 + TWF / KH )
          CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
          REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
        ELSE
          ALFA = 0.0
          CEND ( PNTR ) = POLC( PNTR )
          REMOV( PNTR ) = 0.0
        END IF
      END DO

C...aerosol scavenging and wet deposition

      DO VAR = 1, N_AE_SCAV
        SPC = SPC + 1

        PNTR = CGRID_SCAV_MAP( SPC )

        IF ( INDEX( CGRID_SCAV( SPC ), 'AITKEN' ) .GT. 0 ) THEN
          IF ( INDEX( CGRID_SCAV( SPC ), 'NUM' ) .GT. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA0
          ELSE IF ( INDEX( CGRID_SCAV( SPC ), 'SRF' ) .GT. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA2
          ELSE
            ALFA = CGRID_SCAV_FAC( SPC ) * ALFA3
          END IF
        ELSE
csjr          IF ( INDEX( CGRID_SCAV( SPC ), 'SRF' ) .EQ. 0 ) THEN
            ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH
csjr          ELSE
csjr            ALFA = CGRID_SCAV_FAC( SPC ) * TWOTHIRDS * ONE_OVER_TWASH
csjr          END IF
        END IF

        CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
        REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM

      END DO

C...non-reactive scavenging and wet deposition

      DO VAR = 1, N_NR_SCAV
        SPC = SPC + 1
        PNTR = CGRID_SCAV_MAP( SPC )
        KH = HLCONST( CGRID_SCAV( SPC ), TBARC, .TRUE., HPLUS )
        IF ( KH .GT. 0.0 ) THEN
          ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH / ( 1.0 + TWF / KH )
          CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
          REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
        ELSE
          ALFA = 0.0
          CEND ( PNTR ) = POLC( PNTR )
          REMOV( PNTR ) = 0.0
        END IF
      END DO

C...tracer scavenging and wet deposition

      DO VAR = 1, N_TR_SCAV
        SPC = SPC + 1
        PNTR = CGRID_SCAV_MAP( SPC )
        KH = HLCONST( CGRID_SCAV( SPC ), TBARC, .TRUE., HPLUS )
        IF ( KH .GT. 0.0 ) THEN
          ALFA = CGRID_SCAV_FAC( SPC ) * ONE_OVER_TWASH / ( 1.0 + TWF / KH )
          CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
          REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
        ELSE
          ALFA = 0.0
          CEND ( PNTR ) = POLC( PNTR )
          REMOV( PNTR ) = 0.0
        END IF
      END DO

      RETURN
      END
