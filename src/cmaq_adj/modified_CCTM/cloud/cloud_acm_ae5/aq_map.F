
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/cloud/cloud_acm_ae5/aq_map.F,v 1.1.1.1 2010/06/14 16:03:03 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE AQ_MAP ( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                    CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                    REMOV, REMOVAC, ALFA0, ALFA2, ALFA3 )

C-----------------------------------------------------------------------
C
C  DESCRIPTION:  This is an interface processor between the cloud dynamics
C     module(s) and the aqueous chemistry module.  It uses indices to
C     map species from their order in the array POLC, CEND, REMOV to
C     the species that are expected for AQCHEM, ie GAS and AEROSOL.
C
C  Revision History:
C      No   Date   Who	What
C      -- -------- ---  -----------------------------------------
C       0 01/15/98 sjr  created program
C       1 02/13/98 sjr  modified/generalized the mapping procedure
C         Dec 00   Jeff move CGRID_MAP into f90 module
C       3 06/07/05 sjr  added logic for coarse sulfate
C       4 04/11/08 jtk  added logic for coarse ammonium
C
C  Called by:  RADMCLD and RESCLD
C
C  Calls the following subroutines:  AQCHEM
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
C    ALFA0      real     input  scav coef for aitken aerosol number
C    ALFA2      real     input  scav coef for aitken aerosol surface area
C    ALFA3      real     input  scav coef for aitken aerosol mass
C
C-----------------------------------------------------------------------

      USE CGRID_SPCS   ! CGRID species number and offsets

      IMPLICIT NONE

      INCLUDE SUBST_GC_SPC              ! gas chemistry species table
      INCLUDE SUBST_AE_SPC              ! aerosol species table
      INCLUDE SUBST_NR_SPC              ! non-reactive species table
      INCLUDE SUBST_TR_SPC              ! tracer species table
      INCLUDE SUBST_GC_G2AQ             ! gas to aqueous maping table
      INCLUDE SUBST_AE_A2AQ             ! aerosol to aqueous maping table
      INCLUDE SUBST_NR_N2AQ             ! non-reactive to aqueous maping table
      INCLUDE SUBST_TR_T2AQ             ! tracer to aqueous maping table


      INCLUDE 'AQ_PARAMS.EXT'           ! aqueous chemistry shared parameters

      CHARACTER*120 XMSG                ! Exit status message
      DATA          XMSG / ' ' /

C...........PARAMETERS and their descriptions:

      INTEGER       MXSPCS              ! Number of species in CGRID
      PARAMETER   ( MXSPCS = N_GC_SPCD
     &                     + N_AE_SPC
     &                     + N_NR_SPC
     &                     + N_TR_SPC )

      INTEGER       N_CGRID2AQD         ! # of species used in Aqeuous module
      PARAMETER   ( N_CGRID2AQD = N_GC_G2AQD
     &                          + N_AE_A2AQ
     &                          + N_NR_N2AQ
     &                          + N_TR_T2AQ )

      INTEGER       N_CGRID2AQ          ! # of species used in Aqeuous module
      PARAMETER   ( N_CGRID2AQ = N_GC_G2AQ
     &                         + N_AE_A2AQ
     &                         + N_NR_N2AQ
     &                         + N_TR_T2AQ )

C...........ARGUMENTS and their descriptions

      INTEGER       JDATE               ! current model date, coded YYYYDDD
      INTEGER       JTIME               ! current model time, coded HHMMSS

      REAL          AIRM                ! total air mass (moles/m2) in cloudy air
      REAL          ALFA0               ! scav coef for aitken aerosol number
      REAL          ALFA2               ! scav coef for aitken aerosol sfc area
      REAL          ALFA3               ! scav coef for aitken aerosol mass
      REAL          CTHK1               ! aq chem calc cloud thickness
      REAL          PBARC               ! mean cloud pressure (Pa)
      REAL          PRATE1              ! storm rainfall rate (mm/hr)
      REAL          REMOVAC             ! variable storing H+ deposition
      REAL          TAUCLD              ! cloud lifetime
      REAL          TBARC               ! mean cloud temp (K)
      REAL          WCBAR               ! liq water content of cloud (kg/m3)
      REAL          WTBAR               ! total wat cont (kg/m2) int. thru cld

      REAL          POLC ( MXSPCS )     ! ave vert conc incloud (moles/mole)
      REAL          CEND ( MXSPCS )     ! ending incloud conc (moles/mole)
      REAL          REMOV( MXSPCS )     ! moles/m2 or mm*mol/lit scavenged

C...........LOCAL VARIABLES (scalars) and their descriptions:

      LOGICAL       FIRSTIME            ! flag for first pass thru
      DATA          FIRSTIME / .TRUE. /
      SAVE          FIRSTIME

      CHARACTER*16  PNAME               ! program name
      DATA          PNAME / 'AQ_MAP' /
      SAVE          PNAME
      CHARACTER*16  CGRID2AQ( N_CGRID2AQD )  ! CGRID species used in AQCHEM
      SAVE          CGRID2AQ
      CHARACTER*16  BKGROUND            ! temporary background var string
      CHARACTER*16  BKUNITS             ! temporary background units string

      INTEGER       NSGRGAS( NGAS )     ! # surrogates in CGRID for each gas
      SAVE          NSGRGAS
      INTEGER       NSGRAER( NAER )     ! # surrogates in CGRID for each aerosol
      SAVE          NSGRAER
      INTEGER       IAER                ! aerosol loop counter
      INTEGER       IGAS                ! gas loop counter
      INTEGER       ISGR                ! surrogate loop counter
      INTEGER       PNTR                ! relative pointer variable
      INTEGER       SPC                 ! liquid species loop counter
      INTEGER       VAR                 ! variable loop counter
      INTEGER       LSGR( MXSPCS )      ! temporary pointer array
      INTEGER       LSGRGAS( NGAS, MXSPCS ) ! pointers in CGRID to each gas
      SAVE          LSGRGAS
      INTEGER       LSGRAER( NAER, MXSPCS ) ! pointers in CGRID to each aerosol
      SAVE          LSGRAER
      INTEGER       CGRID2AQ_MAP( N_CGRID2AQD ) ! CGRID map to aqueous species
      SAVE          CGRID2AQ_MAP

      REAL          HPWDEP              ! hydrogen wet dep (mm mol/liter)
      REAL          GAS    ( NGAS )     ! gas phase conc (mol/mol)
      REAL          GASWDEP( NGAS )     ! gas phase wet dep array (mm mol/liter)
      REAL          AEROSOL( NAER )     ! aerosol conc (mol/mol)
      REAL          AERWDEP( NAER )     ! aerosol wet dep array (mm mol/liter)
      REAL          WSGRGAS( NGAS, MXSPCS )  ! weights for surrogate
      REAL          WSGRAER( NAER, MXSPCS )  ! weights for surrogate

C...........EXTERNAL FUNCTIONS and their descriptions:

      INTEGER       INDEXN              ! external func to get species pointers
      INTEGER       TRIMLEN             !  string length, excl. trailing blanks

      EXTERNAL      INDEXN, TRIMLEN

C*********************************************************************
C     begin body of subroutine AQ_MAP

C...INITIALIZATION
C...  event-statistics variables.

      IF ( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

C...first check to make sure that some surrogates were specified for the
C...  aqueous species, otherwise notify the user and return

        IF ( N_CGRID2AQ .LE. 0 ) THEN
          XMSG = 'No surrogates specified for any aqueous species' //
     &           '...aqueous chemistry WILL NOT BE PERFORMED!'
          CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
          RETURN
        END IF

        CALL CGRID_MAP ( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C...prepare indices for aqueous species
C...  load the CGRID to aqueous species pointers for the gases

        SPC = 0
        DO VAR=1, N_GC_G2AQ
          SPC = SPC + 1
          CGRID2AQ( SPC ) = GC_G2AQ( VAR )
          CGRID2AQ_MAP( SPC ) = GC_G2AQ_MAP( VAR ) + GC_STRT - 1
        END DO

C...  load the CGRID to aqueous species pointers for the aerosols

        DO VAR=1, N_AE_A2AQ
          SPC = SPC + 1
          CGRID2AQ( SPC ) = AE_A2AQ( VAR )
          CGRID2AQ_MAP( SPC ) = AE_A2AQ_MAP( VAR ) + AE_STRT - 1
        END DO

C...  load the CGRID to aqueous species pointers for the non-reactives

        DO VAR=1, N_NR_N2AQ
          SPC = SPC + 1
          CGRID2AQ( SPC ) = NR_N2AQ( VAR )
          CGRID2AQ_MAP( SPC ) = NR_N2AQ_MAP( VAR ) + NR_STRT - 1
        END DO

C...  load the CGRID to aqueous species pointers for the tracers

        DO VAR=1, N_TR_T2AQ
          SPC = SPC + 1
          CGRID2AQ( SPC ) = TR_T2AQ( VAR )
          CGRID2AQ_MAP( SPC ) = TR_T2AQ_MAP( VAR ) + TR_STRT - 1
        END DO

C...create the pointers from CGRID to the gaseous species needed by AQCHEM

        DO IGAS=1, NGAS
          NSGRGAS( IGAS ) = INDEXN( SGRGAS( IGAS ), N_CGRID2AQ,
     &                              CGRID2AQ, LSGR )
          IF ( NSGRGAS( IGAS ) .EQ. 0 ) THEN
            WRITE ( BKGROUND, '(1PE8.2)' ) BGNDGAS( IGAS )
            BKUNITS = BUNTSGAS( IGAS )
            XMSG = 'No surrogates specified for aqueous species ' //
     &              SGRGAS( IGAS )( 1:TRIMLEN( SGRGAS( IGAS ) ) ) //
     &             '...Using background value of ' //
     &              BKGROUND( 1:TRIMLEN( BKGROUND ) ) // ' ' //
     &              BKUNITS ( 1:TRIMLEN( BKUNITS  ) )
            CALL M3MESG( XMSG )
          ELSE
            DO ISGR = 1, NSGRGAS( IGAS )
              LSGRGAS( IGAS, ISGR ) = LSGR( ISGR )
            END DO
          END IF
        END DO

C...create the pointers from CGRID to the aerosol species needed by AQCHEM

        DO IAER=1, NAER
          NSGRAER( IAER ) = INDEXN( SGRAER( IAER ), N_CGRID2AQ,
     &                              CGRID2AQ, LSGR )
          IF ( NSGRAER( IAER ) .EQ. 0 ) THEN
            WRITE ( BKGROUND, '(1PE8.2)' ) BGNDAER( IAER )
            BKUNITS = BUNTSAER( IAER )
            XMSG = 'No surrogates specified for aqueous species ' //
     &              SGRAER( IAER )( 1:TRIMLEN( SGRAER( IAER ) ) ) //
     &             '...Using background value of ' //
     &              BKGROUND( 1:TRIMLEN( BKGROUND ) ) // ' ' //
     &              BKUNITS ( 1:TRIMLEN( BKUNITS  ) )
            CALL M3MESG( XMSG )
          ELSE
            DO ISGR = 1, NSGRAER( IAER )
              LSGRAER( IAER, ISGR ) = LSGR( ISGR )
            END DO
          END IF
        END DO

      END IF

C...for subsequent calls, check to make sure some surrogates were
C...  specified, otherwise there is no need to perform aqueous chemistry

      IF ( N_CGRID2AQ .LE. 0 ) THEN
        RETURN
      END IF

C...load gas-phase concentrations

      DO IGAS = 1, NGAS

        GAS( IGAS ) = 0.0

        DO ISGR = 1, NSGRGAS( IGAS )
          PNTR = CGRID2AQ_MAP( LSGRGAS( IGAS, ISGR ) )
          GAS( IGAS ) = GAS( IGAS ) + POLC( PNTR )
        END DO

        IF ( GAS( IGAS ) .GT. 0.0 ) THEN
          DO ISGR = 1, NSGRGAS( IGAS )
            PNTR = CGRID2AQ_MAP( LSGRGAS( IGAS, ISGR ) )
            WSGRGAS( IGAS, ISGR ) = POLC( PNTR ) / GAS( IGAS )
          END DO
        ELSE
          DO ISGR = 1, NSGRGAS( IGAS )
            WSGRGAS( IGAS, ISGR ) = 1.0 / FLOAT( NSGRGAS( IGAS ) )
          END DO
        END IF

C...set background values for gases if no surrogates were specified

        IF ( NSGRGAS( IGAS ) .EQ. 0 ) THEN
          GAS( IGAS ) = BGNDGAS( IGAS ) * 1.0E-6
        END IF

      END DO

C...load aerosol concentrations

      DO IAER = 1, NAER

        AEROSOL( IAER ) = 0.0

        DO ISGR = 1, NSGRAER( IAER )
          PNTR = CGRID2AQ_MAP( LSGRAER( IAER, ISGR ) )
          AEROSOL( IAER ) = AEROSOL( IAER ) + POLC( PNTR )
        END DO

        IF ( AEROSOL( IAER ) .GT. 0.0 ) THEN
          DO ISGR = 1, NSGRAER( IAER )
            PNTR = CGRID2AQ_MAP( LSGRAER( IAER, ISGR ) )
            WSGRAER( IAER, ISGR ) = POLC( PNTR ) / AEROSOL( IAER )
          END DO
        ELSE
          DO ISGR = 1, NSGRAER( IAER )
            WSGRAER( IAER, ISGR ) = 1.0 / FLOAT( NSGRAER( IAER ) )
          END DO
        END IF

C...set background values for aerosols if no surrogates were specified

        IF ( NSGRAER( IAER ) .EQ. 0 ) THEN
          IF ( INDEX( BUNTSAER( IAER ), 'ug/m3' ) .GT. 0 ) THEN
            AEROSOL( IAER ) = BGNDAER( IAER ) * 1.0E-6 * CTHK1
     &                      / ( SGRAERMW( IAER ) * AIRM )
          ELSE
            AEROSOL( IAER ) = BGNDAER( IAER ) * CTHK1 / AIRM
          END IF
        END IF

      END DO

C...now perform aqueous-phase chemistry calculations

      CALL AQCHEM ( JDATE, JTIME, TBARC,  PBARC, TAUCLD, PRATE1,
     &              WCBAR, WTBAR, AIRM, ALFA0, ALFA2, ALFA3, GAS,
     &              AEROSOL, GASWDEP, AERWDEP, HPWDEP )

C...store the amount of hydrogen deposition

      REMOVAC = HPWDEP

C...since the aqueous chemistry routine doesn't know which surrogate
C...  species were specified by the user, several checks must be made
C...  on the surrogate species to make sure that all of the mass is
C...  conserved and doesn't fall through the cracks.

      DO IGAS = 1, NGAS

C...check for the special case for species H2SO4

        IF ( IGAS .EQ. LH2SO4 ) THEN

C...  load sulfate accumulation aerosol mass into the CEND slot for
C...    sulfuric acid gas if there is no slot for it in CEND

          IF ( NSGRAER( LSO4ACC ) .EQ. 0 ) THEN
            GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LSO4ACC )
            GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LSO4ACC )
            AEROSOL( LSO4ACC ) = 0.0
            AERWDEP( LSO4ACC ) = 0.0

            IF ( NSGRAER( LSO4COR ) .EQ. 0 ) THEN

              GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LSO4COR )
              GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LSO4COR )
              AEROSOL( LSO4COR ) = 0.0
              AERWDEP( LSO4COR ) = 0.0

              IF ( NSGRAER( LSO4AKN ) .EQ. 0 ) THEN
                GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LSO4AKN )
                GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LSO4AKN )
                AEROSOL( LSO4AKN ) = 0.0
                AERWDEP( LSO4AKN ) = 0.0
              END IF

            ELSE

              IF ( NSGRAER( LSO4AKN ) .EQ. 0 ) THEN
                GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LSO4AKN )
                GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LSO4AKN )
                AEROSOL( LSO4AKN ) = 0.0
                AERWDEP( LSO4AKN ) = 0.0
              END IF

            END IF
          END IF

C...check for the special case for species NH3

        ELSE IF ( IGAS .EQ. LNH3 ) THEN

C...  load ammonium accumulation aerosol mass into the CEND slot for
C...   ammonia gas if there is no slot for it in CEND

          IF ( NSGRAER( LNH4ACC ) .EQ. 0 ) THEN
            GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNH4ACC )
            GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNH4ACC )
            AEROSOL( LNH4ACC ) = 0.0
            AERWDEP( LNH4ACC ) = 0.0

C...  load ammonium coarse aerosol mass into the CEND slot for
C...   ammonia gas if there is no slot for it in CEND

            IF ( NSGRAER( LNH4COR ) .EQ. 0 ) THEN 

              GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNH4COR )
              GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNH4COR )
              AEROSOL( LNH4COR ) = 0.0
              AERWDEP( LNH4COR ) = 0.0

C...  load ammonium aitken aerosol mass into the CEND slot for
C...   ammonia gas if there is no slot for it in CEND

              IF ( NSGRAER( LNH4AKN ) .EQ. 0 ) THEN
               GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNH4AKN )
               GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNH4AKN )
               AEROSOL( LNH4AKN ) = 0.0
               AERWDEP( LNH4AKN ) = 0.0
              END IF

            ELSE 

              IF ( NSGRAER( LNH4AKN ) .EQ. 0 ) THEN
               GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNH4AKN )
               GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNH4AKN )
               AEROSOL( LNH4AKN ) = 0.0
               AERWDEP( LNH4AKN ) = 0.0
              END IF

            END IF 
          END IF

C...check for the special case for species HNO3

        ELSE IF ( IGAS .EQ. LHNO3 ) THEN

C...  add-on any NITRATE mass output by AQCHEM which has no slot in array CEND

          IF ( NSGRAER( LNO3ACC ) .EQ. 0 ) THEN
            GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNO3ACC )
            GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNO3ACC )
            AEROSOL( LNO3ACC ) = 0.0
            AERWDEP( LNO3ACC ) = 0.0

            IF ( NSGRAER( LNO3COR ) .EQ. 0 ) THEN

              GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNO3COR )
              GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNO3COR )
              AEROSOL( LNO3COR ) = 0.0
              AERWDEP( LNO3COR ) = 0.0

              IF ( NSGRAER( LNO3AKN ) .EQ. 0 ) THEN
                GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNO3AKN )
                GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNO3AKN )
                AEROSOL( LNO3AKN ) = 0.0
                AERWDEP( LNO3AKN ) = 0.0
              END IF

            ELSE

              IF ( NSGRAER( LNO3AKN ) .EQ. 0 ) THEN
                GAS    ( IGAS    ) = GAS    ( IGAS ) + AEROSOL( LNO3AKN )
                GASWDEP( IGAS    ) = GASWDEP( IGAS ) + AERWDEP( LNO3AKN )
                AEROSOL( LNO3AKN ) = 0.0
                AERWDEP( LNO3AKN ) = 0.0
              END IF

            END IF
          END IF

        END IF
      END DO

      DO IAER = 1, NAER

        IF ( IAER .EQ. LSO4AKN ) THEN

C...  load sulfuric acid gas and sulfate accumulation aerosol mass into
C...    the CEND slot for sulfate aitken aerosol if there are no slots
C...    for either of those species in CEND

          IF ( ( NSGRGAS( LH2SO4  ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LSO4ACC ) .EQ. 0 ) ) THEN
            AEROSOL( IAER ) = AEROSOL( IAER    )
     &                      + GAS    ( LH2SO4  )
     &                      + AEROSOL( LSO4ACC )
            AERWDEP( IAER ) = AERWDEP( IAER    )
     &                      + GASWDEP( LH2SO4  )
     &                      + AERWDEP( LSO4ACC )
            GAS    ( LH2SO4  ) = 0.0
            GASWDEP( LH2SO4  ) = 0.0
            AEROSOL( LSO4ACC ) = 0.0
            AERWDEP( LSO4ACC ) = 0.0

            IF ( NSGRAER( LSO4COR ) .EQ. 0 ) THEN
              AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LSO4COR )
              AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LSO4COR )
              AEROSOL( LSO4COR ) = 0.0
              AERWDEP( LSO4COR ) = 0.0
            END IF

          END IF

        ELSE IF ( IAER .EQ. LSO4ACC ) THEN

C...  load sulfuric acid gas mass into the CEND slot for sulfate accumulation
C...    aerosol if there is no slot for it in CEND

          IF ( NSGRGAS( LH2SO4 ) .EQ. 0 ) THEN
            AEROSOL( IAER   ) = AEROSOL( IAER ) + GAS    ( LH2SO4 )
            AERWDEP( IAER   ) = AERWDEP( IAER ) + GASWDEP( LH2SO4 )
            GAS    ( LH2SO4 ) = 0.0
            GASWDEP( LH2SO4 ) = 0.0
          END IF

C...  load sulfate aitken aerosol mass into the CEND slot for
C...   sulfate accumulation aerosol if there is no slot for it in CEND

          IF ( NSGRAER( LSO4AKN ) .EQ. 0 ) THEN
            AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LSO4AKN )
            AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LSO4AKN )
            AEROSOL( LSO4AKN ) = 0.0
            AERWDEP( LSO4AKN ) = 0.0
          END IF

C...  load sulfate coarse aerosol mass into the CEND slot for
C...   sulfate accumulation aerosol if there is no slot for it in CEND

          IF ( NSGRAER( LSO4COR ) .EQ. 0 ) THEN
            AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LSO4COR )
            AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LSO4COR )
            AEROSOL( LSO4COR ) = 0.0
            AERWDEP( LSO4COR ) = 0.0
          END IF

        ELSE IF ( IAER .EQ. LSO4COR ) THEN

C...  check for H2SO4, SO4ACC, and SO4AKN in cgrid...and as
C...    a last resort, put all of the sulfate mass into SO4COR

          IF ( ( NSGRGAS( LH2SO4  ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LSO4ACC ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LSO4AKN ) .EQ. 0 ) ) THEN
            AEROSOL( IAER ) = AEROSOL( IAER    )
     &                      + GAS    ( LH2SO4  )
     &                      + AEROSOL( LSO4AKN )
     &                      + AEROSOL( LSO4ACC )
            AERWDEP( IAER ) = AERWDEP( IAER    )
     &                      + GASWDEP( LH2SO4  )
     &                      + AERWDEP( LSO4AKN )
     &                      + AERWDEP( LSO4ACC )
            GAS    ( LH2SO4  ) = 0.0
            GASWDEP( LH2SO4  ) = 0.0
            AEROSOL( LSO4AKN ) = 0.0
            AERWDEP( LSO4AKN ) = 0.0
            AEROSOL( LSO4ACC ) = 0.0
            AERWDEP( LSO4ACC ) = 0.0
          END IF

        ELSE IF ( IAER .EQ. LNH4AKN ) THEN

C...  load ammonia gas and ammonium accumulation aerosol mass into
C...    the CEND slot for ammonium aitken aerosol if there are no slots
C...    for either of those species in CEND

          IF ( ( NSGRGAS( LNH3    ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LNH4ACC ) .EQ. 0 ) ) THEN
            AEROSOL( IAER ) = AEROSOL( IAER    )
     &                      + GAS    ( LNH3    )
     &                      + AEROSOL( LNH4ACC )
            AERWDEP( IAER ) = AERWDEP( IAER    )
     &                      + GASWDEP( LNH3    )
     &                      + AERWDEP( LNH4ACC )
            GAS    ( LNH3    ) = 0.0
            GASWDEP( LNH3    ) = 0.0
            AEROSOL( LNH4ACC ) = 0.0
            AERWDEP( LNH4ACC ) = 0.0

            IF ( NSGRAER( LNH4COR ) .EQ. 0 ) THEN 
              AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LNH4COR )
              AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LNH4COR )
              AEROSOL( LNH4COR ) = 0.0
              AERWDEP( LNH4COR ) = 0.0
            END IF

          END IF

        ELSE IF ( IAER .EQ. LNH4ACC ) THEN

C...  load ammonia gas mass into the CEND slot for ammonium accumulation
C...    aerosol if there is no slot for it in CEND

          IF ( NSGRGAS( LNH3 ) .EQ. 0 ) THEN
            AEROSOL( IAER ) = AEROSOL( IAER ) + GAS    ( LNH3 )
            AERWDEP( IAER ) = AERWDEP( IAER ) + GASWDEP( LNH3 )
            GAS    ( LNH3 ) = 0.0
            GASWDEP( LNH3 ) = 0.0
          END IF

C...  load ammonium aitken aerosol mass into the CEND slot for
C...   ammonium accumulation aerosol if there is no slot for it in CEND

          IF ( NSGRAER( LNH4AKN ) .EQ. 0 ) THEN
            AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LNH4AKN )
            AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LNH4AKN )
            AEROSOL( LNH4AKN ) = 0.0
            AERWDEP( LNH4AKN ) = 0.0
          END IF

C...  load ammonium coarse aerosol mass into the CEND slot for
C...   ammonium accumulation aerosol if there is no slot for it in CEND 

          IF ( NSGRAER( LNH4COR ) .EQ. 0 ) THEN 
            AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LNH4COR )
            AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LNH4COR )
            AEROSOL( LNH4COR ) = 0.0
            AERWDEP( LNH4COR ) = 0.0
          END IF

        ELSE IF ( IAER .EQ. LNH4COR ) THEN 

C...  check for NH3, NH4ACC, and NH4AKN in cgrid...and as
C...    a last resort, put all of the ammonium mass into NH4COR

          IF ( ( NSGRGAS( LNH3   ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LNH4ACC ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LNH4AKN ) .EQ. 0 ) ) THEN
            AEROSOL( IAER ) = AEROSOL( IAER    )
     &                      + GAS    ( LNH3   )
     &                      + AEROSOL( LNH4AKN )
     &                      + AEROSOL( LNH4ACC )
            AERWDEP( IAER ) = AERWDEP( IAER    )
     &                      + GASWDEP( LNH3   )
     &                      + AERWDEP( LNH4AKN )
     &                      + AERWDEP( LNH4ACC )
            GAS    ( LNH3   ) = 0.0
            GASWDEP( LNH3   ) = 0.0
            AEROSOL( LNH4AKN ) = 0.0
            AERWDEP( LNH4AKN ) = 0.0
            AEROSOL( LNH4ACC ) = 0.0
            AERWDEP( LNH4ACC ) = 0.0
          END IF

        ELSE IF ( IAER .EQ. LNO3AKN ) THEN

C...  check for HNO3 in cgrid...and if there then add it mass to HNO3

          IF ( ( NSGRGAS( LHNO3   ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LNO3ACC ) .EQ. 0 ) ) THEN

            AEROSOL( IAER ) = AEROSOL( IAER    )
     &                      + GAS    ( LHNO3   )
     &                      + AEROSOL( LNO3ACC )
            AERWDEP( IAER ) = AERWDEP( IAER    )
     &                      + GASWDEP( LHNO3   )
     &                      + AERWDEP( LNO3ACC )
            GAS    ( LHNO3    ) = 0.0
            GASWDEP( LHNO3    ) = 0.0
            AEROSOL( LNO3ACC  ) = 0.0
            AERWDEP( LNO3ACC  ) = 0.0

            IF ( NSGRAER( LNO3COR ) .EQ. 0 ) THEN
              AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LNO3COR )
              AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LNO3COR )
              AEROSOL( LNO3COR ) = 0.0
              AERWDEP( LNO3COR ) = 0.0
            END IF

          END IF

        ELSE IF ( IAER .EQ. LNO3ACC ) THEN

C...  check for HNO3 in cgrid...and if there then add it mass to HNO3

          IF ( NSGRGAS( LHNO3 ) .EQ. 0 ) THEN
            AEROSOL( IAER  ) = AEROSOL( IAER ) + GAS    ( LHNO3 )
            AERWDEP( IAER  ) = AERWDEP( IAER ) + GASWDEP( LHNO3 )
            GAS    ( LHNO3 ) = 0.0
            GASWDEP( LHNO3 ) = 0.0
          END IF

C...  check for NO3AKN in cgrid...and if there then add it mass to HNO3

          IF ( NSGRAER( LNO3AKN ) .EQ. 0 ) THEN
            AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LNO3AKN )
            AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LNO3AKN )
            AEROSOL( LNO3AKN ) = 0.0
            AERWDEP( LNO3AKN ) = 0.0
          END IF

C...  check for NO3COR in cgrid...and if there then add it mass to HNO3

          IF ( NSGRAER( LNO3COR ) .EQ. 0 ) THEN
            AEROSOL( IAER    ) = AEROSOL( IAER ) + AEROSOL( LNO3COR )
            AERWDEP( IAER    ) = AERWDEP( IAER ) + AERWDEP( LNO3COR )
            AEROSOL( LNO3COR ) = 0.0
            AERWDEP( LNO3COR ) = 0.0
          END IF

        ELSE IF ( IAER .EQ. LNO3COR ) THEN

C...  check for HNO3, NO3ACC, and NO3AKN in cgrid...and as
C...    a last resort, put all of the nitrate mass into NO3COR

          IF ( ( NSGRGAS( LHNO3   ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LNO3ACC ) .EQ. 0 ) .AND.
     &         ( NSGRAER( LNO3AKN ) .EQ. 0 ) ) THEN
            AEROSOL( IAER ) = AEROSOL( IAER    )
     &                      + GAS    ( LHNO3   )
     &                      + AEROSOL( LNO3AKN )
     &                      + AEROSOL( LNO3ACC )
            AERWDEP( IAER ) = AERWDEP( IAER    )
     &                      + GASWDEP( LHNO3   )
     &                      + AERWDEP( LNO3AKN )
     &                      + AERWDEP( LNO3ACC )
            GAS    ( LHNO3   ) = 0.0
            GASWDEP( LHNO3   ) = 0.0
            AEROSOL( LNO3AKN ) = 0.0
            AERWDEP( LNO3AKN ) = 0.0
            AEROSOL( LNO3ACC ) = 0.0
            AERWDEP( LNO3ACC ) = 0.0
          END IF

        END IF

      END DO

C...Now, re-apportion mass back into cend/remov (cgrid-type) array

      DO IGAS = 1, NGAS
        DO ISGR = 1, NSGRGAS( IGAS )
          PNTR = CGRID2AQ_MAP( LSGRGAS( IGAS, ISGR ) )
          CEND ( PNTR ) = GAS    ( IGAS ) * WSGRGAS( IGAS, ISGR )
          REMOV( PNTR ) = GASWDEP( IGAS ) * WSGRGAS( IGAS, ISGR )
        END DO
      END DO

      DO IAER = 1, NAER
        DO ISGR = 1, NSGRAER( IAER )
          PNTR = CGRID2AQ_MAP( LSGRAER( IAER, ISGR ) )
          CEND ( PNTR ) = AEROSOL( IAER ) * WSGRAER( IAER, ISGR )
          REMOV( PNTR ) = AERWDEP( IAER ) * WSGRAER( IAER, ISGR )
        END DO
      END DO

      RETURN
      END
