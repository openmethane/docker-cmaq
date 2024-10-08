
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhg/PMEM_DEFN.F,v 1.1.1.1 2010/06/14 16:03:08 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE PMEM_DEFN

C set up PM species mapping from emissions file(s) for use in GET_AERO_EMIS

C Revision History:
C        24 Aug 2007 J.Young: initial
C        17 Sept 2009 B.Hutzell: enables different PM emissions
C                                needed based on the mechanism name
C-----------------------------------------------------------------------

      IMPLICIT NONE

C Number of PM species in emission input file
      INTEGER, SAVE      :: N_EMIS_PM     
      INTEGER, PARAMETER :: N_EMIS_PM_CAP = 6
      INTEGER, PARAMETER :: N_EMIS_PM_TX  = 20
      INTEGER, PARAMETER :: N_EMIS_PM_HG  = 1
      INTEGER, PARAMETER :: N_EMIS_PM_HAP = N_EMIS_PM_TX + N_EMIS_PM_HG


C Aerosol (PM) species names in emission input files
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: PM_EMIS( : )

C Descriptive variables from the emission input file
      CHARACTER( 16 ), SAVE               :: UNITSCKPM    ! units of 1st PM species
      CHARACTER( 16 ), ALLOCATABLE, SAVE  :: UNITSPM( : ) ! units of all PM species

C Indices in the PM_EM_CAP array
      INTEGER, SAVE :: VPOC  = 1   ! primary fine organic carbon
      INTEGER, SAVE :: VPSO4 = 2   ! primary fine sulfate
      INTEGER, SAVE :: VPNO3 = 3   ! primary fine nitrate
      INTEGER, SAVE :: VPEC  = 4   ! primary fine elemental carbon
      INTEGER, SAVE :: VPMF  = 5   ! unspeciated fine aerosol
      INTEGER, SAVE :: VPMC  = 6   ! primary coarse aerosol

C PM variables in emissions file
      CHARACTER( 16 ), SAVE :: EMIS_PM_HAP( N_EMIS_PM_HAP )   ! in emission input files
      INTEGER, SAVE         :: EMIS_MAP(    N_EMIS_PM_HAP )   ! indices emissions 

!BACKWARD COMPATIBLE with POA
C       DATA EMIS_PM_HAP( 1 ) / 'DIESEL_POC'    /   
! 
       DATA EMIS_PM_HAP( 2 ) / 'DIESEL_PSO4'   /
       DATA EMIS_PM_HAP( 3 ) / 'DIESEL_PNO3'   /
       DATA EMIS_PM_HAP( 4 ) / 'DIESEL_PEC'    /
       DATA EMIS_PM_HAP( 5 ) / 'DIESEL_PMFINE' /
       DATA EMIS_PM_HAP( 6 ) / 'DIESEL_PMC'    /
       DATA EMIS_PM_HAP( 7 ) / 'BERYLLIUM_F'   /
       DATA EMIS_PM_HAP( 8 ) / 'BERYLLIUM_C'   /
       DATA EMIS_PM_HAP( 9 ) / 'CHROMHEX_F'    /
       DATA EMIS_PM_HAP( 10) / 'CHROMHEX_C'    /
       DATA EMIS_PM_HAP( 11) / 'CHROMTRI_F'    /
       DATA EMIS_PM_HAP( 12) / 'CHROMTRI_C'    /
       DATA EMIS_PM_HAP( 13) / 'LEAD_F'        /
       DATA EMIS_PM_HAP( 14) / 'LEAD_C'        /
       DATA EMIS_PM_HAP( 15) / 'MANGANESE_F'   /
       DATA EMIS_PM_HAP( 16) / 'MANGANESE_C'   /
       DATA EMIS_PM_HAP( 17) / 'NICKEL_F'      /
       DATA EMIS_PM_HAP( 18) / 'NICKEL_C'      /
       DATA EMIS_PM_HAP( 19) / 'CADMIUM_F'     /
       DATA EMIS_PM_HAP( 20) / 'CADMIUM_C'     /
       DATA EMIS_PM_HAP( 21) / 'PHGI'          /

      INTEGER, SAVE :: INDX10 = 0  ! flag when PMC is calc'd from PM10

C-----------------------------------------------------------------------

      CONTAINS

         SUBROUTINE GET_PMEM_MAP ( JDATE, JTIME )

            USE HGRD_DEFN

            INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
!           INCLUDE SUBST_IODECL    ! I/O definitions and declarations
#include            SUBST_IODECL    # I/O definitions and declarations
            INCLUDE SUBST_IOFDESC   ! file header data structure
            INCLUDE SUBST_FILES_ID  ! file name parameters
            INCLUDE SUBST_RXCMMN    ! to get mech name

C Arguments:
            INTEGER, INTENT( IN ) :: JDATE, JTIME

C External Functions:
            INTEGER, EXTERNAL :: INDEX1, TRIMLEN, SETUP_LOGDEV

C Local Variables:
            CHARACTER( 16 ), SAVE :: PNAME = 'GET_PMEM_MAP'
            CHARACTER( 16 ) :: VNAME
            CHARACTER( 180 ) :: XMSG = ' '

!           LOGICAL, SAVE :: FIRSTIME = .TRUE.
            INTEGER, SAVE :: LOGDEV
            LOGICAL, SAVE :: AETOX_VRSN
            LOGICAL, SAVE :: AE_MERCURY

            INTEGER INDX, L, V, IOS
            INTEGER ALLOCSTAT
            INTEGER START, FINISH

C-----------------------------------------------------------------------

            LOGDEV = SETUP_LOGDEV()


            FLAG_HAPS: SELECT CASE( MECHNAME )
                CASE('CB05CLTX_AE4_AQ','CB05CLTX_AE5_AQ')
                     AETOX_VRSN = .TRUE.
                     AE_MERCURY = .FALSE.
                     N_EMIS_PM  = N_EMIS_PM_CAP
     &                          + N_EMIS_PM_TX
                CASE('SAPRC99TX3_AE4_AQ','SAPRC99TX3_AE5_AQ')
                     AETOX_VRSN = .TRUE.
                     AE_MERCURY = .FALSE.
                     N_EMIS_PM  = N_EMIS_PM_CAP
     &                          + N_EMIS_PM_TX
                CASE('CB05TXHG_AE4_AQ','CB05TXHG_AE5_AQ')
                     AETOX_VRSN = .TRUE.
                     AE_MERCURY = .TRUE.
                     N_EMIS_PM  = N_EMIS_PM_CAP
     &                          + N_EMIS_PM_HAP
                CASE('CB05HG_AE4_AQ','CB05HG_AE5_AQ')
                     AETOX_VRSN = .FALSE.
                     AE_MERCURY = .TRUE.
                     N_EMIS_PM  = N_EMIS_PM_CAP
     &                          + N_EMIS_PM_HG
                CASE DEFAULT
                     AETOX_VRSN = .FALSE.
                     AE_MERCURY = .FALSE.
                     N_EMIS_PM  = N_EMIS_PM_CAP
            END SELECT FLAG_HAPS

            ALLOCATE ( PM_EMIS( N_EMIS_PM ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** PM_EM memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( UNITSPM( N_EMIS_PM ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** UNITSPM memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

C Open the gridded emissions file, which contains gas, aerosol, and non-reactive
C species

            IF ( .NOT. OPEN3( EMIS_1, FSREAD3, PNAME ) ) THEN
               XMSG = 'Could not open '// EMIS_1 // ' file'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            IF ( .NOT. DESC3( EMIS_1 ) ) THEN
               XMSG = 'Could not get '// EMIS_1 // ' file description'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

C Search emissions file for certain species names. For each species found,
C   1. Assign a value to the index variables in the PM_EM array
C      (e.g., VPOC, VPSO4, etc.)
C   2. Populate the PM_EMIS array with the emitted species names
C      Note: species names are hardcoded in the present version
C   3. Populate the UNITSPM array based on the units in which the
C      emissions data are input (e.g., 'G/S', 'KG/H')
C If a species is not found, print error message and halt.

            VNAME ='POC'
            INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
            IF ( INDX .LE. 0 ) THEN
               VNAME ='POA'
               INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
               IF ( INDX .NE. 0 ) THEN
                  XMSG = 'Emissions file contains POA instead of POC'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               ELSE
                  XMSG = 'Could not find '
     &                 // VNAME( 1:TRIMLEN( VNAME ) )
     &                 // ' in emissions file'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
            END IF
            PM_EMIS( VPOC ) = VNAME
            UNITSCKPM = UNITS3D( INDX )
            UNITSPM( VPOC ) = UNITS3D( INDX )

            VNAME ='PSO4'
            INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
            IF ( INDX .NE. 0 ) THEN
               PM_EMIS( VPSO4 ) = VNAME
               UNITSPM( VPSO4 ) = UNITS3D( INDX )
            ELSE
               XMSG = 'Could not find '
     &              // VNAME( 1:TRIMLEN( VNAME ) )
     &              // ' in emissions file'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'PNO3'
            INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
            IF ( INDX .NE. 0 ) THEN
               PM_EMIS( VPNO3 ) = VNAME
               UNITSPM( VPNO3 ) = UNITS3D( INDX )
            ELSE
               XMSG = 'Could not find '
     &              // VNAME( 1:TRIMLEN( VNAME ) )
     &              // ' in emissions file'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'PEC'
            INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
            IF ( INDX .NE. 0 ) THEN
               PM_EMIS( VPEC ) = VNAME
               UNITSPM( VPEC ) = UNITS3D( INDX )
            ELSE
               XMSG = 'Could not find '
     &              // VNAME( 1:TRIMLEN( VNAME ) )
     &              // ' in emissions file'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'PMFINE'
            INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
            IF ( INDX .NE. 0 ) THEN
               PM_EMIS( VPMF ) = VNAME
               UNITSPM( VPMF ) = UNITS3D( INDX )
            ELSE
               XMSG = 'Could not find '
     &              // VNAME( 1:TRIMLEN( VNAME ) )
     &              // ' in emissions file'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

C Compatibility with emission files that do not list 'PMC' explicitly:
C   1. If 'PMC' variable exists, use those data for coarse PM emissions.
C   2. Otherwise, look for the 'PM10' variable.  If found, assign INDX10
C      a non-zero value to flag the fact that PM10 was read instead of PMC.
C   3. If 'PM10' is not found either, print error message and halt.

            VNAME = 'PMC'
            INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
            IF ( INDX .EQ. 0 ) THEN   ! try another
               XMSG = 'Could not find PMC in emissions file.'
     &             // '  Using PM10 - sum(PM2.5) instead.'
               CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               VNAME = 'PM10'
               INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
               INDX10 = INDX
            END IF

            IF ( INDX .NE. 0 ) THEN
               PM_EMIS( VPMC ) = VNAME
               UNITSPM( VPMC ) = UNITS3D( INDX )
            ELSE
               XMSG = 'Could not find PMC nor PM10 in emissions file'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF
  
            START  =  N_EMIS_PM_CAP + 1

            IF( AETOX_VRSN )THEN

              FINISH = START + N_EMIS_PM_TX - 1

              V = START
              VNAME = 'DIESEL_POC'
              INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
              IF ( INDX .EQ. 0 ) THEN
                  XMSG = 'FAILED to find ' // VNAME //
     &                   ' in emissions file. Attempting ' //
     &                   ' to find DIESEL_POA '
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
                  VNAME = 'DIESEL_POA'
                  INDX  = INDEX1( VNAME, NVARS3D, VNAME3D )
              ENDIF
              IF ( INDX .NE. 0 ) THEN
                   UNITSPM( V ) = UNITS3D( INDX )
                   PM_EMIS( V ) = VNAME
                   EMIS_PM_HAP(V - N_EMIS_PM_CAP) = VNAME
              ELSE
                   XMSG = 'Could not find '
     &                  // VNAME( 1:LEN_TRIM( VNAME ) )
     &                 // ' in aerosol table'
                   CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
              ENDIF
                           
              DO V = START+1, FINISH
                VNAME = EMIS_PM_HAP( V - N_EMIS_PM_CAP )
                INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
                IF ( INDX .NE. 0 ) THEN
                   PM_EMIS( V ) = VNAME
                   UNITSPM( V ) = UNITS3D( INDX )
                ELSE
                   XMSG = 'Could not find '
     &                 // VNAME( 1:TRIMLEN( VNAME ) )
     &                 // ' in emissions file'
                   CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
                ENDIF
              ENDDO

              START = FINISH + 1

            ENDIF


            IF( AE_MERCURY )THEN
 
              FINISH = START + N_EMIS_PM_HG - 1
              L  = 0

              DO V = START, FINISH
                L = L + 1
                VNAME = EMIS_PM_HAP( L + N_EMIS_PM_TX )
                INDX = INDEX1( VNAME, NVARS3D, VNAME3D )
                IF ( INDX .NE. 0 ) THEN
                   PM_EMIS( V ) = VNAME
                   UNITSPM( V ) = UNITS3D( INDX )
                ELSE
                   XMSG = 'Could not find '
     &                 // VNAME( 1:TRIMLEN( VNAME ) )
     &                 // ' in emissions file'
                   CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
                ENDIF
              ENDDO

            ENDIF

C Save the number of aerosol species read from the emissions file and write their
C names to the log file

            WRITE( LOGDEV, '( /5X, I3, 2X, A )' ) N_EMIS_PM,
     &                         'PM species in emissions file:'
            DO V = 1, N_EMIS_PM
               WRITE( LOGDEV, '( 15X, A )' ) PM_EMIS( V )
            END DO

            WRITE( LOGDEV,* ) ' '

C Confirm that all aerosol species in the emission input file have the same units.
C  If not, print error message and halt.

            DO V = 1, N_EMIS_PM
               IF ( UNITSPM( V ) .NE. UNITSCKPM ) THEN
                     XMSG = 'PM Units not uniform in PM_EMIS species ' //
     &                      ' Required Units are '    // UNITSCKPM //
     &                      ' Emissions '   // PM_EMIS( V ) // 
     &                      ' has units '   // UNITSPM( V )
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            END DO

            RETURN

         END SUBROUTINE GET_PMEM_MAP

      END MODULE PMEM_DEFN

