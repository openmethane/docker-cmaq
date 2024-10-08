
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/PMEM_DEFN.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE PMEM_DEFN

C set up PM species mapping from emissions file(s) for use in GET_AERO_EMIS

C Revision History:
C        24 Aug 2007 J.Young: initial
C-----------------------------------------------------------------------

      IMPLICIT NONE

C Number of PM species in emission input file
      INTEGER, PARAMETER :: N_EMIS_PM = 6

C Aerosol (PM) species names in emission input files
      CHARACTER( 16 ), SAVE :: PM_EMIS( N_EMIS_PM )

C Descriptive variables from the emission input file
      CHARACTER( 16 ) :: UNITSCKPM            ! units of 1st PM species
      CHARACTER( 16 ) :: UNITSPM( N_EMIS_PM ) ! units of all PM species

C Indices in the PM_EM array
      INTEGER, SAVE :: VPOC  = 1   ! primary fine organic carbon
      INTEGER, SAVE :: VPSO4 = 2   ! primary fine sulfate
      INTEGER, SAVE :: VPNO3 = 3   ! primary fine nitrate
      INTEGER, SAVE :: VPEC  = 4   ! primary fine elemental carbon
      INTEGER, SAVE :: VPMF  = 5   ! unspeciated fine aerosol
      INTEGER, SAVE :: VPMC  = 6   ! primary coarse aerosol

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

C Arguments:
            INTEGER, INTENT( IN ) :: JDATE, JTIME

C External Functions:
            INTEGER, EXTERNAL :: INDEX1, TRIMLEN, SETUP_LOGDEV

C Local Variables:
            CHARACTER( 16 ), SAVE :: PNAME = 'GET_PMEM_MAP'
            CHARACTER( 16 ) :: VNAME
            CHARACTER( 96 ) :: XMSG = ' '

!           LOGICAL, SAVE :: FIRSTIME = .TRUE.
            INTEGER, SAVE :: LOGDEV

            INTEGER INDX, V, IOS

C-----------------------------------------------------------------------

            LOGDEV = SETUP_LOGDEV()

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
                  XMSG = 'PM Units not uniform on ' // EMIS_1
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            END DO

            RETURN

         END SUBROUTINE GET_PMEM_MAP

      END MODULE PMEM_DEFN

