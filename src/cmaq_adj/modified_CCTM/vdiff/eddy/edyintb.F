 
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/eddy/edyintb.F,v 1.1.1.1 2010/06/14 16:03:10 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)edyintb.F	1.2 /project/mod3/CMAQ/src/vdiff/eddy/SCCS/s.edyintb.F 16 Jun 1997 15:47:48

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE EDYINTB ( EDDYV, DT, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Calculates vertical eddy diffusivity using integration method
C   Estimates max. admissible diffusion time step for
C   Crank-Nicholson method (with some safety factor)
C   Integrated Version of Heinz Formula - by Daewon Byun
C     ZL, ZU; low and upper limit of integration
C     ZSL   ; surface layer thickness
C     PBL   ; boundary layer thickness
C     Summary of Cases
C     Case 1 : 0 <= ZL < ZU  < ZSL < PBL (A,B,C)
C     Case 2 : 0 <= ZL < ZSL < ZU  < PBL (A,B,C) + (D,E)
C     Case 3 : 0 <= ZL < ZSL < PBL < ZU  (A,B,C) + (D,E) + (F)
C     Case 4 : 0 < ZSL < ZL  < ZU  < PBL (D,E)
C     Case 5 : 0 < ZSL < ZL  < PBL < ZU  (D,E)   + (F)
C     Case 6 : 0 < ZSL < PBL < ZL  < ZU  (F)
C     Where (A): stable surface layer formula
C           (B): neutral surface layer formula
C           (C): unstable surface layer formula
C           (D): stable pbl formula
C           (E): unstable mixed layer formula
C           (F): formula for free atmosphere
 
C Preconditions:

C Subroutines and Functions Called:
C   INIT3, M3EXIT, SEC2TIME, TIME2SEC, INDEX1

C Revision history:
C   NO.   DATE     WHO    WHAT
C          Sep 97  Jeff  generalize, get vert. data from COORD.EXT
C       12 Mar 96  Jeff  for free trop. modifications
C       24 Feb 96  Jeff  for upgraded Models-3 I/O API
C   2   8/16/94  Dongming Hwang  CM template
C   1   7/1/93     CJC   Now uses INTERP3()
C   2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   3 Mar 1999 by Dave Wong at LM: extend diagnostic output to parallel
C                                  implementation
C     Aug 1999 by Daewon Byun: found that Kz in the PBL min is not set
C                          Reset min(Kz) to be KZ0UT
C   21 Nov 00 J.Young: use Dave Wong's f90 stenex GLOBAL_MIN_DATA
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX
C   30 Aug 02 J.Young: fix ZINT - use ZH
C   27 Jan 03 J.Young: restore ZINT, fix UWIND, VWIND dimensions
C   15 Apr 04 J.Pleim: use urban fraction to set min Kz
C   28 Jan 05 J.Young: optional kz0uT to Min Kz
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   16 Feb 06 J.Young: cap kzmin at 500 m
C   07 Mar 06 J.Young: trap *all* eddyv's at kzmin; improve logging max
C                      eddyv's
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex

      IMPLICIT NONE
 
C Includes:
 
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
!     INCLUDE SUBST_COORD_ID    ! coord. and domain definitions (req IOPARMS)
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Arguments:

!     REAL         EDDYV ( NCOLS,NROWS,NLAYS ) ! eddy diffusivity (m**2/s)
!     REAL         DT    ( NCOLS,NROWS )       ! computed diffusion time step
      REAL      :: EDDYV ( :,:,: )     ! eddy diffusivity (m**2/s)
      REAL      :: DT    ( :,: )       ! computed diffusion time step
      INTEGER      JDATE        ! current model date , coded YYYYDDD
      INTEGER      JTIME        ! current model time , coded HHMMSS
      INTEGER      TSTEP        ! sciproc sync. step (chem)

C Parameters:

C von Karman constant
      REAL, PARAMETER :: KARM = 0.4
      REAL, PARAMETER :: KARMC = KARM / 0.74
C asymptotic mixing length
      REAL, PARAMETER :: XLAMB = 1600.0
C default eddy diffusivity (m**2/sec) KZ0 for Upper Troposphere
      REAL, PARAMETER :: KZ0UT = 1.0
!     REAL, PARAMETER :: KZL = 0.1                 ! lowest KZ
      REAL, PARAMETER :: KZL = 0.5                 ! lowest KZ
      REAL, PARAMETER :: KZU = 2.0                 ! highest KZ
      REAL, PARAMETER :: ETA = 0.5

C lower limit for absolute inverse Monin-Obukhov Length
      REAL, PARAMETER :: MINMOLI = 1.0E-04 ! 1/10,000

!     REAL, PARAMETER :: MXFTEDDY = 10.0 ! m**2/sec

C External Functions not previously declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, INDEX1
      LOGICAL, EXTERNAL :: ENVYN

!                                123456789012345678901234567890
      CHARACTER( 30 ) :: MSG1 = ' Error interpolating variable '

C File Variables:

      REAL          PBL  ( NCOLS,NROWS )           ! pbl height (m)
      REAL          USTAR( NCOLS,NROWS )           ! friction velocity
      REAL          MOLI ( NCOLS,NROWS )           ! inverse Monin-Obukhov Len
      REAL          WSTAR( NCOLS,NROWS )           ! convective velocity scale
      REAL          ZH   ( NCOLS,NROWS,NLAYS )     ! mid-layer elevation
      REAL          ZF   ( NCOLS,NROWS,0:NLAYS )   ! full layer elevation
      REAL          TA   ( NCOLS,NROWS,NLAYS )     ! temperature (K)
      REAL          QV   ( NCOLS,NROWS,NLAYS )     ! water vapor mixing ratio
      REAL          PRES ( NCOLS,NROWS,NLAYS )     ! pressure

      REAL, ALLOCATABLE, SAVE :: MSFX2( :,: )      ! Squared map scale factors

      LOGICAL, SAVE :: MINKZ
!     REAL          PURB( NCOLS,NROWS )            ! percent urban
      REAL, ALLOCATABLE, SAVE ::PURB( :,: )        ! percent urban
      REAL          UFRAC                          ! urban fraction
      REAL, ALLOCATABLE, SAVE :: KZMIN( :,:,: )    ! minimum Kz (m**2/s)
      REAL          KZM                            ! local KZMIN
      LOGICAL, ALLOCATABLE, SAVE :: KZLAY( :,:,: ) ! minimum Kz applied
      REAL, PARAMETER :: KZMAXL = 500.0            ! upper limit for min Kz (m)

      REAL          UWIND( NCOLS+1,NROWS+1,NLAYS ) ! x-direction winds
      REAL          VWIND( NCOLS+1,NROWS+1,NLAYS ) ! y-direction winds
!     REAL, ALLOCATABLE, SAVE :: UWIND( :,:,: )
!     REAL, ALLOCATABLE, SAVE :: VWIND( :,:,: )
      REAL, ALLOCATABLE, SAVE :: UVBUF( :,:,: )    ! U, V read buffer
      INTEGER, SAVE :: MCOLS, MROWS                ! for allocating

C Local variables:

      LOGICAL,SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 16 ) :: PNAME = 'EDYINTB'
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 16 ) :: UNITSCK
      CHARACTER( 120 ) :: XMSG = ' '

      REAL, SAVE :: P0          ! 1000 mb reference pressure
      REAL, SAVE :: CONVPA      ! Pressure conversion factor file units to Pa

      INTEGER      ASTAT
      INTEGER      MDATE, MTIME, STEP
      INTEGER      C, R, L, V
      INTEGER      ICASE( NCOLS,NROWS,NLAYS-1 ), ICS
      INTEGER      KFTB ( NCOLS,NROWS )    ! free troposphere height
      LOGICAL      FREETROP
      LOGICAL      NEVER

      REAL         DTSEC
      REAL         TV                      ! virtual temperature (K)
      REAL         DZL                     ! Z(L+1)-Z(L)
      REAL         WW2                     ! (wind speed)**2
      REAL         WS2                     ! (wind shear)**2
      REAL         RIB                     ! Bulk Richardson Number
      REAL         RIC                     ! critical Richardson Number
      REAL         RL, RU, ZL, ZU
      REAL         HEAD, ARG1, BETA
      REAL         THETAV( NCOLS,NROWS,NLAYS )    ! potential temp
      REAL         ZINT  ( NCOLS,NROWS,NLAYS )
      REAL         DZINT ( NCOLS,NROWS )
      REAL         EDV                     ! local eddyv

      INTEGER      GXOFF, GYOFF            ! global origin offset from file
C for INTERPX
      INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2
      INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER, SAVE :: STRTCOLMD3, ENDCOLMD3, STRTROWMD3, ENDROWMD3

      INTEGER      MCOL                   ! these don't need to be initialized
      INTEGER      MROW
      INTEGER      MLVL
      INTEGER      MTROP                  ! free troposphere flag
      REAL         MTH1                   ! pot. temp. in layer L
      REAL         MTH2                   ! pot. temp. in layer L+1
      REAL         MRIB                   ! bulk Richardson Number
      REAL         MWS                    ! wind shear (/sec)
      REAL         MEDDYV                 ! max eddy diffusivity (m**2/sec)
      REAL         DT0

      INTEGER, SAVE :: LOGDEV

      LOGICAL STABLE

C Statement Functions:

      REAL        FKA, FKC, FKD, FKE, FREE

C Statement Function Parameters and Arguments:

      REAL, PARAMETER :: THIRD = 1.0 / 3.0
      REAL, PARAMETER :: FIFTH = 1.0 / 5.0

      REAL       ALPHA, AS, AS2, AS3, X1, X2, BETA1
      REAL       Z1, Z2, ZK, ZPBL, ZSL, ZWST

C Function Definitions:

      FKA( Z1, Z2, BETA1 ) =
     &                       ( Z2 - Z1 
     &                     - LOG( ( BETA1 * Z2 + 1.0 )
     &                          / ( BETA1 * Z1 + 1.0 ) )
     &                     / BETA1 ) / BETA1

      FKC( Z1, Z2, ALPHA ) =
     &     ( ( 3.0 * ALPHA * Z2 - 2.0 ) * SQRT( 1.0 + ALPHA * Z2 )**3
     &     - ( 3.0 * ALPHA * Z1 - 2.0 ) * SQRT( 1.0 + ALPHA * Z1 )**3
     &     ) / ( 7.5 * ALPHA * ALPHA )

      FKD( X1, X2, AS2, AS3 ) = AS3 * (
     &                  AS2 * FIFTH * ( X1 ** 5 - X2 ** 5 )
     &           + ( AS2 - 1.0 ) * (
     &                        THIRD * ( X1 ** 3 - X2 ** 3 )
     &           +                      X1      - X2
     &               - 0.5 * ( LOG( (1.0 + X1)
     &                            / (1.0 - X1) )
     &                       - LOG( (1.0 + X2)
     &                            / (1.0 - X2) )
     &                       )
     &                             )
     &                                )

      FKE( Z1, Z2, ZK, ZWST, ZPBL ) =
     &   ZK * ZWST * ( Z2 * Z2 * ( 0.5 - THIRD * Z2 / ZPBL )
     &               - Z1 * Z1 * ( 0.5 - THIRD * Z1 / ZPBL ) )

      FREE( Z1, Z2, ZK ) = ZK * (Z2 - Z1)

C-----------------------------------------------------------------------

      IF ( FIRSTIME )  THEN
         FIRSTIME  =  .FALSE.
         LOGDEV = INIT3 ()

C Open the met files

         IF ( .NOT. OPEN3( MET_CRO_3D, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open '// MET_CRO_3D // ' file'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
                          
         IF ( .NOT. DESC3( MET_CRO_3D ) ) THEN
            XMSG = 'Could not get ' // MET_CRO_3D // ' file description'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF         !  error abort if if desc3() failed

         V = INDEX1( 'PRES', NVARS3D, VNAME3D )
         IF ( V .NE. 0 ) THEN
            UNITSCK = UNITS3D( V )
            ELSE
            XMSG = 'Could not get variable PRES from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

         IF ( UNITSCK .EQ. 'PASCAL' .OR. UNITSCK .EQ. 'pascal' .OR.
     &        UNITSCK .EQ. 'Pascal' .OR. UNITSCK .EQ. 'PA'     .OR.
     &        UNITSCK .EQ. 'pa'     .OR. UNITSCK .EQ. 'Pa' ) THEN
            CONVPA = 1.0
            P0 = 100000.0
            ELSE IF ( UNITSCK .EQ. 'MILLIBAR' .OR. UNITSCK .EQ. 'millibar' .OR.
     &                UNITSCK .EQ. 'Millibar' .OR. UNITSCK .EQ. 'MB'       .OR.
     &                UNITSCK .EQ. 'mb'       .OR. UNITSCK .EQ. 'Mb' ) THEN
            CONVPA = 1.0E-02
            P0 = 1000.0
            ELSE IF ( UNITSCK .EQ. 'CENTIBAR' .OR. UNITSCK .EQ. 'centibar' .OR.
     &                UNITSCK .EQ. 'Centibar' .OR. UNITSCK .EQ. 'CB'       .OR.
     &                UNITSCK .EQ. 'cb'       .OR. UNITSCK .EQ. 'Cb' ) THEN
            CONVPA = 1.0E-03
            P0 = 100.0
            ELSE
            XMSG = 'Units incorrect on ' // MET_CRO_3D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

!        IF ( .NOT. OPEN3( MET_DOT_3D, FSREAD3, PNAME ) ) THEN
!           XMSG = 'Could not open '// MET_DOT_3D // ' file'
!           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!           END IF
                          
!        IF ( .NOT. OPEN3( GRID_CRO_2D, FSREAD3, PNAME ) ) THEN
!           XMSG = 'Could not open '// GRID_CRO_2D // ' file'
!           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!           END IF

         ALLOCATE ( MSFX2( NCOLS,NROWS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating MSFX2'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         CALL SUBHFILE ( GRID_CRO_2D, GXOFF, GYOFF,
     &                   STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2 )

         VNAME = 'MSFX2'
         IF ( .NOT. INTERPX( GRID_CRO_2D, VNAME, PNAME,
!    &                       1,NCOLS, 1,NROWS, 1,1,
     &                       STRTCOLGC2,ENDCOLGC2, STRTROWGC2,ENDROWGC2, 1,1,
     &                       JDATE, JTIME, MSFX2 ) ) THEN
            XMSG = MSG1 // VNAME // ' from ' // GRID_CRO_2D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
 
         MINKZ = .TRUE.   ! default
         MINKZ = ENVYN( 'KZMIN', 'Kz min on flag', MINKZ, ASTAT )
         IF ( ASTAT .NE. 0 ) WRITE( LOGDEV,'(5X, A)' ) 'Kz min on flag'
         IF ( ASTAT .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            ELSE IF ( ASTAT .EQ. -1 ) THEN
            XMSG = 'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV,'(5X, A)' ) XMSG
            ELSE IF ( ASTAT .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV,'(5X, A)' ) XMSG
            END IF

         IF ( MINKZ ) THEN

            ALLOCATE ( PURB( NCOLS,NROWS ), STAT = ASTAT )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating PURB'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            VNAME = 'PURB'
            IF ( .NOT. INTERPX( GRID_CRO_2D, VNAME, PNAME,
     &                          STRTCOLGC2,ENDCOLGC2, STRTROWGC2,ENDROWGC2, 1,1,
     &                          JDATE, JTIME, PURB ) ) THEN
               XMSG = 'Either make the data available from MCIP'
               WRITE( LOGDEV,'(/5X, A)' ) XMSG
               XMSG = 'or set the env var KZMIN to F or N,'
               WRITE( LOGDEV,'( 5X, A)' ) XMSG
               XMSG = 'in which case you will revert back to the'
               WRITE( LOGDEV,'( 5X, A)' ) XMSG
               XMSG = 'previous version of subroutine edyintb using Kz0UT'
               WRITE( LOGDEV,'( 5X, A, /)' ) XMSG
               XMSG = ' '
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            ALLOCATE ( KZLAY( NCOLS,NROWS,NLAYS ), STAT = ASTAT )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating KZLAY'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            ELSE IF ( .NOT. MINKZ ) THEN
            XMSG = 'This run uses Kz0UT, *NOT* KZMIN in subroutine edyintb.'
            WRITE( LOGDEV,'(/5X, A, /)' ) XMSG

            END IF   ! MINKZ

         ALLOCATE ( KZMIN( NCOLS,NROWS,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating KZMIN'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         CALL SUBHFILE ( MET_CRO_2D, GXOFF, GYOFF,
     &                   STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2 )
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )
         CALL SUBHFILE ( MET_DOT_3D, GXOFF, GYOFF,
     &                   STRTCOLMD3, ENDCOLMD3, STRTROWMD3, ENDROWMD3 )

         MCOLS = ENDCOLMD3 - STRTCOLMD3 + 1
         MROWS = ENDROWMD3 - STRTROWMD3 + 1

         ALLOCATE ( UVBUF( MCOLS,MROWS,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating UVBUF'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         END IF          !  if firstime

C Interpolate time dependent layered input variables

      MDATE  = JDATE
      MTIME  = JTIME
      STEP   = TIME2SEC( TSTEP )
      CALL NEXTIME( MDATE, MTIME, SEC2TIME( STEP / 2 ) )

      DTSEC = FLOAT( STEP )

      VNAME = 'UWIND'
      IF ( .NOT. INTERPX( MET_DOT_3D, VNAME, PNAME,
     &                    STRTCOLMD3,ENDCOLMD3, STRTROWMD3,ENDROWMD3, 1,NLAYS,
     &                    MDATE, MTIME, UVBUF ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_DOT_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      UWIND = 0.0
      DO L = 1, NLAYS
         DO R = 1, MROWS
            DO C = 1, MCOLS
               UWIND( C,R,L ) = UVBUF( C,R,L )
               END DO
            END DO
         END DO

      VNAME = 'VWIND'
      IF ( .NOT. INTERPX( MET_DOT_3D, VNAME, PNAME,
     &                    STRTCOLMD3,ENDCOLMD3, STRTROWMD3,ENDROWMD3, 1,NLAYS,
     &                    MDATE, MTIME, UVBUF ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_DOT_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VWIND = 0.0
      DO L = 1, NLAYS
         DO R = 1, MROWS
            DO C = 1, MCOLS
               VWIND( C,R,L ) = UVBUF( C,R,L )
               END DO
            END DO
         END DO

      VNAME = 'PBL'
      IF ( .NOT. INTERPX( MET_CRO_2D, VNAME, PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    MDATE, MTIME, PBL ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_2D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'USTAR'
      IF ( .NOT. INTERPX( MET_CRO_2D, VNAME, PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    MDATE, MTIME, USTAR ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_2D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'MOLI'
      IF ( .NOT. INTERPX( MET_CRO_2D, VNAME, PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    MDATE, MTIME, MOLI ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_2D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'WSTAR'
      IF ( .NOT. INTERPX( MET_CRO_2D, VNAME, PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    MDATE, MTIME, WSTAR ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_2D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'TA'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, TA ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'QV'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, QV ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'PRES'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, PRES ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'ZF'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, ZF ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Move 3rd dimension slabbed data from INTERP3 into proper order
C ( Using ZF both as a read buffer and an argument variable.)
 
      IF ( MINKZ ) THEN
         KZLAY = .FALSE.
         DO L = NLAYS, 1, -1
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  ZF( C,R,L ) = ZF( C,R,L-1 )
                  IF ( ZF( C,R,L ) .LE. KZMAXL ) KZLAY( C,R,L ) = .TRUE.
                  END DO
               END DO
            END DO
         ELSE
         DO L = NLAYS, 1, -1
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  ZF( C,R,L ) = ZF( C,R,L-1 )
                  END DO
               END DO
            END DO
         END IF
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            ZF( C,R,0 ) = 0.0
            END DO
         END DO

      IF ( MINKZ ) THEN
         KZMIN = KZL
         DO L = 1, NLAYS
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  IF ( KZLAY( C,R,L ) ) THEN
                     UFRAC = 0.01 * PURB( C,R )
                     KZMIN( C,R,L ) = KZL + ( KZU - KZL ) * UFRAC
                     END IF
                  END DO
               END DO
            END DO
         ELSE
         KZMIN = KZ0UT
         END IF

      VNAME = 'ZH'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, ZH ) ) THEN
         XMSG = MSG1 // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
C bound MOLI away from (-e-04, e-04)
            MOLI( C,R ) = SIGN ( 
     &                         MAX (
     &                             ABS ( MOLI( C,R ) ),
     &                             MINMOLI
     &                             ),
     &                         MOLI( C,R )
     &                         )
            TV              = TA( C,R,1 ) * ( 1.0 + 0.6 * QV( C,R,1 ) )
            THETAV( C,R,1 ) = TV * ( P0 / PRES( C,R,1 ) ) ** 0.286
            DZINT( C,R )  = ETA * ZF( C,R,1 )
            ZINT( C,R,1 ) = DZINT( C,R )
            KFTB( C,R ) = NLAYS - 1
            END DO
         END DO

      DO L = 2, NLAYS
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               TV              = TA( C,R,L ) * ( 1.0 + 0.6 * QV( C,R,L ) )
               THETAV( C,R,L ) = TV * ( P0 / PRES( C,R,L ) ) ** 0.286
               ZINT ( C,R,L ) =  ZF( C,R,L-1 ) + DZINT( C,R )
               DZINT( C,R )   =  ZF( C,R,L ) - ZINT( C,R,L )
               END DO
            END DO
         END DO

C Determine layer indices, RADM layer cases, KFTB

      DO 188 L = 1, NLAYS-1
      DO 177 R = 1, MY_NROWS
      DO 166 C = 1, MY_NCOLS

C Surface layer height is 10% of pbl, but at most 50 meters
         ZSL  = MIN ( 50.0, 0.1 * PBL( C,R ) )

         IF ( ZINT( C,R,L ) .LT. ZSL ) THEN
            IF ( ZINT( C,R,L+1 ) .LT. ZSL ) THEN
               ICASE( C,R,L ) = 1
               ELSE IF ( ZINT( C,R,L+1 ) .LT. PBL( C,R ) ) THEN
               ICASE( C,R,L ) = 2
               ELSE
               ICASE( C,R,L ) = 3
               END IF
            ELSE IF ( ZINT( C,R,L ) .LT. PBL( C,R ) ) THEN
               IF ( ZINT( C,R,L+1 ) .LT. PBL( C,R ) ) THEN
               ICASE( C,R,L ) = 4
               ELSE
               ICASE( C,R,L ) = 5
               END IF
            ELSE
            ICASE( C,R,L ) = 6
            KFTB( C,R ) = MIN ( KFTB( C,R ), L )
            END IF

166      CONTINUE
177      CONTINUE
188      CONTINUE

      NEVER = .TRUE.
      DT0 = DTSEC

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DT( C,R ) = DTSEC
            END DO
         END DO

C get ghost values for wind fields in case of free trop.

      CALL SUBST_COMM ( UWIND, DSPL_N0_E1_S0_W0, DRCN_E )
      CALL SUBST_COMM ( VWIND, DSPL_N1_E0_S0_W0, DRCN_N )

      MEDDYV = 0.0
      DO 233 L = 1, NLAYS - 1
      DO 222 R = 1, MY_NROWS
      DO 211 C = 1, MY_NCOLS
 
         ARG1 = KARMC * USTAR( C,R )
         ZSL  = MIN ( 50.0, 0.1 * PBL( C,R ) )
         ICS  = ICASE( C,R,L )
         KZM  = KZMIN( C,R,L )

C Depending on the case, use appropriate functions

         IF ( L .LT. KFTB( C,R ) ) THEN
            FREETROP = .FALSE.

            ZL = ZINT( C,R,L )
            ZU = ZINT( C,R,L+1 )

            IF ( MOLI( C,R ) .GT. 0.0 ) THEN      !  stable case:

               STABLE = .TRUE.

               BETA1 = 6.345 * MOLI( C,R )
               BETA  = PBL( C,R ) * BETA1   ! 4.7*pbl/(0.74*l)
               HEAD  = 2.0 * ARG1 * PBL( C,R ) * PBL( C,R ) / BETA
               AS2   = ( 1.0 + BETA ) / BETA
               AS    = SQRT( AS2 )
               AS3   = AS * AS2

               IF ( ICS .EQ. 1 ) THEN

                  EDV = ARG1 * FKA( ZL, ZU, BETA1 ) / ( ZU - ZL )
                  EDV = MAX ( KZM, EDV )

                  ELSE IF ( ICS .EQ. 2 ) THEN

                  RL = SQRT( MAX ( 1.0E-5, 1.0 - ZSL / PBL( C,R ) ) ) / AS
                  RU = SQRT( 1.0 - ZU / PBL( C,R ) ) / AS
                  EDV = ( ARG1 * FKA( ZL, ZSL, BETA1 )
     &                +   HEAD * FKD( RL, RU, AS2, AS3 ) ) / ( ZU - ZL )
                  EDV = MAX ( KZM, EDV )

                  ELSE IF ( ICS .EQ. 3 ) THEN

                  RL = SQRT( MAX ( 1.0E-5, 1.0 - ZSL / PBL( C,R ) ) ) / AS
                  RU = 0.0
                  EDV = ( ARG1 * FKA( ZL, ZSL, BETA1 )
     &                +   HEAD * FKD( RL, RU, AS2, AS3 )
     &                +   FREE( PBL( C,R ), ZU, KZM ) ) / ( ZU - ZL )
                  EDV = MAX ( KZM, EDV )

                  ELSE IF ( ICS .EQ. 4 ) THEN

                  RL = SQRT( 1.0 - ZL / PBL( C,R ) ) / AS
                  RU = SQRT( 1.0 - ZU / PBL( C,R ) ) / AS
                  EDV = HEAD * FKD( RL, RU, AS2, AS3 ) / ( ZU - ZL )
                  EDV = MAX ( KZM, EDV )

                  ELSE IF ( ICS .EQ. 5 ) THEN

                  RL = SQRT( 1.0 - ZL / PBL( C,R ) ) / AS
                  RU = 0.0
                  EDV = ( HEAD * FKD( RL, RU, AS2, AS3 )
     &                +   FREE( PBL( C,R ), ZU, KZM ) ) / ( ZU - ZL )
                  EDV = MAX ( KZM, EDV )

                  ELSE		!  case 6

                  EDV = KZM

                  END IF

               ELSE                      ! unstable case:  MOL < 0

               STABLE = .FALSE.

               ALPHA = -9.0 * MOLI( C,R )

               IF ( ICS .EQ. 1 ) THEN

                  EDV = ARG1 * FKC( ZL, ZU, ALPHA ) / ( ZU - ZL )

                  ELSE IF ( ICS .EQ. 2 ) THEN

                  EDV = ( ARG1 * FKC( ZL, ZSL, ALPHA )
     &                +   FKE( ZSL, ZU, KARM, WSTAR( C,R ), PBL( C,R ) ) )
     &                / ( ZU - ZL )

                  ELSE IF ( ICS .EQ. 3 ) THEN

                  EDV = ( ARG1 * FKC( ZL, ZSL, ALPHA )
     &                +   FKE( ZSL, PBL( C,R ), KARM, WSTAR( C,R ), PBL( C,R ) )
     &                +   FREE( PBL( C,R ), ZU, KZM ) ) / ( ZU - ZL )

                  ELSE IF ( ICS .EQ. 4 ) THEN

                  EDV = FKE( ZL, ZU, KARM, WSTAR( C,R ), PBL( C,R ) ) 
     &                / ( ZU - ZL )

                  ELSE IF ( ICS .EQ. 5 ) THEN

                  EDV = ( FKE( ZL, PBL( C,R ), KARM, WSTAR( C,R ), PBL( C,R ) )
     &                +   FREE( PBL( C,R ), ZU, KZM ) ) / ( ZU - ZL )

                  ELSE	!  case 6

                  EDV = KZM

                  END IF

               END IF    !  if stable case, or else unstable case

cccccc error section begin ccccccccccccccccccccccccccccccccccccccccccccc
               IF ( EDV .LE. 0.0 ) THEN
                  WRITE( LOGDEV,* ) ' '
                  WRITE( LOGDEV,* ) '   *** Eddy Diffusivity .le. 0 ***'
                  WRITE( LOGDEV,* ) ' Stable: ', STABLE
                  WRITE( LOGDEV,* ) ' C,R,L, ICASE, EDDYV = ',
     &                          C,R,L, ICS, EDV
                  WRITE( LOGDEV,* ) ' ZL, ZU, PBL, MOL, AS, HEAD = ',
     &                          ZL, ZU, PBL( C,R ), 1.0/MOLI( C,R ), AS, HEAD
                  IF ( STABLE ) THEN
                     IF ( ICS .GE. 2 .AND.
     &                    ICS .LT. 6 ) THEN
                        WRITE( LOGDEV,* ) ' RL, RU = ', RL, RU
                        END IF
                     IF ( ICS .LE. 3 ) THEN
                        WRITE( LOGDEV,* ) ' FKA = ', FKA( ZL, ZU, BETA1 )
                        END IF
                     IF ( ICS .GE. 3 .AND.
     &                    ICS .LE. 5 ) THEN
                        WRITE( LOGDEV,* ) ' FKD = ', FKD( RL, RU, AS2, AS3 )
                        WRITE( LOGDEV,* ) ' FKD term1 = ',
     &                                      AS2 * FIFTH * (RL**5 - RU**5)
                        WRITE( LOGDEV,* ) ' FKD subterm2 = ',
     &                                      THIRD * (RL**3-RU**3)
                        WRITE( LOGDEV,* ) ' FKD subterm3 = ', (RL-RU)
                        WRITE( LOGDEV,* ) ' FKD ln arg''s = ',
     &                           (1.0+RL)/(1.0-RL), '    ', (1.0+RU)/(1.0-RU)
                        WRITE( LOGDEV,* ) ' sum FKD ln terms = ',
     &                       LOG( (1.0+RL)/(1.0-RL) ) - LOG( (1.0+RU)/(1.0-RU) )
                        END IF
                     ELSE
                     IF ( ICS .GE. 2 .AND.
     &                    ICS .LE. 5 ) THEN
                        WRITE( LOGDEV,* ) ' Wstar = ', WSTAR( C,R )
                        END IF
                     IF ( ICS .EQ. 1 ) THEN
                        WRITE( LOGDEV,* ) ' FKC = ', FKC( ZL, ZU, ALPHA )
                        END IF
                     IF ( ICS .GT. 1 .AND.
     &                    ICS .LE. 3 ) THEN
                        WRITE( LOGDEV,* ) ' FKC = ', FKC( ZL, ZSL, ALPHA )
                        END IF
                     IF ( ICS .EQ. 2 ) THEN
                        WRITE( LOGDEV,* ) ' FKE = ',
     &                  FKE( ZSL, ZU, KARM, WSTAR( C,R ), PBL( C,R ) )
                        END IF
                     IF ( ICS .EQ. 3 ) THEN
                        WRITE( LOGDEV,* ) ' FKE = ',
     &                  FKE( ZSL, PBL( C,R ), KARM, WSTAR( C,R ), PBL( C,R ) )
                        END IF
                     IF ( ICS .EQ. 4 ) THEN
                        WRITE( LOGDEV,* ) ' FKE = ',
     &                  FKE( ZL, ZU, KARM, WSTAR( C,R ), PBL( C,R ) )
                        END IF
                     IF ( ICS .EQ. 5 ) THEN
                        WRITE( LOGDEV,* ) ' FKE = ',
     &                  FKE( ZL, PBL( C,R ), KARM, WSTAR( C,R ),  PBL( C,R ) )
                        END IF
                     END IF
                  IF ( ICS .EQ. 3 .OR.
     &                 ICS .EQ. 5 ) THEN
                     WRITE( LOGDEV,* ) ' FREE = ', 
     &                                   FREE( PBL( C,R ), ZU, KZM )
                     END IF
                  XMSG = ' '
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
                  END IF
cccccc error section end ccccccccccccccccccccccccccccccccccccccccccccccc
 
C For upper layers:

            ELSE	!  L >= KFTB
            FREETROP = .TRUE.

            DZL = ZH( C,R,L+1 ) - ZH( C,R,L )
            RIC = 0.257 * DZL ** 0.175
            WW2 = 0.25 * MSFX2( C,R )      ! component-wise wind shear
     &          * ( ( UWIND( C+1,R,  L+1 ) - UWIND( C+1,R  ,L  )
     &              + UWIND( C,  R,  L+1 ) - UWIND( C,  R  ,L  ) ) ** 2
     &          +   ( VWIND( C,  R+1,L+1 ) - VWIND( C,  R+1,L )
     &              + VWIND( C,  R,  L+1 ) - VWIND( C,  R,  L  ) ) ** 2 )
            WS2 = WW2 / ( DZL * DZL ) + 1.0E-9
 
            RIB =   2.0 * GRAV * ( THETAV( C,R,L+1 ) - THETAV( C,R,L ) )
     &          / ( DZL * WS2  * ( THETAV( C,R,L+1 ) + THETAV( C,R,L ) ) )
            RIB = MAX ( RIB, 0.0 )

            IF ( ( RIB - RIC ) .GE. 0.0 ) THEN
               EDV = KZM
               ELSE
               EDV = KZM + XLAMB * SQRT( WS2 ) * ( 1.0 - RIB / RIC )
               END IF

cccccc error section begin ccccccccccccccccccccccccccccccccccccccccccccc
            IF ( EDV .LE. 0.0 ) THEN
               WRITE( LOGDEV,* ) '   *** Eddy Diffusivity .le. 0 ***'
               WRITE( LOGDEV,* ) ' C,R,L, EDDYV = ', C,R,L, EDV
               WRITE( LOGDEV,* ) ' DZL, RIC, WW2, WS2, RIB, HEAD = ',
     &                             DZL, RIC, WW2, WS2, RIB, HEAD
               XMSG = ' '
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
cccccc error section end ccccccccccccccccccccccccccccccccccccccccccccccc
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!           EDV = MAX ( EDV, MXFTEDDY )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            END IF	!  if L < KFTB, or not

         DT( C,R ) = MIN ( DT( C,R ),
     &                   0.75 * ( ZF( C,R,L )   - ZF( C,R,L-1 ) )
     &                        * ( ZH( C,R,L+1 ) - ZH( C,R,L ) )
     &                        / EDV )
 
!        IF ( DT( C,R ) .LT. DT0 ) THEN
         IF ( EDV .GT. MEDDYV ) THEN
C Capture the col, row, lvl, and EDDYV for the global min DT
!           DT0 = DT( C,R )
            MCOL = C
            MROW = R
            MLVL = L
            MEDDYV = EDV
            NEVER = .FALSE.
            MTROP = 0
            IF ( FREETROP ) THEN
               MTROP = 1
               MTH1 = THETAV( C,R,L )
               MTH2 = THETAV( C,R,L+1 )
               MRIB = RIB
               MWS  = SQRT ( WS2 )
               END IF
            END IF

         EDDYV( C,R,L ) = MAX ( KZM, EDV )

211      CONTINUE       !  end loop on columns
222      CONTINUE       !  end loop on rows
233      CONTINUE       !  end loop on levels

!     CALL SUBST_GLOBAL_MIN_DATA ( DT0, MCOL, MROW, MLVL, MEDDYV, NEVER,
!    &                             MTROP, FREETROP, MTH1, MTH2, MRIB, MWS )

      IF ( NEVER ) THEN
!        DT0 = DT( 1,1 )
         MCOL = 1
         MROW = 1
         MLVL = 1
         MEDDYV = EDDYV( 1,1,1 )
         MTROP = 0
         END IF
      WRITE( LOGDEV,* ) '    '
      WRITE( LOGDEV,1001 ) MEDDYV, MCOL, MROW, MLVL
1001  FORMAT(/ 5X, 'Maximum eddy diffusivity of:', 1PG13.5,
     &         1X, '(m**2/sec)'
     &       / 5X, 'at col, row, layer:', I4, 2(', ', I3) )
      IF ( MTROP .EQ. 1 ) THEN
         WRITE( LOGDEV,1003 ) MWS, MRIB, MTH1, MTH2, ZF( MCOL,MROW,MLVL )
1003     FORMAT(  5X, 'corresponding to a free tropospheric wind shear of:',
     &            1PG11.3,  1X, '(/sec),'
     &          /28X, 'a bulk Richardson Number of:', 1PG13.5, ','
     &          / 5X, 'and pot. temps. in layer and layer+1:', 2( 1PG13.5 ) 
     &          / 5X, 'Layer height above ground: ', 1PG13.5, ' meters' )
         END IF
      WRITE( LOGDEV,* ) ' '

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            EDDYV( C,R,NLAYS ) = 0.0
            END DO
         END DO

      RETURN
      END
