
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/init/yamo_ddm3d/initscen.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)initscen.F        1.4 /project/mod3/CMAQ/src/init/init/SCCS/s.initscen.F 25 Jul 1997 14:32:52

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef sens
      SUBROUTINE INITSCEN ( CGRID, SENGRID, STDATE, STTIME, TSTEP,
     &                      NSTEPS, NPMAX, IPT, IPARM, ISPLIT,
     &                      IAMOUNT, ILAYER, NDAYS,
     &                      IDATE, ITIME, IREGION, IRXN, SEN_PAR)
#else
      SUBROUTINE INITSCEN ( CGRID, STDATE, STTIME, TSTEP, NSTEPS )
#endif

C-----------------------------------------------------------------------
C Function:
C   Initialize simulation time period and time stepping constants for
C   core model driver
C   Environment variable can reference a previous CONC file to use as
C   initial data.
C   Write initial conc data as step "0" on output conc file
 
C Preconditions:
 
C Subroutines and functions called:
C   INIT3, M3EXIT, OPEN3, CLOSE3, DESC3, ENVINT, TIME2SEC, HHMMSS
 
C Revision history:
C   prototype 6/92 by CJC for "concept model"
C   revised   2/93 by CJC for initial M3 Prototype Alpha
C   revised   6/14/94 by Dongming Hwang for configuration management prototype
C   revised   8/30/94 by Dongming Hwang: If CTM_CONC_1 file doesn't exist,
C             use header info from CHEMINIT3D file.
 
C Interim Version ... (Jeff)

C   revised  6/12/97 by Jerry Gipson: Get ICs by species name, by surrogate
C            name, or zero
C   Jeff - Aug 97 - fixed problems, cleaned up
C   Jeff - Dec 97 - add CMIN
C   Jeff - Dec 97 - put in aerosol sulfate inititalization
C   Jeff - Feb 98 - close init cond files after reading
C   2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   revised 10/7/99 by Shawn Roselle: added surface area species to
C            aerosol species types
C   Jeff - Dec 00 - check if append, split out opconc and load_cgrid
C                 - move CGRID_MAP into f90 module
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; assumed shape arrays
C   17 Mar 03 D.Wong: move barrier to avoid race conditions
C   28 Aug 03 J.Young: following Zion Wang at CERT, only pe 0 closes
C   30 May 05 J.Young: add call to load RHOJ into CGRID
C   15 Nov 08 S.Napelenok DDM-3D implementation for CMAQv4.7
C-----------------------------------------------------------------------

      USE HGRD_DEFN             ! horizontal domain specifications
!     USE CGRID_SPCS            ! CGRID species number and offsets

      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex
#ifdef sens
      USE VGRD_DEFN             ! vertical layer specifications              
#endif

      IMPLICIT NONE

C Includes:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
!     INCLUDE SUBST_GC_ICBC     ! gas chem ic/bc surrogate names and map table
!     INCLUDE SUBST_AE_ICBC     ! aerosol ic/bc surrogate names and map table
!     INCLUDE SUBST_NR_ICBC     ! non-react ic/bc surrogate names and map table
!     INCLUDE SUBST_TR_ICBC     ! tracer ic/bc surrogate names and map table
!     INCLUDE SUBST_GC_CONC     ! gas chem conc file species and map table
!     INCLUDE SUBST_AE_CONC     ! aerosol conc file species and map table
!     INCLUDE SUBST_GC_G2AE     ! gas chem aerosol species and map table
!     INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
!     INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID    ! coord. and domain definitions (req IOPARMS)
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

!     REAL      :: CGRID( :,:,:,: )  ! for initial CONC
      REAL, POINTER :: CGRID( :,:,:,: )  ! for initial CONC
      INTEGER      STDATE       ! starting date,    format YYYYDDD
      INTEGER      STTIME       ! starting time,    format HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
      INTEGER      NSTEPS       ! run duration, as number of output time steps

C Parameters:

      REAL, PARAMETER :: CMIN = 1.0E-30

C External Functions (not already declared by IODECL3.EXT):

      INTEGER, EXTERNAL :: ENVINT       !  get environment variable as integer
      INTEGER      STATUS               !  ENVINT status
      INTEGER, EXTERNAL :: TIME2SEC     !  converts HHMMSS to raw seconds
      INTEGER, EXTERNAL :: TRIMLEN      !  string length, excl. trailing blanks

      CHARACTER( 10 ), EXTERNAL :: HHMMSS  !  converts to string "HH:MM:SS"

C Local Variables

C  environment variable start date
      CHARACTER( 16 ) :: CTM_STDATE = 'CTM_STDATE'
C  environment variable start time
      CHARACTER( 16 ) :: CTM_STTIME = 'CTM_STTIME'
C  environment variable run duration
      CHARACTER( 16 ) :: CTM_RUNLEN = 'CTM_RUNLEN'
C  environment variable output time step
      CHARACTER( 16 ) :: CTM_TSTEP = 'CTM_TSTEP'
C driver program name
      CHARACTER( 16 ) :: CTM_PROGNAME = 'CTM_PROGNAME'
C  environment variable description
      CHARACTER( 80 ) :: VARDESC

      CHARACTER( 16 ) :: PNAME = 'INITSCEN'
      CHARACTER( 16 ) :: FNAME
      CHARACTER( 16 ) :: PROGNAME
      CHARACTER( 96 ) :: MSG = ' '

      INTEGER      RUNLEN            ! run duration, HHMMSS
      INTEGER      STEPSECS          ! seconds per time step
      INTEGER      TOTSECS           ! run duration seconds
      INTEGER      LOGDEV            ! FORTRAN unit number for log file
!     INTEGER      C, R, L, SPC      ! loop counters

#ifdef sens
c  VARIABLES SPECIFIC TO SENSITIVITY
!     REAL SENGRID( NCOLS,NROWS,NLAYS,NPMAX,* ) ! sensitivities
      REAL, POINTER :: SENGRID( :,:,:,:,: )  ! for initial sens
      INTEGER  NPMAX                    ! Number of sens parameters
      INTEGER  P                        ! Loop counter over NPMAX
      INTEGER, POINTER :: IPT( : )
!     INTEGER IPT(NPMAX)                ! Type of sensitivity
                                        ! 1 = INIT; 2 = BOUN; 3 = EMIS
                                        ! 4 = HIGH; 5 = RATE
      INTEGER, POINTER :: IPARM( :,: )
!     INTEGER IPARM( NPMAX, MXSPCS )    ! set to 1 if sens to that species
                                        ! set to 0 otherwise
      INTEGER, POINTER :: ISPLIT( :,: )
!     INTEGER ISPLIT( NPMAX, 6 )        ! 1 = total; 2 = area; 3 = bio
                                        ! 4 = mobile; 5 = point; 6 = non-rd
      REAL, POINTER :: IAMOUNT( :,:,: )
!      REAL IAMOUNT ( NPMAX, MXSPCS, 25 )  ! emission rate during the hour
      INTEGER, POINTER :: ILAYER( :,: )
!      INTEGER ILAYER ( NPMAX, NLAYS )   ! 1 if that layer's emis included, 
                                        ! 0 else
                                        ! only used if emis amt specified,
                                        ! otherwise set to 1 for all layers

      ! IRXN( NPMAX, NRXNS ) = 1 if sens to that reaction rate;
      !                      = 0 otherwise
      INTEGER, POINTER :: IRXN( :,: )

      INTEGER  NDAYS
      INTEGER, POINTER :: IDATE( :,: )
!      INTEGER IDATE ( NPMAX, NDAYS )   ! date(s) on which emission or boundary
                                        ! condition should occur
      INTEGER, POINTER :: ITIME( :,: )
!      INTEGER ITIME ( NPMAX, 2 )        ! begin time & end time for each sens
                                        ! default is 0, 240000
                                        ! Note: to span midnight, begin time
                                        ! may be larger than endtime
!      REAL IREGION( SENCOLS, SENROWS, NLAYS, NPMAX )
      REAL, POINTER ::  IREGION( :,:,:,: )
      CHARACTER( 8 ), POINTER :: SEN_PAR ( : )
!      CHARACTER( 8 ) :: SEN_PAR(NPMAX)  ! names for sens parameter
      INTEGER      C, R, L, SPC      ! loop counters

      INTEGER       MXSPCS              ! Number of species in CGRID
      PARAMETER   ( MXSPCS = N_GC_SPCD
     &                     + N_AE_SPC
     &                     + N_NR_SPC
     &                     + N_TR_SPC )

      LOGICAL, EXTERNAL      :: ENVYN                    ! get environment variable as boolean
      CHARACTER( 16 ), SAVE  :: DDM3D_RST = 'DDM3D_RST'  ! switch for sens restart file
      LOGICAL, SAVE          :: RST

#endif

      INTERFACE
         SUBROUTINE LOAD_CGRID ( FNAME, STDATE, STTIME, SPC_CAT, CMIN, CGRID )
            IMPLICIT NONE
            CHARACTER( 16 ), INTENT( IN ) :: FNAME
            INTEGER, INTENT( IN )         :: STDATE, STTIME
            CHARACTER(  2 ), INTENT( IN ) :: SPC_CAT
            REAL,    INTENT( IN )         :: CMIN
            REAL, POINTER                 :: CGRID( :,:,:,: )
         END SUBROUTINE LOAD_CGRID
         SUBROUTINE OPCONC ( CGRID, STDATE, STTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER                 :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )         :: STDATE, STTIME, TSTEP
         END SUBROUTINE OPCONC
#ifdef sens
         SUBROUTINE LOAD_SENGRID ( FNAME, STDATE, STTIME, SPC_CAT,
     &                             SENGRID, NPMAX, SEN_PAR )
            IMPLICIT NONE
            CHARACTER( 16 ), INTENT( IN ) :: FNAME
            INTEGER, INTENT( IN )         :: STDATE, STTIME
            CHARACTER(  2 ), INTENT( IN ) :: SPC_CAT
            REAL, POINTER                 :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )         :: NPMAX
            CHARACTER( 8 ), POINTER  :: SEN_PAR( : )
         END SUBROUTINE LOAD_SENGRID
         SUBROUTINE SINPUT ( STDATE, STTIME, NPMAX, IPT, IPARM,
     &                ISPLIT, IAMOUNT, ILAYER, NDAYS, IDATE,
     &                ITIME, IREGION, IRXN, SEN_PAR )
            IMPLICIT NONE
            INTEGER, INTENT( IN )         :: STDATE, STTIME
            INTEGER, INTENT( IN )         :: NPMAX
            INTEGER, POINTER              :: IPT( : )
            INTEGER, POINTER        :: IPARM( :, : )
            INTEGER, POINTER        :: ISPLIT( :, : )
            REAL, POINTER           :: IAMOUNT( :, :, : )
            INTEGER, POINTER        :: ILAYER( :, : )
            INTEGER, INTENT( IN )   :: NDAYS
            INTEGER, POINTER        :: IDATE( :, : )
            INTEGER, POINTER        :: ITIME( :, : )
            REAL, POINTER           :: IREGION( :, :, :, : )
            INTEGER, POINTER        :: IRXN( :,: )
            CHARACTER( 8 ), POINTER :: SEN_PAR( : )
         END SUBROUTINE SINPUT
#endif
      END INTERFACE

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

      VARDESC = 'Main Program Name'
      CALL ENVSTR( CTM_PROGNAME, VARDESC, 'DRIVER', PROGNAME, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 'Environment variable set, but empty...Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            END IF

!     WRITE( LOGDEV, '( 5X, A)') 'Program ' // PROGNAME( 1:TRIMLEN( PROGNAME ))

!     WRITE( LOGDEV, '( 5X, A)') 'INPUT VALUES:'

      STDATE = 1995192        ! default
      VARDESC = 'Scenario Starting Date (YYYYDDD)'
      STDATE = ENVINT( CTM_STDATE, VARDESC, STDATE, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 'Environment variable set, but empty...Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            END IF

      STTIME = 000000         ! default
      VARDESC = 'Scenario Starting Time (HHMMSS)'
      STTIME = ENVINT( CTM_STTIME, VARDESC, STTIME, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 'Environment variable set, but empty...Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            END IF

      RUNLEN = 480000         ! default
      VARDESC = 'Scenario Run Duration (HHMMSS)'
      RUNLEN = ENVINT( CTM_RUNLEN, VARDESC, RUNLEN, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 'Environment variable set, but empty...Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            END IF

      TSTEP( 1 )  = 010000         ! default
      VARDESC = 'Scenario Output Time Step (HHMMSS)'
      TSTEP( 1 ) = ENVINT( CTM_TSTEP, VARDESC, TSTEP( 1 ), STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 'Environment variable set, but empty...Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
            END IF

#ifdef sens
         RST = .TRUE.         ! default
         VARDESC = 'Flag for using sensitivity restart files'
         RST = ENVYN( DDM3D_RST, VARDESC, RST, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 
     &          'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, STTIME
         END IF
#endif

C Calculate number of output time steps for this model run:

      IF ( RUNLEN .LT. 1000000 ) THEN
         TOTSECS  = TIME2SEC( RUNLEN )
         ELSE
         RUNLEN = RUNLEN - 1000000
         TOTSECS  = TIME2SEC( RUNLEN )
         TOTSECS  = TOTSECS + 360000
         END IF
      STEPSECS = TIME2SEC( TSTEP( 1 )  )

      IF ( MOD( TOTSECS, STEPSECS ) .EQ. 0 ) THEN
         NSTEPS = TOTSECS / STEPSECS
         ELSE
         MSG = 'Output time step ' // HHMMSS( TSTEP( 1 ) ) //
     &         ' does not divide duration ' // HHMMSS( RUNLEN )
         CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT1 )
         END IF

#ifdef sens
C This subroutine reads in the names, types and number of sens parameters
         CALL SINPUT( STDATE, STTIME, NPMAX, IPT, IPARM, ISPLIT,
     &                IAMOUNT, ILAYER, NDAYS, IDATE, ITIME, IREGION,
     &                IRXN, SEN_PAR )
#endif

C Initialize the CGRID array

!     DO SPC = 1, NSPCSD   ! load CGRID with this many variables
!        DO L = 1, NLAYS
!           DO R = 1, MY_NROWS
!              DO C = 1, MY_NCOLS
!                 CGRID( C,R,L,SPC ) = CMIN
!                 END DO
!              END DO
!           END DO
!        END DO

      CGRID = CMIN

C Load CGRID

      IF ( N_GC_SPC .GT. 0 ) THEN

         FNAME = INIT_GASC_1

         CALL LOAD_CGRID ( FNAME, STDATE, STTIME, 'GC', CMIN, CGRID )

         END IF  ! N_GC_SPC > 0

C load rhoj

      FNAME = MET_CRO_3D

      CALL LOAD_CGRID ( FNAME, STDATE, STTIME, 'RJ', CMIN, CGRID )

C Get aerosols IC's. Note: INIT_AERO_1 may be setenv'ed to INIT_GASC_1
 
      IF ( N_AE_SPC .GT. 0 ) THEN

         FNAME = INIT_AERO_1
 
         CALL LOAD_CGRID ( FNAME, STDATE, STTIME, 'AE', CMIN, CGRID )

         END IF  ! N_AE_SPC > 0

C Get non-reactives IC's. Note: INIT_NONR_1 may be setenv'ed to INIT_GASC_1

      IF ( N_NR_SPC .GT. 0 ) THEN

         FNAME = INIT_NONR_1

         CALL LOAD_CGRID ( FNAME, STDATE, STTIME, 'NR', CMIN, CGRID )

         END IF  ! N_NR_SPC > 0

C Get tracer IC's. Note: INIT_TRAC_1 may be setenv'ed to INIT_GASC_1

      IF ( N_TR_SPC .GT. 0 ) THEN

         FNAME = INIT_TRAC_1

         CALL LOAD_CGRID ( FNAME, STDATE, STTIME, 'TR', CMIN, CGRID )

         END IF  ! N_TR_SPC > 0

#ifdef sens
C Begin Sensitivity Calls
C SENGRID is initially set to CGRID if we're doing sens to initial conditions
C Otherwise, set to 0
      IF ( RST ) THEN
C Load SENGRID

         IF ( N_GC_SPC .GT. 0 ) THEN
 
            FNAME = INIT_GASC_S

            CALL LOAD_SENGRID ( FNAME, STDATE, STTIME, 'GC', SENGRID,
     &                          NPMAX, SEN_PAR )

            END IF  ! N_GC_SPC > 0

C DSC: Is it necessary to load RJ into SENGRID???

            FNAME = MET_CRO_3D

            CALL LOAD_SENGRID ( FNAME, STDATE, STTIME, 'RJ', SENGRID,
     &                          NPMAX, SEN_PAR )

C Get aerosols IC's. Note: INIT_AERO_S may be setenv'ed to INIT_GASC_S

         IF ( N_AE_SPC .GT. 0 ) THEN
  
            FNAME = INIT_AERO_S

            CALL LOAD_SENGRID ( FNAME, STDATE, STTIME, 'AE', SENGRID,
     &                          NPMAX, SEN_PAR )

            END IF  ! N_AE_SPC > 0

C Get non-reactives IC's. Note: INIT_NONR_S may be setenv'ed to INIT_GASC_S

         IF ( N_NR_SPC .GT. 0 ) THEN

            FNAME = INIT_NONR_S

            CALL LOAD_SENGRID ( FNAME, STDATE, STTIME, 'NR', SENGRID,
     &                          NPMAX, SEN_PAR )

            END IF  ! N_NR_SPC > 0

C Get tracer IC's. Note: INIT_TRAC_1 may be setenv'ed to INIT_GASC_1
!
!         IF ( N_TR_SPC .GT. 0 ) THEN
! 
!            FNAME = INIT_TRAC_S
!
!            CALL LOAD_SENGRID ( FNAME, STDATE, STTIME, 'TR', SENGRID,
!     &                          NPMAX, SEN_PAR )
!
!            END IF  ! N_TR_SPC > 0

      ELSE

       DO P = 1, NPMAX
         ! Sens to Initial Conditions
         IF ( IPT(P) .EQ. 1 ) THEN
           DO SPC = 1, MXSPCS
             DO L = 1, NLAYS
               DO R = 1, NROWS
                 DO C = 1, NCOLS
                   SENGRID( C, R, L, P, SPC) =
     &                   CGRID( C, R, L, SPC ) *
     &                   REAL( IPARM( P, SPC ) ) *
     &                   IREGION( C, R, L, P )
                 END DO
               END DO
             END DO
           END DO
         ELSE
           DO SPC = 1, MXSPCS
             DO L = 1, NLAYS
               DO R = 1, NROWS
                 DO C = 1, NCOLS
                   SENGRID(C,R,L,P,SPC) = 0.0
                 END DO
               END DO
             END DO
           END DO
         END IF
       END DO

      END IF
#endif

C Test opening existing conc file for update

      IF ( .NOT. OPEN3( CTM_CONC_1, FSRDWR3, PNAME ) ) THEN

         MSG = 'Could not open ' // CTM_CONC_1 // ' file for update - '
     &        // 'try to open new'
         CALL M3MESG( MSG )

C Open the file new

         CALL SUBST_BARRIER

         CALL OPCONC ( CGRID, STDATE, STTIME, TSTEP( 1 ) )

         END IF

!     IF ( MYPE .EQ. 0 ) THEN

         IF ( .NOT. CLOSE3( CTM_CONC_1 ) ) THEN
            MSG = 'Could not close ' // CTM_CONC_1
!           CALL M3EXIT( PNAME, STDATE, STTIME, MSG, XSTAT1 )
            CALL M3WARN( PNAME, STDATE, STTIME, MSG )
            END IF

!        END IF

      CALL SUBST_BARRIER

      RETURN

      END
