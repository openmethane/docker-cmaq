
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/init/init/flcheck.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE FLCHECK ( STDATE, STTIME, TSTEP )

C----------------------------------------------------------------------
C Function: To check file header data for consistency with COORD.EXT
C            and runlength; only gridded and boundary files are checked.
 
C Preconditions: None
 
C Key Subroutines/Functions Called:
 
C Revision History: Prototype created by Jerry Gipson, September, 1997
C                   Modified Jun, 1998 by Jerry Gipson to remove PING
C                   and MEPSE files from being checked
C                   Jeff - Nov 2000 - ck only input files for start date/time
C   30 Mar 01 J.Young: elim reference to grid
C----------------------------------------------------------------------

      IMPLICIT NONE 

C..INCLUDES:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
      
C..ARGUMENTS:
      INTEGER STDATE           ! Current date (YYYYDDD)
      INTEGER STTIME           ! Current time (HHMMSS)
      INTEGER TSTEP            ! Current time step (HHMMSS)

C..PARAMETERS:

C..EXTERNAL FUNCTIONS:
      INTEGER ENVINT    ! Gets integer value of environment variable 
      LOGICAL ENVYN     ! Gets logical value of environment variable

C..SAVED LOCAL VARIABLES: None
      
C..SCRATCH LOCAL VARIABLES:
      CHARACTER*80 MSG    ! Message
      CHARACTER*16 PNAME  ! Program name
      
      INTEGER JDATE     ! Time step date (YYYYDDD)
      INTEGER JTIME     ! Time step time (HHMMSS)
      INTEGER N         ! Loop index
      INTEGER RUNLEN    ! Length of run
      INTEGER STATUS    ! Status of returned value

      LOGICAL LERROFF   ! Flag to stop run if errors found
      LOGICAL LERROR    ! Flag to indicate error for 1 or more files
      LOGICAL LOK       ! Flag to indicate error for 1 file

C----------------------------------------------------------------------

      DATA PNAME/'FLCHECK'/

      LERROFF = .TRUE.
      LERROFF = ENVYN( 'FL_ERR_STOP', 'Flag for stopping on errors',
     &           LERROFF, STATUS )  

      LERROR = .FALSE.

      RUNLEN = 480000         ! default
      RUNLEN = ENVINT( 'CTM_RUNLEN', 'Scenario Run Duration (HHMMSS)',
     &                 RUNLEN, STATUS )


      DO N = 1, N_FILES      
        
C skip Ping and Mepse files
         IF ( INDEX( LGCL_NAME( N ), 'PING'  ) .EQ. 0 .AND. 
     &        INDEX( LGCL_NAME( N ), 'MEPSE' ) .EQ. 0 ) THEN

            LOK = OPEN3( LGCL_NAME( N ), FSREAD3, PNAME )

            IF ( LOK ) THEN

               IF ( DESC3( LGCL_NAME( N ) ) ) THEN

C only check gridded or boundary files
                  IF ( FTYPE3D .EQ. GRDDED3 .OR. FTYPE3D .EQ. BNDARY3 ) THEN

                     MSG = 'Checking header data for file: ' // LGCL_NAME( N )
                     CALL M3MESG( MSG )

                     CALL GRDCHECK ( LGCL_NAME( N ), LOK )

c check time steps for time dependent files
                     IF ( TSTEP3D .GT. 0 ) THEN
                        JDATE = STDATE
                        JTIME = STTIME

                        IF ( .NOT. CHECK3( LGCL_NAME( N ), ALLVAR3,
     &                                     JDATE, JTIME ) ) THEN
                           LOK = .FALSE.
                           MSG = 'Starting time not on file'
                           CALL M3MESG( MSG )
                        END IF
 
c skip checking ending time for output files
                        IF ( IN_MODE( N ) ) THEN
                           CALL NEXTIME ( JDATE, JTIME, RUNLEN )
                           IF ( .NOT. CHECK3( LGCL_NAME( N ), ALLVAR3,
     &                                        JDATE, JTIME ) ) THEN
                              LOK = .FALSE.
                              MSG = 'Ending time not on file'
                              CALL M3MESG( MSG )
                           END IF
                        END IF

                     END IF   ! if TSTEP3D .GT. 0

                     IF ( .NOT. LOK ) LERROR = .TRUE.

                  END IF   ! only check gridded or boundary files
               ELSE
               LOK = .FALSE.
               MSG = 'DESC3 failed for ' // LGCL_NAME( N )
               CALL M3MESG( MSG )
               END IF   ! DESC3 successful
            END IF   ! successful open

!           IF ( LOK ) THEN
               IF( .NOT. CLOSE3( LGCL_NAME( N ) ) ) THEN
                  MSG = '  >>> COULD NOT CLOSE ' // LGCL_NAME( N )
                  CALL M3WARN ( PNAME, STDATE, STTIME, MSG )
                  END IF
!              END IF

         END IF  ! skip Ping and Mepse files

      END DO

      IF ( LERROR ) THEN
         MSG = 'Inconsistent header data on input files'
         IF (LERROFF ) THEN       
            CALL M3EXIT ( PNAME, STDATE, STTIME, MSG, XSTAT2 )
            ELSE
            CALL M3WARN ( PNAME, STDATE, STTIME, MSG )
            RETURN
            END IF
         ELSE
         MSG = '   No inconsistent header data found on input files'
         CALL M3MESG( MSG )
         END IF

      END
