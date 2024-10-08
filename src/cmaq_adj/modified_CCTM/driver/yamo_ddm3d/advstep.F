
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_ddm3d/advstep.F,v 1.1.1.1 2010/06/14 16:03:04 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ADVSTEP( JDATE, JTIME, TSTEP, ASTEP, NREPS )

C-----------------------------------------------------------------------
C Function:
C    Use Courant-number conditions to determine the (time-splitting) model
C    base synchronization time step TSTEP(2) and number NREPS of sync time
C    steps per output time step TSTEP(1). The base synchronization time step
C    is considered to be the Courant condition safe advection time step.
C    Minimum returned sync time step = 900 sec
C    Determines a Courant condition safe advection step in TSTEP(3) with
C    respect to a user-defined minimum sync time step [TSTEP(2)].

C    June 2003 J.Young
C    Restructure to produce layer dependent advection step
C    change arguments:
C       TSTEP( 1 ) = output timestep (same) - input
C       TSTEP( 2 ) = synchronization timestep (same) - output
C       eliminate TSTEP( 3 ). Add new argument ASTEP( NLAYS ) - output

C Preconditions:

C Subroutines and functions called:
C    SEC2TIME, TIME2SEC, DESC3, M3EXIT, HVELOC, ENVREAL

C Revision history:
C    prototype 2/93 by CJC
C    configuration management prototype 6/94 by Dongming Hwang

C    Apr 15, 1996 by M. Talat Odman at NCSC:  the number of model time
C    steps is made even to assure symmetry of transport operators within
C    each output time step.

C    Jan 20, 1996 by M. Talat Odman at NCSC:  now applies the CFL condition
C    independently in all three contravariant directions and
C    finds minimum time step.

C    22 Apr 97 Jeff: only for horizontal advection
C                    keep Talat's even step mod for symmetric sciproc
C     2 Sep 97 Jeff: calc. separate sync and adv time steps, TSTEP(3)
C    12 Feb 98 Jeff: make sync step = chem step
C    27 Jun 98 Jeff: allow for constant vel field (time-independent file)

C    2 October, 1998 by Al Bourgeois at LM: parallel implementation
C
C    8 Jan 99 by David Wong at LM: GLOBAL_MAX function call is added
C
C    1/19/99 David Wong
C      -- add a loop_index call
C      -- change loop index ending point to avoid accessing invalid region.
C         (reason to do this is to prevent using boundary data from PINTERP,
C          which sets pseudo-boundary data to 0)
C
C    9 Jun 99 by David Wong at LM: change GLOBAL_MAX to GLOBAL_RMAX
C    21 Nov 00 J.Young: GLOBAL_RMAX -> Dave Wong's f90 stenex GLOBAL_MAX
C    Jeff - Feb 01 - module LOCALFILE to set CCTM IOAPI file defn's
C                  - assumed shape arrays
C   30 Mar 01 J.Young: Use HGRD_DEFN; replace INTERP3 with INTERPX
C    9 Sep 03 J.Young: return layer dependent advection tstep, ASTEP
C   27 Aug 03 J.Young: symmetric processing is now deprecated
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   26 Jul 06 J.Young: add WRF sigma VGTYP
C   27 Apr 07 J.Young: Using Peter Percell's (UH) two timestep max wind
C   20 Jan 09 J.Young: In addition to the CFL criterion, use horizontal divergence
C                      to establish safe advection time step
C   28 Oct 09 J.Young: Revert back from Percell's algorithm to Talat's fix for
C                      multiple adv steps per sync step
C   16 Mar 10 J.Young: Fix hdiv bug for MPI by communicating ghost cell data using
C                      new hveloc to add hdiv and replace hcontvel
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE SUBST_MODULES         ! stenex
!     USE SUBST_GLOBAL_MAX_MODULE    ! stenex
!     USE SUBST_UTIL_MODULE          ! stenex

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:
      INTEGER, INTENT( IN ) :: JDATE ! current model simulation date (YYYYDDD)
      INTEGER, INTENT( IN ) :: JTIME ! current model simulation time (HHMMSS)
      INTEGER, INTENT( INOUT ) :: TSTEP( 2 ) ! time step vector (HHMMSS)
                                ! TSTEP( 1 ) = local output step - INTENT( IN )
                                ! TSTEP( 2 ) = sciproc sync. step (chem) - INTENT( OUT )
      INTEGER, INTENT( OUT ) :: ASTEP( : )   ! layer advection sub-steps per sync step
                                             ! (HHMMSS)
      INTEGER, INTENT( OUT ) :: NREPS        ! no. of sync steps per output step

C Parameters:
      REAL, PARAMETER :: CC = 0.75 ! maximum Courant number allowed
!     REAL, PARAMETER :: CC = 0.99 ! maximum Courant number allowed

C External Functions (not already declared by IODECL3.EXT):
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, SECSDIFF, ENVINT
      LOGICAL, EXTERNAL :: ENVYN, CURRSTEP
      REAL,    EXTERNAL :: ENVREAL

C Local Variables:
      INTEGER, SAVE :: LOGDEV       ! unit number for log device

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: MAXSYNC       ! force max TSTEP(2) (sec)
      INTEGER, SAVE :: MINSYNC       ! force min TSTEP(2) (sec)
      REAL,    SAVE :: SIGST         ! sigma_sync_top value
      INTEGER, SAVE :: ADVLAYR       ! adv=sync at least up to this level

      INTEGER, SAVE :: WSTEP = 0     ! wind file interpolation time step
      INTEGER, SAVE :: SDATE, STIME  ! wind file start date and time
      INTEGER, SAVE :: FDATE, FTIME  ! wind file current date and time
      INTEGER, SAVE :: FSTEP         ! wind file time step
      INTEGER, SAVE :: OSTEP         ! output time step (sec)

      REAL,    SAVE :: IDX1, IDX2    ! 1/dx1, 1/dx2
      REAL,    SAVE :: HDIV_LIM = 1.0  ! cutoff for max horizontal divergence step mod

      INTEGER :: MINSECS             ! effective MINSYNC (sec) that divides TSTEP(1)
 
      INTEGER :: C, R, L, K  ! loop counters
      INTEGER :: EDATE, ETIME
      INTEGER :: SYNC, SUB, REP      ! loop variables
      INTEGER :: ADV, TADV           ! test sub adv step
      INTEGER :: NADVS( NLAYS )      ! no. of advection sub-steps per sync step
      REAL    :: RSYNC
      REAL    :: MXUOVDX             ! max ( WIND / DS ) all layers 1, ADVLAYR
      LOGICAL :: ADJFLG

      REAL    :: UWIND( NCOLS+1,NROWS+1,NLAYS )
      REAL    :: VWIND( NCOLS+1,NROWS+1,NLAYS )
      REAL    :: HDIV( NCOLS,NROWS,NLAYS )
                                                    ! MAX ( WIND / DS ) ...
      REAL, ALLOCATABLE, SAVE :: WIND_IDX_BEG( : )  ! at beginning of output step
      REAL, ALLOCATABLE, SAVE :: WIND_IDX_END( : )  ! at end of output step
      REAL, ALLOCATABLE, SAVE :: WIND_IDX    ( : )  ! for whole output step
!     REAL,              SAVE :: WIND_IDX    ( NLAYS ) Can't SAVE an automatic array!
                                                    ! MAX ( HDIV ) ...
      REAL, ALLOCATABLE, SAVE :: MX_HDIV_BEG( : )   ! at beginning of output step
      REAL, ALLOCATABLE, SAVE :: MX_HDIV_END( : )   ! at end of output step
      REAL, ALLOCATABLE, SAVE :: MX_HDIV    ( : )   ! for whole output step

      CHARACTER( 16 ) :: PNAME = 'ADVSTEP'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER :: ALLOCSTAT
      INTEGER :: STATUS         ! ENVINT status

      CHARACTER( 80 ) :: VARDESC  ! environment variable description

C environment variable max sync step
      CHARACTER( 16 ) :: CTM_MAXSYNC = 'CTM_MAXSYNC'
C environment variable min sync step
      CHARACTER( 16 ) :: CTM_MINSYNC = 'CTM_MINSYNC'
C environment variable adv layer
      CHARACTER( 16 ) :: SIGMA_SYNC_TOP = 'SIGMA_SYNC_TOP'
C environment variable max horiz. div.
      CHARACTER( 16 ) :: ADV_HDIV_LIM = 'ADV_HDIV_LIM'

      INTERFACE
        SUBROUTINE HVELOC( FDATE, FTIME, FSTEP, UWIND, VWIND, HDIV )
            USE GRID_CONF
            USE SUBST_MODULES
            IMPLICIT NONE
            INTEGER, INTENT( IN )  :: FDATE, FTIME, FSTEP
            REAL,    INTENT( OUT ) :: UWIND( :,:,: ), VWIND( :,:,: )
            REAL,    INTENT( OUT ) :: HDIV( :,:,: )
         END SUBROUTINE HVELOC
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

         MAXSYNC = 720          ! default
         VARDESC = 'Maximum Synchronization Time Step (sec)'
         MAXSYNC = ENVINT( CTM_MAXSYNC, VARDESC, MAXSYNC, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV,'(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

         MINSYNC = 60           ! default
         VARDESC = 'Minimum Synchronization Time Step (sec)'
         MINSYNC = ENVINT( CTM_MINSYNC, VARDESC, MINSYNC, STATUS )
         IF ( STATUS .NE. 0 ) WRITE(LOGDEV,'(5X, A)') VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

         SIGST = 0.7            ! default
         VARDESC = 'Minimum layer limit for which adv = sync'
         SIGST = ENVREAL( SIGMA_SYNC_TOP, VARDESC, SIGST, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV,'(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

         IF ( VGTYP_GD .NE. VGSGPN3 .AND.
     &        VGTYP_GD .NE. VGSGPH3 .AND.
     &        VGTYP_GD .NE. VGWRFEM ) THEN
            XMSG = 'Wrong vertical coordinate type'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

         IF ( SIGST .GT. VGLVS_GD( 2 ) .OR.
     &        SIGST .LT. VGLVS_GD( NLAYS + 1 ) ) THEN
            WRITE(XMSG, '(A, 1PE12.3)') 'SIGMA_SYNC_TOP incorrect', SIGST
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

         IF ( SIGST .EQ. VGLVS_GD( NLAYS + 1 ) ) THEN
            L = NLAYS
         ELSE
            DO K = 1, NLAYS
               IF ( SIGST .GE. VGLVS_GD( K+1 ) ) THEN
                  IF ( SIGST - VGLVS_GD( K+1 ) .LE.
     &                 VGLVS_GD( K ) - SIGST )  THEN
                     L = K
                  ELSE
                     L = K - 1
                  END IF
                  EXIT
               END IF
            END DO
         END IF

         ADVLAYR = L
         WRITE( LOGDEV,92005 ) ADVLAYR

C Open wind field file and get header data
         IF ( .NOT. OPEN3( MET_DOT_3D, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open ' // MET_DOT_3D // ' file'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT0 )
         END IF

         IF ( .NOT. DESC3( MET_DOT_3D ) ) THEN
            XMSG = 'Could not get ' // MET_DOT_3D // ' file description'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         SDATE = SDATE3D
         STIME = STIME3D
         FSTEP = TSTEP3D

C Check file data against COORD.EXT

         IF ( XCELL3D .NE. XCELL_GD .OR. YCELL3D .NE. YCELL_GD ) THEN
            XMSG = 'File grid sizes do not match CTM domain definition'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
         END IF

C Get cell size in meters (from COORD.EXT)

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            IDX2 = 1.0 / ( DG2M * YCELL_GD )
            IDX1 = 1.0 / ( DG2M * XCELL_GD
     &           * COS( PI180*( YORIG_GD + YCELL_GD*FLOAT( NROWS/2 ) ) ) )
         ELSE
            IDX1 = 1.0 / XCELL_GD
            IDX2 = 1.0 / YCELL_GD
         END IF

C Get the time of starting wind field: do not interpolate
         IF ( .NOT. CURRSTEP( JDATE, JTIME, SDATE, STIME, FSTEP,
     &                        FDATE, FTIME) ) THEN
            XMSG = 'Cannot get step-starting date and time'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
         END IF

         ALLOCATE( WIND_IDX_BEG( NLAYS ),
     &             WIND_IDX_END( NLAYS ),
     &             WIND_IDX    ( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating WIND_IDX*'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE( MX_HDIV_BEG( NLAYS ),
     &             MX_HDIV_END( NLAYS ),
     &             MX_HDIV    ( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating MX_HDIV*'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         CALL HVELOC( FDATE, FTIME, WSTEP, UWIND, VWIND, HDIV )

         DO L = 1, NLAYS
            WIND_IDX_BEG( L ) = 0.0
            MX_HDIV_BEG( L ) = 0.0
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  WIND_IDX_BEG( L ) = MAX( WIND_IDX_BEG( L ),
     &                                     ABS( UWIND( C,R,L ) ) * IDX1,
     &                                     ABS( VWIND( C,R,L ) ) * IDX2 )
                  MX_HDIV_BEG( L ) = MAX( MX_HDIV_BEG( L ), HDIV( C,R,L ) )
               END DO
            END DO
            WIND_IDX_BEG( L ) = SUBST_GLOBAL_MAX( WIND_IDX_BEG( L ) )
            MX_HDIV_BEG( L )  = SUBST_GLOBAL_MAX( MX_HDIV_BEG( L ) )
         END DO

         VARDESC = 'Maximum horiz. div. limit for which adv = sync'
         HDIV_LIM = ENVREAL( ADV_HDIV_LIM, VARDESC, HDIV_LIM, STATUS )

      END IF   ! FIRSTIME

C Compute the least number of equal time steps that satisfy the Courant
C condition (force TSTEP(2) to be no greater than MAXSYNC):

      OSTEP = TIME2SEC( TSTEP( 1 ) )
      MAXSYNC = MIN( OSTEP, MAXSYNC )

C Make sure MINSYNC <= output time step
      IF ( MINSYNC .GT. OSTEP ) THEN
         WRITE( LOGDEV,92009 ) OSTEP, MINSYNC
         XMSG = ' '
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
      END IF

C Make sure MINSYNC <= MAXSYNC
      IF ( MINSYNC .GT. MAXSYNC ) THEN
         WRITE( LOGDEV,92010 ) MAXSYNC, MINSYNC
         XMSG = ' '
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
      END IF

C Determine MINSECS, the smallest time step (sec) between MINSYNC and MAXSYNC
C that divides the output step
      IF ( MINSYNC .LE. 1 ) THEN
         MINSECS = 1
      ELSE
         DO MINSECS = MINSYNC, MAXSYNC
            IF ( MOD( OSTEP, MINSECS ) .EQ. 0 ) THEN
               IF ( MINSECS .NE. MINSYNC )
     &              WRITE( LOGDEV,92011 ) JDATE, JTIME, MINSYNC, MINSECS
               EXIT
            ELSE IF ( MINSECS .LT. MAXSYNC ) THEN
               CYCLE
            END IF
C A satisfactory MINSECS could not be found
            WRITE( XMSG,92012 ) MINSYNC, OSTEP
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
         END DO
      END IF

      WRITE( LOGDEV,92013 ) JDATE, JTIME

C Establish ending time for this Courant number calculation:
      EDATE = JDATE
      ETIME = JTIME
      CALL NEXTIME( EDATE, ETIME, TSTEP( 1 ) )

      IF ( SECSDIFF( FDATE, FTIME, EDATE, ETIME ) .GT. 0 )
C Get the next wind field: do not interpolate
     &   CALL NEXTIME( FDATE, FTIME, FSTEP )

      CALL HVELOC( FDATE, FTIME, WSTEP, UWIND, VWIND, HDIV )

      DO L = 1, NLAYS
         WIND_IDX_END( L ) = 0.0
         MX_HDIV_END( L ) = 0.0
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               WIND_IDX_END( L ) = MAX( WIND_IDX_END( L ),
     &                                  ABS( UWIND( C,R,L ) ) * IDX1,
     &                                  ABS( VWIND( C,R,L ) ) * IDX2 )
               MX_HDIV_END( L ) = MAX( MX_HDIV_END( L ), HDIV( C,R,L ) )
            END DO
         END DO

!        if ( jtime .gt. 075600 ) then
D        if ( jtime .gt. 005600 ) then
D           write( logdev,* ) 'local MX_HDIV_END: ', l, mx_hdiv_end( l )
D        end if

         WIND_IDX_END( L ) = SUBST_GLOBAL_MAX( WIND_IDX_END( L ) )
         MX_HDIV_END( L )  = SUBST_GLOBAL_MAX( MX_HDIV_END( L ) )

!        if ( jtime .gt. 075600 ) then
D        if ( jtime .gt. 005600 ) then
D           write( logdev,* ) 'global MX_HDIV_END: ', l, mx_hdiv_end( l )
D        end if

      END DO

      WIND_IDX = MAX( WIND_IDX_BEG, WIND_IDX_END ) ! array
      WIND_IDX_BEG = WIND_IDX_END                  ! save for next call
      MX_HDIV = MAX( MX_HDIV_BEG, MX_HDIV_END )    ! array
      MX_HDIV_BEG = MX_HDIV_END                    ! save for next call

C Determine the synchronization timestep, TSTEP(2), and NREPS, the corresponding
C number of sync steps per output timestep

C Find the largest timestep that satisfies the Courant condition in all of the
C lower ADVLAYR layers. If this is smaller than MINSECS, then set the sync step
C to MINSECS. Otherwise, set the sync step to the largest time step that both
C satisfies the Courant condition and is an integer divisor of the output step
C (existence of such a sync step is guaranteed because, at this point, MINSECS
C satisfies both conditions).
      MXUOVDX = MAXVAL( WIND_IDX( 1:ADVLAYR ) )
      ADJFLG = .FALSE.
      DO REP = 1, OSTEP
         IF ( MOD( OSTEP, REP ) .EQ. 0 ) THEN   ! make TSTEP(2) divide TSTEP(1)
            SYNC = OSTEP / REP
            IF ( SYNC .LE. MAXSYNC ) THEN   ! enforce max TSTEP(2)
               ADV = SYNC
               IF ( MXUOVDX * REAL( ADV ) .LT. CC ) THEN   ! CC satisfied
                  IF ( SYNC .GE. MINSECS ) THEN   ! force min TSTEP(2)
                     NREPS = REP
                     GO TO 301
                  ELSE   ! multiple adv steps per sync step (= MINSECS)
                     SYNC = MINSECS
                     NREPS = OSTEP / MINSECS
                     DO SUB = 1, MINSECS
                        IF ( MOD( MINSECS, SUB ) .EQ. 0 ) THEN  ! make ADV divide TSTEP(2)
                           TADV = MINSECS / SUB
                           IF ( TADV .LE. ADV ) THEN
                              ADV = TADV
                              ADJFLG = .TRUE.
                              GO TO 301
                           END IF
                        END IF
                     END DO
                     GO TO 291   ! No solution
                  END IF
               END IF    ! if Courant condition satisfied
            END IF    ! if SYNC .le. MAXSYNC
         END IF    ! if REP divides STEP evenly
      END DO    ! REP loop

291   continue

C If you get here: could not determine satisfactory advection time step.

      WRITE( LOGDEV,94011 ) TSTEP( 1 ), MXUOVDX
      XMSG = ' '
      CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )

301   CONTINUE

      TSTEP( 2 ) = SEC2TIME( SYNC )
      DO L = 1, ADVLAYR
         ASTEP( L ) = SEC2TIME( ADV )
         NADVS( L ) = SYNC / ADV
      END DO

      IF ( ADJFLG ) THEN
         WRITE( LOGDEV,92019 ) MINSECS, TSTEP( 2 ), NREPS
         ELSE
         WRITE( LOGDEV,92021 ) TSTEP( 2 ), NREPS
      END IF

      DO L = ADVLAYR + 1, NLAYS
         DO REP = 1, SYNC
            IF ( MOD( SYNC, REP ) .EQ. 0 ) THEN
               ADV = SYNC / REP
               IF ( WIND_IDX( L ) * REAL( ADV ) .LT. CC ) THEN  ! Layer CC satisfied
                  NADVS( L ) = REP
                  ASTEP( L ) = SEC2TIME( ADV )
                  GO TO 401
               END IF   ! if courant condition satisfied
            END IF   !  if rep divides step
         END DO

C If you get here: could not determine satisfactory advection time step.

         WRITE( LOGDEV,94013 ) TSTEP( 1 ), L, WIND_IDX( L )
         XMSG = ' '
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )

401      CONTINUE

      END DO   ! layer loop

C For each layer ensure that Max(HDIV)*ADV < 1.

      WRITE( LOGDEV,92025 )
      DO L = NLAYS, 1, -1   ! start at the top, purely for the 92027, 92029 write
         ADV = TIME2SEC( ASTEP( L ) )
         IF ( REAL( ADV ) * MX_HDIV( L ) .LT. HDIV_LIM ) THEN
            WRITE( LOGDEV,92027 ) L, ASTEP( L ), NADVS( L )
            CYCLE
         ELSE
            SUB = ADV
            DO WHILE ( SUB .GT. 0 )
               SUB = SUB - 1
               IF ( MOD( SYNC, SUB ) .EQ. 0 ) THEN
                  ASTEP( L ) = SEC2TIME( SUB )
                  NADVS( L ) = SYNC / SUB
                  GO TO 501
               END IF
            END DO
C If you get here: could not determine satisfactory advection time step.

            WRITE( LOGDEV,94013 ) TSTEP( 1 ), L, WIND_IDX( L )
            XMSG = ' '
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )

501         CONTINUE
            WRITE( LOGDEV,92029 ) L, ASTEP( L ), NADVS( L )
         END IF
      END DO

      RETURN

C------------------  FORMAT  STATEMENTS   ------------------------------

92005 FORMAT(   5X, 'Top layer thru which sync step determined:', I3 )

92009 FORMAT( / 5X, 'Output time step:', I8,
     &          1X, 'less than minimum synchronization step:', I8 )

92010 FORMAT( / 5X, 'Maximum synchronization step:', I8,
     &          1X, 'less than minimum synchronization step:', I8 )

92011 FORMAT( / 5X, 'From ADVSTEP - date/time: ', I8, '/', I6.6
     &        / 5X, 'Minimum Synchronization Step adjusted from:', I8,
     &          1X, 'to:', I8 )

92012 FORMAT( / 5X, 'Starting from:', I6, ',',
     &          1X, 'could not determine minimum step that divides TSTEP',
     &          1X, 'for model step:', I7.6, ' HHMMSS' )

92013 FORMAT( / 5X, 'From ADVSTEP - date/time: ', I8, '/', I6.6 )

92019 FORMAT( / 5X, 'Synchronization step adjusted up to mimimum (SEC):', I7
     &        /46X, '(HHMMSS): ', I6.6
     &        / 5X, 'Number of Synchronization steps:', I5 / )

92021 FORMAT( / 5X, 'Computed synchronization step (HHMMSS): ', I6.6
     &        / 5X, 'Number of Synchronization steps:', I5 / )

92025 FORMAT( / 5X, 'Layer', 3X, 'Advection', 3X, 'per Sync'
     &        /11X, 'Step (HHMMSS)', 2X, 'Step' )

92027 FORMAT(   5X, I4, 6X, I6.6, 6X, I2 )
92029 FORMAT(   5X, I4, 6X, I6.6, 6X, I2, 2X, '<- Hdiv adjusted' )

94011 FORMAT( / 5X, 'Could not determine Courant-condition safe sync step',
     &          1X, 'for model step:', I7.6, ' HHMMSS',
     &        / 5X, '(Max vel)/(dX) =', 1PE10.3)

94013 FORMAT( / 5X, 'Could not determine Courant-condition safe sync step',
     &          1X, 'for model step:', I7.6, ' HHMMSS',
     &          1X, 'in layer:', I3
     &        / 5X, '(Max vel)/(dX) =', 1PE10.3)

      END
