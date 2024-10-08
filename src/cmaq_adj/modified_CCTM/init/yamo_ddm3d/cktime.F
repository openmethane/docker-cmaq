
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/init/yamo_ddm3d/cktime.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE CKTIME(IDATE,JDATE,STDATE,ITIME,JTIME,TIMEFLAG,P )

C-----------------------------------------------------------------------
C Function:
C   CKTIME determines whether the current model time, JTIME, is
C   within the time range ITIME specified in sensinput.dat.
C   Returns TIMEFLAG as TRUE or FALSE.
C   Called by advection and diffusion routines in sensitivity analysis.
C Revision History:
C   Created July 02 Daniel Cohan
C 15th March 2006: SRTonse, LBL, if time period straddles midnight
C then both dates need to be set as true and only the hours connected
C to the straddling period are set with TIMEFLAG TRUE
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
 
C   list of active dates IDATE(Sens Par #,NDAYS). 1 if active 0 if not
      INTEGER,POINTER:: IDATE ( :,: )
      INTEGER,INTENT(IN):: JDATE          !current date
      INTEGER,INTENT(IN):: STDATE       ! starting date,    format YYYYDDD
      INTEGER,INTENT(IN)::  JTIME       ! current model time, coded HHMMSS
C  begin & end times ITIME(Sens Par #,1:2) 1 is beg, 2 is end
C  default is 0, 240000
C  Note: to span midnight, begin time may be larger than endtime
      INTEGER,POINTER:: ITIME (:,: )
      INTEGER,INTENT(IN)::  P   !index of current sens par in calling routine
      LOGICAL,INTENT(OUT) :: TIMEFLAG
C local variables
      INTEGER:: I, NAD
      LOGICAL,SAVE :: FIRSTTIME = .TRUE.
C flollowing block of vars for polling the Environment variables 
C to calculate NDAYS
      INTEGER      LOGDEV       ! FORTRAN unit number for log file
C  environment variable start time
      CHARACTER( 16 ) :: CTM_STTIME = 'CTM_STTIME'
C  environment variable run duration
      CHARACTER( 16 ) :: CTM_RUNLEN = 'CTM_RUNLEN'
      CHARACTER( 16 ) :: PNAME = 'CKTIME'
      INTEGER RUNLEN, BGNTIME
      INTEGER, EXTERNAL :: ENVINT       !  get environment variable as integer
      INTEGER      STATUS               !  ENVINT status
      CHARACTER( 80 ) :: VARDESC
      CHARACTER( 96 ) :: XMSG = ' '
      INTEGER,SAVE::  NDAYS


C-----------------------------------------------------------------------
      IF_FIRST: IF(FIRSTTIME)THEN
         FIRSTTIME = .FALSE.
         CALL M3MESG(PNAME//' first entry')
C Use the begin time and run length to determine how many
C dates are in the model run, for use in dimensioning IDATE
         BGNTIME = 000000       ! default
         VARDESC = 'Scenario Starting Time (HHMMSS)'
         BGNTIME = ENVINT( CTM_STTIME, VARDESC, BGNTIME, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, BGNTIME, XMSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Envirnmt variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
         END IF

         RUNLEN = 480000        ! default
         VARDESC = 'Scenario Run Duration (HHMMSS)'
         RUNLEN = ENVINT( CTM_RUNLEN, VARDESC, RUNLEN, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, BGNTIME, XMSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Envirnmt variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
         END IF

         NDAYS = 1 + FLOOR( REAL( BGNTIME + RUNLEN ) / 240000.0 )
      END IF IF_FIRST

      TIMEFLAG = .FALSE.

C time period does not straddle midnight
      IF ( ( JTIME .GE. ITIME ( P,1 )) .AND. (JTIME .LE. ITIME(P,2 ))) THEN
        TIMEFLAG = .TRUE.

C time period straddles midnight
      ELSE IF ( ITIME ( P,1 ) .GT. ITIME ( P,2 ) ) THEN
C count number of active dates set
        NAD = 0
        DO I = 1,NDAYS
           IF(IDATE(P,I) .EQ. 1)NAD = NAD + 1
        END DO
C if all dates set true do not bother to test on date, test on time only
        IF(NAD .EQ. NDAYS)THEN
           TIMEFLAG = 
     &     (JTIME .GE. ITIME(P,1) .AND. JTIME .LE. 240000) .OR. 
     &     (JTIME .LE. ITIME(P,2) .AND. JTIME .GE. 0 ) 
        ELSE
C test on date, and only pick times that straddle 2 active dates
           I =  1 + JDATE - STDATE     !date index in IDATE array
           TIMEFLAG =
     &     (JTIME.GE.ITIME(P,1) .AND. JTIME .LE. 240000
     &     .AND. IDATE(P,I).EQ.1 .AND. IDATE(P,I+1).EQ.1)
     &     .OR.
     &     (JTIME.LE.ITIME(P,2) .AND. JTIME .GE. 0
     &     .AND. IDATE(P,I-1).EQ.1 .AND. IDATE(P,I).EQ.1)
        END IF
      END IF

      RETURN
      END
