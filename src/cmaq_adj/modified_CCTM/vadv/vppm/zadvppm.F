
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vadv/vppm/zadvppm.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ZADV ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Advection in the vertical, x3-direction:
C   The process time step is set equal to TSTEP
 
C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times. Process time step divides TSTEP exactly
C   CGRID in transport units: SQRT{DET[metric tensor]}*concentration (Mass/Vol)
 
C Subroutines and functions called:
C   TIME2SEC
 
C Revision history:
C   02/19/93 by M. Talat Odman  at NCSC
C   05/17/93 by Carlie J. Coats at NCSC:  now uses INTERP3()
C   06/14/94 by Dongming Hwang at NCSC: 
C              include statement and subroutine name template
C   10/15/95 by M. Talat Odman at NCSC: generalized coordinates

C   Sep 97 Jeff
C   Aug 98 Jeff better Courant condition tstep limit

C    David Wong, Sep. 1998
C      -- parallelized the code

C    15 Dec 00 J.Young: move CGRID_MAP into f90 module
C                       GLOBAL_RSUM -> Dave Wong's f90 stenex GLOBAL_SUM
C                       GLOBAL_ISUM -> Dave Wong's f90 stenex GLOBAL_SUM

C    28 Jul 01 J.Young: allocatable arrays ...
C                       Since F90 does not preserve dummy argument array
C                       indices, the 3rd dimension of WHAT has been changed
C                       from 0:NLAYS to 1:NLAYS+1 for the sake of vcontvel

C    03 Sep 01 David Wong
C      -- inserted F90 DEALLOCATE statement for NX3

C    31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE SUBST_MODULES              ! stenex
!     USE SUBST_GLOBAL_SUM_MODULE    ! stenex

      IMPLICIT NONE

C Includes:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID    ! coordinate & domain definitions (req IOPARMS)
 
C Arguments:
 
!     REAL        CGRID ( NCOLS,NROWS,NLAYS,*) ! concentrations
!     REAL     :: CGRID( :,:,:,: )
      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER     JDATE         ! current model date, coded YYYYDDD
      INTEGER     JTIME         ! current model time, coded HHMMSS
      INTEGER     TSTEP( 2 )    ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
 
C External Functions not declared in IODECL3.EXT:
 
      INTEGER, EXTERNAL :: TIME2SEC
 
C Parameters:
 
C Advected species dimension
 
      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV
     &                                + N_AE_ADV
     &                                + N_NR_ADV
     &                                + N_TR_ADV
     &                                + 1       ! for advecting air

C File Variables:

!     REAL        WHAT( NCOLS,NROWS,0:NLAYS )   ! x3-component CX-velocity
      REAL, ALLOCATABLE :: WHAT( :,:,: )        ! x3-component CX-velocity

C Local variables:
 
      CHARACTER( 16 ) :: PNAME = 'ZADVPPM'
 
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      REAL          CON( NLAYS )         ! Concentrations in a vertical column
      REAL          VEL( NLAYS )         ! Velocities in a vertical column
!     REAL, SAVE :: DS ( NLAYS )         ! dx3 (dimensionless in sigma coord.)
      REAL, ALLOCATABLE, SAVE :: DS ( : ) ! dx3 (dimensionless in sigma coord.)
      REAL          DTSEC                ! model time step in seconds
 
      INTEGER, SAVE :: ADV_MAP( N_SPC_ADV ) ! global adv map to CGRID

      INTEGER      COL, ROW, LVL, SPC, VAR   ! loop counters
      INTEGER      A2C

      CHARACTER( 96 ) :: XMSG = ' '

C for Courant stability
 
      REAL, PARAMETER :: CFL = 0.75     ! maximum allowable Courant number
!     REAL        DTX3( NCOLS,NROWS )   ! Courant safe time step
      REAL, ALLOCATABLE, SAVE :: DTX3( :,: )  ! Courant safe time step
      REAL        FF                    ! vertical velocity component / DS
!     INTEGER     NX3( NCOLS,NROWS )    ! Number of iterations per TSTEP(2)
      INTEGER, ALLOCATABLE, SAVE :: NX3( :,: ) ! Number of iterations per TSTEP(2)
      INTEGER     IX3, VSTEP            ! loop and intermediate variables
      INTEGER     ALLOCSTAT
 
!     INTEGER, SAVE :: LOGDEV

!     INTEGER     TNX3                  ! stats acummulator
!     REAL        TDTX3                 ! stats acummulator
!     INTEGER     MXNX3                 ! max no. of integration steps
!     INTEGER     MNNX3                 ! min no. of integration steps
!     REAL        MXDTX3                ! max integration time step
!     REAL        MNDTX3                ! min integration time step
!     REAL        ANX3                  ! grid average no. of integration steps
!     REAL        ADTX3                 ! grid average integration time step

C Required interface for allocatable array dummy arguments

      INTERFACE
         SUBROUTINE VCONTVEL ( FDATE, FTIME, TSTEP, WHAT )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: FDATE, FTIME, TSTEP
            REAL, INTENT( OUT )         :: WHAT( :,:,: )
         END SUBROUTINE VCONTVEL
      END INTERFACE
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

!        LOGDEV = INIT3 ()

C Get dx3 from COORD include file

         ALLOCATE ( DS( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DS'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         DO LVL = 1, NLAYS
            DS ( LVL ) = ABS ( X3FACE_GD( LVL ) - X3FACE_GD( LVL-1 ) )
            END DO

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
 
C Create global map to CGRID
 
         SPC = 0
         DO VAR = 1, N_GC_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = GC_STRT - 1 + GC_ADV_MAP( VAR )
            END DO
         DO VAR = 1, N_AE_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = AE_STRT - 1 + AE_ADV_MAP( VAR )
            END DO
         DO VAR = 1, N_NR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = NR_STRT - 1 + NR_ADV_MAP( VAR )
            END DO
         DO VAR = 1, N_TR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = TR_STRT - 1 + TR_ADV_MAP( VAR )
            END DO
 
         ADV_MAP( N_SPC_ADV ) = N_GC_SPCD

         ALLOCATE ( DTX3( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DTX3'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         ALLOCATE ( NX3( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating NX3'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         END IF                    ! if firstime

      ALLOCATE ( WHAT( MY_NCOLS,MY_NROWS,1:NLAYS+1 ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'Failure allocating WHAT'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Get the contravarient vertical velocity d(x3)/dt for the middle
C of sync time step, TSTEP(2)

      CALL VCONTVEL ( JDATE, JTIME, TSTEP( 2 ), WHAT )

C Process time step is TSTEP(2).
C Time-stepped gridded computation for Z-direction advection.

      DTSEC  = FLOAT( TIME2SEC( TSTEP( 2 ) ) ) ! process time step (seconds)

C vertical velocities are at face centers, positive upward.
C No boundary conditions are needed because VEL(1) = VEL(NLAYS+1) = 0

!     TNX3 = 0
!     TDTX3 = 0.0
!     MXNX3 = -1
!     MNNX3 = STEP + 1
!     MXDTX3 = 0.0
!     MNDTX3 = FLOAT( STEP ) + 1

      DO ROW = 1, MY_NROWS
         DO COL = 1, MY_NCOLS
 
            DO LVL = 1, NLAYS
               VEL( LVL ) = WHAT( COL,ROW,LVL )
               END DO
 
            FF = MAX (
     &               ABS( VEL( 1 ) / DS( 1 ) ),
     &               ABS( VEL( 2 ) / DS( 1 ) )
     &               )
            FF = MAX ( FF,
     &               ABS( VEL( NLAYS-1 ) / DS( NLAYS ) ),
     &               ABS( VEL( NLAYS )   / DS( NLAYS ) )
     &               )
            DO LVL = 2, NLAYS-1
               FF = MAX ( FF,
     &                  ABS( VEL( LVL-1 ) / DS( LVL ) ),
     &                  ABS( VEL( LVL )   / DS( LVL ) ),
     &                  ABS( VEL( LVL+1 ) / DS( LVL ) )
     &                  )
               END DO
 
            DO IX3 = 1, DTSEC
               VSTEP = DTSEC / IX3
               IF ( FF * FLOAT( VSTEP ) .LT. CFL ) THEN
                  NX3( COL,ROW ) = IX3
                  DTX3( COL,ROW ) = VSTEP
                  GO TO 111
                  END IF
               END DO

C unable to get vertical advection step
            WRITE( XMSG,94010 ) TSTEP( 2 ), FF
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
94010       FORMAT( 'Could not compute vertical advection step for MODEL STEP:',
     &               I7, ' HHMMSS', 2X,
     &              'Max vel / dS =', 1PE10.3)

111         CONTINUE

C accumulate step stats
!           TNX3 = TNX3 + NX3( COL,ROW )
!           TDTX3 = TDTX3 + DTX3( COL,ROW )
!           MXNX3 = MAX( NX3( COL,ROW ), MXNX3 )
!           MNNX3 = MIN( NX3( COL,ROW ), MNNX3 )
!           MXDTX3 = MAX( DTX3( COL,ROW ), MXDTX3)
!           MNDTX3 = MIN( DTX3( COL,ROW ), MNDTX3)
 
            END DO
         END DO

      DO SPC = 1, N_SPC_ADV
         A2C = ADV_MAP(SPC)

         DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS

               DO LVL = 1, NLAYS
                  VEL( LVL ) = WHAT( COL,ROW,LVL )
                  CON( LVL ) = CGRID( COL,ROW,LVL,A2C )
                  END DO

               DO IX3 = 1, NX3( COL,ROW )
                  CALL VPPM ( NLAYS, CON, VEL, DTX3( COL,ROW ), DS )
                  END DO

               DO LVL = 1, NLAYS
                  CGRID( COL,ROW,LVL,A2C ) = CON( LVL )
                  END DO

               END DO
            END DO

         END DO

      DEALLOCATE ( WHAT )

!     ANX3 = FLOAT( SUBST_GLOBAL_SUM( TNX3 ) ) / FLOAT( GL_NCOLS*GL_NROWS )
!     ADTX3 = SUBST_GLOBAL_SUM( TDTX3 ) / FLOAT( GL_NCOLS*GL_NROWS )
!     WRITE( LOGDEV,92010 ) ANX3, MXNX3, MNNX3,
!    &                      ADTX3, MXDTX3, MNDTX3
92010 FORMAT ( / 5X, 'Grid Avg, Max, Min - No. of vert. integration steps:',
     &         F7.2, 2(I4)
     &         / 5X, 'Grid Avg, Max, Min - integration time step:',
     &         3(1PE10.3) )

      RETURN
      END
