
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hdiff/multiscale_adj_bwd/rho_j.F,v 1.1 2010/11/17 21:06:36 sergey Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE RHO_J ( JDATE, JTIME, TSTEP, RHOJI )
      
C-----------------------------------------------------------------------
C Function:
C   Get Air Density X SqRDMT for computational grid including boundary
C   where SqRDMT = Sq. Root [det ( metric tensor )]
C   = Jacobian / (map scale factor)**2
 
C Preconditions:
C   Can be used only for conformal map coordinates in the horizontal.
C   Dates and times should be represented YYYYDDD:HHMMSS.
 
C Subroutines and functions called:
C    INTERP3, M3EXIT, TIME2SEC, SEC2TIME, NEXTIME, CHECK3
      
C Revision history:
C    5 Nov 97 Jeff

C    Sep. 1998 David Wong
C      -- parallelized the code
C      -- removed the indirect index reference

C    6 Aug 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with
C                      INTERPX and INTERPB; allocatable arrays ...
C                      Since F90 does not preserve dummy argument array
C                      indices, RHOJI( 1:NCOLS+2,, ) is copied into local array
C                      RHOJI( 0:NCOLS+1,, ).
C                      The caller of RHO_J dimensions the actual argument,
C                      as RHOJ( 0:NCOLS+1,, ).
C                      NOTE: NTHIK must be = 1
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------
      
      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Includes:
      
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! horizontal dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structuer
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
 
C Parameters

C Arguments:
      
      INTEGER     JDATE        ! current model date, coded YYYYDDD
      INTEGER     JTIME        ! current model time, coded HHMMSS
      INTEGER     TSTEP( 2 )   ! time step vector (HHMMSS)
                               ! TSTEP(1) = local output step
                               ! TSTEP(2) = sciproc sync. step (chem)
!     REAL        RHOJ( 0:NCOLS+1,0:NROWS+1,NLAYS )  ! Jacobian * air density
      REAL     :: RHOJI( :,:,: )
      
C Parameters:

!     INTEGER, PARAMETER :: NTHIN = NTHIK - 1

C file variables:
      
!     REAL        DENSJ_BUF( NCOLS,NROWS,NLAYS ) ! Jacobian * air density
!     REAL        DENSJ_BND( NBNDY,NLAYS )    ! boundary Jacobian * air density
      REAL, ALLOCATABLE :: DENSJ_BUF( :,:,: ) ! Jacobian * air density
      REAL, ALLOCATABLE :: DENSJ_BND( :,: )    ! boundary Jacobian * air density
 
C External Functions (not already declared by IODECL3.EXT):
 
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, TRIMLEN
#ifdef parallel
      LOGICAL, EXTERNAL :: INTERPB
#endif
      
C local variables:
      
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 16 ) :: PNAME = 'RHO_J'
      CHARACTER( 16 ) :: AMSG
      CHARACTER( 96 ) :: XMSG = ' '
 
      REAL, ALLOCATABLE :: RHOJ( :,:,: )

      INTEGER, SAVE :: MLAYS

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
       
      INTEGER     ALLOCSTAT
      INTEGER     COL, ROW, LVL ! column, row, layer indices
      INTEGER     MDATE         ! mid-advection date
      INTEGER     MTIME         ! mid-advection time
      INTEGER     STEP          ! advection time step in seconds
      INTEGER, SAVE :: LDATE( 2 )    ! last date for data on file
      INTEGER, SAVE :: LTIME( 2 )    ! last time for data on file
      LOGICAL     REVERT        ! recover last time step if true
 
      INTEGER COUNT

      INTEGER   GXOFF, GYOFF            ! global origin offset from file
      LOGICAL, SAVE :: WINDOW = .FALSE. ! posit same file and global
                                        ! processing domain
      INTEGER, SAVE :: NCOLSDENS, NROWSDENS       ! local for DENSJ_BUF
C for INTERPX
      INTEGER, SAVE :: STRTCOL,   ENDCOL,   STRTROW,   ENDROW
      INTEGER       :: STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC

      INTEGER, SAVE :: LCB, HCB, LRB, HRB ! convenience pointers
 
C-----------------------------------------------------------------------
 
      IF ( FIRSTIME ) THEN
 
         FIRSTIME = .FALSE.
 
         LCB = 1 - NTHIK
         HCB = MY_NCOLS + NTHIK
         LRB = 1 - NTHIK
         HRB = MY_NROWS + NTHIK
         MLAYS = SIZE ( RHOJI,3 )

         CALL LSTEPF( MET_CRO_3D, LDATE( 1 ), LTIME( 1 ) )
         CALL LSTEPF( MET_BDY_3D, LDATE( 2 ), LTIME( 2 ) )

         LDATE( 1 ) = MIN( LDATE( 1 ), LDATE( 2 ) )
         LTIME( 1 ) = SEC2TIME( MIN(
     &                              TIME2SEC( LTIME( 1 ) ),
     &                              TIME2SEC( LTIME( 2 ) )
     &                              ) )

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC )

         NCOLSDENS = ENDCOLMC - STRTCOLMC + 1
         NROWSDENS = ENDROWMC - STRTROWMC + 1

         IF ( NCOLSDENS .NE. MY_NCOLS .OR.
     &        NROWSDENS .NE. MY_NROWS ) THEN
            WRITE( XMSG,'( A, 4I8 )' ) 'Local Columns or Rows incorrect',
     &      NCOLSDENS, MY_NCOLS, NROWSDENS, MY_NROWS
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

C currently not implemented: case where only one origin component matches file's
         IF ( GXOFF .NE. 0 .AND. GYOFF .NE. 0 ) THEN
            WINDOW = .TRUE.       ! windowing from file
            STRTCOL = STRTCOLMC - NTHIK
            ENDCOL  = ENDCOLMC  + NTHIK
            STRTROW = STRTROWMC - NTHIK
            ENDROW  = ENDROWMC  + NTHIK
            ELSE
            STRTCOL = STRTCOLMC
            ENDCOL  = ENDCOLMC
            STRTROW = STRTROWMC
            ENDROW  = ENDROWMC
            END IF

         END IF                    ! if firstime
 
      MDATE  = JDATE
      MTIME  = JTIME
      STEP   = TIME2SEC( TSTEP( 2 ) )
      CALL NEXTIME( MDATE, MTIME, SEC2TIME( STEP / 2 ) )

      IF ( MDATE .LT. LDATE( 1 ) ) THEN
         REVERT = .FALSE.
         ELSE IF ( MDATE .EQ. LDATE( 1 ) ) THEN
            IF ( MTIME .LE. LTIME( 1 ) ) THEN
               REVERT = .FALSE.
               ELSE
               REVERT = .TRUE.
               END IF
         ELSE   ! MDATE .GT. LDATE
         REVERT = .TRUE.
         END IF

      IF ( REVERT ) THEN
         XMSG = 'Current scenario interpolation step not available in both '
     &        // MET_CRO_3D(1:TRIMLEN( MET_CRO_3D ) ) // 'and '
     &        // MET_BDY_3D(1:TRIMLEN( MET_BDY_3D ) )
         CALL M3MESG( XMSG )
!        CALL NEXTIME( MDATE, MTIME, -SEC2TIME( STEP / 2 ) )
         WRITE( AMSG,'( 2I8 )' ) LDATE( 1 ), LTIME( 1 )
         XMSG = 'Using data for last file step: ' // AMSG
         CALL M3MESG( XMSG )
         MDATE = LDATE( 1 )
         MTIME = LTIME( 1 )
         END IF
 
C Interpolate Jacobian X Air Density

      ALLOCATE ( RHOJ( LCB:HCB,LRB:HRB,MLAYS ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'Failure allocating RHOJ'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      RHOJ = -9999999.9

      IF ( WINDOW ) THEN

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                       STRTCOL,ENDCOL, STRTROW,ENDROW, 1,MLAYS,
     &                       MDATE, MTIME, RHOJ ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         ELSE ! need to extend data from bndy file

         ALLOCATE ( DENSJ_BUF( NCOLSDENS,NROWSDENS,NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DENSJ_BUF'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                       STRTCOL,ENDCOL, STRTROW,ENDROW, 1,MLAYS,
     &                       MDATE, MTIME, DENSJ_BUF ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         ALLOCATE ( DENSJ_BND( NBNDY,MLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DENSJ_BND'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPB ( MET_BDY_3D, VNAME, PNAME,
     &                        MDATE, MTIME, NBNDY*MLAYS,
     &                        DENSJ_BND ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_BDY_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

C Load core of RHOJ array

         DO LVL = 1, MLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  RHOJ( COL,ROW,LVL ) = DENSJ_BUF( COL,ROW,LVL )
                  END DO
               END DO
            END DO

C Fill in DENSJ array for boundaries

         DO LVL = 1, MLAYS
            COUNT = 0
            DO ROW = 1-NTHIK, 0
               DO COL = 1, MY_NCOLS+NTHIK
                  COUNT = COUNT + 1
                  RHOJ( COL,ROW,LVL ) = DENSJ_BND( COUNT,LVL )  ! South
               END DO
            END DO
            DO ROW = 1, MY_NROWS+NTHIK
               DO COL = MY_NCOLS+1, MY_NCOLS+NTHIK
                  COUNT = COUNT + 1
                  RHOJ( COL,ROW,LVL ) = DENSJ_BND( COUNT,LVL )  ! East
               END DO
            END DO
            DO ROW = MY_NROWS+1, MY_NROWS+NTHIK
               DO COL = 1-NTHIK, MY_NCOLS
                  COUNT = COUNT + 1
                  RHOJ( COL,ROW,LVL ) = DENSJ_BND( COUNT,LVL )  ! North
               END DO
            END DO
            DO ROW = 1-NTHIK, MY_NROWS
               DO COL = 1-NTHIK, 0
                  COUNT = COUNT + 1
                  RHOJ( COL,ROW,LVL ) = DENSJ_BND( COUNT,LVL )  ! West
               END DO
            END DO
         END DO

         DEALLOCATE ( DENSJ_BUF )
         DEALLOCATE ( DENSJ_BND )

         END IF   ! WINDOW

C Adjust for argument offset (f90 dummy arguments assumed 1-based)

      DO LVL = 1, MLAYS
         DO ROW = 0, NROWS + 1
            DO COL = 0, NCOLS + 1
               RHOJI( COL+1,ROW+1,LVL ) = RHOJ( COL,ROW,LVL )
               END DO
            END DO
         END DO

      DEALLOCATE ( RHOJ )

      RETURN
      END
