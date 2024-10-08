
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/adjcon/denrate/adjadv.F,v 1.1.1.1 2010/06/14 16:02:58 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ADJADV ( CGRID, JDATE, JTIME, TSTEP )
 
C-----------------------------------------------------------------------
C Function:
C   Adjusts concentration fields for mass balance inconsistencies
C   LOAD_RHOJ does the initial loading of Air Density X SqRDMT into RHOJ
C   Where SqRDMT = Sq. Root [det ( metric tensor )]
C                = Jacobian / (map scale factor)**2
C   ADJADV adjusts CGRID(*,*,*,1:NSPCS)
C Adjustment flag:
C  IADJ = 0  no adjustment
C  IADJ = 1  UAM method (calibration gas)
C  IADJ = 2  SAQM method (density adj)
C  IADJ = 3  anelastic approximation
C  IADJ = 4  CMAQ method #1 (used in "semi-demo" testing)
C  IADJ = 5  CMAQ method #2
C  Default   CMAQ method #3 - Byun's method (the best from tracer tests)
     
C Preconditions:
C   CGRID slot for advected RHOJ is at the end of the gas chem sub-array,
C   pointed to by the N_GC_SPCD parameter in GC_SPC table
C   LOAD_RHOJ is called initially to load RHOJ into CGRID
C   ADJADV   is called after each (total, x-y-z) advection step.
      
C Subroutines and Functions Called:
C   SEC2TIME, TIME2SEC, INIT3
      
C Revision History:
C   March 07, 1997 by D.W. Byun & Jeff Young
C   Adapted for Models-3 Interim hydrostatic version
C   Jeff - 19 Aug 97 for targetted
C   Daewon - Dec 97 - alternative methods
C   Jeff - 4 Feb 98 - deal with end of scenario for Default method
C   Jeff - 12 Feb 98 - deal with sync step = chem step
C   2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   Jeff - Dec 00 - move CGRID_MAP into f90 module
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX
C   13 Aug 03 J.Young: back to sync step from adv step
C   07 Dec 04 J.Young: dyn alloc - Use VGRD_DEFN
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   23 Jun 05 J.Young: remove LOAD_RHOJ entry - now done in LOAD_CGRID
C-----------------------------------------------------------------------
      
      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE

C Include Files:
      
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_FILES_ID    ! file name parameters
!     INCLUDE SUBST_COORD_ID    ! coord. and domain definitions (req IOPARMS)
 
C Arguments:
 
!     REAL        CGRID( NCOLS,NROWS,NLAYS,* ) ! concentrations
!     REAL     :: CGRID( :,:,:,: )
      REAL, POINTER :: CGRID( :,:,:,: )

      INTEGER     JDATE        ! current model date , coded YYYYDDD
      INTEGER     JTIME        ! current model time , coded HHMMSS
      INTEGER     TSTEP( 2 )   ! time step vector (HHMMSS)
                               ! TSTEP(1) = local output step
                               ! TSTEP(2) = sciproc sync. step (chem)

C External Functions:
 
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, TRIMLEN
       
C Layered gridded variables (cross)
      
      REAL         RHOJ  ( NCOLS,NROWS,NLAYS )   ! RhoJ (Kg/m**3)
      REAL         RHOJ1 ( NCOLS,NROWS,NLAYS )   ! RhOJ at JTIME
      REAL         ADJFAC( NCOLS,NROWS,NLAYS )   ! adjustment factor
      REAL         DENS  ( NCOLS,NROWS,NLAYS )   ! air density
      
C Parameters:

      REAL, PARAMETER :: TESTVAL = 1.5

      INTEGER, PARAMETER :: IADJ = 6

C local variables:

      REAL         DRATE
C  environment variable for adjcon option
      CHARACTER( 16 ) :: CTM_ADJOPT = 'CTM_ADJOPT'
      CHARACTER( 28 ) :: MSCA = ' mass consistency adjustment'
      
      INTEGER      COL, ROW, LVL, SPC        ! loop counters
      INTEGER      MDATE                     ! mid-advection date
      INTEGER      MTIME                     ! mid-advection time
      INTEGER      STEP                      ! advection time step in seconds
      INTEGER, SAVE :: LDATE                     ! last date for data on file
      INTEGER, SAVE :: LTIME                     ! last time for data on file
      INTEGER      ASPC                      ! species pointer to adjust CGRID
      LOGICAL      REVERT                    ! recover last time step if true
      LOGICAL      FRSTWRT                   ! debug write control
 
      INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER( 16 ) :: PNAME = 'ADJADV'
      CHARACTER( 16 ) :: AMSG
      CHARACTER( 96 ) :: XMSG = ' '

C-----------------------------------------------------------------------
 
      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

         IF ( IADJ .EQ. 0 ) THEN

            WRITE( LOGDEV,* ) '    No' // MSCA

            ELSE IF ( IADJ .EQ. 1 ) THEN

            WRITE( LOGDEV,* ) '    UAM-type calibration gas' // MSCA

            ELSE IF ( IADJ .EQ. 2 ) THEN

            WRITE( LOGDEV,* ) '    SAQM-type density' // MSCA

            ELSE IF ( IADJ .EQ. 3 ) THEN

            WRITE( LOGDEV,* ) '    Anelastic approximation' // MSCA

            CALL LSTEPF( MET_CRO_3D, LDATE, LTIME )
 
            ELSE IF ( IADJ .EQ. 4 ) THEN

            WRITE( LOGDEV,* ) '    CMAQ method 1' // MSCA

            CALL LSTEPF( MET_CRO_3D, LDATE, LTIME )
 
            ELSE IF ( IADJ .EQ. 5 ) THEN

            WRITE( LOGDEV,* ) '    CMAQ method 2' // MSCA

            CALL LSTEPF( MET_CRO_3D, LDATE, LTIME )
 
            ELSE

            WRITE( LOGDEV,* ) '    CMAQ method 3' // MSCA

            CALL LSTEPF( MET_CRO_3D, LDATE, LTIME )
 
            END IF

C Get CGRID offsets

         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C get file local domain offsets

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         END IF    ! FIRSTIME

      IF ( IADJ .EQ. 3 .OR. IADJ .EQ. 4 ) THEN

C Interpolate rho*J at start of step

         IF ( .NOT. INTERPX( MET_CRO_3D, 'DENSA_J', PNAME,
!    &              1,NCOLS, 1,NROWS, 1,NLAYS,
     &              STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &              JDATE, JTIME, RHOJ1 ) ) THEN
            XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

c        CALL DJINTERP ( JDATE, JTIME, TSTEP(2), RHOJ1 ) ! for hydrostatic

         END IF

      IF ( IADJ .GE. 3 ) THEN

C Interpolate rho*J at step
 
         MDATE  = JDATE
         MTIME  = JTIME
         STEP   = TIME2SEC( TSTEP( 2 ) )
!        STEP   = TIME2SEC( TSTEP( 3 ) )
         CALL NEXTIME( MDATE, MTIME, TSTEP( 2 ) )
!        CALL NEXTIME( MDATE, MTIME, TSTEP( 3 ) )
 
         IF ( MDATE .LT. LDATE ) THEN
            REVERT = .FALSE.
            ELSE IF ( MDATE .EQ. LDATE ) THEN
               IF ( MTIME .LE. LTIME ) THEN
                  REVERT = .FALSE.
                  ELSE
                  REVERT = .TRUE.
                  END IF
            ELSE   ! MDATE .GT. LDATE
            REVERT = .TRUE.
            END IF
 
         IF ( REVERT ) THEN
            XMSG = 'Current scenario interpolation step not available in '
     &           // MET_CRO_3D(1:TRIMLEN( MET_CRO_3D ) )
            CALL M3MESG( XMSG )
!           CALL NEXTIME( MDATE, MTIME, -SEC2TIME( STEP / 2 ) )
            WRITE( AMSG,'( 2I8 )' ) LDATE, LTIME
            XMSG = 'Using data for last file step: ' // AMSG
            CALL M3MESG( XMSG )
            MDATE = LDATE
            MTIME = LTIME
            END IF

         IF ( .NOT. INTERPX( MET_CRO_3D, 'DENSA_J', PNAME,
!    &              1,NCOLS, 1,NROWS, 1,NLAYS,
     &              STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &              MDATE, MTIME, RHOJ ) ) THEN
            XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

!        CALL DJINTERP ( MDATE, MTIME, TSTEP( 2 ), RHOJ ) ! for hydrostatic

         END IF

C Compute adjustment factor

      ASPC = GC_STRT - 1 + N_GC_SPCD
 
      IF ( IADJ .EQ. 0 ) THEN

         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  ADJFAC( COL,ROW,LVL ) = 1.0
                  END DO
               END DO
            END DO

         ELSE IF ( IADJ .EQ. 1 ) THEN

         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  ADJFAC( COL,ROW,LVL ) = 1.0 / CGRID( COL,ROW,LVL,ASPC )
                  END DO
               END DO
            END DO

         ELSE IF ( IADJ .EQ. 2 ) THEN

C Interpolate DENS at step
 
         MDATE  = JDATE
         MTIME  = JTIME
         STEP   = TIME2SEC( TSTEP( 2 ) )
         CALL NEXTIME( MDATE, MTIME, TSTEP( 2 ) )
!        STEP   = TIME2SEC( TSTEP( 3 ) )
!        CALL NEXTIME( MDATE, MTIME, TSTEP( 3 ) )
 
         IF ( .NOT. INTERPX( MET_CRO_3D, 'DENS', PNAME,
!    &              1,NCOLS, 1,NROWS, 1,NLAYS,
     &              STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &              MDATE, MTIME, DENS ) ) THEN
            XMSG = 'Could not read DENS from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  ADJFAC( COL,ROW,LVL ) = DENS( COL,ROW,LVL )
     &                                  / CGRID( COL,ROW,LVL,ASPC )
                  END DO
               END DO
            END DO

         ELSE IF ( IADJ .EQ. 3 ) THEN

         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  DRATE  = 1.0 - CGRID( COL,ROW,LVL,ASPC )
     &                   / RHOJ1( COL,ROW,LVL ) 
                  ADJFAC( COL,ROW,LVL ) = EXP( DRATE )
                  END DO
               END DO
            END DO

         ELSE IF ( IADJ .EQ. 4 ) THEN

         FRSTWRT = .TRUE.
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  DRATE = 1.0 - CGRID( COL,ROW,LVL,ASPC )
     &                  / RHOJ( COL,ROW,LVL )
     &                  + 1.0 - RHOJ1( COL,ROW,LVL )
     &                  / RHOJ( COL,ROW,LVL )
                  ADJFAC( COL,ROW,LVL ) = EXP( DRATE )
                  IF ( ADJFAC( COL,ROW,LVL ) .GT. TESTVAL ) THEN
                     IF ( FRSTWRT ) THEN
                        FRSTWRT = .FALSE.
                        WRITE( LOGDEV, 1001) TESTVAL
1001                    FORMAT( /5X, 'AdjFactor exceeds ', F10.5,
     &                           2X, 'Col,Row,Lvl,AdjFactor:' )
                        END IF
                     WRITE( LOGDEV,1003 ) COL, ROW, LVL, ADJFAC( COL,ROW,LVL )
1003                 FORMAT( 5X, 3I6.3, 2X, 1PE15.7)
                     END IF
                  END DO
               END DO
            END DO

         ELSE IF ( IADJ .EQ. 5 ) THEN

         FRSTWRT = .TRUE.
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  DRATE = 1.0 - CGRID( COL,ROW,LVL,ASPC )
     &                  / RHOJ( COL,ROW,LVL ) 
                  ADJFAC( COL,ROW,LVL ) = EXP( DRATE )
                  IF ( ADJFAC( COL,ROW,LVL ) .GT. TESTVAL ) THEN
                     IF ( FRSTWRT ) THEN
                        FRSTWRT = .FALSE.
                        WRITE( LOGDEV, 1001 ) TESTVAL
                        END IF
                     WRITE( LOGDEV,1003 ) COL, ROW, LVL, ADJFAC( COL,ROW,LVL )
                     END IF
                  END DO
               END DO
            END DO

         ELSE

         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  ADJFAC( COL,ROW,LVL ) = RHOJ( COL,ROW,LVL )
     &                                  / CGRID( COL,ROW,LVL,ASPC )
                  END DO
               END DO
            END DO

         END IF  ! adjcon choice

C Adjust concentrations

      DO SPC = 1, N_GC_ADV
         ASPC = GC_ADV_MAP( SPC ) + GC_STRT - 1
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = CGRID( COL,ROW,LVL,ASPC )
     &                                      * ADJFAC( COL,ROW,LVL )
                  END DO
               END DO
            END DO
         END DO
 
      DO SPC = 1, N_AE_ADV
         ASPC = AE_ADV_MAP( SPC ) + AE_STRT - 1
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = CGRID( COL,ROW,LVL,ASPC )
     &                                      * ADJFAC( COL,ROW,LVL )
                  END DO
               END DO
            END DO
         END DO
 
      DO SPC = 1, N_NR_ADV
         ASPC = NR_ADV_MAP( SPC ) + NR_STRT - 1
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = CGRID( COL,ROW,LVL,ASPC )
     &                                      * ADJFAC( COL,ROW,LVL )
                  END DO
               END DO
            END DO
         END DO
 
!     DO SPC = 1, N_TR_ADV
C** testing, exclude adjustment for the last advected tracer (TRAC01ADV)
      DO SPC = 1, N_TR_ADV-1
         ASPC = TR_ADV_MAP( SPC ) + TR_STRT - 1
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = CGRID( COL,ROW,LVL,ASPC )
     &                                      * ADJFAC( COL,ROW,LVL )
                  END DO
               END DO
            END DO
         END DO

C Load CGRID(*,*,*,N_GC_SPCD) with RHOJ
 
      ASPC = GC_STRT - 1 + N_GC_SPCD

      IF ( IADJ .EQ. 0 ) THEN             ! do nothing

         ELSE IF ( IADJ .EQ. 1 ) THEN     ! calibration gas -- 1.0 
 
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = 1.0
                  END DO
               END DO
            END DO

         ELSE IF ( IADJ .EQ. 2 ) THEN     ! SAQM method #1
 
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = DENS( COL,ROW,LVL )
                  END DO
               END DO
            END DO

         ELSE                             ! all others
 
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  CGRID( COL,ROW,LVL,ASPC ) = RHOJ( COL,ROW,LVL )
                  END DO
               END DO
            END DO

         END IF  ! adjcon choice
  
      RETURN

      END
