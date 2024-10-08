#ifdef sens
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/eddy_ddm3d/s_opemis.F,v 1.1.1.1 2010/06/14 16:03:10 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE S_OPEMIS ( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM,
     &                        EMISLYRS, EMISFILE )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C S_OPEMIS:
C   This subroutine, based on opemis.F, is used to read in emission
C files that have been split into area, biogenic, mobile, non-road,
C and point emissions
C
C   Written by D. Cohan, 2002
C  7 Dec 04 - J.Young: vert dyn alloc - Use VGRD_DEFN (opemis.F)
C  Jul 06  - D.Cohan: adjust to reflect above changes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      USE VGRD_DEFN             ! vertical layer specifications

      IMPLICIT NONE

!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
! ythu removed 11202006
!     INCLUDE SUBST_GC_EMIS_B   ! gas chem emis surrogate names and map table
!     INCLUDE SUBST_GC_EMIS_M   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
!     INCLUDE SUBST_NR_EMIS_B   ! gas chem emis surrogate names and map table
!     INCLUDE SUBST_NR_EMIS_M   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names and map table
!     INCLUDE SUBST_EMLYRS_ID ! emissions layers parameter
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations

C Arguments:

      INTEGER      JDATE      ! current model date, coded YYYYDDD
      INTEGER      JTIME      ! current model time, coded HHMMSS
      INTEGER      NEMIS      ! no. of gas chem emissions species for vdiff
      LOGICAL      EM_TRAC    ! are there tracer emissions?
      REAL         CONVEM     ! conversion for emissions rates
      INTEGER      EMISLYRS   ! no. of emissions layers on file

C External Functions not previously declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: INDEX1, TRIMLEN

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'S_OPEMIS'
      CHARACTER( 96 ) :: XMSG
      CHARACTER( 16 ) :: UNITSCK

      CHARACTER( 16 ) :: EMISFILE

      LOGICAL ::   WRFLG = .FALSE.
      INTEGER      LOGDEV
      INTEGER      V, N, S, N2        ! induction variables

C-----------------------------------------------------------------------
 
      LOGDEV = INIT3 ()

C Following lines commented out ... 
C Sensitivity does not yet handle tracer emissions

C Open the tracer emissions file

!      IF ( N_TR_EMIS .GT. 0 ) THEN
!
!         IF ( .NOT. OPEN3( EMIS_TRAC_1, FSREAD3, PNAME ) ) THEN
!
!            XMSG = 'Could not open '// EMIS_TRAC_1 // ' file'
!            CALL M3MESG( XMSG )
!            EM_TRAC = .FALSE.
!
!            ELSE
!
!            IF ( .NOT. DESC3( EMIS_TRAC_1 ) ) THEN
!               XMSG = 'Could not get '// EMIS_TRAC_1 // ' file description'
!               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!               END IF
! 
!            V = INDEX1( TR_EMIS( 1 ), NVARS3D, VNAME3D )
!            IF ( V .NE. 0 ) THEN
!               UNITSCK = UNITS3D( V )
!               ELSE
!               XMSG = 'Emissions species '
!     &              // TR_EMIS( 1 )( 1:TRIMLEN( TR_EMIS( 1 ) ) )
!     &              // ' not found on ' // EMISFILE
!               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!               END IF
!
!            DO N = 2, N_TR_EMIS
!               V = INDEX1( TR_EMIS( N ), NVARS3D, VNAME3D )
!               IF ( V .NE. 0 ) THEN
!                  IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
!                     XMSG = 'Units not uniform on ' // EMIS_TRAC_1
!                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!                     END IF
!                  ELSE
!                  XMSG = 'Emissions species '
!     &                 // TR_EMIS( N )( 1:TRIMLEN( TR_EMIS( 1 ) ) )
!     &                 // ' not found on ' // EMIS_TRAC_1
!                  CALL M3MESG( XMSG )
!                  END IF
!               END DO
!
!            END IF
!
!         END IF   ! tracer emissions

!     IF ( NEMIS + N_AE_EMIS + N_NR_EMIS + N_TR_EMIS .GT. 0 ) THEN
      IF ( NEMIS + N_NR_EMIS .GT. 0 ) THEN

C Open the emissions file (for gas chem, aerosols and non-reactive species)

         IF ( .NOT. OPEN3( EMISFILE, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open '// EMISFILE // ' file'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         IF ( .NOT. DESC3( EMISFILE ) ) THEN
            XMSG = 'Could not get '// EMISFILE // ' file description'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

         ! Area, point source, and non-road emissions files have same
         ! emitted gas-phase species as the total file
         IF ( ( EMISFILE .EQ. EMIS_A ) .OR.
     &        ( EMISFILE .EQ. EMIS_B ) .OR.
     &        ( EMISFILE .EQ. EMIS_M ) .OR.
     &        ( EMISFILE .EQ. EMIS_P ) .OR.
     &        ( EMISFILE .EQ. EMIS_N ) ) THEN
           DO N = 1, N_GC_EMIS
             V = INDEX1( GC_EMIS( N ), NVARS3D, VNAME3D )
             IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D( V )
               S = N + 1
               GO TO 101
             ELSE
                 XMSG = 'Emissions species '
     &             // GC_EMIS( N )( 1:TRIMLEN( GC_EMIS( N ) ) )
     &             // ' not found on ' // EMISFILE
                 CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
c sln 12/15/09 Check NR_EMIS for the chance that there is only NR_species in the sector emissions file
                 IF ( N .EQ. N_GC_EMIS ) THEN
                    DO N2 = 1, N_NR_EMIS
                       V = INDEX1( NR_EMIS( N2 ), NVARS3D, VNAME3D )
                       IF ( V .NE. 0 ) THEN
                          UNITSCK = UNITS3D( V )
                          S = 1
                          GO TO 101
                       ELSE
                          XMSG = 'Emissions species '
     &                    // NR_EMIS( N2 )( 1:TRIMLEN( NR_EMIS( N2 ) ) )
     &                    // ' not found on ' // EMISFILE
                          CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
                       END IF
                    END DO
                 ENDIF
             END IF
           END DO
           XMSG = ' No emissions species ' // ' found on ' // EMISFILE
           CALL M3WARN( PNAME, JDATE, JTIME, XMSG )

101        CONTINUE

           DO N = S, N_GC_EMIS
              V = INDEX1( GC_EMIS( N ), NVARS3D, VNAME3D )
              IF ( V .NE. 0 ) THEN
                IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
                  XMSG = 'Units not uniform on ' // EMISFILE
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                END IF
              ELSE
                XMSG = 'Emissions species '
     &              // GC_EMIS( N )( 1:TRIMLEN( GC_EMIS( N ) ) )
     &              // ' not found on ' // EMISFILE
                CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
              END IF
           END DO

           DO N = 1, N_NR_EMIS
             V = INDEX1( NR_EMIS( N ), NVARS3D, VNAME3D )
             IF ( V .NE. 0 ) THEN
               IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
                  XMSG = 'Units not uniform on ' // EMISFILE
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
             ELSE
               XMSG = 'Emissions species '
     &              // NR_EMIS( N )( 1:TRIMLEN( NR_EMIS( N ) ) )
     &              // ' not found on ' // EMISFILE
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
             END IF
           END DO

C ythu removed 11202006
C         ! Biogenic Emissions
C         ELSE IF ( EMISFILE .EQ. EMIS_B ) THEN
C           DO N = 1, N_GC_EMIS_B
C             V = INDEX1( GC_EMIS_B( N ), NVARS3D, VNAME3D )
C             IF ( V .NE. 0 ) THEN
C               UNITSCK = UNITS3D( V )
C               S = N + 1
C               GO TO 102
C             ELSE
C               XMSG = 'Emissions species '
C     &           // GC_EMIS_B( N )( 1:TRIMLEN( GC_EMIS_B( N ) ) )
C     &           // ' not found on ' // EMISFILE
C               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C             END IF
C           END DO
C           XMSG = ' No emissions species ' // ' found on ' // EMISFILE
C           CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C
C102        CONTINUE
C
C           DO N = S, N_GC_EMIS_B
C              V = INDEX1( GC_EMIS_B( N ), NVARS3D, VNAME3D )
C              IF ( V .NE. 0 ) THEN
C                IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
C                  XMSG = 'Units not uniform on ' // EMISFILE
C                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
C                END IF
C              ELSE
C                XMSG = 'Emissions species '
C     &              // GC_EMIS_B( N )( 1:TRIMLEN( GC_EMIS_B( N ) ) )
C     &              // ' not found on ' // EMISFILE
C                CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C              END IF
C           END DO
C
CC There are not currently biogenic NR emissions
C!          DO N = 1, N_NR_EMIS_B
C!            V = INDEX1( NR_EMIS_B( N ), NVARS3D, VNAME3D )
C!             IF ( V .NE. 0 ) THEN
C!               IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
C!                  XMSG = 'Units not uniform on ' // EMISFILE
C!                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
C!               END IF
C!             ELSE
C!               XMSG = 'Emissions species '
C!    &              // NR_EMIS_B( N )( 1:TRIMLEN( NR_EMIS_B( N ) ) )
C!    &              // ' not found on ' // EMISFILE
C!               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C!             END IF
C!           END DO
C
C         ! Mobile emissions
C         ELSE IF ( EMISFILE .EQ. EMIS_M ) THEN
C           DO N = 1, N_GC_EMIS_M
C             V = INDEX1( GC_EMIS_M( N ), NVARS3D, VNAME3D )
C             IF ( V .NE. 0 ) THEN
C               UNITSCK = UNITS3D( V )
C               S = N + 1
C               GO TO 103
C             ELSE
C               XMSG = 'Emissions species '
C     &           // GC_EMIS_M( N )( 1:TRIMLEN( GC_EMIS_M( N ) ) )
C     &           // ' not found on ' // EMISFILE
C               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C             END IF
C           END DO
C           XMSG = ' No emissions species ' // ' found on ' // EMISFILE
C           CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C
C103        CONTINUE
C
C           DO N = S, N_GC_EMIS_M
C              V = INDEX1( GC_EMIS_M( N ), NVARS3D, VNAME3D )
C              IF ( V .NE. 0 ) THEN
C                IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
C                  XMSG = 'Units not uniform on ' // EMISFILE
C                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
C                END IF
C              ELSE
C                XMSG = 'Emissions species '
C     &              // GC_EMIS_M( N )( 1:TRIMLEN( GC_EMIS_M( N ) ) )
C     &              // ' not found on ' // EMISFILE
C                CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C              END IF
C           END DO
C
CC There are not currently mobile NR emissions
C!          DO N = 1, N_NR_EMIS_M
C!            V = INDEX1( NR_EMIS_M( N ), NVARS3D, VNAME3D )
C!            IF ( V .NE. 0 ) THEN
C!              IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
C!                 XMSG = 'Units not uniform on ' // EMISFILE
C!                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
C!              END IF
C!            ELSE
C!              XMSG = 'Emissions species '
C!    &              // NR_EMIS_M( N )( 1:TRIMLEN( NR_EMIS_M( N ) ) )
C!    &              // ' not found on ' // EMISFILE
C!              CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
C!            END IF
C!          END DO
C
         ELSE
           XMSG = 'Emissions file not recognized.'
           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

!        DO N = 1, N_AE_EMIS
!           V = INDEX1( AE_EMIS( N ), NVARS3D, VNAME3D )
!           IF ( V .NE. 0 ) THEN
!              IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
!                 XMSG = 'Units not uniform on ' // EMISFILE
!                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!                 END IF
!              ELSE
!              XMSG = 'Emissions species '
!    &              // AE_EMIS( N )( 1:TRIMLEN( AE_EMIS( N ) ) )
!    &              // ' not found on ' // EMISFILE
!              CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
!              END IF
!           END DO


C Assume units uniform across gas and non-reactive species classes

         CALL UPCASE ( UNITSCK )

!        IF ( UNITSCK .EQ. 'G/S' .OR. UNITSCK .EQ. 'g/s' ) THEN
         IF ( UNITSCK .EQ. 'MOLES/S'   .OR.
     &        UNITSCK .EQ. 'MOLE/S'    .OR.
     &        UNITSCK .EQ. 'MOL/S'     .OR.
     &        UNITSCK .EQ. 'MOLES/SEC' .OR.
     &        UNITSCK .EQ. 'MOLE/SEC'  .OR.
     &        UNITSCK .EQ. 'MOL/SEC' ) THEN
            CONVEM = 1.0E-03  ! assuming gram-moles
!           ELSE IF ( UNITSCK .EQ. 'KG/HR' .OR. UNITSCK .EQ. 'kg/hr' ) THEN
!           CONVEM = 1.0 / 3600.0
            ELSE
            XMSG = 'Units incorrect on ' // EMISFILE
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

         END IF  ! (gas chem or non-react emissions)

!      EMISLYRS = NLAYS3D

!      IF ( NEMIS .GT. 0 ) THEN
!         WRFLG = .TRUE.
!         WRITE( LOGDEV,1001 )
!1001     FORMAT( / 10X, 'Gas Chemistry Emissions Processing in',
!     &              1X, 'Vertical diffusion ...' )
!         END IF

!     IF ( N_AE_EMIS .GT. 0 ) THEN
!        WRFLG = .TRUE.
!        WRITE( LOGDEV,1003 )
!1003     FORMAT( / 10X, 'Aerosol Emissions Processing in',
!     &              1X, 'Vertical diffusion ...' )
!        END IF

!      IF ( N_NR_EMIS .GT. 0 ) THEN
!         WRFLG = .TRUE.
!         WRITE( LOGDEV,1005 )
!1005     FORMAT( / 10X, 'Non-reactives Emissions Processing in',
!     &              1X, 'Vertical diffusion ...' )
!         END IF

!      IF ( N_TR_EMIS .GT. 0 .AND. EM_TRAC ) THEN
!         WRFLG = .TRUE.
!         WRITE( LOGDEV,1007 )
!1007     FORMAT( / 10X, 'Tracer Emissions Processing in',
!     &              1X, 'Vertical diffusion ...' )
!         END IF

!      IF ( WRFLG ) THEN
!         WRITE( LOGDEV,1009 ) EMISLYRS, NLAYS
!1009     FORMAT( / 10X, 'Number of Emissions Layers:         ', I3
!     &           / 10X, 'out of total Number of Model Layers:', I3 )
!         END IF

      RETURN
      END

#else
      SUBROUTINE S_OPEMIS ( JDATE )

      INTEGER JDATE 

       RETURN
      END


#endif
