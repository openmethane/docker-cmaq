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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/eddy_ddm3d/s_rdemis.F,v 1.1.1.1 2010/06/14 16:03:10 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
C
C S_RDEMIS_GC:
C   This subroutine, based on rdemis.F, reads in emissions files
C that have been split into area, biogenic, point, non-road and mobile.
C 
C   Written by Daniel Cohan, 2002, based on rdemis.F
C   Jul 2006: updated for following changes in rdemis.F:
C        31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module

C 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE S_RDEMIS_GC ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, EMIS,
     &                        EMISFILE )

C read & interpolate gas chem emissions (NSPC_EMIS is zero if EMISVD = .FALSE.)

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
 
!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
! ythu removed 11202006
!     INCLUDE SUBST_GC_EMIS_B ! biogenic
!     INCLUDE SUBST_GC_EMIS_M ! mobile
!     INCLUDE SUBST_AE_EMIS   ! aerosol emis surrogate names and map table
      INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
!     INCLUDE SUBST_NR_EMIS_B ! biogenic
!     INCLUDE SUBST_NR_EMIS_M ! mobile
      INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names and map table
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID  ! coord. and domain definitions (req IOPARMS)

C arguments:

      INTEGER      MDATE, MTIME              ! internal simulation date&time
      INTEGER      NSPC_EMIS                 ! no. of emis species in gas chem
      INTEGER      EMISLYRS                  ! no. of file emissions layers
      LOGICAL      EM_TRAC                   ! are there tracer emissions?
      REAL      :: EMIS( :,:,:,: )           ! emis file input

      INTEGER, EXTERNAL :: TRIMLEN

C Local variables:

!     REAL, ALLOCATABLE :: BUFF( :,:,: ) ! read buffer
      REAL BUFF( NCOLS, NROWS, EMISLYRS )
      INTEGER      ALLOCSTAT, ELAYS
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER( 16 ), SAVE :: PNAME = 'S_RDEMIS'
      CHARACTER( 120 ) :: XMSG = ' '
      INTEGER      C, R, L, V                ! loop induction variables

      INTEGER      GXOFF, GYOFF              ! global origin offset from file
      INTEGER, SAVE :: STARTCOL, ENDCOL, STARTROW, ENDROW
      CHARACTER( 16 ) :: EMISFILE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
C Get domain decomp info from the emissions file
         CALL SUBHFILE ( EMISFILE, GXOFF, GYOFF,
     &                   STARTCOL, ENDCOL, STARTROW, ENDROW )
      END IF   ! FIRSTIME

!     ALLOCATE ( BUFF( MY_NCOLS,MY_NROWS,EMISLYRS ), STAT = ALLOCSTAT )
!     IF ( ALLOCSTAT .NE. 0 ) THEN
!         XMSG = 'BUFF memory allocation failed'
!         CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
!     END IF

         
C     DO V = 1, NSPC_EMIS

         ! Area, pt, and non-rd emissions have same species as overall
         IF ( ( EMISFILE .EQ. EMIS_A ) .OR. 
     &        ( EMISFILE .EQ. EMIS_B ) .OR.
     &        ( EMISFILE .EQ. EMIS_M ) .OR.
     &        ( EMISFILE .EQ. EMIS_P ) .OR.
     &        ( EMISFILE .EQ. EMIS_N ) ) THEN

         IF ( .NOT. DESC3( EMISFILE ) ) THEN
            XMSG = 'Could not get '// EMISFILE // ' file description'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
            END IF

         ELAYS = MIN(NLAYS3D, EMISLYRS)  

      DO V = 1, NSPC_EMIS
           IF ( .NOT. INTERPX( EMISFILE, GC_EMIS( V ), PNAME,
     &                       STARTCOL,ENDCOL, STARTROW,ENDROW, 1,ELAYS,
!    &                       MDATE, MTIME, EMIS( 1,1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
             BUFF = 0.0
             XMSG = 'Could not read  '
     &           // GC_EMIS( V )( 1:TRIMLEN( GC_EMIS( V ) ) )
     &           // ' from ' // EMISFILE
!            CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
           END IF

           DO R = 1, MY_NROWS
             DO C = 1, MY_NCOLS
               DO L = 1, ELAYS   
                  EMIS( V,L,C,R ) = GC_EMIS_FAC( V ) * BUFF( C,R,L ) 
               END DO
             END DO
           END DO

      ENDDO
C ythu removed 11202006
C         ! Biogenic Emissions
C         ELSE IF ( EMISFILE .EQ. EMIS_B ) THEN
C           IF ( .NOT. INTERPX( EMISFILE, GC_EMIS_B( V ), PNAME,
C     &                       STARTCOL,ENDCOL, STARTROW,ENDROW, 1,1,
C     &                       MDATE, MTIME, BUFF ) ) THEN
C             BUFF = 0.0
C             XMSG = 'Could not read  '
C     &           // GC_EMIS_B( V )( 1:TRIMLEN( GC_EMIS_B( V ) ) )
C     &           // ' from ' // EMISFILE
C             CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
C           END IF
C
C           DO R = 1, MY_NROWS
C             DO C = 1, MY_NCOLS
C               EMIS( V,1,C,R ) = GC_EMIS_FAC_B( V ) * BUFF( C,R,1 )
C             END DO
C           END DO
C
C         ! Mobile Emissions
C         ELSE IF ( EMISFILE .EQ. EMIS_M ) THEN
C           IF ( .NOT. INTERPX( EMISFILE, GC_EMIS_M( V ), PNAME,
C     &                       STARTCOL,ENDCOL, STARTROW,ENDROW, 1,1,
C     &                       MDATE, MTIME, BUFF ) ) THEN
C             BUFF = 0.0
C             XMSG = 'Could not read  '
C     &           // GC_EMIS_M( V )( 1:TRIMLEN( GC_EMIS_M( V ) ) )
C     &           // ' from ' // EMISFILE
C             CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
C           END IF
C
C          DO R = 1, MY_NROWS
C            DO C = 1, MY_NCOLS
C               EMIS( V,1,C,R ) = GC_EMIS_FAC_M( V ) * BUFF( C,R,1 )
C            END DO
C          END DO
C
        ELSE
          XMSG = 'Emissions file type not recognized.'
          CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
        END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        WRITE( *,* ) '  GC EMIS: ', V, '  ', GC_EMIS( N )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     END DO

      RETURN

      ENTRY S_RDEMIS_NR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, EMIS,
     &                    EMISFILE )

C     DO V = 1, NSPC_EMIS

         ! Area, pt, and non-rd emissions have same species as overall
         IF ( ( EMISFILE .EQ. EMIS_A ) .OR.
     &        ( EMISFILE .EQ. EMIS_B ) .OR.
     &        ( EMISFILE .EQ. EMIS_M ) .OR.
     &        ( EMISFILE .EQ. EMIS_P ) .OR. 
     &        ( EMISFILE .EQ. EMIS_N ) ) THEN

         IF ( .NOT. DESC3( EMISFILE ) ) THEN
            XMSG = 'Could not get '// EMISFILE // ' file description'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
            END IF

         ELAYS = MIN(NLAYS3D, EMISLYRS)

      DO V = 1, NSPC_EMIS 
           IF ( .NOT. INTERPX( EMISFILE, NR_EMIS( V ), PNAME,
     &                       STARTCOL,ENDCOL, STARTROW,ENDROW, 1, ELAYS,
!    &                       MDATE, MTIME, EMIS( 1,1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
             BUFF = 0.0
             XMSG = 'Could not read  '
     &           // NR_EMIS( V )( 1:TRIMLEN( NR_EMIS( V ) ) )
     &           // ' from ' // EMISFILE
C            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
!            CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
           END IF

           DO R = 1, MY_NROWS
             DO C = 1, MY_NCOLS
               DO L = 1, ELAYS 
                  EMIS( V,L,C,R ) = NR_EMIS_FAC( V ) * BUFF( C,R,L )
               END DO
             END DO
           END DO
      ENDDO
C ythu removed 11202006
C         ! Biogenic Emissions
C         ELSE IF ( EMISFILE .EQ. EMIS_B ) THEN
CC There are not currently biogenic NR emissions
C            CONTINUE
C
C         ! Mobile emissions
C         ELSE IF ( EMISFILE .EQ. EMIS_M ) THEN
CC There are not currently mobile NR emissions
C            CONTINUE
C
        ELSE
          XMSG = 'Emissions file type not recognized.'
          CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
        END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        WRITE( *,* ) '  NR EMIS: ', V, '  ', NR_EMIS( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     END DO

      RETURN

      END

#else

      SUBROUTINE S_RDEMIS_GC ( MDATE )
 
      INTEGER MDATE

      RETURN
      END
#endif
