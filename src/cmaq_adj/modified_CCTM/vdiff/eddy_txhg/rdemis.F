
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/eddy_txhg/rdemis.F,v 1.1.1.1 2010/06/14 16:03:11 sjr Exp $
 
C what( 1 ) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE RDEMIS_GC ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, EMIS )

C read & interpolate gas chem emissions (NSPC_EMIS is zero if EMISVD = .FALSE.)

C-----------------------------------------------------------------------
C 1/12/99 David Wong at LM:
C   -- changed NROWS to MY_NROWS and NCOLS to MY_NCOLS in the loop structure 
C   -- added HGRID.EXT in the include file list
C   -- removed NCOLS and NROWS in the argument list
C   30 Aug 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX;
C                      allocatable VDEMIS
C   20 Aug 03 J.Young: fix error loading buff into emis when no input emissions
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
 
!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
!     INCLUDE SUBST_AE_EMIS   ! aerosol emis surrogate names and map table
      INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
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
!     REAL         EMIS( NCOLS,NROWS,NLAYS,* )  ! emis file input
      REAL      :: EMIS( :,:,:,: )           ! emis file input

      INTEGER, EXTERNAL :: TRIMLEN

C Local variables:

      REAL, ALLOCATABLE, SAVE :: BUFF( :,:,: ) ! read buffer
      INTEGER      ALLOCSTAT
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER( 16 ), SAVE :: PNAME = 'RDEMIS'
      CHARACTER( 120 ) :: XMSG = ' '
      INTEGER      C, R, L, V                ! loop induction variables

      INTEGER      GXOFF, GYOFF              ! global origin offset from file
      INTEGER, SAVE :: STARTCOL, ENDCOL, STARTROW, ENDROW

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
         ALLOCATE ( BUFF( MY_NCOLS,MY_NROWS,EMISLYRS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'BUFF memory allocation failed'
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
C Get domain decomp info from the emissions file
         CALL SUBHFILE ( EMIS_1, GXOFF, GYOFF,
     &                   STARTCOL, ENDCOL, STARTROW, ENDROW )
         END IF   ! FIRSTIME
         
      DO V = 1, NSPC_EMIS

C                                    vvvvvvv - emis surrogate names
         IF ( .NOT. INTERPX( EMIS_1, GC_EMIS( V ), PNAME,
     &                       STARTCOL,ENDCOL, STARTROW,ENDROW, 1,EMISLYRS,
!    &                       MDATE, MTIME, EMIS( 1,1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
            BUFF = 0.0
            XMSG = 'Could not read  '
     &           // GC_EMIS( V )( 1:TRIMLEN( GC_EMIS( V ) ) )
     &           // ' from ' // EMIS_1
            CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
            END IF

C "Stores read an entire cache line, modify the target, then write back the
C  entire line. Thus, non-consecutive stores are worse than non-consecutive
C  loads."
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, EMISLYRS
                  EMIS( V,L,C,R ) = GC_EMIS_FAC( V ) * BUFF( C,R,L )
                  END DO
               END DO
            END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        WRITE( *,* ) '  GC EMIS: ', V, '  ', GC_EMIS( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO

      RETURN

      ENTRY RDEMIS_NR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, EMIS )

      DO V = 1, NSPC_EMIS

C                                    vvvvvvv - emis surrogate names
         IF ( .NOT. INTERPX( EMIS_1, NR_EMIS( V ), PNAME,
     &                       STARTCOL,ENDCOL, STARTROW,ENDROW, 1,EMISLYRS,
!    &                       MDATE, MTIME, EMIS( 1,1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
            XMSG = 'Could not read  '
     &           // NR_EMIS( V )( 1:TRIMLEN( NR_EMIS( V ) ) )
     &           // ' from ' // EMIS_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, EMISLYRS
                  EMIS( V,L,C,R ) = NR_EMIS_FAC( V ) * BUFF( C,R,L )
                  END DO
               END DO
            END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        WRITE( *,* ) '  NR EMIS: ', V, '  ', NR_EMIS( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO

      RETURN

      ENTRY RDEMIS_TR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, EMIS )

      DO V = 1, NSPC_EMIS

C                                            vvvvvvv - emis surrogate names
            IF ( .NOT. INTERPX( EMIS_TRAC_1, TR_EMIS( V ), PNAME,
     &                          STARTCOL,ENDCOL, STARTROW,ENDROW, 1,EMISLYRS,
!    &                          MDATE, MTIME, EMIS( 1,1,1,V ) ) ) THEN
     &                          MDATE, MTIME, BUFF ) ) THEN
            XMSG = 'Could not read  '
     &           // TR_EMIS( V )( 1:TRIMLEN( TR_EMIS( V ) ) )
     &           // ' from ' // EMIS_TRAC_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, EMISLYRS
                  EMIS( V,L,C,R ) = TR_EMIS_FAC( V ) * BUFF( C,R,L )
                  END DO
               END DO
            END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        WRITE( *,* ) '  TR EMIS: ', V, '  ', TR_EMIS( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO

      RETURN
      END
