
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/init/yamo_ddm3d/load_sengrid.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LOAD_SENGRID ( FNAME, JDATE, JTIME, SPC_CAT, SENGRID, 
     &                          NPMAX, SEN_PAR )

C-----------------------------------------------------------------------
C Function:
C   Initialize the model SENGRID array from file data

C Revision history:
C   Created by Daniel Cohan, Nov. 2002, based on load_cgrid
C   Modified by D. Cohan, Jan. 2004, for compatibility with Intel compiler.
C   July 2006: Modified by D. Cohan to reflect the following in load_cgrid:
C       20 Nov 03 J.Young: enable loading RHOJ
C       31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_CONC     ! gas chem conc file species and map table
      INCLUDE SUBST_AE_CONC     ! aerosol conc file species and map table
      INCLUDE SUBST_GC_ICBC     ! gas chem ic/bc surrogate names and map table
      INCLUDE SUBST_AE_ICBC     ! aerosol ic/bc surrogate names and map table
      INCLUDE SUBST_NR_ICBC     ! non-react ic/bc surrogate names and map table
      INCLUDE SUBST_TR_ICBC     ! tracer ic/bc surrogate names and map table
      INCLUDE SUBST_GC_G2AE     ! gas chem aerosol species and map table
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      CHARACTER( 16 ) :: FNAME
      INTEGER      JDATE
      INTEGER      JTIME
      CHARACTER( 2 ) ::  SPC_CAT
      REAL, POINTER :: SENGRID( :,:,:,:,: )  ! for initial SENS
      INTEGER NPMAX
      CHARACTER( 8 ), POINTER :: SEN_PAR( : )  ! names for sens parameters

      INTEGER SENNUM
      CHARACTER( 7 ) :: SPCNAME
      CHARACTER( 16) :: SENNAME
      INTEGER P
 
C Parameters:

      INTEGER, PARAMETER :: MXSPC = N_GC_SPCD + N_AE_SPC
     &                            + N_NR_SPC + N_TR_SPC

C File variables:

      REAL         DENS( NCOLS,NROWS,NLAYS )       ! air density (kg/m^3)
      REAL         RHOJ( NCOLS,NROWS,NLAYS ) ! air density X Jacobian (kg/m^2)

C External Functions (not already declared by IODECL3.EXT):

      INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.
      INTEGER, EXTERNAL :: INDEX1       !  looks up name in table.
      INTEGER, EXTERNAL :: TRIMLEN      !  string length, excl. trailing blanks

C Local Variables


      INTEGER      GXOFF, GYOFF               ! global origin offset from file

C for INTERPX
!     INTEGER :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C for XTRACT3
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER          STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI
      REAL         DBUFF( NCOLS,NROWS,NLAYS )

      INTEGER      SPC_STRT
      INTEGER      N_SPCS
      INTEGER      INDX( MXSPC * NPMAX )      ! Variable and surrogate indices
      INTEGER      NDX                        ! loop copy of INDX
      INTEGER      C, R, L, SPC, V            ! loop counters
      INTEGER      ASPC                       ! CGRID RHOJ pointer


      CHARACTER( 16 ) :: PNAME = 'LOAD_SENGRID'
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 40 ) :: CHWARN = 'Domain extents difft from model for '

      CHARACTER( 16 ) :: VNAME

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )
C IOFDESC common now loaded with MET_CRO_3D header
         END IF

      IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
         XMSG = 'Could not open ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
 
      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // FNAME // ' file description'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
 
      IF ( GL_NCOLS .NE. NCOLS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NCOLS: ', GL_NCOLS, '  NCOLS3D: ', NCOLS3D
         END IF
 
      IF ( GL_NROWS .NE. NROWS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NROWS: ', GL_NROWS, '  NROWS3D: ', NROWS3D
         END IF
 
      IF ( NLAYS .NE. NLAYS3D ) THEN
         XMSG = 'Wrong number of layers in ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Get CGRID offsets

      CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C Get INDX
      IF ( SPC_CAT .EQ. 'GC' ) THEN
         SPC_STRT = GC_STRT
         N_SPCS = N_GC_SPC
         DO SPC = 1, N_SPCS
            SPCNAME = GC_SPC( SPC )( 1:7 )
            DO P = 1, NPMAX
               SENNUM = ( SPC - 1 ) * NPMAX + P
               SENNAME = SPCNAME(1:trimlen(SPCNAME)) // '_' //
     &                         SEN_PAR( P )
               NDX = INDEX1( SENNAME, NVARS3D, VNAME3D )
               
               INDX( SENNUM ) = NDX
               IF ( NDX .EQ. 0 ) THEN
                 XMSG = 'No IC found for sensitivity ' //
     &                SENNAME  //
     &                ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                ' Set to zero.'
                 CALL M3MESG ( XMSG )
               END IF
            END DO
         END DO
      ELSE IF ( SPC_CAT .EQ. 'AE' ) THEN
         SPC_STRT = AE_STRT
         N_SPCS = N_AE_SPC
         DO SPC = 1, N_SPCS
            SPCNAME = AE_SPC( SPC )( 1:7 )
            DO P = 1, NPMAX
              SENNUM = ( SPC + SPC_STRT - 2 ) * NPMAX + P
              SENNAME = SPCNAME(1:trimlen(SPCNAME)) // '_' //
     &                         SEN_PAR( P )
              NDX = INDEX1( SENNAME, NVARS3D, VNAME3D )
              INDX( SENNUM ) = NDX
              IF ( NDX .EQ. 0 ) THEN
                XMSG = 'No IC found for sensitivity ' //
     &                  SENNAME  //
     &                 ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                 ' Set to zero.'
                CALL M3MESG ( XMSG )
              END IF
            END DO
         END DO
      ELSE IF ( SPC_CAT .EQ. 'NR' ) THEN
         SPC_STRT = NR_STRT
         N_SPCS = N_NR_SPC
         DO SPC = 1, N_SPCS
           SPCNAME = NR_SPC( SPC )( 1:7 )
           DO P = 1, NPMAX
              SENNUM = ( SPC + SPC_STRT - 2 ) * NPMAX + P
              SENNAME = SPCNAME(1:trimlen(SPCNAME)) // '_' //
     &                         SEN_PAR( P )
              NDX = INDEX1( SENNAME, NVARS3D, VNAME3D )
              INDX( SENNUM ) = NDX
              IF ( NDX .EQ. 0 ) THEN
                XMSG = 'No IC found for sensitivity ' //
     &                  SENNAME  //
     &                 ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                 ' Set to zero.'
                CALL M3MESG ( XMSG )
              END IF
           END DO
         END DO
      ELSE IF ( SPC_CAT .EQ. 'TR' ) THEN
         SPC_STRT = TR_STRT
         N_SPCS = N_TR_SPC
         DO SPC = 1, N_SPCS
           DO P = 1, NPMAX
              SENNUM = ( SPC + SPC_STRT - 2 ) * NPMAX + P
              SENNAME = SPCNAME(1:trimlen(SPCNAME)) // '_' //
     &                         SEN_PAR( P )
              NDX = INDEX1( SENNAME, NVARS3D, VNAME3D )
              INDX( SENNUM ) = NDX
              IF ( NDX .EQ. 0 ) THEN
                XMSG = 'No IC found for sensitivity ' //
     &                  SENNAME  //
     &                 ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                 ' Set to zero.'
                CALL M3MESG ( XMSG )
              END IF
           END DO
         END DO
      ! DSC: I'm not clear if RJ would ever come up in sens and whether
      !      following lines are appropriate ???
      ELSE IF ( SPC_CAT .EQ. 'RJ' ) THEN
         N_SPCS = 0
!        WRITE( LOGDEV,'(/ 5X, A)' ) 'loading Density*Jacobian into SENGRID'
      ELSE
         XMSG = 'Species categories incorrect for SENGRID '
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
      END IF
        
C Read into SENGRID

      CALL SUBHFILE ( FNAME, GXOFF, GYOFF,
     &                STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI )
C IOFDESC common now loaded with FNAME header

      DO SPC = 1, N_SPCS
         V = SPC_STRT - 1 + SPC
         DO P = 1, NPMAX
           SENNUM = ( SPC + SPC_STRT - 2 ) * NPMAX + P
           NDX = INDX( SENNUM )

           IF ( NDX .GT. 0 ) THEN
             IF ( .NOT. XTRACT3( FNAME, VNAME3D ( NDX ),
     &           1,NLAYS, STRTROWINI,ENDROWINI, STRTCOLINI,ENDCOLINI,
     &           JDATE, JTIME, DBUFF ) ) THEN
               XMSG = 'Could not read ' //
     &                VNAME3D( NDX )( 1:TRIMLEN( VNAME3D( NDX ) ) ) //
     &                ' from ' // FNAME
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
             END IF

             DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     SENGRID( C,R,L,P,V ) = DBUFF( C,R,L )
                  END DO
               END DO
             END DO
   
           ELSE
             DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     SENGRID( C,R,L,P,V ) = 0.0
                  END DO
               END DO
             END DO

           END IF   ! INDX .GT. 0

         END DO

      END DO

      IF ( N_SPCS .NE. 0 ) WRITE( LOGDEV,'(/ 5X, A)' )
     &                            SPC_CAT // ' loaded into SENGRID'

      ! DSC: Not clear about applicability of RHOJ code to DDM ???
      IF ( SPC_CAT .EQ. 'RJ' ) THEN

C Load RHOJ for transport and mixing ratio advection adjustment

         VNAME = 'DENSA_J'
         IF ( .NOT. XTRACT3( MET_CRO_3D, VNAME,
     &              1,NLAYS, STRTROWMC3,ENDROWMC3, STRTCOLMC3,ENDCOLMC3,
     &              JDATE, JTIME, RHOJ ) ) THEN
             XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         ASPC = GC_STRT - 1 + N_GC_SPCD
         DO P = 1, NPMAX
           DO L = 1, NLAYS
             DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  SENGRID( C,R,L,P,ASPC ) = RHOJ( C,R,L )
               END DO
             END DO
           END DO
         END DO

         WRITE( LOGDEV,'(/ 5X,A)' ) 'Density*Jacobian loaded to SENGRID'

      END IF   ! End RHOJ Code

C Close the file

      IF ( .NOT. CLOSE3( FNAME ) ) THEN
         XMSG = 'Could not close ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      RETURN
      END
