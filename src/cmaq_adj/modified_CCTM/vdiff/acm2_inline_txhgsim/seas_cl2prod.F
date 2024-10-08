
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhgsim/seas_cl2prod.F,v 1.1.1.1 2010/06/14 16:03:09 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SEAS_CL2PROD ( MDATE, MTIME, NSPC_EMIS, EMIS )

C added constant Cl2 source to EMIS array based on ocean and surf zone fraction

C-----------------------------------------------------------------------
C  24/03/08 - WTH: start initial version 
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
 
!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC    ! gas species tabel
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_CONST     ! constants

!     INCLUDE SUBST_COORD_ID  ! coord. and domain definitions (req IOPARMS)

C arguments:

      INTEGER      MDATE, MTIME              ! internal simulation date&time
      INTEGER      NSPC_EMIS                 ! no. of emis species in gas chem
      LOGICAL      EM_TRAC                   ! are there tracer emissions?

      REAL      :: EMIS( :,:,: )             ! emis file input

      INTEGER, EXTERNAL :: TRIMLEN
      INTEGER, EXTERNAL :: INDEX1

C Local variables:

      REAL, ALLOCATABLE, SAVE :: OCEAN( :,: )     ! fractional seawater cover
      REAL, ALLOCATABLE, SAVE :: SZONE( :,: )     ! fractional surf-zone cover

      REAL, ALLOCATABLE, SAVE :: FRACTION_LAYER( :,: )  ! fraction of layer with Cl2 production


      REAL, ALLOCATABLE, SAVE :: DENS( :, :, : )  ! Air density [ kg/m**-3 ]

      INTEGER      ALLOCSTAT
      INTEGER, SAVE :: LOGDEV

      INTEGER, SAVE :: INDEX_CL2   ! array index of Cl2 in EMIS array

      CHARACTER( 16 ), SAVE :: PNAME = 'SEAS_CL2PROD'
      CHARACTER( 120 ) :: XMSG = ' '
      INTEGER      C, R, L, N, V             ! loop induction and scratch variables

      INTEGER      GXOFF, GYOFF              ! global origin offset from file

      INTEGER, SAVE :: STARTCOL_O, ENDCOL_O, STARTROW_O, ENDROW_O
      INTEGER, SAVE :: STRTCOL_M, ENDCOL_M, STRTROW_M, ENDROW_M


      CHARACTER( 16 ) :: VNAME            ! variable name

! Observed Cl2 production on based 100 pptV per day reported in Spicer et al (1998)

      REAL, PARAMETER :: OBSER_CL2PROD  = 4.8115E-12 ! (moles/(m^2*s))
      REAL, PARAMETER :: ZONE_CL2PROD   = 100.0      ! height of Cl2 production zone, m
      REAL, PARAMETER :: DENS_293K      = 1.204      ! air mass density at 293.15 K(kg/m^3)

      REAL, SAVE      :: FACT_CL2PROD     !  conversion factor for Cl2 production [(moles/s)*(m^3/Kg)]

      LOGICAL, SAVE   :: FIRSTIME = .TRUE.

C-----------------------------------------------------------------------

      XMSG = 'Simulation uses CL2 emissions over oceans for '
     &       // 'mercury photochemistry.'
      CALL M3MESG( XMSG )


      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()


         ALLOCATE ( OCEAN( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'memory allocation for OCEAN failed'
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         ENDIF

         ALLOCATE ( SZONE( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'memory allocation for SZONE failed'
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         ENDIF

         ALLOCATE ( DENS( MY_NCOLS,MY_NROWS,NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'memory allocation for DENS failed'
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         ENDIF

c *** Open the ocean file, which contains the ocean and surf-zone fractions

         IF ( .NOT. OPEN3( OCEAN_1, FSREAD3, PNAME ) ) THEN
               XMSG = 'Could not open '// OCEAN_1 // ' file'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         ENDIF
               
         IF ( .NOT. DESC3( OCEAN_1 ) ) THEN
               XMSG = 'Could not get '// OCEAN_1 // ' file description'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
         END IF

C Get domain decomp info from the OCEANs file

         CALL SUBHFILE ( OCEAN_1, GXOFF, GYOFF,
     &                   STARTCOL_O, ENDCOL_O, STARTROW_O, ENDROW_O )


C Open the met file for mass density

         IF ( .NOT. OPEN3( MET_CRO_3D, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open  MET_CRO_3D  file '
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOL_M, ENDCOL_M, STRTROW_M, ENDROW_M )


         VNAME = 'CL2'
         N = INDEX1( VNAME, N_GC_SPC, GC_SPC )
         IF ( N .LT. 1 ) THEN
            XMSG = 
     &        'Could not find ' // VNAME // 'in gas species table'
        
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT3 )
         ENDIF

         INDEX_CL2 = -1
         DO L = 1, N_GC_EMIS
            IF( N .EQ. GC_EMIS_MAP( L ) )THEN
                INDEX_CL2 = L
                EXIT
             ENDIF
         ENDDO

         IF( INDEX_CL2 .LT. 1 )THEN
            XMSG = 
     &        'Could not find index of ' // VNAME // 
     &        'in gas emissions map'
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT3 )
         ELSE
            XMSG =  'Adding to emissions species ' 
     &           //  GC_EMIS( INDEX_CL2 ) 
     &           // 'for model species ' // GC_SPC( N )
         ENDIF
       
         FACT_CL2PROD  = ( XCELL_GD * YCELL_GD)
     &                 *  OBSER_CL2PROD/DENS_293K

      ENDIF   ! FIRSTIME

C *** air density (kg/m3)

      VNAME = 'DENS'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOL_M, ENDCOL_M, STRTROW_M, ENDROW_M, 
     &                    1, NLAYS, MDATE, MTIME, DENS ) ) THEN
        XMSG = 'Could not interpolate '// VNAME // ' from MET_CRO_3D '
        CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      ENDIF

      VNAME = 'OPEN'
      IF ( .NOT. INTERPX( OCEAN_1, VNAME, PNAME,
     &                    STARTCOL_O, ENDCOL_O, STARTROW_O, ENDROW_O,
     &                    1, 1, MDATE, MTIME, OCEAN ) ) THEN
          XMSG = 'Could not interpolate OPEN from ' // OCEAN_1
          CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      ENDIF
         

      VNAME = 'SURF'
      IF ( .NOT. INTERPX( OCEAN_1, VNAME, PNAME,
     &                    STARTCOL_O, ENDCOL_O, STARTROW_O, ENDROW_O,
     &                    1, 1, MDATE, MTIME, SZONE ) ) THEN
          XMSG = 'Could not interpolate OPEN from ' // OCEAN_1
          CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      ENDIF

C "Stores read an entire cache line, modify the target, then write back the
C  entire line. Thus, non-consecutive stores are worse than non-consecutive
C  loads."

!W              r = 24
!W              c =  4
!W              WRITE(LOGDEV,100)PNAME,':At ',r,c,GC_EMIS(INDEX_CL2 ),
!W     &               EMIS(  INDEX_CL2, C, R )

!W              r = 26
!W              c = 27
!W              WRITE(LOGDEV,100)PNAME,':At ',r,c,GC_EMIS(INDEX_CL2 ),
!W     &               EMIS(  INDEX_CL2, C, R )


         EMIS = 0.0

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS


            EMIS(  INDEX_CL2, C, R ) = EMIS(  INDEX_CL2, C, R )
     &                               + FACT_CL2PROD
     &                               * DENS(C, R, 1)*OCEAN(C, R)


            ENDDO
         ENDDO

         
!W              r = 24
!W              c =  4
!W              WRITE(LOGDEV,100)PNAME,': At ',r,c,GC_EMIS(INDEX_CL2 ),
!W     &               EMIS(  INDEX_CL2, C, R ),OCEAN(C, R)
!W
!W              r = 26
!W              c = 27
!W              WRITE(LOGDEV,100)PNAME,': At ',r,c,GC_EMIS(INDEX_CL2 ),
!W     &               EMIS(  INDEX_CL2, C, R ),OCEAN(C, R)


      XMSG = 'CL2 emissions over oceans estimated for '
     &       // 'mercury photochemistry.'
      CALL M3MESG( XMSG )

100   FORMAT(A16,A5,2(I3,1X),A16,1X,3(ES12.4,1X))
      RETURN
      END
