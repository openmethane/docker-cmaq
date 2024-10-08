
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline/openlayout.F,v 1.1.1.1 2010/06/14 16:03:07 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPENLAYOUT( SDATE, STIME, TSTEP, EMLAYS, VGLVSXG,
     &                       NPTGRPS, MY_NSRC, LNAME )

C-----------------------------------------------------------------------
 
C  Description:
C    Opens the diagnostic layer fractions file(s).
 
C  Preconditions:
 
C  Subroutines and Functions Called:
C       I/O API 

C-----------------------------------------------------------------------
C Modified from:

C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C             System
C File: @(#)$Id: openlayout.F,v 1.1.1.1 2010/06/14 16:03:07 sjr Exp $
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C See file COPYRIGHT for conditions of use.
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C smoke@emc.mcnc.org
C Pathname: $Source: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline/openlayout.F,v $
C Last updated: $Date: 2010/06/14 16:03:07 $ 
 
C-----------------------------------------------------------------------
 
      USE HGRD_DEFN, ONLY: MYPE

      IMPLICIT NONE
 
C Includes:
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
!include      SUBST_IODECL    # I/O definitions and declarations
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC   ! file header data structure

C External Functions:
      CHARACTER( 16 ), EXTERNAL :: PROMPTMFILE
      INTEGER,         EXTERNAL :: SETUP_LOGDEV, TRIMLEN

C Local Parameters:

C Arguments:
      INTEGER,        INTENT( IN )  :: SDATE    ! Julian start date
      INTEGER,        INTENT( IN )  :: STIME    ! start time (HHMMSS)
      INTEGER,        INTENT( IN )  :: TSTEP    ! time step (HHMMSS)
      INTEGER,        INTENT( IN )  :: EMLAYS   ! number of emissions layers
      REAL,           INTENT( IN )  :: VGLVSXG( 0:MXLAYS3 ) !  vertical coord values
      INTEGER,        INTENT( IN )  :: NPTGRPS  ! number of output files
      INTEGER,        INTENT( IN )  :: MY_NSRC( NPTGRPS ) ! number of point sources
      CHARACTER( * ), INTENT( OUT ) :: LNAME( NPTGRPS )   ! layer fractions logical file name

C Local variables:

      INTEGER J, N
      CHARACTER( 256 ) :: XMSG      !  buffer for M3EXIT() messages
      CHARACTER(  16 ) :: PNAME = 'OPENLAYOUT'   ! procedure name

      CHARACTER(   3 ) :: CMYPE
      CHARACTER(   4 ) :: POSTSTR = '.ncf'
      CHARACTER(  16 ) :: PLAY_BASE
      CHARACTER(  96 ) :: PLAY_EQ
      CHARACTER(  96 ) :: EQNAME
      CHARACTER( 256 ) :: BASE

      INTEGER :: LOGDEV, IOS

C-----------------------------------------------------------------------

      LOGDEV = SETUP_LOGDEV()

C Set up and open diagnostic file(s)

      FTYPE3D = GRDDED3
      SDATE3D = SDATE
      STIME3D = STIME
      TSTEP3D = TSTEP

      NROWS3D = EMLAYS
      NLAYS3D = 1
      NTHIK3D = 1
      NVARS3D = 2
      GDTYP3D = IMISS3
      P_ALP3D = AMISS3
      P_BET3D = AMISS3
      P_GAM3D = AMISS3
      XCENT3D = AMISS3
      YCENT3D = AMISS3
      XORIG3D = AMISS3
      YORIG3D = AMISS3
      XCELL3D = AMISS3
      YCELL3D = AMISS3
      VGTYP3D = IMISS3
      VGTOP3D = AMISS3
      GDNAM3D = ' '

      VNAME3D = ' '     ! array
      VTYPE3D = 0       ! array
      UNITS3D = ' '     ! array
      VDESC3D = ' '     ! array
      VNAME3D( 1 ) = 'LFRAC'
      VTYPE3D( 1 ) = M3REAL
      UNITS3D( 1 ) = 'none'
      VDESC3D( 1 ) = 'Fraction of plume emitted into layer'
      VNAME3D( 2 ) = 'ZPLUME'
      VTYPE3D( 2 ) = M3REAL
      UNITS3D( 2 ) = 'none'
      VDESC3D( 2 ) = 'plume elevation above ground'

      VGLVS3D = 0.0     ! array
      J = LBOUND( VGLVS3D, 1 )
      VGLVS3D( J:J+EMLAYS ) = VGLVSXG( 0:EMLAYS )  ! array
      VGTYP3D = -1
      VGTOP3D = 0.0

!     UPNAM3D = 'FAKE_GRIDDED'   ! Mandatory!
      UPNAM3D = 'DISTR_GRIDDED'   ! Mandatory!

      FDESC3D = ' '  ! array
      FDESC3D( 1 ) = 'By-source hourly plume rise layer fractions'
      FDESC3D( 2 ) = '/FROM/ '    // PNAME
      FDESC3D( 3 ) = '/VERSION/ Inline'

      WRITE( CMYPE,'(I3.3)' ) MYPE
      DO N = 1, NPTGRPS

         IF ( MY_NSRC( N ) .GT. 0 ) THEN
            NCOLS3D = MY_NSRC( N )

            WRITE( PLAY_BASE,'( "PLAY_BASE", I2.2, "_" )' ) N
            XMSG = 'Logical name for LAYER FRACTIONS MATRIX'
            CALL ENVSTR( PLAY_BASE, XMSG, ' ', BASE, IOS )
            IF ( IOS .NE. 0 ) THEN
               XMSG = 'Error retrieving Layer Fractions File base name env var'
               CALL M3EXIT( PNAME, SDATE, STIME, XMSG, XSTAT1 )
            END IF

            PLAY_EQ = BASE( 1:TRIMLEN( BASE ) ) // '_' // CMYPE // POSTSTR
            WRITE( LNAME( N ),'( "PLAY_", I2.2 )' ) N
            CALL SETENVVAR ( LNAME( N ), PLAY_EQ )
            IF ( .NOT. OPEN3( LNAME( N ), FSNEW3, UPNAM3D ) ) THEN
               XMSG = 'Could not open '
     &              // FDESC3D( 1 )( 1:TRIMLEN( FDESC3D( 1 ) ) ) // ' file'
               CALL M3EXIT( PNAME, SDATE, STIME, XMSG, XSTAT1 )
            END IF

         ELSE

            WRITE( LOGDEV,1009 ) N

         END IF

      END DO

      RETURN

1009  FORMAT( /5X, 'No sources in this subdomain for file sector', I2
     &        /5X, 'Therefore no plume layer fraction file created' )

      END SUBROUTINE OPENLAYOUT
