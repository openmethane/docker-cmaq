
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhg/oppt3d_diag.F,v 1.1.1.1 2010/06/14 16:03:09 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPPT3D_DIAG( PT3DNAME, SDATE, STIME, TSTEP, EMLAYS,
     &                        N_SPC_EM, STKEM_NAME,
     &                        N_SPC_PM, STKPM_NAME )

C-----------------------------------------------------------------------
 
C  Description:
C    Opens the output 3d point source emissions diagnostic file
 
C  Preconditions:
 
C  Subroutines and Functions Called:
C       I/O API 

C Revision History:
C      4 Sep 2007 J.Young: initial
C     29 Jan 2010 D. Wong: Eliminate potential race condition at the MYPE = 0
C                          to open the diagnostic file. Eliminate the BARRIER
C                          and the CLOSE3 and OPEN3.
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
 
C Includes:
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
!include      SUBST_IODECL    # I/O definitions and declarations
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC   ! file header data structure

C External Functions:

C Local Parameters:

C Arguments:
      CHARACTER( 16 ), INTENT( IN ) :: PT3DNAME ! diagnostic file logical name
      INTEGER,         INTENT( IN ) :: SDATE    ! Julian start date
      INTEGER,         INTENT( IN ) :: STIME    ! start time (HHMMSS)
      INTEGER,         INTENT( IN ) :: TSTEP    ! time step (HHMMSS)
      INTEGER,         INTENT( IN ) :: EMLAYS   ! number of emissions layers
!        REAL,         INTENT( IN ) :: VGLVSXG( 0:MXLAYS3 ) !  vertical coord values
      INTEGER,         INTENT( IN ) :: N_SPC_EM ! number of point source gas species
      CHARACTER( 16 ), INTENT( IN ) :: STKEM_NAME( N_SPC_EM ) ! non-PM species names
      INTEGER,         INTENT( IN ) :: N_SPC_PM ! number of point source PM species
      CHARACTER( 16 ), INTENT( IN ) :: STKPM_NAME( N_SPC_PM ) ! PM species names

C Local variables:

      INTEGER L, V, IOS
      CHARACTER( 256 ) :: XMSG      ! buffer for M3EXIT() messages
      CHARACTER( 256 ) :: RET_VAL   ! env var returned value
      CHARACTER(  16 ) :: PNAME = 'OPPT3DDIAG'   ! procedure name

C-----------------------------------------------------------------------

C Set output file characteristics based on GRID_CONF module and open the 3D emis file

      FTYPE3D = GRDDED3
      SDATE3D = SDATE
      STIME3D = STIME
      TSTEP3D = TSTEP
      NCOLS3D = GL_NCOLS
      NROWS3D = GL_NROWS
      NLAYS3D = EMLAYS
      NTHIK3D = 1
      NVARS3D = N_SPC_EM + N_SPC_PM
      GDTYP3D = GDTYP_GD
      P_ALP3D = P_ALP_GD
      P_BET3D = P_BET_GD
      P_GAM3D = P_GAM_GD
      XORIG3D = XORIG_GD
      YORIG3D = YORIG_GD
      XCENT3D = XCENT_GD
      YCENT3D = YCENT_GD
      XCELL3D = XCELL_GD
      YCELL3D = YCELL_GD
      VGTYP3D = VGTYP_GD
      VGTOP3D = VGTOP_GD
      GDNAM3D = GRID_NAME  ! from HGRD_DEFN

!     VGLVS3D = 0.0     ! array
!     L = LBOUND( VGLVS3D, 1 )
!     VGLVS3D( L:L+EMLAYS ) = VGLVSXG( 0:EMLAYS )  ! array

      DO L = 1, NLAYS3D + 1
         VGLVS3D( L ) = VGLVS_GD( L )
      END DO

      VNAME3D = ' '     ! array
      DO V = 1, N_SPC_EM
         VNAME3D( V ) = STKEM_NAME( V )
      END DO
      L = N_SPC_EM
      DO V = 1, N_SPC_PM
         VNAME3D( V+L ) = STKPM_NAME( V )
      END DO

      VTYPE3D = 0       ! array
      VTYPE3D( 1:NVARS3D ) = M3REAL

      UNITS3D = ' '     ! array
      DO V = 1, N_SPC_EM
         UNITS3D( V ) = 'mol/s'
      END DO
      L = N_SPC_EM
      DO V = 1, N_SPC_PM
         UNITS3D( V+L ) = 'g/s'
      END DO

      VDESC3D = ' '     ! array
      DO V = 1, N_SPC_EM
         VDESC3D( V ) = 'Model species: ' // TRIM( STKEM_NAME( V ) )
      END DO
      L = N_SPC_EM
      DO V = 1, N_SPC_PM
         VDESC3D( V+L) = 'Primary particulate species: ' // TRIM( STKPM_NAME( V ) )
      END DO

      FDESC3D = ' '  ! array
      FDESC3D( 1 ) = 'gridded hourly point source layered emissions'
      FDESC3D( 2 ) = '/From/ ' // PNAME
      FDESC3D( 3 ) = '/Version/ Inline'

!     MESG = 'Enter logical name for Point Source 3D Emissions File'
!     PT3DNAME = PROMPTMFILE( MESG, FSNEW3, 'PT3DEMIS', PNAME )
!     PT3DNAME = 'PT3DEMIS_DIAG'

!     XMSG = 'logical name for Point Source 3D Emissions File'
!     CALL ENVSTR( PT3DNAME, XMSG, ' ', RET_VAL, IOS )
!     IF ( IOS .NE. 0 ) THEN
!     XMSG = 'Error retrieving Point Source 3D Emissions File name env var'
!     CALL M3EXIT( PNAME, SDATE, STIME, XMSG, XSTAT1 )
!     END IF

C create header

      IF ( MYPE .EQ. 0 ) THEN   ! open new
         IF ( .NOT. OPEN3( PT3DNAME, FSNEW3, PNAME ) ) THEN
            XMSG = 'Could not open ' // TRIM( PT3DNAME )
            CALL M3EXIT( PNAME, SDATE, STIME, XMSG, XSTAT1 )
         END IF
      END IF

      RETURN

      END SUBROUTINE OPPT3D_DIAG

