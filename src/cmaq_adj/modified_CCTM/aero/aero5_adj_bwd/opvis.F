
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_adj_bwd/opvis.F,v 1.1 2010/11/18 14:51:26 sergey Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPVIS ( JDATE, JTIME, TSTEP )

C Revision history
C   ??? Frank Binkowski
C   8 Sep 01 J.Young: dyn alloc - Use HGRD_DEFN
C   03 Sep 01 David Wong: for new pario
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
 
!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID  ! coordinate and domain definitions (req IOPARMS)

C...Arguments:

      INTEGER      JDATE      ! current model date, coded YYYYDDD
      INTEGER      JTIME      ! current model time, coded HHMMSS
      INTEGER      TSTEP      ! output time step

C...Local variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'OPVIS'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER      L          ! loop induction variables

C-----------------------------------------------------------------------

C Try to open existing file for update

      IF ( .NOT. OPEN3( CTM_VIS_1, FSRDWR3, PNAME ) ) THEN

         XMSG = 'Could not open ' // CTM_VIS_1 // ' file for update - '
     &        // 'try to open new'
         CALL M3MESG( XMSG )

C Set output file characteristics based on COORD.EXT and open the dry dep file

         FTYPE3D = GRDDED3
         SDATE3D = JDATE
         STIME3D = JTIME
         TSTEP3D = TSTEP
         CALL NEXTIME ( SDATE3D, STIME3D, TSTEP3D ) !  start the next hour

         NVARS3D = 4
         NCOLS3D = GL_NCOLS
         NROWS3D = GL_NROWS
         NLAYS3D =     1
         NTHIK3D =     1
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
!        VGTPUN3D = VGTPUN_GD   ! currently, not defined

         DO L = 1, NLAYS3D + 1
            VGLVS3D( L ) = VGLVS_GD( L )
            END DO

!        GDNAM3D = GDNAME_GD
         GDNAM3D = GRID_NAME  ! from HGRD_DEFN

         VTYPE3D( 1 ) = M3REAL
         VNAME3D( 1 ) = 'DCV_Mie'
         UNITS3D( 1 ) = 'deciview'
         VDESC3D( 1 ) = 'hourly visual range'

         VTYPE3D( 2 ) = M3REAL
         VNAME3D( 2 ) = 'EXT_Mie'
         UNITS3D( 2 ) = '1/km'
         VDESC3D( 2 ) = 'hourly extinction coefficient'
            
         VTYPE3D( 3 ) = M3REAL
         VNAME3D( 3 ) = 'DCV_Recon'
         UNITS3D( 3 ) = 'deciview'
         VDESC3D( 3 ) = 'hourly visual range'

         VTYPE3D( 4 ) = M3REAL
         VNAME3D( 4 ) = 'EXT_Recon'
         UNITS3D( 4 ) = '1/km'
         VDESC3D( 4 ) = 'hourly extinction coefficient'

         FDESC3D( 1 ) = 'hourly 1-layer cross-point Visibility Information'
         DO L = 2, MXDESC3
            FDESC3D( L ) = ' '
            END DO

C Open visibility file new

         IF ( .NOT. OPEN3( CTM_VIS_1, FSNEW3, PNAME ) ) THEN
            XMSG = 'Could not create '// CTM_VIS_1 // ' file'
            CALL M3EXIT( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
            END IF

      END IF

      RETURN
      END
