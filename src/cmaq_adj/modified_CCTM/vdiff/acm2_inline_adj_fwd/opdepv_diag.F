
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/opdepv_diag.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPDEPV_DIAG ( JDATE, JTIME, TSTEP,
     &                         N_GDEPV_NAMES, GDEPV_NAMES,
     &                         N_ADEPV_NAMES, ADEPV_NAMES, 
     &                         N_SDEPV_NAMES, SDEPV_NAMES)

C 2006-Jun-30 - created by W Hutzell
C 2006-Dec-05 - modified by J.Young
C 2010-Jan-29 - removed unneccesary CLOSE3 statement for CTM_DEPV_DIAG file
C               to eliminate potential MPI race condition (D.Wong)

      USE HGRD_DEFN           ! horizontal grid specifications
      USE DEPVVARS

      IMPLICIT NONE
 
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations

C Arguments:

      INTEGER         JDATE            ! current model date, coded YYYYDDD
      INTEGER         JTIME            ! current model time, coded HHMMSS
      INTEGER         TSTEP            ! output time step

C number and names of gas deposition velocities
      INTEGER            N_GDEPV_NAMES
      CHARACTER( 16 ) :: GDEPV_NAMES( : )

C number and names of aerosol deposition velocities
      INTEGER            N_ADEPV_NAMES
      CHARACTER( 16 ) :: ADEPV_NAMES( : )

C number and names of (SOA+aerosol) species deposition velocities
      INTEGER            N_SDEPV_NAMES
      CHARACTER( 16 ) :: SDEPV_NAMES( : )

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'OPDDEP_DIAG'
      CHARACTER( 96 ) :: XMSG = ' '
      INTEGER         :: LOGDEV

      INTEGER      V, N, L     ! loop induction variables

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

C Try to open existing file for update

      IF ( .NOT. OPEN3( CTM_DEPV_DIAG, FSRDWR3, PNAME ) ) THEN

         XMSG = 'Could not open ' // CTM_DEPV_DIAG // ' file for update - '
     &        // 'try to open new'
         CALL M3MESG( XMSG )

C Set output file characteristics based on conc file

         IF ( .NOT. OPEN3( CTM_CONC_1, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open ' // TRIM( CTM_CONC_1 )
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         IF ( .NOT. DESC3( CTM_CONC_1  ) ) THEN
            XMSG = 'Could not get ' // TRIM( CTM_CONC_1 ) // ' file description'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         IF ( .NOT. CLOSE3( CTM_CONC_1 ) ) THEN
            XMSG = 'Could not close ' // TRIM( CTM_CONC_1 )
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         FDESC3D = ' '   ! array

         L = LEN_TRIM( DEPV_METHOD )
         FDESC3D( 1 ) = 'Instanteous values of dry deposition velocities:'
     &               // ' the gas phase velocities are' 
         FDESC3D( 2 ) = 'estimated by the '
     &               //  DEPV_METHOD( 1:L ) // ' algorithm.'

         FTYPE3D = GRDDED3

         SDATE3D = JDATE 
         STIME3D = JTIME 
         TSTEP3D = TSTEP 

!        WRITE( LOGDEV,* ) PNAME, ': TSTEP_DEPV_FILE = ', TSTEP_DEPV_FILE

         NVARS3D = N_GDEPV_NAMES + N_ADEPV_NAMES + N_SDEPV_NAMES
!        NCOLS3D = GL_NCOLS
!        NROWS3D = GL_NROWS
         NLAYS3D =     1
!        NTHIK3D =     1
!        GDTYP3D = GDTYP_GD
!        P_ALP3D = P_ALP_GD
!        P_BET3D = P_BET_GD 
!        P_GAM3D = P_GAM_GD
!        XORIG3D = XORIG_GD
!        YORIG3D = YORIG_GD
!        XCENT3D = XCENT_GD
!        YCENT3D = YCENT_GD
!        XCELL3D = XCELL_GD
!        YCELL3D = YCELL_GD
!        VGTYP3D = VGTYP_GD
!        VGTOP3D = VGTOP_GD
!!       VGTPUN3D = VGTPUN_GD ! currently, not defined
!        DO L = 1, NLAYS3D + 1
!           VGLVS3D( L ) = VGLVS_GD( L )
!        END DO
!        GDNAM3D = GDNAME_GD
!        GDNAM3D = GRID_NAME  ! from HGRD_DEFN

         N = 0

         DO V = 1, N_GDEPV_NAMES
            N = N + 1
            VTYPE3D( V ) = M3REAL
            VNAME3D( V ) = GDEPV_NAMES( V )
            UNITS3D( V ) = 'cm/sec'
            VDESC3D( V ) = 'gas phase dry deposition velocity'
         END DO

         N = N_GDEPV_NAMES

         DO V = 1, N_ADEPV_NAMES
            N = N + 1
            VTYPE3D( N ) = M3REAL
            VNAME3D( N ) = ADEPV_NAMES( V )
            UNITS3D( N ) = 'cm/sec'
            VDESC3D( N ) = 'aerosol phase dry deposition velocity'
         END DO

         N = N_GDEPV_NAMES + N_ADEPV_NAMES

         DO V = 1, N_SDEPV_NAMES
            N = N + 1
            VTYPE3D( N ) = M3REAL
            VNAME3D( N ) = SDEPV_NAMES( V )
            UNITS3D( N ) = 'cm/sec'
            VDESC3D( N ) = 'deposition velocity for (GAS-AEROSOL) '
     &                  // ' semi-volatile compound'
         END DO

C Open diagnostic file, then close it for subsequent open by all processors

         IF ( .NOT. OPEN3( CTM_DEPV_DIAG, FSNEW3, PNAME ) ) THEN
            XMSG = 'Could not create '// TRIM( CTM_DRY_DEP_1 ) // ' file'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END IF

      RETURN
      END
