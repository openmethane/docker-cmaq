
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhg/opddep.F,v 1.1.1.1 2010/06/14 16:03:09 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPDDEP ( JDATE, JTIME, TSTEP, N_SPC_DDEP )

C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN
C   03 Sep 01 David Wong: for new pario
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
 
!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_DDEP   ! gas chem dry dep species and map table
      INCLUDE SUBST_AE_DDEP   ! aerosol dry dep species and map table
      INCLUDE SUBST_NR_DDEP   ! non-react dry dep species and map table
      INCLUDE SUBST_TR_DDEP   ! tracer dry dep species and map table
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID  ! coordinate and domain definitions (req IOPARMS)

C Arguments:

      INTEGER      JDATE      ! current model date, coded YYYYDDD
      INTEGER      JTIME      ! current model time, coded HHMMSS
      INTEGER      TSTEP      ! output time step
      INTEGER      N_SPC_DDEP ! no. of dry dep species

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'OPDDEP'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER      V, N, L     ! loop induction variables
      INTEGER      STRT, FINI  ! loop induction variables

C-----------------------------------------------------------------------

C Try to open existing file for update

      IF ( .NOT. OPEN3( CTM_DRY_DEP_1, FSRDWR3, PNAME ) ) THEN

         XMSG = 'Could not open ' // CTM_DRY_DEP_1 // ' file for update - '
     &        // 'try to open new'
         CALL M3MESG( XMSG )

C Set output file characteristics based on COORD.EXT and open the dry dep file

         FTYPE3D = GRDDED3
         SDATE3D = JDATE
         STIME3D = JTIME
         TSTEP3D = TSTEP
         CALL NEXTIME( SDATE3D, STIME3D, TSTEP3D ) !  start the next hour

         NVARS3D = N_SPC_DDEP
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
!        VGTPUN3D = VGTPUN_GD ! currently, not defined
         DO L = 1, NLAYS3D + 1
            VGLVS3D( L ) = VGLVS_GD( L )
            END DO
!        GDNAM3D = GDNAME_GD
         GDNAM3D = GRID_NAME  ! from HGRD_DEFN

         N = 0
         STRT = 1
         FINI = N_GC_DDEP
         DO V = STRT, FINI
            N = N + 1
            VTYPE3D( V ) = M3REAL
            VNAME3D( V ) = GC_DDEP( N )
            UNITS3D( V ) = 'kg/hectare'
            VDESC3D( V ) = 'hourly dry deposition values'
            END DO

         N = 0
         STRT = N_GC_DDEP + 1
         FINI = N_GC_DDEP + N_AE_DDEP
         DO V = STRT, FINI
            N = N + 1
            VTYPE3D( V ) = M3REAL
            VNAME3D( V ) = AE_DDEP( N )
            IF ( AE_DDEP( N )( 1:3 ) .EQ. 'NUM' ) THEN
               UNITS3D( V ) = 'number/hectare'
               ELSE IF ( AE_DDEP( N )( 1:3 ) .EQ. 'SRF' ) THEN
               UNITS3D( V ) = 'm**2/hectare'
               ELSE
               UNITS3D( V ) = 'kg/hectare'
               END IF
            VDESC3D( V ) = 'hourly dry deposition values'
            END DO

         N = 0
         STRT = N_GC_DDEP + N_AE_DDEP + 1
         FINI = N_GC_DDEP + N_AE_DDEP + N_NR_DDEP
         DO V = STRT, FINI
            N = N + 1
            VTYPE3D( V ) = M3REAL
            VNAME3D( V ) = NR_DDEP( N )
            UNITS3D( V ) = 'kg/hectare'
            VDESC3D( V ) = 'hourly dry deposition values'
            END DO

         N = 0
         STRT = N_GC_DDEP + N_AE_DDEP + N_NR_DDEP + 1
         FINI = N_GC_DDEP + N_AE_DDEP + N_NR_DDEP + N_TR_DDEP
         DO V = STRT, FINI
            N = N + 1
            VTYPE3D( V ) = M3REAL
            VNAME3D( V ) = TR_DDEP( N )
            UNITS3D( V ) = '----'
            VDESC3D( V ) = 'hourly dry deposition values'
            END DO

         FDESC3D( 1 ) = 'hourly 1-layer cross-point RADM dry deposition data'
         DO L = 2, MXDESC3
            FDESC3D( L ) = ' '
            END DO

C    Open dry deposition file

         IF ( .NOT. OPEN3( CTM_DRY_DEP_1, FSNEW3, PNAME ) ) THEN
            XMSG = 'Could not create '// CTM_DRY_DEP_1 // ' file'
            CALL M3EXIT( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
            END IF

         END IF

      RETURN
      END
