
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero4/opdiam.F,v 1.1.1.1 2010/06/14 16:02:58 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPDIAM ( JDATE, JTIME, TSTEP )

C Revision history
C   ??? Frank Binkowski
C   8 Sep 01 J.Young: dyn alloc - Use HGRD_DEFN
C   03 Sep 01 David Wong: for new pario
C   13 May 04 P.Bhave: added RH to species list; removed M0 and dry M2
C                      conc's from species list
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   18 Jul 05 P.Bhave: Added mass fractions of each mode < 2.5um to output list
C   06 Apr 06 P.Bhave: Added GAMMA_N2O5 to output list
C   25 May 06 P.Bhave: Changed units of all DG variables from m to um, as 
C                      suggested by Dr. Bill Hutzell

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

      CHARACTER( 16 ), SAVE :: PNAME = 'OPDIAM'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER      L          ! loop induction variables

C-----------------------------------------------------------------------

C Try to open existing file for update

      IF ( .NOT. OPEN3( CTM_DIAM_1, FSRDWR3, PNAME ) ) THEN

         XMSG = 'Could not open ' // CTM_DIAM_1 // ' file for update - '
     &        // 'try to open new'
         CALL M3MESG( XMSG )

C Set output file characteristics based on COORD.EXT and open the aerosol
C diagnostic file

         FTYPE3D = GRDDED3
         SDATE3D = JDATE
         STIME3D = JTIME
         TSTEP3D = TSTEP
         CALL NEXTIME ( SDATE3D, STIME3D, TSTEP3D ) !  start the next hour

         NVARS3D = 19
         NCOLS3D = GL_NCOLS
         NROWS3D = GL_NROWS
         NLAYS3D = NLAYS
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

         VTYPE3D( 1 ) = M3REAL
         VNAME3D( 1 ) = 'STDEVATKN'
         UNITS3D( 1 ) = ' '
         VDESC3D( 1 ) = 'Aitken mode standard deviation (dry and wet)'

         VTYPE3D( 2 ) = M3REAL
         VNAME3D( 2 ) = 'STDEVACC'
         UNITS3D( 2 ) = ' '
         VDESC3D( 2 ) = 'Accumulation mode standard deviation (dry and wet)'
            
         VTYPE3D( 3 ) = M3REAL
         VNAME3D( 3 ) = 'DGATKN_DRY'
         UNITS3D( 3 ) = 'um'
         VDESC3D( 3 ) = 'Aitken mode mean diameter (dry)'

         VTYPE3D( 4 ) = M3REAL
         VNAME3D( 4 ) = 'DGACC_DRY'
         UNITS3D( 4 ) = 'um'
         VDESC3D( 4 ) = 'Accumulation mode mean diameter (dry)'
         
         VTYPE3D( 5 ) = M3REAL
         VNAME3D( 5 ) = 'DGATKN_WET'
         UNITS3D( 5 ) = 'um'
         VDESC3D( 5 ) = 'Aitken mode mean diameter (wet)'

         VTYPE3D( 6 ) = M3REAL
         VNAME3D( 6 ) = 'DGACC_WET'
         UNITS3D( 6 ) = 'um'
         VDESC3D( 6 ) = 'Accumulation mode mean diameter (wet)'

         VTYPE3D( 7 ) = M3REAL
         VNAME3D( 7 ) = 'DGCOR'
         UNITS3D( 7 ) = 'um'
         VDESC3D( 7 ) = 'Coarse mode mean diameter'

         VTYPE3D( 8 ) = M3REAL
         VNAME3D( 8 ) = 'M2ATKN_WET'
         UNITS3D( 8 ) = 'm**2/m**3'
         VDESC3D( 8 ) = 'Aitken mode 2nd moment (wet)'

         VTYPE3D( 9 ) = M3REAL
         VNAME3D( 9 ) = 'M2ACC_WET'
         UNITS3D( 9 ) = 'm**2/m**3'
         VDESC3D( 9 ) = 'Accumulation mode 2nd moment (wet)'

         VTYPE3D( 10 ) = M3REAL
         VNAME3D( 10 ) = 'M3ATKN_DRY'
         UNITS3D( 10 ) = 'm**3/m**3'
         VDESC3D( 10 ) = 'Aitken mode 3rd moment (dry)'

         VTYPE3D( 11 ) = M3REAL
         VNAME3D( 11 ) = 'M3ACC_DRY'
         UNITS3D( 11 ) = 'm**3/m**3'
         VDESC3D( 11 ) = 'Accumulation mode 3rd moment (dry)'

         VTYPE3D( 12 ) = M3REAL
         VNAME3D( 12 ) = 'M3ATKN_WET'
         UNITS3D( 12 ) = 'm**3/m**3'
         VDESC3D( 12 ) = 'Aitken mode 3rd moment (wet)'

         VTYPE3D( 13 ) = M3REAL
         VNAME3D( 13 ) = 'M3ACC_WET'
         UNITS3D( 13 ) = 'm**3/m**3'
         VDESC3D( 13 ) = 'Accumulation mode 3rd moment (wet)'

         VTYPE3D( 14 ) = M3REAL
         VNAME3D( 14 ) = 'M3COR'
         UNITS3D( 14 ) = 'm**3/m**3'
         VDESC3D( 14 ) = 'Coarse mode 3rd moment'

         VTYPE3D( 15 ) = M3REAL
         VNAME3D( 15 ) = 'RH'
         UNITS3D( 15 ) = ' '
         VDESC3D( 15 ) = 'Relative humidity'

         VTYPE3D( 16 ) = M3REAL
         VNAME3D( 16 ) = 'PM25AT'
         UNITS3D( 16 ) = ' '
         VDESC3D( 16 ) = 'Fine fraction of Aitken mode'

         VTYPE3D( 17 ) = M3REAL
         VNAME3D( 17 ) = 'PM25AC'
         UNITS3D( 17 ) = ' '
         VDESC3D( 17 ) = 'Fine fraction of accumulation mode'

         VTYPE3D( 18 ) = M3REAL
         VNAME3D( 18 ) = 'PM25CO'
         UNITS3D( 18 ) = ' '
         VDESC3D( 18 ) = 'Fine fraction of coarse mode'

         VTYPE3D( 19 ) = M3REAL
         VNAME3D( 19 ) = 'GAMMA_N2O5'
         UNITS3D( 19 ) = ' '
         VDESC3D( 19 ) = 'N2O5 Heterogeneous reaction probability'

         FDESC3D( 1 ) = 'aerosol distribution parameters'
         DO L = 2, MXDESC3
            FDESC3D( L ) = ' '
            END DO

C Open the aerosol diagnostic file

         IF ( .NOT. OPEN3( CTM_DIAM_1, FSNEW3, PNAME ) ) THEN
            XMSG = 'Could not create '// CTM_DIAM_1 // ' file'
            CALL M3EXIT ( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
            END IF

         END IF

      RETURN
      END
