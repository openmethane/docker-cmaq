
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_fwd/VGRD_DEFN.F,v 1.1 2010/11/17 19:22:55 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE VGRD_DEFN

C Define the vertical layer structure (replaces COORD.EXT and VGRD.EXT)
C Revision History: Jeff Young: 6 Dec 04: created
C                              26 Jul 06: comments in line with PARMS3.EXT
C.......................................................................

      IMPLICIT NONE

      INTEGER, SAVE :: NLAYS = 0

      REAL, ALLOCATABLE, SAVE :: VGLVS_GD( : )  ! vert layer surface values
      REAL, ALLOCATABLE, SAVE :: X3FACE_GD( : ) ! mono. incr. vert coord values
                                                ! = 1 - VGLVS_GD
      CHARACTER( 16 ), SAVE :: GDNAME_GD

C The definitions and declarations in this F90 Module follow those
C given in the I/O-API include file FDESC3C.EXT and replace similar
C definitions. This F90 Module is compatible with FDESC3C.EXT.

C VGTYP_GD:
C The vertical grid type:
C  1: VGSGPH3 : hydrostatic sigma-P
C  2: VGSGPN3 : non-h sigma-P
C  3: VGSIGZ3 : sigma-Z
C  4: VGPRES3 : pressure (mb)
C  5: VGZVAL3 : Z (m) (above sea lvl)
C  6: VGHVAL3 : H (m) (above ground)
C  7: VGWRFEM : WRF mass-core sigma
C  8: VGWRFNM : WRF NMM
C  9: IMISS3 for vertical coordinates not stored in VGLVSD
C            (e.g., temporally or spatially changing vertical coordinates)

      INTEGER, SAVE :: VGTYP_GD

C VGTPUN_GD:
C The units of the vertical coordinate top.

      CHARACTER( 16 ), SAVE :: VGTPUN_GD

C VGTOP_GD:
C The value for the model top used in the definition of the sigma coordinate
C systems in the VGTPUN_GD units.

      REAL, SAVE :: VGTOP_GD

C For sigma-P, the relationship between pressure levels P and sigma-P is
C given by:
C    sigma-P = ( P - VGTOP_GD ) / (P_srf - VGTOP_GD ),
C    where P_srf is the surface pressure.

C VGLVUN_GD:
C The units of the vertical coordinate surface values,

      CHARACTER( 16 ) :: VGLVUN_GD

C-----------------------------------------------------------------------

      CONTAINS

         FUNCTION VGRD_INIT ( ) RESULT ( SUCCESS )

            IMPLICIT NONE

            INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
            INCLUDE SUBST_IOFDESC     ! file header data structure
            INCLUDE SUBST_IODECL      ! I/O definitions and declarations

            LOGICAL :: SUCCESS

!           INTEGER, SAVE :: LOGDEV
            LOGICAL, SAVE :: FIRSTIME = .TRUE.
            CHARACTER( 96 ) :: XMSG = ' '
            CHARACTER( 16 ) :: LAYER_FILE = 'LAYER_FILE'

            INTEGER L, ALLOCSTAT

C This function is expected to be called only once - at startup

            IF ( FIRSTIME ) THEN
               FIRSTIME = .FALSE.
!              LOGDEV = INIT3()
               SUCCESS = .TRUE.

C get vert layer definition from LAYER_FILE

               IF ( .NOT. OPEN3( LAYER_FILE, FSREAD3, 'VGRD_INIT' ) ) THEN
                  XMSG = 'Could not open '// LAYER_FILE
                  CALL M3WARN ( 'VGRD_INIT', 0, 0, XMSG )
                  SUCCESS = .FALSE.; RETURN
                  END IF

               IF ( .NOT. DESC3( LAYER_FILE ) ) THEN
                  XMSG = 'Could not get ' // LAYER_FILE // ' file description'
                  CALL M3WARN ( 'VGRD_INIT', 0, 0, XMSG )
                  SUCCESS = .FALSE.; RETURN
                  END IF

               GDNAME_GD = GDNAM3D

               NLAYS = NLAYS3D

               VGTYP_GD  = VGTYP3D
               VGTOP_GD  = VGTOP3D
               VGTPUN_GD = 'Pa'
               VGLVUN_GD = 'none'

               ALLOCATE ( VGLVS_GD( NLAYS + 1 ), STAT = ALLOCSTAT )
               ALLOCATE ( X3FACE_GD( 0:NLAYS ), STAT = ALLOCSTAT )
               IF ( ALLOCSTAT .NE. 0 ) THEN
                  XMSG = 'Failure allocating VGLVS_GD, X3FACE_GD'
                  CALL M3WARN ( 'VGRD_INIT', 0, 0, XMSG )
                  SUCCESS = .FALSE.; RETURN
                  END IF

               DO L = 1, NLAYS + 1
                  VGLVS_GD( L ) = VGLVS3D( L )
                  END DO

C set layer layer thickenesses

               DO L = 0, NLAYS
                  X3FACE_GD( L ) = 1.0 - VGLVS_GD( L + 1 )
                  END DO

            ELSE
            XMSG = 'Vertical layer structure already defined'
            CALL M3WARN ( 'VGRD_INIT', 0, 0, XMSG )
            SUCCESS = .FALSE.; RETURN

            END IF   ! FIRSTIME

         RETURN
         END FUNCTION VGRD_INIT

      END MODULE VGRD_DEFN

C  21   VGLVS_GD   X3FACE_GD
C  
C   1     1.0       0.0
C   2     0.995     0.005
C   3     0.99      0.01
C   4     0.985     0.015
C   5     0.98      0.02
C   6     0.97      0.03
C   7     0.96      0.04
C   8     0.945     0.055
C   9     0.93      0.07
C  10     0.91      0.09
C  11     0.89      0.11
C  12     0.865     0.135
C  13     0.84      0.16
C  14     0.78      0.22
C  15     0.7       0.3
C  16     0.6       0.4
C  17     0.5       0.5
C  18     0.4       0.6
C  19     0.3       0.7
C  20     0.2       0.8
C  21     0.1       0.9
C  22     0.0       1.0
