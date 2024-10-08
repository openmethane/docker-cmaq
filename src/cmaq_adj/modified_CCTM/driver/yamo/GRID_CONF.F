
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo/GRID_CONF.F,v 1.1.1.1 2010/06/14 16:03:04 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE GRID_CONF

C Define the horizontal grid and vertical layer structure by inheriting
C HGRD_DEFN and VGRD_DEFN. In addition, set up the computational blocking
C parameters.
C VGRD_DEFN replaces COORD.EXT and VGRD.EXT
C Revision History: Jeff Young 31 Jen 05: created
C.......................................................................

      USE HGRD_DEFN             ! horizontal domain specifications
      USE VGRD_DEFN             ! vertical layer specifications

      IMPLICIT NONE

!     INTEGER, PARAMETER :: BLKSIZE = 500
      INTEGER, PARAMETER :: BLKSIZE = 50
      INTEGER, SAVE :: MXCELLS
      INTEGER, SAVE :: MXBLKS

      CONTAINS

         FUNCTION GRID_INIT ( NPROCS, MYID ) RESULT ( SUCCESS )

         INTEGER :: NPROCS
         INTEGER :: MYID
         LOGICAL :: SUCCESS

         INTEGER, SAVE :: LOGDEV
         LOGICAL, SAVE :: FIRSTIME = .TRUE.
         CHARACTER( 96 ) :: MSG = ' '

         INTEGER, EXTERNAL :: SETUP_LOGDEV

C-----------------------------------------------------------------------

C This function is expected to be called only once - at startup

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            LOGDEV = SETUP_LOGDEV()
            SUCCESS = .TRUE.

C Set up horizontal domain and calculate processor-to-subdomain maps
            IF ( .NOT. HGRD_INIT ( NPROCS, MYID ) ) THEN
               MSG = '*** Failure defining horizontal domain'
               CALL M3WARN ( 'GRID_INIT', 0, 0, MSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

C define vertical layer structure
            IF ( .NOT. VGRD_INIT() ) THEN
               MSG = '*** Failure defining vertical structure'
               CALL M3WARN ( 'GRID_INIT', 0, 0, MSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            MXCELLS = NCOLS * NROWS * NLAYS
            MXBLKS  = 1 + ( MXCELLS - 1 ) / BLKSIZE

            ELSE
            MSG = 'Domain configuration already defined'
            CALL M3WARN ( 'GRID_INIT', 0, 0, MSG )
            SUCCESS = .FALSE.; RETURN

            END IF   ! FIRSTIME

         RETURN
         END FUNCTION GRID_INIT

      END MODULE GRID_CONF


