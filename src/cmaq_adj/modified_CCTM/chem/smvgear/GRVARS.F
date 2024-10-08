
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/smvgear/GRVARS.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE GRVARS

C smvgear variables that are dependent on NCOLS, NROWS

C Preconditions: GRID_INIT() called from PAR_INIT, which is called from DRIVER

C Revision History: J.Young 17 Aug 01: create
C                   J.Young 31 Jan 05: get MXCELLS from dyn alloc horizontal
C                   & vertical domain specifications module (GRID_CONF)
C                   29 Jul 05 WTH: added variable used by degrade routines.
C.......................................................................

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Column index of ordered cells
      INTEGER, ALLOCATABLE, SAVE :: CCOL( : )
C Row index for ordered cells
      INTEGER, ALLOCATABLE, SAVE :: CROW( : )
C Layer index of ordered cells
      INTEGER, ALLOCATABLE, SAVE :: CLEV( : )
C Cell number offset for each block
      INTEGER, ALLOCATABLE, SAVE :: BLKCNO( : )
C Number of cells in each block
      INTEGER, ALLOCATABLE, SAVE :: BLKLEN( : )

C Original cell number of ordered cell index
      INTEGER, ALLOCATABLE, SAVE :: NORDCELL( : )
C Estimated stiffness of each cell
      REAL( 8 ), ALLOCATABLE, SAVE :: ERRMX2( : )

      LOGICAL, SAVE :: CALL_DEG = .FALSE.    ! WTH: SWITCH for calling DEGRADE routine

      CONTAINS
         FUNCTION CELLVAR_ALLOC () RESULT ( SUCCESS )

         INCLUDE SUBST_IODECL      ! I/O definitions and declarations

         LOGICAL :: SUCCESS
         LOGICAL, SAVE :: FIRSTIME = .TRUE.
         INTEGER :: ALST
         CHARACTER( 96 ) :: XMSG = ' '

C This function is expected to be called only once - at startup

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            SUCCESS = .TRUE.

            ALLOCATE ( CCOL( MXCELLS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** CCOL Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ALLOCATE ( CROW( MXCELLS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** CROW Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ALLOCATE ( CLEV( MXCELLS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** CLEV Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ALLOCATE ( BLKCNO( MXBLKS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** BLKCNO Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ALLOCATE ( BLKLEN( MXBLKS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** BLKLEN Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ALLOCATE ( NORDCELL( MXCELLS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** NORDCELL Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ALLOCATE ( ERRMX2( MXCELLS ), STAT = ALST )
            IF ( ALST .NE. 0 ) THEN
               XMSG = '*** ERRMX2 Memory allocation failed'
               CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
               SUCCESS = .FALSE.
               RETURN
               END IF

            ELSE   ! called more than once

            XMSG = 'Horizontal domain dependent variables already allocated'
            CALL M3WARN ( 'CELLVAR_ALLOC', 0, 0, XMSG )
            SUCCESS = .FALSE.; RETURN

            END IF   ! FIRSTIME

         RETURN
         END FUNCTION CELLVAR_ALLOC

      END MODULE GRVARS
