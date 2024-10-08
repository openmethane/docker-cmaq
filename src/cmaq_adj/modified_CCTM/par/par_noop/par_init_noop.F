
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/par/par_noop/par_init_noop.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE PAR_INIT( COLROW, NSPCS, CLOCK, IERROR )
C.....................................................................
C
C  PURPOSE:    No-operation version of parallel initialization, for
C              running in serial configuration. Local processor
C              dimensions are set to the one-processor values.
C              
C
C  REVISION HISTORY: 
C       Original version 10/1998 by Al Bourgeois for merging parallel
C                and serial codes.
C       Modified 01/1998 by Al Bourgeois to fix subscripts in arrays 
C                ROWSX_PE and COLSX_PE. 
C       Modified 11/2001 by Jeff Young for dynamic alloc version
C       Modified 31/01/05 J.Young: dyn alloc - establish both horizontal &
C                vertical domain specifications in one module
C                16 Jan 06 J.Young: possible PA grid init
C
C
C  ARGUMENT LIST DESCRIPTION:
C  IN:
C     CHARACTER( 2 ), INTENT( IN ) :: COLROW ! columns/rows arg list check
C     INTEGER NSPCS             ! Number of species in CGRID.
C
C  OUT:
C     REAL    CLOCK             ! Wall-clock time (sec) at MPI initialization.
C     INTEGER IERROR            ! Error code.
C
C
C  LOCAL VARIABLE DESCRIPTION:  see below
C
C  CALLS:  None 
C
C........................................................................

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE PAGRD_DEFN            ! PA horiz domain specs

      IMPLICIT  NONE

C Include Files

      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PACTL_ID    ! PA control parameters

C Arguments

      CHARACTER( 2 ), INTENT( IN ) :: COLROW ! columns/rows arg list check
      INTEGER, INTENT( IN )        :: NSPCS  ! Number of species in CGRID
      REAL, INTENT( OUT ) :: CLOCK ! Wall-clock time (sec) at MPI initialization
      INTEGER, INTENT( OUT )       :: IERROR ! Error code

C Local Variables

      INTEGER NPROCS
      CHARACTER( 96 ) :: MSG

C Set error code to zero.
      IERROR = 0

      NPROCS = 1
      MYPE = 0
      CLOCK = 0

C Set up horizontal domain and calculate processor-to-subdomain maps
C and define vertical layer structure (in module GRID_CONF)

      IF ( .NOT. GRID_INIT ( NPROCS, MYPE ) ) THEN
         MSG = '*** Failure defining domain configuration'
         CALL M3WARN ( 'PAR_INIT', 0, 0, MSG )
         IERROR = 1
         RETURN
         END IF

C Set up horizontal domain and calculate processor-to-subdomain maps for
C process analysis, if required

      IF ( LIPR .OR. LIRR ) THEN
         IF ( .NOT. PAGRD_INIT( NPROCS, MYPE ) ) THEN
            MSG = '*** Failure defining PA domain configuration'
            CALL M3WARN ( 'PAR_INIT', 0, 0, MSG )
            IERROR = 1
            RETURN
            END IF
         END IF

      RETURN
      END
