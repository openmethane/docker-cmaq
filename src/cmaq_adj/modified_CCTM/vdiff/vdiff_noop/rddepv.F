**********************************************************************
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/rddepv.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE RDDEPV ( JDATE, JTIME, SGSPEC, DEPV )

C-----------------------------------------------------------------------
C Function:
C   Get deposition velocities for vertical diffusion
 
C Preconditions:
 
C Subroutines and functions called:
C INTERP3, M3EXIT, AERO_DEPV
 
C Revision History:
C   Jeff - Sep 97 - targetted, call subroutine for aerosol species
C   Jeff - Dec 97 - add MDEPV
C   2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   Jeff - May 2000 - correct DEPV to "real world" values
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX
C   15 Mar 03 J.Young: swap depv and mdepv dimensions
C   21 Jul 03 J.Gipson: revised to compute avg DEPV for aero3 NR CG/SOA species
C   22 Aug 03 J.Young: aero in ppm units, eliminate MDEPV
C   23 Oct 03 J.Young: change depv from flux to m/s
C   14 Sep 04 J.Gipson: generalized avg DEPV calcs for aero3 SGTOT species
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   07 Jun 05 P.Bhave: added logic for 'AE4' mechanism
C-----------------------------------------------------------------------

      USE HGRD_DEFN           ! horizontal grid specifications

      IMPLICIT NONE

      INCLUDE SUBST_GC_DEPV   ! gas chem dep vel surrogate names and map table
      INCLUDE SUBST_AE_DEPV   ! aerosol dep vel surrogate names and map table
      INCLUDE SUBST_NR_DEPV   ! non-react dep vel surrogate names and map table
      INCLUDE SUBST_TR_DEPV   ! tracer dep vel surrogate names and map table
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block (for MECHNAME)

C Arguments:

      INTEGER, INTENT( IN )  :: JDATE, JTIME ! internal simulation date&time, tstep
      LOGICAL, INTENT( IN )  :: SGSPEC( : )  ! Flag for SGTOT spec in NR_DEPV list
      REAL,    INTENT( OUT ) :: DEPV( :,:,: )

C Paramters:

C External Functions:

C Local variables:

      REAL, ALLOCATABLE, SAVE :: DBUFF( :,: )       ! read buffer

      LOGICAL, SAVE :: LAERO3_4 = .FALSE.   ! Flag for aero3 or aero4

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 16 ), SAVE :: PNAME = 'RDDEPV'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER ALLOCSTAT
      INTEGER C, R, V, N                      ! loop induction variables
      INTEGER STRT, FINI                      ! loop induction variables

      INTEGER GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2

C-----------------------------------------------------------------------

      IF ( FIRSTIME )  THEN

         FIRSTIME  =  .FALSE.

         IF ( INDEX( MECHNAME, 'AE3' ) .NE. 0 .OR.
     &        INDEX( MECHNAME, 'AE4' ) .NE. 0 ) LAERO3_4 = .TRUE.

         ALLOCATE ( DBUFF( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating BUFF'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         CALL SUBHFILE ( MET_CRO_2D, GXOFF, GYOFF,
     &                   STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2 )

      END IF   ! Firstime

      N = 0
      STRT = 1
      FINI = N_GC_DEPV
      DO V = STRT, FINI
         N = N + 1

C                                        vvvvvvv - dep vel surrogate name
         IF ( .NOT. INTERPX( MET_CRO_2D, GC_DEPV( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                       JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not read  ' // GC_DEPV( N ) // ' from ' // MET_CRO_2D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DEPV( V,C,R ) = GC_DEPV_FAC( N ) * DBUFF( C,R ) 
            END DO
         END DO

      END DO

      N = 0
      STRT = N_GC_DEPV + N_AE_DEPV + 1
      FINI = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV
      DO V = STRT, FINI
         N = N + 1

         IF ( LAERO3_4 .AND. SGSPEC( N ) ) CYCLE

C                                        vvvvvvv - dep vel surrogate name
         IF ( .NOT. INTERPX( MET_CRO_2D, NR_DEPV( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                       JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not read  ' // NR_DEPV( N ) // ' from ' // MET_CRO_2D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DEPV( V,C,R ) = NR_DEPV_FAC( N ) * DBUFF( C,R )
            END DO
         END DO

      END DO

      N = 0
      STRT = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV + 1
      FINI = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV + N_TR_DEPV 
      DO V = STRT, FINI
         N = N + 1

C                                         vvvvvvv - dep vel surrogate name
         IF ( .NOT. INTERPX( DEPV_TRAC_1, TR_DEPV( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                       JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not read  ' // TR_DEPV( N ) // ' from ' // DEPV_TRAC_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DEPV( V,C,R ) = TR_DEPV_FAC( N ) * DBUFF( C,R ) 
            END DO
         END DO

      END DO

      RETURN
      END
