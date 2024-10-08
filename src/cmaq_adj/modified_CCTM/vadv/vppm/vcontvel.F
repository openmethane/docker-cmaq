
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vadv/vppm/vcontvel.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VCONTVEL SUBST_GRID_ID ( JDATE, JTIME, TSTEP, WHAT )
      
C-----------------------------------------------------------------------
C Function:
C    This subroutine computes the contravariant x3 velocity component
 
C Preconditions:
C    This routine can be used only for conformal map coordinates 
C    in the horizontal.
C    Dates and times should be represented YYYYDDD:HHMMSS.
 
C Subroutines and functions called:
C    INTERPX, M3EXIT, TIME2SEC, SEC2TIME, NEXTIME
      
C Revision history:
C    January 30, 1996 by Clint L. Ingram at NCSC: created for
C    RADM-coordinates

C    22 Apr 97 Jeff:
C     7 Aug 97 Jeff: for NTHIK = 1

C    Sep. 1998 David Wong
C      -- parallelized the code 

C    28 Jul 01 J.Young: Use HGRD_DEFN; replace INTERP3 with INTERPX;
C                       allocatable arrays ...
C                       Since F90 does not preserve dummy argument array
C                       indices, the 3rd dimension of WHAT has been changed
C                       from 0:NLAYS to 1:NLAYS+1
C    31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module
C-----------------------------------------------------------------------
      
      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE
      
C Includes:
      
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structuer
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
 
C Parameters

C Arguments:
      
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
!     REAL         WHAT( NCOLS,NROWS,0:NLAYS ) ! Contravariant xi-velocity 
      REAL      :: WHAT( :,:,: ) ! Contravariant xi-velocity 
      INTEGER      TSTEP        ! time step (HHMMSS)
      
C file variables:
      
      REAL         DENSJ( NCOLS,NROWS,NLAYS )  ! Jacobian * air density
 
C External functions:
 
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC
      
C local variables:
      
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      
      INTEGER      COL, ROW, LVL  ! loop indices
      INTEGER      MDATE          ! mid-time step date
      INTEGER      MTIME          ! mid-time step time
      INTEGER      STEP           ! advection time step in seconds
      REAL         DJ             ! temporary Jacobian * air density
      INTEGER      ALLOCSTAT
 
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 16 ) :: PNAME = 'VCONTVEL'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER    GXOFF, GYOFF ! global origin offset from file

C for INTERPX
      INTEGER, SAVE :: STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC
 
C-----------------------------------------------------------------------
 
      IF ( FIRSTIME ) THEN
 
         FIRSTIME = .FALSE.
 
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC )

         END IF                    ! if firstime
 
      MDATE  = JDATE
      MTIME  = JTIME
      STEP   = TIME2SEC( TSTEP )
      CALL NEXTIME( MDATE, MTIME, SEC2TIME( STEP / 2 ) )
 
C Interpolate time dependent input variables:
 
      VNAME = 'DENSA_J'
      IF ( .NOT. INTERPX ( MET_CRO_3D, VNAME, PNAME,
     &                     STRTCOLMC,ENDCOLMC, STRTROWMC,ENDROWMC, 1,NLAYS,
     &                     MDATE, MTIME, DENSJ ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

      VNAME = 'WHAT_JD'     ! x3 component of CX-vel * Jacobian * air density
      IF ( .NOT. INTERPX ( MET_CRO_3D, VNAME, PNAME,
     &                     STRTCOLMC,ENDCOLMC, STRTROWMC,ENDROWMC, 1,NLAYS,
     &                     MDATE, MTIME, WHAT ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_DOT_3D
         CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF
 
C Move 3rd dimension slabbed data from INTERP3 into proper order
C ( Using WHAT both as an argument variable and a read buffer)

      DO LVL = NLAYS, 1, -1
         DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
               WHAT( COL,ROW,LVL+1 ) = WHAT( COL,ROW,LVL ) 
               END DO
            END DO
         END DO

C Obtain flux point values of Jacobian * air density and retrieve
C contravariant velocities (Note: the recursive writes into the WHAT
C array are data-independent.)

      DO LVL = 1, NLAYS-1
         DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
               DJ = 0.5 * ( DENSJ( COL,ROW,LVL ) + DENSJ( COL,ROW,LVL+1 ) )
               WHAT( COL,ROW,LVL+1 ) = WHAT( COL,ROW,LVL+1 ) / DJ
               END DO
            END DO
         END DO

      DO ROW = 1, MY_NROWS
         DO COL = 1, MY_NCOLS
            WHAT( COL,ROW,1 ) = 0.0
            END DO
         END DO

      RETURN
      END
