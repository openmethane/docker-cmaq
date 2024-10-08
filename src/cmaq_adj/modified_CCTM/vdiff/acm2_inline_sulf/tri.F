
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_sulf/tri.F,v 1.1.1.1 2010/06/14 16:03:08 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)tridiag.F	1.1 /project/mod3/CMAQ/src/vdiff/eddy/SCCS/s.tridiag.F 03 Jun 1997 11:45:09

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE TRI ( L, D, U, B, X )
C-----------------------------------------------------------------------

C  FUNCTION:
C    Solves tridiagonal system by Thomas algorithm.  Algorithm fails
C    ( M3ERR ) if first pivot is zero.  In that case, rewrite the
C    equation as a set of order KMAX-1, with X(2) trivially eliminated.
C The associated tri-diagonal system is stored in 3 arrays
C   D : diagonal
C   L : sub-diagonal
C   U : super-diagonal
C   B : right hand side function
C   X : return solution from tridiagonal solver

C     [ D(1) U(1) 0    0    0 ...       0     ]
C     [ L(2) D(2) U(2) 0    0 ...       .     ]
C     [ 0    L(3) D(3) U(3) 0 ...       .     ]
C     [ .       .     .     .           .     ] X(i) = B(i)
C     [ .             .     .     .     0     ]
C     [ .                   .     .     .     ]
C     [ 0                           L(n) D(n) ]

C   where n = NLAYS

C  PRECONDITIONS REQUIRED:
C    Dimensionality set up in terms of NLAYS from SUBST_VGRD_ID

C  SUBROUTINES AND FUNCTIONS CALLED:

C  REVISION HISTORY:
C    NO.   DATE     WHO      WHAT
C    __    ____     ___      ____
C    5     Apr 06   JEP   adapted tridiag for use with ACM2
C    4     Aug 96    yoj  cleaner
C    3     8/16/94   XKX  configuration management include statements
C    2     3/15/92   CJC  For use in Models-3 LCM.
C    1     10/19/89  JKV  converted for use on IBM
C    0      3/89     BDX  Initial version
C-----------------------------------------------------------------------

      USE VGRD_DEFN           ! vertical layer specifications

      IMPLICIT NONE

C Includes:

!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table

      INTEGER     N_SPC_DIFF               ! global diffusion species
      PARAMETER ( N_SPC_DIFF = N_GC_DIFF   ! 
     &                       + N_AE_DIFF
     &                       + N_NR_DIFF
     &                       + N_TR_DIFF)

C Arguments:

      REAL        L( NLAYS )               ! subdiagonal
      REAL        D( NLAYS )               ! diagonal
      REAL        U( NLAYS )               ! superdiagonal
      REAL        B( N_SPC_DIFF,NLAYS )    ! R.H. side
      REAL        X( N_SPC_DIFF,NLAYS )    ! solution

C Local Variables:

      REAL        GAM( NLAYS )
      REAL        BET
      INTEGER     V, K

C Decomposition and forward substitution:
      BET = 1.0 / D( 1 )
      DO V = 1, N_SPC_DIFF
         X( V,1 ) = BET * B(V,1 )
      END DO

      DO K = 2, NLAYS
         GAM( K ) = BET * U( K-1 )
         BET = 1.0 / ( D( K ) - L( K ) * GAM( K ) )
	 DO V = 1, N_SPC_DIFF
            X( V, K ) = BET * ( B( V,K ) - L( K ) * X( V,K-1 ) )
	 END DO
      END DO

C Back-substitution:

      DO K = NLAYS - 1, 1, -1
         DO V = 1, N_SPC_DIFF
            X( V,K ) = X( V,K ) - GAM( K+1 ) * X( V,K+1 )
         END DO
      END DO

      RETURN
      END
