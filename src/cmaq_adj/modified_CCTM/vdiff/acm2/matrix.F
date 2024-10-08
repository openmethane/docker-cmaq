
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2/matrix.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE MATRIX ( A, B, C, D, E, X )

C---------------------------------------------------------
C-- Bordered band diagonal matrix solver for ACM2

C-- ACM2 Matrix is in this form:
C   B1 E1
C   A2 B2 E2
C   A3 C3 B3 E3
C   A4    C4 B4 E4
C   A5       C5 B5 E5
C   A6          C6 B6

C--Upper Matrix is
C  U11 U12
C      U22 U23
C          U33 U34
C              U44 U45
C                  U55 U56
C                      U66

C--Lower Matrix is:
C  1
C L21  1
C L31 L32  1
C L41 L42 L43  1
C L51 L52 L53 L54  1
C L61 L62 L63 L64 L65 1
C---------------------------------------------------------

      USE VGRD_DEFN           ! vertical layer specifications

      IMPLICIT NONE

C Includes:

!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table

      INTEGER     N_SPC_DIFF               ! global diffusion species
      PARAMETER ( N_SPC_DIFF = N_GC_DIFF   ! = N_GC_SPC + 1
     &                       + N_AE_DIFF
     &                       + N_NR_DIFF
     &                       + N_TR_DIFF)

C...Arguments

!     INTEGER, INTENT(IN)   :: KL
!     INTEGER, INTENT(IN)   :: NSP
      REAL A( NLAYS )
      REAL B( NLAYS )
      REAL C( NLAYS )
      REAL E( NLAYS )
      REAL D( N_SPC_DIFF,NLAYS )
      REAL X( N_SPC_DIFF,NLAYS )

C...Locals

      REAL Y( N_SPC_DIFF,NLAYS )
      REAL L( NLAYS,NLAYS )
      REAL U( NLAYS )
      REAL UP1( NLAYS )
      REAL RU( NLAYS )
      REAL DD, DD1, YSUM
      INTEGER I, J, V

C-- Define Upper and Lower matrices

      L( 1,1 ) = 1.0
      U( 1 ) = B( 1 )
      RU( 1 ) = 1.0 / B( 1 )

      DO I = 2, NLAYS
         L( I,I ) = 1.0
         L( I,1 ) = A( I ) / B( 1 )
         UP1( I-1 ) = E( I-1 )
      END DO

      DO I = 3, NLAYS
         DO J = 2, I - 2
            DD = B( J ) - L( J,J-1 ) * E( J-1 )
            L( I,J ) = - L( I,J-1 ) * E( J-1 ) / DD
         END DO
         J = I - 1
         DD = B( J ) - L( J,J-1 ) * E( J-1 )
         L( I,J ) = ( C( I ) - L( I,J-1 ) * E( J-1 ) ) / DD
      END DO

      DO I = 2, NLAYS
         U( I ) = B( I ) - L( I,I-1 ) * E( I-1 )
         RU( I ) = 1.0 / U( I )
      END DO

C-- Forward sub for Ly=d

      DO V = 1, N_SPC_DIFF
         Y( V,1 ) = D( V,1 )
         DO I = 2, NLAYS
            YSUM = D( V,I )
            DO J = 1, I-1
               YSUM = YSUM - L( I,J ) * Y( V,J )
            END DO
            Y( V,I ) = YSUM
         END DO
      END DO

C-- Back sub for Ux=y

      DO V= 1, N_SPC_DIFF
         X( V,NLAYS ) = Y( V,NLAYS ) * RU( NLAYS )
      END DO

      DO I = NLAYS - 1, 1, -1
         DD = RU( I )
         DD1 = UP1( I )
         DO V = 1, N_SPC_DIFF
            X( V,I ) = ( Y( V,I ) - DD1 * X( V,I+1 ) ) * DD
         END DO
      END DO

      RETURN
      END

