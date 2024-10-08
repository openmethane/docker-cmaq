
C***************************************************************************
C  Significant portions of Models-3/CMAQ software were developed by        *
C  Government employees and under a United States Government contract.     *
C  Portions of the software were also based on information from non-       *
C  Federal sources, including software developed by research institutions  *
C  through jointly funded cooperative agreements. These research institu-  *
C  tions have given the Government permission to use, prepare derivative   *
C  works, and distribute copies of their work to the public within the     *
C  Models-3/CMAQ software release and to permit others to do so. EPA       *
C  therefore grants similar permissions for use of Models-3/CMAQ software, *
C  but users are requested to provide copies of derivative works to the    *
C  Government without re-strictions as to use by others.  Users are        *
C  responsible for acquiring their own copies of commercial software       *
C  associated with the Models-3/CMAQ release and are also responsible      *
C  to those vendors for complying with any of the vendors' copyright and   *
C  license restrictions. In particular users must obtain a Runtime license *
C  for Orbix from IONA Technologies for each CPU used in Models-3/CMAQ     *
C  applications.                                                           *
C                                                                          *
C  Portions of I/O API, PAVE, and the model builder are Copyrighted        *
C  1993-1997 by MCNC--North Carolina Supercomputing Center and are         *
C  used with their permissions subject to the above restrictions.          *
C***************************************************************************

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/eddy_ddm3d/s_tridiag.F,v 1.1.1.1 2010/06/14 16:03:10 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)tridiag.F	1.1 /project/mod3/CMAQ/src/vdiff/eddy/SCCS/s.tridiag.F 03 Jun 1997 11:45:09

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE STRIDIAG SUBST_GRID_ID ( L, D, U, B, X, NPMAX )
C-----------------------------------------------------------------------

C  FUNCTION:
C    Modified from TRIDIAG, to handle sensitivities. Main difference
C    is that it loops over the sens parameters.
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
C    6     Jul 2006  dsc  vert dyn alloc - use VGRD_DEFN
C    5     Jan 2002  dsc  sensitivities
C    4     Aug 96    yoj  cleaner
C    3     8/16/94   XKX  configuration management include statements
C    2     3/15/92   CJC  For use in Models-3 LCM.
C    1     10/19/89  JKV  converted for use on IBM
C    0      3/89     BDX  Initial version
C-----------------------------------------------------------------------

      USE VGRD_DEFN           ! vertical layer specifications

      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE

C Includes:

!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table
      INCLUDE SUBST_GC_SPC    ! gas chemistry species table
      INCLUDE SUBST_AE_SPC    ! aerosol species table
      INCLUDE SUBST_NR_SPC    ! non-reactive species table
      INCLUDE SUBST_TR_SPC    ! tracer species table


      INTEGER     N_SPC_DIFF               ! global diffusion species
      PARAMETER ( N_SPC_DIFF = N_GC_DIFFD  ! = N_GC_SPC + 1
     &                       + N_AE_DIFF
     &                       + N_NR_DIFF
     &                       + N_TR_DIFF)

C Arguments:


      REAL        L( NLAYS )            ! subdiagonal
      REAL        D( N_SPC_DIFF,NLAYS )      ! diagonal
      REAL        U( NLAYS )            ! superdiagonal
      REAL        B( NPMAX,N_SPC_DIFF,NLAYS )      ! R.H. side
      REAL        X( NPMAX,N_SPC_DIFF,NLAYS )      ! solution

C Local Variables:

      REAL        GAM( N_SPC_DIFF,NLAYS )
      REAL        BET( N_SPC_DIFF )
      INTEGER     V, K
      INTEGER     P, NPMAX

C Decomposition and forward substitution:

      DO V = 1, N_SPC_DIFF
         BET( V ) = 1.0 / D( V,1 )
         DO P = 1, NPMAX
           X( P,V,1 ) = BET( V ) * B( P,V,1 )
         END DO
      END DO

      DO K = 2, NLAYS
         DO V = 1, N_SPC_DIFF
            GAM( V,K ) = BET( V ) * U( K-1 )
            BET( V ) = 1.0 / ( D( V,K ) - L( K ) * GAM( V,K ) )
            DO P = 1, NPMAX
              X( P,V,K ) = 
     &           BET( V ) * ( B( P,V,K ) - L( K ) * X( P,V,K-1 ) )
            END DO
         END DO
      END DO

C Back-substitution:

      DO K = NLAYS - 1, 1, -1
         DO V = 1, N_SPC_DIFF
            DO P = 1, NPMAX
              X( P,V,K ) = X( P,V,K ) - GAM( V,K+1 ) * X( P,V,K+1 )
            END DO
         END DO
      END DO

      RETURN
      END
