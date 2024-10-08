
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vadv/vppm/vppm.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VPPM SUBST_GRID_ID ( NI, CON, VEL, DT, DS )
      
C----------------------------------------------------------------------
C Function      
C   This is the one-dimensional implementation of piecewise parabolic
C   method.  Variable grid spacing is allowed. The scheme is positive
C   definite and monotonic. It is conservative, and causes small
C   numerical diffusion.
      
C   A piecewise continuous parabola is used as the intepolation polynomial.
C   The slope of the parabola at cell edges are computed from a cumulative
C   function of the advected quantity.  These slopes are further modified
C   so that the interpolation function is monotone. For more detailed
C   information see:
      
C   Colella, P., and P. L. Woodward, (1984), "The Piecewise Parabolic
C   Method (PPM) for Gas-Dynamical Simulations," J. Comput. Phys. 54,
C   174-201.
      
C   The concentrations at boundary cells (i.e., at 1 and NI) are not
C   computed here.  They should be updated according to the boundary
C   conditions.
      
C   The following definitions are used:
     
C              |---------------> Positive direction
C     
C  -->|Boundary|<----------------Main Grid----------------->|Boundary|<--
C     
C     |---><---|---><---|       ~|---><---|~       |---><---|---><---|
C       CON(0)   CON(1)            CON(i)            CON(n)  CON(n+1)
C     
C     VEL(1)-->|        VEL(i)-->|        |-->VEL(i+1)      |-->VEL(n+1)
C    
C      FP(0)-->|       FP(i-1)-->|        |-->FP(i)         |-->FP(n)
C     
C      FM(1)<--|         FM(i)<--|        |<--FM(i+1)       |<--FM(n+1)
C    
C                             -->| DS(i)  |<--
      
C----------------------------------------------------------------------
      
C Revision History:
      
C   20 April, 1993 by M. Talat Odman at NCSC: 
C   Created based on Colella and Woodward (1984)
      
C   15 Sept., 1993 by Daewon Byun at EPA:
C   Original code obtained from Phillip Colella at Berkeley
      
C   29 Nov.,  1993 by M. Talat Odman at NCSC:
C   Found no difference from original code
      
C   05 Oct.,  1993 by M. Talat Odman at NCSC:
C   Modified for EDSS archive, made discontinuity capturing an option
      
C   Sep 97 - Jeff
C   Aug 98 - Jeff - optimize for mesh coefficients

C   06/16/04 by Peter Percell & Daewon Byun at UH-IMAQS:
C     - Fixed bug in using fluxes in non-uniform grids to update concentrations

C   07 Dec 04 J.Young: vert dyn alloc - Use VGRD_DEFN
C----------------------------------------------------------------------
      
      USE VGRD_DEFN             ! vertical layer specifications

      IMPLICIT NONE
      
C Includes:
      
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
 
C Arguments:
 
      INTEGER     NI            ! number of zones (cells)
      REAL        CON( NI )     ! concentrations in a vertical column
      REAL        VEL( NI )     ! velocities at zone (cell) boundaries
      REAL        DS ( NI )     ! distance between zone (cell) boundaries
      REAL        DT            ! time step
      
C Parameters:
      
C     Flag for discontinuty capturing (steepening)

      LOGICAL, PARAMETER :: STEEPEN = .TRUE.  ! Flag for discontinuty capturing
                                              ! (steepening)
      REAL, PARAMETER :: ETA1 = 20.0
      REAL, PARAMETER :: ETA2 = 1.0 / ETA1
      REAL, PARAMETER :: EPS = 0.01

      REAL, PARAMETER :: TWO3RDS = 2.0 / 3.0
      
C Local variables:
      
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER ALLOCSTAT
      CHARACTER( 120 ) :: XMSG = ' '

!     INTEGER, PARAMETER :: NDIMS = NLAYS + 1

      REAL A                    ! temp lattice var.
      REAL B                    ! temp lattice var.
      REAL C                    ! temp lattice var.

!     REAL ALPHA ( NDIMS )      ! temp lattice var.
!     REAL BETA  ( NDIMS )      ! temp lattice var.
!     REAL MU    ( NDIMS )      ! lattice var. for CM
!     SAVE MU

      REAL ALPHA ( NLAYS+1 )    ! temp lattice var.
      REAL BETA  ( NLAYS+1 )    ! temp lattice var.
!     SAVE MU
      REAL, ALLOCATABLE, SAVE :: MU( : )  ! lattice var. for CM

!     REAL NU    ( NDIMS )      ! lattice var. for CM
!     SAVE NU
!     REAL LAMBDA( NDIMS )      ! lattice var. for CM
!     SAVE LAMBDA
!     REAL CHI   ( NDIMS )      ! lattice var. for DC
!     SAVE CHI
!     REAL PSI   ( NDIMS )      ! lattice var. for DC
!     SAVE PSI
!     REAL ZETA  ( NDIMS )      ! lattice var. for ETABAR
!     SAVE ZETA
!     REAL SIGMA ( NDIMS )      ! lattice var. for D2C
!     SAVE SIGMA
!     REAL TAU   ( NDIMS )      ! lattice var. for D2C
!     SAVE TAU

      REAL, ALLOCATABLE, SAVE :: NU    ( : )  ! lattice var. for CM
      REAL, ALLOCATABLE, SAVE :: LAMBDA( : )  ! lattice var. for CM
      REAL, ALLOCATABLE, SAVE :: CHI   ( : )  ! lattice var. for DC
      REAL, ALLOCATABLE, SAVE :: PSI   ( : )  ! lattice var. for DC
      REAL, ALLOCATABLE, SAVE :: ZETA  ( : )  ! lattice var. for ETABAR
      REAL, ALLOCATABLE, SAVE :: SIGMA ( : )  ! lattice var. for D2C
      REAL, ALLOCATABLE, SAVE :: TAU   ( : )  ! lattice var. for D2C

!     REAL FM  ( 0:NDIMS )      ! outflux from left or bottom of cell
!     REAL FP  ( 0:NDIMS )      ! outflux from right or top of cell
!     REAL CM  ( 0:NDIMS )      ! zone R.H. trial intercept
!     REAL CL    ( NDIMS )      ! zone L.H. intercept
!     REAL CR    ( NDIMS )      ! zone R.H. intercept
!     REAL DC    ( NDIMS )      ! CR - CL
!     REAL C6    ( NDIMS )      ! coefficient of second-order term
!     REAL D2C   ( NDIMS )      ! second derivative
!     REAL ETA   ( NDIMS )      ! discontinuity homotopy function
!     REAL ETABAR( NDIMS )      ! third to first order derivative ratio
!     REAL CLD   ( NDIMS )      ! zone L.H. intercept with discontinuity
!     REAL CRD   ( NDIMS )      ! zone R.H. intercept with discontinuity

      REAL FM  ( 0:NLAYS+1 )    ! outflux from left or bottom of cell
      REAL FP  ( 0:NLAYS+1 )    ! outflux from right or top of cell
      REAL CM  ( 0:NLAYS+1 )    ! zone R.H. trial intercept
      REAL CL    ( NLAYS+1 )    ! zone L.H. intercept
      REAL CR    ( NLAYS+1 )    ! zone R.H. intercept
      REAL DC    ( NLAYS+1 )    ! CR - CL
      REAL C6    ( NLAYS+1 )    ! coefficient of second-order term
      REAL D2C   ( NLAYS+1 )    ! second derivative
      REAL ETA   ( NLAYS+1 )    ! discontinuity homotopy function
      REAL ETABAR( NLAYS+1 )    ! third to first order derivative ratio
      REAL CLD   ( NLAYS+1 )    ! zone L.H. intercept with discontinuity
      REAL CRD   ( NLAYS+1 )    ! zone R.H. intercept with discontinuity

      REAL X                    ! Courant number
      REAL Y
      
      INTEGER I                 ! loop index

      INTEGER, SAVE :: LOGDEV

C----------------------------------------------------------------------
      
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

         ALLOCATE ( MU    ( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( NU    ( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( LAMBDA( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( CHI   ( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( PSI   ( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( ZETA  ( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( SIGMA ( NLAYS+1 ), STAT = ALLOCSTAT )
         ALLOCATE ( TAU   ( NLAYS+1 ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating '
     &           // 'MU, NU, LAMBDA, CHI, PSI, ZETA, SIGMA, or TAU'
            CALL M3EXIT( 'VPPM', 0, 0, XMSG, XSTAT1 )
            END IF

         WRITE( LOGDEV,* ) '   VPPM Steepen: ', STEEPEN

         DO I = 2, NI-1
            ALPHA( I ) = DS( I ) + DS( I+1 )
            BETA( I ) = DS( I-1 ) + DS( I )
            C = DS( I ) / ( BETA( I ) + DS( I + 1 ) )
            CHI( I ) = C * ( DS( I-1 ) + BETA( I ) ) / ALPHA( I )
            PSI( I ) = C * ( ALPHA( I ) + DS( I+1 ) ) / BETA( I )
            END DO
         DO I = 2, NI-2
            A = DS( I ) / ALPHA( I )
            B = 2.0 * DS( I+1 ) / ALPHA( I )
            C = 1.0 / ( DS( I-1 ) + ALPHA( I ) + DS( I+2 ) )
            MU( I ) = C * DS( I )
     &              * ( DS( I-1 ) + DS( I ) )   / ( DS( I )   + ALPHA( I ) )
            NU( I ) = C * DS( I+1 )
     &              * ( DS( I+1 ) + DS( I+2 ) ) / ( DS( I+1 ) + ALPHA( I ) )
            LAMBDA( I ) = A + MU( I ) * B - 2.0 * NU( I ) * A
            END DO

         IF ( STEEPEN ) THEN
            DO I = 2, NI-1
               C = 1.0 / ( BETA( I ) + DS( I + 1 ) )
               SIGMA( I ) = C / ALPHA( I )
               TAU( I )   = C / BETA( I )
               ZETA( I ) = 0.25 * ( ALPHA( I ) * ALPHA( I )
     &                   -          ALPHA( I ) * BETA( I )
     &                   +          BETA( I ) * BETA( I ) )       
               END DO
            END IF

         END IF   ! FIRSTIME

C Set all fluxes to zero. Either positive or negative flux will
C remain zero depending on the sign of the velocity.
      
      DO I = 0, NI+1
         FM( I ) = 0.0
         FP( I ) = 0.0
         END DO

C Zeroth order polynomial at the boundary cells
      
      CM( 1 )    = CON( 1 )
      CM( NI+1 ) = CON( NI )
      
C First order polynomial at the next cells, no monotonicity constraint needed
      
      CM( 2 )  = ( DS( 1 ) * CON( 2 ) + DS( 2 ) * CON( 1 ) )
     &         / ( DS( 1 ) + DS( 2 ) )
      CM( NI ) = ( DS( NI-1 ) * CON( NI ) + DS( NI ) * CON( NI-1 ) )
     &         / ( DS( NI-1 ) + DS( NI ) )
      
C Second order polynomial inside the domain
      
      DO 101 I = 2, NI-1
      
C Compute average slope in the i'th zone
      
         DC( I ) = CHI( I ) * ( CON( I+1 ) - CON( I ) )
     &           + PSI( I ) * ( CON( I )   - CON( I-1 ) )      ! Equation (1.7)
      
C Guarantee that CM lies between CON(I) and CON(I+1) - monotonicity constraint
      
         IF ( ( CON( I+1 ) - CON( I ) ) * ( CON( I ) - CON( I-1 ) )
     &       .GT. 0.0 ) THEN
            DC( I ) = SIGN( 1.0, DC( I ) ) * MIN(
     &                                      ABS( DC( I ) ),
     &                                2.0 * ABS( CON( I+1 ) - CON( I ) ),
     &                                2.0 * ABS( CON( I ) - CON( I-1 ) ) )
            ELSE
            DC( I ) = 0.0
            END IF                                    ! Equation (1.8)

101      CONTINUE

      DO I = 2, NI-2                                ! Equation (1.6)
         CM( I+1 ) = CON( I ) + LAMBDA( I ) * ( CON( I+1 ) - CON( I ) )
     &             - MU( I ) * DC( I+1 ) + NU( I ) * DC( I )
         END DO
      
C Initialize variables for discontinuty capturing. This is necessary
C even if discontinuity capturing is deactivated.
      
      DO I = 1, NI         ! these reduce Eqn. 1.15 
         ETA( I ) = 0.0
         CLD( I ) = CON( I )
         CRD( I ) = CON( I )
         END DO

      IF ( STEEPEN ) THEN
 
C Finite diff. approximation to 2nd derivative as in Equation (1.17)
 
         DO I = 2, NI-1  
            D2C( I ) = SIGMA( I ) * ( CON( I+1 ) - CON( I ) )
     &               - TAU( I )   * ( CON( I )   - CON( I-1 ) )
            END DO
 
C No discontinuity detection near the boundary: cells 1, 2, NI-1, NI
 
         DO 201 I = 3, NI-2  
 
C Compute etabars in Equation (1.16)
 
            IF ( ( - D2C( I+1 ) * D2C( I-1 ) .GT. 0.0 ) .AND.
     &           (
     &           ABS( CON( I+1 ) - CON( I-1 ) ) 
     &           - EPS * MIN( ABS( CON( I+1 ) ), ABS( CON( I-1 ) )
     &           )
     &           .GT. 0.0 ) ) THEN        ! 2nd derivative changes sign
               ETABAR( I ) = - ZETA( I ) * ( D2C( I+1 ) - D2C( I-1 ) )
     &                     / ( CON( I+1 ) - CON( I-1 ) )
               ELSE
               ETABAR( I ) = 0.0
               END IF
 
C Eqn. (1.16)
 
            ETA( I ) = MAX( 0.0, MIN( ETA1 * ( ETABAR( I ) - ETA2 ), 1.0 ) ) 
 
C Equation (1.14)
 
            CRD( I ) = CON( I+1 ) - 0.5 * DC( I+1 )
            CLD( I ) = CON( I-1 ) + 0.5 * DC( I-1 )

201         CONTINUE

            END IF                    ! if STEEPEN
      
C Generate piecewise parabolic distributions
      
      DO 301 I = 1, NI

         CR( I ) = CM( I+1 )
     &           + ETA( I ) * ( CRD( I ) - CM( I+1 ) ) ! Equation (1.15)
         CL( I ) = CM( I )
     &           + ETA( I ) * ( CLD( I ) - CM( I ) )
 
C Monotonicity
 
         IF ( ( CR( I ) - CON( I ) ) * ( CON( I ) - CL( I ) ) .GT. 0.0 ) THEN

            DC( I ) = CR( I ) - CL( I )  ! Temporary computation of DC and C6
            C6( I ) = 6.0 * ( CON( I ) - 0.5 * ( CL( I ) + CR( I ) ) )

C overshoot cases

            IF ( DC( I ) * C6( I ) .GT. DC( I ) * DC( I ) ) THEN
               CL( I ) = 3.0 * CON( I ) - 2.0 * CR( I )
               ELSE IF ( -DC( I ) * DC( I ) .GT. DC( I ) * C6( I ) ) THEN
               CR( I ) = 3.0 * CON( I ) - 2.0 * CL( I )
               END IF

            ELSE                   ! Local extremum: Interpolation  
                                   ! function is set to be a constant
            CL( I ) = CON( I )
            CR( I ) = CON( I )

            END IF

         DC( I ) = CR( I ) - CL( I )      ! Equation (1.5)
         C6( I ) = 6.0 * ( CON( I ) - 0.5 * ( CL( I ) + CR( I ) ) )

301      CONTINUE
      
C Compute integrals of fluxes over the time step from the parabolic
C distribution as in Equation (1.12)
      
      DO I = 1, NI

C mass leaving interval I at bottom
C = length of segment leaving * average concentration in that segment

         IF ( VEL( I ) .LT. 0.0 ) THEN
            Y = -VEL( I ) * DT
            X = Y / DS( I )
            FM( I ) = Y * ( CL( I ) + 0.5 * X
     &              * ( DC( I ) + C6( I ) * ( 1.0 - TWO3RDS * X ) ) )
            END IF

         END DO
      
C mass leaving interval I at top

      DO I = 1, NI-1

         IF ( VEL( I+1 ) .GT. 0.0 ) THEN
            Y = VEL( I+1 ) * DT
            X = Y / DS( I )
            FP( I ) = Y * ( CR( I ) - 0.5 * X
     &              * ( DC( I ) - C6( I ) * ( 1.0 - TWO3RDS * X ) ) )
            END IF

         END DO
      
C Compute fluxes from boundary cells assuming uniform distribution
      
c      IF ( VEL( 1 ) .GT. 0.0 ) THEN   !!! Vel(1) = 0
c         X = VEL( 1 ) * DT / DS( 0 )
c         FP( 0 ) = X * CON( 0 )
c         END IF

c      IF ( VEL( NI+1 ) .LT. 0.0 ) THEN   !!! Vel(NI+1) = 0
c         X = -VEL( NI+1 ) * DT / DS( NI+1 )
c         FM( NI+1 ) = X * CON( NI+1 )
c         END IF
      
C Update concentrations as in Equation (1.13)
      
      DO I = 1, NI
         CON( I ) = CON( I )
     &            + ( FP( I-1 ) - FP( I ) + FM( I+1 ) - FM( I ) ) / DS( I )
         END DO
      
      RETURN
      END
