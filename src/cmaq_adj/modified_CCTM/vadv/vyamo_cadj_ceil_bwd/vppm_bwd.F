
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vadv/vyamo_cadj_fwd/vppm_bwd.F,v 1.1 2010/11/17 20:56:40 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VPPM_BWD ( NI, DT, DS, VEL, CON, FLX )

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
C     VEL(1)<->|        VEL(i)<->|        |<->VEL(i+1)      |<->VEL(n+1)
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
C   08 May 09 J.Young: dimension CON with species; eliminate "STEEPEN" option

C   18 Nov 09 J.Young: combine PPM velocity adjustment and vertically advected
C                      concentrations on one code
C   Aug 2010 by Peter Percell at UH-IMAQS: modified for adjoint/4dvar
C----------------------------------------------------------------------

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

C Arguments:

      INTEGER, INTENT( IN ) :: NI ! number of zones (cells)
      REAL, INTENT( IN )    :: DT ! time step
      REAL, INTENT( IN )    :: DS ( NI ) ! distance between zone (cell) boundaries
      REAL, INTENT( INOUT ) :: VEL( NI+1 ) ! velocities at zone (cell) boundaries
      REAL, INTENT( INOUT ) :: CON( NI ) ! concs in a vertical column
      REAL, INTENT( IN ), OPTIONAL :: FLX( NI+1 ) ! fluxes at zone (cell) boundaries

C Parameters:

      REAL, PARAMETER :: TWO3RDS = 2.0 / 3.0
      REAL, PARAMETER :: EPSF = 0.001

C Local variables:

      REAL    :: FM ( 1:NI+1 )  ! outflux from left or bottom of cell
      REAL    :: FP ( 0:NI )    ! outflux from right or top of cell
      REAL    :: CR ( 1:NI )    ! zone R.H. intercept
      REAL    :: CL ( 1:NI )    ! zone L.H. intercept
      REAL    :: DC ( 1:NI )    ! CR - CL
      REAL    :: C6 ( 1:NI )    ! coefficient of second-order term

      REAL X                    ! Courant number
      REAL Y                    ! removed zone slab
      REAL FDN, FUP             ! upstream donor cell versions of fm, fp
      INTEGER ICNT              ! no. of times fdn, fup exceeds fm, fp error range
      
      CHARACTER( 120 ) :: XMSG = ' '

      INTEGER I, S              ! loop index
      INTEGER, PARAMETER :: MAXCNT = 50

C----------------------------------------------------------------------

      CALL PPM ( NI, DS, CON, CR, CL, DC, C6 )

C Update velocities so that flux of RHOJ advected with the PPM algorithm matches
C the flux advected with the upstream donor cell method

      IF ( PRESENT(FLX) ) THEN
         FM( 1:NI+1 ) = 0.0
         FP( 0:NI ) = 0.0

         DO I = 1, NI
            IF ( VEL( I ) .LT. 0.0 ) THEN
               FDN = -FLX( I ) * DT
               ICNT = 0
 66            CONTINUE
               Y = -VEL( I ) * DT
               X = Y / DS( I )
               FM( I ) = Y * ( CL( I ) + 0.5 * X
     &              * ( DC( I ) + C6( I ) * ( 1.0 - TWO3RDS * X ) ) )
               IF ( ABS( FM( I ) - FDN ) .GT. EPSF * FDN ) THEN
                  VEL( I ) = VEL( I ) * SQRT( FDN / FM( I ) )
                  ICNT = ICNT + 1
                  IF ( ICNT .GT. MAXCNT ) THEN
                     XMSG = 'max iterations exceeded in vppm_bwd at 66'
                     CALL M3EXIT( 'VPPM_BWD', 0, 0, XMSG, XSTAT1 )
                  END IF
                  GO TO 66
               END IF
            END IF
            IF ( VEL( I+1 ) .GT. 0.0 ) THEN
               FUP = FLX( I+1 ) * DT
               ICNT = 0
 77            CONTINUE
               Y = VEL( I+1 ) * DT
               X = Y / DS( I )
               FP( I ) = Y * ( CR( I ) - 0.5 * X
     &              * ( DC( I ) - C6( I ) * ( 1.0 - TWO3RDS * X ) ) )
               IF ( ABS( FP( I ) - FUP ) .GT. EPSF * FUP ) THEN
                  VEL( I+1 ) = VEL( I+1 ) * SQRT( FUP / FP( I ) )
                  ICNT = ICNT + 1
                  IF ( ICNT .GT. MAXCNT ) THEN
                     XMSG = 'max iterations exceeded in vppm_bwd at 77'
                     CALL M3EXIT( 'VPPM_BWD', 0, 0, XMSG, XSTAT1 )
                  END IF
                  GO TO 77
               END IF
            END IF
         END DO
      END IF                    ! PRESENT(FLX)

c set all fluxes to zero. either positive or negative flux will
c remain zero depending on the sign of the velocity.

c fm: function for mass leaving interval i at lower face (i-1/2)
c = length of segment leaving * integral average concentration in that segment:
c   length of segment leaving = y = -v(i)dt
c   segment integral ave. conc. = cl(i) + 1/2[c6(i)+dc(i)]dx -1/3[c6(i)]dx**2
c fp: function for mass leaving interval i at upper face (i+1/2)
c = length of segment leaving * integral average concentration in that segment:
c   length of segment leaving = y = v(i+1)dt
c   segment integral ave. conc. = cr(i) + 1/2[c6(i)-dc(i)]dx -1/3[c6(i)]dx**2

      FM( 1:NI+1 ) = 0.0
      FP( 0:NI ) = 0.0

C Update concentrations of species

      DO I = 1, NI

         IF ( VEL( I ) .LT. 0.0 ) THEN
            Y = -VEL( I ) * DT
            X = Y / DS( I )
            FM( I ) = Y * ( CL( I ) + 0.5 * X
     &           * ( DC( I ) + C6( I ) * ( 1.0 - TWO3RDS * X ) ) )
         END IF

         IF ( VEL( I+1 ) .GT. 0.0 ) THEN
            Y = VEL( I+1 ) * DT
            X = Y / DS( I )
            FP( I ) = Y * ( CR( I ) - 0.5 * X
     &           * ( DC( I ) - C6( I ) * ( 1.0 - TWO3RDS * X ) ) )
         END IF

      END DO

c compute fluxes for top face

      I = NI+1
      IF ( VEL( I ) .LT. 0.0 ) THEN
         IF ( PRESENT(FLX) ) THEN
            Y = -VEL( I ) * DT
            FM( I ) = Y * CON( I-1 )
         ELSE
            FM( I ) = 0.
         END IF
      END IF

      DO I = 1, NI
         CON( I ) = CON( I )
     &        + ( FP( I-1 ) - FP( I ) + FM( I+1 ) - FM( I ) ) / DS( I )
      END DO

      RETURN
      END SUBROUTINE VPPM_BWD

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PPM ( NI, DS, CN, CR, CL, DC, C6 )

C get ppm coefficients CR, CL, DC, and C6

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

C arguments:

      INTEGER, INTENT(    IN ) :: NI ! number of zones (cells)
      REAL,    INTENT(    IN ) :: DS( NI ) ! distance between zone (cell) boundaries
      REAL,    INTENT(    IN ) :: CN( NI ) ! concentrations in a vertical column
      REAL,    INTENT( INOUT ) :: CR( NI ) ! zone r.h. intercept
      REAL,    INTENT( INOUT ) :: CL( NI ) ! zone l.h. intercept
      REAL,    INTENT( INOUT ) :: DC( NI ) ! CR - CL
      REAL,    INTENT( INOUT ) :: C6( NI ) ! coefficient of second-order term

C local variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER ALLOCSTAT
      CHARACTER( 120 ) :: XMSG = ' '

      REAL A, B, C, BETA        ! temp lattice vars.

      REAL, ALLOCATABLE, SAVE :: ALPHA ( : ) ! temp lattice var.
      REAL, ALLOCATABLE, SAVE :: CHI   ( : ) ! lattice var. for dc
      REAL, ALLOCATABLE, SAVE :: PSI   ( : ) ! lattice var. for dc
      REAL, ALLOCATABLE, SAVE :: MU    ( : ) ! lattice var. for cm
      REAL, ALLOCATABLE, SAVE :: NU    ( : ) ! lattice var. for cm
      REAL, ALLOCATABLE, SAVE :: LAMBDA( : ) ! lattice var. for cm

      REAL :: CM ( 1:NI+1 )                  ! zone r.h. trial intercept

      INTEGER I

C----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         ALLOCATE ( ALPHA ( 2:NI-1 ),
     &        CHI   ( 2:NI-1 ),
     &        PSI   ( 2:NI-1 ),
     &        MU    ( 2:NI-2 ),
     &        NU    ( 2:NI-2 ),
     &        LAMBDA( 2:NI-2 ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating ALPHA, CHI, PSI, MU, NU, LAMBDA'
            CALL M3EXIT( 'VPPM_BWD', 0, 0, XMSG, XSTAT1 )
         END IF

         DO I = 2, NI-1
            ALPHA( I ) = DS( I ) + DS( I+1 )
            BETA = DS( I-1 ) + DS( I )
            C = DS( I ) / ( BETA + DS( I+1 ) )
            CHI( I ) = C * ( DS( I-1 ) + BETA ) / ALPHA( I )
            PSI( I ) = C * ( ALPHA( I ) + DS( I+1 ) ) / BETA
         END DO
         DO I = 2, NI-2
            A = DS( I ) / ALPHA( I )
            B = 2.0 * DS( I+1 ) / ALPHA( I )
            C = 1.0 / ( DS( I-1 ) + ALPHA( I ) + DS( I+2 ) )
            MU( I ) = C * DS( I )
     &           * ( DS( I-1 ) + DS( I ) )   / ( DS( I )   + ALPHA( I ) )
            NU( I ) = C * DS( I+1 )
     &           * ( DS( I+1 ) + DS( I+2 ) ) / ( DS( I+1 ) + ALPHA( I ) )
            LAMBDA( I ) = A + MU( I ) * B - 2.0 * NU( I ) * A
         END DO
      END IF                    ! FIRSTIME

C zeroth order polynomial at the boundary cells
C first order polynomial at the next cells, no monotonicity constraint needed

      CM( 1 )    = CN( 1 )
      CM( 2 )    = ( DS( 1 ) * CN( 2 ) + DS( 2 ) * CN( 1 ) )
     &     / ( DS( 1 ) + DS( 2 ) )
      CM( NI+1 ) = CN( NI )
      CM( NI )   = ( DS( NI-1 ) * CN( NI ) + DS( NI ) * CN( NI-1 ) )
     &     / ( DS( NI-1 ) + DS( NI ) )

C second order polynomial inside the domain

      DO I = 2, NI-1

C compute average slope in zone i

         DC( I ) = CHI( I ) * ( CN( I+1 ) - CN( I ) )
     &        + PSI( I ) * ( CN( I )   - CN( I-1 ) ) ! equation (1.7)

C guarantee that cm lies between con(i) and con(i+1) - monotonicity constraint

         IF ( ( CN( I+1 ) - CN( I ) ) * ( CN( I ) - CN( I-1 ) ) .GT. 0.0 ) THEN
            DC( I ) = SIGN( 1.0, DC( I ) ) * MIN(
     &           ABS( DC( I ) ),
     &           2.0 * ABS( CN( I+1 ) - CN( I ) ),
     &           2.0 * ABS( CN( I ) - CN( I-1 ) ) )
         ELSE
            DC( I ) = 0.0
         END IF                 ! equation (1.8)

      END DO

      DO I = 2, NI-2            ! equation (1.6)
         CM( I+1 ) = CN( I ) + LAMBDA( I ) * ( CN( I+1 ) - CN( I ) )
     &        - MU( I ) * DC( I+1 ) + NU( I ) * DC( I )
      END DO

C generate piecewise parabolic distributions

      DO I = 1, NI

         CR( I ) = CM( I+1 )    ! equation (1.15)
         CL( I ) = CM( I )

C monotonicity

         IF ( ( CR( I ) - CN( I ) ) * ( CN( I ) - CL( I ) ) .GT. 0.0 ) THEN

            DC( I ) = CR( I ) - CL( I ) ! temporary computation of dc and c6
            C6( I ) = 6.0 * ( CN( I ) - 0.5 * ( CL( I ) + CR( I ) ) )

C overshoot cases

            IF ( DC( I ) * C6( I ) .GT. DC( I ) * DC( I ) ) THEN
               CL( I ) = 3.0 * CN( I ) - 2.0 * CR( I )
            ELSE IF ( -DC( I ) * DC( I ) .GT. DC( I ) * C6( I ) ) THEN
               CR( I ) = 3.0 * CN( I ) - 2.0 * CL( I )
            END IF

         ELSE                   ! local extremum: interpolation function is set to be a constant
            CL( I ) = CN( I )
            CR( I ) = CN( I )

         END IF

         DC( I ) = CR( I ) - CL( I ) ! EQUATION (1.5)
         C6( I ) = 6.0 * ( CN( I ) - 0.5 * ( CL( I ) + CR( I ) ) )

      END DO

      RETURN
      END SUBROUTINE PPM
