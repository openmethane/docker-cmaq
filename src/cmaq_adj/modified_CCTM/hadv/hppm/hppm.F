
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/hppm/hppm.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HPPM ( NI, CONI, VEL, DT, DSI, ORI )
      
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

C   Sep 97 Jeff
C   Aug 98 - Jeff - optimize for mesh coefficients      

C   David Wong - Sep. 1998
C     -- parallelized the code
C     -- Expanded the one-level nested loop which involves either with row or
C        column, into a three-level nested loop with layers and species.
C        Corresponding arrays' dimensions were adjusted accordingly
C   Jeff - optimize for mesh coefficients
C
C   David Wong - 1/8/99
C     -- BARRIER is removed
C
C   David Wong - 1/12/99
C     -- inside BNDY_HI_PE conditional code segment, NI is changed to MY_NI
C
C   David Wong - 1/12/99
C     -- change se_loop_index argument list
C     -- add new subroutine call to determine lo and hi boundary processor

C   22 Nov 00 J.Young: PE_COMM2E -> Dave Wong's f90 stenex COMM
C                      PE_COMM3E -> Dave Wong's f90 stenex COMM

C   23 Feb 01 J.Young: allocatable arrays ...
C                      Since F90 does not preserve dummy argument array
C                      indices, CONI( 1:NI+2,, ) is copied into local array
C                      CON( 0:NI+1,, ).
C                      The caller of HPPM dimensions the actual argument,
C                      as CON( -NTHIK+1:MY_NCOLS+NTHIK,, ).

C   3 Sep 01 David Wong
C     -- use "dynamic" data structure instead of F90 ALLOCATE statement to
C        avoid memory fragmentation which eventually leads to not enough
C        contigous memory (F90 bug?)
C   24 Mar 04 G.Hammond: moved all mpi communication to caller

C   06/16/04 by Peter Percell & Daewon Byun at UH-IMAQS: 
C     - Fixed bug in using fluxes in non-uniform grids to update concentrations
C     (Note lines containing !UH)

C   14 Feb 05 J.Young: fix DS dimension bug
C----------------------------------------------------------------------
      
      USE HGRD_DEFN

      USE SUBST_MODULES              ! stenex
!     USE SUBST_COMM_MODULE          ! stenex
!     USE SUBST_UTIL_MODULE          ! stenex

      IMPLICIT NONE

C Includes:
      
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INTEGER, PARAMETER :: NTHIK = 1
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
!     INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
!     INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
!     INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
!     INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

!     INTEGER     N_SPC_ADV
!     INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV
!    &                                + N_AE_ADV
!    &                                + N_NR_ADV
!    &                                + N_TR_ADV
!    &                                + 1          ! for advecting RHO*SqRDMT
 
C Arguments:
 
      INTEGER     NI              ! number of zones (cells) along the stride
!     REAL        CON(  0:NI+1,MLAYS,N_SPC_ADV ) ! conc's along the stride
!     REAL        VEL(    NI+1,MLAYS) ! velocities at zone (cell) boundaries
!     REAL        DS ( -1:NI+1 )  ! distance between zone (cell) boundaries
      REAL     :: CONI( :,: )     ! conc's along the stride
      REAL     :: VEL( : )        ! velocities at zone (cell) boundaries
      REAL        DT              ! time step
      REAL     :: DSI ( : )       ! distance between zone (cell) boundaries
      CHARACTER   ORI             ! orientation of advection ('C'-x or 'R'-y)

C Parameters:
      
C Flag for discontinuty capturing (steepening)
#ifdef STEEPEN_FLAG
      LOGICAL, PARAMETER :: STEEPEN = .TRUE.
#else
      LOGICAL, PARAMETER :: STEEPEN = .FALSE.
#endif

      REAL, PARAMETER :: ETA1 = 20.0
      REAL, PARAMETER :: ETA2 = 0.05
      REAL, PARAMETER :: EPS = 0.01
      
      REAL, PARAMETER :: TWO3RDS = 2.0 / 3.0

C Local variables:

C orientation flag (ORI)
      CHARACTER, SAVE :: FIRSTORI = ' '
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      
!     INTEGER, PARAMETER :: NDIMS = MAX (NCOLS, NROWS)
!     PARAMETER ( NDIMS = NCOLS + NROWS )
      INTEGER, SAVE :: NSPCS
      INTEGER, SAVE :: MNR

!     REAL, ALLOCATABLE :: CON ( :,:,: ) ! conc's along the stride (local)
!     REAL, ALLOCATABLE :: DS ( : ) ! distance between cell boundaries (local)
#ifdef parallel
!*begin change by snl
      REAL :: CON ( -2:NI+3,SIZE( CONI,2 ) )
      REAL :: DS ( -3:NI+3 ) ! distance between cell boundaries (local)
!*end change by snl
#else
      REAL :: CON ( 0:NI+1,SIZE( CONI,2 ) )
      REAL :: DS ( -1:NI+1 )
#endif

      REAL A                               ! temp lattice var.
      REAL B                               ! temp lattice var.
      REAL C                               ! temp lattice var.
      REAL D                               ! temp lattice var.
!     REAL ALPHA ( NDIMS,NLAYS,N_SPC_ADV ) ! temp lattice var.
!     REAL BETA  ( NDIMS,NLAYS,N_SPC_ADV ) ! temp lattice var.
!     REAL GAMMA ( NDIMS,NLAYS,N_SPC_ADV ) ! temp lattice var.

      REAL :: GAMMA                          ! temp lattice var.
!*begin change by snl
!     REAL :: ALPHA( NI )
!     REAL :: BETA ( NI )    
      REAL, ALLOCATABLE, SAVE :: ALPHA ( : )
      REAL, ALLOCATABLE, SAVE :: BETA  ( : )
!*end change by snl
      REAL, ALLOCATABLE, SAVE :: MU    ( : ) ! lattice var. for CM
      REAL, ALLOCATABLE, SAVE :: NU    ( : ) ! lattice var. for CM
      REAL, ALLOCATABLE, SAVE :: LAMBDA( : ) ! lattice var. for CM
      REAL, ALLOCATABLE, SAVE :: CHI   ( : ) ! lattice var. for DC
      REAL, ALLOCATABLE, SAVE :: PSI   ( : ) ! lattice var. for DC
      REAL, ALLOCATABLE, SAVE :: ZETA  ( : ) ! lattice var. for ETABAR
      REAL, ALLOCATABLE, SAVE :: SIGMA ( : ) ! lattice var. for D2C
      REAL, ALLOCATABLE, SAVE :: TAU   ( : ) ! lattice var. for D2C

!     REAL FM  ( 0:NDIMS,NLAYS,N_SPC_ADV ) ! outflux from left or bottom of cell
!     REAL FP  ( 0:NDIMS,NLAYS,N_SPC_ADV ) ! outflux from right or top of cell
!     REAL CM  ( 0:NDIMS,NLAYS,N_SPC_ADV ) ! zone R.H. trial intercept
!     REAL CL    ( NDIMS,NLAYS,N_SPC_ADV ) ! zone L.H. intercept
!     REAL CR    ( NDIMS,NLAYS,N_SPC_ADV ) ! zone R.H. intercept
!     REAL DC  ( 0:NDIMS,NLAYS,N_SPC_ADV ) ! CR - CL
!     REAL C6    ( NDIMS,NLAYS,N_SPC_ADV ) ! coefficient of second-order term
!     REAL D2C   ( NDIMS,NLAYS,N_SPC_ADV ) ! second derivative
!     REAL ETA   ( NDIMS,NLAYS,N_SPC_ADV ) ! discontinuity homotopy function
!     REAL ETABAR( NDIMS,NLAYS,N_SPC_ADV ) ! 3rd to 1st order derivative ratio
!     REAL CLD   ( NDIMS,NLAYS,N_SPC_ADV ) ! zone L.H.intercept w/ discontinuity
!     REAL CRD   ( NDIMS,NLAYS,N_SPC_ADV ) ! zone R.H.intercept w/ discontinuity

C outflux from left or bottom of cell
      REAL :: FM    (    NI+1,SIZE( CONI,2 ) )
C outflux from right or top of cell
      REAL :: FP    (  0:NI,  SIZE( CONI,2 ) )
C zone R.H. trial intercept
#ifdef parallel
!*begin change by snl
      REAL :: CM    (  0:NI+2,SIZE( CONI,2 ) )
C zone L.H. intercept
      REAL :: CL    (  0:NI+1,SIZE( CONI,2 ) )
C zone R.H. intercept
      REAL :: CR    (  0:NI+1,SIZE( CONI,2 ) )
C CR - CL
      REAL :: DC    ( -1:NI+2,SIZE( CONI,2 ) )
C coefficient of second-order term
      REAL :: C6    (  0:NI+1,SIZE( CONI,2 ) )
C second derivative
      REAL :: D2C   ( -1:NI+2,SIZE( CONI,2 ) )
C discontinuity homotopy function
      REAL :: ETA   (  0:NI+1,SIZE( CONI,2 ) )
C 3rd to 1st order derivative ratio
      REAL :: ETABAR(  0:NI+1,SIZE( CONI,2 ) )
C zone L.H. intercept w/ discontinuity
      REAL :: CLD   (  0:NI+1,SIZE( CONI,2 ) )
C zone R.H. intercept w/ discontinuity
      REAL :: CRD   (  0:NI+1,SIZE( CONI,2 ) )
!*end change by snl
#else
      REAL :: CM    (    NI+1,SIZE( CONI,2 ) )
      REAL :: CL    (    NI,  SIZE( CONI,2 ) )
      REAL :: CR    (    NI,  SIZE( CONI,2 ) )
      REAL :: DC    (  0:NI+1,SIZE( CONI,2 ) )
      REAL :: C6    (    NI,  SIZE( CONI,2 ) )
      REAL :: D2C   (  0:NI+1,SIZE( CONI,2 ) )
      REAL :: ETA   (    NI,  SIZE( CONI,2 ) )
      REAL :: ETABAR(    NI,  SIZE( CONI,2 ) )
      REAL :: CLD   (    NI,  SIZE( CONI,2 ) )
      REAL :: CRD   (    NI,  SIZE( CONI,2 ) )
#endif

      CHARACTER( 96 ) :: XMSG = ' '

      REAL X                    ! Courant number
      REAL Y
      INTEGER ALLOCSTAT
      
      INTEGER I                 ! loop index
      INTEGER L                 ! layer index
      INTEGER S                 ! species index

#ifdef parallel
!*begin change by snl
      INTEGER :: START, END
!*end change by snl
#else
      INTEGER, SAVE :: MY_NI
      INTEGER, SAVE :: START1, END1
      INTEGER, SAVE :: START2, END2
      INTEGER, SAVE :: START3, END3
#endif
      LOGICAL, SAVE :: BNDY_LO_PE, BNDY_HI_PE

!     integer, save :: logdev

C----------------------------------------------------------------------

C get local DS
#ifdef parallel
!*begin change by snl
      DO I = -2, NI + 3
         DS( I ) = DSI( I+3 )
         END DO
      DS( -3 ) = DS( -2 )
!*end change by snl
#else
      DO I = -1, NI + 1
         DS( I ) = DSI( I+2 ) 
         END DO
#endif

      IF ( ORI .NE. FIRSTORI ) THEN
         FIRSTORI = ORI

         NSPCS = SIZE ( CONI,2 )
         MNR   = MAX( NCOLS, NROWS )

#ifndef parallel
!*begin change by snl
C get starting and ending loop indices
         CALL SUBST_LOOP_INDEX ( ORI, 2, NI, -1, MY_NI, START1, END1 )
         CALL SUBST_LOOP_INDEX ( ORI, 3, NI, -1, MY_NI, START2, END2 )
         CALL SUBST_LOOP_INDEX ( ORI, 3, NI, -2, MY_NI, START3, END3 )
!*end change by snl
#endif

         CALL SUBST_HI_LO_BND_PE ( ORI, BNDY_LO_PE, BNDY_HI_PE )

!        WRITE( *,* ) '   HPPM Steepen: ', STEEPEN

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
!           logdev = init3()
#ifdef parallel
!*begin change by snl
            ALLOCATE ( ALPHA ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( BETA  ( -1:MNR+2 ), STAT = ALLOCSTAT )
!*end change by snl
            ALLOCATE ( MU    ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( NU    ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( LAMBDA( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( CHI   ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( PSI   ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( ZETA  ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( SIGMA ( -1:MNR+2 ), STAT = ALLOCSTAT )
            ALLOCATE ( TAU   ( -1:MNR+2 ), STAT = ALLOCSTAT )
#else
!*begin change by snl
            ALLOCATE ( ALPHA ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( BETA  ( MNR ), STAT = ALLOCSTAT )
!*end change by snl
            ALLOCATE ( MU    ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( NU    ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( LAMBDA( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( CHI   ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( PSI   ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( ZETA  ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( SIGMA ( MNR ), STAT = ALLOCSTAT )
            ALLOCATE ( TAU   ( MNR ), STAT = ALLOCSTAT )
#endif
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating lattice variable(s)'
               CALL M3EXIT ( 'HPPM', 0, 0, XMSG, XSTAT1 )
               END IF
            END IF

#ifdef parallel
!*begin change by snl
         DO I = -1, NI + 2
!*end change by snl
#else
         DO I = START1, END1
#endif
            ALPHA( I ) = DS( I )   + DS( I+1 )
            BETA( I )  = DS( I-1 ) + DS( I )
            GAMMA      = DS( I-2 ) + DS( I-1 )
            D = DS( I ) / ( BETA( I ) + DS( I+1 ) )
            CHI( I ) = D * ( DS( I-1 ) + BETA( I ) ) / ALPHA( I )
            PSI( I ) = D * ( ALPHA( I ) + DS( I+1 ) ) / BETA( I )
            A = DS( I-1 ) / BETA( I )
            B = 2.0 * DS( I ) / BETA( I )
            C = 1.0 / ( ALPHA( I ) + GAMMA )
            MU( I ) = C * DS( I-1 ) * GAMMA / ( DS( I-1 ) + BETA( I ) )
            NU( I ) = C * DS( I ) * ALPHA( I ) / ( DS( I ) + BETA( I ) )
            LAMBDA( I ) = A + MU( I ) * B - 2.0 * NU( I ) * A
            END DO

C No need to update these arrays with off-PE data

         IF ( STEEPEN ) THEN
#ifdef parallel
!*begin change by snl
            DO I = -1, NI + 2
!*end change by snl
#else
            DO I = START1, END1
#endif
               C = 1.0 / ( BETA( I ) + DS( I+1 ) )
               SIGMA( I ) = C / ALPHA( I )
               TAU( I )   = C / BETA( I )
               ZETA( I ) = 0.25 * ( ALPHA( I ) * ALPHA( I )
     &                   -          ALPHA( I ) * BETA( I )
     &                   +          BETA( I ) * BETA( I ) )
               END DO
            END IF

         END IF   ! FIRSTORI

C get local con

#ifdef parallel
!*begin change by snl
      DO I = -2, NI + 3
         CON( I,: ) = CONI( I+3,: )
         END DO
!*end change by snl
#else
      DO I = 0, NI + 1
         CON( I,: ) = CONI( I+1,: )
         END DO
#endif

C Set all fluxes to zero. Either positive or negative flux will
C remain zero depending on the sign of the velocity.
      
      DO S = 1, NSPCS
         FP( 0,S ) = 0.0
         DO I = 1, NI
            FM( I,S ) = 0.0
            FP( I,S ) = 0.0
            END DO
         FM( NI+1,S ) = 0.0
         END DO

#ifndef parallel
!*begin change by snl
C If PE near bottom or left boundary...
C Zeroth order polynomial at the boundary cells
C First order polynomial at the next cells, no monotonicity constraint needed

      IF ( BNDY_LO_PE ) THEN
         DO S = 1, NSPCS
            CM( 1,S ) = CON( 1,S )
            CM( 2,S ) = ( DS( 1 ) * CON( 2,S )
     &                +   DS( 2 ) * CON( 1,S ) )
     &                / ( DS( 1 ) + DS( 2 ) )
            END DO
         END IF

C If PE near top or right boundary...
C Zeroth order polynomial at the boundary cells
C First order polynomial at the next cells, no monotonicity constraint needed

      IF ( BNDY_HI_PE ) THEN
         DO S = 1, NSPCS
            CM( MY_NI+1,S ) = CON( MY_NI,S )
            CM( MY_NI,S )   = ( DS( MY_NI-1 ) * CON( MY_NI,S )
     &                      +   DS( MY_NI )   * CON( MY_NI-1,S ) )
     &                      / ( DS( MY_NI-1 ) + DS( MY_NI ) )
            END DO
         END IF

C Update the remaining conc array with off-PE data, if nec.

      IF ( ORI .EQ. 'R' ) THEN
         CALL SUBST_COMM ( CON, DSPL_N1_E0_S1_W0, DRCN_N_S, 1, '1 0' )
         ELSE
         CALL SUBST_COMM ( CON, DSPL_N0_E1_S0_W1, DRCN_E_W, 1, '1 0' )
         END IF
!*end change by snl
#endif

C Second order polynomial inside the domain
      
      DO S = 1, NSPCS
#ifdef parallel
!*begin change by snl
         DO I = -1, NI + 2
!*end change by snl
#else
         DO I = START1, END1
#endif
      
C Compute average slope in the i'th zone
      
C Equation (1.7)

            DC( I,S ) = CHI( I ) * ( CON( I+1,S ) - CON( I,S ) )
     &                + PSI( I ) * ( CON( I,S )   - CON( I-1,S ) )
      
C Guarantee that CM lies between CON(I) and CON(I+1) - monotonicity constraint
      
            IF ( ( CON( I+1,S ) - CON( I,S ) )
     &         * ( CON( I,S )   - CON( I-1,S ) ) .GT. 0.0 ) THEN
               DC( I,S ) = SIGN( 1.0, DC( I,S ) )
     &                   * MIN(      ABS( DC( I,S ) ),
     &                         2.0 * ABS( CON( I+1,S ) - CON( I,S ) ),
     &                         2.0 * ABS( CON( I,S ) - CON( I-1,S ) ) )
               ELSE
               DC( I,S ) = 0.0
               END IF                                    ! Equation (1.8)

            END DO   ! I
         END DO   ! S

#ifndef parallel
!*begin change by snl
C Update the DC array with off-PE data, if nec.

      IF ( ORI .EQ. 'R' ) THEN
         CALL SUBST_COMM ( DC, DSPL_N1_E0_S1_W0, DRCN_N_S, 1, '1 0' )
         ELSE
         CALL SUBST_COMM ( DC, DSPL_N0_E1_S0_W1, DRCN_E_W, 1, '1 0' )
         END IF
!*end change by snl
#endif

      DO S = 1, NSPCS
#ifdef parallel
!*begin change by snl
         DO I = 0, NI + 2     ! Equation (1.6)
!*end change by snl
#else
         DO I = START2, END2     ! Equation (1.6)
#endif
            CM( I,S ) = CON( I-1,S )
     &                + LAMBDA( I ) * ( CON( I,S ) - CON( I-1,S ) )
     &                - MU( I ) * DC( I,S ) + NU( I ) * DC( I-1,S )
            END DO

C Initialize variables for discontinuty capturing. This is necessary
C even if discontinuity capturing is deactivated.
      
#ifdef parallel
!*begin change by snl
         DO I = 0, NI + 1         ! these reduce Eqn. 1.15 
!*end change by snl
#else
         DO I = 1, NI             ! these reduce Eqn. 1.15 
#endif
            ETA( I,S ) = 0.0
            CLD( I,S ) = CON( I,S )
            CRD( I,S ) = CON( I,S )
            END DO

         END DO   ! S

#ifdef parallel
!*begin change by snl
C Update the CM array with off-PE data, if nec.

C If PE near bottom or left boundary...
C Zeroth order polynomial at the boundary cells
C First order polynomial at the next cells, no monotonicity constraint needed

      IF ( BNDY_LO_PE ) THEN
         DO S = 1, NSPCS
            CM( 1,S ) = CON( 1,S )
            CM( 2,S ) = ( DS( 1 ) * CON( 2,S )
     &                +   DS( 2 ) * CON( 1,S ) )
     &                / ( DS( 1 ) + DS( 2 ) )
            END DO
         END IF

C If PE near top or right boundary...
C Zeroth order polynomial at the boundary cells
C First order polynomial at the next cells, no monotonicity constraint needed

      IF ( BNDY_HI_PE ) THEN
         DO S = 1, NSPCS
            CM( NI+1,S ) = CON( NI,S )
            CM( NI,S )   = ( DS( NI-1 ) * CON( NI,S )
     &                      +   DS( NI )   * CON( NI-1,S ) )
     &                      / ( DS( NI-1 ) + DS( NI ) )
            END DO
         END IF
!*end change by snl
#else
C Update the CM array with off-PE data, if nec.

      IF ( ORI .EQ. 'R' ) THEN
!        CALL SUBST_COMM ( CM, DSPL_N1_E0_S0_W0, DRCN_N, 1, '1 0' )
         CALL SUBST_COMM ( CM, DSPL_N1_E0_S0_W0, DRCN_N, 1 )
         ELSE
!        CALL SUBST_COMM ( CM, DSPL_N0_E1_S0_W0, DRCN_E, 1, '1 0' )
         CALL SUBST_COMM ( CM, DSPL_N0_E1_S0_W0, DRCN_E, 1 )
         END IF
#endif
      
      IF ( STEEPEN ) THEN
 
C Finite diff. approximation to 2nd derivative as in Equation (1.17)

         DO S = 1, NSPCS

#ifdef parallel
!*begin change by snl
            DO I = -1, NI + 2
!*end change by snl
#else
            DO I = START1, END1
#endif
               D2C( I,S ) = SIGMA( I ) * ( CON( I+1,S ) - CON( I,S ) )
     &                    -   TAU( I ) * ( CON( I,S )   - CON( I-1,S ) )
               END DO

            END DO

#ifdef parallel
!*begin change by snl
C Update the D2C array with off-PE data, if nec.

C No discontinuity detection near the boundary: cells 1, 2, NI-1, NI
 
         IF ( BNDY_LO_PE ) THEN
            START = 3
            ELSE
            START = 0
            END IF
         IF ( BNDY_HI_PE ) THEN
            END = NI - 2
            ELSE
            END = NI + 1
            END IF

         DO S = 1, NSPCS
            DO I = START, END
!*end change by snl
#else
C Update the D2C array with off-PE data, if nec.

         IF ( ORI .EQ. 'R' ) THEN
!           CALL SUBST_COMM ( D2C, DSPL_N1_E0_S1_W0, DRCN_N_S, 1 )
            CALL SUBST_COMM ( D2C, DSPL_N1_E0_S1_W0, DRCN_N_S, 1, '1 0' )
            ELSE
!           CALL SUBST_COMM ( D2C, DSPL_N0_E1_S0_W1, DRCN_E_W, 1 )
            CALL SUBST_COMM ( D2C, DSPL_N0_E1_S0_W1, DRCN_E_W, 1, '1 0' )
            END IF

C No discontinuity detection near the boundary: cells 1, 2, NI-1, NI
 
         DO S = 1, NSPCS
            DO I = START3, END3
#endif

C Compute etabars in Equation (1.16)
 
               IF ( ( - D2C( I+1,S ) * D2C( I-1,S ) .GT. 0.0 ) .AND.
     &              (
     &              ABS( CON( I+1,S ) - CON( I-1,S ) )
     &              - EPS * MIN( ABS( CON( I+1,S ) ), ABS( CON( I-1,S ) )
     &              )
     &              .GT. 0.0 ) ) THEN        ! 2nd derivative changes sign
                  ETABAR( I,S ) = - ZETA( I )
     &                          * ( D2C( I+1,S ) - D2C( I-1,S ) )
     &                          / ( CON( I+1,S ) - CON( I-1,S ) )
                  ELSE
                  ETABAR( I,S ) = 0.0
                  END IF
 
C Equation (1.16)
 
               ETA( I,S ) = MAX( 0.0,
     &                      MIN( ETA1 * ( ETABAR( I,S ) - ETA2 ), 1.0 ) ) 
 
C Equation (1.14)
 
               CRD( I,S ) = CON( I+1,S ) - 0.5 * DC( I+1,S )
               CLD( I,S ) = CON( I-1,S ) + 0.5 * DC( I-1,S )

               END DO   ! I
            END DO   ! S

         END IF                    ! if STEEPEN

#ifndef parallel
!*begin change by snl
C Update the VEL array with off-PE data, if nec.

      IF ( ORI .EQ. 'R' ) THEN
         CALL SUBST_COMM ( VEL, DSPL_N1_E0_S0_W0, DRCN_N )
         ELSE
         CALL SUBST_COMM ( VEL, DSPL_N0_E1_S0_W0, DRCN_E )
         END IF
!*end change by snl
#endif

C Generate piecewise parabolic distributions

      DO S = 1, NSPCS

#ifdef parallel
!*begin change by snl
         DO I = 0, NI + 1
!*end change by snl
#else
         DO I = 1, MY_NI
#endif

C Equation (1.15)

            CR( I,S ) = CM( I+1,S )
     &                + ETA( I,S ) * ( CRD( I,S ) - CM( I+1,S ) )
            CL( I,S ) = CM( I,S )
     &                + ETA( I,S ) * ( CLD( I,S ) - CM( I,S ) )
 
C Monotonicity
 
            IF ( ( CR( I,S ) - CON( I,S ) )
     &        * ( CON( I,S ) - CL( I,S ) ) .GT. 0.0 ) THEN

C Temporary computation of DC and C6

               DC( I,S ) = CR( I,S ) - CL( I,S )
               C6( I,S ) = 6.0 * ( CON( I,S )
     &                   - 0.5 * ( CL( I,S ) + CR( I,S ) ) )

C overshoot cases

               IF ( DC( I,S ) * C6( I,S ) .GT.
     &              DC( I,S ) * DC( I,S ) ) THEN
                  CL( I,S ) = 3.0 * CON( I,S ) - 2.0 * CR( I,S )
                  ELSE IF ( -DC( I,S ) * DC( I,S ) .GT.
     &                       DC( I,S ) * C6( I,S ) ) THEN
                  CR( I,S ) = 3.0 * CON( I,S ) - 2.0 * CL( I,S )
                  END IF

               ELSE                   ! Local extremum: Interpolation  
                                      ! function is set to be a constant
               CL( I,S ) = CON( I,S )
               CR( I,S ) = CON( I,S )

               END IF

            DC( I,S ) = CR( I,S ) - CL( I,S )      ! Equation (1.5)
            C6( I,S ) = 6.0 * ( CON( I,S )
     &                - 0.5 * ( CL( I,S ) + CR( I,S ) ) )

            END DO   ! I

C Compute fluxes from the parabolic distribution as in Equation (1.12)

#ifdef parallel
!*begin change by snl
         I = 0
         IF ( VEL( I+1 ) .GT. 0.0 ) THEN
!UH         X = VEL( I+1 ) * DT / DS( I )
            Y = VEL( I+1 ) * DT                                   !UH
            X = Y / DS( I )                                       !UH
!UH         FP( I,S ) = X * ( CR( I,S ) - 0.5 * X * ( DC( I,S )
            FP( I,S ) = Y * ( CR( I,S ) - 0.5 * X * ( DC( I,S )   !UH
     &                - C6( I,S ) * ( 1.0 - TWO3RDS * X ) ) )
            END IF
      
         DO I = 1, NI
!*end change by snl
#else
         DO I = 1, MY_NI
#endif

C mass leaving interval I at left end
C = length of segment leaving * average concentration in that segment

            IF ( VEL( I ) .LT. 0.0 ) THEN
!UH            X = -VEL( I ) * DT / DS( I )
               Y = -VEL( I ) * DT                                    !UH
               X = Y / DS( I )                                       !UH
!UH            FM( I,S ) = X * ( CL( I,S ) + 0.5 * X * ( DC( I,S )
               FM( I,S ) = Y * ( CL( I,S ) + 0.5 * X * ( DC( I,S )   !UH
     &                   + C6( I,S ) * ( 1.0 - TWO3RDS * X ) ) )
               END IF

C mass leaving interval I at right end

            IF ( VEL( I+1 ) .GT. 0.0 ) THEN
!UH            X = VEL( I+1 ) * DT / DS( I )
               Y = VEL( I+1 ) * DT                                   !UH
               X = Y / DS( I )                                       !UH
!UH            FP( I,S ) = X * ( CR( I,S ) - 0.5 * X * ( DC( I,S )
               FP( I,S ) = Y * ( CR( I,S ) - 0.5 * X * ( DC( I,S )   !UH
     &                   - C6( I,S ) * ( 1.0 - TWO3RDS * X ) ) )
               END IF

            END DO   ! I

#ifdef parallel
!*begin change by snl
         I = NI + 1
         IF ( VEL( I ) .LT. 0.0 ) THEN
!UH         X = -VEL( I ) * DT / DS( I )
            Y = -VEL( I ) * DT                                       !UH
            X = Y / DS( I )                                          !UH
!UH         FM( I,S ) = X * ( CL( I,S ) + 0.5 * X * ( DC( I,S )
            FM( I,S ) = Y * ( CL( I,S ) + 0.5 * X * ( DC( I,S )      !UH
     &                + C6( I,S ) * ( 1.0 - TWO3RDS * X ) ) )
            END IF
!*end change by snl
#endif

         END DO   ! S

C Compute fluxes from boundary cells assuming uniform distribution
      
C If PE near top or left boundary...

      IF ( BNDY_LO_PE ) THEN
         IF ( VEL( 1 ) .GT. 0.0 ) THEN
            DO S = 1, NSPCS
!UH            X = VEL( 1 ) * DT / DS( 0 )
               Y = VEL( 1 ) * DT                                     !UH
!UH            FP( 0,S ) = X * CON( 0,S )
               FP( 0,S ) = Y * CON( 0,S )                            !UH
               END DO
            END IF
         END IF

C If PE near bottom or right boundary...

      IF ( BNDY_HI_PE ) THEN
         IF ( VEL( NI+1 ) .LT. 0.0 ) THEN
            DO S = 1, NSPCS
!UH            X = -VEL( NI+1 ) * DT / DS( NI+1 )
               Y = -VEL( NI+1 ) * DT                                 !UH
!UH            FM( NI+1,S ) = X * CON( NI+1,S )
               FM( NI+1,S ) = Y * CON( NI+1,S )                      !UH
               END DO
            END IF
         END IF

#ifndef parallel
!*begin change by snl
C Update the FP and FM arrays with off-PE data, if nec.

      IF ( ORI .EQ. 'R' ) THEN
         CALL SUBST_COMM ( FP, DSPL_N0_E0_S1_W0, DRCN_S, 1, '1 0' )
!        CALL SUBST_COMM ( FM, DSPL_N1_E0_S0_W0, DRCN_N, 1, '1 0' )
         CALL SUBST_COMM ( FM, DSPL_N1_E0_S0_W0, DRCN_N, 1 )
         ELSE
         CALL SUBST_COMM ( FP, DSPL_N0_E0_S0_W1, DRCN_W, 1, '1 0' )
!        CALL SUBST_COMM ( FM, DSPL_N0_E1_S0_W0, DRCN_E, 1, '1 0' )
         CALL SUBST_COMM ( FM, DSPL_N0_E1_S0_W0, DRCN_E, 1 )
         END IF
!*end change by snl
#endif

C Update concentrations as in Equation (1.13)

      DO S = 1, NSPCS
      
#ifdef parallel
!*begin change by snl
         DO I = 1, NI
!*end change by snl
#else
         DO I = 1, MY_NI
#endif
            CON( I,S ) = CON( I,S )
!UH  &                 + FP( I-1,S ) - FP( I,S )
!UH  &                 + FM( I+1,S ) - FM( I,S )
     &                 + ( FP( I-1,S ) - FP( I,S )               !UH
     &                 +   FM( I+1,S ) - FM( I,S ) ) / DS( I )   !UH
            END DO
      
         END DO

#ifdef parallel
!*begin change by snl
      DO I = 1, NI
         CONI( I+3,: ) = CON( I,: )
         END DO
!*end change by snl
#else
      DO I = 0, NI + 1
         CONI( I+1,: ) = CON( I,: )
         END DO
#endif

      RETURN
      END
