
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/hppm/yadvppm.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE YADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )

C-----------------------------------------------------------------------
C Function:
C   Advection in the horizontal plane; x2-direction:
C   The process time step is set equal to TSTEP(2). Boundary concentrations
C   are coupled in RDBCON with SqRDMT = Sq. Root [det ( metric tensor )]
C   = Jacobian / (map scale factor)**2
C   where Air Density X SqRDMT is loaded into last BCON slot for advection.
      
C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.  All boundary input variables have the
C   same boundary perimeter structure with a thickness of 1
C   CGRID in transport units: SQRT{DET[metric tensor]}*concentration (Mass/Vol)
      
C Subroutines and functions called:
C   INTERP3, M3EXIT, TIME2SEC, ZFDBC, HPPM
 
C Revision history:
C   02/23/93 by M. Talat Odman  at NCSC: created
C   01/17/94 by Dongming Hwang at NCSC:   adds species mapping array A2C
C   09/08/94 by Dongming Hwang at NCSC:  ports from sequoia into
C   archive directory and putin CM template
      
C   11/08/94 by M. Talat Odman at NCSC:
C   1. Generalized module for a boundary of thickness NTHIK
C   2. Still uses the zero slope boundary condition at the outflow
C   boundaries.  When NTHIK is larger than one, we do not know how 
C   to treat the outflow boundaries, yet.
      
C   15 Dec. 1994,  by M. Talat Odman at NCSC:
C     1. Now, the boundary conditions are read in at the beginning 
C        of the time step (JDATE, JTIME) because the Bott scheme
C        is explicit. Previously, they were being read in at the
C        middle of the time step (MDATE, MTIME).  This proved to
C        create problems during testing of advection in presence
C        of chemistry.
 
C   26 June 1995, by M. Talat Odman at NCSC:
C   The Nested Grid Version
       
C   7 August 1995,  by Clint L. Ingram at NCSC:
C   The Generalized Coordinate Version

C   9 July 1996, by M. Talat Odman at NCSC:
C   Zero Slope Boundary Condition (ZSBC) of Pleim (1991) does not work with
C   NTHIK > 1, when there is divergence in the first cell near boundary.
C   When the boundary velocity is outward but very small, this condition
C   may lead to a very large boundary (outflow) concentration.  Since we
C   use the 4th order area-preserving polynomials to describe the subgrid
C   distribution of concentrations near the boundary when NTHIK > 1, this
C   eventually may affect inter-cell fluxes. So, ZSBC was taken out.

C   Jeff - Dec, 1966 - Aug, 1977
C   Put ZSBC (really Pleim's no flux divergence outflow - ZFDBC) back in
C   Gen. coord.

C   David Wong - Sep. 1998
C     -- parallelized the code
C     -- Converted the four-level nested loop which calls HPPM into one-level
C        loop with respect to column dimenion. Corresponding arrays' dimensions
C        were adjusted accordingly

C    1/19/99 David Wong
C      -- add a loop_index call
C      -- change loop index ending point to avoid accessing invalid region.
C         (reason to do this is to prevent using boundary data from PINTERP,
C          which sets pseudo-boundary data to 0)

C    3/18/99 David Wong
C      -- replaced NROWS with MY_NROWS is the DS calculation loop

C    15 Dec 00 J.Young: move CGRID_MAP into f90 module
C                       Dave Wong's f90 stenex, PE_COMM1 -> COMM

C    23 Feb 01 J.Young: allocatable arrays

C    03 Sep 01 David Wong: changed back from allocatable arrays
C    20 Mar 02 J.Young: bcon for NTHIK = 1 only
C     3 Apr 02 J.Young: remove XSTAT.EXT - now in PARMS3.EXT
C     5 Apr 02 J.Young: change copyright header
C    23 Jun 03 J.Young: for layer dependent advection tstep
C    25 Mar 04 G.Hammond: moved all communication out of HPPM to this level;
C                         using "swap_sandia" communication; update only local
C                         values in the CGRID array within a time step,
C                         discarding previous ghost values.
C    31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module
C-----------------------------------------------------------------------
      
      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE SUBST_MODULES         ! stenex
!     USE SUBST_COMM_MODULE     ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex

#ifdef parallel
!*begin added by snl
      USE SWAP_SANDIA
!*end added by snl
#endif

      IMPLICIT NONE
      
C Includes:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
!     INCLUDE SUBST_AE_SPC      ! aerosol species table
!     INCLUDE SUBST_NR_SPC      ! non-reactive species table
!     INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID    ! coordinate & domain definitions (req IOPARMS)
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Arguments:
      
!     REAL        CGRID( NCOLS,NROWS,NLAYS,* ) ! concentrations
!     REAL     :: CGRID( :,:,:,: )
      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER     JDATE         ! current model date, coded YYYYDDD
      INTEGER     JTIME         ! current model time, coded HHMMSS
      INTEGER     TSTEP( 2 )    ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
!     INTEGER     ASTEP( NLAYS )  ! layer advection time step
      INTEGER     ASTEP( : )    ! layer advection time step

C External Functions not declared in IODECL3.EXT:
      
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC
      REAL,    EXTERNAL :: ZFDBC
      
C Parameters:

C Advected species dimension

      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV
     &                                + N_AE_ADV
     &                                + N_NR_ADV
     &                                + N_TR_ADV
     &                                + 1        ! for advecting RHO*SqRDMT
 
C File Variables:
 
!     REAL         VHAT( NCOLS+1,NROWS+1,NLAYS ) ! x2-component CX-velocity
!     REAL         BCON( NBNDY,NLAYS,N_SPC_ADV ) ! boundary concentrations
      REAL         VHAT( NCOLS+1,NROWS+1 )       ! x2-component CX-velocity
      REAL         BCON( NBNDY,N_SPC_ADV )       ! boundary concentrations
!     REAL, ALLOCATABLE :: VHAT( :,:,: )         ! x2-component CX-velocity
!     REAL, ALLOCATABLE :: BCON( :,:,: )         ! boundary concentrations

C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'YADVPPM'
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 16 ) :: UORV = 'X2VEL'          ! flag for HCONTVEL
      CHARACTER( 96 ) :: XMSG = ' '

      REAL, SAVE :: DX2                           ! dx1 (meters)

! Concentrations along a column
#ifdef parallel
!*begin change by snl
      REAL         CON( -2:NROWS+3,N_SPC_ADV )
      REAL         DS(  -2:NROWS+3 )              ! ds = DX2
!*end change by snl
#else
      REAL         CON(  0:NROWS+1,N_SPC_ADV )
      REAL         DS(  -1:NROWS+1 )
#endif
      REAL         VEL( NROWS+1 )        ! Velocities along a column
!     REAL, ALLOCATABLE :: CON( :,:,: )
!     REAL, ALLOCATABLE :: VEL( :,: )
!     REAL, ALLOCATABLE :: DS( : )               ! ds = DX2
      REAL         DT                            ! dt
      INTEGER      STEP                          ! FIX dt
      INTEGER      DSTEP                         ! dt accumulator
      INTEGER      FDATE                         ! interpolation date
      INTEGER      FTIME                         ! interpolation time
      INTEGER      ALLOCSTAT

      INTEGER, SAVE :: ADV_MAP( N_SPC_ADV )      ! global adv map to CGRID

      INTEGER      COL, ROW, LVL, SPC, VAR       ! loop counters
      INTEGER      A2C

      INTEGER MY_TEMP
      INTEGER, SAVE :: STARTROW, ENDROW
      LOGICAL, SAVE :: BNDY_LO_PE, BNDY_HI_PE

C Statement functions:
      
#ifdef parallel
!*begin added by snl
      integer, parameter :: swap_width = 3
      real bound_north( ncols,swap_width,n_spc_adv )
      real bound_south( ncols,swap_width,n_spc_adv )
      real buf_ns( swap_width*ncols*n_spc_adv )
      integer north_row, m
!*end added by snl
#else
      integer, parameter :: swap_width = 1
#endif

      INTEGER, SAVE :: SFX    ! fixed parameter for southern boundary
      INTEGER, SAVE :: NFX    ! fixed parameter for northern boundary

      REAL    BCCN            ! boundary concentrations

      INTEGER BFX             ! dummy positional parameter
      INTEGER CC              ! column index
      INTEGER SS              ! species index

      BCCN ( BFX, CC, SS ) = BCON( BFX + CC, SS )

C Required interface for allocatable array dummy arguments

      INTERFACE
         SUBROUTINE RDBCON ( FDATE, FTIME, TSTEP, LVL, BCON )
            IMPLICIT NONE
            INTEGER,         INTENT( IN )     :: FDATE, FTIME, TSTEP, LVL
            REAL,            INTENT( OUT )    :: BCON( :,: )
         END SUBROUTINE RDBCON
         SUBROUTINE HCONTVEL ( FDATE, FTIME, TSTEP, LVL, UORV, UHAT )
            IMPLICIT NONE
            INTEGER,         INTENT( IN )     :: FDATE, FTIME, TSTEP, LVL
            CHARACTER( 16 ), INTENT( IN )     :: UORV
            REAL,            INTENT( OUT )    :: UHAT( :,: )
         END SUBROUTINE HCONTVEL
         SUBROUTINE HPPM ( NI, CON, VEL, DT, DS, ORI )
            IMPLICIT NONE
            INTEGER,         INTENT( IN )     :: NI
            REAL,            INTENT( IN OUT ) :: CON( :,: )
            REAL,            INTENT( IN )     :: VEL( : )
            REAL,            INTENT( IN )     :: DT
            REAL,            INTENT( IN )     :: DS ( : )
            CHARACTER,       INTENT( IN )     :: ORI
         END SUBROUTINE HPPM

      END INTERFACE
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         SFX = 0
         NFX = MY_NCOLS + MY_NROWS + 3

C Get dx2 from COORD include file

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX2 = DG2M * YCELL_GD   ! in m.
            ELSE
            DX2 = YCELL_GD          ! in m.
            END IF

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
 
C Create global map to CGRID
 
         SPC = 0
         DO VAR = 1, N_GC_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = GC_STRT - 1 + GC_ADV_MAP( VAR )
            END DO
         DO VAR = 1, N_AE_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = AE_STRT - 1 + AE_ADV_MAP( VAR )
            END DO
         DO VAR = 1, N_NR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = NR_STRT - 1 + NR_ADV_MAP( VAR )
            END DO
         DO VAR = 1, N_TR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = TR_STRT - 1 + TR_ADV_MAP( VAR )
            END DO
 
         ADV_MAP( N_SPC_ADV ) = N_GC_SPCD

         CALL SUBST_LOOP_INDEX ( 'R', 1, MY_NROWS, 1, MY_TEMP,
     &                           STARTROW, ENDROW )

         CALL SUBST_HI_LO_BND_PE ( 'R', BNDY_LO_PE, BNDY_HI_PE )

         END IF                    ! if firstime

#ifdef parallel
!*begin added by snl
!     ALLOCATE ( DS( -2:MY_NROWS+3 ), STAT = ALLOCSTAT )
      DO ROW = -2, MY_NROWS + 3
         DS ( ROW ) = DX2
         END DO
      call swap2d( ds( 1 ),ds( my_nrows+1 ), swap_width, 1, 1, NORTH )
      call swap2d( ds( my_nrows-2 ), ds( -2 ), swap_width, 1, 1, SOUTH )
!*end added by snl
#else
      DO ROW = -1, MY_NROWS + 1
         DS ( ROW ) = DX2
         END DO

      CALL SUBST_COMM ( DS, DSPL_N1_E0_S2_W0, DRCN_N_S, '1 -1' )
#endif

      DO 301 LVL = 1, NLAYS

         STEP = TIME2SEC( ASTEP( LVL ) )
         DT = FLOAT( STEP )
         DSTEP = STEP
         FDATE = JDATE
         FTIME = JTIME

101      CONTINUE

C Get the boundary concentration values (converted to transport units)

         CALL RDBCON ( FDATE, FTIME, ASTEP( LVL ), LVL, BCON )

C Get the contravariant x2 velocity component

         CALL HCONTVEL ( FDATE, FTIME, ASTEP( LVL ), LVL, UORV, VHAT )

#ifdef parallel
!*begin added by snl
         call swap2d( vhat( 1,1 ), vhat( 1,nrows+1 ), ncols+1, 1, ncols+1, NORTH )

         m = 0
         north_row = nrows - swap_width
         do spc = 1, n_spc_adv
            a2c = adv_map( spc )
            do row = 1, swap_width
               do col = 1, ncols
                  bound_south( col,row,spc ) = cgrid( col,row,lvl,a2c )
                  bound_north( col,row,spc ) = cgrid( col,north_row+row,lvl,a2c )
                  m = m + 1
                  buf_ns(m) = bound_north(col,row,spc)
               enddo
            enddo
         enddo

! swap3d interface ...
!     swap3d( 1st value of array to be sent,
!             1st value of array to be received,
!             number of values to send in 1st dimension of array,
!             number of values to send in 2nd dimension of array,
!             number of values to send in 3rd dimension of array,
!             size of 1st dimension in arrays,  
!             size of 2nd dimension in arrays,
!             direction to receive from: NSEW)

         call swap3d( bound_south, bound_north, ncols, swap_width, n_spc_adv,
     &                ncols, swap_width, NORTH )
         call swap3d( buf_ns, bound_south, ncols, swap_width, n_spc_adv,
     &                ncols, swap_width, SOUTH )
!*end added by snl
#endif

C Do the gridded computation for horizontal advection

         DO 233 COL = 1, MY_NCOLS    ! x2-direction

            DO ROW = 1, MY_NROWS+1
               VEL( ROW ) = VHAT( COL,ROW )
               END DO

            DO SPC = 1, N_SPC_ADV

               A2C = ADV_MAP( SPC )
               DO ROW = 1, MY_NROWS
                  CON( ROW,SPC ) = CGRID( COL,ROW,LVL,A2C )
                  END DO

#ifdef parallel
!*begin added by snl
               do row = 1, swap_width
                  con( row-swap_width,spc ) = bound_south( col,row,spc )
                  con( nrows+row,spc ) = bound_north( col,row,spc )
               enddo
!*end added by snl
#endif

C South boundary

               IF ( BNDY_LO_PE ) THEN
                  IF ( VEL( 1 ) .LT. 0.0 ) THEN          ! outflow
                     CON( 1-swap_width:0,SPC) = 
     &                              ZFDBC ( CON( 1,SPC ), CON( 2,SPC ),
     &                                      VEL( 1 ), VEL( 2 ) )
                     ELSE    ! inflow
                     CON( 1-swap_width:0,SPC ) = BCCN( SFX,COL,SPC )
                     END IF
                  END IF

C North boundary

               IF ( BNDY_HI_PE ) THEN
                  IF ( VEL( MY_NROWS+1 ) .GT. 0.0) THEN     ! outflow
                     CON( MY_NROWS+1:MY_NROWS+swap_width,SPC ) = 
     &                                       ZFDBC ( CON( MY_NROWS,SPC ), 
     &                                               CON( MY_NROWS-1,SPC ),
     &                                               VEL( MY_NROWS+1 ), 
     &                                               VEL( MY_NROWS ) )
                     ELSE    ! inflow
                     CON( MY_NROWS+1:MY_NROWS+swap_width,SPC ) = 
     &                                                   BCCN( NFX,COL,SPC )
                     END IF
                  END IF

               END DO

C PPM scheme

            CALL HPPM ( NROWS, CON, VEL, DT, DS, 'R' )

            DO SPC = 1, N_SPC_ADV
               A2C = ADV_MAP(SPC)
               DO ROW = 1, MY_NROWS
                  CGRID( COL,ROW,LVL,A2C ) = CON( ROW,SPC )
                  END DO
               END DO

233         CONTINUE

         DSTEP = DSTEP + STEP
         IF ( DSTEP .LE. TIME2SEC( TSTEP( 2 ) ) ) THEN
            CALL NEXTIME( FDATE, FTIME, SEC2TIME( STEP ) )
            GO TO 101
            END IF

301      CONTINUE

      RETURN
      END
