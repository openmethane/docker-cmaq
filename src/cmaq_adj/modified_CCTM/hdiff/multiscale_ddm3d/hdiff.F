
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hdiff/multiscale_ddm3d/hdiff.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef sens
      SUBROUTINE HDIFF ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
#else
      SUBROUTINE HDIFF ( CGRID, JDATE, JTIME, TSTEP )
#endif

C-----------------------------------------------------------------------
C Function:
C   Horizontal diffusion with constant eddy diffusivity - gen. coord.
C   The process time step is set equal to TSTEP(2). Boundary concentrations
C   are set using a Dirichlet (no flux) condition
      
C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.  All boundary input variables have the
C   same boundary perimeter structure with a thickness of 1
C   CGRID in ppm units or micro-g/m**3, #/m**3 for aerosols
      
C Subroutines and functions called:
C   INIT3, TIME2SEC, SEC2TIME, CGRID_MAP, NEXTIME, RHO_J, 
C   HCDIFF3D
 
C Revision history:
C   Jeff - 5 Nov 97, 1 Jan 98
C   DWB  - 1 Feb 98, use simple B/C (no conc gradient at domain boundary)

C   David Wong Sep. 1998
C     -- parallelized the code
C     -- removed the intermediate constant CRHOJ_Q and placed the answer of
C        the calculation directly into CGRID. Removed the next immediate
C        loop completely.

C   David Wong 1/19/99
C      -- add a loop_index call
C      -- change loop index ending point to avoid accessing invalid region.
C         (reason to do this is to prevent using boundary data from PINTERP,
C          which sets pseudo-boundary data to 0)
 
C   Daewon Byun 10/10/2000
C      -- generalized 3d horizontal diffusivity module
C      -- accomdates 3d hdiff values

C    15 Dec 00 J.Young: PE_COMM3 -> Dave Wong's f90 stenex COMM 
C     6 Aug 01 J.Young: Use HGRD_DEFN
C    25 Mar 04 G.Hammond: RK11/RK22 ghost cell updates moved outside main loop;
C                         use explicit boundary arrays for CGRID ghost cells;
C                         use SNL's "swap3d".
C    31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module
C
C    Jul 06 D.Cohan: incorporated DDM sensitivity, but NOT for parallel
C    Sep 06 D.Cohan and S.Tonse: incorporated parallel capability for DDM sens
C    Nov 08 S.Napelenok: DDM-3D implementation for CMAQv4.7
C    17 Dec 09 J.Young: fix CONC initializing error for sub-cycling timesteps (DO 344)
C                       reported by Talat Odman and Yongtao Hu of GA tech.
C     9 Feb 10 S.Napelenok: added DDM-3D to this version 
C-----------------------------------------------------------------------
      
      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex
!     USE SUBST_COMM_MODULE     ! stenex
      
#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_DIFF     ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF     ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF     ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF     ! tracer diffusion species and map table
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Arguments:
      
      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER     JDATE         ! current model date, coded YYYYDDD
      INTEGER     JTIME         ! current model time, coded HHMMSS
      INTEGER     TSTEP( 2 )    ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C External Functions not declared in IODECL3.EXT:
      
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, GETEFILE
      
C Parameters:

C Advected species dimension

      INTEGER, PARAMETER :: N_SPC_DIFF = N_GC_DIFF
     &                                 + N_AE_DIFF
     &                                 + N_NR_DIFF
     &                                 + N_TR_DIFF
!    &                                 + 1       ! diff RHO_J
 
C File Variables:
 
      REAL          CONC  ( 0:NCOLS+1,0:NROWS+1 )   ! conc working array
      REAL          RHOJ  ( 0:NCOLS+1,0:NROWS+1,NLAYS ) ! density X Jacobian
      CHARACTER( 8 ), SAVE :: COMMSTR                 ! for both CONC and RHOJ
      REAL          RK11  (   NCOLS+1,NROWS+1,NLAYS ) ! initially used as RHOJ
                              ! at x1 cell face, then finally as 11 eddy diff. factor
      REAL          RK22  (   NCOLS+1,NROWS+1,NLAYS ) ! initially used as RHOJ
                              ! at x2 cell face, then finally as 22 eddy diff. factor
      REAL          K11BAR3D ( NCOLS+1,NROWS+1,NLAYS ) ! ave. Cx11 eddy diff
      REAL          K22BAR3D ( NCOLS+1,NROWS+1,NLAYS ) ! ave. Cx22 eddy diff
      REAL          DT                          ! diffusion time step
      REAL          CRHOJ_Q                     ! intermediate, coupled conc.

C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'HDIFF'
      
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      REAL          DX1                         ! dx1 (meters)
      REAL          DX2                         ! dx2 (meters)
      REAL, SAVE :: RDX1S                       ! reciprocal dx1*dx1
      REAL, SAVE :: RDX2S                       ! reciprocal dx2*dx2
      
      REAL          DTDX1S                      ! dt/dx1**2
      REAL          DTDX2S                      ! dt/dx2**2
      REAL          DTSEC                       ! model time step in seconds
      INTEGER       NSTEPS                      ! diffusion time steps
      INTEGER       STEP                        ! FIX dt
      INTEGER       FDATE                       ! interpolation date
      INTEGER       FTIME                       ! interpolation time

#ifdef sens
      REAL, POINTER :: SENGRID( :,:,:,:,: )      ! sensitivities (COLS,ROWS,LAYS,P,SPC)
      REAL         SENSWORK  ( 0:NCOLS+1,0:NROWS+1 )   ! sens working array
      INTEGER      P                           ! loop counter for sensitivity
      INTEGER      NPMAX                       ! # sensitivity parameters
#endif


#ifdef parallel
      INTEGER      M
      REAL         HALO_SOUTH( NCOLS,NLAYS,N_SPC_DIFF )
      REAL         HALO_NORTH( NCOLS,NLAYS,N_SPC_DIFF )
      REAL         HALO_WEST ( NROWS,NLAYS,N_SPC_DIFF )
      REAL         HALO_EAST ( NROWS,NLAYS,N_SPC_DIFF )
      REAL         BUF_NS    ( NCOLS*NLAYS*N_SPC_DIFF )
      REAL         BUF_EW    ( NROWS*NLAYS*N_SPC_DIFF )

#ifdef sens
      REAL         S_HALO_SOUTH( NCOLS,NLAYS,NPMAX,N_SPC_DIFF )
      REAL         S_HALO_NORTH( NCOLS,NLAYS,NPMAX,N_SPC_DIFF )
      REAL         S_HALO_WEST ( NROWS,NLAYS,NPMAX,N_SPC_DIFF )
      REAL         S_HALO_EAST ( NROWS,NLAYS,NPMAX,N_SPC_DIFF )
      REAL         S_BUF_NS    ( NCOLS*NLAYS*NPMAX*N_SPC_DIFF )
      REAL         S_BUF_EW    ( NROWS*NLAYS*NPMAX*N_SPC_DIFF )
#endif

#endif

      INTEGER, SAVE :: DIFF_MAP( N_SPC_DIFF )   ! global diff map to CGRID
      INTEGER, SAVE :: LOGDEV

      INTEGER      C, R, L, S, V, N            ! loop counters
      INTEGER      D2C

      CHARACTER( 96 ) :: XMSG = ' '
     
      INTEGER MY_TEMP
      INTEGER, SAVE :: STARTROW, ENDROW
      INTEGER, SAVE :: STARTCOL, ENDCOL

      INTERFACE
         SUBROUTINE RHO_J ( JDATE, JTIME, TSTEP, RHOJ )
            IMPLICIT NONE
            INTEGER, INTENT( IN )  :: JDATE, JTIME, TSTEP( 2 )
            REAL,    INTENT( OUT ) :: RHOJ( :,:,: )
         END SUBROUTINE RHO_J
         SUBROUTINE HCDIFF3D ( JDATE, JTIME, K11BAR, K22BAR, DT )
            IMPLICIT NONE
            INTEGER, INTENT( IN )  :: JDATE, JTIME
            REAL,    INTENT( OUT ) :: K11BAR( :,:,: ), K22BAR( :,:,: )
            REAL,    INTENT( OUT ) :: DT
         END SUBROUTINE HCDIFF3D
      END INTERFACE
 
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

C Get dx1 from COORD include file

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD
     &          * COS( PI180*( YORIG_GD + YCELL_GD*FLOAT( GL_NROWS/2 ))) ! in m.
            DX2 = DG2M * YCELL_GD   ! in m.
         ELSE
            DX1 = XCELL_GD          ! in m.
            DX2 = YCELL_GD          ! in m.
         END IF

         RDX1S = 1.0 / ( DX1 * DX1 )
         RDX2S = 1.0 / ( DX2 * DX2 )

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
 
C Create global map to CGRID
 
         S = 0
         DO V = 1, N_GC_DIFF
            S = S + 1
            DIFF_MAP( S ) = GC_STRT - 1 + GC_DIFF_MAP( V )
         END DO
         DO V = 1, N_AE_DIFF
            S = S + 1
            DIFF_MAP( S ) = AE_STRT - 1 + AE_DIFF_MAP( V )
         END DO
         DO V = 1, N_NR_DIFF
            S = S + 1
            DIFF_MAP( S ) = NR_STRT - 1 + NR_DIFF_MAP( V )
         END DO
         DO V = 1, N_TR_DIFF
            S = S + 1
            DIFF_MAP( S ) = TR_STRT - 1 + TR_DIFF_MAP( V )
         END DO
 
C open output ASCII file
 
         CALL SUBST_LOOP_INDEX ( 'C', 1, NCOLS, 1, MY_TEMP, STARTCOL, ENDCOL )
         CALL SUBST_LOOP_INDEX ( 'R', 1, NROWS, 1, MY_TEMP, STARTROW, ENDROW )

         WRITE( COMMSTR,'(4I2)' )  1, 0, 2, 0

      END IF                    ! if firstime
                                     
      DTSEC = FLOAT( TIME2SEC( TSTEP( 2 ) ) )
      FDATE = JDATE
      FTIME = JTIME
 
C Get the computational grid ( rho X Jacobian ) for this step

      CALL RHO_J ( FDATE, FTIME, TSTEP, RHOJ )

#ifdef parallel
C swap3d interface ... 
C     swap3d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             number of values to send in 3rd dimension of array,
C             size of 1st dimension in arrays,  
C             size of 2nd dimension in arrays,
C             direction to receive from: NSEW)

      CALL SWAP3D( RHOJ( 1,NROWS,1 ), RHOJ( 1,0,1 ), NCOLS, 1, NLAYS,
     &             NCOLS+2, NROWS+2, SOUTH )
      CALL SWAP3D( RHOJ( NCOLS,1,1 ), RHOJ( 0,1,1 ), 1, NROWS, NLAYS,
     &             NCOLS+2, NROWS+2, WEST )

#else
      CALL SUBST_COMM ( RHOJ, DSPL_N0_E0_S1_W1, DRCN_S_W, COMMSTR )
#endif

C get face values for RHOJ (assumes dx1 = dx2)

      DO L = 1, NLAYS
         DO R = STARTROW, ENDROW        !  DO R = 1, NROWS + 1
            DO C = STARTCOL, ENDCOL     !     DO C = 1, NCOLS + 1
               RK11( C,R,L ) = 0.5 * ( RHOJ( C,R,L ) + RHOJ( C-1,R,  L ) )
               RK22( C,R,L ) = 0.5 * ( RHOJ( C,R,L ) + RHOJ( C,  R-1,L ) )
            END DO
         END DO
      END DO

C Do the gridded computation for horizontal diffusion

C Get the contravariant eddy diffusivities

      CALL HCDIFF3D ( FDATE, FTIME, K11BAR3D, K22BAR3D, DT )

C get number of steps based on eddy time 
 
      NSTEPS = INT ( DTSEC / DT ) + 1
      DT = DTSEC / FLOAT( NSTEPS )
 
      WRITE( LOGDEV,1005 ) DT, NSTEPS

      DTDX1S = DT * RDX1S
      DTDX2S = DT * RDX2S

      DO L = 1, NLAYS
         DO R = STARTROW, ENDROW        !  DO R = 1, NROWS + 1
            DO C = STARTCOL, ENDCOL     !     DO C = 1, NCOLS + 1
               RK11( C,R,L ) = RK11( C,R,L ) * K11BAR3D( C,R,L )
               RK22( C,R,L ) = RK22( C,R,L ) * K22BAR3D( C,R,L )
            END DO
         END DO
      END DO

#ifdef parallel
      CALL SWAP3D( RK11( 1,1,1 ), RK11( NCOLS+1,1,1 ), 1, NROWS, NLAYS,
     &             NCOLS+1, NROWS+1, EAST )
      CALL SWAP3D( RK22( 1,1,1 ), RK22( 1,NROWS+1,1 ), NCOLS, 1, NLAYS,
     &             NCOLS+1, NROWS+1, NORTH )

C Create and fill four 3d boundary arrays with cgrid/rhoj
C    send to SOUTH, recv from NORTH, fill bound_north
C    send to NORTH, recv from SOUTH, fill bound_south
C    send to EAST, recv from WEST, fill bound_west
C    send to WEST, recv from EAST, fill bound_east

      M = 0
      DO S = 1, N_SPC_DIFF
         D2C = DIFF_MAP( S )
         DO L = 1, NLAYS
            DO C = 1, NCOLS
               HALO_SOUTH( C,L,S ) = CGRID( C,1,L,D2C ) / RHOJ( C,1,L )
               HALO_NORTH( C,L,S ) = CGRID( C,NROWS,L,D2C ) / RHOJ( C,NROWS,L )
               M = M + 1
               BUF_NS( M ) = HALO_NORTH( C,L,S )
            END DO
         END DO
      END DO

      CALL SWAP3D( HALO_SOUTH, HALO_NORTH, NCOLS, NLAYS, N_SPC_DIFF,
     &             NCOLS, NLAYS, NORTH )
      CALL SWAP3D( BUF_NS,     HALO_SOUTH, NCOLS, NLAYS, N_SPC_DIFF,
     &             NCOLS, NLAYS, SOUTH )


#ifdef sens
      M = 0
      DO S = 1, n_spc_diff
         D2C = DIFF_MAP( S )
         DO P = 1, NPMAX
           DO L = 1, NLAYS
             DO C = 1, NCOLS
                S_HALO_SOUTH( C,L,P,S ) = SENGRID( C,1,L,P,D2C ) / RHOJ( C,1,L )
                S_HALO_NORTH( C,L,P,S ) = SENGRID( C,NROWS,L,P,D2C ) / RHOJ( C,NROWS,L )
                M = M + 1
                S_BUF_NS( M ) = S_HALO_NORTH( C,L,P,S )
             END DO
           END DO
         END DO
      END DO

      CALL SWAP4D( S_HALO_SOUTH, S_HALO_NORTH, NCOLS, NLAYS,
     &             NPMAX, N_SPC_DIFF, NCOLS, NLAYS, NPMAX, NORTH )
      CALL SWAP4D( S_BUF_NS,     S_HALO_SOUTH, NCOLS, NLAYS,
     &             NPMAX, N_SPC_DIFF, NCOLS, NLAYS, NPMAX, SOUTH )

#endif

      M = 0
      DO S = 1, N_SPC_DIFF
         D2C = DIFF_MAP( S )
         DO L = 1, NLAYS
            DO R = 1, NROWS
               HALO_WEST( R,L,S ) = CGRID( 1,R,L,D2C ) / RHOJ( 1,R,L )
               HALO_EAST( R,L,S ) = CGRID( NCOLS,R,L,D2C ) / RHOJ( NCOLS,R,L )
               M = M + 1
               BUF_EW( M ) = HALO_EAST( R,L,S )
            END DO
         END DO
      END DO

      CALL SWAP3D( HALO_WEST, HALO_EAST, NROWS, NLAYS, N_SPC_DIFF,
     &             NROWS, NLAYS, EAST )
      CALL SWAP3D( BUF_EW,    HALO_WEST, NROWS, NLAYS, N_SPC_DIFF,
     &             NROWS, NLAYS, WEST )

#ifdef sens

      M = 0
      DO S = 1, N_SPC_DIFF
         D2C = DIFF_MAP( S )
         DO P = 1, NPMAX
           DO L = 1, NLAYS
             DO R = 1, NROWS
               S_HALO_WEST( R,L,P,S ) = SENGRID( 1,R,L,P,D2C ) / RHOJ( 1,R,L )
               S_HALO_EAST( R,L,P,S ) = SENGRID( NCOLS,R,L,P,D2C ) / RHOJ( NCOLS,R,L )
               M = M + 1
               S_BUF_EW( M ) = S_HALO_EAST( R,L,P,S )
             END DO
           END DO
         END DO
      END DO

      CALL SWAP4D( S_HALO_WEST, S_HALO_EAST, NROWS, NLAYS,
     &             NPMAX, N_SPC_DIFF, NROWS, NLAYS, NPMAX, EAST )
      CALL SWAP4D( S_BUF_EW,    S_HALO_WEST, NROWS, NLAYS,
     &             NPMAX, N_SPC_DIFF, NROWS, NLAYS, NPMAX, WEST )

#endif

#endif

C Loop over species, layers, nsteps

      DO 366 S = 1, N_SPC_DIFF
         D2C = DIFF_MAP( S )

         DO 355 L = 1, NLAYS

            DO 344 N = 1, NSTEPS

C Load working array (CGRID is coupled, CONC is mixing ratio)

               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CONC( C,R ) = CGRID( C,R,L,D2C ) / RHOJ( C,R,L )
                  END DO
               END DO

#ifdef parallel
C Fill 4 boundaries of working array CONC: south, north, west, east

               DO C = 1, NCOLS
                  CONC( C,0 )       = HALO_SOUTH( C,L,S )
                  CONC( C,NROWS+1 ) = HALO_NORTH( C,L,S )
               END DO

               DO R = 1, NROWS
                  CONC( 0,R )       = HALO_WEST( R,L,S )
                  CONC( NCOLS+1,R ) = HALO_EAST( R,L,S )
               END DO
#else

C South boundary

               R = 1
               DO C = 1, MY_NCOLS
                  CONC( C,R-1 ) = CONC( C,R )
               END DO

C North boundary

               R = MY_NROWS
               DO C = 1, MY_NCOLS
                  CONC( C,R+1 ) = CONC( C,R )
               END DO

C West boundary

               C = 1
               DO R = 1, MY_NROWS
                  CONC( C-1,R ) = CONC( C,R )
               END DO

C East boundary

               C = MY_NCOLS
               DO R = 1, MY_NROWS
                  CONC( C+1,R ) = CONC( C,R )
               END DO

               CALL SUBST_COMM ( CONC, DSPL_N1_E1_S1_W1, DRCN_N_E_S_W, COMMSTR )
               CALL SUBST_COMM ( RK11, DSPL_N0_E1_S0_W0, DRCN_E )
               CALL SUBST_COMM ( RK22, DSPL_N1_E0_S0_W0, DRCN_N )
#endif

C Update CGRID

               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS

                     CGRID( C,R,L,D2C ) = RHOJ( C,R,L ) * CONC( C,R )
     &                                  + DTDX1S
     &                                  * ( RK11( C+1,R,L )
     &                                  * ( CONC( C+1,R ) - CONC( C,R ) )
     &                                  - RK11( C,R,L )
     &                                    * ( CONC( C,R )   - CONC( C-1,R ) ) )
     &                                  + DTDX2S
     &                                  * ( RK22( C,R+1,L )
     &                                    * ( CONC( C,R+1 ) - CONC( C,R ) )
     &                                    - RK22( C,R,L )
     &                                    * ( CONC( C,R )   - CONC( C,R-1 ) ) )

                  END DO
               END DO

344         CONTINUE

355      CONTINUE
366   CONTINUE

#ifdef sens
c BEGIN SENSITIVITY CALLS
c SENGRID is treated identically to the way CGRID was treated above
      DO 466 S = 1, N_SPC_DIFF
        D2C = DIFF_MAP( S )

        DO 460 P = 1, NPMAX

          DO 455 L = 1, NLAYS

            DO 444 N = 1, NSTEPS

C Load working array (SENGRID is coupled, SENSWORK is mixing ratio)

            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  SENSWORK( C,R ) = SENGRID( C,R,L,P,D2C ) / RHOJ( C,R,L )
               END DO
            END DO

#ifdef parallel

C Fill 4 boundaries of working array SENSWORK: south, north, west, east

               DO C = 1, NCOLS
                  SENSWORK( C,0 ) = S_HALO_SOUTH( C,L,P,S )
                  SENSWORK( C,NROWS+1 ) = S_HALO_NORTH( C,L,P,S )
               END DO

               DO R = 1, NROWS
                  SENSWORK( 0,R ) = S_HALO_WEST( R,L,P,S )
                  SENSWORK( NCOLS+1,R ) = S_HALO_EAST( R,L,P,S )
               END DO

#else

c    South boundary

               R = 1
               DO C = 1, MY_NCOLS
                  SENSWORK( C,R-1 ) = SENSWORK( C,R )
               END DO

c       North boundary

               R = MY_NROWS
               DO C = 1, MY_NCOLS
                  SENSWORK( C,R+1 ) = SENSWORK( C,R )
               END DO

c       West boundary

               C = 1
               DO R = 1, MY_NROWS
                  SENSWORK( C-1,R ) = SENSWORK( C,R )
               END DO

c       East boundary

               C = MY_NCOLS
               DO R = 1, MY_NROWS
                  SENSWORK( C+1,R ) = SENSWORK( C,R )
               END DO

               CALL SUBST_COMM (SENSWORK,DSPL_N1_E1_S1_W1,DRCN_N_E_S_W,COMMSTR)

! The following sharp_endif refers to the ifdef parallel within the sens code
#endif


               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS

                     SENGRID( C,R,L,P,D2C ) = RHOJ( C,R,L ) * SENSWORK( C,R)
     &                               + DTDX1S
     &                               * ( RK11( C+1,R,L )
     &                               * ( SENSWORK( C+1,R ) - SENSWORK( C,R ) )
     &                               - RK11( C,R,L )
     &                               * ( SENSWORK( C,R )   - SENSWORK( C-1,R )))
     &                               + DTDX2S
     &                               * ( RK22( C,R+1,L )
     &                               * ( SENSWORK( C,R+1 ) - SENSWORK( C,R ) )
     &                               - RK22( C,R,L )
     &                               * ( SENSWORK( C,R )   - SENSWORK( C,R-1 )))


                     END DO
                  END DO

444            CONTINUE


455         CONTINUE

460        CONTINUE

466      CONTINUE
! END SENSITIVITY CALLS
#endif

      RETURN

1001  FORMAT( 5X, 'Negative concentrations reset to', 1PE11.3 )
1003  FORMAT( 1X, 4I7, 9X, 1PE11.3)
1005  FORMAT( / 5X, 'H-eddy DT & integration steps: ', 1PE15.7, I8 )

      END
