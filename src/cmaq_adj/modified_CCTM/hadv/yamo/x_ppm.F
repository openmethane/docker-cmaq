
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo/x_ppm.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE X_PPM ( CGRID, FDATE, FTIME, TSTEP, LVL, BCON )

C-----------------------------------------------------------------------
C Function:
C   Piecewise Parabolic Method advection in the X-direction

C Preconditions:

C Subroutines and functions called:

C Revision history:
C  28 Jun 2004: Jeff Young

C   1 Nov 06: Jeff Young - Following Glenn Hammond, moved all communication
C   out of HPPM to this level; using "swap_sandia" communication; update only
C   local values in the CGRID array within a time step, discarding previous
C   ghost values.
C   11 May 2009: Jeff Young: Simplify - assume constant cell widths, DS( i )
C   11 May 2010: Jeff Young: New hppm fix for PGI compiler by David Wong
C-----------------------------------------------------------------------

      USE HGRD_DEFN             ! horizontal domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE SUBST_MODULES         ! stenex
!     USE SUBST_COMM_MODULE     ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex

#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Arguments:

      REAL,    POINTER      :: CGRID( :,:,:,: )
      INTEGER, INTENT( IN ) :: FDATE         ! current model date, coded YYYYDDD
      INTEGER, INTENT( IN ) :: FTIME         ! current model time, coded HHMMSS
      INTEGER, INTENT( IN ) :: TSTEP         ! time step (HHMMSS)
      INTEGER, INTENT( IN ) :: LVL           ! layer
!     REAL,    INTENT( IN ) :: BCON( NBNDY,1: ) ! boundary concentrations
      REAL,    INTENT( IN ) :: BCON( :,: )      ! boundary concentrations

C External Functions not declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC
      REAL,    EXTERNAL :: ZFDBC

C Parameters:

C Advected species dimension

      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV
     &                                + N_AE_ADV
     &                                + N_NR_ADV
     &                                + N_TR_ADV
     &                                + 1       ! for advecting RHO*SqRDMT

#ifdef parallel
      INTEGER, PARAMETER :: SWP = 3
#else
      INTEGER, PARAMETER :: SWP = 1
#endif

C File Variables:

      REAL         UHAT( NCOLS+1,NROWS+1 )       ! x1-component CX-velocity

C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'X_PPM'
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      CHARACTER( 96 ) :: XMSG = ' '

      integer, save :: logdev

      REAL,    SAVE :: DX1                      ! dx1 (meters)
      INTEGER, SAVE :: ASPC                     ! RHOJ index in CGRID

      REAL, ALLOCATABLE, SAVE :: VELX( : ),     ! Velocities along a row
     &                           CONX( :,: )    ! Conc's along a row

      REAL          DT                          ! TSTEP in sec
      INTEGER       ALLOCSTAT

      INTEGER, SAVE :: ADV_MAP( N_SPC_ADV )     ! global adv map to CGRID

      CHARACTER( 16 ) :: X1VEL = 'X1VEL'

      INTEGER      COL, ROW, SPC, VAR           ! loop counters
      INTEGER      A2C

      INTEGER MY_TEMP
!     INTEGER, SAVE :: STARTCOL, ENDCOL
      LOGICAL, SAVE :: BNDY_PE_LOX, BNDY_PE_HIX

#ifdef parallel
      INTEGER EAST_COL, M
!     REAL    HALO_EAST( SWP,NROWS,N_SPC_ADV )
!     REAL    HALO_WEST( SWP,NROWS,N_SPC_ADV )
!     REAL    BUF_EW   ( SWP*NROWS*N_SPC_ADV )
      REAL, ALLOCATABLE, SAVE :: HALO_EAST( :,:,: )
      REAL, ALLOCATABLE, SAVE :: HALO_WEST( :,:,: )
      REAL, ALLOCATABLE, SAVE :: BUF_EW( : )
#endif

C Statement functions:

      INTEGER, SAVE :: EFX    ! fixed parameter for eastern boundary
      INTEGER, SAVE :: WFX    ! fixed parameter for western boundary

      REAL    BCCN            ! boundary concentrations stmt fn

      INTEGER BFX             ! dummy positional parameter
      INTEGER CR              ! row or column index
      INTEGER SS              ! species index

      BCCN ( BFX, CR, SS ) = BCON( BFX + CR, SS )

C Required interface for allocatable array dummy arguments

      INTERFACE
         SUBROUTINE HCONTVEL( FDATE, FTIME, TSTEP, LVL, UORV, UHAT )
            IMPLICIT NONE
            INTEGER, INTENT( IN )         :: FDATE, FTIME, TSTEP, LVL
            CHARACTER( 16 ), INTENT( IN ) :: UORV
            REAL,    INTENT( OUT )        :: UHAT( :,: )
         END SUBROUTINE HCONTVEL
         SUBROUTINE HPPM ( NI, NJ, CON, VEL, DT, DS, ORI )
            IMPLICIT NONE
#ifdef parallel
            INTEGER, PARAMETER         :: SWP = 3
#else
            INTEGER, PARAMETER         :: SWP = 1
#endif
            INTEGER,   INTENT( IN )    :: NI, NJ
            REAL,      INTENT( INOUT ) :: CON( 1-SWP:,1: )
            REAL,      INTENT( IN )    :: VEL( : )
            REAL,      INTENT( IN )    :: DT
            REAL,      INTENT( IN )    :: DS
            CHARACTER, INTENT( IN )    :: ORI
         END SUBROUTINE HPPM
      END INTERFACE
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         logdev = init3 ()

         EFX = MY_NCOLS + 1
         WFX = 2 * MY_NCOLS + MY_NROWS + 4

C Get dx1 from HGRD_DEFN module

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD
     &          * COS( PI180*( YORIG_GD + YCELL_GD*FLOAT( GL_NROWS/2 ))) ! in m.
         ELSE
            DX1 = XCELL_GD        ! in m.
         END IF

         ALLOCATE ( CONX( 1-SWP:MY_NCOLS+SWP,N_SPC_ADV ),
     &              VELX( MY_NCOLS+1 ), STAT = ALLOCSTAT ) ! Vel along a row
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating VELX, or CONX'
            CALL M3EXIT ( PNAME, FDATE, FTIME, XMSG, XSTAT1 )
         END IF

#ifdef parallel
         ALLOCATE ( HALO_EAST( SWP,MY_NROWS,N_SPC_ADV ),
     &              HALO_WEST( SWP,MY_NROWS,N_SPC_ADV ),
     &              BUF_EW   ( SWP*MY_NROWS*N_SPC_ADV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating HALO_EAST, HALO_WEST, or BUF_EW'
            CALL M3EXIT ( PNAME, FDATE, FTIME, XMSG, XSTAT1 )
         END IF
         HALO_EAST = 0.0   ! array
         HALO_WEST = 0.0   ! array
         BUF_EW    = 0.0   ! array
#endif

C Get CGRID offsets

         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

         ASPC = GC_STRT - 1 + N_GC_SPCD

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

         ADV_MAP( N_SPC_ADV ) = ASPC

!        CALL SUBST_LOOP_INDEX ( 'C', 1, MY_NCOLS, 1, MY_TEMP,
!    &                           STARTCOL, ENDCOL )

!        write( logdev,* ) 'mype, my_ncols, startcol, endcol: ',
!    &                      mype, my_ncols, startcol, endcol

         CALL SUBST_HI_LO_BND_PE ( 'C', BNDY_PE_LOX, BNDY_PE_HIX )

         END IF                    ! if firstime

      DT = FLOAT ( TIME2SEC ( TSTEP ) )

C Do the computation for x advection

C Get the contravariant x1 velocity component

      CALL HCONTVEL ( FDATE, FTIME, TSTEP, LVL, X1VEL, UHAT )

#ifdef parallel
C SWAPnD interface (same array in both procs)
C                                                     SWAP1D SWAP2D SWAP3D
C 1st value of array to be sent,                         X     X      X
C 1st value of array to be received,                     X     X      X
C number of values to send in 1st dimension of array,    X     X      X
C number of values to send in 2nd dimension of array,          X      X
C number of values to send in 3rd dimension of array,                 X
C size of 1st dimension in array,                              X      X
C size of 2nd dimension in array,                                     X
C direction to receive from: NORTH,SOUTH,EAST,WEST       X     X      X

      CALL SWAP2D( UHAT( 1,1 ), UHAT( NCOLS+1,1 ), 1, NROWS+1, NCOLS+1, EAST )
      M = 0
      EAST_COL = MY_NCOLS - SWP
      DO SPC = 1, N_SPC_ADV
         A2C = ADV_MAP( SPC )
         DO ROW = 1, MY_NROWS
            DO COL = 1, SWP
               HALO_WEST( COL,ROW,SPC ) = CGRID( COL,ROW,LVL,A2C )
               HALO_EAST( COL,ROW,SPC ) = CGRID( EAST_COL+COL,ROW,LVL,A2C )
               M = M + 1
               BUF_EW( M ) = HALO_EAST( COL,ROW,SPC )
            END DO
         END DO
      END DO

      CALL SWAP3D( HALO_WEST, HALO_EAST, SWP, MY_NROWS, N_SPC_ADV,
     &             SWP, MY_NROWS, EAST )
      CALL SWAP3D( BUF_EW,    HALO_WEST, SWP, MY_NROWS, N_SPC_ADV,
     &             SWP, MY_NROWS, WEST )
#endif

      DO 233 ROW = 1, MY_NROWS

!        DO COL = STARTCOL, ENDCOL         !    DO COL = 1, MY_NCOLS+1
         DO COL = 1, MY_NCOLS+1
            VELX( COL ) = UHAT( COL,ROW )
         END DO

         DO SPC = 1, N_SPC_ADV

            A2C = ADV_MAP( SPC )
            DO COL = 1, MY_NCOLS
               CONX( COL,SPC ) = CGRID( COL,ROW,LVL,A2C )
            END DO

#ifdef parallel
            DO COL = 1, SWP
               CONX( COL-SWP,SPC )      = HALO_WEST( COL,ROW,SPC )
               CONX( MY_NCOLS+COL,SPC ) = HALO_EAST( COL,ROW,SPC )
            END DO
#endif

C West boundary

            IF ( BNDY_PE_LOX ) THEN
               IF ( VELX( 1 ) .LT. 0.0 ) THEN          ! outflow
                  CONX( 1-SWP:0,SPC ) =
     &               ZFDBC ( CONX( 1,SPC ), CONX( 2,SPC ),
     &                       VELX( 1 ),     VELX( 2 ) )
               ELSE    ! inflow
                  CONX( 1-SWP:0,SPC ) = BCCN( WFX,ROW,SPC )
               END IF
            END IF

C East boundary

            IF ( BNDY_PE_HIX ) THEN
               IF ( VELX( MY_NCOLS+1 ) .GT. 0.0 ) THEN     ! outflow
                  CONX( MY_NCOLS+1:MY_NCOLS+SWP,SPC ) =
     &               ZFDBC ( CONX( MY_NCOLS,SPC ), CONX( MY_NCOLS-1,SPC ),
     &                       VELX( MY_NCOLS+1 ),   VELX( MY_NCOLS ) )
               ELSE    ! inflow
                  CONX( MY_NCOLS+1:MY_NCOLS+SWP,SPC ) = BCCN( EFX,ROW,SPC )
               END IF
            END IF

         END DO

C PPM scheme

         CALL HPPM ( MY_NCOLS, MY_NROWS, CONX, VELX, DT, DX1, 'C' )

         DO SPC = 1, N_SPC_ADV
            A2C = ADV_MAP( SPC )
            DO COL = 1, MY_NCOLS
               CGRID( COL,ROW,LVL,A2C ) = CONX( COL,SPC )
!              if ( conx( col,spc ) .le. 0.0 )
!    &            write( logdev,2019 ) ftime, col, row, lvl, spc, conx(col,spc )
            END DO
         END DO

!        do col = 1, my_ncols
!           if ( conx( col,n_spc_adv ) .le. 0.0 ) then
!              write( logdev,2019 ) ftime, col, row, lvl,
!    &                              conx( col,n_spc_adv )
!           end if
!        end do
!019  format( 'x_ppm# date, time, col, row, lvl, conx: ',
!    &         2I7, 3I4, 1pe12.3 )
2019  format( 'x_ppm# time, c, r, l, s, conx: ',
     &         I7.6, 4I4, 1pe12.3 )

233   CONTINUE

      RETURN
      END
