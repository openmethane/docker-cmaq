
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo/y_yamo.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE Y_YAMO ( CGRID, FDATE, FTIME, TSTEP, LVL, BCON, YTRHOJ )

C-----------------------------------------------------------------------
C Function:
C   Advection in the horizontal plane, Y-direction first
      
C Preconditions:
      
C Subroutines and functions called:
 
C Revision history:
C  19 Jan 2004: Jeff Young

C   1 Nov 06: Jeff Young - Following Glenn Hammond, moved all communication
C   out of HPPM to this level; using "swap_sandia" communication; update only
C   local values in the CGRID array within a time step, discarding previous
C   ghost values.
C   7 Jul 2007: Jeff Young - declare YTRHOJ dimensions properly in order to retain
C                            them as declared in HADV
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
      
      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER, INTENT( IN ) :: FDATE         ! current model date, coded YYYYDDD
      INTEGER, INTENT( IN ) :: FTIME         ! current model time, coded HHMMSS
      INTEGER, INTENT( IN ) :: TSTEP         ! time step (HHMMSS)
      INTEGER, INTENT( IN ) :: LVL           ! layer
!     REAL,    INTENT( IN ) :: BCON( NBNDY,* ) ! boundary concentrations
      REAL,    INTENT( IN ) :: BCON( :,: )     ! boundary concentrations
!     REAL,    INTENT( IN ) :: YTRHOJ( NCOLS,0:NROWS )
      REAL,    INTENT( IN ) :: YTRHOJ( 1:,0: ) ! met RHOJ

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
 
      REAL         VHAT( NCOLS+1,NROWS+1 )       ! x1-component CX-velocity

C Local Variables:

      REAL         TRRHOJ( 0:NROWS )

      CHARACTER( 16 ) :: PNAME = 'Y_YAMO'
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      CHARACTER( 96 ) :: XMSG = ' '

      integer, save :: logdev

      REAL,    SAVE :: DX2                      ! dx1 (meters)
      INTEGER, SAVE :: ASPC                     ! RHOJ index in CGRID
      REAL          ADJFAC                      ! Yamo's adjustment - jp 11/02
      REAL          ADJFAC_L

      REAL, ALLOCATABLE, SAVE :: VELY( : ),     ! Velocities along a column
     &                           CONY( :,: )    ! Conc's along a column

      REAL          DT                          ! TSTEP in sec
      INTEGER       ALLOCSTAT

      INTEGER, SAVE :: ADV_MAP( N_SPC_ADV )     ! global adv map to CGRID

      CHARACTER( 16 ) :: X2VEL = 'X2VEL'

      INTEGER      COL, ROW, SPC, VAR           ! loop counters
      INTEGER      A2C

      INTEGER MY_TEMP
!     INTEGER, SAVE :: STARTROW, ENDROW, STARTROW2, ENDROW2, ENDROW2_LIM
      INTEGER, SAVE :: STARTROW2, ENDROW2, ENDROW2_LIM
      LOGICAL, SAVE :: BNDY_PE_LOY, BNDY_PE_HIY

#ifdef parallel
      INTEGER NORTH_ROW, M
!     REAL    HALO_NORTH( NCOLS,SWP,N_SPC_ADV )
!     REAL    HALO_SOUTH( NCOLS,SWP,N_SPC_ADV )
!     REAL    BUF_NS    ( SWP*NCOLS*N_SPC_ADV )
      REAL, ALLOCATABLE, SAVE :: HALO_NORTH( :,:,: )
      REAL, ALLOCATABLE, SAVE :: HALO_SOUTH( :,:,: )
      REAL, ALLOCATABLE, SAVE :: BUF_NS( : )
#endif

C Statement functions:

      INTEGER, SAVE :: SFX    ! fixed parameter for southern boundary
      INTEGER, SAVE :: NFX    ! fixed parameter for northern boundary

      REAL    BCCN            ! boundary concentrations stmt fn

      INTEGER BFX             ! dummy positional parameter
      INTEGER CR              ! row or column index
      INTEGER SS              ! species index

      BCCN ( BFX, CR, SS ) = BCON( BFX + CR, SS )

C Required interface for allocatable array dummy arguments

      INTERFACE
         SUBROUTINE HCONTVEL ( FDATE, FTIME, TSTEP, LVL, UORV, UHAT )
            IMPLICIT NONE
            INTEGER,         INTENT( IN )     :: FDATE, FTIME, TSTEP, LVL
            CHARACTER( 16 ), INTENT( IN )     :: UORV
            REAL,            INTENT( OUT )    :: UHAT( :,: )
         END SUBROUTINE HCONTVEL
         SUBROUTINE HPPM ( NI, NJ, CON, VEL, DT, DS, ORI ) 
            IMPLICIT NONE
#ifdef parallel
            INTEGER, PARAMETER                :: SWP = 3
#else
            INTEGER, PARAMETER                :: SWP = 1
#endif
            INTEGER,         INTENT( IN )     :: NI, NJ
            REAL,            INTENT( IN OUT ) :: CON( 1-SWP:,1: )
            REAL,            INTENT( IN )     :: VEL( : )
            REAL,            INTENT( IN )     :: DT
            REAL,            INTENT( IN )     :: DS
            CHARACTER,       INTENT( IN )     :: ORI
         END SUBROUTINE HPPM
      END INTERFACE
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         logdev = init3 ()

         SFX = 0
         NFX = MY_NCOLS + MY_NROWS + 3

C Get dx2 from HGRD_DEFN module

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX2 = DG2M * YCELL_GD   ! in m.
         ELSE
            DX2 = YCELL_GD          ! in m.
         END IF

         ALLOCATE ( CONY( 1-SWP:MY_NROWS+SWP,N_SPC_ADV ),
     &              VELY( MY_NROWS+1 ), STAT = ALLOCSTAT ) ! Vel along a col
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating VELY, or CONY'
            CALL M3EXIT ( PNAME, FDATE, FTIME, XMSG, XSTAT1 )
         END IF

#ifdef parallel
         ALLOCATE ( HALO_NORTH( MY_NCOLS,SWP,N_SPC_ADV ),
     &              HALO_SOUTH( MY_NCOLS,SWP,N_SPC_ADV ),
     &              BUF_NS    ( SWP*MY_NCOLS*N_SPC_ADV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating HALO_NORTH, HALO_SOUTH, or BUF_NS'
            CALL M3EXIT ( PNAME, FDATE, FTIME, XMSG, XSTAT1 )
         END IF
         HALO_NORTH = 0.0   ! array
         HALO_SOUTH = 0.0   ! array
         BUF_NS     = 0.0   ! array
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

!        CALL SUBST_LOOP_INDEX ( 'R', 1, MY_NROWS, 1, MY_TEMP,
!    &                           STARTROW, ENDROW )

         CALL SUBST_LOOP_INDEX ( 'R', 2, MY_NROWS, 1, MY_TEMP,
     &                           STARTROW2, ENDROW2 )

         CALL SUBST_HI_LO_BND_PE ( 'R', BNDY_PE_LOY, BNDY_PE_HIY )

         IF ( BNDY_PE_HIY ) THEN
            ENDROW2_LIM = ENDROW2
         ELSE
            ENDROW2_LIM = ENDROW2 + 1
         END IF

!        write( logdev,2005 ) bndy_pe_loy, bndy_pe_hiy,
!    &                        startrow, endrow, startrow2, endrow2, endrow2_lim
!        call flush_
2005  format( '@@ #y_yamo> :', 2l6, 2x, 5i6.3 )

      END IF                    ! if firstime

      DT = FLOAT ( TIME2SEC ( TSTEP ) )

C Do the computation for y advection

C Get the contravariant x2 velocity component

      CALL HCONTVEL ( FDATE, FTIME, TSTEP, LVL, X2VEL, VHAT )

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

      CALL SWAP2D( VHAT( 1,1 ), VHAT( 1,MY_NROWS+1 ), MY_NCOLS+1, 1,
     &             MY_NCOLS+1, NORTH )
      M = 0
      NORTH_ROW = MY_NROWS - SWP
      DO SPC = 1, N_SPC_ADV
         A2C = ADV_MAP( SPC )
         DO ROW = 1, SWP
            DO COL = 1, MY_NCOLS
               HALO_SOUTH( COL,ROW,SPC ) = CGRID( COL,ROW,LVL,A2C )
               HALO_NORTH( COL,ROW,SPC ) = CGRID( COL,NORTH_ROW+ROW,LVL,A2C )
               M = M + 1
               BUF_NS( M ) = HALO_NORTH( COL,ROW,SPC )
            END DO
         END DO
      END DO

      CALL SWAP3D( HALO_SOUTH, HALO_NORTH, MY_NCOLS, SWP, N_SPC_ADV,
     &             MY_NCOLS, SWP, NORTH )
      CALL SWAP3D( BUF_NS,     HALO_SOUTH, MY_NCOLS, SWP, N_SPC_ADV,
     &             MY_NCOLS, SWP, SOUTH )
#endif

      DO 253 COL = 1, MY_NCOLS

!        DO ROW = STARTROW, ENDROW      ! DO ROW = 1, MY_NROWS+1
         DO ROW = 1, MY_NROWS+1
            VELY( ROW ) = VHAT( COL,ROW )
         END DO

         TRRHOJ( 1:MY_NROWS ) = CGRID( COL,1:MY_NROWS,LVL,ASPC )

#ifdef parallel
         TRRHOJ( 0 ) = HALO_SOUTH( COL,SWP,N_SPC_ADV )
#else
         TRRHOJ( 0 ) = TRRHOJ( 1 )
#endif

         IF ( BNDY_PE_LOY ) THEN
            ADJFAC = YTRHOJ( COL,1 ) / TRRHOJ( 1 )
            IF ( VELY( 1 ) .LT. 0.0 ) VELY( 1 ) = VELY( 1 ) * ADJFAC
         END IF

         DO ROW = STARTROW2, ENDROW2

            IF ( ROW .EQ. STARTROW2 ) THEN
               ADJFAC_L = YTRHOJ( COL,ROW-1 ) / TRRHOJ( ROW-1 )
            ELSE
               ADJFAC_L = ADJFAC
            END IF

            IF ( ROW .LT. ENDROW2_LIM ) ADJFAC = YTRHOJ( COL,ROW )
     &                                         / TRRHOJ( ROW )

            IF ( VELY( ROW ) .LT. 0.0 ) THEN
               VELY( ROW ) = VELY( ROW ) * ADJFAC
            ELSE IF ( VELY( ROW ) .GT. 0.0 ) THEN
               VELY( ROW ) = VELY( ROW ) * ADJFAC_L
            END IF

!           if ( vely( row ) .lt. 0.0 ) then
!              if ( adjfac .lt. 0.70 .or. adjfac .gt. 1.36 ) then
!                 write( logdev,2009 ) col, row, lvl, adjfac, vely(row)
!                 end if  !
!              else
!              if ( adjfac_l .lt. 0.70 .or. adjfac_l .gt. 1.36 ) then
!                 write( logdev,2011 ) col, row, lvl, adjfac_l, vely(row)
!                 end if  !
!              end if
2009           format(  '#y_yamo> col, row, lvl, adjfac,   vely: ',
     &                3i3, 1x, 2f7.2 )
2011           format(  '#y_yamo> col, row, lvl, adjfac_l, vely: ',
     &                3i3, 1x, 2f7.2 )

         END DO

#ifdef parallel
         CALL SWAP1D( VELY( 1 ), VELY( MY_NROWS+1 ), 1, NORTH )
#endif

         DO SPC = 1, N_SPC_ADV

            A2C = ADV_MAP( SPC )
            DO ROW = 1, MY_NROWS
               CONY( ROW,SPC ) = CGRID( COL,ROW,LVL,A2C )
            END DO

#ifdef parallel
            DO ROW = 1, SWP
               CONY( ROW-SWP,SPC )      = HALO_SOUTH( COL,ROW,SPC )
               CONY( MY_NROWS+ROW,SPC ) = HALO_NORTH( COL,ROW,SPC )
            END DO
#endif

C South boundary

            IF ( BNDY_PE_LOY ) THEN
               IF ( VELY( 1 ) .LT. 0.0 ) THEN          ! outflow
                  CONY( 1-SWP:0,SPC) =
     &               ZFDBC ( CONY( 1,SPC ), CONY( 2,SPC ),
     &                       VELY( 1 ),     VELY( 2 ) )
               ELSE    ! inflow
                  CONY( 1-SWP:0,SPC ) = BCCN( SFX,COL,SPC )
               END IF
            END IF

C North boundary

            IF ( BNDY_PE_HIY ) THEN
               IF ( VELY( MY_NROWS+1 ) .GT. 0.0 ) THEN     ! outflow
                  CONY( MY_NROWS+1:MY_NROWS+SWP,SPC ) =
     &               ZFDBC ( CONY( MY_NROWS,SPC ), CONY( MY_NROWS-1,SPC ),
     &                       VELY( MY_NROWS+1 ),   VELY( MY_NROWS ) )
               ELSE    ! inflow
                  CONY( MY_NROWS+1:MY_NROWS+SWP,SPC ) = BCCN( NFX,COL,SPC )
               END IF
            END IF

         END DO

!        if ( lvl .eq. 1 .and. col .eq. 22 ) then
!           write( logdev,* ) '@@ y_yamo# time, row, vely'
!           write( logdev,* ) '@@ y_yamo# time, row, trrhoj'
!           do row = 20, my_nrows+1
!           do row = 0, my_nrows
!              write( logdev,2017 ) ftime, row, vely( row )
!              write( logdev,2017 ) ftime, row, trrhoj( row )
!           end do
!        end if
2017  format( '@@', 5x, i7.6, i4, 1pe12.3 )

C PPM scheme

         CALL HPPM ( MY_NROWS, MY_NCOLS, CONY, VELY, DT, DX2, 'R' )

         DO SPC = 1, N_SPC_ADV
            A2C = ADV_MAP( SPC )
            DO ROW = 1, MY_NROWS
               CGRID( COL,ROW,LVL,A2C ) = CONY( ROW,SPC )
!              if ( cony( row,spc ) .le. 0.0 )
!    &            write( logdev,2019 ) ftime, col, row, lvl, spc, cony(row,spc )
            END DO
         END DO

2019  format( 'y_yamo# time, c, r, l, s, cony: ',
     &         I7.6, 4I4, 1pe12.3 )

253   CONTINUE

      RETURN
      END
