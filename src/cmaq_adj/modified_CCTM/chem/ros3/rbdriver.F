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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ros3/rbdriver.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


       SUBROUTINE CHEM( CONC, JDATE, JTIME, TSTEP )

C**********************************************************************
C
C  FUNCTION: To control gas phase chemistry calculations performed by
C            the vectorized Rosenbrock solver
C
C  PRECONDITIONS: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: RBINIT
C                                    RBSPARSE
C                                    PHOT
C                                    CALCKS
C                                    RBSOLVER
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, August, 2004
C                    Based on the solver described by Sandu et al
C                    ( Atm. Env., Vol. 31, No. 20, 1997 ) and included
C                    in the Kinetic PreProcessor ( see for example 
C                    Sandu et al., At, Env., Vol. 37, 5097-5114, 
C                    2003). This code also incorporates efficiency
C                    concepts originally developed by M. Jacobson
C                    for SMVGEAR (Atm. Env., Vol 28, No 2, 1994).
C                    Adapted from Subroutine CHEM in CMAQ SMVGEAR
C
C                    31 Jan 05 J.Young: dyn alloc - establish both horizontal
C                    & vertical domain specifications in one module (GRID_CONF)
C                    29 Jul 05     WTH: Added IF blocks that call degrade 
C                                       routines if MECHNAME contains 'TX' 
C                                       substring.
C**********************************************************************

      USE SUBST_MODULES
      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE RBDATA                ! ROS3 solver data
      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE 

C..INCLUDES:

!     INCLUDE SUBST_VGRD_ID     ! Vertical grid data
      INCLUDE SUBST_IOPARMS     ! Io/api parameters
      INCLUDE SUBST_IOFDESC     ! Io/api file descriptions
#include      SUBST_IODECL      ! Io/api declarations
      INCLUDE SUBST_FILES_ID    ! CMAQ files
!     INCLUDE SUBST_COORD_ID    ! Coordinate and domain defs (req IOPARMS)
      INCLUDE SUBST_CONST       ! CMAQ constants
      INCLUDE SUBST_GC_SPC      ! Gas chem species names and MWs
      INCLUDE SUBST_RXCMMN      ! Mechanism reaction common block
      INCLUDE SUBST_PACTL_ID    ! Process analysis control parameters
 
     
C..ARGUMENTS:

      REAL, POINTER :: CONC( :,:,:,: )  ! Concentrations

      INTEGER JDATE                     ! Current date (YYYYDDD)
      INTEGER JTIME                     ! Current time (HHMMSS)
      INTEGER TSTEP( 2 )                ! Time step vector (HHMMSS)


C..PARAMETERS:

      INTEGER, PARAMETER :: IZERO = 0                ! Integer zero

      REAL, PARAMETER    :: CONCMIN = 1.0E-30        ! Minimum conc
      REAL, PARAMETER    :: CONCOFM = 1.0E+06        ! conc. of M = 1E+06 ppm
      REAL, PARAMETER    :: PA2ATM  = 1.0 / STDATMPA ! Pascal to atm conv fac


C..EXTERNAL FUNCTIONS:

      INTEGER SEC2TIME                 ! Returns time interval from seconds
      INTEGER TIME2SEC                 ! Returns seconds in time interval


C..SAVED LOCAL VARIABLES:

      LOGICAL, SAVE :: LFIRST = .TRUE. ! Flag for first call to this subroutine
      LOGICAL, SAVE :: LIRRBLK         ! Flag for IRR to be done for block

      INTEGER, SAVE :: LOGDEV          ! Unit number of output log
      INTEGER, SAVE :: NOXYZ           ! Total number of grid cells

      REAL,    SAVE :: AIRFC           ! Factor to convert gms air to ppm
      REAL,    SAVE :: MAOMV           ! Mol Wt of air over Mol Wt of water
!     REAL,    SAVE :: CONV( NLAYS )   ! Volumetric emis conv facs

      
C..SCRATCH LOCAL VARIABLES:

      CHARACTER(  16 ) :: PNAME = 'RBDRIVER' ! Routine name
      CHARACTER(  16 ) :: VNAME              ! Name of I/O API data variable
      CHARACTER( 144 ) :: MSG                ! Message text
     
      INTEGER C, R, L, S      ! Loop indices
      INTEGER ALLOCSTAT       ! Allocate status code
      INTEGER OFFSET          ! Starting cell number of a block
      INTEGER NCSP            ! Mech no: 1=gas/day 2=gas/night
      INTEGER BLK             ! Loop index for block of cells
      INTEGER CELLNUM         ! Cell number 
      INTEGER COL             ! Column index
      INTEGER IPAR            ! Pointer for cell sort routine
      INTEGER IRVAL           ! Pointer for cell sort routine
      INTEGER IRXN            ! Reaction number
      INTEGER ISP             ! Species index
      INTEGER ISPOLD          ! Species number in original order
      INTEGER ISPNEW          ! Species number in new sorted order 
      INTEGER ITMSTEP         ! Chemistry integration interval (sec)   
      INTEGER JPAR            ! Pointer for cell sort routine
      INTEGER JREORD          ! Index holder for sort routine
      INTEGER LEV             ! Layer index
      INTEGER LVAL            ! Pointer for cell sort routine
      INTEGER MIDDATE         ! Date at time step midpoint
      INTEGER MIDTIME         ! Time at time step midpoint
      INTEGER NCELL           ! Index for number of cells
      INTEGER NDARK           ! Number of layer 1 cells in darkness
      INTEGER NIRRCLS         ! No. of cells in block for IRR
      INTEGER NPH             ! Index for number of phot. rxns in PHOT
      INTEGER NRX             ! Index for number of reactions
      INTEGER ROW             ! Row index
      INTEGER SPC             ! Species loop index
      INTEGER VAR             ! Variable number on I/O API file

      INTEGER IRRCELL( BLKSIZE )               ! Cell No. of an IRR cell

C for INTERPX
      INTEGER     GXOFF, GYOFF                 ! global origin offset from file
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER, SAVE :: STRTCOLEM1, ENDCOLEM1, STRTROWEM1, ENDROWEM1
 
      REAL CHEMSTEP                            ! Chem integration interval (min)
      REAL VALLOW                              ! Value holder for sort routine
      
      REAL DENSA_J( NCOLS, NROWS, NLAYS )      ! Cell density (Kg/m**3)
      REAL DENS(    NCOLS, NROWS, NLAYS )      ! Cell density (Kg/m**3)
      REAL PRES(    NCOLS, NROWS, NLAYS )      ! Cell pressure (Pa)
      REAL QV(      NCOLS, NROWS, NLAYS )      ! Cell water vapor (Kg/Kg air)
      REAL TA(      NCOLS, NROWS, NLAYS )      ! Cell temperature (K)

      REAL RJBLK( BLKSIZE, NPHOTAB )           ! J-values for each cell in 
                                               ! a block

      REAL RJ( NCOLS, NROWS, NLAYS, NPHOTAB )  ! J-values for each cell      



#ifdef rbstats

      CHARACTER( 16 ), SAVE              :: CTM_RBSTATS_1 = 'CTM_RBSTATS_1' 
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: VSTATS( : )        !

      INTEGER, SAVE                      :: WSTEP = 0     
      INTEGER, ALLOCATABLE, SAVE         :: STAT_SUM( : , : , : , : )
      INTEGER  EDATE, ETIME

      REAL                               :: STATOUT( NCOLS, NROWS, NLAYS )

#endif


      INTERFACE
!        SUBROUTINE PHOT ( MDATE, MTIME, JDATE, JTIME, NDARK, RJ )
!           IMPLICIT NONE
!           INTEGER, INTENT( IN )  :: MDATE, MTIME, JDATE, JTIME
!           INTEGER, INTENT( OUT ) :: NDARK
!           REAL,    INTENT( OUT ) :: RJ( :,:,:,: )
!        END SUBROUTINE PHOT
         SUBROUTINE PA_IRR_CKBLK ( NUMCELLS, LIRRBLK, OFFSET,
     &                             CCOL, CROW, CLEV, NORDCELL,
     &                             NIRRCLS, IRRCELL )
            IMPLICIT NONE
            INTEGER, INTENT( IN )   :: NUMCELLS, OFFSET
            LOGICAL, INTENT( OUT )  :: LIRRBLK
            INTEGER, INTENT( IN )   :: CCOL( : ), CROW( : ), CLEV( : )
            INTEGER, INTENT( IN )   :: NORDCELL( : )
            INTEGER, INTENT( IN )   :: NIRRCLS
            INTEGER, INTENT( OUT )  :: IRRCELL( : )
         END SUBROUTINE PA_IRR_CKBLK
         SUBROUTINE PA_IRR_BLKENDC ( OFFSET, 
     &                               CCOL, CROW, CLEV, NORDCELL,
     &                               NIRRCLS, IRRCELL )
            IMPLICIT NONE
            INTEGER, INTENT( IN )   :: OFFSET
            INTEGER, INTENT( IN )   :: CCOL( : ), CROW( : ), CLEV( : )
            INTEGER, INTENT( IN )   :: NORDCELL( : )
            INTEGER, INTENT( IN )   :: NIRRCLS
            INTEGER, INTENT( IN )   :: IRRCELL( : )
         END SUBROUTINE PA_IRR_BLKENDC
         SUBROUTINE RBSOLVER ( LOGDEV, JDATE, JTIME, CHEMSTEP, NCSP,
     &                         LIRRFLAG, NIRRCLS, IRRCELL )
!           USE HGRD_DEFN
            USE GRID_CONF
            IMPLICIT NONE
            INTEGER, INTENT( IN )   :: LOGDEV, JDATE, JTIME
            REAL, INTENT( IN )      :: CHEMSTEP
            INTEGER, INTENT( IN )   :: NCSP
            LOGICAL, INTENT( IN )   :: LIRRFLAG
            INTEGER, INTENT( IN )   :: NIRRCLS
            INTEGER, INTENT( IN )   :: IRRCELL( BLKSIZE )
         END SUBROUTINE RBSOLVER
      END INTERFACE

C**********************************************************************

      IF ( N_GC_SPC .EQ. 0 ) RETURN

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  On first call, call routines to set-up for Gear solver and 
c  set-up to do emissions here if that option is invoked
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         LOGDEV = INIT3( )

         IF ( .NOT. CELLVAR_ALLOC() ) THEN
            MSG = 'Failure allocating variables dependent on horizontal extents'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
         ENDIF

         CALL RBINIT
         CALL RBSPARSE( LOGDEV )
         NOXYZ = MY_NCOLS * MY_NROWS * NLAYS
         MAOMV =  MWAIR / MWWAT


C...Initialize and report data

         WRITE( LOGDEV, 92020 ) NOXYZ, BLKSIZE, NBLKS, BLKLEN( 1 ), 
     &                          BLKLEN( NBLKS )

         WRITE( LOGDEV, 92040 ) GLBL_RTOL, GLBL_ATOL

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

c..Open file for solver stats if requested
#ifdef rbstats
         ALLOCATE( VSTATS( 3 ) )
         VSTATS( 1 ) = 'N_STRT_FAILS'
         VSTATS( 2 ) = 'N_FAILS'
         VSTATS( 3 ) = 'N_STEPS'

         IF( MYPE .EQ. 0 ) THEN

            IF( .NOT. OPEN3( CTM_CONC_1, FSREAD3, PNAME ) ) THEN
               MSG = 'Could not open ' // CTM_CONC_1 // ' file for readonly'
               CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT2 )
            ENDIF

            IF( .NOT. DESC3( CTM_CONC_1 ) ) THEN
               MSG = 'Could not get description of concentration file ' //
     &             CTM_CONC_1
               CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT2 )
            ENDIF

            EDATE = JDATE
            ETIME = JTIME
            CALL NEXTIME( EDATE, ETIME, TSTEP( 1 ) )

            SDATE3D = EDATE
            STIME3D = ETIME
            NVARS3D = 3
            NCOLS3D = GL_NCOLS
            NROWS3D = GL_NROWS
            NLAYS3D = NLAYS
            VNAME3D( 1 ) = 'N_STRT_FAILS'
            VNAME3D( 2 ) = 'N_FAILS'
            VNAME3D( 3 ) = 'N_STEPS'
            VDESC3D( 1 ) = 'Number of fails at start'
            VDESC3D( 2 ) = 'Number of step fails'
            VDESC3D( 3 ) = 'Number of steps'
            UNITS3D( 1 ) = 'None'
            UNITS3D( 2 ) = 'None'
            UNITS3D( 3 ) = 'None'
            VTYPE3D( 1 ) = M3REAL
            VTYPE3D( 2 ) = M3REAL
            VTYPE3D( 3 ) = M3REAL
            IF ( .NOT. OPEN3( CTM_RBSTATS_1, FSNEW3, PNAME ) ) THEN
               MSG = 'Could not create '// CTM_RBSTATS_1 // ' file'
               CALL M3EXIT( PNAME, SDATE3D, STIME3D, MSG, XSTAT2 )
            ENDIF

         ENDIF

         ALLOCATE( STAT_SUM( MY_NCOLS, MY_NROWS, NLAYS, 3 ) )

         STAT_SUM = 0

#endif

C Determine whether DEGRADE rountines are needed.

         IF( INDEX( MECHNAME, 'TX' ) .NE. 0 ) THEN
            CALL_DEG = .TRUE.
            WRITE( LOGDEV, * ) 'TX DEGRADE ROUTINES USED'
            WRITE( LOGDEV, * ) 'MECHNAME contains TX substring'
         ENDIF

C setup degradation array

         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
         ALLOCATE( Y_DEGRADE( BLKSIZE, NSPCSD ) )

      ENDIF      ! First call

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Start of integration driver after first call
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NIRRCLS = 0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Set date and time to center of time step, get necessary physical 
C  data, and get photolysis rates
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MIDDATE = JDATE
      MIDTIME = JTIME
      ITMSTEP = TIME2SEC( TSTEP( 2 ) )
      CHEMSTEP = FLOAT( ITMSTEP ) / 60.0
      CALL NEXTIME( MIDDATE, MIDTIME, SEC2TIME( ITMSTEP / 2 ) )
 
C.. Get ambient temperature in K

      VNAME = 'TA' 
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, TA ) ) THEN
         MSG = 'Could not read TA from MET_CRO_3D'
         CALL M3EXIT( PNAME, MIDDATE, MIDTIME, MSG, XSTAT1 )
      ENDIF
      
C.. Get specific humidity in Kg H2O / Kg air

      VNAME = 'QV'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, QV ) ) THEN
         MSG = 'Could not read QV from MET_CRO_3D'
         CALL M3EXIT( PNAME, MIDDATE, MIDTIME, MSG, XSTAT1 )
      ENDIF 
      
C.. Get pressure in Pascals

      VNAME = 'PRES'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, PRES ) ) THEN
         MSG = 'Could not read PRES from MET_CRO_3D'
         CALL M3EXIT ( PNAME, MIDDATE, MIDTIME, MSG, XSTAT1 )
      ENDIF

C.. Get density in Kg/m^3

      VNAME = 'DENS'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, DENS ) ) THEN
         MSG = 'Could not read DENS from MET_CRO_3D'
         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF

C.. Get photolysis rates in /min

      CALL PHOT ( MIDDATE, MIDTIME, JDATE, JTIME, NDARK, RJ )                    

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set flag for reordering of cells and put cells in sequential  
c  order initially
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LORDERING = .TRUE.
      IF( .NOT. LREORDER .OR. NBLKS .EQ. 1 ) LORDERING = .FALSE.
      DO NCELL = 1, NOXYZ
         NORDCELL( NCELL ) = NCELL
      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Top of loop over blocks. This loop will be done once if
C  no reordering, twice if reordering is required
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  100 CONTINUE

      DO 500 BLK = 1, NBLKS
         BLKID = BLK
         NUMCELLS = BLKLEN( BLK )
         OFFSET = BLKCNO( BLK )
         IF( .NOT. LORDERING .AND. LIRR ) THEN
             LIRRBLK = .FALSE.
             CALL PA_IRR_CKBLK ( NUMCELLS, LIRRBLK, OFFSET,
     &                           CCOL, CROW, CLEV, NORDCELL, NIRRCLS,
     &                           IRRCELL )
             CALL PA_IRR_BLKSTRT( NUMCELLS )
         ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Put the grid cell physical data in the block arrays, converting
C  pressure to atmospheres, water vapor to ppm, emissions to ppm/min 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         DO NCELL = 1, NUMCELLS
            CELLNUM = NORDCELL( OFFSET + NCELL )
            COL = CCOL( CELLNUM )
            ROW = CROW( CELLNUM )
            LEV = CLEV( CELLNUM )
            BLKTEMP( NCELL )  = TA( COL, ROW, LEV )
            BLKDENS( NCELL )  = DENS( COL, ROW, LEV )
            BLKPRES( NCELL )  = PA2ATM * PRES( COL, ROW, LEV )
            BLKCH2O( NCELL )  = MAX ( QV( COL, ROW, LEV ) * MAOMV * 
     &                          CONCOFM, 0.0 )          
         ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Put the grid cell concentrations in the block arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO ISP = 1, ISCHANG( NCS )
            ISPNEW = IOLD2NEW( ISP, NCS )
            DO NCELL = 1, NUMCELLS
               CELLNUM = NORDCELL( OFFSET + NCELL )
               COL = CCOL( CELLNUM )
               ROW = CROW( CELLNUM )
               LEV = CLEV( CELLNUM )
               Y( NCELL, ISPNEW ) = MAX( CONC( COL, ROW, LEV, ISP ), CONCMIN )
            ENDDO                 
         ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C   Get photolytic and thermal rate constants & call solver
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         LSUNLIGHT = .FALSE.

         DO NCELL = 1, NUMCELLS
            CELLNUM = NORDCELL( OFFSET + NCELL )
            COL = CCOL( CELLNUM )
            ROW = CROW( CELLNUM )
            LEV = CLEV( CELLNUM )
            DO NPH = 1, NPHOTAB
               RJBLK( NCELL, NPH ) = RJ( COL, ROW, LEV, NPH )
               IF( RJBLK( NCELL, NPH ) .GT. 0.0 ) LSUNLIGHT = .TRUE.
            ENDDO                         
         ENDDO
         
         CALL CALCKS( NPHOTAB, RJBLK )

         IF( LSUNLIGHT ) THEN
            NCSP = NCS
         ELSE
            NCSP = NCS + 1
         ENDIF

C  update cell concentrations for degradation routines

         IF( CALL_DEG ) THEN

            Y_DEGRADE = 0.0
            DO ISP = 1, NSPCSD
               DO NCELL = 1, NUMCELLS
                  CELLNUM = NORDCELL( OFFSET + NCELL )
                  COL = CCOL( CELLNUM )
                  ROW = CROW( CELLNUM )
                  LEV = CLEV( CELLNUM )
                  Y_DEGRADE( NCELL, ISP ) = MAX( CONC( COL, ROW, LEV, ISP ),
     &                                           CONCMIN )
               ENDDO
            ENDDO

C initialize degradation routines

            CALL INIT_DEGRADE( Y_DEGRADE, BLKTEMP, BLKDENS, RJBLK,
     &                         NUMCELLS,   JDATE,   JTIME, LOGDEV, BLKID )

         ENDIF

#ifdef rbstats

         NSTEPS = 0
         NFAILS = 0
         N_BAD_STARTS = 0

#endif


         CALL RBSOLVER(  LOGDEV, JDATE, JTIME, CHEMSTEP, NCSP,
     &                   LIRRBLK, NIRRCLS, IRRCELL ) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  If not ordering cells, save performance statistics and
C  store updated concentrations.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF( .NOT. LORDERING ) THEN

C..Update concentrations
            DO ISP = 1, ISCHANG( NCS )
               ISPOLD  = INEW2OLD( ISP, NCS )
               DO NCELL = 1, NUMCELLS
                  CELLNUM = NORDCELL( OFFSET + NCELL )
                  ROW = CROW( CELLNUM )
                  COL = CCOL( CELLNUM )
                  LEV = CLEV( CELLNUM )
                  CONC( COL, ROW, LEV, ISPOLD ) = Y( NCELL, ISP ) 
               ENDDO
            ENDDO

            IF( CALL_DEG ) THEN

C  Update degradation array with species treated by Rosenbach solver

               DO ISP = 1, ISCHANG( NCS )
                  ISPOLD  = INEW2OLD( ISP, NCS )
                  DO NCELL = 1, NUMCELLS
                     Y_DEGRADE( NCELL, ISPOLD ) = Y( NCELL, ISP )
                  ENDDO
               ENDDO

C  Update CGRID based on the degradation routines

               CALL FINAL_DEGRADE( Y_DEGRADE )
               DO ISP = NR_STRT, (TR_STRT-1)
                  DO NCELL = 1, NUMCELLS
                     CELLNUM = NORDCELL( OFFSET + NCELL )
                     ROW = CROW( CELLNUM )
                     COL = CCOL( CELLNUM )
                     LEV = CLEV( CELLNUM )
                     CONC( COL, ROW, LEV, ISP ) = Y_DEGRADE( NCELL, ISP )
                 ENDDO
               ENDDO
           ENDIF            !WTH

#ifdef rbstats

            DO NCELL = 1, NUMCELLS
               CELLNUM = NORDCELL( OFFSET + NCELL )
               ROW = CROW( CELLNUM )
               COL = CCOL( CELLNUM )
               LEV = CLEV( CELLNUM )
               STAT_SUM( COL, ROW, LEV, 1 ) = STAT_SUM( COL, ROW, LEV, 1 ) +
     &            N_BAD_STARTS
               STAT_SUM( COL, ROW, LEV, 2 ) = STAT_SUM( COL, ROW, LEV, 2 ) +
     &            NFAILS
               STAT_SUM( COL, ROW, LEV, 3 ) = STAT_SUM( COL, ROW, LEV, 3 ) +
     &            NSTEPS
            ENDDO


#endif
                      
            IF( LIRRBLK ) CALL PA_IRR_BLKENDC ( OFFSET, CCOL, CROW, CLEV,
     &                                          NORDCELL, NIRRCLS, IRRCELL )

         ENDIF

  500 CONTINUE

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  End of block loop; reorder cells if necessary and go back do the  
C  block loop again.  Taken from Jacobson 1994.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF( LORDERING ) THEN
         LORDERING = .FALSE.     
         LVAL = NOXYZ * 0.5 + 1
         IRVAL = NOXYZ
  600    IF( LVAL .GT. 1 ) THEN
            LVAL = LVAL - 1
            VALLOW = ERRMX2( LVAL )
            JREORD = NORDCELL( LVAL )
         ELSE
            VALLOW = ERRMX2( IRVAL )
            JREORD = NORDCELL( IRVAL )
            ERRMX2( IRVAL ) = ERRMX2( 1 )
            NORDCELL( IRVAL ) = NORDCELL( 1 )
            IRVAL = IRVAL - 1
            IF( IRVAL.EQ.1 ) THEN
               ERRMX2( IRVAL ) = VALLOW
               NORDCELL( IRVAL ) = JREORD
               GO TO 100
            ENDIF
         ENDIF
         IPAR = LVAL
         JPAR = LVAL + LVAL
  650    IF( JPAR .LE. IRVAL ) THEN
            IF( JPAR .LT. IRVAL ) THEN
               IF( ERRMX2( JPAR ) .LT. ERRMX2( JPAR + 1 ) )
     &                 JPAR = JPAR + 1
            ENDIF
            IF( VALLOW .LT. ERRMX2( JPAR )) THEN
               ERRMX2( IPAR ) = ERRMX2( JPAR )
               NORDCELL( IPAR ) = NORDCELL( JPAR )
               IPAR = JPAR
               JPAR = JPAR + JPAR
            ELSE
               JPAR = IRVAL + 1
            ENDIF
            GO TO 650
         ENDIF
         ERRMX2( IPAR ) = VALLOW
         NORDCELL( IPAR ) = JREORD
         GO TO 600
      ENDIF
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Output performance statistics if required and return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef rbstats

      WSTEP = WSTEP + TIME2SEC( TSTEP( 2 ) )
      EDATE = JDATE
      ETIME = JTIME
      CALL NEXTIME( EDATE, ETIME, TSTEP( 2 ) )
      IF( WSTEP .GE. TIME2SEC( TSTEP( 1 ) ) ) THEN

         WSTEP = 0

         DO S = 1, 3
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  DO L = 1, NLAYS
                     STATOUT( C, R, L ) = INT( STAT_SUM( C, R, L, S ) +
     &                       0.00001 )
                  ENDDO
               ENDDO
            ENDDO

            IF( .NOT. WRITE3( CTM_RBSTATS_1, VSTATS( S ),
     &            EDATE, ETIME, STATOUT ) ) THEN
               XMSG = 'Could not write ' // CTM_RBSTATS_1 // ' file'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            ENDIF  

         ENDDO

         STAT_SUM = 0.0

      ENDIF

#endif

      RETURN
      
C*********************** FORMAT STATEMENTS ****************************
92020 FORMAT( / 10X, 'Chemistry Solver Blocking Parameters ... ',
     &        / 10X, 'Domain Size (number of cells):             ', I10
     &        / 10X, 'Dimensioning Block Size (number of cells): ', I10
     &        / 10X, 'Number of Blocks:        ', I10
     &        / 10X, 'Size of General Blocks:  ', I10
     &        / 10X, 'Size of Last Block:      ', I10 )
92040 FORMAT( / 10X, 'Rosenbrock Chemistry Solver Error Control ',
     &               'Parameters ...',
     &        / 10X, 'RTOL : ', 1PE12.3,
     &        / 10X, 'ATOL : ', 1PE12.3, ' ppm' )

      END
                            

