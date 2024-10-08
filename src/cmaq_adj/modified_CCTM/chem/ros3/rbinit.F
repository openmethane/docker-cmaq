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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ros3/rbinit.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE RBINIT
C***********************************************************************
C
C  FUNCTION: To initialize species tolerances, allocate arrays, and 
C            define grid structure
C
C  PRECONDITIONS: None
C
C  RETURN VALUES: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, August, 2004
C                    31 Jan 05 J.Young: dyn alloc - establish both horizontal
C                    & vertical domain specifications in one module (GRID_CONF)
C                    Get BLKSIZE from module GRID_CONF
C                    29 Jul 05 WTH: allocate variables used by degrade
C                                   routines.
C                      
C***********************************************************************
      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE RBDATA                   ! Rosenbrock solver data

      IMPLICIT NONE

C.....INCLUDES:
      INCLUDE SUBST_IOPARMS        ! Io/api parameters
      INCLUDE SUBST_IOFDESC        ! Io/api file descriptions
      INCLUDE SUBST_IODECL         ! Io/api declarations
!     INCLUDE SUBST_VGRD_ID        ! Vertical grid data
      INCLUDE SUBST_GC_SPC         ! Gas chem species names and MWs
      INCLUDE SUBST_RXCMMN         ! Mechanism reaction common block

C.....ARGUMENTS: NONE
 

C.....PARAMETERS:

      CHARACTER( 16 ), PARAMETER   :: PNAME = 'RBINIT'    ! Program name
      REAL, PARAMETER              :: DEF_ATOL = 1.0E-07  ! Default ATOL
      REAL, PARAMETER              :: DEF_RTOL = 1.0E-03  ! Default RTOL


C.....EXTERNAL FUNCTIONS:
      REAL ENVREAL                 ! Gets real value of environment variable


C.....SAVED VARIABLES: NONE

 
C.....LOCAL VARIABLES: 
      CHARACTER( 132 ) :: MSG      ! Log message
      CHARACTER(  80 ) :: VARDESC  ! Description of environment variable 

      INTEGER N                    ! Loop index

      INTEGER COL                  ! Column number index
      INTEGER IAVGSIZE             ! Average number of cells per block
      INTEGER IDUMMY               ! Dummy integer variable 
      INTEGER IERRST               ! Ouput error status number
      INTEGER LEV                  ! Level number index
      INTEGER OFFSET               ! Pointer for start cell number in a block
      INTEGER NBLK                 ! Block number index
      INTEGER NCOUNT               ! Counter for number of cells for grid
      INTEGER NOXYZ                ! Total number of cells for grid
      INTEGER ROW                  ! Row number index
      INTEGER STATUS               ! Status code for functions

      REAL    DEFTOL               ! Default tolerance value


C***********************************************************************


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..Initialize vars & allocate arrays used in sparse matrix treatment
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      LREORDER = .TRUE.

      N_SPEC = N_GC_SPC
      N_RXNS = NRXNS

      MXRR = 3 * MXRCT
      MXRP = 3 * MXPRD

      MXCOUNT1 = N_GC_SPC * MAXGL3 * 3
      MXCOUNT2 = N_GC_SPC * MAXGL3 * 3

      ALLOCATE( NKUSERAT( NRXNS, NCS2 ) )
      ALLOCATE( NDERIVL(  NRXNS, NCS2 ) )
      ALLOCATE( NDERIVP(  NRXNS, NCS2 ) )
      ALLOCATE( IRM2( NRXNS, MXRCT + MXPRD, NCS2 ) )
      ALLOCATE( ICOEFF( NRXNS, MXRP, NCS2 ) )

      ALLOCATE( JARRAYPT( N_GC_SPC, N_GC_SPC, NCS2 ) )
      ALLOCATE( JARRL( NRXNS, MXRR, NCS2 ) )
      ALLOCATE( JARRP( NRXNS, MXRP, NCS2 ) )
      ALLOCATE( JLIAL( NRXNS, MXRR, NCS2 ) )
      ALLOCATE( JPIAL( NRXNS, MXRP, NCS2 ) )

      ALLOCATE( INEW2OLD( N_GC_SPC,  NCS ) )
      ALLOCATE( IOLD2NEW( N_GC_SPC,  NCS ) )

      ALLOCATE( JZEROA( MXARRAY ) )
      ALLOCATE( JZEROB( MXARRAY ) )

      ALLOCATE( JZLO( NCS2 ) )
      ALLOCATE( IDEC1LO(  N_GC_SPC, NCS2 ) )
      ALLOCATE( IDEC1HI(  N_GC_SPC, NCS2 ) )

      ALLOCATE( IJDECA( MXCOUNT2 ) )
      ALLOCATE( IJDECB( MXCOUNT2 ) )
      ALLOCATE( IKDECA( MXCOUNT2 ) )
      ALLOCATE( IKDECB( MXCOUNT2 ) )
      ALLOCATE( KJDECA( MXCOUNT2 ) )
      ALLOCATE( KJDECB( MXCOUNT2 ) )

      ALLOCATE( JHIZ1( N_GC_SPC, NCS2 ) )
      ALLOCATE( JHIZ2( N_GC_SPC, NCS2 ) )
      ALLOCATE( KZLO1( N_GC_SPC, NCS2 ) )
      ALLOCATE( KZLO2( N_GC_SPC, NCS2 ) )
      ALLOCATE( KZHI0( N_GC_SPC, NCS2 ) )
      ALLOCATE( KZHI1( N_GC_SPC, NCS2 ) )
      ALLOCATE( KZERO( MXARRAY,  NCS2 ) )
      ALLOCATE( KZILCH( N_GC_SPC, NCS2 ) )
      ALLOCATE( MZHI0(  N_GC_SPC, NCS2 ) )
      ALLOCATE( MZHI1(  N_GC_SPC, NCS2 ) )
      ALLOCATE( MZILCH( N_GC_SPC, NCS2 ) )
      ALLOCATE( MZLO1(  N_GC_SPC, NCS2 ) )
      ALLOCATE( MZLO2(  N_GC_SPC, NCS2 ) )

      ALLOCATE( VDIAG( BLKSIZE, N_GC_SPC ) )

      ALLOCATE( CC2( BLKSIZE, 0 : MXARRAY ) )

c..cell & solver data
      ALLOCATE( BLKTEMP( BLKSIZE ) )        
      ALLOCATE( BLKPRES( BLKSIZE ) )        
      ALLOCATE( BLKCH2O( BLKSIZE ) )       
      ALLOCATE( BLKDENS( BLKSIZE ) )

      ALLOCATE( RKI( BLKSIZE, N_RXNS ) )    
      ALLOCATE( ATOL( N_SPEC ) )            
      ALLOCATE( RTOL( N_SPEC ) )            
      ALLOCATE( Y(   BLKSIZE, N_SPEC ) )    


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set convergence tolerances for each species; currently uses
c  one set of tolerances for all species
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      VARDESC = 'Relative tolerance for Rosenbrock Solver'
      DEFTOL = DEF_RTOL
      GLBL_RTOL = ENVREAL( 'RB_RTOL', VARDESC, DEFTOL, STATUS )
      IF ( STATUS .EQ. 1 ) THEN
         MSG = 'Environment variable RB_RTOL improperly formatted'
         CALL M3EXIT( PNAME, 0, 0, MSG, 0 )
      ELSEIF ( STATUS .LT. 0 ) THEN
         GLBL_RTOL = DEF_RTOL
      ENDIF

      VARDESC = 'Absolute tolerance for Rosenbrock Solver'
      DEFTOL = DEF_ATOL
      GLBL_ATOL = ENVREAL( 'RB_ATOL', VARDESC, DEFTOL, STATUS )
      IF ( STATUS .EQ. 1 ) THEN
         MSG = 'Environment variable RB_ATOL improperly formatted'
         CALL M3EXIT( PNAME, 0, 0, MSG, 0 )
      ELSEIF ( STATUS .LT. 0 ) THEN
         GLBL_ATOL = DEF_ATOL
      ENDIF
      

      DO N = 1, N_GC_SPC
         ATOL( N ) = GLBL_ATOL
         RTOL( N ) = GLBL_RTOL
      ENDDO


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get number of cells in grid and store i,j,k indices of cells in
c  sequential order
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NCOUNT = 0
      DO 120 LEV = 1, NLAYS
         DO 100 COL = 1, MY_NCOLS
            DO 80 ROW = 1, MY_NROWS
               NCOUNT = NCOUNT + 1
               CROW( NCOUNT ) = ROW
               CCOL( NCOUNT ) = COL
               CLEV( NCOUNT ) = LEV
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Define block structure for grid; stop if maxblks exceeded
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NOXYZ = NCOUNT
      NBLKS = 1 + ( NOXYZ - 1 ) / BLKSIZE
      IF( NBLKS .GT. MXBLKS ) THEN
         WRITE( MSG, 92020 ) NBLKS, MXBLKS
         IDUMMY = 0
         IERRST = 2
         CALL M3EXIT( PNAME, IDUMMY, IDUMMY, MSG, IERRST )
      ENDIF

      IAVGSIZE = 1 + ( NOXYZ - 1 ) / NBLKS
      IAVGSIZE = MIN( IAVGSIZE, BLKSIZE )
      OFFSET = 0

      DO 200 NBLK = 1, NBLKS
         IF( NBLK .LT. NBLKS ) THEN
            BLKCNO( NBLK ) = OFFSET
            BLKLEN( NBLK ) = IAVGSIZE
            OFFSET = OFFSET + IAVGSIZE
         ELSE
            BLKCNO( NBLK ) = OFFSET
            BLKLEN( NBLK ) = NOXYZ - ( ( NBLK-1 ) * IAVGSIZE )
         ENDIF
  200 CONTINUE

      RETURN
      
C********************** FORMAT STATEMENTS ******************************      
92020 FORMAT( 1X, 'ERROR: MAXIMUM NUMBER OF BLOCKS EXCEEDED',
     &            ' FOR GRID', 'NBLKS=', I3, 1X, ' MAXBLKS=',
     &            I3, '  CHANGE GRPARMS.EXT' )

      END
