
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhg/rdemis_pm_hg.F,v 1.1.1.1 2010/06/14 16:03:09 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C ///////////////////////////////////////////////////////////////////////////
C  SUBROUTINE RDEMIS_PM_HG reads aerosol emissions from gridded input file and
C   converts into molar-mixing-ratio units, as required by the vertical 
C   diffusion routines
C
C  DEPENDENT UPON:  NONE
C
C  REVISION HISTORY:
C
C     17 Sept 2009 B.Hutzell: modified routine to allow in-line plume rise
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE RDEMIS_PM_HG ( JDATE, JTIME, EMISLYRS, NSPC_EMIS, RJACM, VDEMIS_AE )

      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE PT3D_DEFN, ONLY: PT3DEMIS, PMEMIS_PT ! from plume rise
      USE PMEM_DEFN                            ! defines PM emissions species
      USE PTMAP, ONLY: N_SPC_PTPM, PTPM_MAP, STKPM_NAME    ! defines emissions file(s) PM mapping
       
      IMPLICIT NONE
      
     
      INCLUDE SUBST_RXCMMN    ! to get mech name
      INCLUDE SUBST_AE_SPC    ! aerosol species names and molecular weights
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names and map table
      INCLUDE SUBST_CONST     ! physical and mathematical constants
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include      SUBST_IODECL      # I/O definitions and declarations
!        INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOFDESC   ! file header data structure


C *** ARGUMENTS
     
       INTEGER JDATE           ! current model date, coded YYYYDDD
       INTEGER JTIME           ! current model time, coded HHMMSS
C      INTEGER TSTEP( 2 )      ! time step vector (HHMMSS)
C Note:  TSTEP(1) = local output step
C Note:  TSTEP(2) = sciproc sync. step (chem)

      INTEGER EMISLYRS        ! number of vertical layers on emissions file
      INTEGER NSPC_EMIS
      
      REAL    RJACM( :,:,: ) ! recip of mid-layer Jacobian [1/m]
      REAL :: VDEMIS_AE( :,:,:,: )       ! aerosol emissions 

C *** LOCAL VARIABLES
      CHARACTER( 16 ), SAVE :: PNAME = 'RDEMIS_PM_HG'
C     Aerosol version name
      CHARACTER( 16 ), SAVE :: AE_VRSN
C     States whether the mechanism is for Mercuric Aerosols
      LOGICAL,         SAVE :: HG_VRSN


      REAL  DX1, DX2                               ! grid-cell width and length [ m ]
      REAL  GRDAREA                                ! grid area [m2]
      REAL, ALLOCATABLE, SAVE :: GRDHGT( : )       ! grid height [sigma]
      REAL, ALLOCATABLE, SAVE :: GRDVOL( : )       ! grid volume [m2*sigma]
      REAL, ALLOCATABLE, SAVE :: CONVEM_AE( : )    ! into [ug/m2/sec]


C Array dimensions
      INTEGER, PARAMETER :: N_EMIS_SPECIES = 1  ! maximum number of PM species
                                                ! in emission input file and point emissions
C PM variables in emissions file
      CHARACTER( 16 ), SAVE :: EMIS_SPECIES( N_EMIS_SPECIES )   ! in emission input files
      INTEGER, SAVE         :: SEMIS_MAP(     N_EMIS_SPECIES )   ! indices emissions 
      INTEGER, SAVE         :: PTEMIS_MAP(   N_EMIS_SPECIES )   ! indices emissions 
      INTEGER, SAVE         :: PTEMIS_INDX(  N_EMIS_SPECIES )   ! point emissions species
      CHARACTER( 16 )       :: UNITSAE(      N_EMIS_SPECIES )   ! units of emissions

      DATA EMIS_SPECIES( 1 ) / 'PHGI'    /

C     Emission rate of all aerosol species interpolated to current time
       REAL, ALLOCATABLE, SAVE :: EMBUFF(:, :,:,: )  ! in all grid cells
       REAL, ALLOCATABLE, SAVE :: PM_EM (:)           ! in one grid cell

C     Speciation factors for coarse mode emissions, from paragraph #15 of
C     Binkowski & Roselle (2003)
       REAL, PARAMETER :: FAC_DUST     = 0.90   ! (fugitive dust)/PMC
       REAL, PARAMETER :: FAC_OTHER    = 0.10  ! (non-fugitive dust)/PMC

C     Factors for splitting primary carbon emissions into Aitken and
C     accumulation modes, from paragraph #12 of Binkowski & Roselle (2003)
         REAL, PARAMETER :: FACEM25_ACC  = 0.999  ! accumulation mode
         REAL, PARAMETER :: FACEM25_ATKN = 0.001 ! Aitken mode

C     Number and names of toxic PM species and surrogates in emissions file
         INTEGER,       PARAMETER    :: N_AETOX_SPC = 2
         CHARACTER(16),      SAVE    :: AE_TOX_SPC(   N_AETOX_SPC )
         CHARACTER(16),     SAVE     :: AE_TOX_SUR( N_AETOX_SPC )
         REAL,              SAVE     :: TOX_EMIS_FAC(  N_AETOX_SPC )

         DATA AE_TOX_SPC(  1 ), AE_TOX_SUR(  1 ) / 'APHGI', 'PHGI' /
         DATA AE_TOX_SPC(  2 ), AE_TOX_SUR(  2 ) / 'APHGJ', 'PHGI' /
C         DATA AE_TOX_SPC(  3 ), AE_TOX_SUR(  3 ) / 'APHGK', 'PHGK' /
         
         DATA TOX_EMIS_FAC(  1) /      FACEM25_ATKN / 
         DATA TOX_EMIS_FAC(  2) /      FACEM25_ACC  /
C         DATA TOX_EMIS_FAC(  3) /      1.0 /

C local indicies, molecular weights, and emission pointers of toxic species
         INTEGER,       ALLOCATABLE, SAVE ::     N_AETOX( : )
         REAL,          ALLOCATABLE, SAVE ::   AE_TOX_MW( : )
         INTEGER      , ALLOCATABLE, SAVE ::  TOX_EMIS_MAP( : )
 
C     Descriptive variables from the emission input file
C         CHARACTER( 16 ) :: UNITSCK                   ! units of 1st aero species
         INTEGER, SAVE   :: NAESPCEMIS                ! number of input aero spcs

C     Factor for converting aerosol emissions from input units ...
         REAL  CONVEM_AE_MASS                       ! into [ug/sec]
         REAL  GSFAC                                ! into [ug/m3/sec]
         REAL  CONVM     ! conversion factor for mass emissions [m3/mol]
         REAL  CONVN     ! conversion factor for number emissions [1e6*m3]
         REAL  CONVS     ! conversion factor for surface area emissions [m3/mol]

C     Domain decomposition info from emission and meteorology files
         INTEGER, SAVE ::   STARTCOL,    ENDCOL,   STARTROW, ENDROW
         INTEGER, SAVE :: STRTCOL_O1, ENDCOL_O1, STRTROW_O1, ENDROW_O1
         INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C     Variables interpolated from the meteorological input files
         REAL DENS( NCOLS,NROWS,NLAYS )        ! air density [ kg/m3 ]

         REAL EMIS_TOX                         ! diesel PM emissions

C     Geometric Constants
         REAL,      PARAMETER :: F6DPI   = 6.0 / PI
         REAL,      PARAMETER :: F6DPIM9 = 1.0E-9 * F6DPI

C     Variables for converting emission rates into molar-mixing-ratio units
         REAL,      PARAMETER :: GPKG  = 1.0E+03     ! g/kg
         REAL,      PARAMETER :: MGPG  = 1.0E+06     ! ug/g
         REAL,      PARAMETER :: RAVO  = 1.0 / AVO   ! reciprocal of Avogadro
!        REAL,      PARAMETER :: MWAIR = 28.9628     ! molar mass of dry air [g/mol]
                                                     ! assuming 78.06% N2, 21% O2, 
                                                     ! and 0.943% Ar 
                                                     ! (Source : Hobbs, 1995) pp 69-70

C     Miscellaneous variables
         INTEGER, SAVE     :: LOGDEV
         INTEGER, EXTERNAL :: INDEX1
         INTEGER, EXTERNAL :: TRIMLEN

         CHARACTER( 96 )   :: XMSG = ' '
         CHARACTER( 16 )   :: VNAME       ! temp var for species names

         INTEGER  GXOFF, GYOFF            ! origin offset
         INTEGER C,R,L,N,V                ! Loop indices
         INTEGER INDX
         INTEGER ALLOCSTAT

         LOGICAL, EXTERNAL :: ENVYN       ! get environment variable as boolean
         CHARACTER( 80 ) :: VARDESC       ! environment variable description
         INTEGER  STATUS                  ! ENV... status

         LOGICAL, SAVE :: FIRSTIME = .TRUE.

C ----------------------------------------------------------------------

C        VDEMIS_AE = 0.0

         IF ( FIRSTIME ) THEN
 
   
            FIRSTIME = .FALSE.
            LOGDEV = INIT3()


   
            IF ( INDEX ( MECHNAME, 'AE3' ) .GT. 0 ) THEN
               AE_VRSN  = 'AE3'
            ELSE IF ( INDEX ( MECHNAME, 'AE4' ) .GT. 0 ) THEN
               AE_VRSN  = 'AE4'
            ELSE IF ( INDEX ( MECHNAME, 'AE5' ) .GT. 0 ) THEN
              AE_VRSN  = 'AE5'
            ELSE
               XMSG = 'This version of the emission processing code '
     &             // 'can only be used with the AE3, AE4, AE5 aerosol '
     &             // 'mechanisms.'
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF ! check on MECHNAME


            IF (INDEX(MECHNAME,'CB05TXHG') .GT. 0 
     &                 .OR. INDEX(MECHNAME,'CB05HG') .GT. 0) THEN
               HG_VRSN  = .TRUE.
               WRITE(LOGDEV,*)'Mechanism with mercuric PM species '
            ELSE 
               HG_VRSN  = .FALSE.
               RETURN
            END IF ! check for HG mechanism 

            ALLOCATE(     N_AETOX(N_AETOX_SPC) )
            ALLOCATE(   AE_TOX_MW(N_AETOX_SPC) )
            ALLOCATE(TOX_EMIS_MAP(N_AETOX_SPC) )

            N_AETOX = 0
            DO V = 1, N_AETOX_SPC
               VNAME = AE_TOX_SPC( V )
               N_AETOX(V) = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N_AETOX(V) .EQ. 0 ) THEN
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
            ENDDO


C *** Save array of MWs
C Debugging statements
             XMSG = 'Below species belong to the toxic AE group.'
C            WRITE(LOGDEV,*)XMSG
             DO V = 1, N_AETOX_SPC
                N = AE_EMIS_MAP( N_AETOX(V) )
                AE_TOX_MW( V ) = AE_MOLWT( N )
                INDX = N_AETOX( V )
C Debugging statements
D               WRITE(LOGDEV,'(3(I3,1X),A16,2(1X,E12.4))')V,N,
D    &                       INDX,AE_SPC( N ),AE_TOX_MW(V),
D    &                       AE_EMIS_FAC( INDX )
             END DO

C *** Open the gridded emissions file, which contains gas, aerosol, 
C     and non-reactive species

            IF ( .NOT. OPEN3( EMIS_1, FSREAD3, PNAME ) ) THEN
               XMSG = 'Could not open '// EMIS_1 // ' file'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            IF ( .NOT. DESC3( EMIS_1 ) ) THEN
               XMSG = 'Could not get '// EMIS_1 // ' file description'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

C *** Search emission file for certain species names.  For each species found, 
C     1. Populate the AEMIS array with the emitted species names
C        note: species names are hardcoded in the present version
C     2. Populate the UNITSAE array based on the units in which the 
C        emissions data are input (e.g., 'G/S', 'KG/H')
C
C     If a species is not found, print error message and halt.


            CALL M3MESG ( 'define aetox emissions pointers ' )

            SEMIS_MAP = 0
            DO V = 1, N_EMIS_SPECIES
               VNAME = EMIS_SPECIES( V )
               INDX = INDEX1( VNAME, N_EMIS_PM, PM_EMIS )
               IF ( INDX .NE. 0 ) THEN
                  SEMIS_MAP( V ) = INDX
                  UNITSAE( V ) = UNITSPM( INDX )
               ELSE
                  XMSG = 'Could not find '
     &                // VNAME( 1:LEN_TRIM( VNAME ) )
     &                // ' in PM_EMIS array'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
            ENDDO
D            DO V = 1, N_EMIS_SPECIES
D              WRITE(LOGDEV,1049)V,EMIS_SPECIES(V),V,SEMIS_MAP(V)
D            ENDDO
 
           PTEMIS_MAP = 0
           IF ( PT3DEMIS ) THEN
            DO V = 1, N_EMIS_SPECIES
               VNAME = EMIS_SPECIES( V )
               INDX = INDEX1( VNAME, N_SPC_PTPM, STKPM_NAME )
               IF ( INDX .NE. 0 ) THEN
                  PTEMIS_MAP( V ) = PTPM_MAP(INDX)
                  PTEMIS_INDX(V ) = INDX
                      XMSG = PNAME // ':Found '
     &                    // VNAME( 1:LEN_TRIM( VNAME ) )
     &                    // ' in STKPM_NAME array named as '
     &                    // STKPM_NAME(INDX)
                      CALL M3MESG ( XMSG )
               ELSE
                  XMSG = 'Could not find '
     &                // VNAME( 1:LEN_TRIM( VNAME ) )
     &                // ' in STKPM_NAME array'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG)
               END IF
            ENDDO
D            DO V = 1, N_SPC_PTPM
D               WRITE( LOGDEV,1043 ) V, STKPM_NAME( V ), PTPM_MAP( V )
D            ENDDO
D            DO V = 1, N_EMIS_SPECIES
D              INDX = PTEMIS_MAP( V )
D              IF(INDX .GT. 0)THEN
D                 WRITE(LOGDEV,1050)V,EMIS_SPECIES(V),V,PTEMIS_MAP( V )
D              ELSE
D                 WRITE(LOGDEV,1050)V,EMIS_SPECIES(V),INDX,INDX
D	      ENDIF
D            ENDDO
           ENDIF

           TOX_EMIS_MAP = 0
           DO V = 1, N_AETOX_SPC
              VNAME = AE_TOX_SUR( V )
              TOX_EMIS_MAP( V ) = INDEX1( VNAME, N_EMIS_SPECIES, EMIS_SPECIES )
               IF ( TOX_EMIS_MAP( V ) .EQ. 0 ) THEN
                  XMSG = 'Could not find ' // VNAME // 'in EMIS file'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
           ENDDO

            
C *** Save the number of aerosol species read from the emission file and
C     write their names to the log file

            NAESPCEMIS = N_EMIS_SPECIES
            VNAME      = ' UNITS '
            WRITE( LOGDEV, '( /5X, A )' ) 'Mercuric PM species in EMIS_1 file:'
            DO V = 1, NAESPCEMIS
C               WRITE( LOGDEV, '( /5X, 3(1X,A) )' ) EMIS_SPECIES( V ),VNAME,
C     &                                             UNITSAE( V )
               WRITE( LOGDEV, '( /5X, 3(1X,A) )' ) PM_EMIS( INDX ),VNAME,
     &                                             UNITSPM( INDX )
            END DO

C *** Allocate memory for PM_EM, EMBUFF, GRDHGT, GRDVOL, CONVEM_AE, and SSBF

            ALLOCATE ( PM_EM( NAESPCEMIS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** PM_EM memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            ALLOCATE ( EMBUFF( MY_NCOLS,MY_NROWS,EMISLYRS,NAESPCEMIS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** EMBUFF memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( CONVEM_AE( EMISLYRS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** CONVEM_AE memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            ALLOCATE ( GRDHGT( EMISLYRS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** GRDHGT memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            ALLOCATE ( GRDVOL( EMISLYRS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = '*** GRDVOL memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            

C *** Calculate length and width of each grid cell
C     note: crude estimate is made for LAT/LONG coordinate systems

            IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
               DX1 = DG2M * XCELL_GD ! in m
               DX2 = DG2M * YCELL_GD
     &             * COS( PI180*( YORIG_GD + YCELL_GD
     &             * FLOAT( GL_NROWS/2 ) ) ) ! in m
            ELSE
               DX1 = XCELL_GD        ! in m
               DX2 = YCELL_GD        ! in m
            END IF

C *** Calculate height of grid cell in each layer in sigma coordinates
C     Multiply by grid area [m2] to obtain grid volume

            GRDAREA = DX1 * DX2
            DO L = 1, EMISLYRS
               GRDHGT( L ) = X3FACE_GD( L ) - X3FACE_GD( L-1 )
               GRDVOL( L ) = GRDHGT( L ) * GRDAREA
            END DO

C *** Confirm that all aerosol species in the emission input file have 
C     the same units.  If not, print error message and halt.

            DO N = 1, NAESPCEMIS
               V = SEMIS_MAP( N )
               IF ( UNITSPM( V ) .NE. UNITSCKPM ) THEN
                     XMSG = 'PM Units not uniform in PM_EMIS species ' //
     &                      ' Required Units are '    // UNITSCKPM //
     &                      ' Emissions '   // PM_EMIS( V ) // 
     &                      ' has units '   // UNITSPM( V )
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            END DO

C *** Calculate scaling factor for converting aerosol emissions from
C     their input units to [ug/s] and then to [ug/m2/s] using layer-
C     specific grid-cell volume

            IF ( UNITSCKPM .EQ. 'G/S' .OR.
     &           UNITSCKPM .EQ. 'g/s' ) THEN
               CONVEM_AE_MASS = MGPG                  ! (g/s) -> (ug/s)
            ELSE IF ( UNITSCKPM .EQ. 'KG/HR' .OR.
     &                UNITSCKPM .EQ. 'kg/hr' ) THEN
               CONVEM_AE_MASS = GPKG * MGPG / 3600.0  ! (kg/hr) -> (ug/s)
            ELSE
               XMSG = 'Units incorrect on ' // EMIS_1
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

            DO L = 1, EMISLYRS
               CONVEM_AE( L ) = CONVEM_AE_MASS / GRDVOL( L )
            END DO

C *** Get domain decomposition info from input files

            CALL SUBHFILE ( EMIS_1, GXOFF, GYOFF,
     &                      STARTCOL, ENDCOL, STARTROW, ENDROW )

            CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                      STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

 
         END IF    ! FIRSTIME

C ----------------------------------------------------------------------

         IF( .NOT. HG_VRSN )RETURN

C *** Read aerosol emission rates from file and interpolate to the current
C     time.  Store result in EMBUFF array.


         EMBUFF = 0.0   ! clear array

         IF ( PT3DEMIS ) THEN

             DO N = 1, N_EMIS_SPECIES
               V = SEMIS_MAP( N )
D               WRITE(LOGDEV,*)'Reading surface emissions for ',PM_EMIS( V ),' at ',
D    &                         JDATE, JTIME
               IF ( .NOT. INTERPX( EMIS_1, PM_EMIS( V ), PNAME,
     &                             STARTCOL,ENDCOL, STARTROW,ENDROW, 1,1,
     &                             JDATE, JTIME, EMBUFF( 1,1,1,N ) ) ) THEN
                  XMSG = 'Could not read '
     &                 // PM_EMIS( N )( 1:TRIMLEN( PM_EMIS( N ) ) )
     &                 // ' from ' // EMIS_1
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG  )
               END IF
            END DO

!            DO N = 1, N_SPC_PTPM     ! merge inline point src PM species
             DO N = 1, N_EMIS_SPECIES 
               V = PTEMIS_MAP( N )
               IF(V .LT. 1)CYCLE
               INDX = PTEMIS_INDX( N )
D               XMSG =  'Adding point emissions for ' // STKPM_NAME(INDX)
D    &              // ' to ' // PM_EMIS( SEMIS_MAP( N ) ) 
D    &              // ' emissions'
D               CALL M3MESG( XMSG )
               DO L = 1, EMISLYRS
                  DO R = 1, MY_NROWS
                     DO C = 1, MY_NCOLS
                     EMBUFF( C,R,L,N ) = EMBUFF( C,R,L,N ) + PMEMIS_PT( C,R,L,V )
                     END DO
                  END DO
               END DO
            END DO


         ELSE

             DO N = 1, N_EMIS_SPECIES
                V = SEMIS_MAP( N )
D                WRITE(LOGDEV,*)'Reading emissions for ',PM_EMIS( V ),' at ',
D    &                         JDATE, JTIME
                IF ( .NOT. INTERPX( EMIS_1, PM_EMIS( V ), PNAME,
     &                             STARTCOL,ENDCOL, STARTROW,ENDROW, 1,EMISLYRS,
     &                             JDATE, JTIME, EMBUFF( 1,1,1,N ) ) ) THEN
                  XMSG = 'Could not read '
     &                // PM_EMIS( V )( 1:LEN_TRIM( PM_EMIS( V ) ) )
     &                // ' from ' // EMIS_1
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG  )
                END IF
             ENDDO

         END IF

C *** Read air density [ kg/m3 ] and store results.

         IF ( .NOT. INTERPX( MET_CRO_3D, 'DENS', PNAME,
     &                       STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                       1,NLAYS, JDATE, JTIME, DENS ) ) THEN
            XMSG = 'Could not interpolate DENS from ' // MET_CRO_3D
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF


C *** LOOP OVER EACH GRID CELL

         DO L = 1, EMISLYRS
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS

C *** Calculate scaling factor for converting mass emissions into [ ug/m3/s ]
C     note: RJACM converts grid heights from sigma coordinates to meters
C     Also calculate scaling factors for converting to molar-mixing-ratio units

                  GSFAC = CONVEM_AE( L )  * RJACM( C,R,L )
                  CONVM = MWAIR / GPKG / DENS( C,R,L )  !  [m3/mol]
                  CONVN = CONVM * RAVO * 1.0E+06        !  10^6 * [m3]
                  CONVS = CONVM                         !  [m3/mol]

C *** Convert emission rates ug/m3/s and store into a 1D array

                  DO N = 1, NAESPCEMIS
                     PM_EM( N ) = GSFAC*EMBUFF( C,R,L,N )
                  END DO

C Compute and store toxic PM emissions

                  DO N = 1, N_AETOX_SPC

                      EMIS_TOX = 0.0

                      V    = TOX_EMIS_MAP( N )

                      IF(V .LT. 1 .OR. V .GT. NAESPCEMIS)CYCLE
                      
                      INDX      = N_AETOX( N )
                      EMIS_TOX  = PM_EM( V )

                      VDEMIS_AE( INDX, L,C,R ) =  CONVM
     &                                         * EMIS_TOX
     &                                         * AE_EMIS_FAC( INDX ) 
     &                                         * TOX_EMIS_FAC( N )
     &                                         / AE_TOX_MW( N )

                  ENDDO


               END DO   ! loop on MY_NCOLS
            END DO   ! loop on MY_NROWS
         END DO   ! loop on EMISLYRS

C Debugging statements
D         write(logdev,*)'After RDEMIS_AETOX'
1043     FORMAT( 5X, I3, 1X,'STKPM_NAME, PTPM_MAP = ', A16, I3 )
1049     FORMAT(1X,'EMIS_SPECIES( ',I3,' ) = ',A16,' SEMIS_MAP( ',I3,' ) = ',I3)
1050     FORMAT(1X,'EMIS_SPECIES( ',I3,' ) = ',A16,'PTEMIS_MAP( ',I3,' ) = ',I3)
        

         RETURN
         
         END SUBROUTINE RDEMIS_PM_HG

      
