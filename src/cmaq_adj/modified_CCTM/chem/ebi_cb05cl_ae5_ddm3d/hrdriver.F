
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_cb05cl_ae5_ddm3d/hrdriver.F,v 1.1.1.1 2010/06/14 16:03:01 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

#ifdef sens
      SUBROUTINE CHEM ( CGRID, SENGRID, JDATE, JTIME,
     &                  TSTEP, NPMAX, IPT, IPARM, STDATE,
     &                  IDATE, ITIME, IREGION, IRXN )
#else
      SUBROUTINE CHEM ( CGRID, JDATE, JTIME, TSTEP )
#endif

C**********************************************************************
C
C  FUNCTION: Driver subroutine for Euler Backward Iterative solver
C
C  PRECONDITIONS: For the CB05CL family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  HRINIT, PHOT, HRCALCKS, HRSOLVER
C
C  REVISION HISTORY: Created by EBI solver program, Jan. 31, 2008
C                       Based on the algorithm in "Test of Two Numerical
C                       Schemes for Use in Atmospheric Transport-Chemistry
C                       Models", O. Hertel, R. Berkowicz, J. Christensen,
C                       and O. Hov, Atm Env., Vol. 27A, No. 16, 1993.
C                       Original MEBI code developed by Ho-Chun Huang,
C                       SUNY, Albany -- "On the performance of numerical
C                       solvers for a chemistry submodel in three-dimensional
C                       air quality models 1. Box model simulations",
C                       H. Huang and J.S. Chang, JGR, Vol 106, No. D17, 2001.
C                       This version replaces Huang and Chang use of numerical
C                       solutions with analytical solutions derived in
C                       Hertel et al.
C
C Aug 08    S.Napelenok:  DDM-3D implementation for CMAQv4.7 based on previous work by
C                         D.Cohan
C
C**********************************************************************

      USE HGRD_DEFN             ! horizontal domain specifications
      USE VGRD_DEFN             ! vertical layer specifications
      USE HRDATA

      IMPLICIT NONE

C..Includes:
      INCLUDE SUBST_IOPARMS   ! Io/api parameters
      INCLUDE SUBST_IOFDESC   ! Io/api file descriptions
      INCLUDE SUBST_IODECL    ! Io/api declarations
      INCLUDE SUBST_FILES_ID  ! CMAQ files
      INCLUDE SUBST_CONST     ! CMAQ constants
      INCLUDE SUBST_GC_SPC    ! Gas chem species names and MWs
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block
      INCLUDE SUBST_GC_EMIS   ! Gas chem emissions name and mapping tables

#ifdef emis_chem
      INCLUDE SUBST_EMPR_CH   ! Emissions processing in chem
#else
      INCLUDE SUBST_EMPR_VD   ! Emissions processing in vdif
#endif

      INCLUDE SUBST_PACTL_ID  ! Process analysis control parameters

#ifdef sens
      INCLUDE SUBST_AE_SPC    ! aerosol species table
      INCLUDE SUBST_NR_SPC    ! non-reactive species table
      INCLUDE SUBST_TR_SPC    ! tracer species table
#endif

C..Arguments:
      REAL, POINTER :: CGRID( :,:,:,: )  ! Species concentrations
      INTEGER JDATE           ! Current date (YYYYDDD)
      INTEGER JTIME           ! Current time (HHMMSS)
      INTEGER TSTEP( 2 )      ! Time step vector (HHMMSS)

C..Parameters:
      REAL, PARAMETER :: PA2ATM = 1.0 / STDATMPA    ! Pascal to atm conv fac

C..External Functions:
      INTEGER INDEX1          ! Looks up name in a list
      INTEGER JUNIT           ! Gets logical device number
      INTEGER SEC2TIME        ! Returns time interval from seconds
      INTEGER TIME2SEC        ! Returns seconds in time interval

C..Saved Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'HRDRIVER'     ! Program name

      INTEGER, SAVE :: ISTFL            ! Unit no. of iteration stat output file
      LOGICAL, SAVE :: LFIRST = .TRUE.  ! Flag for first call to this subroutine

      REAL, SAVE :: MAOMV               ! Mol Wt of air over Mol Wt of water

C..Scratch Local Variables:
      CHARACTER( 132 ) :: MSG           ! Message text
      CHARACTER(  16 ) :: VNAME         ! Name of I/O API data variable

      INTEGER C, E, L, R, S   ! Loop indices

      INTEGER AVGEBI          ! Average no. of EBI iterations
      INTEGER DELT_SEC        ! EBI max time step in seconds
      INTEGER ESP             ! Loop index for emissions species
      INTEGER ITMSTEP         ! Chemistry integration interval (sec)
      INTEGER LEV             ! Layer index
      INTEGER MIDDATE         ! Date at time step midpoint
      INTEGER MIDTIME         ! Time at time step midpoint
      INTEGER MNEBI           ! Min no. of EBI iterations
      INTEGER MXEBI           ! Max no. of EBI iterations
      INTEGER NDARK           ! Number of layer 1 cells in darkness
      INTEGER NPH             ! Index for number of phot. rxns in PHOT
      INTEGER SPC             ! Species loop index
      INTEGER STATUS          ! Status code
      INTEGER VAR             ! Variable number on I/O API file

      LOGICAL LSUNLIGHT       ! Flag for sunlight

      REAL ATMPRES            ! Cell pressure
      REAL CHEMSTEP           ! Chemistry integration interval (min)
      REAL H2O                ! Cell H2O mixing ratio (ppmV)
      REAL SUMEBI             ! Sum of EBI iterations
      REAL TEMP               ! Cell Temperature

      REAL PRES(    NCOLS, NROWS, NLAYS )        ! Cell pressure (Pa)
      REAL QV(      NCOLS, NROWS, NLAYS )        ! Cell water vapor (Kg/Kg air)
      REAL TA(      NCOLS, NROWS, NLAYS )        ! Cell temperature (K)
      REAL RJIN( NPHOTAB )                       ! J-values for a cell
      REAL RJ( NCOLS, NROWS, NLAYS, NPHOTAB )    ! J-values for each cell

      INTEGER     GXOFF, GYOFF            ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

#ifdef sens
!     REAL         SENGRID( NCOLS,NROWS,NLAYS,NPMAX,* ) ! sensitivities
      REAL, POINTER :: SENGRID( :,:,:,:,: )  ! sensitivities
      INTEGER NPMAX         ! number of sens parameters
      INTEGER I, J, IP, RXN
      INTEGER, PARAMETER :: NTXNT = N_GC_SPC * N_GC_SPC
      DOUBLE PRECISION  W( N_GC_SPC )
      REAL*8 SRK( NRXNS )              ! rate constants
      LOGICAL RXNFLAG( NPMAX )
      INTEGER DATENUM, HOURNUM

C Begin code specific to DDM-AERO
      REAL S_SO2I( NPMAX )
      INTEGER NP

      INTEGER N, NSO2, NSULAER, NSULF
C End code specific to DDM-AERO

C Variables used for hddm-3d
      CHARACTER( 16 ), SAVE  :: DDM3D_HIGH = 'DDM3D_HIGH'
      LOGICAL, SAVE          :: HIGH
      LOGICAL, EXTERNAL      :: ENVYN         ! get environment variable as boolean
      CHARACTER( 80 )        :: VARDESC       ! environment variable description

      REAL*8 SRK2 ( NRXNS )            ! rate constants
      LOGICAL ORDER1 ( NRXNS )       ! true if order 1; else, false
      INTEGER HIPARM( 2 )            ! index numbers of the 1st order sens
                                     ! parameters to which we're taking
                                     ! 2nd order sensitivity
      REAL PDT2( N_GC_SPC, N_GC_SPC )  ! Used for 2nd order call of JAC
      REAL SMID( NPMAX, N_GC_SPC )     ! SENGRID in middle of timestep
      REAL SMIDJAC( N_GC_SPC )         ! SMID for the relevant 1st order
                                       ! sensitivity parameter
      REAL*8 RK ( NRXNS )
      INTEGER HITMP1
      INTEGER HITMP2



      INTEGER, POINTER :: IPT( : )
!      INTEGER IPT( NPMAX )            ! Type of sensitivity
                                       ! 1 = INIT; 2 = BOUN; 3 = EMIS; 
                                       ! 4 = HIGH; 5 = RATE
      INTEGER, POINTER :: IPARM( :,: )
!      INTEGER IPARM( NPMAX, MXSPCS )   ! set to 1 if sens to that species
                                        ! set to 0 otherwise

      INTEGER, POINTER :: IDATE( :,: )
!      INTEGER IDATE ( NPMAX, NDAYS ) ! set to 1 if that date is included

      INTEGER      STDATE       ! starting date,    format YYYYDDD

      INTEGER, POINTER :: ITIME( :,: )
!      INTEGER ITIME ( NPMAX, 2 )        ! begin time & end time for each sens
                                        ! default is 0, 240000
                                        ! Note: to span midnight, begin time
                                        ! may be larger than endtime

      REAL, POINTER :: IREGION( :,:,:,: )
!      REAL IREGION( SENCOLS, SENROWS, NLAYS, NPMAX )
      REAL IREGTEMP                  ! Holds relevant value of IREGION

      ! IRXN( NPMAX, NRXNS ) = 1 if sens to that reaction rate;
      !                      = 0 otherwise
      INTEGER, POINTER :: IRXN( :,: )

      REAL RXNMULT                   ! Holds product of concs of reactants
      REAL KSTEP                     ! Holds k times timestep(in min)
      INTEGER NRCT                   ! Counter over reactants
      INTEGER NPROD                  ! Counter over products

      REAL A( N_GC_SPC, N_GC_SPC )   ! An upper triangular matrix and the
                                     ! multipliers used to obtain it
                                     ! (See s_lu.F)
      REAL A1( N_GC_SPC, N_GC_SPC )
      REAL PDT( N_GC_SPC, N_GC_SPC )
      REAL PRD( N_GC_SPC )
      REAL PRD_RATE( N_GC_SPC )      ! Contribution to PRD from rxn rate sens
      REAL SOLD( N_GC_SPC )
      REAL SUM
      INTEGER INFO                   ! see s_lu.F
      INTEGER IPVT ( N_GC_SPC )      ! an integer vector of pivot indices.

#endif


      INTERFACE
#ifdef sens
         SUBROUTINE CKTIME(IDATE,JDATE,STDATE,ITIME,JTIME,TIMEFLAG,P )
           IMPLICIT NONE
           INTEGER,POINTER:: IDATE ( :,: )
           INTEGER,INTENT(IN):: JDATE
           INTEGER,INTENT(IN):: STDATE
           INTEGER,INTENT(IN)::  JTIME
           INTEGER,POINTER:: ITIME (:,: )
           INTEGER,INTENT(IN)::  P
           LOGICAL,INTENT(OUT) :: TIMEFLAG
         END SUBROUTINE CKTIME
#endif
      END INTERFACE

C**********************************************************************

      IF( N_GC_SPC .EQ. 0 ) RETURN

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  On first call, call routines to set-up for EBI solver
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF( LFIRST ) THEN

         LOGDEV = INIT3( )

         IF( MECHNAME .NE. 'CB05CL_AE5_AQ' ) THEN
             MSG ='This version of the EBI solver can only be used with'
     &            // ' the CB05CL_AE5_AQ chemical mechanism'
             CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
         ENDIF

         IF( INDEX( MECHNAME, 'AE' ) .NE. 0 ) THEN
           L_AE_VRSN = .TRUE.
         ELSE
           L_AE_VRSN = .FALSE.
         ENDIF

         IF( INDEX( MECHNAME, 'AQ' ) .NE. 0 ) THEN
           L_AQ_VRSN = .TRUE.
         ELSE
           L_AQ_VRSN = .FALSE.
         ENDIF

         IF( LIRR ) THEN
            MSG = 'IRR Analysis not allowed with EBI solver'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
         ENDIF

         CALL HRINIT

         ITMSTEP = TIME2SEC( TSTEP( 2 ) )
         CHEMSTEP = FLOAT( ITMSTEP ) / 60.0
         WRITE( LOGDEV, 92000 ) CHEMSTEP, DELTAT

         WRITE( LOGDEV, 92020 )
         DO SPC = 1, N_GC_SPC
            WRITE( LOGDEV, 92040 ) GC_SPC( SPC ), RTOL( SPC )
         ENDDO

         MAOMV =  MWAIR / MWWAT

c..If emissions processing requested stop
         IF( EMISCH ) THEN

            MSG = 'ERROR: EBI solver not configured to '//
     &            'process emissions in chemistry'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )

         ENDIF   ! End if doing emissions


#ifdef hrstats
         ISTFL = JUNIT()
         OPEN( UNIT=ISTFL, FILE='iterstat.dat' )
         WRITE( ISTFL, 94020 )
#endif

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         LFIRST = .FALSE.

#ifdef sens
         ALLOCATE( YCDDM (N_GC_SPC) )

         HIGH = .FALSE.         ! default
         VARDESC='Flag for computing higher order sensitivities'
         HIGH = ENVYN( DDM3D_HIGH, VARDESC, HIGH, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            MSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG =
     &          'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, JTIME
         END IF
#endif



      ENDIF      ! First time

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Set date and time to center of time step, get necessary physical
C  data, and get photolysis rates
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MIDDATE = JDATE
      MIDTIME = JTIME
      ITMSTEP = TIME2SEC( TSTEP( 2 ) )
      CHEMSTEP = FLOAT( ITMSTEP ) / 60.0D+00
      CALL NEXTIME( MIDDATE, MIDTIME, SEC2TIME( ITMSTEP / 2 ) )

C.. Compute number of time step loops and step size for EBI solver
      DELT_SEC = DELTAT * 60.0 + 0.1
      IF( DELT_SEC .GE. ITMSTEP ) THEN
         N_EBI_STEPS = 1
         EBI_TMSTEP = FLOAT( ITMSTEP ) / 60.0
      ELSE
         IF( MOD( ITMSTEP, DELT_SEC ) .EQ. 0 ) THEN
            N_EBI_STEPS = ITMSTEP / DELT_SEC
         ELSE
            N_EBI_STEPS = ITMSTEP / DELT_SEC + 1
         ENDIF
         EBI_TMSTEP =  FLOAT( ITMSTEP ) / FLOAT( N_EBI_STEPS ) / 60.0
      ENDIF

      N_INR_STEPS = 1


C.. Get ambient temperature in K

      VNAME = 'TA'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, TA ) ) THEN
         MSG = 'Could not read TA from MET_CRO_3D'
         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF

C.. Get specific humidity in Kg H2O / Kg air
      VNAME = 'QV'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, QV ) ) THEN
         MSG = 'Could not read QV from MET_CRO_3D'
         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF

C.. Get pressure in Pascals
      VNAME = 'PRES'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, PRES ) ) THEN
         MSG = 'Could not read PRES from MET_CRO_3D'
         CALL M3EXIT ( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF

C.. Get photolysis rates in /min
      CALL PHOT ( MIDDATE, MIDTIME, JDATE, JTIME, NDARK, RJ )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Top of loop over cells
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef hrstats
      MNEBI = 1000
      MXEBI = 0
      SUMEBI = 0.0
#endif

#ifdef sens

C Begin code necessary for reaction rate sensitivities

      ! Set the date and hour counters used in sensitivity calls
      DATENUM = 1 + JDATE - STDATE

      DO IP = 1, NPMAX
        IF ( IPT( IP ) .EQ. 5 ) THEN
          ! Rxnflag set to true iff ipt=5 and time, date within bounds
          CALL CKTIME ( IDATE, JDATE, STDATE, ITIME ,
     &                      JTIME, RXNFLAG( IP) , IP)
          IF ( IDATE( IP, DATENUM ) .NE. 1 ) RXNFLAG( IP ) = .FALSE.
        ELSE
          RXNFLAG( IP ) = .FALSE.
        END IF
      END DO

C End code necessary for reaction rate sensitivities

c Begin code specific to DDM-AERO
      VNAME = 'SO2'
      N = INDEX1( VNAME, N_GC_SPC, GC_SPC )
      IF ( N .NE. 0 ) THEN
           NSO2 = N
      ELSE
           MSG = 'Could not find ' // VNAME // 'in gas table'
           CALL M3EXIT ( PNAME, JDATE, JTIME, MSG, XSTAT3 )
      END IF

      VNAME = 'SULRXN'
      N = INDEX1( VNAME, N_GC_SPC, GC_SPC )
      IF ( N .NE. 0 ) THEN
           NSULAER = N
      ELSE
           MSG = 'Could not find ' // VNAME // 'in gas table'
           CALL M3EXIT ( PNAME, JDATE, JTIME, MSG, XSTAT3 )
      END IF

      VNAME = 'SULF'
      N = INDEX1( VNAME, N_GC_SPC, GC_SPC )
      IF ( N .NE. 0 ) THEN
           NSULF = N
      ELSE
           MSG = 'Could not find ' // VNAME // 'in gas table'
           CALL M3EXIT ( PNAME, JDATE, JTIME, MSG, XSTAT3 )
      END IF

c End code specific to DDM-AERO

      IF ( HIGH ) THEN
        DO RXN = 1, NRXNS
          IF( NREACT( RXN ) .EQ. 1 ) THEN
            ORDER1( RXN ) = .TRUE.
          ELSE
            ORDER1( RXN ) = .FALSE.
          END IF
        END DO
      END IF
#endif

      DO L = 1, NLAYS
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS

#ifdef sens
c             DO NP = 1, NPMAX
c                S_SO2I(NP) = SENGRID(C,R,L,NP,NSO2)
c                SENGRID(C,R,L,NP,NSULAER) = 0.0
c             ENDDO
#endif

c..Load ICs
               DO S = 1, N_GC_SPC
                  YC( S ) = MAX( CGRID( C, R, L, S ), 1.0E-30 )
               ENDDO


c..Set physical quantities
               TEMP = TA( C, R, L )
               ATMPRES = PA2ATM * PRES( C, R, L )
               H2O  = MAX ( QV( C, R, L ) * MAOMV *  1.0E+06, 0.0 )


c..Get rate constants
               LSUNLIGHT = .FALSE.
               DO NPH = 1, NPHOTAB
                  RJIN( NPH ) = RJ( C, R, L, NPH )
                  IF( RJ( C, R, L, NPH ) .GT. 0.0 ) LSUNLIGHT = .TRUE.
               ENDDO

               CALL HRCALCKS( NPHOTAB, LSUNLIGHT, RJIN, TEMP,
     &                        ATMPRES, H2O, RKI )


c..Call EBI solver
               N_EBI_IT = 0

               CALL HRSOLVER( JDATE, JTIME, C, R, L )

#ifdef hrstats
               MXEBI  = MAX( MXEBI, N_EBI_IT )
               MNEBI  = MIN( MNEBI, N_EBI_IT )
               SUMEBI = SUMEBI + FLOAT( N_EBI_IT )
#endif


c..Update concentration array
               DO S = 1, N_GC_SPC
                  CGRID( C, R, L, S ) = YC( S )
               ENDDO

#ifdef sens

          DO RXN = 1, NRXNS
            SRK( RXN ) = RKI( RXN )

            IF ( HIGH ) THEN

              IF ( ORDER1 (RXN ) ) THEN
                SRK2( RXN ) = 0.0
              ELSE
                SRK2( RXN ) = RKI( RXN )
              END IF
            END IF

          END DO

          ! Create the Jacobian matrix
          ! YCDDM is set in mid-timestep of hrsolver.F

          CALL JAC ( YCDDM, PDT, SRK )

          DO 433 J = 1, N_GC_SPC
            DO 434 I = 1, N_GC_SPC
              A( I, J ) = 0.0
              A1( I, J ) = 0.0
              A( I, J ) = -0.5 * CHEMSTEP * PDT( I, J )
              IF ( I .EQ. J ) THEN
                  A( I, J ) = 1.0 + A( I, J )
              END IF
              A1( I, J ) = 0.5 * CHEMSTEP * PDT( I, J )
              IF ( I .EQ. J ) THEN
                  A1( I, J ) = 1.0 + A1( I, J )
              END IF
434         CONTINUE
433       CONTINUE


C Factor matrix A by Gaussian elimination
          CALL SGEFA(A, N_GC_SPC, N_GC_SPC, IPVT,INFO)


          DO 495 IP = 1, NPMAX

            IF ( IPT( IP ) .NE. 4 ) THEN

            DO S = 1, N_GC_SPC
              PRD( S ) = 0.0
              PRD_RATE( S ) = 0.0
              SOLD(S) = SENGRID( C, R, L, IP,S )
              IF (ABS(SOLD( S )) .LT. 1.e-25 ) THEN
                SOLD(S) = 0.
              END IF
            END DO


C Begin code specific to reaction rate sensitivities
              ! RXNFLAG is true if IPT = 5 and time, date within bounds
              IF ( RXNFLAG( IP ) ) THEN

                IREGTEMP = IREGION( C, R, L, IP )

                DO RXN = 1, NRXNS

                  ! This checks if it's a reaction in the sens parameter
                  IF ( IRXN( IP, RXN ) ) THEN

                    ! KSTEP = k * timestep(in min)
                    KSTEP = SRK( RXN ) * CHEMSTEP

                    ! RXNMULT is the product of the concs of the reactants
                    ! Note that the first 3 slots of IRR are for reactants,
                    ! and slots 4- are for products
                    IF ( NREACT( RXN ) .EQ. 1 ) THEN
                      RXNMULT = KSTEP
     &                        * YCDDM( IRR( RXN, 1 ) )
                    ELSE IF ( NREACT( RXN ) .EQ. 2 ) THEN
                      RXNMULT = KSTEP
     &                        * YCDDM( IRR( RXN, 1 ) )
     &                        * YCDDM( IRR( RXN, 2 ) )
                    ELSE IF ( NREACT( RXN ) .EQ. 3 ) THEN
                      RXNMULT = KSTEP
     &                        * YCDDM( IRR( RXN, 1 ) )
     &                        * YCDDM( IRR( RXN, 2 ) )
     &                        * YCDDM( IRR( RXN, 3 ) )
                    ELSE
                      MSG = 'NREACT out of expected bounds of 1-3.'
                      CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
                    END IF

                    ! Loop over the reactants
                    DO NRCT = 1, NREACT( RXN )
                      ! Subtract RXNMULT from PRD_RATE for reactant species
                      PRD_RATE( IRR( RXN, NRCT ) ) =
     &                       PRD_RATE( IRR( RXN, NRCT ) )
     &                     - RXNMULT * IREGTEMP
                    END DO

                    ! Loop over the products
                    DO NPROD = 1, NPRDCT( RXN )
                      ! Add RXNMULT to PRD_RATE for product species
                      ! The +3 reflects that slots >=4 of IRR are for products
                      ! SC are the stoichiometric coefficients of products
                      !    and do not need the +3 (see RXDT.EXT)
                      PRD_RATE( IRR( RXN, 3+NPROD ) ) =
     &                    PRD_RATE( IRR( RXN, 3+NPROD ) )
     &                + ( RXNMULT * SC( RXN, NPROD ) * IREGTEMP )
                    END DO

                  END IF

                END DO  ! RXN

              END IF  ! RXNFLAG

C End code specific to reaction rate sensitivities

            DO S = 1, N_GC_SPC
              SUM = 0.0
              DO J = 1, N_GC_SPC
                SUM = SUM + A1( S, J ) * SOLD( J )
              END DO
C edits by A.Digar
              PRD( S ) = SUM + PRD_RATE( S )
C end edit
            END DO

            CALL SGESL( A, N_GC_SPC, N_GC_SPC, IPVT, PRD, 0 )

            DO S = 1, N_GC_SPC
              IF ( ABS ( PRD ( S ) ) .LT. 1.e-25 ) THEN

                IF ( HIGH ) THEN
                  SMID( IP, S ) =
     &                0.5 * SENGRID( C, R, L, IP, S )
                END IF

                SENGRID( C, R, L, IP, S ) = 0.
              ELSE

                IF ( HIGH ) THEN

                  ! SMID is the average of SENGRID before
                  ! and after chemistry
                  SMID( IP, S ) = 0.5 *
     &                      ( SENGRID( C, R, L, IP, S ) + PRD( S) )
                END IF

                SENGRID( C, R, L, IP, S ) = PRD( S )
              END IF
            END DO

            ELSE

              ! 2nd order sensitivity
              HIPARM( 1 ) = 0
              HIPARM( 2 ) = 0
              DO J = 1, IP - 1
                IF ( IPARM( IP, J ) .EQ. 1 ) THEN
                  HIPARM( 1 ) = J
                ELSE IF ( IPARM( IP, J ) .EQ. 2 ) THEN
                  HIPARM( 2 ) = J
                ELSE IF ( IPARM( IP, J ) .EQ. 3 ) THEN
                  HIPARM( 1 ) = J
                  HIPARM( 2 ) = J
                END IF
              END DO

              DO S = 1, N_GC_SPC
                SMIDJAC( S ) = SMID( HIPARM( 1 ), S )
              END DO
C Added by A.Digar
              DO S = 1, N_GC_SPC
                PRD( S ) = 0.0
                SOLD(S) = SENGRID( C, R, L,IP,S )
                IF (ABS(SOLD( S )) .LT. 1.e-25 ) THEN
                  SOLD(S) = 0.
                END IF
              END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Begin code specific to high-order sensitivity with one/more 1st order term/s
C being reaction rate sensitivity
C added by A.Digar
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

              HITMP1 = HIPARM( 1 )
              HITMP2 = HIPARM( 2 )

              DO N = 1, 2 ! loop for two 1st-order sens parameters

                 IF ( ( IPT( HITMP1 ) .EQ. 5 ) .AND. ( RXNFLAG( HITMP1 ) ) ) THEN
                      ! check for rate constant sens, date & time

                  IREGTEMP = IREGION ( C, R, L, HITMP1 )

                  DO RXN = 1, NRXNS
                
                   ! keeping the rate terms that contain only the
                   ! rate constants of interest and setting rest to zero

                   RK( RXN ) = SRK( RXN ) * IRXN( HITMP1, RXN )

                  END DO

                   ! Jacobian for first-order,
                   ! called with sensitivities and
                   ! rxn rates with 1st order rxns effective  

                   CALL JAC ( YCDDM, PDT, RK )

                   DO S = 1, N_GC_SPC
                     DO J = 1, N_GC_SPC
                       PRD( S ) = PRD( S ) + CHEMSTEP * PDT( S,J ) * SMID( HITMP2,J ) * IREGTEMP
                     END DO
                   END DO

                   IF ( IPT( HITMP1 ) .eq. IPT( HITMP2 ) ) THEN
                       PRD = 2.0 * PRD
                       EXIT
                   ENDIF

                 ENDIF

                 HITMP1 = HIPARM( 2 )
                 HITMP2 = HIPARM( 1 )

              END DO

C End of modification by A.Digar

              ! Jacobian for higher-order,
              ! called with sensitivities and
              ! rxn rates with 1st order rxns removed
              CALL JAC ( SMIDJAC, PDT2, SRK2 )

C removed by A.Digar
c              DO S = 1, N_GC_SPC
c                PRD( S ) = 0.0
c                SOLD(S) = SENGRID( C, R, L,IP,S )
c                IF (ABS(SOLD( S )) .LT. 1.e-25 ) THEN
c                  SOLD(S) = 0.
c                END IF
c              END DO
C end of edits

              DO S = 1, N_GC_SPC
                SUM = 0.0
                DO J = 1, N_GC_SPC
                  SUM = SUM + A1( S, J ) * SOLD( J )
     &              + CHEMSTEP * PDT2( S,J ) * SMID( HIPARM( 2 ),J )
                END DO
C edits by A.Digar
                PRD( S ) = SUM + PRD( S )
C end of edits
              END DO

              CALL SGESL( A, N_GC_SPC, N_GC_SPC, IPVT, PRD, 0 )

              DO S = 1, N_GC_SPC
                IF ( ABS ( PRD ( S ) ) .LT. 1.e-25 ) THEN
                  SENGRID( C, R, L, IP, S ) = 0.
                ELSE
                  SENGRID( C, R, L, IP, S ) = PRD( S )
                END IF
              END DO

            END IF

495         CONTINUE    ! End of do loop over IP


C Begin code specific to DDM-AERO
c           DO NP = 1, NPMAX
c            SENGRID(C,R,L,NP,NSULAER)=S_SO2I(NP)-SENGRID(C,R,L,NP,NSO2)
c            SENGRID(C,R,L,NP,NSULF)=SENGRID(C,R,L,NP,NSULF)
c    &                              +SENGRID(C,R,L,NP,NSULAER)
c           ENDDO
c End code specific to DDM-AERO
#endif

            ENDDO
         ENDDO
      ENDDO

#ifdef hrstats
      AVGEBI = SUMEBI / FLOAT( NCOLS * NROWS * NLAYS )
      WRITE( ISTFL, 94040 ) JDATE, JTIME, MNEBI, AVGEBI, MXEBI
#endif

      RETURN

C*********************** FORMAT STATEMENTS ****************************

92000 FORMAT( / 10X, 'Euler Backward Iterative Parameters -'
     &        / 10X, 'Chemistry Integration Time Interval (min):', F12.4,
     &        / 10X, 'EBI maximum time step (min):              ', F12.4 )

92020 FORMAT( //10X, 'Species convergence tolerances:' )

92040 FORMAT(   10X, A16, 2X, 1PE12.2 )

92060 FORMAT( / 10X, 'Emissions Processing in Chemistry ...'
     &        / 10X, 'Number of Emissions Layers:         ', I3
     &        / 10X, 'out of total Number of Model Layers:', I3 )


94020 FORMAT( 'DATE      TIME ', 'MNEBI AVEBI MXEBI' )

94040 FORMAT( I7, 1X, I6, 1X, 3( I5, 1X ) )
      END
