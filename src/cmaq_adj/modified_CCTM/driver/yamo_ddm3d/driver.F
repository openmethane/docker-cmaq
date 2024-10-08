
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_ddm3d/driver.F,v 1.1.1.1 2010/06/14 16:03:04 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)driver.F  1.2 /project/mod3/CMAQ/src/driver/ctm/SCCS/s.driver.F 07 Jul 1997 09:20:18

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      PROGRAM  DRIVER

C-----------------------------------------------------------------------
C Function:
C    CMAQ CTM driver
 
C Preconditions:
C    Initialized file CONCFILE for output; completed
C    files HISTORY containing initial conditions, SPCCONST for
C    conversion of concentration field from computational units
C    to output units.
 
C Subroutines and functions called:
C    INITSCEN, ADVSTEP, M3EXIT, WRITE3
C    science processes SCIPROC, PA_OUTPUT
 
C Revision History:
C    prototype 6/92 by CJC for proof-of-concept
C    Revised   2/93 by CJC for initial LCM Prototype Alpha
 
C    31 August, 1995 by M. Talat Odman at NCSC: special version for one 
C    single grid
 
C    16 April 1995 by M. Talat Odman at NCSC: write (or rewrite if restart)
C    initial conditions to the output file

C    Jeff
C    18 June 98 - put load of mechanism common here because of ping/ping_noop,
C    chem/chem_noop options

C    2 October, 1998 by Al Bourgeois at LM: parallel implementation
C    Jeff - Dec 00 - move CGRID_MAP into f90 module, re-order cols/rows in
C    call to PAR_INIT
C    Jeff - Jul 01 - enable integral average conc data
C    Sep 01  J. Young        Dyn Alloc - Use HGRD_DEFN

C    3 Sep 01 David wong
C      -- removed M3IO SHUT3 call which is done in PAR_TERM
C      -- removed SET_CTMFILE call

C   23 Jun 03 J.Young: for layer dependent advection tstep
C   18 Aug 03 J. Pleim - move vdiff before advection
C   07 Dec 04 J.Young: for layer dyn alloc - Use VGRD_DEFN
C   30 May 05 J.Young: mass-conserving advection (yamo)
C   20 Jan 06 J.Young: add circular buffer CGRID state file
C   24 May 06 J.Young: par_init/pio_init col/row order check
C    6 Sep 06 J.Young: one-write cgrid file; SGRID in module
C   15 Nov 08 S.Napelenok: DDM-3D implementation for CMAQv4.7 
C-----------------------------------------------------------------------

      USE PCGRID_DEFN           ! inherits HGRD_DEFN and CGRID_SPCS
      USE VGRD_DEFN             ! vertical layer specifications
      USE STD_CONC              ! standard CONC
      USE AVG_CONC              ! integral average CONC
      USE WVEL_DEFN             ! derived vertical velocity component
#ifdef sens
      USE S_PCGRID_DEFN
#endif

      IMPLICIT NONE

C Include Files:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! I/O definitions and declarations
      INCLUDE SUBST_PACTL_ID    ! PA control parameters

C Load the mechanism COMMON ... (next two INCLUDE files)
      INCLUDE SUBST_RXCMMN      ! chemical mechamism reactions COMMON
      INCLUDE SUBST_RXDATA      ! chemical mechamism reactions DATA

#ifdef sens
      INCLUDE SUBST_GC_SPC              ! gas chemistry species table
      INCLUDE SUBST_AE_SPC              ! aerosol species table
      INCLUDE SUBST_NR_SPC              ! non-reactive species table
      INCLUDE SUBST_TR_SPC              ! tracer species table
#endif

C External Functions (not already declared by IODECL3.EXT):

      LOGICAL, EXTERNAL :: ENVYN ! get environment variable as boolean

C Local variables:

      INTEGER      LOGDEV       ! FORTRAN unit number for log file
      INTEGER      STDATE       ! starting date,    format YYYYDDD
      INTEGER      STTIME       ! starting time,    format HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
      INTEGER, ALLOCATABLE, SAVE :: ASTEP( : )
      INTEGER      NSTEPS       ! run duration: number of output time steps
      INTEGER      NREPS        ! number of model time steps per output step
      INTEGER      ISTEP        ! current output time step number
      INTEGER      IREP         ! model step number within this output step
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      C, R, L, K, S, V     ! loop induction variables
      INTEGER      ALLOCSTAT

      CHARACTER(  2 ) :: COLROW = 'CR'  ! col/row arg list order
      CHARACTER( 16 ) :: PNAME = 'DRIVER'
      CHARACTER( 96 ) :: XMSG = ' '

      REAL, POINTER     :: CGRID( :,:,:,: )
      REAL, ALLOCATABLE :: AGRID( :,:,:,: )
      REAL    DIVFAC      ! trapezoidal average factor
      INTEGER A_NLYS

      INTEGER PAR_ERR     ! Error code from parallel initialization. 
      REAL    CLOCK       ! Wall-clock time (sec) at initialization.

#ifdef sens
C BEGIN SENSITIVITY VARIABLES
      REAL, POINTER :: SENGRID( :,:,:,:,: )
      REAL, ALLOCATABLE :: SAGRID( :,:,:,:,: )
      INTEGER  P
      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_IPT( : )
      INTEGER, POINTER :: IPT (:)
!     INTEGER IPT( NPMAX )              ! Type of sensitivity
                                        ! 1 = INIT; 2 = BOUN;
                                        ! 3 = EMIS; 4 = HIGH;
                                        ! 5 = RATE
      INTEGER       MXSPCS              ! Number of species in CGRID
      PARAMETER   ( MXSPCS = N_GC_SPCD
     &                     + N_AE_SPC
     &                     + N_NR_SPC
     &                     + N_TR_SPC )

      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_IPARM( :,: )
      INTEGER, POINTER :: IPARM ( :,: )
!     INTEGER IPARM( NPMAX, MXSPCS )    ! set to 1 if sens to that species
                                        ! set to 0 otherwise
      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_ISPLIT( :,: )
      INTEGER, POINTER :: ISPLIT( :,: )
!     INTEGER ISPLIT( NPMAX, 6 )        ! 1 = total; 2 = area; 3 = bio
                                        ! 4 = mobile; 5 = point; 6 = non-rd
      REAL, ALLOCATABLE, SAVE, TARGET :: TGT_IAMOUNT( :,:,: )
      REAL, POINTER :: IAMOUNT( :,:,: )
!     REAL IAMOUNT ( NPMAX, MXSPCS, 25 )  ! emission rate during the hour
      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_ILAYER( :,: )
      INTEGER, POINTER :: ILAYER( :,: )
!     INTEGER ILAYER ( NPMAX, NLAYS )   ! 1 if that layer's emis included,
                                        ! 0 else
                                        ! only used if emis amount specified,
                                        ! otherwise set to 1 for all layers

      ! IRXN( NPMAX, NRXNS ) = 1 if sens to that reaction rate; 
      !                      = 0 otherwise 
      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_IRXN( :,: )
      INTEGER, POINTER :: IRXN( :,: )

      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_IDATE( :,: )
      INTEGER, POINTER :: IDATE( :,: )
!     INTEGER IDATE ( NPMAX, NDAYS )    ! date(s) on which emission or boundary
                                        ! condition should occur
      INTEGER, ALLOCATABLE, SAVE, TARGET :: TGT_ITIME( :,: )
      INTEGER, POINTER :: ITIME( :,: )
!     INTEGER ITIME ( NPMAX, 2 )        ! begin time & end time for each sens
                                        ! default is 0, 240000
                                        ! Note: to span midnight, begin time
                                        ! may be larger than endtime
      REAL, ALLOCATABLE, SAVE, TARGET :: TGT_IREGION( :,:,:,:)
!     REAL IREGION( SENCOLS, SENROWS, NLAYS, NPMAX )
      REAL, POINTER :: IREGION( :,:,:,: )
      CHARACTER( 8 ), ALLOCATABLE, SAVE, TARGET :: TGT_SENPAR( : )
      CHARACTER( 8 ), POINTER :: SEN_PAR( : )
!     CHARACTER( 8 ) :: SEN_PAR(NPMAX)  ! names for sens parameter
C max number of sensitivity parameters
      CHARACTER( 16 ) :: CTM_NPMAX = 'CTM_NPMAX'
C  environment variable start time
      CHARACTER( 16 ) :: CTM_STTIME = 'CTM_STTIME'
C  environment variable run duration
      CHARACTER( 16 ) :: CTM_RUNLEN = 'CTM_RUNLEN'
      INTEGER RUNLEN, BGNTIME
      INTEGER, EXTERNAL :: ENVINT       !  get environment variable as integer
      INTEGER      STATUS               !  ENVINT status
      CHARACTER( 80 ) :: VARDESC
      INTEGER      NPMAX, NDAYS
C END SENSITIVITY VARIABLES
#endif

      INTERFACE
#ifdef sens
         SUBROUTINE INITSCEN
     &                      ( CGRID, SENGRID, STDATE, STTIME, TSTEP,
     &                        NSTEPS, NPMAX, IPT, IPARM, ISPLIT,
     &                        IAMOUNT, ILAYER, NDAYS, IDATE, ITIME,
     &                        IREGION, IRXN, SEN_PAR )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( OUT )    :: STDATE, STTIME
            INTEGER, INTENT( OUT )    :: TSTEP( 2 )
            INTEGER, INTENT( OUT )    :: NSTEPS
            INTEGER, INTENT( IN )     :: NPMAX
            INTEGER, POINTER    :: IPT( : )
            INTEGER, POINTER    :: IPARM( :,: )
            INTEGER, POINTER    :: ISPLIT( :,: )
            REAL, POINTER       :: IAMOUNT( :,:,: )
            INTEGER, POINTER    :: ILAYER( :,: )
            INTEGER, INTENT( IN ) :: NDAYS
            INTEGER, POINTER    :: IDATE( :,: )
            INTEGER, POINTER    :: ITIME( :,: )
            REAL, POINTER        :: IREGION( :,:,:,: )
            INTEGER, POINTER     :: IRXN( :,: )
            CHARACTER ( * ), POINTER :: SEN_PAR( : )
         END SUBROUTINE INITSCEN
#else
         SUBROUTINE INITSCEN ( CGRID, STDATE, STTIME, TSTEP, NSTEPS )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( OUT )    :: STDATE, STTIME, TSTEP( 2 )
            INTEGER, INTENT( OUT )    :: NSTEPS
         END SUBROUTINE INITSCEN
#endif
         SUBROUTINE ADVSTEP ( JDATE, JTIME, TSTEP, ASTEP, NREPS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( INOUT )  :: TSTEP( 2 )
            INTEGER, INTENT( OUT )    :: ASTEP( : )
            INTEGER, INTENT( OUT )    :: NREPS
         END SUBROUTINE ADVSTEP
         SUBROUTINE CKSUMMER ( PRNAME, CGRID, JDATE, JTIME )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PRNAME
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
         END SUBROUTINE CKSUMMER
         SUBROUTINE PA_INIT ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME, TSTEP( 2 )
         END SUBROUTINE PA_INIT
#ifdef sens
         SUBROUTINE SCIPROC ( CGRID, SENGRID, JDATE, JTIME,
     &                        TSTEP, ASTEP, NPMAX, IPT, IPARM, ISPLIT,
     &                        IAMOUNT, ILAYER, STDATE, IDATE,
     &                        ITIME, IREGION, IRXN, SEN_PAR )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN OUT ) :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: ASTEP( : )
            INTEGER, INTENT( IN )     :: NPMAX
            INTEGER, POINTER     :: IPT( : )
            INTEGER, POINTER     :: IPARM( :,: )
            INTEGER, POINTER     :: ISPLIT( :,: )
            REAL, POINTER        :: IAMOUNT( :,:,: )
            INTEGER, POINTER     :: ILAYER( :,: )
            INTEGER, INTENT( IN )     :: STDATE
            INTEGER, POINTER     :: IDATE( :,: )
            INTEGER, POINTER     :: ITIME( :,: )
            REAL, POINTER        :: IREGION( :,:,:,: )
            INTEGER, POINTER     :: IRXN( :,: )
            CHARACTER ( * ), POINTER :: SEN_PAR( : )
         END SUBROUTINE SCIPROC
#else
         SUBROUTINE SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN OUT ) :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 ), ASTEP( : )
         END SUBROUTINE SCIPROC
#endif
         SUBROUTINE WR_ACONC ( AGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL,    INTENT( IN )     :: AGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME, TSTEP
         END SUBROUTINE WR_ACONC
         SUBROUTINE WR_CGRID ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME, TSTEP
         END SUBROUTINE WR_CGRID
#ifdef sens
         SUBROUTINE WR_ASENS ( SAGRID, JDATE, JTIME, TSTEP,
     &                         NPMAX, SEN_PAR )
            IMPLICIT NONE
            REAL,    INTENT( IN )     :: SAGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP
            INTEGER, INTENT( IN )     :: NPMAX
            CHARACTER ( * ), POINTER :: SEN_PAR( : )
         END SUBROUTINE WR_ASENS
         SUBROUTINE WR_SENGRID( SENGRID, JDATE, JTIME, TSTEP, NPMAX,
     &                          SEN_PAR )
            IMPLICIT NONE
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
            CHARACTER ( * ), POINTER  :: SEN_PAR( : )
         END SUBROUTINE WR_SENGRID
#endif
         SUBROUTINE PA_OUTPUT ( CGRID, JDATE, JTIME )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
         END SUBROUTINE PA_OUTPUT
      END INTERFACE

C-----------------------------------------------------------------------

      JDATE = 0; JTIME = 0

C Get number of species, and starting indices for CGRID array.

      CALL CGRID_MAP ( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C Start up processor communication and define horizontal domain decomposition 
C and vertical layer structure

      CALL PAR_INIT( COLROW, NSPCSD, CLOCK, PAR_ERR )

      IF ( PAR_ERR .NE. 0 ) THEN
         XMSG = 'Error in PAR_INIT'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, 1 )
         END IF

      LOGDEV = INIT3 ()

      IF ( NSPCSD .GT. MXVARS3 ) THEN
         WRITE( XMSG,'(5X, A, I5, A)' ) 'The number of variables,', NSPCSD,
     &   ' to be written to the State CGRID File'
         WRITE( LOGDEV, '( A )' ) XMSG
         WRITE( XMSG,'(5X, A, I5)' ) 'exceeds the I/O-API limit:', MXVARS3
         WRITE( LOGDEV, '( A )' ) XMSG
         XMSG = 'Recompile with an I/O-API lib having a larger MXVARS3'
         WRITE( LOGDEV, '( 5X, A )' ) XMSG
         CALL M3EXIT( PNAME, JDATE, JTIME, ' ', XSTAT1 )
      END IF

C Initialize PCGRID

      IF ( .NOT. PCGRID_INIT () ) THEN
         XMSG = 'Failure defining horizontal domain'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2  )
         END IF
      
      CGRID => PCGRID( 1:MY_NCOLS,1:MY_NROWS,:,: )   ! required for PinG

#ifdef sens
      NPMAX = 1    ! default
      VARDESC = 'Max number of sens parameters.'
      NPMAX = ENVINT( CTM_NPMAX, VARDESC, NPMAX, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Envt variable set, but empty...Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, 0
            ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, 0
            END IF

C Use the begin time and run length to determine how many
C dates are in the model run, for use in dimensioning IDATE
      BGNTIME = 000000         ! default
      VARDESC = 'Scenario Starting Time (HHMMSS)'
      BGNTIME = ENVINT( CTM_STTIME, VARDESC, BGNTIME, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, BGNTIME, XMSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Envirnmt variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
            END IF

      RUNLEN = 480000         ! default
      VARDESC = 'Scenario Run Duration (HHMMSS)'
      RUNLEN = ENVINT( CTM_RUNLEN, VARDESC, RUNLEN, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, STDATE, BGNTIME, XMSG, XSTAT2 )
            ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Envirnmt variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, BGNTIME
            END IF

      NDAYS = 1 + FLOOR( REAL( BGNTIME + RUNLEN ) / 240000.0 )

C Initialize S_PCGRID

      IF ( .NOT. S_PCGRID_IN( NPMAX ) ) THEN
         XMSG = 'Failure defining horizontal sens domain'
         CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT2  )
         END IF

      SENGRID => S_PCGRID(1:MY_NCOLS,1:MY_NROWS,:,:,:)
#endif

C Initalize CONC definitions (in STD_CONC F90 module)

      CALL CONC_DEFN ()

C Get avg CONC definitions, species and layer pointers (in AVG_CONC F90 module)

      CALL A_CONC_DEFN ()
      A_NLYS = ACONC_ELEV - ACONC_BLEV + 1

C Initialize optional derived vertical velocity writes to conc file

      IF ( .NOT. WVEL_INIT () ) THEN
         XMSG = 'Failure initializing derived vertical velocity writes'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2  )
         END IF

#ifdef sens
      ! Allocate IREGION
      ALLOCATE   ( TGT_IREGION( NCOLS,NROWS,NLAYS,NPMAX ),
     &             STAT = ALLOCSTAT ) 
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_IREGION'
               CALL M3EXIT ( 'TGT_IREGION', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_IREGION = 0.0
      IREGION => TGT_IREGION


      ! Allocate IPT
      ALLOCATE   ( TGT_IPT( NPMAX ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_IPT'
               CALL M3EXIT ( 'TGT_IPT', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_IPT = 0
      IPT => TGT_IPT


      ! Allocate IPARM
      ALLOCATE   ( TGT_IPARM( NPMAX, MXSPCS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_IPARM'
               CALL M3EXIT ( 'TGT_IPARM', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_IPARM = 0
      IPARM => TGT_IPARM


      ! Allocate ISPLIT
      ALLOCATE   ( TGT_ISPLIT( NPMAX, 6 ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_ISPLIT'
               CALL M3EXIT ( 'TGT_ISPLIT', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_ISPLIT = 0
      ISPLIT => TGT_ISPLIT


      ! Allocate IAMOUNT
      ALLOCATE   ( TGT_IAMOUNT( NPMAX, MXSPCS, 25 ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_IAMOUNT'
               CALL M3EXIT ( 'TGT_IAMOUNT', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_IAMOUNT = 0.0
      IAMOUNT => TGT_IAMOUNT

      ! Allocate ILAYER
      ALLOCATE   ( TGT_ILAYER( NPMAX, NLAYS), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_ILAYER'
               CALL M3EXIT ( 'TGT_ILAYER', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_ILAYER = 0
      ILAYER => TGT_ILAYER


      ! Allocate IRXN
      ALLOCATE   ( TGT_IRXN( NPMAX,NRXNS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_IRXN'
               CALL M3EXIT ( 'TGT_IRXN', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_IRXN = 0
      IRXN => TGT_IRXN


      ! Allocate IDATE
      ALLOCATE   ( TGT_IDATE( NPMAX, NDAYS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_IDATE'
               CALL M3EXIT ( 'TGT_IDATE', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_IDATE = 0
      IDATE => TGT_IDATE


      ! Allocate ITIME
      ALLOCATE   ( TGT_ITIME( NPMAX, 2 ),
     &           STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_ITIME'
               CALL M3EXIT ( 'TGT_ITIME', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_ITIME = 0
      ITIME => TGT_ITIME

      ! Allocate SENPAR
      ALLOCATE   ( TGT_SENPAR( NPMAX ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TGT_SENPAR'
               CALL M3EXIT ( 'TGT_SENPAR', 0, 0, XMSG, XSTAT2 )
               END IF
            TGT_SENPAR = ' '
      SEN_PAR => TGT_SENPAR
#endif

C Initialize conc field: Copy IC's to CONC file as step 0
C Convention: the input file concentration units are always ppmV.

#ifdef sens
      CALL INITSCEN ( CGRID, SENGRID, STDATE, STTIME,
     &                TSTEP, NSTEPS, NPMAX, IPT, IPARM,
     &                ISPLIT, IAMOUNT, ILAYER, NDAYS, IDATE,
     &                ITIME, IREGION, IRXN, SEN_PAR)
#else
      CALL INITSCEN ( CGRID, STDATE, STTIME, TSTEP, NSTEPS )
#endif
      CALL CKSUMMER ( 'INITSCEN', CGRID, STDATE, STTIME )

C Verify input file header consistency with COORD.EXT and run duration

      CALL FLCHECK ( STDATE, STTIME, TSTEP( 1 ) )

      JDATE = STDATE; JTIME = STTIME

      IF ( LIPR .OR. LIRR ) CALL PA_INIT ( CGRID, STDATE, STTIME, TSTEP )

      ALLOCATE ( AGRID( MY_NCOLS,MY_NROWS,A_NLYS,N_ASPCS ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'AGRID memory allocation failed'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      AGRID = 0.0

#ifdef sens
      ALLOCATE ( SAGRID( MY_NCOLS,MY_NROWS,A_NLYS,NPMAX,N_ASPCS ),
     &           STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'SAGRID memory allocation failed'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      SAGRID = 0.0
#endif

      ALLOCATE ( ASTEP( NLAYS ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'ASTEP memory allocation failed'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Main processing loop:

      DO 201 ISTEP = 1, NSTEPS   ! output time step loop

C Get synchronization and advection time steps, TSTEP(2), ASTEP(L) and NREPS

         CALL ADVSTEP ( JDATE, JTIME, TSTEP, ASTEP, NREPS )

         DO V = 1, N_ASPCS
            S = AVG_CONC_MAP( V )
            L = 0
            DO K = ACONC_BLEV, ACONC_ELEV
               L = L + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     AGRID( C,R,L,V ) = CGRID( C,R,K,S )
#ifdef sens
                     DO P = 1, NPMAX
                       SAGRID( C,R,L,P,V ) = SENGRID( C,R,K,P,S )
                       END DO
#endif
                     END DO
                  END DO
               END DO
            END DO

C science process sequence:
         
         DIVFAC = 0.5 / FLOAT( NREPS )

         DO IREP = 1, NREPS - 1

#ifdef sens
            CALL SCIPROC (CGRID, SENGRID, JDATE, JTIME, TSTEP, ASTEP,
     &                      NPMAX, IPT, IPARM, ISPLIT,
     &                      IAMOUNT, ILAYER, STDATE, IDATE,
     &                      ITIME, IREGION, IRXN, SEN_PAR )
#else
            CALL SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
#endif

            DO V = 1, N_ASPCS
               S = AVG_CONC_MAP( V )
               L = 0
               DO K = ACONC_BLEV, ACONC_ELEV
                  L = L + 1
                  DO R = 1, MY_NROWS
                     DO C = 1, MY_NCOLS
                        AGRID( C,R,L,V ) = AGRID( C,R,L,V )
     &                                   + 2.0 * CGRID( C,R,K,S )
#ifdef sens
                        DO P = 1, NPMAX
                           SAGRID( C,R,L,P,V ) = SAGRID( C,R,L,P,V )
     &                                   + 2.0 * SENGRID( C,R,K,P,S )
                           END DO  ! P
#endif
                        END DO
                     END DO
                  END DO
               END DO

            END DO

#ifdef sens
         CALL SCIPROC ( CGRID, SENGRID, JDATE, JTIME, TSTEP, ASTEP,
     &                    NPMAX, IPT, IPARM, ISPLIT,
     &                    IAMOUNT, ILAYER, STDATE, IDATE,
     &                    ITIME, IREGION, IRXN, SEN_PAR )
#else
         CALL SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
#endif

         DO V = 1, N_ASPCS
            S = AVG_CONC_MAP( V )
            L = 0
            DO K = ACONC_BLEV, ACONC_ELEV
               L = L + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     AGRID( C,R,L,V ) = DIVFAC * ( AGRID( C,R,L,V )
     &                                +            CGRID( C,R,K,S ) )
#ifdef sens
                     DO P = 1, NPMAX
                        SAGRID( C,R,L,P,V ) = DIVFAC *
     &                    ( SAGRID( C,R,L,P,V ) + SENGRID( C,R,K,P,S ) )
                        END DO
#endif
                     END DO
                  END DO
               END DO
            END DO

         DO V = 1, N_CSPCS
            S = CONC_MAP( V )
            L = 0
            DO K = CONC_BLEV, CONC_ELEV
               L = L + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     SGRID( C,R,L,V ) = CGRID( C,R,K,S )
                     END DO
                  END DO
               END DO
            END DO

C write conc fields

         CALL WR_CONC ( JDATE, JTIME, TSTEP( 1 ) )

         CALL WR_ACONC ( AGRID, JDATE, JTIME, TSTEP( 1 ) )

         IF ( LIPR .OR. LIRR ) CALL PA_OUTPUT ( CGRID, JDATE, JTIME )

#ifdef sens
C Write hourly-averaged sensitivity file
         CALL WR_ASENS ( SAGRID, JDATE, JTIME, TSTEP( 1 ) ,
     &                   NPMAX, SEN_PAR)
#endif

201      CONTINUE

      CALL WR_CGRID ( CGRID, JDATE, JTIME, TSTEP( 1 ) )

#ifdef sens
      CALL WR_SENGRID ( SENGRID, JDATE, JTIME, TSTEP, NPMAX, SEN_PAR )
#endif

C Clean up communications.

      CALL PAR_TERM ( PNAME, JDATE, JTIME, CLOCK )

      END
