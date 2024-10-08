
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/ctm/driver.F,v 1.1.1.1 2010/06/14 16:03:04 sjr Exp $ 

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
C   20 Jan 06 J.Young: add circular buffer CGRID state file
C   24 May 06 J.Young: par_init/pio_init col/row order check
C    6 Sep 06 J.Young: one-write cgrid file; SGRID in module
C-----------------------------------------------------------------------

      USE PCGRID_DEFN           ! inherits HGRD_DEFN and CGRID_SPCS
      USE VGRD_DEFN             ! vertical layer specifications
      USE STD_CONC              ! standard CONC
      USE AVG_CONC              ! integral average CONC

      IMPLICIT NONE

C Include Files:

!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! I/O definitions and declarations
      INCLUDE SUBST_PACTL_ID    ! PA control parameters

C Load the mechanism COMMON ... (next two INCLUDE files)
      INCLUDE SUBST_RXCMMN      ! chemical mechamism reactions COMMON
      INCLUDE SUBST_RXDATA      ! chemical mechamism reactions DATA

C External Functions (not already declared by IODECL3.EXT):

      LOGICAL, EXTERNAL :: ENVYN ! get environment variable as boolean

C Local variables:

      INTEGER      LOGDEV       ! FORTRAN unit number for log file
      INTEGER      STDATE       ! starting date,    format YYYYDDD
      INTEGER      STTIME       ! starting time,    format HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
!     INTEGER      ASTEP( NLAYS )
!     INTEGER, ALLOCATABLE :: ASTEP( : )
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
!     REAL, ALLOCATABLE :: SGRID( :,:,:,: )
      REAL, ALLOCATABLE :: AGRID( :,:,:,: )
!     INTEGER C_NLYS
      REAL    DIVFAC      ! trapezoidal average factor
      INTEGER A_NLYS

      INTEGER PAR_ERR     ! Error code from parallel initialization. 
      REAL    CLOCK       ! Wall-clock time (sec) at initialization.

      INTERFACE
         SUBROUTINE INITSCEN ( CGRID, STDATE, STTIME, TSTEP, NSTEPS )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( OUT )    :: STDATE, STTIME, TSTEP( 2 )
            INTEGER, INTENT( OUT )    :: NSTEPS
         END SUBROUTINE INITSCEN
         SUBROUTINE ADVSTEP ( JDATE, JTIME, TSTEP, ASTEP, NREPS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
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
         SUBROUTINE SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN OUT ) :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 ), ASTEP( : )
         END SUBROUTINE SCIPROC
!        SUBROUTINE WR_CONC ( SGRID, JDATE, JTIME, TSTEP )
!           IMPLICIT NONE
!           REAL,    INTENT( IN )     :: SGRID( :,:,:,: )
!           INTEGER, INTENT( IN )     :: JDATE, JTIME, TSTEP
!        END SUBROUTINE WR_CONC
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
!     CALL PAR_INIT( PNAME, NSPCSD, CLOCK, PAR_ERR )

      IF ( PAR_ERR .NE. 0 ) THEN
         XMSG = 'Error in PAR_INIT'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, 1 )
         END IF

      LOGDEV = INIT3 ()

C Initialize PCGRID

      IF ( .NOT. PCGRID_INIT () ) THEN
         XMSG = 'Failure defining horizontal domain'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2  )
         END IF
      
      CGRID => PCGRID( 1:MY_NCOLS,1:MY_NROWS,:,: )   ! required for PinG

C Initalize CONC definitions (in STD_CONC F90 module)

      CALL CONC_DEFN ()
!     C_NLYS = CONC_ELEV - CONC_BLEV + 1

C Get avg CONC definitions, species and layer pointers (in AVG_CONC F90 module)

      CALL A_CONC_DEFN ()
      A_NLYS = ACONC_ELEV - ACONC_BLEV + 1

C Initialize conc field: Copy IC's to CONC file as step 0
C Convention: the input file concentration units are always ppmV.

      CALL INITSCEN ( CGRID, STDATE, STTIME, TSTEP, NSTEPS )
      CALL CKSUMMER ( 'INITSCEN', CGRID, STDATE, STTIME )

C Verify input file header consistency with COORD.EXT and run duration

      CALL FLCHECK ( STDATE, STTIME, TSTEP( 1 ) )

      JDATE = STDATE; JTIME = STTIME

      IF ( LIPR .OR. LIRR ) CALL PA_INIT ( CGRID, STDATE, STTIME, TSTEP )

!     write( logdev,* ) 'driver - A_NLYS,N_ASPCS : ', mype, a_nlys, n_aspcs 

!     ALLOCATE ( SGRID( MY_NCOLS,MY_NROWS,C_NLYS,N_CSPCS ), STAT = ALLOCSTAT )
!     IF ( ALLOCSTAT .NE. 0 ) THEN
!        XMSG = 'SGRID memory allocation failed'
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!        END IF

!     SGRID = 0.0

      ALLOCATE ( AGRID( MY_NCOLS,MY_NROWS,A_NLYS,N_ASPCS ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'AGRID memory allocation failed'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      AGRID = 0.0

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
                     END DO
                  END DO
               END DO
            END DO

C science process sequence:
         
         DO IREP = 1, NREPS - 1

            CALL SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )

            DO V = 1, N_ASPCS
               S = AVG_CONC_MAP( V )
               L = 0
               DO K = ACONC_BLEV, ACONC_ELEV
                  L = L + 1
                  DO R = 1, MY_NROWS
                     DO C = 1, MY_NCOLS
                        AGRID( C,R,L,V ) = AGRID( C,R,L,V )
     &                                   + 2.0 * CGRID( C,R,K,S )
                        END DO
                     END DO
                  END DO
               END DO

            END DO

         CALL SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )

         DIVFAC = 0.5 / FLOAT( NREPS )

         DO V = 1, N_ASPCS
            S = AVG_CONC_MAP( V )
            L = 0
            DO K = ACONC_BLEV, ACONC_ELEV
               L = L + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     AGRID( C,R,L,V ) = DIVFAC * ( AGRID( C,R,L,V )
     &                                +            CGRID( C,R,K,S ) )
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

!        CALL WR_CONC ( SGRID, JDATE, JTIME, TSTEP( 1 ) )
         CALL WR_CONC ( JDATE, JTIME, TSTEP( 1 ) )

         CALL WR_ACONC ( AGRID, JDATE, JTIME, TSTEP( 1 ) )

!        CALL WR_CGRID ( CGRID, JDATE, JTIME, TSTEP( 1 ) )

         IF ( LIPR .OR. LIRR ) CALL PA_OUTPUT ( CGRID, JDATE, JTIME )

201      CONTINUE

!     CALL WR_CGRID ( CGRID, JDATE, JTIME, NSTEPS*TSTEP( 1 ) )
      CALL WR_CGRID ( CGRID, JDATE, JTIME, TSTEP( 1 ) )

C Clean up communications.

      CALL PAR_TERM ( PNAME, JDATE, JTIME, CLOCK )

      END
