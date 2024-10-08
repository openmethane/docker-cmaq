
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/procan/pa/pa_init.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $ 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)pa_init.F	1.1 /project/mod3/CMAQ/src/procan/pa/SCCS/s.pa_init.F 07 Jul 1997 10:20:52

      SUBROUTINE PA_INIT( CGRID, SDATE, STIME, TSTEP )
 
C**********************************************************************
C
C  FUNCTION: To create the Process Analysis ouput files and initialize
C            Process Analysis arrays           
C
C  PRECONDITIONS: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: SUBST_PA_INDEX, SUBST_IRR_INDEX,
C                                    PA_MKHDR, SUBST_4D_DATA_COPY
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, July, 1996
C                    allow env var for file names Jeff, Dec, 1996
C                    Changed CGOUT dimension parameter from NSPCSD to
C                       NSPCS -- Jerry Gipson, May, 1997.
C                    Modified May, 1997 by Jerry Gipson to be consistent
C                      with beta CTM
C                    Modified June 1997 by Jerry Gipson to change CGOUT
C                      to CGRID       
C                    Modified Sept, 1997 by Jerry Gipson to be consistent
C                      with targeted CTM
C                    Modified Jun, 1998 by Jerry Gipson to add PA domain
C                      error check
C                    Modified 1/19/99 by David Wong at LM:
C                      -- add four include files because of new PA_CMN.EXT
C                      -- modify PA_INDEX function call parameter list
C                      -- add DATA_COPY function call to redistribute PA grid
C                    Modified 2/26/99 by David Wong at LM:
C                      -- remove SUBST_AE_SPC, SUBST_NR_SPC, SUBST_TR_SPC,
C                         three .EXT files
C                      -- call SUBST_IRR_INDEX to determine row, column, and
C                         level loop indexes which are associated with IRR
C                         calculations
C                      -- renamed DATA_COPY function name
C                      -- use ifdef statement to distinguish parallel
C                         implementation of IRR calculation which does not 
C                         shift to the origin
C                    Modified 4/12/00 by Jerry Gipson to correct units for
C                      aerosol IPRs writen to IO/API file
C                    Modified 1 Jul 2000 by Jeff Young to allow for file
C                      updating - added PA_MKHDR
C                    Modified 22 Nov 00 by J.Young: Dave Wong's f90 stenex 
C                      SUBGRID_INDEX and DATA_COPY
C    6 Oct 01 J.Young: dyn alloc - generic CGRID declaration; Use PAGRD_DEFN,
C                     which uses HGRD_DEFN; new pario
C                    Modified 10 Oct 01 by David Wong
C                      -- let other processors open PA output files after
C                         processor has created them.
C                      -- used a new formula to compute number of IPR and IRR 
C                         output files
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module (GRID_CONF)
C   16 Jan 06 J.Young: move PAGRD_INIT to par_init
C**********************************************************************

      USE GRID_CONF             ! horizontal & vertical domain specifications
!     USE HGRD_DEFN             ! horizontal domain specifications
!     USE VGRD_DEFN             ! vertical layer specifications
      USE PAGRD_DEFN            ! PA horiz domain specs
      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex
!     USE SUBST_DATA_COPY_MODULE   ! stenex

      IMPLICIT NONE 

C..INCLUDES:
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
!     INCLUDE SUBST_BLKPRM      ! cell blocking parameters

      INCLUDE SUBST_PACTL_ID    ! PA control parameters
      INCLUDE SUBST_PACMN_ID    ! PA common block definitions
      INCLUDE SUBST_PADAT_ID    ! PA data statements
     

C..ARGUMENTS:
      INTEGER SDATE    !  starting date,    format YYYYDDD
      INTEGER STIME    !  starting time,    format HHMMSS
      INTEGER TSTEP( 2 )   ! time step vector (HHMMSS)
                           ! TSTEP(1) = local output step
                           ! TSTEP(2) = sciproc sync. step (chem)

!#ifdef parallel
!     INCLUDE SUBST_GC_SPC      ! gas chemistry species table (dim./name/molwt)
!     INCLUDE SUBST_AE_SPC      ! aerosol species table
!     INCLUDE SUBST_NR_SPC      ! non-reactive species table
!     INCLUDE SUBST_TR_SPC      ! tracer species table

!     INTEGER, PARAMETER :: ALLSPC = N_GC_SPCD + N_AE_SPC + N_NR_SPC + N_TR_SPC
!     REAL :: CGRID( NCOLS,NROWS,NLAYS,ALLSPC )
!     REAL :: CGRID( :,:,:,: )
      REAL, POINTER :: CGRID( :,:,:,: )
!#else
!     REAL CGRID( NCOLS,NROWS,NLAYS,* )
!#endif

C..PARAMETERS:

C..EXTERNAL FUNCTIONS: None

C..SAVED LOCAL VARIABLES: None

C..SCRATCH LOCAL VARIABLES:
      CHARACTER( 16 ) ::  ENV_DFLT  ! Environment variable default value
      CHARACTER( 80 ) ::  ENV_DESC  ! Environment variab;e description
      CHARACTER( 16 ) ::  OUTFNAME  ! Assigned IPR or IRR output file name
      CHARACTER( 16 ) ::  PNAME = 'PA_INIT' ! Routine Name
      CHARACTER( 256 ) :: RET_VAL   ! Returned value of environment variable

      LOGICAL LSTOP     ! Flag to stop because a PA file not assigned
 
      INTEGER C         ! Loop index for columns
      INTEGER ICG       ! Index for species in cgout array
      INTEGER L         ! Loop index for layers
      INTEGER NIRR      ! Loop index for IRR outputs
      INTEGER NFL       ! Loop index for no. of files
      INTEGER NOUT      ! Loop index for process monitoring outputs
      INTEGER NS        ! Loop index for number cgout concs saved
      INTEGER NUMFLS    ! Number of PA or IRR output files needed
      INTEGER PC        ! Index for PA output column
      INTEGER PL        ! Index for PA output level
      INTEGER PR        ! Index for PA output row
      INTEGER R         ! Loop index for rows
      INTEGER STATUS    ! Status Code

C**********************************************************************

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check process analysis domain for consistency with CTM domain
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C move to PAGRD_DEFN.F module

!     IF( BEGCOL .LT. 1 .OR. ENDCOL .GT. GL_NCOLS .OR.
!    &    BEGROW .LT. 1 .OR. ENDROW .GT. GL_NROWS .OR.
!    &    BEGLEV .LT. 1 .OR. ENDLEV .GT. NLAYS ) THEN
!         CALL M3EXIT( PNAME, SDATE, STIME,
!    &        'Process Analysis domain inconsistent with CTM domain',
!    &         2 )
!     END IF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Get local processor's computational domain bounds
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C move to PAGRD_DEFN.F module

C Determine IPR subgrid mapping with data redistribution
!     CALL SUBST_SUBGRID_INDEX
!    &     ( BEGCOL, ENDCOL, BEGROW, ENDROW, BEGLEV, ENDLEV,
!    &       MY_BEGCOL, MY_ENDCOL, MY_BEGROW, MY_ENDROW,
!    &       MY_BEGLEV, MY_ENDLEV, MY_PACOLS, MY_PAROWS )

C Determine IRR subgrid mapping without data redistribution
!     CALL SUBST_SUBGRID_INDEX
!    &     ( BEGCOL, ENDCOL, BEGROW, ENDROW, BEGLEV, ENDLEV,
!    &       MY_IRR_BEGCOL, MY_IRR_ENDCOL, MY_IRR_BEGROW,
!    &       MY_IRR_ENDROW, MY_IRR_BEGLEV, MY_IRR_ENDLEV, 'N' )

C move to par_init

!     IF ( .NOT. PAGRD_INIT() ) THEN
!        CALL M3EXIT( PNAME, SDATE, STIME,
!    &                'Could not initialize PA horiz grid defintions',
!    &                2 )
!        END IF

      LSTOP = .FALSE.         

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Open the Integrated Process Rates output file(s) if needed
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF ( LIPR ) THEN

C..compute the number of files needed

         NUMFLS = (NIPRVAR - 1) / MXVARS3 + 1

!        IF( MOD( NIPRVAR, MXVARS3 ) .EQ. 0 ) THEN
!           NUMFLS = NIPRVAR / MXVARS3
!        ELSE
!           NUMFLS = ( NIPRVAR / MXVARS3 ) + 1
!        END IF

C..open each one 
         DO NFL = 1, NUMFLS

            IF( NFL .EQ. 1 ) THEN
               OUTFNAME = CTM_IPR_1
            ELSEIF( NFL .EQ. 2 ) THEN
               OUTFNAME = CTM_IPR_2
            ELSEIF( NFL .EQ. 3 ) THEN
               OUTFNAME = CTM_IPR_3
            ELSEIF( NFL .GT. 3 ) THEN
               CALL M3EXIT(PNAME, SDATE, STIME,
     &            'Maximum number of IPR output files exceeded' , 2 )
            END IF

            ENV_DESC = 'IPR Output file ' // OUTFNAME 
            ENV_DFLT = ' '       
            CALL ENVSTR( OUTFNAME, ENV_DESC, ENV_DFLT, RET_VAL, STATUS)
            IF( STATUS .NE. 0 ) CALL M3EXIT( PNAME, SDATE, STIME,
     &             'IPR output file ' // OUTFNAME // ' not assigned', 2 )

!           CALL SUBST_BARRIER

C..try to open existing file for update
            IF ( .NOT. OPEN3( OUTFNAME, FSRDWR3, PNAME ) ) THEN
               CALL M3MESG( 'Could not open ' // OUTFNAME //
     &                      ' file for update - try to open new' )

               IF ( MYPE .EQ. 0 ) THEN   ! open new

                  CALL PA_MKHDR ( 'IPR', NFL, OUTFNAME, SDATE, STIME, TSTEP )

                  IF ( OPEN3( OUTFNAME, FSNEW3, PNAME ) ) THEN
                     CALL M3MESG( 'Opened Integrated Reaction Rate Output'
     &                           // ' File ' //  OUTFNAME )
                  ELSE
                     CALL M3EXIT( PNAME, SDATE, STIME,
     &                   'Could not open Integrated Process Rate Output File ' 
     &                   // OUTFNAME, 2 )
                  END IF
               END IF

!           CALL SUBST_BARRIER

!              IF ( MYPE .NE. 0 ) THEN
!                 IF ( .NOT. OPEN3( OUTFNAME, FSRDWR3, PNAME ) ) THEN
!                    CALL M3EXIT( PNAME, SDATE, STIME,
!    &                   'Could not open Integrated Process Rate Output File ' 
!    &                   // OUTFNAME, 2 )
!                 END IF
!              END IF

            END IF

         END DO

      CALL SUBST_BARRIER

      END IF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Open the Integrated Reaction Rates output file(s) if needed
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF ( LIRR ) THEN

C..compute the number of files needed

         NUMFLS = (NIRRVAR - 1) / MXVARS3 + 1

!        IF( MOD( NIRRVAR, MXVARS3 ) .EQ. 0 ) THEN
!           NUMFLS = NIRRVAR / MXVARS3
!        ELSE
!           NUMFLS = ( NIRRVAR / MXVARS3 ) + 1
!        END IF

C..open each one
         DO NFL = 1, NUMFLS

            IF( NFL .EQ. 1 ) THEN
               OUTFNAME = CTM_IRR_1
            ELSEIF( NFL .EQ. 2 ) THEN
               OUTFNAME = CTM_IRR_2
            ELSEIF( NFL .EQ. 3 ) THEN
               OUTFNAME = CTM_IRR_3
            ELSEIF( NFL .GT. 3 ) THEN
               CALL M3EXIT(PNAME, SDATE, STIME,
     &            'Maximum number of IRR output files exceeded' , 2 )
            END IF

            ENV_DESC = 'IRR Output file ' // OUTFNAME 
            ENV_DFLT = ' '       
            CALL ENVSTR( OUTFNAME, ENV_DESC, ENV_DFLT, RET_VAL, STATUS)
            IF( STATUS .NE. 0 ) CALL M3EXIT( PNAME, SDATE, STIME,
     &             'IRR output file ' // OUTFNAME // ' not assigned', 2 )

!           CALL SUBST_BARRIER

C..try to open existing file for update
            IF ( .NOT. OPEN3( OUTFNAME, FSRDWR3, PNAME ) ) THEN
               CALL M3MESG( 'Could not open ' // OUTFNAME //
     &                      ' file for update - try to open new' )

               IF ( MYPE .EQ. 0 ) THEN   ! open new

                  CALL PA_MKHDR ( 'IRR', NFL, OUTFNAME, SDATE, STIME, TSTEP )

                  IF( OPEN3( OUTFNAME, FSNEW3, PNAME ) ) THEN
                     CALL M3MESG( 'Opened Integrated Reaction Rate Output'
     &                           // ' File ' //  OUTFNAME )
                  ELSE
                     CALL M3EXIT( PNAME, SDATE, STIME,
     &                    'Could not open Integrated Reaction Rate Output'
     &                    // ' File ' // OUTFNAME, 2 )
                  END IF
               END IF

!              CALL SUBST_BARRIER

!              IF ( MYPE .NE. 0 ) THEN
!                 IF ( .NOT. OPEN3( OUTFNAME, FSRDWR3, PNAME ) ) THEN
!                    CALL M3EXIT( PNAME, SDATE, STIME,
!    &                    'Could not open Integrated Reaction Rate Output'
!    &                    // ' File ' // OUTFNAME, 2 )
!                 END IF
!              END IF

            END IF

         END DO

      CALL SUBST_BARRIER

      END IF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Initialize arrays
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO  NOUT = 1, NIPRVAR
         DO L = 1, PALEVS
            DO R = 1, MY_PAROWS
               DO C = 1, MY_PACOLS
                  DELC( C, R, L, NOUT ) = 0.0
               END DO
            END DO
         END DO
      END DO

      DO NS = 1, NCSAVE
         ICG = SV2GRD( NS )
#ifdef parallel
         CALL SUBST_DATA_COPY (CGRID, CSAV, ICG, NS)
#else
         DO L = MY_BEGLEV, MY_ENDLEV
            PL = L - MY_BEGLEV + 1
            DO R = MY_BEGROW, MY_ENDROW
               PR = R - MY_BEGROW + 1
               DO C = MY_BEGCOL, MY_ENDCOL
                  PC = C - MY_BEGCOL + 1
                  CSAV( PC, PR, PL, NS ) = CGRID ( C, R, L, ICG ) 
               END DO
           END DO
         END DO
#endif
      END DO

      IF ( LIRR ) THEN

         DO NIRR = 1, NIRRVAR

#ifdef parallel
            DO L = MY_IRR_BEGLEV, MY_IRR_ENDLEV
               DO R = MY_IRR_BEGROW, MY_IRR_ENDROW
                  DO C = MY_IRR_BEGCOL, MY_IRR_ENDCOL
                     IRROUT( C, R, L, NIRR ) = 0.0
                     END DO
               END DO
            END DO
#else
            DO L = 1, PALEVS
               DO R = 1, MY_PAROWS
                  DO C = 1, MY_PACOLS
                     IRROUT( C, R, L, NIRR ) = 0.0
                  END DO
               END DO
            END DO
#endif

         END DO

      END IF

      END
