
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_saprc99_ae5_ddm3d/hrsolver.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE HRSOLVER( JDATE, JTIME, C, R, L )


C**********************************************************************
C
C  FUNCTION: EBI solver
C
C  PRECONDITIONS: For the SAPRC99 family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  HRRATES, HRG1, HRG2, HRG3,
C                                     HRG4, HRPRODLOSS
C
C  REVISION HISTORY: Created by EBI solver program, Sept. 10, 2008
C                     Revised hrsolver.F by Golam Sarwar, December, 2007
C                             rearranged a statement suggested by CARB
C                             to reduce the chance of code crash; it does
C                             not change results
C**********************************************************************
      USE HRDATA

      IMPLICIT NONE

C..INCLUDES:
      INCLUDE SUBST_GC_SPC    ! Gas chem species names and MWs

C..ARGUMENTS:
      INTEGER JDATE           ! Current date (YYYYDDD)
      INTEGER JTIME           ! Current time (HHMMSS)
      INTEGER C, R, L         ! Cell col, row, lev

C..PARAMETERS:
      INTEGER, PARAMETER :: MXBKUPS = 5  ! Max no. of back-ups allowed
      INTEGER, PARAMETER :: STAT = 1     ! Status code

      REAL, PARAMETER :: EPSLON = 1.0E-30     ! Small number
      REAL, PARAMETER :: MAXPRED = 1.0E+03    ! Upper limit on predicted conc
      REAL, PARAMETER :: ZERO = 0.0               ! zero

C..EXTERNAL FUNCTIONS:
      INTEGER JUNIT


C..SAVED LOCAL VARIABLES:
      CHARACTER( 16 ), SAVE ::  PNAME = 'HRSOLVER'      ! Program name


C..SCRATCH LOCAL VARIABLES:

      CHARACTER( 132 ) :: MSG           ! Message text

      INTEGER CELLNO          ! Cell no. fo debug output
      INTEGER ITER            ! Loop index for Backward Euler iterations
      INTEGER S               ! Loop index for species
      INTEGER NEBI            ! Loop index for time steps
      INTEGER NINR            ! No. of inner time steps
      INTEGER N               ! Loop index
      INTEGER EBI             ! Loop index
      INTEGER NBKUPS          ! No. of times time step reduced

      LOGICAL LEBI_CONV             ! Flag for EBI convergence
      LOGICAL LEBISPFL( N_GC_SPC )  ! Flag for EBI species
      LOGICAL MXFL                  ! hit MAXPRED flag

      REAL DTC              ! Time step to take
      REAL FXDLOSS          ! Total loss due to negative stoichiometry
      REAL VARLOSS          ! Loss excluding negative stoichiometry


#ifdef hrdebug
      CHARACTER*8  NOTE       ! Convergence fail note

      INTEGER COL             ! Column to generate deboug output for
      INTEGER ROW             ! Row to generate deboug output for
      INTEGER LEV             ! Level to generate deboug output for
      INTEGER DBGOUT          ! Output unit for debu outpt

      LOGICAL LDEBUG          ! Debug output flag
      LOGICAL, SAVE  :: LOPEN = .FALSE.
#endif

#ifdef sens
      INTEGER N_EBI_MID    ! the midpoint ebi step; half of N_EBI_STEPS
      LOGICAL ODD_STEPS    ! true if N_EBI_STEPS is odd
      LOGICAL DDM_CHECK( N_GC_SPC)     ! check YCDDM set for all species
#endif


C**********************************************************************




c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
      COL = 0
      ROW = 0
      LEV = 0
      IF( C .EQ. COL .AND. R .EQ. ROW .AND. L .EQ. LEV ) THEN
c      IF( JTIME .EQ. 160000 ) THEN
         LDEBUG = .TRUE.
      ELSE
         LDEBUG = .FALSE.
      ENDIF

      IF( LDEBUG ) THEN
           IF( .NOT. LOPEN ) THEN
              DBGOUT = JUNIT()
              OPEN( UNIT = DBGOUT, FILE = 'debug.out' )
              LOPEN = .TRUE.
           ENDIF

           WRITE( DBGOUT, '( A, 2I4, I3, 1X, I7, 1X, I6 ) ' )
     &             'Debug output for col/row/lev/date/time:',
     &              C, R, L, JDATE, JTIME
           WRITE( DBGOUT, '( A, F7.2) ' )
     &             'EBI_TMSTEP = ', EBI_TMSTEP
           WRITE( DBGOUT, '( A )' ) 'Starting concs and rate constants'
           DO N = 1, N_SPEC
             WRITE( DBGOUT,  '( A, I3, 1X, A, 1X, 1PE13.5 )' )
     &                     'SP ',N, GC_SPC( N ), YC( N )
           ENDDO
           DO N = 1, N_RXNS
             WRITE( DBGOUT, '( A, I3, 1X, 1PE13.5 )' )
     &                     'RKI ', N, RKI( N )
           ENDDO
      ENDIF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++

#ifdef sens
      DO S = 1, N_SPEC
        DDM_CHECK( S ) = .FALSE.
      END DO
      N_EBI_MID = CEILING( N_EBI_STEPS / 2.0 )
      IF ( MOD( N_EBI_STEPS, 2 ) .EQ. 1 ) THEN
        ODD_STEPS = .TRUE.
      ELSE
        ODD_STEPS = .FALSE.
      END IF
#endif


      N_EBI_IT = 0

      DO 3000 NEBI = 1, N_EBI_STEPS    ! outer EBI time-tep loop

         DTC = EBI_TMSTEP
         NBKUPS = 0
         N_INR_STEPS = 1

 100     CONTINUE                        !  Restart location

         DO 2000 NINR = 1, N_INR_STEPS   ! No. of time steps for back-up

            DO S = 1, N_SPEC             ! Set ICs for EBI iterations
               YC0( S ) = YC( S )
            ENDDO

            DO 1000 ITER = 1, NEBITER    ! EBI iteration loop

               N_EBI_IT = N_EBI_IT + 1

               CALL HRRATES

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
               IF( LDEBUG ) THEN
                  WRITE( DBGOUT, '( A, I5 )' ) 'ITER NO ', ITER
                  WRITE( DBGOUT, '( A, F12.5 )' )
     &               ' DTC=', DTC

                  IF( ITER .EQ. 1 ) THEN
                     WRITE( DBGOUT, '( A )' ) 'Starting reaction rates'
                     DO N = 1, N_RXNS
                        WRITE( DBGOUT, '( A, I3, 1X, 1PE13.5 )' )
     &                        'RXRAT ', N, RXRAT( N )
                     ENDDO
                  ENDIF
               ENDIF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Call routines to compute concentrations of groups 1-4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               CALL HRG1( DTC )

               CALL HRG2( DTC )

               CALL HRG3( DTC )

               CALL HRG4( DTC )


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do the Euler backward method
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               CALL HRPRODLOSS

               DO EBI = 1, N_EBISP
                  S = EBISP( EBI )

                  YCP( S ) = YC( S )*(( YC0( S ) + PROD( S ) * DTC ) /
     &                               ( YC ( S ) + LOSS( S ) * DTC ) )


               ENDDO




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check for convergence
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               LEBI_CONV = .TRUE.
               MXFL = .FALSE.
               DO S = 1, N_SPEC
                  LEBISPFL( S ) = .FALSE.
                  YCP( S ) = MAX( EPSLON, YCP( S ) )
                  IF( ABS( YC( S ) - YCP( S ) ) .GT. RTOL( S ) *
     &               ( YC( S ) + YCP( S ) ) ) THEN
                     LEBI_CONV = .FALSE.
                     LEBISPFL( S ) = .TRUE.
                  ENDIF
c..if predictions growing too large, treat as a convergence failure
                  IF( YCP( S ) .GT. MAXPRED ) then
                     MXFL = .TRUE.
                     GO TO 1010
                  END IF
#ifdef sens
C YCDDM is the value at the middle time step 
                  IF ( N_EBI_MID .EQ. NEBI ) THEN
                    DDM_CHECK ( S ) = .TRUE.
                    IF ( ODD_STEPS ) THEN
                      YCDDM( S ) = 0.5 * ( YC( S ) + YCP( S ) )
                    ELSE
                      YCDDM( S ) = YCP( S )
                    END IF

                    IF ( YCDDM( S ) .LT. 1.e-25) YCDDM( S ) = 0.

                  END IF
#endif

                  YC( S ) = YCP( S )
               ENDDO

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef hrdebug
               IF( LDEBUG ) THEN
                  WRITE( DBGOUT, '( A, I5 )' ) 'Concs after ITER= ', ITER
                  DO S = 1, N_SPEC

                     IF( LEBISPFL( S ) ) THEN
                        NOTE = 'CONV FAIL'
                     ELSE
                        NOTE = '         '
                     ENDIF

                     WRITE( DBGOUT, '( I3, 1X, A, 1PE13.5, 1X, A )' )
     &                            S, GC_SPC( S ), YC( S ), NOTE
                  ENDDO
                  IF( LEBI_CONV ) WRITE( DBGOUT, '( A )' )
     &                 '****Convergence achieved'
               ENDIF
#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++


               IF( LEBI_CONV ) GO TO 2000

 1000       CONTINUE

c...Convergence failure section; cut the inner time step in half &
c.....start inner loop over unless max backups exceeded

 1010       CONTINUE

            NBKUPS = NBKUPS + 1

            IF( NBKUPS .LE. MXBKUPS ) THEN

               IF ( MXFL ) THEN
                  WRITE( LOGDEV, 92010 ) C, R, L, TRIM( GC_SPC( S ) ), NBKUPS
               ELSE
                  WRITE( LOGDEV, 92000 ) C, R, L, NBKUPS
               END IF

               DO S = 1, N_SPEC
                  YC( S ) = YC0( s )
               ENDDO

               DTC = 0.5 * DTC

               N_INR_STEPS = 2 ** NBKUPS

               GO TO 100

            ELSE

               WRITE( LOGDEV, 92040 ) C, R, L

               WRITE( LOGDEV, 92060 )
               DO S = 1, N_SPEC
                  IF( LEBISPFL( S ) ) WRITE( LOGDEV, 92080 ) GC_SPC( S )
               ENDDO

               MSG='ERROR: Stopping because of EBI convergence failures'
               CALL M3EXIT( PNAME, JDATE, JTIME, MSG, STAT )

            ENDIF

 2000    CONTINUE

 3000 CONTINUE

#ifdef sens
      DO S = 1, N_SPEC
        IF ( .NOT. DDM_CHECK( S ) ) THEN
          MSG = 'ERROR: HRSOLVER failed to assign YCDDM to ' //
     &          GC_SPC( S )
          CALL M3EXIT( PNAME, JDATE, JTIME, MSG, STAT )
        END IF
      END DO
#endif


      RETURN


92000 FORMAT( 'WARNING: EBI Euler convergence failure' /
     &        '         Reducing EBI time step because of ',
     &         'convergence failure for ' /
     &        '         Cell (', I3, ', ', I3, ', ', I3, ')' ,
     &        '  Back-up number', I2 )

92010 FORMAT( 'WARNING: EBI Euler convergence failure' /
     &        '         Reducing EBI time step because of ',
     &         'MAXPRED convergence failure for ' /
     &        '         Cell (', I3, ', ', I3, ', ', I3, ')' ,
     &        ' and species ', A,
     &        '  Back-up number', I2 )

92040 FORMAT( 'ERROR: Max number of EBI time step reductions exceeded'
     &      / '      Convergence failure for cell (', I3, ', ', I3,
     &                ', ', I3, ')' )

92060 FORMAT( '      Convergence failure for the following species:' )

92080 FORMAT( 10X, A )

      END
