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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ros3/rbsolver.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE RBSOLVER( LOGDEV, JDATE, JTIME, CHEMSTEP, NCSP,
     &                      LIRRFLAG, NIRRCLS, IRRCELL )


C**********************************************************************
C
C  FUNCTION: ROS3 solver 
C
C  PRECONDITIONS: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: 
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, August, 2004
C                    Based on the solver described by Sandu et al
C                    ( Atm. Env., Vol. 31, No. 20, 1997 ) and included
C                    in the Kinetic PreProcessor-KPP ( see for example 
C                    Sandu et al., At, Env., Vol. 37, 5097-5114, 
C                    2003). This code also incorporates efficiency
C                    concepts originally developed by M. Jacobson
C                    for SMVGEAR (Atm. Env., Vol 28, No 2, 1994)
C                                          
C                    31 Jan 05 J.Young: get BLKSIZE from dyn alloc horizontal
C                    & vertical domain specifications module (GRID_CONF)
C                    29 Jul 05     WTH: Added IF blocks that call degrade 
C                                       routines if CALL_DEG is true, i.e.,
C                                       if MECHNAME contains 'TX' substring.
C**********************************************************************
      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE RBDATA                  ! ROS3 solver data

      IMPLICIT NONE 

C..INCLUDES:
      INCLUDE SUBST_GC_SPC        ! Gas chem species names and MWs


C..ARGUMENTS:
      INTEGER LOGDEV              ! Unit number of output log
      INTEGER JDATE               ! Current date (YYYYDDD)
      INTEGER JTIME               ! Current time (HHMMSS)
      INTEGER NCSP                ! Index of chem mech to use
                                  ! 1=gas/day, 2=gas/night
      INTEGER NIRRCLS             ! No. of cells in block for IRR
      INTEGER IRRCELL( BLKSIZE )  ! Cell No. of an IRR cell

      LOGICAL LIRRFLAG            ! Flag for IRR calculations

      REAL    CHEMSTEP            ! Chem integration interval (min)


C..PARAMETERS:

c..ROS3 solver parameters - from KPP
      REAL( 8 ), PARAMETER :: GAM =  0.43586652150845899941601945119356D+00
      REAL( 8 ), PARAMETER :: C21 = -0.10156171083877702091975600115545D+01
      REAL( 8 ), PARAMETER :: C31 =  0.40759956452537699824805835358067D+01
      REAL( 8 ), PARAMETER :: C32 =  0.92076794298330791242156818474003D+01
      REAL( 8 ), PARAMETER :: B1  =  0.10000000000000000000000000000000D+01
      REAL( 8 ), PARAMETER :: B2  =  0.61697947043828245592553615689730D+01
      REAL( 8 ), PARAMETER :: B3  = -0.42772256543218573326238373806514D+00
      REAL( 8 ), PARAMETER :: D1  =  0.50000000000000000000000000000000D+00
      REAL( 8 ), PARAMETER :: D2  = -0.29079558716805469821718236208017D+01
      REAL( 8 ), PARAMETER :: D3  =  0.22354069897811569627360909276199D+00
      REAL( 8 ), PARAMETER :: A21 =  1.0D+00
      REAL( 8 ), PARAMETER :: A31 =  1.0D+00
      REAL( 8 ), PARAMETER :: A32 =  0.0D+00
      REAL( 8 ), PARAMETER :: G1  =  0.43586652150845899941601945119356D+00
      REAL( 8 ), PARAMETER :: G2  =  0.24291996454816804366592249683314D+00
      REAL( 8 ), PARAMETER :: G3  =  0.21851380027664058511513169485832D+01
      REAL( 8 ), PARAMETER :: GROW = 1.0D+00 / 3.00D+00

      REAL( 8 ), PARAMETER :: DTSTART  =  0.5D-01  ! Starting time step (min)
      REAL( 8 ), PARAMETER :: DTMIN    =  1.0D-08  ! Min time step
      REAL( 8 ), PARAMETER :: DTMAX    =  1.0D+01  ! Max time step
      REAL( 8 ), PARAMETER :: UROUND   =  1.0D-18  ! Roundoff parameter

      REAL( 8 ), PARAMETER :: FACMAX   =  1.0D+01  ! Max time step factor
      REAL( 8 ), PARAMETER :: FACMIN   =  1.0D-01  ! Min time step factor
      REAL( 8 ), PARAMETER :: FACONE   =  1.0D+00  ! Time step fac of 1.0

      REAL( 8 ), PARAMETER :: CONMIN   =  1.0D-30  ! Min conc


C..EXTERNAL FUNCTIONS:
      INTEGER JUNIT                                ! Gets unit number of file


C..SAVED LOCAL VARIABLES:
      CHARACTER( 16 ), SAVE :: PNAME = 'RBSOLVER'  ! Program name

     
C..SCRATCH LOCAL VARIABLES:
      INTEGER I, J, N, JSPC       ! Loop indices
      INTEGER IDIAGBEG            ! Index of diagonal start in Jac array 
      INTEGER OFFSET              ! Cell offset for blcock
      INTEGER NCELL               ! Cell loop index
      INTEGER DUMMY( BLKSIZE )    ! Dummy array for IRR call
      INTEGER NRX                 ! Loop index for reactions
      INTEGER ISPOLD              ! Species index for old array order
      INTEGER NCALL_DEGRADE       ! WTH

      LOGICAL LPASS               ! Flag for convergence achieved

      REAL DTCELL( BLKSIZE )             ! Time step for each cell for IRR
      REAL CIRR(   BLKSIZE, N_GC_SPC )   ! Species concs for IRR analysis
      REAL RKIRR(  BLKSIZE, N_RXNS )     ! Rate constants for IRR analysis

      REAL( 8 )  D
      REAL( 8 )  DT               ! Time step
      REAL( 8 )  DTINV            ! Inverse of time step
      REAL( 8 )  DTFAC            ! Time step scale factor
      REAL( 8 )  GDTINV           ! Inverse of gamma x time step
      REAL( 8 )  TNOW             ! Elapsed time at start of integration step
      REAL( 8 )  TEND             ! Elapsed time at end of integration step
      REAL( 8 )  DT_DEGRADE       ! WTH: Time step for degradation routine
      REAL( 8 )  YTOL             ! Species tolerance
      REAL( 8 )  RNSPEC           ! Recipricol of # of species
      REAL( 8 )  X1, X2           ! Temp ROS3 variables
      REAL( 8 )  ERRYMAX          ! Cell/species stiffness estimate
      REAL( 8 )  YLOWEPS          ! Tolerance ratio used in stiffness calc
      REAL( 8 )  MAXERR           ! Max of cell error estimates

      REAL( 8 )  ERR(  BLKSIZE )           ! Error est for each cell
      REAL( 8 )  YDOT( BLKSIZE, N_SPEC )   ! dc/dt array

c.....ROS3 intermediate variables
      REAL( 8 )  K1( BLKSIZE, N_SPEC )    
      REAL( 8 )  K2( BLKSIZE, N_SPEC )
      REAL( 8 )  K3( BLKSIZE, N_SPEC )
      REAL( 8 )  K4( BLKSIZE, N_SPEC )

      REAL( 8 )  YP( BLKSIZE, N_SPEC )     ! Predicted conc

#ifdef rbdebug

      INTEGER COL                          ! Column for debug output 
      INTEGER ROW                          ! Row for debug output
      INTEGER LEV                          ! Level for debug output
      INTEGER DBGOUT                       ! Output unit for debug output

      LOGICAL LDEBUG                       ! Debug output flag
      LOGICAL, SAVE  :: LOPEN = .FALSE.    ! Flag for debug file opened

#endif


C**********************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The s-stage Rosenbrock method solves the following equations
C                        s
C         Y(n+1) = Yn + SUM( Bi * Ki) 
C                       i=1
C     where
C                          i-1                          i
C         Ki = H * F[ Yn + SUM( Aij * Kj ) ] + H * J * SUM( GAMMAij * Kj ) 
C                          j=1                         j=1
C
C See Sandu et al. for details and the values of Bi, GAMMAij, Aij, etc.
C
C For computational efficiency, the equations are re-arranged as  
C follows (e.g., see Press, Numerical Recipes, Chap 16 on ODEs )
C
C     Gi = SUM( GAMMAij * Kj ) + GAMMA * Ki  i=1,...s
C
C    [ I / GAMMA H - J ] G1 = F[ Yn ]
C    [ I / GAMMA H - J ] G2 = F[ Yn + A21 * G1 ] + ( C21 * G1 ) / H
C    [ I / GAMMA H - J ] G3 = F[ Yn + A31 * G1 + A32 * G2 ] + 
C                              ( C31 * G1 + C32 G2 ) / H   
C 
C The code below sequentially calculates the Gi, and then computes
C Y(n+1) via the first formula.  Note that the parameter values will be
C different form those shown in Sandu et al. because the code computes
C Gi instead of Ki
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      RNSPEC = 1.0 / FLOAT( N_SPEC )


c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef rbdebug 

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
     &             'CHEMTMSTEP = ', CHEMSTEP
           WRITE( DBGOUT, '( A )' ) 'Starting concs and rate constants'
           DO N = 1, N_SPEC
             WRITE( DBGOUT,  '( A, I3, 1X, A, 1X, 1PE13.5 )' )
     &                     'SP ',N, GC_SPC( N ), Y( N )
           ENDDO
           DO N = 1, N_RXNS
             WRITE( DBGOUT, '( A, I3, 1X, 1PE13.5 )' )
     &                     'RKI ', N, RKI( N )
           ENDDO
      ENDIF

#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++


      ISCHAN = ISCHANG( NCS )

      IDIAGBEG = IARRAY( NCSP) - ISCHAN + 1

      DT = MAX( DTMIN, DTSTART )

      TNOW = 0.0D+00

      NCALL_DEGRADE = 0

      LPASS = .FALSE.

      OFFSET = BLKCNO( BLKID )

      IF( LIRRFLAG ) THEN

         DO NRX = 1, N_RXNS
            DO NCELL = 1, NIRRCLS
               RKIRR( NCELL, NRX ) = RKI( IRRCELL( NCELL ), NRX )
            ENDDO
         ENDDO

         DO NCELL = 1, NUMCELLS
            DTCELL( NCELL ) = 0.0
         ENDDO
         DO I = 1, ISCHAN
            ISPOLD = INEW2OLD( I, NCS )
            DO NCELL = 1, NIRRCLS
               CIRR( NCELL, ISPOLD ) = Y( IRRCELL( NCELL ), I )
            ENDDO
         ENDDO
         CALL PA_IRR ( .TRUE., .FALSE., RKIRR, CIRR, DTCELL, 
     &                    NIRRCLS, DUMMY )

      ENDIF


  100 CONTINUE


      CALL RBFEVAL( NCSP, Y, YDOT )


      IF( LORDERING ) THEN


         DO JSPC = 1, ISCHAN  
            YLOWEPS = ATOL( JSPC ) / ( MIN( RTOL( JSPC ), 0.003D+00 ) )
            DO NCELL = 1, NUMCELLS
               ERRYMAX  = YDOT( NCELL, JSPC ) / 
     &                   ( Y( NCELL, JSPC ) + YLOWEPS )
               ERRMX2( OFFSET + NCELL ) = ERRMX2( OFFSET + NCELL ) +
     &                                    ERRYMAX * ERRYMAX
            ENDDO
         ENDDO

         RETURN

      ENDIF

      TEND = TNOW + DT 

      IF( TEND .GT. CHEMSTEP ) THEN
         DT = CHEMSTEP - TNOW
         TEND = CHEMSTEP
      ENDIF    

      DTINV = 1.0D+00 / DT

      GDTINV = DTINV / GAM
      
      CALL RBJACOB( NCSP, Y )

      DO J = IDIAGBEG, IARRAY( NCSP )
         DO NCELL = 1, NUMCELLS
           CC2( NCELL, J ) = CC2( NCELL, J ) - GDTINV
         ENDDO
      ENDDO


      CALL RBDECOMP( NCSP )

c..stage 1
      DO N = 1, N_SPEC
         DO NCELL = 1, NUMCELLS
            K1( NCELL, N ) = -YDOT( NCELL, N )
         ENDDO
      ENDDO      

      CALL RBSOLVE( NCSP, K1 )

c..stage 2
      DO N = 1, N_SPEC
         DO NCELL = 1, NUMCELLS
            YP( NCELL, N ) = Y( NCELL, N ) + K1( NCELL, N )
         ENDDO
      ENDDO      

      CALL RBFEVAL( NCSP, YP, YDOT )


      X1 = C21 * DTINV
      DO N = 1, N_SPEC
         DO NCELL = 1, NUMCELLS
            K2( NCELL, N ) = -YDOT( NCELL, N ) - X1 * K1( NCELL, N )
         ENDDO
      ENDDO      

      CALL RBSOLVE( NCSP, K2 )


c..stage 3
      X1 = C31 * DTINV
      X2 = C32 * DTINV
      DO N = 1, N_SPEC
         DO NCELL = 1, NUMCELLS
            K3( NCELL, N ) = -YDOT( NCELL, N ) - X1 * K1( NCELL, N ) - 
     &                       X2 * K2( NCELL, N )
         ENDDO
      ENDDO

      CALL RBSOLVE( NCSP, K3 )


c..final solution
      DO N = 1, N_SPEC
         DO NCELL = 1, NUMCELLS
            YP( NCELL, N ) = Y( NCELL, N ) + B1 * K1( NCELL, N ) + 
     &                       B2 * K2( NCELL, N ) + B3 * K3( NCELL, N )
         ENDDO
      ENDDO


c..Estimate error
      ERR = 0.0D+00
      DO N = 1, N_SPEC
         DO NCELL = 1, NUMCELLS
            YTOL = ATOL( N ) + RTOL( N ) * ABS( YP( NCELL, N ) )
            ERR( NCELL ) = ERR( NCELL ) + ( ( D1 * K1( NCELL, N ) + 
     &          D2 * K2( NCELL, N ) + D3 * K3( NCELL, N ) ) / YTOL ) ** 2
        ENDDO
      ENDDO


      MAXERR = 0.0
      DO NCELL = 1, NUMCELLS
         MAXERR = MAX( MAXERR, UROUND, SQRT( ERR( NCELL ) * RNSPEC ) )
      ENDDO

      DTFAC = 0.9D+00 / MAXERR ** GROW

      IF( MAXERR .LE. 1.0D+00 ) THEN

#ifdef rbstats
         NSTEPS = NSTEPS + 1
#endif


         DO N = 1, N_GC_SPC
            DO NCELL = 1, NUMCELLS
               Y( NCELL, N ) = MAX( YP( NCELL, N ), CONMIN )
            ENDDO
         ENDDO

         IF( LIRRFLAG ) THEN
            DO NCELL = 1, NUMCELLS
               DTCELL( NCELL ) = DT
            ENDDO
            DO I = 1, ISCHAN
               ISPOLD = INEW2OLD( I, NCS )
               DO NCELL = 1, NIRRCLS
                  CIRR( NCELL, ISPOLD ) = Y( IRRCELL( NCELL ), I )
               ENDDO
            ENDDO
            CALL PA_IRR ( .FALSE., .FALSE., RKIRR, CIRR, DTCELL, 
     &                       NIRRCLS, DUMMY )
         ENDIF


         IF( CALL_DEG ) THEN  !:WTH applying degradation algorithm

            DT_DEGRADE = 60.0D0 * ( TEND - TNOW )

            DO I = 1, ISCHAN
               ISPOLD = INEW2OLD( I, NCS )
               DO NCELL = 1, NUMCELLS
                  Y_DEGRADE( NCELL, ISPOLD ) = Y( NCELL, I )
               ENDDO
            ENDDO

            NCALL_DEGRADE = NCALL_DEGRADE + 1

            CALL DEGRADE( Y_DEGRADE, DT_DEGRADE, JDATE, JTIME, BLKID )

         ENDIF            !:WTH

         TNOW = TEND

         IF( LPASS ) THEN
            DTFAC = MAX( FACMIN, MIN( DTFAC, FACMAX ) )
         ELSE
            DTFAC = MAX( FACMIN, MIN( DTFAC, FACONE ) )
         ENDIF 

         DT = MIN( DTMAX, MAX( DTMIN, DTFAC * DT ) )

         LPASS = .TRUE.

      ELSE

#ifdef rbstats
         IF( NFAILS .EQ. 0 .AND. TNOW .EQ. 0.0 ) N_BAD_STARTS = N_BAD_STARTS + 1
         NFAILS = NFAILS + 1
#endif


c        DTFAC = MAX( FACMIN, MIN( DTFAC, FACONE ) )

         DT = FACMIN * DT  

         LPASS = .FALSE.

         IF( DT .LT. DTMIN ) THEN

            WRITE( LOGDEV, 92100 ) JDATE, JTIME

            STOP

         ENDIF

      ENDIF

      IF( TNOW .LT. CHEMSTEP ) GO TO 100

      RETURN

  999 CONTINUE



92100 FORMAT( '      Convergence failure ', 
     &       '  JDATE = ', I7, '  JTIME = ' , I6 )


      END

