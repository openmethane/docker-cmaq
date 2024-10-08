
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ros3/degrade.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE DEGRADE( CBLK, DT, JDATE, JTIME, BLKID )
C**********************************************************************
C
C Function: Calculate changes in gas species based on a exponential decay.
C           The decay rate sums losses from processes in DEGRADE_DATA.
C
C CALLED BY: RBSOLVER or GRSMVGEAR
C
C WARNING: THIS ROUTINE ASSUMES SIMPLE AND LINEAR TRANSFORMATIONS FROM
C          ATMOSPHERIC CHEMISTRY.
C
C Species being degraded are governed by the equation,
C     dx/dt = -b*x, where b is the sum of N loss rates
C
C IT DOES NOT SOLVE A SYSTEM OF ODE's AS IN SMVGEAR, ROS3, and EBI SOLVERS.
C
C  REVISION HISTORY:  07/29/05 : B.Hutzell - Initial version
C
C**********************************************************************

      USE DEGRADE_SETUP_TOX

      IMPLICIT NONE

C.....ARGUMENTS:

      REAL(8)      CBLK( BLKSIZE, N_GC_SPC )   ! array holding species concentrations
      REAL(8)      TSTEP                       ! time step for integrations [sec]
      INTEGER      JDATE                       ! current model date , coded YYYYDDD
      INTEGER      JTIME                       ! current model time , coded HHMMSS
      INTEGER      BLKID                       ! ID number for the BLK

C.....PARAMETERS:

      CHARACTER(16), PARAMETER :: PNAME = ' DEGRADE    '  ! name of routine calling I/OAPI

      INTEGER, PARAMETER :: LOCAL_DT = 3     ! minimum time step, mili-seconds

      REAL(8), PARAMETER :: CONMIN = 1.0D-30 ! concentration lower limit

C.....LOCAL VARIABLES:

      CHARACTER(16)  ::  VNAME                   ! variable name
      CHARACTER(120) ::  XMSG

      LOGICAL, SAVE  ::  FIRSTCALL  = .TRUE.

      INTEGER        :: TIME_SECONDS                ! TIME, sec
      INTEGER        :: I_RXT, I_RAD, J_RAD, I_PROD ! indices
      INTEGER        :: I, J, K, L, I_CELL          ! loop counters
      INTEGER, SAVE  :: I_SIZE                      ! scratch

      REAL(8)        ::  DT                                ! degradation time step, sec
      REAL(8)        ::  TRANS    ( BLKSIZE )              ! molecules/cm^3 transferred to products
      REAL(8)        ::  NET_RATE ( BLKSIZE )              ! net rate of transfer   [sec^-1]
      REAL(8)        ::  NET_LIFE ( BLKSIZE )              ! lifetime based on net transfer rate  [sec]
      REAL(8)        ::  LOSS_RATE( BLKSIZE, N_PROCESSES ) ! individual loss rates  [sec^-1]

C***********************************************************************

      IF ( FIRSTCALL ) THEN  ! initialize maps
         I_SIZE = SIZE( CURR_CONC, 2 )
         FIRSTCALL = .FALSE.
      ENDIF

C..Initialize concentrations changes

      DELT_CONC = 0.0D0

C..Quality Control on time step

      TSTEP = DT
      BLOCK_A : IF ( TSTEP < 0.0D0 ) THEN
         XMSG = ' Time step has negative value. '
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
      ENDIF BLOCK_A

C..update current concentrations

      DO J = 1, N_GC_SPC
         DO I_CELL = 1, NUM_CELLS
            CURR_CONC( I_CELL, J ) = CBLK( I_CELL, J )
         ENDDO
      ENDDO

C..Loop over each reactant

      LOOP_REACT: DO I = 1, N_REACT

         I_RXT = RXTANT_MAP( I )

         LOSS_RATE = 0.0D0
         NET_RATE  = 0.0D0
         NET_LIFE  = 0.0D0

         LOOP_UNIRATE: DO J = UNI_START, UNI_STOP
            DO I_CELL = 1, NUM_CELLS
               LOSS_RATE( I_CELL, J ) = RATE_CONST( I_CELL, I, J )
            ENDDO
         ENDDO LOOP_UNIRATE

         L = 0

         LOOP_BIRATE: DO J = BI_START, BI_STOP
            L = L + 1
            I_RAD = RAD_MAP( I, L )
            IF ( I_RAD < 1 ) CYCLE   ! radical species is undefined
            IF ( I_RAD > 900 ) THEN
               DO I_CELL = 1, NUM_CELLS
                  LOSS_RATE( I_CELL, J ) = RATE_CONST( I_CELL, I, J )
     &                                  * NUMB_DENS( I_CELL )
               ENDDO
            ELSE
               DO I_CELL = 1, NUM_CELLS
                  LOSS_RATE( I_CELL, J ) = 0.5D0 * RATE_CONST( I_CELL, I, J )
     &                                   * ( PREV_CONC( I_CELL, I_RAD )
     &                                   +   CURR_CONC( I_CELL, I_RAD ) )
               ENDDO
            ENDIF
         ENDDO LOOP_BIRATE

         L = 0

         LOOP_TRIRATE: DO J = TRI_START, TRI_STOP
            L = L + 1
            I_RAD = RAD2_MAP( I, L, 1 )
            J_RAD = RAD2_MAP( I, L, 2 )
            IF ( I_RAD < 1 .OR. J_RAD < 1 ) CYCLE   ! radical species are undefined
            DO I_CELL = 1, NUM_CELLS
               LOSS_RATE( I_CELL, J ) =  0.5D0 * RATE_CONST( I_CELL, I, J )
     &                                * ( PREV_CONC( I_CELL, I_RAD )
     &                                *   PREV_CONC( I_CELL, J_RAD )
     &                                +   CURR_CONC( I_CELL, I_RAD )
     &                                *   CURR_CONC( I_CELL, J_RAD ) )
            ENDDO
         ENDDO LOOP_TRIRATE

         L = 0
         LOOP_PHOTORATE: DO J = PHOTO_START, PHOTO_STOP
            L = L + 1
            DO I_CELL = 1, NUM_CELLS
               LOSS_RATE( I_CELL, J ) = RATE_CONST( I_CELL, I, J )
            ENDDO
         ENDDO LOOP_PHOTORATE

         LOOP_RATE :  DO J = 1, N_PROCESSES
            DO I_CELL = 1, NUM_CELLS
               NET_RATE( I_CELL ) = NET_RATE( I_CELL )
     &                            + LOSS_RATE( I_CELL, J )
            ENDDO
         ENDDO LOOP_RATE

         LOOP_LIFE: DO I_CELL = 1, NUM_CELLS
C           IF ( NET_RATE( I_CELL ) * DT .LE. 1.0D-20 ) CYCLE
            NET_LIFE( I_CELL ) = 1.0D0 / NET_RATE( I_CELL )
            TRANS( I_CELL ) = CURR_CONC( I_CELL, I_RXT )
     &                      * ( 1.0D0
     &                      -  DBLE( EXP( -REAL( NET_RATE( I_CELL )
     &                      *  DT ) ) ) )
     &                      - CONMIN
            DELT_CONC( I_CELL, I_RXT ) = -MAX( TRANS( I_CELL ), 0.0D0 )
         ENDDO LOOP_LIFE

         LOOP_PROD: DO I_CELL = 1, NUM_CELLS
            IF ( DELT_CONC( I_CELL, I_RXT ) <  -CONMIN ) THEN
               DO J = 1, N_PROCESSES
                  I_PROD = PROD_MAP( I, J )
                  IF ( I_PROD < 1 ) CYCLE    !  No specified product
                  DELT_CONC( I_CELL, I_PROD ) = -LOSS_RATE( I_CELL, J )
     &                                        *  NET_LIFE ( I_CELL )
     &                                        *  RATE_YIELD( I, J )
     &                                        *  DELT_CONC( I_CELL, I_RXT )
               ENDDO
            ENDIF
         ENDDO LOOP_PROD

      ENDDO LOOP_REACT

C..update previous concentrations

      DO J = 1, NSPCSD
         DO I_CELL = 1, NUM_CELLS
            PREV_CONC( I_CELL, J ) = CURR_CONC( I_CELL, J )
         ENDDO
      ENDDO

C..update current concentrations

      DO J = 1, NSPCSD
         DO I_CELL = 1, NUM_CELLS
            CURR_CONC( I_CELL, J ) = PREV_CONC( I_CELL, J )
     &                             + DELT_CONC( I_CELL, J )
         ENDDO
      ENDDO

      RETURN
      END
