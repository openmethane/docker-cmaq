
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_cb05txhg_ae5/init_degrade.F,v 1.1.1.1 2010/06/14 16:03:01 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE INIT_DEGRADE( CBLK, TCELL, DCELL, PHOTO_CELL,
     &                         JDATE, JTIME, LOGDEV )
C**********************************************************************
C
C  FUNCTION:  Initialize arrays used by degrade routines then load
C             CBLK concentration needed in degrade routines.
C
C  CALLED BY: HRDRIVER
C
C  REVISION HISTORY:  07/29/05 : B.Hutzell - Initial version
C
C**********************************************************************

      USE DEGRADE_SETUP_TOX

      IMPLICIT NONE

      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block

C.....ARGUMENTS:

      REAL        :: CBLK(*)       !  species concentration in cell
      REAL        :: TCELL         !  cell temperature  [ k ]
      REAL        :: DCELL         !  cell air density  [ kg/m^3 ]
      REAL        :: PHOTO_CELL(*) !  Photolysis table for cell [1/s]

      INTEGER     :: JDATE         ! current model date , coded YYYYDDD
      INTEGER     :: JTIME         ! current model time , coded HHMMSS
      INTEGER     :: LOGDEV        ! unit number for the log

C.....LOCAL VARIABLES:

      CHARACTER( 144 )        :: XMSG                   ! Message text
      CHARACTER( 16  ), SAVE  :: PNAME = 'INIT_DEGRADE' ! Routine name

      REAL(8), SAVE ::  CONV_M2N      ! factor to convert ppm times mass density in [kg/m^3]
                                      ! into number density in [molecules/cm^3]
      REAL(8)       ::  CONV_FACT     ! conversion factor from ppm to molecules/cm^3
      REAL(8)       ::  INV_TEMP      ! reciprocal of temperature, [K^-1]

      INTEGER       :: I, J, K        ! loop counters

      LOGICAL, SAVE ::  FIRSTCALL
      DATA              FIRSTCALL / .TRUE. /

C**********************************************************************

      IF ( FIRSTCALL ) THEN  ! initialize constants and allocate arrays

         CALL DEGRADE_DATA()

         IF(N_REACT < 1)THEN
           XMSG = 'UNABLE TO SET DEGRADE MODULE BASED ON MECHNAME '
     &         // 'WHICH EQUALS ' // MECHNAME
           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         ENDIF

         CALL DEGRADE_MAP( JDATE, JTIME, LOGDEV )

         CONV_M2N = 1.0D-9 * AVO / MWAIR       ! AVO and MWAIR defined in CONST.EXT

         ALLOCATE( PREV_CONC( NSPCSD ) )
         ALLOCATE( CURR_CONC( NSPCSD ) )
         ALLOCATE( DELT_CONC( NSPCSD ) )

         FIRSTCALL = .FALSE.

      ENDIF

C..initialize concentrations and their changes

      PREV_CONC( 1:NSPCSD ) = CBLK( 1:NSPCSD )
      CURR_CONC( 1:NSPCSD ) = CBLK( 1:NSPCSD )
      DELT_CONC( 1:NSPCSD ) = 0.0D0

      NUMB_DENS = 1.0E-3 * DCELL / MWAIR
      TEMP = TCELL

      RATE_CONST = 0.0D0
      CONV_FACT = CONV_M2N * DCELL
      INV_TEMP = 1.0D0 / TEMP

      LOOP_REACT: DO I = 1, N_REACT ! calculated rate constants

         LOOP_UNIRATE: DO J = 1, N_UNI_LOSS

            RATE_CONST( I, J ) = UNIRATE( I, J )
     &                         * TEMP**UNI_TEXP( I, J )
     &                         * EXP( -UNI_ACT( I, J )
     &                         * INV_TEMP )

         ENDDO LOOP_UNIRATE

         LOOP_BIRATE: DO J = 1, N_BI_LOSS

            RATE_CONST( I, J+UNI_STOP ) = BIRATE( I, J )
     &                                  * TEMP**BI_TEXP( I, J )
     &                                  * EXP( -BI_ACT( I, J )
     &                                  * INV_TEMP )
     &                                  * CONV_FACT

         ENDDO LOOP_BIRATE

         LOOP_TRIRATE: DO J = 1, N_TRI_LOSS

            RATE_CONST( I, J+BI_STOP ) = TRIRATE( I, J )
     &                                 * TEMP**TRI_TEXP( I, J )
     &                                 * EXP( -TRI_ACT( I, J )
     &                                 * INV_TEMP )
     &                                 * CONV_FACT
     &                                 * CONV_FACT

         ENDDO LOOP_TRIRATE

         LOOP_PHOTORATE: DO J = 1, N_PHOTO_LOSS

            K = PHOTO_MAP( I, J )
            IF ( K < 1 ) CYCLE
            RATE_CONST( I, J+TRI_STOP ) = A_PHOTO( I, J )
     &                                  * PHOTO_CELL( K )

         ENDDO LOOP_PHOTORATE

      ENDDO LOOP_REACT

      RETURN

      END SUBROUTINE INIT_DEGRADE
