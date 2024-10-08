
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

C RCS file, release, date & time of last delta, author, state, [locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_txhg/aeroproc_trace.F,v 1.1.1.1 2010/06/14 16:03:00 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


       SUBROUTINE AEROPROC_TRACE( COAG_DECAY_I, MERGE_I, MERGE_J )

**********************************************************************
C
C  FUNCTION:  ADJUST aerosol trace species in the aitken and accumulation modes
C             based on coagulation rates and and mode merging in aero_proc.
C
C  REVISION HISTORY: 01/31/06 : B.Hutzell - Initial version
C
C**********************************************************************

      USE AE_TRACE

      IMPLICIT NONE

C ....inputs

      REAL*8  :: COAG_DECAY_I
      REAL    :: MERGE_I
      REAL    :: MERGE_J

C ....local
   
      REAL, PARAMETER  :: CONMIN = 1.0E-30 ! floor on concentrations
      REAL             :: TOTAL_MASS       ! mass in aitken and accumulation modes
      INTEGER          :: I, J, K          ! loop counters
      INTEGER          :: I_AITK     
      INTEGER          :: I_ACCU

      IF(N_TRACE .LT. 1)RETURN

      COAG_MERG : DO I = 1, N_PAIR_IJ

         I_AITK = TRACE_PAIR_IJ( I, 1 )
         I_ACCU = TRACE_PAIR_IJ( I, 2 )

         TOTAL_MASS = AE_TRACE_CONC( I_AITK ) + AE_TRACE_CONC( I_ACCU )

         AE_TRACE_CONC( I_AITK ) = AE_TRACE_CONC( I_AITK ) * COAG_DECAY_I
         AE_TRACE_CONC( I_AITK ) = MAX( CONMIN,  AE_TRACE_CONC( I_AITK ) )
         AE_TRACE_CONC( I_ACCU ) = MAX( CONMIN,  TOTAL_MASS - AE_TRACE_CONC( I_AITK ) )

          IF ( (MERGE_I + MERGE_J) .LE. 0.0 ) CYCLE

         AE_TRACE_CONC( I_ACCU ) = AE_TRACE_CONC( I_ACCU )
     &                           + AE_TRACE_CONC( I_AITK ) * MERGE_J

         AE_TRACE_CONC( I_AITK ) = AE_TRACE_CONC( I_AITK ) * MERGE_I

         AE_TRACE_CONC( I_AITK ) = MAX( CONMIN,  AE_TRACE_CONC( I_AITK ) )
         AE_TRACE_CONC( I_ACCU ) = MAX( CONMIN,  AE_TRACE_CONC( I_ACCU ) )

      END DO COAG_MERG

      RETURN

      END SUBROUTINE AEROPROC_TRACE
          
