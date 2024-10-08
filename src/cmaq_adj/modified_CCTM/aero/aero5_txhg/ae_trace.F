
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_txhg/ae_trace.F,v 1.1.1.1 2010/06/14 16:03:00 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      MODULE AE_TRACE
C**********************************************************************
C
C  FUNCTION:  Define arrays that identify aerosol species not in the 
C             original version of the mechanims. These sepcies are 
C             assumed to contribute insignifcant mass to each aerosol
C             model.
C
C  REVISION HISTORY: 01/01/06 : B.Hutzell - Initial version
C
C**********************************************************************

      USE CGRID_SPCS   ! CGRID species number and offsets

      IMPLICIT NONE


C..declare and define variables used by maps and data

      CHARACTER(16), PARAMETER :: BLANK = ' '  ! default value for characters

      INTEGER              :: N_TRACE_I        ! number of Aitken trace species 
      INTEGER              :: N_TRACE_J        ! number of Accum. trace species 
      INTEGER              :: N_TRACE_K        ! number of Coarse trace species 
      INTEGER              :: N_TRACE          ! total number of aerosol trace species

      INTEGER, ALLOCATABLE :: AE_TRACE_MAP( : )  ! indicies map of Aitken trace species
      
      INTEGER, ALLOCATABLE :: AE_TRACE_I( : )  ! indicies map of Aitken trace species
      INTEGER, ALLOCATABLE :: AE_TRACE_J( : )  ! indicies map of Accum. trace species
      INTEGER, ALLOCATABLE :: AE_TRACE_C( : )  ! indicies map of Coarse trace species

      REAL,          ALLOCATABLE :: AE_TRACE_CONC( : )     ! mass concentrations of trace species, ug/m^3
      CHARACTER(16), ALLOCATABLE :: AE_TRACE_NAME( : )     ! names of trace species

      INTEGER              :: N_PAIR_IJ	                   ! number of substances in both I and J modes
      INTEGER, ALLOCATABLE :: TRACE_PAIR_IJ( :, : )        ! indicies for a substance in I and J modes 
      INTEGER              :: N_PAIR_IK                    ! number of substances in both I and K modes
      INTEGER, ALLOCATABLE :: TRACE_PAIR_IK( :, : )        ! indicies for a substance in I and K modes 
      INTEGER              :: N_PAIR_JK                    ! number of substances in both J and K modes
      INTEGER, ALLOCATABLE :: TRACE_PAIR_JK( :, : )        ! indicies for a substance in J and K modes 

      INTEGER                    :: N_TRIO                 ! number of substances in all modes
      INTEGER,       ALLOCATABLE :: AE_TRACE_TRIO( :, : )  ! indicies for a substance in each mode
      CHARACTER(16), ALLOCATABLE :: TRACE_TRIO( : )        ! name of trio substance

C**********************************************************************


      END MODULE AE_TRACE

