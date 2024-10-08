
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_saprc99tx3_ae5/hrdata_mod.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      MODULE  HRDATA

C*************************************************************************
C
C  FUNCTION:  Mechanism & solver data for EBI solver
C
C  PRECONDITIONS: For the SAPRC99TX3 family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, May 30, 2008
C
C*************************************************************************

c..EBI solver fixed parameters
      INTEGER, PARAMETER ::  NEBITER = 1000  ! No. of iterations for EBI

      REAL, PARAMETER    ::  DELTAT = 2.5D+00     ! EBI time step

c..Mechanism specific variables
      INTEGER   N_SPEC       ! No. of species in mechanism
      INTEGER   N_RXNS       ! No. of reactions in mechanism
      INTEGER   N_EBISP      ! No. of species solved by EBI
      INTEGER   NING1        ! No. of species in group 1
      INTEGER   NING2        ! No. of species in group 2

c..Control flags
      LOGICAL   L_AE_VRSN    ! Flag for aerosol version of mech
      LOGICAL   L_AQ_VRSN    ! Flag for aqueous chemistry version of mech

      LOGICAL, SAVE :: CALL_DEG = .TRUE.           ! WTH: SWITCH for calling DEGRADE routine

c..Miscellaneous variables
      INTEGER   LOGDEV       ! Unit number of output log
      INTEGER   N_EBI_IT     ! No. of iterations for EBI
      INTEGER   N_EBI_STEPS  ! No. of time steps for EBI
      INTEGER   N_INR_STEPS  ! No. of inner time steps for EBI
      INTEGER   MODEL_SPECIES ! WTH: number of species in CGRID

      REAL    EBI_TMSTEP   ! Time step for EBI loops (min)

c...Allocatable arrays
      INTEGER, ALLOCATABLE :: EBISP( : )         ! Index of EBI species

      REAL*8, ALLOCATABLE ::  RKI(   : )     ! Rate constants
      REAL, ALLOCATABLE ::  RXRAT( : )     ! Reaction rates
      REAL, ALLOCATABLE ::  RTOL(  : )     ! Species tolerances
      REAL, ALLOCATABLE ::  YC(    : )     ! Species concentrations
      REAL, ALLOCATABLE ::  YCCELL( :)     ! WTH: ALL Species in cgrid
      REAL, ALLOCATABLE ::  YC0(   : )     ! Species concentrations
      REAL, ALLOCATABLE ::  YCP(   : )     ! Species concentrations
      REAL, ALLOCATABLE ::  PROD(  : )     ! Prod of species
      REAL, ALLOCATABLE ::  LOSS(  : )     ! Loss of species


c..Species indices
      INTEGER   NO2
      INTEGER   NO
      INTEGER   O3P
      INTEGER   O3
      INTEGER   NO3
      INTEGER   N2O5
      INTEGER   HNO3
      INTEGER   O1D2
      INTEGER   HO
      INTEGER   HONO
      INTEGER   HO2
      INTEGER   CO
      INTEGER   HNO4
      INTEGER   HO2H
      INTEGER   SO2
      INTEGER   SULF
      INTEGER   SULRXN
      INTEGER   C_O2
      INTEGER   HCHO
      INTEGER   COOH
      INTEGER   MEOH
      INTEGER   RO2_R
      INTEGER   ROOH
      INTEGER   R2O2
      INTEGER   RO2_N
      INTEGER   RNO3
      INTEGER   MEK
      INTEGER   PROD2
      INTEGER   CCO_O2
      INTEGER   PAN
      INTEGER   CCO_OOH
      INTEGER   CCO_OH
      INTEGER   RCO_O2
      INTEGER   PAN2
      INTEGER   CCHO
      INTEGER   RCO_OOH
      INTEGER   RCO_OH
      INTEGER   BZCO_O2
      INTEGER   PBZN
      INTEGER   BZ_O
      INTEGER   MA_RCO3
      INTEGER   MA_PAN
      INTEGER   TBU_O
      INTEGER   ACET
      INTEGER   NPHE
      INTEGER   PHEN
      INTEGER   BZNO2_O
      INTEGER   HOCOO
      INTEGER   HCOOH
      INTEGER   RCHO
      INTEGER   GLY
      INTEGER   MGLY
      INTEGER   BACL
      INTEGER   CRES
      INTEGER   BALD
      INTEGER   METHACRO
      INTEGER   MVK
      INTEGER   ISOPROD
      INTEGER   CCHO_SUR
      INTEGER   DCB1
      INTEGER   DCB2
      INTEGER   DCB3
      INTEGER   ETHENE
      INTEGER   ISOPRENE
      INTEGER   ISOPRXN
      INTEGER   TRP1
      INTEGER   TRPRXN
      INTEGER   ALK1
      INTEGER   ALK2
      INTEGER   ALK3
      INTEGER   ALK4
      INTEGER   ALK5
      INTEGER   ALK5RXN
      INTEGER   ARO1
      INTEGER   ARO1RO2
      INTEGER   TOLNRXN
      INTEGER   TOLHRXN
      INTEGER   ARO2
      INTEGER   ARO2RO2
      INTEGER   XYLNRXN
      INTEGER   XYLHRXN
      INTEGER   OLE1
      INTEGER   OLE2
      INTEGER   SESQ
      INTEGER   SESQRXN
      INTEGER   BENZENE
      INTEGER   BENZRO2
      INTEGER   BNZNRXN
      INTEGER   BNZHRXN
      INTEGER   ACROLEIN
      INTEGER   BUTADIENE13
      INTEGER   FORM_PRIMARY
      INTEGER   ALD2_PRIMARY
      INTEGER   ACROLEIN_PRIMARY

      END MODULE HRDATA
