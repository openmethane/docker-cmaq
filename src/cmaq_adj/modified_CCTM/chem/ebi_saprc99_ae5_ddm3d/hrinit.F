
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_saprc99_ae5_ddm3d/hrinit.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE HRINIT
C***********************************************************************
C
C  FUNCTION: To initialize species tolerances, arrays, and indices
C
C  PRECONDITIONS: For the SAPRC99 family of mechanisms
C
C  RETURN VALUES: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, Sept. 10, 2008
C
C***********************************************************************

      USE HRDATA

      IMPLICIT NONE

C.....INCLUDES:
      INCLUDE SUBST_GC_SPC    ! Gas chem species names and MWs
      INCLUDE SUBST_GC_EMIS   ! Gas chem emissions name and mapping tables
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block

C.....ARGUMENTS: NONE

C.....PARAMETERS: NONE


C.....EXTERNAL FUNCTIONS:
      INTEGER  FINDEX         ! Finds location of a number in a list

C.....SAVED VARIABLES:
      CHARACTER( 16 ), SAVE  ::  PNAME = 'HRINIT'   ! Program name


C.....LOCAL VARIABLES:
      CHARACTER*132 MSG       ! Log message

      INTEGER IND             ! Species index
      INTEGER N               ! Loop index

C***********************************************************************

      N_SPEC = N_GC_SPC
      N_RXNS = NRXNS

      ALLOCATE( RKI( NRXNS ) )
      ALLOCATE( RXRAT( NRXNS ) )
      ALLOCATE( RTOL( N_SPEC) )
      ALLOCATE( YC(   N_SPEC) )
      ALLOCATE( YC0(  N_SPEC) )
      ALLOCATE( YCP(  N_SPEC) )
      ALLOCATE( PROD( N_SPEC) )
      ALLOCATE( LOSS( N_SPEC) )

      NING1 = 4
      NING2 = 4

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set species indices and pointers
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NO2       =   1
      NO        =   2
      O3P       =   3
      O3        =   4
      NO3       =   5
      N2O5      =   6
      HNO3      =   7
      O1D2      =   8
      HO        =   9
      HONO      =  10
      HO2       =  11
      CO        =  12
      HNO4      =  13
      HO2H      =  14
      SO2       =  15
      SULF      =  16
      SULRXN    =  17
      C_O2      =  18
      HCHO      =  19
      COOH      =  20
      MEOH      =  21
      RO2_R     =  22
      ROOH      =  23
      R2O2      =  24
      RO2_N     =  25
      RNO3      =  26
      MEK       =  27
      PROD2     =  28
      CCO_O2    =  29
      PAN       =  30
      CCO_OOH   =  31
      CCO_OH    =  32
      RCO_O2    =  33
      PAN2      =  34
      CCHO      =  35
      RCO_OOH   =  36
      RCO_OH    =  37
      BZCO_O2   =  38
      PBZN      =  39
      BZ_O      =  40
      MA_RCO3   =  41
      MA_PAN    =  42
      TBU_O     =  43
      ACET      =  44
      NPHE      =  45
      PHEN      =  46
      BZNO2_O   =  47
      HOCOO     =  48
      HCOOH     =  49
      RCHO      =  50
      GLY       =  51
      MGLY      =  52
      BACL      =  53
      CRES      =  54
      BALD      =  55
      METHACRO  =  56
      MVK       =  57
      ISOPROD   =  58
      DCB1      =  59
      DCB2      =  60
      DCB3      =  61
      ETHENE    =  62
      ISOPRENE  =  63
      ISOPRXN   =  64
      TRP1      =  65
      TRPRXN    =  66
      ALK1      =  67
      ALK2      =  68
      ALK3      =  69
      ALK4      =  70
      ALK5      =  71
      ALK5RXN   =  72
      ARO1      =  73
      ARO1RO2   =  74
      TOLNRXN   =  75
      TOLHRXN   =  76
      ARO2      =  77
      ARO2RO2   =  78
      XYLNRXN   =  79
      XYLHRXN   =  80
      BENZENE   =  81
      BENZRO2   =  82
      BNZNRXN   =  83
      BNZHRXN   =  84
      OLE1      =  85
      OLE2      =  86
      SESQ      =  87
      SESQRXN   =  88

      N_EBISP  =  75
      ALLOCATE( EBISP( N_EBISP ) ) 

      EBISP(   1 ) = HNO3
      EBISP(   2 ) = CO
      EBISP(   3 ) = HO2H
      EBISP(   4 ) = SO2
      EBISP(   5 ) = SULF
      EBISP(   6 ) = SULRXN
      EBISP(   7 ) = C_O2
      EBISP(   8 ) = HCHO
      EBISP(   9 ) = COOH
      EBISP(  10 ) = MEOH
      EBISP(  11 ) = RO2_R
      EBISP(  12 ) = ROOH
      EBISP(  13 ) = R2O2
      EBISP(  14 ) = RO2_N
      EBISP(  15 ) = RNO3
      EBISP(  16 ) = MEK
      EBISP(  17 ) = PROD2
      EBISP(  18 ) = CCO_OOH
      EBISP(  19 ) = CCO_OH
      EBISP(  20 ) = RCO_O2
      EBISP(  21 ) = PAN2
      EBISP(  22 ) = CCHO
      EBISP(  23 ) = RCO_OOH
      EBISP(  24 ) = RCO_OH
      EBISP(  25 ) = BZCO_O2
      EBISP(  26 ) = PBZN
      EBISP(  27 ) = BZ_O
      EBISP(  28 ) = MA_RCO3
      EBISP(  29 ) = MA_PAN
      EBISP(  30 ) = TBU_O
      EBISP(  31 ) = ACET
      EBISP(  32 ) = NPHE
      EBISP(  33 ) = PHEN
      EBISP(  34 ) = BZNO2_O
      EBISP(  35 ) = HOCOO
      EBISP(  36 ) = HCOOH
      EBISP(  37 ) = RCHO
      EBISP(  38 ) = GLY
      EBISP(  39 ) = MGLY
      EBISP(  40 ) = BACL
      EBISP(  41 ) = CRES
      EBISP(  42 ) = BALD
      EBISP(  43 ) = METHACRO
      EBISP(  44 ) = MVK
      EBISP(  45 ) = ISOPROD
      EBISP(  46 ) = DCB1
      EBISP(  47 ) = DCB2
      EBISP(  48 ) = DCB3
      EBISP(  49 ) = ETHENE
      EBISP(  50 ) = ISOPRENE
      EBISP(  51 ) = ISOPRXN
      EBISP(  52 ) = TRP1
      EBISP(  53 ) = TRPRXN
      EBISP(  54 ) = ALK1
      EBISP(  55 ) = ALK2
      EBISP(  56 ) = ALK3
      EBISP(  57 ) = ALK4
      EBISP(  58 ) = ALK5
      EBISP(  59 ) = ALK5RXN
      EBISP(  60 ) = ARO1
      EBISP(  61 ) = ARO1RO2
      EBISP(  62 ) = TOLNRXN
      EBISP(  63 ) = TOLHRXN
      EBISP(  64 ) = ARO2
      EBISP(  65 ) = ARO2RO2
      EBISP(  66 ) = XYLNRXN
      EBISP(  67 ) = XYLHRXN
      EBISP(  68 ) = BENZENE
      EBISP(  69 ) = BENZRO2
      EBISP(  70 ) = BNZNRXN
      EBISP(  71 ) = BNZHRXN
      EBISP(  72 ) = OLE1
      EBISP(  73 ) = OLE2
      EBISP(  74 ) = SESQ
      EBISP(  75 ) = SESQRXN


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set species tolerances
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RTOL( NO2      ) = 1.0E-03
      RTOL( NO       ) = 1.0E-03
      RTOL( O3P      ) = 1.0E+00
      RTOL( O3       ) = 1.0E-03
      RTOL( NO3      ) = 1.0E-03
      RTOL( N2O5     ) = 1.0E-03
      RTOL( HNO3     ) = 1.0E-03
      RTOL( O1D2     ) = 1.0E+00
      RTOL( HO       ) = 1.0E-03
      RTOL( HONO     ) = 1.0E-03
      RTOL( HO2      ) = 1.0E-03
      RTOL( CO       ) = 1.0E-03
      RTOL( HNO4     ) = 1.0E-03
      RTOL( HO2H     ) = 1.0E-03
      RTOL( SO2      ) = 1.0E-03
      RTOL( SULF     ) = 1.0E+00
      RTOL( SULRXN   ) = 1.0E+00
      RTOL( C_O2     ) = 1.0E-03
      RTOL( HCHO     ) = 1.0E-03
      RTOL( COOH     ) = 1.0E-03
      RTOL( MEOH     ) = 1.0E-03
      RTOL( RO2_R    ) = 1.0E-03
      RTOL( ROOH     ) = 1.0E-03
      RTOL( R2O2     ) = 1.0E-03
      RTOL( RO2_N    ) = 1.0E-03
      RTOL( RNO3     ) = 1.0E-03
      RTOL( MEK      ) = 1.0E-03
      RTOL( PROD2    ) = 1.0E-03
      RTOL( CCO_O2   ) = 1.0E-03
      RTOL( PAN      ) = 1.0E-03
      RTOL( CCO_OOH  ) = 1.0E-03
      RTOL( CCO_OH   ) = 1.0E-03
      RTOL( RCO_O2   ) = 1.0E-03
      RTOL( PAN2     ) = 1.0E-03
      RTOL( CCHO     ) = 1.0E-03
      RTOL( RCO_OOH  ) = 1.0E-03
      RTOL( RCO_OH   ) = 1.0E-03
      RTOL( BZCO_O2  ) = 1.0E-03
      RTOL( PBZN     ) = 1.0E-03
      RTOL( BZ_O     ) = 1.0E+00
      RTOL( MA_RCO3  ) = 1.0E-03
      RTOL( MA_PAN   ) = 1.0E-03
      RTOL( TBU_O    ) = 1.0E+00
      RTOL( ACET     ) = 1.0E-03
      RTOL( NPHE     ) = 1.0E-03
      RTOL( PHEN     ) = 1.0E-03
      RTOL( BZNO2_O  ) = 1.0E+00
      RTOL( HOCOO    ) = 1.0E+00
      RTOL( HCOOH    ) = 1.0E-03
      RTOL( RCHO     ) = 1.0E-03
      RTOL( GLY      ) = 1.0E-03
      RTOL( MGLY     ) = 1.0E-03
      RTOL( BACL     ) = 1.0E-03
      RTOL( CRES     ) = 1.0E-03
      RTOL( BALD     ) = 1.0E-03
      RTOL( METHACRO ) = 1.0E-03
      RTOL( MVK      ) = 1.0E-03
      RTOL( ISOPROD  ) = 1.0E-03
      RTOL( DCB1     ) = 1.0E-03
      RTOL( DCB2     ) = 1.0E-03
      RTOL( DCB3     ) = 1.0E-03
      RTOL( ETHENE   ) = 1.0E-03
      RTOL( ISOPRENE ) = 1.0E-03
      RTOL( ISOPRXN  ) = 1.0E+00
      RTOL( TRP1     ) = 1.0E-03
      RTOL( TRPRXN   ) = 1.0E+00
      RTOL( ALK1     ) = 1.0E-03
      RTOL( ALK2     ) = 1.0E-03
      RTOL( ALK3     ) = 1.0E-03
      RTOL( ALK4     ) = 1.0E-03
      RTOL( ALK5     ) = 1.0E-03
      RTOL( ALK5RXN  ) = 1.0E+00
      RTOL( ARO1     ) = 1.0E-03
      RTOL( ARO1RO2  ) = 1.0E-03
      RTOL( TOLNRXN  ) = 1.0E+00
      RTOL( TOLHRXN  ) = 1.0E+00
      RTOL( ARO2     ) = 1.0E-03
      RTOL( ARO2RO2  ) = 1.0E-03
      RTOL( XYLNRXN  ) = 1.0E+00
      RTOL( XYLHRXN  ) = 1.0E+00
      RTOL( BENZENE  ) = 1.0E-03
      RTOL( BENZRO2  ) = 1.0E-03
      RTOL( BNZNRXN  ) = 1.0E+00
      RTOL( BNZHRXN  ) = 1.0E+00
      RTOL( OLE1     ) = 1.0E-03
      RTOL( OLE2     ) = 1.0E-03
      RTOL( SESQ     ) = 1.0E-03
      RTOL( SESQRXN  ) = 1.0E+00


      RETURN

      END
