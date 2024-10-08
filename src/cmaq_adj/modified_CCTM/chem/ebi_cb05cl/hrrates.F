
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_cb05cl/hrrates.F,v 1.1.1.1 2010/06/14 16:03:01 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE HRRATES

C***********************************************************************
C
C  FUNCTION: To calculate rates of reactions
C
C  PRECONDITIONS: For the CB05CL family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, Nov. 6, 2006
C
C***********************************************************************
      USE HRDATA

      IMPLICIT NONE

C.....INCLUDES: NONE


C.....ARGUMENTS: NONE


C.....PARAMETERS: NONE


C.....EXTERNAL FUNCTIONS: NONE


C.....SAVED VARIABLES: NONE
      CHARACTER( 16 ), SAVE  ::  PNAME = 'HRRATES'   ! Program name


C.....LOCAL VARIABLES: None

C***********************************************************************
C..NOTE: RXRAT commented out are not referenced by any routine
!     RXRAT(   1 ) = RKI(   1 ) * YC( NO2 )
!     RXRAT(   2 ) = RKI(   2 ) * YC( O )
!     RXRAT(   3 ) = RKI(   3 ) * YC( O3 ) * YC( NO )
!     RXRAT(   4 ) = RKI(   4 ) * YC( O ) * YC( NO2 )
      RXRAT(   5 ) = RKI(   5 ) * YC( O ) * YC( NO2 )
!     RXRAT(   6 ) = RKI(   6 ) * YC( O ) * YC( NO )
      RXRAT(   7 ) = RKI(   7 ) * YC( NO2 ) * YC( O3 )
      RXRAT(   8 ) = RKI(   8 ) * YC( O3 )
      RXRAT(   9 ) = RKI(   9 ) * YC( O3 )
!     RXRAT(  10 ) = RKI(  10 ) * YC( O1D )
!     RXRAT(  11 ) = RKI(  11 ) * YC( O1D )
!     RXRAT(  12 ) = RKI(  12 ) * YC( O3 ) * YC( OH )
!     RXRAT(  13 ) = RKI(  13 ) * YC( O3 ) * YC( HO2 )
      RXRAT(  14 ) = RKI(  14 ) * YC( NO3 )
      RXRAT(  15 ) = RKI(  15 ) * YC( NO3 )
      RXRAT(  16 ) = RKI(  16 ) * YC( NO3 ) * YC( NO )
      RXRAT(  17 ) = RKI(  17 ) * YC( NO3 ) * YC( NO2 )
!     RXRAT(  18 ) = RKI(  18 ) * YC( NO3 ) * YC( NO2 )
      RXRAT(  19 ) = RKI(  19 ) * YC( N2O5 )
      RXRAT(  20 ) = RKI(  20 ) * YC( N2O5 )
      RXRAT(  21 ) = RKI(  21 ) * YC( N2O5 )
!     RXRAT(  22 ) = RKI(  22 ) * YC( NO ) * YC( NO )
      RXRAT(  23 ) = RKI(  23 ) * YC( NO ) * YC( NO2 )
!     RXRAT(  24 ) = RKI(  24 ) * YC( NO ) * YC( OH )
      RXRAT(  25 ) = RKI(  25 ) * YC( HONO )
      RXRAT(  26 ) = RKI(  26 ) * YC( OH ) * YC( HONO )
      RXRAT(  27 ) = RKI(  27 ) * YC( HONO ) * YC( HONO )
      RXRAT(  28 ) = RKI(  28 ) * YC( NO2 ) * YC( OH )
      RXRAT(  29 ) = RKI(  29 ) * YC( OH ) * YC( HNO3 )
!     RXRAT(  30 ) = RKI(  30 ) * YC( HO2 ) * YC( NO )
!     RXRAT(  31 ) = RKI(  31 ) * YC( HO2 ) * YC( NO2 )
      RXRAT(  32 ) = RKI(  32 ) * YC( PNA )
      RXRAT(  33 ) = RKI(  33 ) * YC( OH ) * YC( PNA )
      RXRAT(  34 ) = RKI(  34 ) * YC( HO2 ) * YC( HO2 )
      RXRAT(  35 ) = RKI(  35 ) * YC( HO2 ) * YC( HO2 )
      RXRAT(  36 ) = RKI(  36 ) * YC( H2O2 )
      RXRAT(  37 ) = RKI(  37 ) * YC( OH ) * YC( H2O2 )
!     RXRAT(  38 ) = RKI(  38 ) * YC( O1D )
!     RXRAT(  39 ) = RKI(  39 ) * YC( OH )
!     RXRAT(  40 ) = RKI(  40 ) * YC( OH ) * YC( O )
      RXRAT(  41 ) = RKI(  41 ) * YC( OH ) * YC( OH )
      RXRAT(  42 ) = RKI(  42 ) * YC( OH ) * YC( OH )
!     RXRAT(  43 ) = RKI(  43 ) * YC( OH ) * YC( HO2 )
!     RXRAT(  44 ) = RKI(  44 ) * YC( HO2 ) * YC( O )
      RXRAT(  45 ) = RKI(  45 ) * YC( H2O2 ) * YC( O )
      RXRAT(  46 ) = RKI(  46 ) * YC( NO3 ) * YC( O )
      RXRAT(  47 ) = RKI(  47 ) * YC( NO3 ) * YC( OH )
      RXRAT(  48 ) = RKI(  48 ) * YC( NO3 ) * YC( HO2 )
      RXRAT(  49 ) = RKI(  49 ) * YC( NO3 ) * YC( O3 )
      RXRAT(  50 ) = RKI(  50 ) * YC( NO3 ) * YC( NO3 )
      RXRAT(  51 ) = RKI(  51 ) * YC( PNA )
      RXRAT(  52 ) = RKI(  52 ) * YC( HNO3 )
      RXRAT(  53 ) = RKI(  53 ) * YC( N2O5 )
      RXRAT(  54 ) = RKI(  54 ) * YC( XO2 ) * YC( NO )
      RXRAT(  55 ) = RKI(  55 ) * YC( XO2N ) * YC( NO )
      RXRAT(  56 ) = RKI(  56 ) * YC( XO2 ) * YC( HO2 )
      RXRAT(  57 ) = RKI(  57 ) * YC( XO2N ) * YC( HO2 )
      RXRAT(  58 ) = RKI(  58 ) * YC( XO2 ) * YC( XO2 )
      RXRAT(  59 ) = RKI(  59 ) * YC( XO2N ) * YC( XO2N )
      RXRAT(  60 ) = RKI(  60 ) * YC( XO2 ) * YC( XO2N )
      RXRAT(  61 ) = RKI(  61 ) * YC( NTR ) * YC( OH )
      RXRAT(  62 ) = RKI(  62 ) * YC( NTR )
      RXRAT(  63 ) = RKI(  63 ) * YC( ROOH ) * YC( OH )
      RXRAT(  64 ) = RKI(  64 ) * YC( ROOH )
      RXRAT(  65 ) = RKI(  65 ) * YC( OH ) * YC( CO )
      RXRAT(  66 ) = RKI(  66 ) * YC( OH )
      RXRAT(  67 ) = RKI(  67 ) * YC( MEO2 ) * YC( NO )
      RXRAT(  68 ) = RKI(  68 ) * YC( MEO2 ) * YC( HO2 )
      RXRAT(  69 ) = RKI(  69 ) * YC( MEO2 ) * YC( MEO2 )
      RXRAT(  70 ) = RKI(  70 ) * YC( MEPX ) * YC( OH )
      RXRAT(  71 ) = RKI(  71 ) * YC( MEPX )
      RXRAT(  72 ) = RKI(  72 ) * YC( MEOH ) * YC( OH )
      RXRAT(  73 ) = RKI(  73 ) * YC( FORM ) * YC( OH )
      RXRAT(  74 ) = RKI(  74 ) * YC( FORM )
      RXRAT(  75 ) = RKI(  75 ) * YC( FORM )
      RXRAT(  76 ) = RKI(  76 ) * YC( FORM ) * YC( O )
      RXRAT(  77 ) = RKI(  77 ) * YC( FORM ) * YC( NO3 )
      RXRAT(  78 ) = RKI(  78 ) * YC( FORM ) * YC( HO2 )
      RXRAT(  79 ) = RKI(  79 ) * YC( HCO3 )
      RXRAT(  80 ) = RKI(  80 ) * YC( HCO3 ) * YC( NO )
      RXRAT(  81 ) = RKI(  81 ) * YC( HCO3 ) * YC( HO2 )
      RXRAT(  82 ) = RKI(  82 ) * YC( FACD ) * YC( OH )
      RXRAT(  83 ) = RKI(  83 ) * YC( ALD2 ) * YC( O )
      RXRAT(  84 ) = RKI(  84 ) * YC( ALD2 ) * YC( OH )
      RXRAT(  85 ) = RKI(  85 ) * YC( ALD2 ) * YC( NO3 )
      RXRAT(  86 ) = RKI(  86 ) * YC( ALD2 )
      RXRAT(  87 ) = RKI(  87 ) * YC( C2O3 ) * YC( NO )
!     RXRAT(  88 ) = RKI(  88 ) * YC( C2O3 ) * YC( NO2 )
      RXRAT(  89 ) = RKI(  89 ) * YC( PAN )
      RXRAT(  90 ) = RKI(  90 ) * YC( PAN )
      RXRAT(  91 ) = RKI(  91 ) * YC( C2O3 ) * YC( HO2 )
      RXRAT(  92 ) = RKI(  92 ) * YC( C2O3 ) * YC( MEO2 )
      RXRAT(  93 ) = RKI(  93 ) * YC( C2O3 ) * YC( XO2 )
      RXRAT(  94 ) = RKI(  94 ) * YC( C2O3 ) * YC( C2O3 )
      RXRAT(  95 ) = RKI(  95 ) * YC( PACD ) * YC( OH )
      RXRAT(  96 ) = RKI(  96 ) * YC( PACD )
      RXRAT(  97 ) = RKI(  97 ) * YC( AACD ) * YC( OH )
      RXRAT(  98 ) = RKI(  98 ) * YC( ALDX ) * YC( O )
      RXRAT(  99 ) = RKI(  99 ) * YC( ALDX ) * YC( OH )
      RXRAT( 100 ) = RKI( 100 ) * YC( ALDX ) * YC( NO3 )
      RXRAT( 101 ) = RKI( 101 ) * YC( ALDX )
      RXRAT( 102 ) = RKI( 102 ) * YC( CXO3 ) * YC( NO )
      RXRAT( 103 ) = RKI( 103 ) * YC( CXO3 ) * YC( NO2 )
      RXRAT( 104 ) = RKI( 104 ) * YC( PANX )
      RXRAT( 105 ) = RKI( 105 ) * YC( PANX )
      RXRAT( 106 ) = RKI( 106 ) * YC( PANX ) * YC( OH )
      RXRAT( 107 ) = RKI( 107 ) * YC( CXO3 ) * YC( HO2 )
      RXRAT( 108 ) = RKI( 108 ) * YC( CXO3 ) * YC( MEO2 )
      RXRAT( 109 ) = RKI( 109 ) * YC( CXO3 ) * YC( XO2 )
      RXRAT( 110 ) = RKI( 110 ) * YC( CXO3 ) * YC( CXO3 )
      RXRAT( 111 ) = RKI( 111 ) * YC( CXO3 ) * YC( C2O3 )
      RXRAT( 112 ) = RKI( 112 ) * YC( PAR ) * YC( OH )
      RXRAT( 113 ) = RKI( 113 ) * YC( ROR )
      RXRAT( 114 ) = RKI( 114 ) * YC( ROR )
      RXRAT( 115 ) = RKI( 115 ) * YC( ROR ) * YC( NO2 )
      RXRAT( 116 ) = RKI( 116 ) * YC( O ) * YC( OLE )
      RXRAT( 117 ) = RKI( 117 ) * YC( OH ) * YC( OLE )
      RXRAT( 118 ) = RKI( 118 ) * YC( O3 ) * YC( OLE )
      RXRAT( 119 ) = RKI( 119 ) * YC( NO3 ) * YC( OLE )
      RXRAT( 120 ) = RKI( 120 ) * YC( O ) * YC( ETH )
      RXRAT( 121 ) = RKI( 121 ) * YC( OH ) * YC( ETH )
      RXRAT( 122 ) = RKI( 122 ) * YC( O3 ) * YC( ETH )
      RXRAT( 123 ) = RKI( 123 ) * YC( NO3 ) * YC( ETH )
      RXRAT( 124 ) = RKI( 124 ) * YC( IOLE ) * YC( O )
      RXRAT( 125 ) = RKI( 125 ) * YC( IOLE ) * YC( OH )
      RXRAT( 126 ) = RKI( 126 ) * YC( IOLE ) * YC( O3 )
      RXRAT( 127 ) = RKI( 127 ) * YC( IOLE ) * YC( NO3 )
      RXRAT( 128 ) = RKI( 128 ) * YC( TOL ) * YC( OH )
      RXRAT( 129 ) = RKI( 129 ) * YC( TO2 ) * YC( NO )
      RXRAT( 130 ) = RKI( 130 ) * YC( TO2 )
      RXRAT( 131 ) = RKI( 131 ) * YC( OH ) * YC( CRES )
      RXRAT( 132 ) = RKI( 132 ) * YC( CRES ) * YC( NO3 )
      RXRAT( 133 ) = RKI( 133 ) * YC( CRO ) * YC( NO2 )
      RXRAT( 134 ) = RKI( 134 ) * YC( CRO ) * YC( HO2 )
      RXRAT( 135 ) = RKI( 135 ) * YC( OPEN )
      RXRAT( 136 ) = RKI( 136 ) * YC( OPEN ) * YC( OH )
      RXRAT( 137 ) = RKI( 137 ) * YC( OPEN ) * YC( O3 )
      RXRAT( 138 ) = RKI( 138 ) * YC( OH ) * YC( XYL )
      RXRAT( 139 ) = RKI( 139 ) * YC( OH ) * YC( MGLY )
      RXRAT( 140 ) = RKI( 140 ) * YC( MGLY )
      RXRAT( 141 ) = RKI( 141 ) * YC( O ) * YC( ISOP )
      RXRAT( 142 ) = RKI( 142 ) * YC( OH ) * YC( ISOP )
      RXRAT( 143 ) = RKI( 143 ) * YC( O3 ) * YC( ISOP )
      RXRAT( 144 ) = RKI( 144 ) * YC( NO3 ) * YC( ISOP )
      RXRAT( 145 ) = RKI( 145 ) * YC( OH ) * YC( ISPD )
      RXRAT( 146 ) = RKI( 146 ) * YC( O3 ) * YC( ISPD )
      RXRAT( 147 ) = RKI( 147 ) * YC( NO3 ) * YC( ISPD )
      RXRAT( 148 ) = RKI( 148 ) * YC( ISPD )
      RXRAT( 149 ) = RKI( 149 ) * YC( TERP ) * YC( O )
      RXRAT( 150 ) = RKI( 150 ) * YC( TERP ) * YC( OH )
      RXRAT( 151 ) = RKI( 151 ) * YC( TERP ) * YC( O3 )
      RXRAT( 152 ) = RKI( 152 ) * YC( TERP ) * YC( NO3 )
      RXRAT( 153 ) = RKI( 153 ) * YC( SO2 ) * YC( OH )
      RXRAT( 154 ) = RKI( 154 ) * YC( OH ) * YC( ETOH )
      RXRAT( 155 ) = RKI( 155 ) * YC( OH ) * YC( ETHA )
      RXRAT( 156 ) = RKI( 156 ) * YC( NO2 ) * YC( ISOP )
      RXRAT( 157 ) = RKI( 157 ) * YC( CL2 )
      RXRAT( 158 ) = RKI( 158 ) * YC( HOCL )
      RXRAT( 159 ) = RKI( 159 ) * YC( CL ) * YC( O3 )
      RXRAT( 160 ) = RKI( 160 ) * YC( CLO ) * YC( CLO )
      RXRAT( 161 ) = RKI( 161 ) * YC( CLO ) * YC( NO )
      RXRAT( 162 ) = RKI( 162 ) * YC( CLO ) * YC( HO2 )
      RXRAT( 163 ) = RKI( 163 ) * YC( OH ) * YC( FMCL )
      RXRAT( 164 ) = RKI( 164 ) * YC( FMCL )
      RXRAT( 165 ) = RKI( 165 ) * YC( CL )
      RXRAT( 166 ) = RKI( 166 ) * YC( CL ) * YC( PAR )
      RXRAT( 167 ) = RKI( 167 ) * YC( CL ) * YC( ETHA )
      RXRAT( 168 ) = RKI( 168 ) * YC( CL ) * YC( ETH )
      RXRAT( 169 ) = RKI( 169 ) * YC( CL ) * YC( OLE )
      RXRAT( 170 ) = RKI( 170 ) * YC( CL ) * YC( IOLE )
      RXRAT( 171 ) = RKI( 171 ) * YC( CL ) * YC( ISOP )
      RXRAT( 172 ) = RKI( 172 ) * YC( CL ) * YC( FORM )
      RXRAT( 173 ) = RKI( 173 ) * YC( CL ) * YC( ALD2 )
      RXRAT( 174 ) = RKI( 174 ) * YC( CL ) * YC( ALDX )
      RXRAT( 175 ) = RKI( 175 ) * YC( CL ) * YC( MEOH )
      RXRAT( 176 ) = RKI( 176 ) * YC( CL ) * YC( ETOH )
      RXRAT( 177 ) = RKI( 177 ) * YC( HCL ) * YC( OH )


      RETURN

      END
