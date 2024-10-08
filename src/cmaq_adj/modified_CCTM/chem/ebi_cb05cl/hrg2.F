
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_cb05cl/hrg2.F,v 1.1.1.1 2010/06/14 16:03:01 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE HRG2( DTC )

C**********************************************************************
C
C  FUNCTION: To solve for the concentration of HO, HO2, HONO, and
C            HNO4 alebraically.
C
C  PRECONDITIONS: For the CB05CL family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, Nov. 6, 2006
C
C**********************************************************************
      USE HRDATA

      IMPLICIT NONE


C..INCLUDES: None


C..ARGUMENTS:
      REAL DTC                      ! Time step

C..PARAMETERS: None


C..EXTERNAL FUNCTIONS: NONE


C..SAVED LOCAL VARIABLES:
      CHARACTER( 16 ), SAVE :: PNAME = 'HRG2'             ! Program name


C..SCRATCH LOCAL VARIABLES:
      REAL    O1D_S                 ! sum of O1D loss frequencies
      REAL    OH_S                  ! stoich coeff for OH from O1D+H2O or H2
      REAL    HO2_S                 ! stoich coeff for HO2 from O1D+H2 rxn
      REAL    R4_19                 ! production of OH from HONO
      REAL    R19_4                 ! production of HONO from OH
      REAL    R4_5                  ! production of OH from HO2
      REAL    R5_4                  ! production of HO2 from OH
      REAL    R5_21                 ! production of HO2 from HNO4
      REAL    R21_5                 ! production of HNO4 from HO2
      REAL    P4, P5, P19, P21      ! species production form other paths
      REAL    L4, L5, L19, L21      ! species loss terms

      REAL    A, B, C               ! coeffs of quadratic eq. for HO2
      REAL    Q                     ! intermediate term

      REAL    T1, T2, T3            ! intermediate terms

      REAL    L21_INV               ! reciprocal of HNO4 loss term

C**********************************************************************


c..stoichiometric coefficient for production of HOX from O3 via O1D
      O1D_S =            RKI(  10 )                   ! O1D=O
     &        +          RKI(  11 )                   ! O1D=2*OH
     &        +          RKI(  38 )                   ! O1D=OH+HO2

      OH_S  =    2.000 * RKI(  11 )                   ! O1D=2*OH
     &        +          RKI(  38 )                   ! O1D=OH+HO2

      OH_S  = OH_S / O1D_S

      HO2_S =            RKI(  38 )                   ! O1D=HO2+OH

      HO2_S  = HO2_S / O1D_S

c..Production of HO from HO2 (r4,5 terms )
      R4_5  =            RKI(  13 ) * YC( O3      )   ! HO2+O3=OH
     &        +          RKI(  30 ) * YC( NO      )   ! HO2+NO=OH+NO2
     &        +          RKI(  44 ) * YC( O       )   ! HO2+O=OH

      R4_5  = R4_5 * DTC

c..Production of HO from HONO (r4,19 terms )
      R4_19 =            RKI(  25 )                   ! HONO=OH+NO

      R4_19 = R4_19 * DTC

c..Remaining HO production
      P4 =   OH_S  * RXRAT(   9 )      ! O1D=2*OH
     &    +  2.000 * RXRAT(  36 )      ! H2O2=2*OH
     &    +          RXRAT(  45 )      ! H2O2+O=OH+HO2
     &    +  0.390 * RXRAT(  51 )      ! PNA=0.39*OH+0.61*HO2+...
     &    +          RXRAT(  52 )      ! HNO3=OH+NO2
     &    +          RXRAT(  64 )      ! ROOH=OH+HO2+0.5*ALD2+0.5*ALDX
     &    +          RXRAT(  71 )      ! MEPX=OH+HO2+FORM
     &    +          RXRAT(  76 )      ! FORM+O=OH+HO2+CO
     &    +          RXRAT(  83 )      ! ALD2+O=OH+C2O3
     &    +          RXRAT(  96 )      ! PACD=OH+MEO2
     &    +          RXRAT(  98 )      ! ALDX+O=OH+CXO3
     &    +  0.100 * RXRAT( 116 )      ! O+OLE=0.1*OH+0.3*HO2+...
     &    +  0.100 * RXRAT( 118 )      ! O3+OLE=0.1*OH+0.44*HO2+...
     &    +  0.300 * RXRAT( 120 )      ! O+ETH=0.3*OH+1.7*HO2+CO+...
     &    +  0.130 * RXRAT( 122 )      ! O3+ETH=0.13*OH+0.13*HO2+...
     &    +  0.500 * RXRAT( 126 )      ! IOLE+O3=0.5*OH+0.5*HO2+...
     &    +  0.080 * RXRAT( 137 )      ! OPEN+O3=0.08*OH+0.76*HO2+...
     &    +  0.266 * RXRAT( 143 )      ! O3+ISOP=0.266*OH+0.066*HO2+...
     &    +  0.268 * RXRAT( 146 )      ! O3+ISPD=0.268*OH+0.154*HO2+...
     &    +  0.570 * RXRAT( 151 )      ! TERP+O3=0.57*OH+0.07*HO2+...
     &    +          RXRAT( 158 )      ! HOCL=OH+CL

      P4 = YC0( OH ) + P4 * DTC

c..Production of HO2 from OH ( r5,4 terms )
      R5_4  =            RKI(  12 ) * YC( O3      )   ! OH+O3=HO2
     &        +          RKI(  37 ) * YC( H2O2    )   ! OH+H2O2=HO2
     &        +          RKI(  39 )                   ! OH=HO2
     &        +          RKI(  40 ) * YC( O       )   ! OH+O=HO2
     &        +          RKI(  47 ) * YC( NO3     )   ! OH+NO3=HO2+NO2
     &        +          RKI(  61 ) * YC( NTR     )   ! OH+NTR=HO2+HNO3+...
     &        +          RKI(  65 ) * YC( CO      )   ! OH+CO=HO2
     &        +  0.300 * RKI(  70 ) * YC( MEPX    )   ! OH+MEPX=0.3*HO2+...
     &        +          RKI(  72 ) * YC( MEOH    )   ! OH+MEOH=HO2+FORM
     &        +          RKI(  73 ) * YC( FORM    )   ! OH+FORM=HO2+CO
     &        +          RKI(  82 ) * YC( FACD    )   ! OH+FACD=HO2
     &        +  0.110 * RKI( 112 ) * YC( PAR     )   ! OH+PAR=0.11*HO2+...
     &        +  0.950 * RKI( 117 ) * YC( OLE     )   ! OH+OLE=0.95*HO2+...
     &        +          RKI( 121 ) * YC( ETH     )   ! OH+ETH=HO2+...
     &        +          RKI( 125 ) * YC( IOLE    )   ! OH+IOLE=HO2+...
     &        +  0.440 * RKI( 128 ) * YC( TOL     )   ! OH+TOL=0.44*HO2+...
     &        +  0.600 * RKI( 131 ) * YC( CRES    )   ! OH+CRES=0.6*HO2+...
     &        +          RKI( 136 ) * YC( OPEN    )   ! OH+OPEN=2*HO2+...
     &        +  0.700 * RKI( 138 ) * YC( XYL     )   ! OH+XYL=0.7*HO2+...
     &        +  0.912 * RKI( 142 ) * YC( ISOP    )   ! OH+ISOP=...
     &        +  0.503 * RKI( 145 ) * YC( ISPD    )   ! OH+ISPD=...
     &        +  0.750 * RKI( 150 ) * YC( TERP    )   ! OH+TERP=0.75*HO2+...
     &        +          RKI( 153 ) * YC( SO2     )   ! OH+SO2=HO2+SULF+...
     &        +          RKI( 154 ) * YC( ETOH    )   ! OH+ETOH=HO2+...
     &        +          RKI( 155 ) * YC( ETHA    )   ! OH+ETHA=HO2+...

      R5_4  = R5_4 * DTC

c..Production of HO2 from HNO4 (r5,21 term )
      R5_21 =            RKI(  32 )                   ! PNA=HO2+NO2
     &        +  0.610 * RKI(  51 )                   ! PNA=0.61*HO2+...

      R5_21 = R5_21 * DTC

c..Remaining HO2 production terms
      P5 =   HO2_S  *RXRAT(   9 )      ! O1D=HO2+OH
     &    +          RXRAT(  45 )      ! H2O2+O=HO2+OH
     &    +          RXRAT(  62 )      ! NTR=HO2+NO2+0.33*FORM+...
     &    +          RXRAT(  64 )      ! ROOH=HO2+OH+0.5*ALD2+0.5*ALDX
     &    +          RXRAT(  67 )      ! MEO2+NO=HO2+FORM+NO2
     &    +  0.740 * RXRAT(  69 )      ! MEO2+MEO2=0.74*HO2+...
     &    +          RXRAT(  71 )      ! MEPX=HO2+OH+FORM
     &    +  2.000 * RXRAT(  74 )      ! FORM=2*HO2+CO
     &    +          RXRAT(  76 )      ! FORM+O=HO2+OH+CO
     &    +          RXRAT(  77 )      ! FORM+NO3=HO2+HNO3+CO
     &    +          RXRAT(  79 )      ! HCO3=HO2+FORM
     &    +          RXRAT(  80 )      ! HCO3+NO=HO2+NO2+FACD
     &    +          RXRAT(  86 )      ! ALD2=HO2+CO+MEO2
     &    +  0.900 * RXRAT(  92 )      ! C2O3+MEO2=0.9*HO2+0.9*MEO2+...
     &    +          RXRAT( 101 )      ! ALDX=HO2+CO+MEO2
     &    +          RXRAT( 102 )      ! CXO3+NO=HO2+NO2+ALD2+XO2
     &    +          RXRAT( 108 )      ! CXO3+MEO2=HO2+0.9*XO2+...
     &    +  2.000 * RXRAT( 110 )      ! CXO3+CXO3=2*HO2+2*XO2+2*ALD2
     &    +          RXRAT( 111 )      ! CXO3+C2O3=HO2+XO2+MEO2+ALD2
     &    +  0.940 * RXRAT( 113 )      ! ROR=0.94*HO2+0.6*ALD2+...
     &    +          RXRAT( 114 )      ! ROR=HO2
     &    +  0.300 * RXRAT( 116 )      ! O+OLE=0.3*HO2+0.1*OH+...
     &    +  0.440 * RXRAT( 118 )      ! O3+OLE=0.44*HO2+0.1*OH+...
     &    +  1.700 * RXRAT( 120 )      ! O+ETH=1.7*HO2+0.3*OH+CO+...
     &    +  0.130 * RXRAT( 122 )      ! O3+ETH=0.13*HO2+0.13*OH+...
     &    +  0.100 * RXRAT( 124 )      ! IOLE+O=0.1*HO2+0.66*ALDX+...
     &    +  0.500 * RXRAT( 126 )      ! IOLE+O3=0.5*HO2+0.5*OH+...
     &    +          RXRAT( 127 )      ! IOLE+NO3=HO2+0.64*ALDX+...
     &    +  0.900 * RXRAT( 129 )      ! TO2+NO=0.9*HO2+0.9*NO2+...
     &    +          RXRAT( 130 )      ! TO2=HO2+CRES
     &    +          RXRAT( 135 )      ! OPEN=HO2+C2O3+CO
     &    +          RXRAT( 136 )      ! OPEN+OH=2*HO2+2*CO+XO2+...
     &    +  0.760 * RXRAT( 137 )      ! OPEN+O3=0.76*HO2+0.08*OH+...
     &    +          RXRAT( 140 )      ! MGLY=HO2+C2O3+CO
     &    +  0.250 * RXRAT( 141 )      ! O+ISOP=0.25*HO2+0.5*FORM+...
     &    +  0.066 * RXRAT( 143 )      ! O3+ISOP=0.066*HO2+0.266*OH+...
     &    +  0.800 * RXRAT( 144 )      ! NO3+ISOP=0.8*HO2+0.8*NTR+...
     &    +  0.154 * RXRAT( 146 )      ! O3+ISPD=0.154*HO2+0.268*OH+...
     &    +  0.925 * RXRAT( 147 )      ! NO3+ISPD=0.925*HO2+...
     &    +  1.033 * RXRAT( 148 )      ! ISPD=1.033*HO2+0.067*ALD2+...
     &    +  0.070 * RXRAT( 151 )      ! TERP+O3=0.07*HO2+0.57*OH+...
     &    +  0.280 * RXRAT( 152 )      ! TERP+NO3=0.28*HO2+0.47*NO2+...
     &    +  0.800 * RXRAT( 156 )      ! NO2+ISOP=0.8*HO2+0.8*NTR+...
     &    +          RXRAT( 164 )      ! FMCL=HO2+CO+CL
     &    +  0.110 * RXRAT( 166 )      ! CL+PAR=0.11*HO2+0.87*XO2+...
     &    +          RXRAT( 167 )      ! CL+ETHA=HO2+0.991*ALD2+...
     &    +          RXRAT( 168 )      ! CL+ETH=HO2+2*XO2+FMCL+FORM
     &    +          RXRAT( 169 )      ! CL+OLE=HO2+0.33*ALD2+...
     &    +          RXRAT( 170 )      ! CL+IOLE=HO2+0.7*FMCL+...
     &    +          RXRAT( 171 )      ! CL+ISOP=HO2+XO2+0.15*HCL+...
     &    +          RXRAT( 172 )      ! CL+FORM=HO2+HCL+CO
     &    +          RXRAT( 175 )      ! CL+MEOH=HO2+HCL+FORM
     &    +          RXRAT( 176 )      ! CL+ETOH=HO2+HCL+ALD2

      P5 = YC0( HO2 ) + P5 * DTC

c..Production of HONO from OH (r19,4 terms )

      R19_4 =   RKI(  24 ) * YC( NO      ) * DTC      ! OH+NO=HONO

c..Remaining HONO production terms
      P19 =   2.000 * RXRAT(  23 )      ! NO+NO2=2*HONO

      P19 = YC0( HONO ) + P19 * DTC

c..Production of PNA from HO2 (r21,5 term )

      R21_5 =   RKI(  31 ) * YC( NO2     ) * DTC      ! HO2+NO2=HNO4

c..Remaining PNA production terms
      P21   =   YC0( PNA )

c..OH loss terms not in R5_4 & R19_4
      L4    =            RKI(  26 ) * YC( HONO    )   ! OH+HONO=NO2
     &        +          RKI(  28 ) * YC( NO2     )   ! OH+NO2=HNO3
     &        +          RKI(  29 ) * YC( HNO3    )   ! OH+HNO3=NO3
     &        +          RKI(  33 ) * YC( PNA     )   ! OH+PNA=NO2
     &        +  2.000 * RKI(  41 ) * YC( OH      )   ! OH+OH=O
     &        +  2.000 * RKI(  42 ) * YC( OH      )   ! OH+OH=H2O2
     &        +          RKI(  43 ) * YC( HO2     )   ! OH+HO2=
     &        +          RKI(  63 ) * YC( ROOH    )   ! OH+ROOH=XO2+...
     &        +          RKI(  66 )                   ! OH=MEO2
     &        +  0.700 * RKI(  70 ) * YC( MEPX    )   ! OH+MEPX=0.3*HO2+...
     &        +          RKI(  84 ) * YC( ALD2    )   ! OH+ALD2=C2O3
     &        +          RKI(  95 ) * YC( PACD    )   ! OH+PACD=C2O3
     &        +          RKI(  97 ) * YC( AACD    )   ! OH+AACD=MEO2
     &        +          RKI(  99 ) * YC( ALDX    )   ! OH+ALDX=CXO3
     &        +          RKI( 106 ) * YC( PANX    )   ! OH+PANX=ALD2+NO2
     &        +  0.890 * RKI( 112 ) * YC( PAR     )   ! OH+PAR=0.11*HO2+...
     &        +  0.050 * RKI( 117 ) * YC( OLE     )   ! OH+OLE=0.95*HO2+...
     &        +  0.560 * RKI( 128 ) * YC( TOL     )   ! OH+TOL=0.44*HO2+...
     &        +  0.400 * RKI( 131 ) * YC( CRES    )   ! OH+CRES=0.6*HO2+...
     &        +  0.300 * RKI( 138 ) * YC( XYL     )   ! OH+XYL=0.7*HO2+...
     &        +          RKI( 139 ) * YC( MGLY    )   ! OH+MGLY=XO2+C2O3
     &        +  0.088 * RKI( 142 ) * YC( ISOP    )   ! OH+ISOP=...
     &        +  0.497 * RKI( 145 ) * YC( ISPD    )   ! OH+ISPD=...
     &        +  0.250 * RKI( 150 ) * YC( TERP    )   ! OH+TERP=0.75*HO2+...
     &        +          RKI( 163 ) * YC( FMCL    )   ! OH+FMCL=CL+CO
     &        +          RKI( 177 ) * YC( HCL     )   ! OH+HCL=CL

      L4    = 1.0 + L4 * DTC + R5_4 + R19_4

c..HO2 loss terms not included in R4_5 & R21_5 (except for HO2+HO2 )
      L5    =            RKI(  43 ) * YC( OH      )   ! HO2+OH=
     &        +          RKI(  48 ) * YC( NO3     )   ! HO2+NO3=HNO3
     &        +          RKI(  56 ) * YC( XO2     )   ! HO2+XO2=ROOH
     &        +          RKI(  57 ) * YC( XO2N    )   ! HO2+XO2N=ROOH
     &        +          RKI(  68 ) * YC( MEO2    )   ! HO2+MEO2=MEPX
     &        +          RKI(  78 ) * YC( FORM    )   ! HO2+FORM=HCO3
     &        +          RKI(  81 ) * YC( HCO3    )   ! HO2+HCO3=MEPX
     &        +          RKI(  91 ) * YC( C2O3    )   ! HO2+C2O3=...
     &        +          RKI( 107 ) * YC( CXO3    )   ! HO2+CXO3=...
     &        +          RKI( 134 ) * YC( CRO     )   ! HO2+CRO=CRES
     &        +          RKI( 162 ) * YC( CLO     )   ! HO2+CLO=HOCL

      L5    = 1.0 + L5 * DTC + R4_5 + R21_5

c..HONO loss terms not included in R4_19
      L19   =            RKI(  26 ) * YC( OH      )   ! HONO+OH=NO2
     &        +  2.000 * RKI(  27 ) * YC( HONO    )   ! HONO+HONO=NO+NO2

      L19   = 1.0 + L19 * DTC + R4_19

c..PNA loss terms not inluded in R5_21
      L21   =            RKI(  33 ) * YC( OH      )   ! PNA+OH=NO2
     &        +  0.390 * RKI(  51 )                   ! PNA=0.61*HO2+...

      L21   = 1.0 + L21 * DTC + R5_21

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Solution section
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c..compute terms used to calculate a,b & c
      L21_INV = 1.0 / L21
      T1 = 1.0 / ( L4 * L19 - R4_19 * R19_4 )
      T2 = R5_4 * T1
      T3 = R5_21 * L21_INV

c..solve quadratic equation for HO2
      A = 2.0 * ( RKI(  34 ) + RKI(  35 ) ) * DTC

      B = L5 - T3 * R21_5 - T2 * R4_5 * L19

      C = P5 + T3 * P21 + T2 * ( P4 * L19 + P19 * R4_19 )

      Q = -0.5 * ( B + SIGN( 1.0, B ) * SQRT( B * B + 4.0 * A * C ) )

      YCP( HO2 ) = MAX( Q / A , -C / Q  )

c..compute remaining species concentrations
      YCP( OH ) = ( ( P4 + R4_5 * YCP( HO2 ) ) * L19 + R4_19 * P19 ) * T1

      YCP( PNA ) = ( P21 + R21_5 * YCP( HO2 ) ) * L21_INV

      YCP( HONO ) = ( P19 + R19_4 * YCP( OH ) ) / L19

      RETURN

      END


