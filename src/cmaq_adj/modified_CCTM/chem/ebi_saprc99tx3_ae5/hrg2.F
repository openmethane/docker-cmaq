
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_saprc99tx3_ae5/hrg2.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE HRG2( DTC )
C**********************************************************************
C
C  FUNCTION: To solve for the concentration of HO, HO2, HONO, and
C            HNO4 alebraically.
C
C  PRECONDITIONS: For the SAPRC99TX3 family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, May 30, 2008
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
      O1D_S =            RKI(  19 )                            ! O1D2=2*HO
     &        +          RKI(  20 )                            ! O1D2=O3P

      OH_S  =    2.000 * RKI(  19 )                            ! O1D2=2*HO

      OH_S  = OH_S / O1D_S


c..Production of HO from HO2 (r4,5 terms )
      R4_5  =            RKI(  31 ) * YC( NO               )   ! HO2+NO=HO+NO2
     &        +          RKI(  36 ) * YC( O3               )   ! HO2+O3=HO
     &        +  0.800 * RKI(  39 ) * YC( NO3              )   ! HO2+NO3=0.8*HO+...

      R4_5  = R4_5 * DTC

c..Production of HO from HONO (r4,19 terms )
      R4_19 =            RKI(  22 )                            ! HONO=HO+NO

      R4_19 = R4_19 * DTC

c..Remaining HO production
      P4 =   OH_S  * RXRAT(  18 )      ! O1D2=2*HO
     &    +          RXRAT(  28 )      ! HNO3=HO+NO2
     &    +  0.390 * RXRAT(  34 )      ! HNO4=0.39*HO+0.61*HO2+...
     &    +  2.000 * RXRAT(  41 )      ! HO2H=2*HO
     &    +          RXRAT( 142 )      ! COOH=HO+HO2+HCHO
     &    +          RXRAT( 144 )      ! ROOH=HO+HO2+RCHO
     &    +  0.208 * RXRAT( 162 )      ! METHACRO+O3=0.208*HO+...
     &    +  0.330 * RXRAT( 165 )      ! METHACRO=0.33*HO+0.34*HO2+...
     &    +  0.164 * RXRAT( 167 )      ! MVK+O3=0.164*HO+0.064*HO2+...
     &    +  0.285 * RXRAT( 171 )      ! ISOPROD+O3=0.285*HO+...
     &    +  0.500 * RXRAT( 179 )      ! DCB1+O3=0.5*HO+1.5*HO2+...
     &    +  0.120 * RXRAT( 186 )      ! ETHENE+O3=0.12*HO+0.12*HO2+...
     &    +  0.266 * RXRAT( 190 )      ! ISOPRENE+O3=0.266*HO+...
     &    +  0.567 * RXRAT( 194 )      ! TRP1+O3=0.567*HO+0.033*HO2+...
     &    +  0.155 * RXRAT( 209 )      ! OLE1+O3=0.155*HO+0.056*HO2+...
     &    +  0.378 * RXRAT( 213 )      ! OLE2+O3=0.378*HO+0.003*HO2+...
     &    +  0.310 * RXRAT( 226 )      ! ACROLEIN+O3=0.31*HO+...
     &    +  0.172 * RXRAT( 229 )      ! ACROLEIN=0.172*HO+1.01*HO2+...
     &    +  0.060 * RXRAT( 231 )      ! BUTADIENE13+O3=0.06*HO+...

      P4 = YC0( HO ) + P4 * DTC

c..Production of HO2 from OH ( r5,4 terms )
      R5_4  =            RKI(  26 ) * YC( NO3              )   ! HO+NO3=HO2+NO2
     &        +          RKI(  29 ) * YC( CO               )   ! HO+CO=HO2
     &        +          RKI(  30 ) * YC( O3               )   ! HO+O3=HO2
     &        +          RKI(  42 ) * YC( HO2H             )   ! HO+HO2H=HO2
     &        +          RKI(  44 ) * YC( SO2              )   ! HO+SO2=HO2+SULF+...
     &        +          RKI(  45 )                            ! HO=HO2
     &        +          RKI( 125 ) * YC( HCHO             )   ! HO+HCHO=HO2+CO
     &        +          RKI( 140 ) * YC( MEOH             )   ! HO+MEOH=HO2+HCHO
     &        +  0.630 * RKI( 147 ) * YC( GLY              )   ! HO+GLY=0.63*HO2+...
     &        +  0.379 * RKI( 174 ) * YC( PROD2            )   ! HO+PROD2=...
     &        +  0.113 * RKI( 176 ) * YC( RNO3             )   ! HO+RNO3=...
     &        +  0.121 * RKI( 198 ) * YC( ALK2             )   ! HO+ALK2=...
     &        +  0.224 * RKI( 202 ) * YC( ARO1             )   ! HO+ARO1=...
     &        +  0.187 * RKI( 205 ) * YC( ARO2             )   ! HO+ARO2=...
     &        +          RKI( 216 ) * YC( HCOOH            )   ! HO+HCOOH=HO2
     &        +  0.236 * RKI( 222 ) * YC( BENZENE          )   ! HO+BENZENE=...

      R5_4  = R5_4 * DTC

c..Production of HO2 from HNO4 (r5,21 term )
      R5_21 =            RKI(  33 )                            ! HNO4=HO2+NO2
     &        +  0.610 * RKI(  34 )                            ! HNO4=0.61*HO2+...

      R5_21 = R5_21 * DTC

c..Remaining HO2 production terms
      P5 =           RXRAT(  23 )      ! HONO=HO2+NO2
     &    +          RXRAT(  46 )      ! C_O2+NO=HO2+HCHO+NO2
     &    +          RXRAT(  48 )      ! C_O2+NO3=HO2+HCHO+NO2
     &    +  2.000 * RXRAT(  50 )      ! C_O2+C_O2=2*HO2+2*HCHO
     &    +          RXRAT(  51 )      ! RO2_R+NO=HO2+NO2
     &    +          RXRAT(  53 )      ! RO2_R+NO3=HO2+NO2
     &    +          RXRAT(  54 )      ! RO2_R+C_O2=HO2+0.75*HCHO+...
     &    +          RXRAT(  55 )      ! RO2_R+RO2_R=HO2
     &    +          RXRAT(  64 )      ! RO2_N+C_O2=HO2+0.25*MEOH+...
     &    +          RXRAT(  65 )      ! RO2_N+NO3=HO2+NO2+MEK
     &    +          RXRAT(  66 )      ! RO2_N+RO2_R=HO2+0.5*MEK+...
     &    +          RXRAT(  68 )      ! RO2_N+RO2_N=HO2+MEK+PROD2
     &    +  2.000 * RXRAT( 123 )      ! HCHO=2*HO2+CO
     &    +          RXRAT( 127 )      ! HOCOO=HO2+HCHO
     &    +          RXRAT( 128 )      ! HOCOO+NO=HO2+NO2+HCOOH
     &    +          RXRAT( 129 )      ! HCHO+NO3=HO2+HNO3+CO
     &    +          RXRAT( 131 )      ! CCHO=HO2+CO+C_O2
     &    +          RXRAT( 134 )      ! RCHO=HO2+RO2_R+CO+CCHO
     &    +          RXRAT( 142 )      ! COOH=HO2+HO+HCHO
     &    +          RXRAT( 144 )      ! ROOH=HO2+HO+RCHO
     &    +  2.000 * RXRAT( 145 )      ! GLY=2*HO2+2*CO
     &    +  0.630 * RXRAT( 148 )      ! GLY+NO3=0.63*HO2+HNO3+...
     &    +          RXRAT( 149 )      ! MGLY=HO2+CO+CCO_O2
     &    +  0.008 * RXRAT( 162 )      ! METHACRO+O3=0.008*HO2+...
     &    +  0.340 * RXRAT( 165 )      ! METHACRO=0.34*HO2+0.33*HO+...
     &    +  0.064 * RXRAT( 167 )      ! MVK+O3=0.064*HO2+0.164*HO+...
     &    +  0.400 * RXRAT( 171 )      ! ISOPROD+O3=0.4*HO2+...
     &    +  1.233 * RXRAT( 173 )      ! ISOPROD=1.233*HO2+...
     &    +  0.341 * RXRAT( 177 )      ! RNO3=0.341*HO2+NO2+...
     &    +  1.500 * RXRAT( 179 )      ! DCB1+O3=1.5*HO2+0.5*HO+...
     &    +  0.500 * RXRAT( 181 )      ! DCB2=0.5*HO2+0.5*CCO_O2+...
     &    +  0.500 * RXRAT( 183 )      ! DCB3=0.5*HO2+0.5*CCO_O2+...
     &    +  0.120 * RXRAT( 186 )      ! ETHENE+O3=0.12*HO2+0.12*HO+...
     &    +  0.500 * RXRAT( 188 )      ! ETHENE+O3P=0.5*HO2+...
     &    +  0.033 * RXRAT( 194 )      ! TRP1+O3=0.033*HO2+0.567*HO+...
     &    +  0.056 * RXRAT( 209 )      ! OLE1+O3=0.056*HO2+0.155*HO+...
     &    +  0.003 * RXRAT( 213 )      ! OLE2+O3=0.003*HO2+0.378*HO+...
     &    +  0.013 * RXRAT( 215 )      ! OLE2+O3P=0.013*HO2+...
     &    +  0.810 * RXRAT( 226 )      ! ACROLEIN+O3=0.81*HO2+...
     &    +  1.010 * RXRAT( 229 )      ! ACROLEIN=1.01*HO2+0.172*HO+...
     &    +  0.060 * RXRAT( 231 )      ! BUTADIENE13+O3=0.06*HO2+...
     &    +  0.250 * RXRAT( 233 )      ! BUTADIENE13+O3P=0.25*HO2+...

      P5 = YC0( HO2 ) + P5 * DTC

c..Production of HONO from OH (r19,4 terms )

      R19_4 =   RKI(  21 ) * YC( NO               ) * DTC      ! OH+NO=HONO

c..Remaining HONO production terms
      P19   = YC0( HONO )

c..Production of HNO4 from HO2 (r21,5 term )

      R21_5 =   RKI(  32 ) * YC( NO2              ) * DTC      ! HO2+NO2=HNO4

c..Remaining HNO4 production terms
      P21   =   YC0( HNO4 )

c..HO loss terms not in R5_4 & R19_4
      L4    =            RKI(  24 ) * YC( HONO             )   ! HO+HONO=NO2
     &        +          RKI(  25 ) * YC( NO2              )   ! HO+NO2=HNO3
     &        +          RKI(  27 ) * YC( HNO3             )   ! HO+HNO3=NO3
     &        +          RKI(  35 ) * YC( HNO4             )   ! HO+HNO4=NO2
     &        +          RKI(  43 ) * YC( HO2              )   ! HO+HO2=
     &        +          RKI( 130 ) * YC( CCHO             )   ! HO+CCHO=CCO_O2
     &        +          RKI( 133 ) * YC( RCHO             )   ! HO+RCHO=...
     &        +          RKI( 136 ) * YC( ACET             )   ! HO+ACET=HCHO+...
     &        +          RKI( 138 ) * YC( MEK              )   ! HO+MEK=...
     &        +  0.650 * RKI( 141 ) * YC( COOH             )   ! HO+COOH=...
     &        +  0.340 * RKI( 143 ) * YC( ROOH             )   ! HO+ROOH=RCHO+...
     &        +  0.370 * RKI( 147 ) * YC( GLY              )   ! HO+GLY=0.63*HO2+...
     &        +          RKI( 150 ) * YC( MGLY             )   ! HO+MGLY=CO+CCO_O2
     &        +          RKI( 153 ) * YC( PHEN             )   ! HO+PHEN=...
     &        +          RKI( 155 ) * YC( CRES             )   ! HO+CRES=...
     &        +          RKI( 158 ) * YC( BALD             )   ! HO+BALD=BZCO_O2
     &        +          RKI( 161 ) * YC( METHACRO         )   ! HO+METHACRO=...
     &        +          RKI( 166 ) * YC( MVK              )   ! HO+MVK=0.3*RO2_R+...
     &        +          RKI( 170 ) * YC( ISOPROD          )   ! HO+ISOPROD=...
     &        +  0.621 * RKI( 174 ) * YC( PROD2            )   ! HO+PROD2=...
     &        +  0.887 * RKI( 176 ) * YC( RNO3             )   ! HO+RNO3=...
     &        +          RKI( 178 ) * YC( DCB1             )   ! HO+DCB1=RCHO+...
     &        +          RKI( 180 ) * YC( DCB2             )   ! HO+DCB2=R2O2+...
     &        +          RKI( 182 ) * YC( DCB3             )   ! HO+DCB3=R2O2+...
     &        +          RKI( 184 )                            ! HO=C_O2
     &        +          RKI( 185 ) * YC( ETHENE           )   ! HO+ETHENE=RO2_R+...
     &        +          RKI( 189 ) * YC( ISOPRENE         )   ! HO+ISOPRENE=...
     &        +          RKI( 193 ) * YC( TRP1             )   ! HO+TRP1=...
     &        +          RKI( 197 ) * YC( ALK1             )   ! HO+ALK1=RO2_R+CCHO
     &        +  0.633 * RKI( 198 ) * YC( ALK2             )   ! HO+ALK2=...
     &        +          RKI( 199 ) * YC( ALK3             )   ! HO+ALK3=...
     &        +          RKI( 200 ) * YC( ALK4             )   ! HO+ALK4=...
     &        +          RKI( 201 ) * YC( ALK5             )   ! HO+ALK5=...
     &        +  0.776 * RKI( 202 ) * YC( ARO1             )   ! HO+ARO1=...
     &        +  0.813 * RKI( 205 ) * YC( ARO2             )   ! HO+ARO2=...
     &        +          RKI( 208 ) * YC( OLE1             )   ! HO+OLE1=...
     &        +          RKI( 212 ) * YC( OLE2             )   ! HO+OLE2=...
     &        +          RKI( 217 ) * YC( CCO_OH           )   ! HO+CCO_OH=...
     &        +          RKI( 218 ) * YC( RCO_OH           )   ! HO+RCO_OH=RO2_R+...
     &        +  0.764 * RKI( 222 ) * YC( BENZENE          )   ! HO+BENZENE=...
     &        +          RKI( 225 ) * YC( ACROLEIN         )   ! HO+ACROLEIN=...
     &        +          RKI( 230 ) * YC( BUTADIENE13      )   ! HO+BUTADIENE13=...

      L4    = 1.0 + L4 * DTC + R5_4 + R19_4

c..HO2 loss terms not included in R4_5 & R21_5 (except for HO2+HO2 )
      L5    =    0.200 * RKI(  39 ) * YC( NO3              )   ! HO2+NO3=0.8*HO+...
     &        +          RKI(  43 ) * YC( HO               )   ! HO2+HO=
     &        +          RKI(  47 ) * YC( C_O2             )   ! HO2+C_O2=COOH
     &        +          RKI(  52 ) * YC( RO2_R            )   ! HO2+RO2_R=ROOH
     &        +          RKI(  63 ) * YC( RO2_N            )   ! HO2+RO2_N=ROOH
     &        +          RKI(  72 ) * YC( CCO_O2           )   ! HO2+CCO_O2=...
     &        +          RKI(  82 ) * YC( RCO_O2           )   ! HO2+RCO_O2=...
     &        +          RKI(  93 ) * YC( BZCO_O2          )   ! HO2+BZCO_O2=...
     &        +          RKI( 105 ) * YC( MA_RCO3          )   ! HO2+MA_RCO3=...
     &        +          RKI( 118 ) * YC( BZ_O             )   ! HO2+BZ_O=PHEN
     &        +          RKI( 121 ) * YC( BZNO2_O          )   ! HO2+BZNO2_O=NPHE
     &        +          RKI( 126 ) * YC( HCHO             )   ! HO2+HCHO=HOCOO

      L5    = 1.0 + L5 * DTC + R4_5 + R21_5

c..HONO loss terms not included in R4_19
      L19   =            RKI(  23 )                            ! HONO=HO2+NO2
     &        +          RKI(  24 ) * YC( HO               )   ! HONO+HO=NO2

      L19   = 1.0 + L19 * DTC + R4_19

c..HNO4 loss terms not inluded in R5_21
      L21   =    0.390 * RKI(  34 )                            ! HNO4=0.61*HO2+...
     &        +          RKI(  35 ) * YC( HO               )   ! HNO4+HO=NO2

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
      A = 2.0 * ( RKI(  37 ) + RKI(  38 ) ) * DTC

      B = L5 - T3 * R21_5 - T2 * R4_5 * L19

      C = P5 + T3 * P21 + T2 * ( P4 * L19 + P19 * R4_19 )

      Q = -0.5 * ( B + SIGN( 1.0, B ) * SQRT( B * B + 4.0 * A * C ) )

      YCP( HO2 ) = MAX( Q / A , -C / Q  )

c..compute remaining species concentrations
      YCP( HO ) = ( ( P4 + R4_5 * YCP( HO2 ) ) * L19 + R4_19 * P19 ) * T1

      YCP( HNO4 ) = ( P21 + R21_5 * YCP( HO2 ) ) * L21_INV

      YCP( HONO ) = ( P19 + R19_4 * YCP( HO ) ) / L19

      RETURN

      END


