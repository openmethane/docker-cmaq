
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_saprc99_ae5/hrg1.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE HRG1( DTC )

C**********************************************************************
C
C  FUNCTION: To solve for the concentration of NO2, NO, O3, and O3P
C            algebraically.
C
C  PRECONDITIONS: For the SAPRC99 family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, Sept. 10, 2008
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
      CHARACTER( 16 ), SAVE  :: PNAME = 'HRG1'   ! Prgram Name


C..SCRATCH LOCAL VARIABLES:
      REAL O1D_S               ! sum of O1D loss frequencies
      REAL O3P_S               ! stoich coeff for O3P from O1D

      REAL R1_2                ! production term for NO from NO2
      REAL R2_1                ! production term for NO2 from NO
      REAL P1, P2, P3, P12     ! production terms for NO, NO2, O3, & O3P
      REAL L1, L2, L3, L12     ! loss terms for NO, NO2, O3, O3P
      REAL L1_INV, L2_INV,
     &     L3_INV, L12_INV     ! inverse of loss terms

      REAL T1, T2, T3, T4, T5  ! intermerdiate terms
      REAL F1, F2, F3          ! intermerdiate terms
      REAL A, B, C             ! coefficients for quadratic equation
      REAL Q, XX, S1, S2       ! intermerdiate terms

      REAL RK1, RK2, RK3       ! rate constants

      REAL PO3                 ! temp variable for O3

C**********************************************************************


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  O1D Section
c    1) sum of the rate constants for all O1D loss reactions
c    2) get fractional yield of O3P from O1D loss
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      O1D_S =            RKI(  19 )                    ! O1D2=2*HO
     &        +          RKI(  20 )                    ! O1D2=O3P

      O3P_S =            RKI(  20 )                    ! O1D2=O3P

      O3P_S  = O3P_S / O1D_S


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  NO Section
c    R1_2 = production of NO from NO2 ( rates of form k[NO2][x] )
c           except NO2+NO3=NO+NO2 (it is treated as if it were NO3=NO )
c    P1 =   remaining NO production terms
c    L1 =   loss of NO (except rxns producing NO2 - they are in R2_1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      R1_2  =            RKI(   1 )                    ! NO2=NO+O3P
     &        +          RKI(   5 ) * YC( O3P      )   ! NO2+O3P=NO
      R1_2  = R1_2 * DTC


      P1 =           RXRAT(  14 )      ! NO2+NO3=NO+NO2
     &    +          RXRAT(  15 )      ! NO3=NO
     &    +          RXRAT(  22 )      ! HONO=NO+HO
      P1    = YC0( NO ) + P1 * DTC


      L1    =            RKI(  21 ) * YC( HO       )   ! NO+HO=HONO
     &        +          RKI(  62 ) * YC( RO2_N    )   ! NO+RO2_N=RNO3
      L1    = 1.0 + L1 * DTC


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  NO2 Section
c    R2_1 = production of NO2 from NO ( rates of form k[NO][x] )
c            a)  NO+O3=NO2 not included
c            b)  NO+NO3=2NO2 ( 1/2 of NO2 formation rate included )
c            c)  NO3+NO2=NO+NO2 is not included for NO2
c    P2 =  remaining NO2 production terms 
c            a)  NO+O3=NO2 not included
c            b)  NO+NO3=2NO2 (1/2 of NO2 formation rate included )
c    L2 = loss of NO2 (except rxns producing NO2 - they are in R1_2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      R2_1  =            RKI(   4 ) * YC( O3P      )   ! NO+O3P=NO2
     &        +          RKI(   9 ) * YC( NO3      )   ! NO+NO3=2*NO2
     &        +  2.000 * RKI(  10 ) * YC( NO       )   ! NO+NO=2*NO2
     &        +          RKI(  31 ) * YC( HO2      )   ! NO+HO2=NO2+HO
     &        +          RKI(  46 ) * YC( C_O2     )   ! NO+C_O2=NO2+HCHO+HO2
     &        +          RKI(  51 ) * YC( RO2_R    )   ! NO+RO2_R=NO2+HO2
     &        +          RKI(  56 ) * YC( R2O2     )   ! NO+R2O2=NO2
     &        +          RKI(  71 ) * YC( CCO_O2   )   ! NO+CCO_O2=NO2+C_O2
     &        +          RKI(  81 ) * YC( RCO_O2   )   ! NO+RCO_O2=NO2+...
     &        +          RKI(  92 ) * YC( BZCO_O2  )   ! NO+BZCO_O2=NO2+...
     &        +          RKI( 104 ) * YC( MA_RCO3  )   ! NO+MA_RCO3=NO2+...
     &        +          RKI( 128 ) * YC( HOCOO    )   ! NO+HOCOO=NO2+...
      R2_1  = R2_1 * DTC


      P2 =           RXRAT(   9 )      ! NO+NO3=2*NO2
     &    +          RXRAT(  12 )      ! N2O5=NO2+NO3
     &    +          RXRAT(  16 )      ! NO3=NO2+O3P
     &    +          RXRAT(  23 )      ! HONO=NO2+HO2
     &    +          RXRAT(  24 )      ! HO+HONO=NO2
     &    +          RXRAT(  26 )      ! HO+NO3=NO2+HO2
     &    +          RXRAT(  28 )      ! HNO3=NO2+HO
     &    +          RXRAT(  33 )      ! HNO4=NO2+HO2
     &    +  0.610 * RXRAT(  34 )      ! HNO4=0.61*NO2+...
     &    +          RXRAT(  35 )      ! HNO4+HO=NO2
     &    +  0.800 * RXRAT(  39 )      ! NO3+HO2=0.8*NO2+...
     &    +  2.000 * RXRAT(  40 )      ! NO3+NO3=2*NO2
     &    +          RXRAT(  48 )      ! C_O2+NO3=NO2+HO2+HCHO
     &    +          RXRAT(  53 )      ! RO2_R+NO3=NO2+HO2
     &    +          RXRAT(  58 )      ! R2O2+NO3=NO2
     &    +          RXRAT(  65 )      ! RO2_N+NO3=NO2+HO2+MEK
     &    +          RXRAT(  70 )      ! PAN=NO2+CCO_O2
     &    +          RXRAT(  73 )      ! CCO_O2+NO3=NO2+C_O2
     &    +          RXRAT(  80 )      ! PAN2=NO2+RCO_O2
     &    +          RXRAT(  83 )      ! RCO_O2+NO3=NO2+CCHO+RO2_R
     &    +          RXRAT(  91 )      ! PBZN=NO2+BZCO_O2
     &    +          RXRAT(  94 )      ! BZCO_O2+NO3=NO2+BZ_O+R2O2
     &    +          RXRAT( 103 )      ! MA_PAN=NO2+MA_RCO3
     &    +          RXRAT( 106 )      ! MA_RCO3+NO3=NO2+HCHO+...
     &    +  0.338 * RXRAT( 176 )      ! RNO3+HO=0.338*NO2+...
     &    +          RXRAT( 177 )      ! RNO3=NO2+0.341*HO2+...
     &    +  0.187 * RXRAT( 191 )      ! ISOPRENE+NO3=...
     &    +  0.474 * RXRAT( 195 )      ! TRP1+NO3=0.474*NO2+...
     &    +  0.391 * RXRAT( 217 )      ! OLE2+NO3=0.391*NO2+...
      P2 = YC0( NO2 ) + P2 * DTC


      L2    =            RKI(   6 ) * YC( O3P      )   ! NO2+O3P=NO3
     &        +          RKI(   8 ) * YC( O3       )   ! NO2+O3=NO3
     &        +          RKI(  11 ) * YC( NO3      )   ! NO2+NO3=N2O5
     &        +          RKI(  25 ) * YC( HO       )   ! NO2+HO=HNO3
     &        +          RKI(  32 ) * YC( HO2      )   ! NO2+HO2=HNO4
     &        +          RKI(  69 ) * YC( CCO_O2   )   ! NO2+CCO_O2=PAN
     &        +          RKI(  79 ) * YC( RCO_O2   )   ! NO2+RCO_O2=PAN2
     &        +          RKI(  90 ) * YC( BZCO_O2  )   ! NO2+BZCO_O2=PBZN
     &        +          RKI( 102 ) * YC( MA_RCO3  )   ! NO2+MA_RCO3=MA_PAN
     &        +          RKI( 115 ) * YC( TBU_O    )   ! NO2+TBU_O=RNO3
     &        +          RKI( 117 ) * YC( BZ_O     )   ! NO2+BZ_O=NPHE
     &        +          RKI( 120 ) * YC( BZNO2_O  )   ! NO2+BZNO2_O=
      L2     = 1.0 + L2 * DTC


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  O3 Section
c    P3 = production of O3 except O+O2=O3
c    L3 =   loss terms for O3 except NO+O3=NO2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      P3 =   0.250 * RXRAT(  72 )      ! CCO_O2+HO2=0.25*O3+...
     &    +  0.250 * RXRAT(  82 )      ! RCO_O2+HO2=0.25*O3+...
     &    +  0.250 * RXRAT(  93 )      ! BZCO_O2+HO2=0.25*O3+...
     &    +  0.250 * RXRAT( 105 )      ! MA_RCO3+HO2=0.25*O3+...
      P3 = YC0( O3 ) + P3 * DTC


      L3    =            RKI(   3 ) * YC( O3P      )   ! O3+O3P=
     &        +          RKI(   8 ) * YC( NO2      )   ! O3+NO2=NO3
     &        +          RKI(  17 )                    ! O3=O3P
     &        +          RKI(  18 )                    ! O3=O1D2
     &        +          RKI(  30 ) * YC( HO       )   ! O3+HO=HO2
     &        +          RKI(  36 ) * YC( HO2      )   ! O3+HO2=HO
     &        +          RKI( 162 ) * YC( METHACRO )   ! O3+METHACRO=...
     &        +          RKI( 167 ) * YC( MVK      )   ! O3+MVK=0.064*HO2+...
     &        +          RKI( 171 ) * YC( ISOPROD  )   ! O3+ISOPROD=...
     &        +          RKI( 179 ) * YC( DCB1     )   ! O3+DCB1=1.5*HO2+...
     &        +          RKI( 186 ) * YC( ETHENE   )   ! O3+ETHENE=...
     &        +          RKI( 190 ) * YC( ISOPRENE )   ! O3+ISOPRENE=...
     &        +          RKI( 194 ) * YC( TRP1     )   ! O3+TRP1=0.567*HO+...
     &        +          RKI( 212 ) * YC( OLE1     )   ! O3+OLE1=0.155*HO+...
     &        +          RKI( 216 ) * YC( OLE2     )   ! O3+OLE2=0.378*HO+...
      L3    = 1.0 + L3 * DTC


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  O3P Section 
c    P12 = production of O3P except NO2+hv=O3P (J1)
c    L12 = loss terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      P12 =           RXRAT(  16 )      ! NO3=O3P+NO2
     &     +          RXRAT(  17 )      ! O3=O3P
     &     +  O3P_S * RXRAT(  18 )      ! O3=O1D2
      P12 = YC0( O3P ) + P12 * DTC


      L12   =            RKI(   2 )                    ! O3P=O3
     &        +          RKI(   3 ) * YC( O3       )   ! O3P+O3=
     &        +          RKI(   4 ) * YC( NO       )   ! O3P+NO=NO2
     &        +          RKI(   5 ) * YC( NO2      )   ! O3P+NO2=NO
     &        +          RKI(   6 ) * YC( NO2      )   ! O3P+NO2=NO3
     &        +          RKI( 164 ) * YC( METHACRO )   ! O3P+METHACRO=RCHO
     &        +          RKI( 168 ) * YC( MVK      )   ! O3P+MVK=0.45*RCHO+...
     &        +          RKI( 188 ) * YC( ETHENE   )   ! O3P+ETHENE=0.5*HO2+...
     &        +          RKI( 192 ) * YC( ISOPRENE )   ! O3P+ISOPRENE=...
     &        +          RKI( 196 ) * YC( TRP1     )   ! O3P+TRP1=0.147*RCHO+...
     &        +          RKI( 214 ) * YC( OLE1     )   ! O3P+OLE1=0.45*RCHO+...
     &        +          RKI( 218 ) * YC( OLE2     )   ! O3P+OLE2=0.013*HO2+...
      L12   = 1.0 + L12 * DTC

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Solution section
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c..compute reciprocal of loss terms
      L1_INV  = 1.0 / L1
      L2_INV  = 1.0 / L2
      L3_INV  = 1.0 / L3
      L12_INV = 1.0 / L12

c..compute specific k*delta t terms
      RK1 = RKI(   1 ) * DTC            ! J1    (NO2+hv=NO+O3P)
      RK2 = RKI(   2 ) * DTC            ! J2    (O3P+O2=O3)
      RK3 = RKI(   7 ) * DTC            ! k1_3  (NO+O3=NO2)

c..compute terms that are used to calulate a,b & c
      T1 = RK1  * L2_INV                ! J1   / ( 1.0 + Lno2 * dt )
      T2 = R1_2 * L2_INV                ! r1,2 / ( 1.0 + Lno2 * dt)
      T3 = R2_1 * L1_INV                ! r2,1 / ( 1.0 + Lno  * dt)
      T4 = RK2  * L12_INV               ! J2   / ( 1.0 + Lo3p * dt )
      T5 = T3   * P1 - T2 * P2          ! T3 * Pno - T2 * Pno2

      F1 = 1.0 + T2 + T3                ! factor in calculating a & b
      F2 = T1 * T4                      ! factor in calculating a & b
      F3 = L3 * L1 + RK3 * P1           ! (1 + Lo3 * dt) (1 + lno * dt )
                                        ! + k1,3 * dt * Pno

      PO3 = P3 + P12 * T4

      A = RK3 * ( F1  - F2 )

      B = F1 * F3 +  RK3 * ( F2 * ( P2 - P1 ) + PO3 +  T5 )

      C = RK3 * P1 * ( PO3 + P2 * F2 ) + F3 * T5

      Q = -0.5 * ( B + SIGN( 1.0, B ) * SQRT( B * B - 4.0 * A * C ) )

      XX = MAX( Q / A , C / Q  )


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Species solutions
c       [NO]   = ( P1 + x ) / ( 1 + L1 )
c       [NO2]  = ( P2 - x ) / ( 1 + L2 )
c       [O3 ]  = ( P3 + Ko3p->O3 ) / (1 + K1,3 * [NO] + L3 )
c       [O3P]  = ( P12 + J1 * [NO2] ) / ( 1 + L12 )
c       [O1D2] = ( Ko3->o1d * [O3] ) / ( O1D_DNM )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      YCP( NO ) = MAX( 0.0, ( P1 + XX ) * L1_INV )

      YCP( NO2 ) = MAX( 0.0, ( P2 - XX ) * L2_INV )

      S1 = P12 + RK1 * YCP( NO2 )

      S2 = T4 * S1

      YCP( O3 ) = ( P3 + S2 ) / ( L3 + RK3 * YCP( NO ) )

      YCP( O3P ) = S1 * L12_INV

      YCP( O1D2 ) = RKI( 18 ) * YCP( O3 ) / O1D_S

      RETURN

      END


