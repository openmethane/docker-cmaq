
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_saprc99_ae5/hrg3.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE HRG3( DTC )

C**********************************************************************
C
C  FUNCTION:  To solve for the concentration of C2O3 and PAN
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

C..INCLUDES:  NONE

C..ARGUMENTS:
      REAL DTC              ! Time step


C..PARAMETERS: NONE


C..EXTERNAL FUNCTIONS: NONE


C..SAVED LOCAL VARIABLES:
      CHARACTER( 16 ), SAVE  :: PNAME = 'HRG3'      ! Program name


C..SCRATCH LOCAL VARIABLES:
      REAL   A, B, C, Q   ! Quadratic equation terms
      REAL   CMN          ! Temp scalar
      REAL   L8           ! Loss of CCO_O2
      REAL   L9           ! Loss of PAN
      REAL   P8           ! Production of CCO_O2

      REAL   K8_8         ! Kcco_o2+cco_o2 * delta t
      REAL   R8_9         ! Kpan-->cco_o2 * delta t
      REAL   R9_8         ! Kcco_o2+no2-->pan * [NO2] * delta t

C**********************************************************************


c..Production of CCO_O2 (except from PAN )
      P8 =           RXRAT( 104 )      ! MA_RCO3+NO=CCO_O2+HCHO+NO2
     &    +          RXRAT( 106 )      ! MA_RCO3+NO3=CCO_O2+HCHO+NO2
     &    +          RXRAT( 112 )      ! MA_RCO3+RCO_O2=CCO_O2+HCHO+...
     &    +          RXRAT( 113 )      ! MA_RCO3+BZCO_O2=CCO_O2+...
     &    +  2.000 * RXRAT( 114 )      ! MA_RCO3+MA_RCO3=2*CCO_O2+...
     &    +          RXRAT( 130 )      ! CCHO+HO=CCO_O2
     &    +          RXRAT( 132 )      ! CCHO+NO3=CCO_O2+HNO3
     &    +          RXRAT( 136 )      ! ACET+HO=CCO_O2+HCHO+R2O2
     &    +          RXRAT( 137 )      ! ACET=CCO_O2+C_O2
     &    +  0.492 * RXRAT( 138 )      ! MEK+HO=0.492*CCO_O2+...
     &    +          RXRAT( 139 )      ! MEK=CCO_O2+CCHO+RO2_R
     &    +          RXRAT( 149 )      ! MGLY=CCO_O2+CO+HO2
     &    +          RXRAT( 150 )      ! MGLY+HO=CCO_O2+CO
     &    +          RXRAT( 151 )      ! MGLY+NO3=CCO_O2+CO+HNO3
     &    +  2.000 * RXRAT( 152 )      ! BACL=2*CCO_O2
     &    +  0.670 * RXRAT( 165 )      ! METHACRO=0.67*CCO_O2+...
     &    +  0.675 * RXRAT( 166 )      ! MVK+HO=0.675*CCO_O2+...
     &    +  0.467 * RXRAT( 173 )      ! ISOPROD=0.467*CCO_O2+...
     &    +  0.029 * RXRAT( 174 )      ! PROD2+HO=0.029*CCO_O2+...
     &    +  0.667 * RXRAT( 175 )      ! PROD2=0.667*CCO_O2+...
     &    +          RXRAT( 180 )      ! DCB2+HO=CCO_O2+RCHO+R2O2
     &    +  0.500 * RXRAT( 181 )      ! DCB2=0.5*CCO_O2+RO2_R+...
     &    +          RXRAT( 182 )      ! DCB3+HO=CCO_O2+RCHO+R2O2
     &    +  0.500 * RXRAT( 183 )      ! DCB3=0.5*CCO_O2+RO2_R+...
     &    +  0.123 * RXRAT( 194 )      ! TRP1+O3=0.123*CCO_O2+...
     &    +  0.011 * RXRAT( 200 )      ! ALK4+HO=0.011*CCO_O2+...
     &    +  0.137 * RXRAT( 216 )      ! OLE2+O3=0.137*CCO_O2+...

c..Loss frequency of CCO_O2 ( not including CCO_O2 + CCO_O2 )
      L8    =            RKI(  69 ) * YC( NO2      )   ! CCO_O2+NO2=PAN
     &        +          RKI(  71 ) * YC( NO       )   ! CCO_O2+NO=C_O2+NO2
     &        +          RKI(  72 ) * YC( HO2      )   ! CCO_O2+HO2=...
     &        +          RKI(  73 ) * YC( NO3      )   ! CCO_O2+NO3=C_O2+NO2
     &        +          RKI(  74 ) * YC( C_O2     )   ! CCO_O2+C_O2=...
     &        +          RKI(  75 ) * YC( RO2_R    )   ! CCO_O2+RO2_R=CCO_OH
     &        +          RKI(  77 ) * YC( RO2_N    )   ! CCO_O2+RO2_N=...
     &        +          RKI(  88 ) * YC( RCO_O2   )   ! CCO_O2+RCO_O2=...
     &        +          RKI(  99 ) * YC( BZCO_O2  )   ! CCO_O2+BZCO_O2=...

c..Loss frequency of PAN
      L9    =            RKI(  70 )                    ! PAN=CCO_O2+NO2

c..K8_8, R8_9, and R9_8 terms
      K8_8  = RKI(  78 ) * DTC

      R8_9  = ( RKI(  70 ) ) * DTC 

      R9_8  = ( RKI(  69 ) * YC( NO2 ) ) * DTC 

c..Solution of quadratic equation to get CCO_O2 & PAN
      CMN = 1.0 + L9 * DTC
      A = 2.0 * K8_8 * CMN
      B = CMN * ( 1.0 + L8 * DTC ) - R8_9 * R9_8
      C = CMN * ( YC0( CCO_O2 ) + P8 * DTC ) +  R8_9 * YC0( PAN )

      Q = -0.5 * ( B + SIGN( 1.0, B ) * SQRT( B * B + 4.0 * A * C ) )

      YCP( CCO_O2 ) = MAX( Q / A , -C / Q  )

      YCP( PAN ) = ( YC0( PAN ) +  R9_8 * YCP( CCO_O2 ) ) / CMN

      RETURN

      END
