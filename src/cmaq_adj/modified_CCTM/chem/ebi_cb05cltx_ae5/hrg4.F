
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ebi_cb05cltx_ae5/hrg4.F,v 1.1.1.1 2010/06/14 16:03:01 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE HRG4( DTC )

C**********************************************************************
C
C  FUNCTION:  To solve for the concentration of NO3 and N2O5
C
C  PRECONDITIONS: For the CB05CLTX family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, May 30, 2008
C
C**********************************************************************
      USE HRDATA

      IMPLICIT NONE


C..INCLUDES: NONE


C..ARGUMENTS:
      REAL DTC             ! Time step


C..PARAMETERS: NONE


C..EXTERNAL FUNCTIONS: NONE


C..SAVED LOCAL VARIABLES:
      CHARACTER( 16 ), SAVE   ::  PNAME = 'HRG4'    ! Program name


C..SCRATCH LOCAL VARIABLES:
      REAL   A, B, C, Q   ! Quadratic equation terms
      REAL   CMN          ! Temp scalar
      REAL   L15          ! Loss of NO3
      REAL   L16          ! Loss of N2O5
      REAL   P15          ! Production of NO3
      REAL   K15_15       ! Kno3+no3 * delta t
      REAL   R15_16       ! Kn2o5-->no3 * delta t
      REAL   R16_15       ! Kno3+no2-->n2o5[NO2] * delta t


C**********************************************************************


c..Production of NO3 (except from N2O5 )
      P15 =           RXRAT(   5 )      ! O+NO2=NO3
     &     +          RXRAT(   7 )      ! NO2+O3=NO3
     &     +          RXRAT(  29 )      ! OH+HNO3=NO3
     &     +  0.390 * RXRAT(  51 )      ! PNA=0.39*NO3+0.61*NO2+...

c..Loss frequency of NO3 ( except NO3 + NO3 if present )
      L15   =            RKI(  14 )                            ! NO3=NO2+O
     &        +          RKI(  15 )                            ! NO3=NO
     &        +          RKI(  16 ) * YC( NO               )   ! NO3+NO=2*NO2
     &        +          RKI(  17 ) * YC( NO2              )   ! NO3+NO2=NO+NO2
     &        +          RKI(  18 ) * YC( NO2              )   ! NO3+NO2=N2O5
     &        +          RKI(  46 ) * YC( O                )   ! NO3+O=NO2
     &        +          RKI(  47 ) * YC( OH               )   ! NO3+OH=HO2+NO2
     &        +          RKI(  48 ) * YC( HO2              )   ! NO3+HO2=HNO3
     &        +          RKI(  49 ) * YC( O3               )   ! NO3+O3=NO2
     &        +          RKI(  77 ) * YC( FORM             )   ! NO3+FORM=HNO3+HO2+CO
     &        +          RKI(  85 ) * YC( ALD2             )   ! NO3+ALD2=C2O3+HNO3
     &        +          RKI( 100 ) * YC( ALDX             )   ! NO3+ALDX=CXO3+HNO3
     &        +          RKI( 119 ) * YC( OLE              )   ! NO3+OLE=NO2+FORM+...
     &        +          RKI( 123 ) * YC( ETH              )   ! NO3+ETH=NO2+XO2+...
     &        +          RKI( 127 ) * YC( IOLE             )   ! NO3+IOLE=...
     &        +          RKI( 132 ) * YC( CRES             )   ! NO3+CRES=CRO+HNO3
     &        +          RKI( 144 ) * YC( ISOP             )   ! NO3+ISOP=...
     &        +          RKI( 147 ) * YC( ISPD             )   ! NO3+ISPD=...
     &        +          RKI( 152 ) * YC( TERP             )   ! NO3+TERP=...

c..Loss frequency of N2O5
      L16   =            RKI(  19 )                            ! N2O5=2*HNO3
     &        +          RKI(  20 )                            ! N2O5=2*HNO3
     &        +          RKI(  21 )                            ! N2O5=NO3+NO2
     &        +          RKI(  53 )                            ! N2O5=NO2+NO3

c..K15_15, R15_16, and R16_15 terms
      K15_15  = RKI(  50 ) * DTC

      R15_16  = ( RKI(  21 )
     &        +   RKI(  53 ) ) * DTC 


      R16_15  = RKI(  18 ) * YC( NO2 ) * DTC

c..Solution of quadratic equation to get NO3 & N2O5
      CMN = 1.0 + L16 * DTC
      A = 2.0 * K15_15 * CMN
      B = CMN * ( 1.0 + L15 * DTC ) - R15_16 * R16_15
      C = CMN * ( YC0( NO3 ) + P15 * DTC ) +  R15_16 * YC0( N2O5 )

      Q = -0.5 * ( B + SIGN( 1.0, B ) * SQRT( B * B + 4.0 * A * C ) )
      YCP( NO3 ) = MAX( Q / A , -C / Q  )
      YCP( N2O5 ) = ( YC0( N2O5 ) + R16_15 * YCP( NO3 ) ) / CMN

      RETURN

      END
