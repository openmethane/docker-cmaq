
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/phot/phot_inline/AERO_DESC.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE AERO_DESC

C-----------------------------------------------------------------------
C
C  FUNCTION:  MODULE AERO_DESC contains aerosol module descriptive
C     information, including species pointers for CBLK, aerosol relevant
C     constants, and methods (e.g. GETNEWPAR2) that operate on the
C     aerosol species
C
C  REVISION  HISTORY:
C       Date   Who          What
C     -------- ---------    -----------------------------------------
C     01/2008  S.Roselle    Adapted from module PHOTMOD12 for AE5
C     05/2008  S.Roselle    Added code for new species ANH4K and SRFCOR;
C                           Also added code for variable coarse mode
C                           standard deviation (no longer fixed at 2.2)
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C***Parameters

      INTEGER, PARAMETER :: NMODE = 3      ! number of aerosol modes (3)

      INTEGER, PARAMETER :: N_AE_MAP = 55  ! # species mapped from CGRID to CBLK
      INTEGER, PARAMETER :: NSPCSDA = N_AE_MAP + 15  ! total # species in CBLK

      INTEGER, PARAMETER :: VSO4AI  = 01   ! CBLK pointer to ASO4I
      INTEGER, PARAMETER :: VNH4AI  = 02   ! CBLK pointer to ANH4I
      INTEGER, PARAMETER :: VNO3AI  = 03   ! CBLK pointer to ANO3I
      INTEGER, PARAMETER :: VORGAI  = 04   ! CBLK pointer to AORGAI  (aero4 only)
      INTEGER, PARAMETER :: VORGPAI = 05   ! CBLK pointer to AORGPAI
      INTEGER, PARAMETER :: VORGBAI = 06   ! CBLK pointer to AORGBI  (aero4 only)
      INTEGER, PARAMETER :: VECI    = 07   ! CBLK pointer to AECI
      INTEGER, PARAMETER :: VP25AI  = 08   ! CBLK pointer to A25I
      INTEGER, PARAMETER :: VH2OAI  = 09   ! CBLK pointer to AH2OI
      INTEGER, PARAMETER :: VNAI    = 10   ! CBLK pointer to ANAI
      INTEGER, PARAMETER :: VCLI    = 11   ! CBLK pointer to ACLI
      INTEGER, PARAMETER :: VSO4AJ  = 12   ! CBLK pointer to ASO4J
      INTEGER, PARAMETER :: VNH4AJ  = 13   ! CBLK pointer to ANH4J
      INTEGER, PARAMETER :: VNO3AJ  = 14   ! CBLK pointer to ANO3J
      INTEGER, PARAMETER :: VORGAJ  = 15   ! CBLK pointer to AORGAJ  (aero4 only)
      INTEGER, PARAMETER :: VORGPAJ = 16   ! CBLK pointer to AORGPAJ
      INTEGER, PARAMETER :: VORGBAJ = 17   ! CBLK pointer to AORGBJ  (aero4 only)
      INTEGER, PARAMETER :: VALKJ   = 18   ! CBLK pointer to ALKJ    (aero5 only)
      INTEGER, PARAMETER :: VXYL1J  = 19   ! CBLK pointer to XYL1J   (aero5 only)
      INTEGER, PARAMETER :: VXYL2J  = 20   ! CBLK pointer to XYL2J   (aero5 only)
      INTEGER, PARAMETER :: VXYL3J  = 21   ! CBLK pointer to XYL3J   (aero5 only)
      INTEGER, PARAMETER :: VTOL1J  = 22   ! CBLK pointer to TOL1J   (aero5 only)
      INTEGER, PARAMETER :: VTOL2J  = 23   ! CBLK pointer to TOL2J   (aero5 only)
      INTEGER, PARAMETER :: VTOL3J  = 24   ! CBLK pointer to TOL3J   (aero5 only)
      INTEGER, PARAMETER :: VBNZ1J  = 25   ! CBLK pointer to BNZ1J   (aero5 only)
      INTEGER, PARAMETER :: VBNZ2J  = 26   ! CBLK pointer to BNZ2J   (aero5 only)
      INTEGER, PARAMETER :: VBNZ3J  = 27   ! CBLK pointer to BNZ3J   (aero5 only)
      INTEGER, PARAMETER :: VTRP1J  = 28   ! CBLK pointer to TRP1J   (aero5 only)
      INTEGER, PARAMETER :: VTRP2J  = 29   ! CBLK pointer to TRP2J   (aero5 only)
      INTEGER, PARAMETER :: VISO1J  = 30   ! CBLK pointer to ISO1J   (aero5 only)
      INTEGER, PARAMETER :: VISO2J  = 31   ! CBLK pointer to ISO2J   (aero5 only)
      INTEGER, PARAMETER :: VISO3J  = 32   ! CBLK pointer to ISO3J   (aero5 only)
      INTEGER, PARAMETER :: VSQTJ   = 33   ! CBLK pointer to SQTJ    (aero5 only)
      INTEGER, PARAMETER :: VOLGAJ  = 34   ! CBLK pointer to OLGAJ   (aero5 only)
      INTEGER, PARAMETER :: VOLGBJ  = 35   ! CBLK pointer to OLGBJ   (aero5 only)
      INTEGER, PARAMETER :: VORGCJ  = 36   ! CBLK pointer to ORGCJ   (aero5 only)
      INTEGER, PARAMETER :: VECJ    = 37   ! CBLK pointer to AECJ
      INTEGER, PARAMETER :: VP25AJ  = 38   ! CBLK pointer to A25J
      INTEGER, PARAMETER :: VH2OAJ  = 39   ! CBLK pointer to AH2OJ
      INTEGER, PARAMETER :: VNAJ    = 40   ! CBLK pointer to ANAJ
      INTEGER, PARAMETER :: VCLJ    = 41   ! CBLK pointer to ACLJ
      INTEGER, PARAMETER :: VANTHA  = 42   ! CBLK pointer to ACORS
      INTEGER, PARAMETER :: VSOILA  = 43   ! CBLK pointer to ASOIL
      INTEGER, PARAMETER :: VNAK    = 44   ! CBLK pointer to ANAK
      INTEGER, PARAMETER :: VCLK    = 45   ! CBLK pointer to ACLK
      INTEGER, PARAMETER :: VSO4K   = 46   ! CBLK pointer to ASO4K
      INTEGER, PARAMETER :: VNH4K   = 47   ! CBLK pointer to ANH4K
      INTEGER, PARAMETER :: VNO3K   = 48   ! CBLK pointer to ANO3K
      INTEGER, PARAMETER :: VH2OK   = 49   ! CBLK pointer to AH2OK
      INTEGER, PARAMETER :: VAT0    = 50   ! CBLK pointer to NUMATKN
      INTEGER, PARAMETER :: VAC0    = 51   ! CBLK pointer to NUMACC
      INTEGER, PARAMETER :: VCOR0   = 52   ! CBLK pointer to NUMCOR
      INTEGER, PARAMETER :: VSURFAT = 53   ! CBLK pointer to SRFATKN
      INTEGER, PARAMETER :: VSURFAC = 54   ! CBLK pointer to SRFACC
      INTEGER, PARAMETER :: VSURFCO = 55   ! CBLK pointer to SRFCOR
      INTEGER, PARAMETER :: VAT2    = 56   ! CBLK pointer to AT2
      INTEGER, PARAMETER :: VAC2    = 57   ! CBLK pointer to AC2
      INTEGER, PARAMETER :: VCOR2   = 58   ! CBLK pointer to AC2
      INTEGER, PARAMETER :: VAT3    = 59   ! CBLK pointer to AT3
      INTEGER, PARAMETER :: VAC3    = 60   ! CBLK pointer to AC3
      INTEGER, PARAMETER :: VCOR3   = 61   ! CBLK pointer to COR3
      INTEGER, PARAMETER :: VDGAT   = 62   ! CBLK pointer to DGAT
      INTEGER, PARAMETER :: VDGAC   = 63   ! CBLK pointer to DGAC
      INTEGER, PARAMETER :: VDGCO   = 64   ! CBLK pointer to DGCO
      INTEGER, PARAMETER :: VLSGAT  = 65   ! CBLK pointer to LSGAT
      INTEGER, PARAMETER :: VLSGAC  = 66   ! CBLK pointer to LSGAC
      INTEGER, PARAMETER :: VLSGCO  = 67   ! CBLK pointer to LSGCO
      INTEGER, PARAMETER :: VOLAT   = 68   ! CBLK pointer to VOLAT
      INTEGER, PARAMETER :: VOLAC   = 69   ! CBLK pointer to VOLAC
      INTEGER, PARAMETER :: VOLCO   = 70   ! CBLK pointer to VOLCO

C***Names of species in CBLK

      CHARACTER(16), SAVE :: CBLK_VARNM( N_AE_MAP )

      DATA CBLK_VARNM( VSO4AI  ) /'ASO4I           ' /
      DATA CBLK_VARNM( VNH4AI  ) /'ANH4I           ' /
      DATA CBLK_VARNM( VNO3AI  ) /'ANO3I           ' /
      DATA CBLK_VARNM( VORGAI  ) /'AORGAI          ' /
      DATA CBLK_VARNM( VORGPAI ) /'AORGPAI         ' /
      DATA CBLK_VARNM( VORGBAI ) /'AORGBI          ' /
      DATA CBLK_VARNM( VECI    ) /'AECI            ' /
      DATA CBLK_VARNM( VP25AI  ) /'A25I            ' /
      DATA CBLK_VARNM( VH2OAI  ) /'AH2OI           ' /
      DATA CBLK_VARNM( VNAI    ) /'ANAI            ' /
      DATA CBLK_VARNM( VCLI    ) /'ACLI            ' /
      DATA CBLK_VARNM( VSO4AJ  ) /'ASO4J           ' /
      DATA CBLK_VARNM( VNH4AJ  ) /'ANH4J           ' /
      DATA CBLK_VARNM( VNO3AJ  ) /'ANO3J           ' /
      DATA CBLK_VARNM( VORGAJ  ) /'AORGAJ          ' /
      DATA CBLK_VARNM( VORGPAJ ) /'AORGPAJ         ' /
      DATA CBLK_VARNM( VORGBAJ ) /'AORGBJ          ' /
      DATA CBLK_VARNM( VALKJ   ) /'AALKJ           ' /
      DATA CBLK_VARNM( VXYL1J  ) /'AXYL1J          ' /
      DATA CBLK_VARNM( VXYL2J  ) /'AXYL2J          ' /
      DATA CBLK_VARNM( VXYL3J  ) /'AXYL3J          ' /
      DATA CBLK_VARNM( VTOL1J  ) /'ATOL1J          ' /
      DATA CBLK_VARNM( VTOL2J  ) /'ATOL2J          ' /
      DATA CBLK_VARNM( VTOL3J  ) /'ATOL3J          ' /
      DATA CBLK_VARNM( VBNZ1J  ) /'ABNZ1J          ' /
      DATA CBLK_VARNM( VBNZ2J  ) /'ABNZ2J          ' /
      DATA CBLK_VARNM( VBNZ3J  ) /'ABNZ3J          ' /
      DATA CBLK_VARNM( VTRP1J  ) /'ATRP1J          ' /
      DATA CBLK_VARNM( VTRP2J  ) /'ATRP2J          ' /
      DATA CBLK_VARNM( VISO1J  ) /'AISO1J          ' /
      DATA CBLK_VARNM( VISO2J  ) /'AISO2J          ' /
      DATA CBLK_VARNM( VISO3J  ) /'AISO3J          ' /
      DATA CBLK_VARNM( VSQTJ   ) /'ASQTJ           ' /
      DATA CBLK_VARNM( VOLGAJ  ) /'AOLGAJ          ' /
      DATA CBLK_VARNM( VOLGBJ  ) /'AOLGBJ          ' /
      DATA CBLK_VARNM( VORGCJ  ) /'AORGCJ          ' /
      DATA CBLK_VARNM( VECJ    ) /'AECJ            ' /
      DATA CBLK_VARNM( VP25AJ  ) /'A25J            ' /
      DATA CBLK_VARNM( VH2OAJ  ) /'AH2OJ           ' /
      DATA CBLK_VARNM( VNAJ    ) /'ANAJ            ' /
      DATA CBLK_VARNM( VCLJ    ) /'ACLJ            ' /
      DATA CBLK_VARNM( VANTHA  ) /'ACORS           ' /
      DATA CBLK_VARNM( VSOILA  ) /'ASOIL           ' /
      DATA CBLK_VARNM( VNAK    ) /'ANAK            ' /
      DATA CBLK_VARNM( VCLK    ) /'ACLK            ' /
      DATA CBLK_VARNM( VSO4K   ) /'ASO4K           ' /
      DATA CBLK_VARNM( VNH4K   ) /'ANH4K           ' /
      DATA CBLK_VARNM( VNO3K   ) /'ANO3K           ' /
      DATA CBLK_VARNM( VH2OK   ) /'AH2OK           ' /
      DATA CBLK_VARNM( VAT0    ) /'NUMATKN         ' /
      DATA CBLK_VARNM( VAC0    ) /'NUMACC          ' /
      DATA CBLK_VARNM( VCOR0   ) /'NUMCOR          ' /
      DATA CBLK_VARNM( VSURFAT ) /'SRFATKN         ' /
      DATA CBLK_VARNM( VSURFAC ) /'SRFACC          ' /
      DATA CBLK_VARNM( VSURFCO ) /'SRFCOR          ' /

      REAL, PARAMETER :: CONMIN  = 1.0E-30 ! min concentration [ug/m**3]

C***(bulk) component densities [ kg/m**3 ]

      REAL, PARAMETER :: RHOSO4  = 1.8E3    ! aerosol sulfate
      REAL, PARAMETER :: RHONH4  = 1.8E3    ! aerosol ammonium
      REAL, PARAMETER :: RHONO3  = 1.8E3    ! aerosol nitrate
      REAL, PARAMETER :: RHOH2O  = 1.0E3    ! aerosol water
      REAL, PARAMETER :: RHOORG  = 2.0E3    ! aerosol organics
      REAL, PARAMETER :: RHOSOIL = 2.6E3    ! aerosol soil dust
      REAL, PARAMETER :: RHOSEAS = 2.2E3    ! marine aerosol
      REAL, PARAMETER :: RHOANTH = 2.2E3    ! anthropogenic aerosol

C***Initial sigma-G for each mode

      REAL, PARAMETER :: SGINIAT = 1.70     ! Aitken mode
      REAL, PARAMETER :: SGINIAC = 2.00     ! accumulation mode
      REAL, PARAMETER :: SGINICO = 2.2      ! coarse mode

C***Inital mean diameters [ m ] for each mode

      REAL, PARAMETER :: DGINIAT = 0.01E-6  ! Aitken mode
      REAL, PARAMETER :: DGINIAC = 0.07E-6  ! accumulation mode
      REAL, PARAMETER :: DGINICO = 1.00E-6  ! coarse mode

      CONTAINS

C//////////////////////////////////////////////////////////////////

      SUBROUTINE GETNEWPAR2 ( NSPCSDA, NLAYS, CBLK, NR, NI )
C-----------------------------------------------------------------------
C  FUNCTION:  This subroutine calculates the volume, geometric mean
C     diameters (Dg) and the natural logs of geometric standard
C     deviations (Sg) of all 3 modes.  The code also calculates the 
C     average modal aerosol refractive indices
C
C     Internally mixed wet aerosols are assumed.
C
C     NOTE
C       2nd moment concentrations (M2) are passed into GETPAR in the
C       CBLK array and are modified within GETPAR only in the event that
C       the Sg value of a given mode has gone outside of the acceptable
C       range (1.05 to 2.50).  The GETPAR calculations implicitly assume
C       that the input value of M2 is consistent with the input value of
C       M3_WET_FLAG.  If, for example, the input M2 value was calculated
C       for a "dry" aerosol and the M3_WET_FLAG is .TRUE., GETPAR would
C       incorrectly adjust the M2 concentrations!
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  NONE
C
C  REVISION HISTORY
C     Modified from GETPAR 10/13/2004
C       Dr. Francis S. Binkowski
C       Carolina Environmental Program
C       University of North Carolina at Chapel Hill
C       email:  frank_binkowski@unc.edu
C     01/29/2008 S.Roselle - Reformatted for inclusion in CMAQ
C     05/21/2008 S.Roselle - Added code for new species ANH4K and SRFCOR;
C                            Also added code for variable coarse mode
C                            standard deviation (no longer fixed at 2.2)
C
C  REFERENCES:
C     1. Binkowski, F.S. Aerosols in Models-3 CMAQ, Chapter 10 of Science
C     Algorithms of the EPA Models-3 Community Multiscale Air Quality
C     (CMAQ) Modeling System, EPA/R-99/030, March 1999.
C     Available at: http://www.epa.gov/asmdnerl/models3
C
C     2. Binkowski, F.S. and S.J. Roselle, Models-3 Community
C     Multiscale Air Quality (CMAQ) model aerosol component 1:
C     Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
C     doi:10.1029/2001JD001409, 2003.
C
C     3. Horvath, H., Size segrated light absorption coefficient for
C     the atmospheric aerosol, Atmos. Environ., Vol. 29, No. 8,
C     pp 875-883, 1995.
C
C     Beyer, K.D. A.R. Ravishankara, and E. R. Lovejoy, Meanurements
C     of H2SO4/H2O and H2SO4/HNO3/H2O solutions, J. Geophys. Res. Vol.
C     101, D9, pp 14519-14524, 1996.
C
C *** coded 02/11/2004 by Dr. Francis S. Binkowski
C     Carolina Environmental Program
C     University of North Carolina at Chapel Hill
C     email: frank_binkowski@unc.edu
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C...include files

      INCLUDE SUBST_CONST      ! physical constants
      INCLUDE SUBST_RXCMMN     ! chemical mechamism reactions COMMON

C***Argument list

      INTEGER, INTENT(IN) :: NSPCSDA         ! number of species in CBLK
      INTEGER, INTENT(IN) :: NLAYS           ! number of layers in CBLK

      REAL, INTENT(INOUT) :: CBLK( NSPCSDA, NLAYS ) ! main array of variables

      REAL, INTENT(OUT)   :: NR( NMODE, NLAYS )  ! real refractive index for 3 modes
      REAL, INTENT(OUT)   :: NI( NMODE, NLAYS )  ! imaginary refactive index for 3 modes

C***Parameters

      CHARACTER(16), PARAMETER :: PNAME = 'GETNEWPAR2'

      REAL,    PARAMETER :: DGMIN   = 1.0E-09        ! min particle diam [m]
      REAL,    PARAMETER :: PIM1    = 1.0 / PI       ! 1/PI
      REAL,    PARAMETER :: PI6     = PI / 6.0       ! PI/6
      REAL,    PARAMETER :: F6DPI   = 6.0 / PI       !  6/PI
      REAL,    PARAMETER :: F6DPIM9 = 1.0E-9 * F6DPI ! 1.0e-9 * 6/PI
      REAL,    PARAMETER :: ONE3    = 1.0 / 3.0      ! 1/3
      REAL(8), PARAMETER :: ONE3D   = 1.0D0 / 3.0D0  ! 1/3
      REAL(8), PARAMETER :: TWO3D   = 2.0D0 * ONE3D  ! 2/3
      REAL(8), PARAMETER :: MINL2SG = 2.38048048D-03 ! min value ln(Sg)**2 (sg=1.05)
      REAL(8), PARAMETER :: MAXL2SG = 8.39588705D-01 ! max value ln(Sg)**2 (sg=2.5)

C***Aerosol version name

      CHARACTER(16), SAVE :: AE_VRSN

      CHARACTER(120) :: XMSG = ' '

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

C***Factors for converting aerosol mass concentration [ ug m**-3] to
C***  to 3rd moment concentration [ m**3 m^-3]

      REAL,    PARAMETER :: SO4FAC  = F6DPIM9 / RHOSO4
      REAL,    PARAMETER :: NH4FAC  = F6DPIM9 / RHONH4
      REAL,    PARAMETER :: H2OFAC  = F6DPIM9 / RHOH2O
      REAL,    PARAMETER :: NO3FAC  = F6DPIM9 / RHONO3
      REAL,    PARAMETER :: ORGFAC  = F6DPIM9 / RHOORG
      REAL,    PARAMETER :: SOILFAC = F6DPIM9 / RHOSOIL
      REAL,    PARAMETER :: SEASFAC = F6DPIM9 / RHOSEAS
      REAL,    PARAMETER :: ANTHFAC = F6DPIM9 / RHOANTH

C***Internal variables

      INTEGER LAY               ! loop variable

C***Geometric mean diameters [m]

      REAL DGATK                ! Aitken mode
      REAL DGACC                ! accumulation mode
      REAL DGCOR                ! coarse mode

C***logarithms of geometric standard deviations

      REAL XXLSGAT              ! Aitken mode
      REAL XXLSGAC              ! accumulation mode
      REAL XXLSGCO              ! coarse mode

      REAL ESAT36, ESAC36, ESC36  ! exp(4.5*L2SG); used in diameter calcs

      REAL(8) XXM0, XXM2, XXM3  ! temporary storage of moment conc's
      REAL(8) XFSUM             ! (ln(M0)+2ln(M3))/3; used in Sg calcs
      REAL(8) LXFM2             ! ln(M2); used in Sg calcs
      REAL(8) L2SGAT, L2SGAC, L2SGCO   ! square of ln(Sg); used in diameter calcs

C***Except as otherwise noted reference values of refractive indices for
C***  aerosol particles values for water soluble, insoluble, sea salt,
C***  and water are from the OPAC Data base.
C***
C***  Hess, M. P. Koepke, and I. Schult, Optical properties of aerosols
C***     and clouds: The software package OPAC, Bulletin of the American
C***     Meteorological Society, Vol 79, No 5, pp 831-844, May 1998.
C***
C***  OPAC is a downloadable data set of optical properties of 10 aerosol
C***  components, 6 water clouds and 3 cirrus clouds at UV, visible and
C***  IR wavelengths:  www.lrz-muenchen.de/~uh234an/www/radaer/opac.htm

C***FSB Water soluble values are modeled on data from Beyer et al. (1996)

      REAL, SAVE :: XNREAL_WS = 1.5  ! water soluble
      REAL, SAVE :: XNIMAG_WS = 0.0

      REAL, SAVE :: XNREAL_IN = 1.53 ! dust like
      REAL, SAVE :: XNIMAG_IN = 0.01 !

      REAL, SAVE :: XNREAL_SS = 1.50 ! sea salt
      REAL, SAVE :: XNIMAG_SS = 0.0

      REAL, SAVE :: XNREAL_H2 = 1.34 ! water
      REAL, SAVE :: XNIMAG_H2 = 0.0

C***FSB The following values for soot like particles are modified from
C***  Table 1 of Horvath (1995). They are size dependent effective
C***  refractive indices fit to data on absorption by carbon

      REAL, SAVE :: XNREAL_ST1 = 1.53
      REAL, SAVE :: XNIMAG_ST1 = 0.14
      REAL, SAVE :: XNREAL_ST2 = 1.53
      REAL, SAVE :: XNIMAG_ST2 = 0.05

C***modal volumes by groups

      REAL VWSAT                ! water soluble Aitken mode
      REAL VINAT                ! insoluble Aitken mode
      REAL VSTAT                ! carbon soot Aitken mode
      REAL VH2AT                ! water Aitken mode

      REAL VWSAC                ! water soluble accumulation mode
      REAL VINAC                ! insolubleaccumulation mode
      REAL VSTAC                ! carbon soot accumulation mode
      REAL VH2AC                ! water accumulation mode

      REAL VWSCO                ! water soluble coarse mode
      REAL VINCO                ! insoluble coarse mode
      REAL VSSCO                ! seasalt coarse mode
      REAL VH2CO                ! water coarse mode
      REAL VSUMAT, VSUMAC, VSUMCO

C------------------ *** Begin solution code *** ----------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         IF ( INDEX ( MECHNAME, 'AE4' ) .GT. 0 ) THEN
            AE_VRSN  = 'AE4'
         ELSE IF ( INDEX ( MECHNAME, 'AE5' ) .GT. 0 ) THEN
            AE_VRSN  = 'AE5'
         ELSE
            XMSG = 'This routine can only be used with the AE4 '
     &          // 'or AE5 aerosol mechanisms.'
            CALL M3EXIT( PNAME, 0, 0, XMSG, 1 )
         END IF ! check on MECHNAME

      END IF    ! FIRSTIME

C***Loop over layers

      DO LAY = 1, NLAYS

C***FSB:  In the following calculations the Pi/6 factor is omitted from
C***  calculating volume from from third moment because only ratios of
C***  volumes are needed when calculating the volume-averaged refractive
C***  indices, and the factor cancels out.

C***Assemble the third moments

C***Aitken mode
C***  Water soluble

         VWSAT = MAX( CONMIN,
     &              ( SO4FAC  * CBLK( VSO4AI, LAY ) +
     &                NH4FAC  * CBLK( VNH4AI, LAY ) +
     &                NO3FAC  * CBLK( VNO3AI, LAY ) +
     &                SEASFAC * CBLK( VNAI,   LAY ) +
     &                SEASFAC * CBLK( VCLI,   LAY ) ) )

C***  insoluble

         IF ( AE_VRSN .EQ. 'AE5' ) THEN

            VINAT = MAX( CONMIN,
     &                 ( ORGFAC  * CBLK( VORGPAI, LAY ) +
     &                   ANTHFAC * CBLK( VP25AI,  LAY ) ) )

         ELSE  ! AE4

            VINAT = MAX( CONMIN,
     &                 ( ORGFAC  * CBLK( VORGPAI, LAY ) +
     &                   ANTHFAC * CBLK( VP25AI,  LAY ) +
     &                   ORGFAC  * CBLK( VORGAI,  LAY ) +
     &                   ORGFAC  * CBLK( VORGBAI, LAY ) ) )

         END IF  ! AE_VRSN

C***  soot like

         VSTAT = MAX( CONMIN,
     &                ANTHFAC * CBLK( VECI, LAY ) )

C***  water

         VH2AT = MAX( CONMIN,
     &                H2OFAC * CBLK( VH2OAI, LAY ) )

         VSUMAT = VWSAT + VINAT + VSTAT + VH2AT

C***accumulation mode
C***  Water soluble

         VWSAC = MAX( CONMIN,
     &              ( SO4FAC  * CBLK( VSO4AJ, LAY ) +
     &                NH4FAC  * CBLK( VNH4AJ, LAY ) +
     &                NO3FAC  * CBLK( VNO3AJ, LAY ) +
     &                SEASFAC * CBLK( VNAJ,   LAY ) +
     &                SEASFAC * CBLK( VCLJ,   LAY ) ) )

C***  insoluble

         IF ( AE_VRSN .EQ. 'AE5' ) THEN

            VINAC = MAX( CONMIN,
     &                 ( ANTHFAC * CBLK( VP25AJ,  LAY ) +
     &                   ORGFAC  * CBLK( VORGPAJ, LAY ) +
     &                   ORGFAC  * CBLK( VALKJ,   LAY ) +
     &                   ORGFAC  * CBLK( VXYL1J,  LAY ) +
     &                   ORGFAC  * CBLK( VXYL2J,  LAY ) +
     &                   ORGFAC  * CBLK( VXYL3J,  LAY ) +
     &                   ORGFAC  * CBLK( VTOL1J,  LAY ) +
     &                   ORGFAC  * CBLK( VTOL2J,  LAY ) +
     &                   ORGFAC  * CBLK( VTOL3J,  LAY ) +
     &                   ORGFAC  * CBLK( VBNZ1J,  LAY ) +
     &                   ORGFAC  * CBLK( VBNZ2J,  LAY ) +
     &                   ORGFAC  * CBLK( VBNZ3J,  LAY ) +
     &                   ORGFAC  * CBLK( VTRP1J,  LAY ) +
     &                   ORGFAC  * CBLK( VTRP2J,  LAY ) +
     &                   ORGFAC  * CBLK( VISO1J,  LAY ) +
     &                   ORGFAC  * CBLK( VISO2J,  LAY ) +
     &                   ORGFAC  * CBLK( VISO3J,  LAY ) +
     &                   ORGFAC  * CBLK( VSQTJ,   LAY ) +
     &                   ORGFAC  * CBLK( VOLGAJ,  LAY ) +
     &                   ORGFAC  * CBLK( VOLGBJ,  LAY ) +
     &                   ORGFAC  * CBLK( VORGCJ,  LAY ) ) )

         ELSE  ! AE4

            VINAC = MAX( CONMIN,
     &                 ( ORGFAC  * CBLK( VORGPAJ, LAY ) +
     &                   ANTHFAC * CBLK( VP25AJ,  LAY ) +
     &                   ORGFAC  * CBLK( VORGAJ,  LAY ) +
     &                   ORGFAC  * CBLK( VORGBAJ, LAY ) ) )

         END IF  ! AE_VRSN

C***  soot like

         VSTAC =  MAX( CONMIN,
     &                 ANTHFAC * CBLK( VECJ, LAY ) )

C*** water

         VH2AC = MAX( CONMIN,
     &                H2OFAC * CBLK( VH2OAJ, LAY ) )

         VSUMAC = VWSAC + VINAC + VSTAC + VH2AC

C***coarse mode
C***  Water soluble

         VWSCO = MAX( CONMIN,
     &              ( SO4FAC * CBLK( VSO4K, LAY ) +
     &                NH4FAC * CBLK( VNH4K, LAY ) +
     &                NO3FAC * CBLK( VNO3K, LAY ) ) )

C***  insoluble

         VINCO = MAX( CONMIN,
     &              ( SOILFAC * CBLK( VSOILA, LAY ) +
     &                ANTHFAC * CBLK( VANTHA, LAY ) ) )

C***  sea salt

         VSSCO = MAX( CONMIN,
     &              ( SEASFAC * CBLK( VNAK, LAY ) +
     &                SEASFAC * CBLK( VCLK, LAY ) ) )

C***  water

         VH2CO = MAX( CONMIN,
     &                H2OFAC * CBLK( VH2OK, LAY ) )

         VSUMCO = VWSCO + VSSCO + VINCO + VH2CO

C***Set aerosol 3rd moment concentrations [ m**3 / m**3 ]

         CBLK( VAT3,  LAY ) = MAX( CONMIN, VSUMAT )
         CBLK( VAC3,  LAY ) = MAX( CONMIN, VSUMAC )
         CBLK( VCOR3, LAY ) = MAX( CONMIN, VSUMCO )

C***set aerosol 2nd moment concentration [ m**2 / m**3 ]

         CBLK ( VAT2,  LAY ) = CBLK( VSURFAT, LAY ) * PIM1
         CBLK ( VAC2,  LAY ) = CBLK( VSURFAC, LAY ) * PIM1
         CBLK ( VCOR2, LAY ) = CBLK( VSURFCO, LAY ) * PIM1

C***FSB calculate average volume-averaged modal refractive indices:

C***real indices
C***  Aitken Mode

         NR( 1, LAY ) = ( VWSAT * XNREAL_WS  +
     &                    VINAT * XNREAL_IN  +
     &                    VSTAT * XNREAL_ST1 +
     &                    VH2AT * XNREAL_H2 ) / VSUMAT

C***  accumulation mode

         NR( 2, LAY ) = ( VWSAC * XNREAL_WS  +
     &                    VINAC * XNREAL_IN  +
     &                    VSTAC * XNREAL_ST2 +
     &                    VH2AC * XNREAL_H2 ) / VSUMAC

C***  coarse mode

         NR( 3, LAY ) = ( VWSCO * XNREAL_WS +
     &                    VINCO * XNREAL_IN +
     &                    VSSCO * XNREAL_SS +
     &                    VH2CO * XNREAL_H2 ) / VSUMCO

C***imaginary indices
C***  Aitken Mode

         NI( 1, LAY ) = ( VWSAT * XNIMAG_WS  +
     &                    VINAT * XNIMAG_IN  +
     &                    VSTAT * XNIMAG_ST1 +
     &                    VH2AT * XNIMAG_H2 ) / VSUMAT

C***  accumulation mode

         NI( 2, LAY ) = ( VWSAC * XNIMAG_WS  +
     &                    VINAC * XNIMAG_IN  +
     &                    VSTAC * XNIMAG_ST2 +
     &                    VH2AC * XNIMAG_H2 ) / VSUMAC

C***  coarse mode

         NI( 3, LAY ) = ( VWSCO * XNIMAG_WS +
     &                    VINCO * XNIMAG_IN +
     &                    VSSCO * XNIMAG_SS +
     &                    VH2CO * XNIMAG_H2 ) / VSUMCO

C***Calculate geometric standard deviations as follows:
C***
C***  ln^2(Sg) = 1/3*ln(M0) + 2/3*ln(M3) - ln(M2)
C***
C***  NOTES:
C***  1. Equation 10-5a of [Binkowski:1999] and Equation 5a of
C***     [Binkowski&Roselle:2003] contain typographical errors.
C***  2. If the square of the logarithm of the geometric standard
C***     deviation is out of an acceptable range, reset this value and
C***     adjust the second moments to be consistent with this value.
C***     In this manner, M2 is artificially increased when Sg exceeds
C***     the maximum limit.  M2 is artificially decreased when Sg falls
C***     below the minimum limit.

C***Aitken Mode:

         XXM0 = CBLK( VAT0, LAY )
         XXM2 = CBLK( VAT2, LAY )
         XXM3 = CBLK( VAT3, LAY )

C***FSB As noted above ther are typographical errors in the published
C***  documentation for the equation defining the logarithm of the
C***  modal geometric standard deviation.

         XFSUM = ONE3D * LOG( XXM0 ) + TWO3D * LOG( XXM3 )

         LXFM2  = LOG( XXM2 )
         L2SGAT = XFSUM - LXFM2  ! this expression is the correct form

         IF ( L2SGAT .LT. MINL2SG ) THEN

            L2SGAT = MINL2SG
            LXFM2 = XFSUM - L2SGAT
            CBLK( VAT2, LAY ) = EXP( LXFM2 )

         ELSE IF ( L2SGAT .GT. MAXL2SG ) THEN

            L2SGAT = MAXL2SG
            LXFM2 = XFSUM - L2SGAT
            CBLK( VAT2, LAY ) = EXP( LXFM2 )

         END IF   ! test on L2SGAT

         ESAT36  = EXP( 4.5 * L2SGAT )
         XXLSGAT = SQRT( L2SGAT )

C***accumulation mode:

         XXM0 = CBLK( VAC0, LAY )
         XXM2 = CBLK( VAC2, LAY )
         XXM3 = CBLK( VAC3, LAY )

         XFSUM = ONE3D * LOG( XXM0 ) +  TWO3D * LOG( XXM3 )

         LXFM2  = LOG( XXM2 )
         L2SGAC = XFSUM - LXFM2  ! this expression is the correct form

         IF ( L2SGAC .LT. MINL2SG ) THEN

            L2SGAC = MINL2SG
            LXFM2 = XFSUM - L2SGAC
            CBLK( VAC2, LAY ) = EXP( LXFM2 )

         ELSE IF ( L2SGAC .GT. MAXL2SG ) THEN

            L2SGAC =  MAXL2SG
            LXFM2 = XFSUM - L2SGAC
            CBLK( VAC2, LAY ) = EXP( LXFM2 )

         END IF   ! test on L2SGAC

         ESAC36  = EXP( 4.5 * L2SGAC )
         XXLSGAC = SQRT( L2SGAC )

C***coarse mode:

         XXM0 = CBLK( VCOR0, LAY )
         XXM2 = CBLK( VCOR2, LAY )
         XXM3 = CBLK( VCOR3, LAY )

         XFSUM = ONE3D * LOG( XXM0 ) +  TWO3D * LOG( XXM3 )

         LXFM2  = LOG( XXM2 )
         L2SGCO = XFSUM - LXFM2  ! this expression is the correct form

         IF ( L2SGCO .LT. MINL2SG ) THEN

            L2SGCO = MINL2SG
            LXFM2 = XFSUM - L2SGCO
            CBLK( VCOR2, LAY ) = EXP( LXFM2 )

         ELSE IF ( L2SGCO .GT. MAXL2SG ) THEN

            L2SGCO =  MAXL2SG
            LXFM2 = XFSUM - L2SGCO
            CBLK( VCOR2, LAY ) = EXP( LXFM2 )

         END IF   ! test on L2SGCO

         ESC36  = EXP( 4.5 * L2SGCO )
         XXLSGCO = SQRT( L2SGCO )

C***Calculate geometric mean diameters [ m ] using Equation 10-5b of
C***  [Binkowski:1999] and Equation 5b of [Binkowski&Roselle:2003]
C***  Dg^3 = M3 / ( M0 * exp(4.5*ln^2(Sg)) )

         DGATK = MAX( DGMIN, ( CBLK( VAT3, LAY ) /
     &              ( CBLK( VAT0, LAY ) * ESAT36 ) ) ** ONE3 )

         DGACC = MAX( DGMIN, ( CBLK( VAC3, LAY ) /
     &              ( CBLK( VAC0, LAY ) * ESAC36 ) ) ** ONE3 )

         DGCOR = MAX( DGMIN, ( CBLK( VCOR3, LAY ) /
     &              ( CBLK( VCOR0, LAY ) * ESC36 ) ) ** ONE3 )

C***Update CBLK values of Dg and Sg

         CBLK( VLSGAT, LAY ) = XXLSGAT
         CBLK( VLSGAC, LAY ) = XXLSGAC
         CBLK( VLSGCO, LAY ) = XXLSGCO
         CBLK( VDGAT,  LAY ) = DGATK
         CBLK( VDGAC,  LAY ) = DGACC
         CBLK( VDGCO,  LAY ) = DGCOR
         CBLK( VOLAT,  LAY ) = PI6 * CBLK( VAT3,  LAY )
         CBLK( VOLAC,  LAY ) = PI6 * CBLK( VAC3,  LAY )
         CBLK( VOLCO,  LAY ) = PI6 * CBLK( VCOR3, LAY )

      END DO   ! loop over layers

      RETURN

      END SUBROUTINE GETNEWPAR2

      END MODULE AERO_DESC
