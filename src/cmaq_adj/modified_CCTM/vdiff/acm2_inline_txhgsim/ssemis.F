
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhgsim/ssemis.F,v 1.1.1.1 2010/06/14 16:03:09 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ///////////////////////////////////////////////////////////////////////////
C  SUBROUTINE SSEMIS calculates the sea-salt emission rates in a grid cell
C   given the fractional grid-cell area covered by open ocean and surf zone
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  NONE
C
C  REVISION HISTORY:
C
C   13 Jun 05 P.Bhave:  first version for public release
C   11 Apr 08 J.Kelly:  added code for (1) emission of coarse surface area,
C                       (2) emission of coarse water, (3) enhanced surf zone
C                       emissions, and (4) variable coarse std. deviation
C
C  REFERENCES:
C    Alduchov and Eskridge, "Improved Magnus Form Approximations of
C          Saturation Vapor Pressure," Jour. of Applied Meteorology, vol. 35,
C          pp 601-609, April, 1996.
C
C    Binkowski F.S., S.J. Roselle. Models-3 Community Multiscale Air Quality 
C          (CMAQ) model aerosol component 1: Model Description.  J. Geophys. 
C          Res., Vol 108, No D6, 4183, doi:10.1029/2001JD001409 (2003).
C
C    Seinfeld, J.H., S.N. Pandis.  "Atmospheric Chemistry and Physics - from
C          air pollution to climate change"  John Wiley & Sons, Inc. (1998).
C
C    Zhang, K.M., E.M. Knipping, A.S. Wexler, P.V. Bhave, G.S. Tonnesen
C          "Size distribution of sea-salt emissions as a function of relative
C          humidity"  Atmospheric Environment, 39, 3373-3379 (2005).
C    Clarke, A.D., S.R. Owens, and J. Zhou "An ultrafine sea-salt flux from
C          breaking waves: Implications for cloud condensation nuclei in the
C          remote marine atmosphere" J. Geophys. Res. (2006)


      SUBROUTINE SSEMIS ( OFRAC, SFRAC, BLKPRS, BLKTA, BLKQV, U10, 
     &                    RLAY1HT, SSOUTM, SSOUTN, SSOUTS )

      USE AERO_EMIS   ! inherits NSSSPC, NSSMOD, RHO*, and indices for SSOUTM

      IMPLICIT NONE

C *** ARGUMENTS

      REAL OFRAC                    ! fractional seawater cover
      REAL SFRAC                    ! fractional surf-zone cover
      REAL BLKPRS                   ! atmospheric pressure [Pa]
      REAL BLKTA                    ! air temperature [K]
      REAL BLKQV                    ! H2O mass mixing ratio [ kg/kg ]
      REAL U10                      ! wind speed at 10m [m/s]
      REAL RLAY1HT                  ! reciprocal of layer-1 height [1/m]
      REAL SSOUTM( NSSSPC,NSSMOD )  ! mass emission rates [ug/m3/s]
      REAL SSOUTN( NSSMOD )         ! number emission rates [1/m3/s]
      REAL SSOUTS( NSSMOD )         ! surface-area emission rates 
                                    !  [m2/m3/s] (omit coarse mode)
         
C *** LOCAL VARIABLES

C     Lognormal parameters fit to the open-ocean flux function
      REAL DGN1( 136 )        ! geometric mean of accum-mode flux [um]
      REAL DGN2( 136 )        ! geometric mean of coarse-mode flux [um]
      REAL SIG1( 136 )        ! geom std deviation of accum-mode flux
      REAL SIG2( 136 )        ! geom std deviation of coarse-mode flux
      REAL FN1 ( 136 )        ! magnitude of accum-mode flux [1/m2/s]
      REAL FN2 ( 136 )        ! magnitude of coarse-mode flux [1/m2/s]
      INTEGER I               ! loop index

C     Open-ocean flux numerically integrated over size range of interest
      REAL VFLUX( 136 )       ! total particle volume [m3/m2/s]

C     Variables for calculating ambient RH as done in the aerosol module
      REAL, PARAMETER :: AA = 610.94          ! coefficients from Alduchov
      REAL, PARAMETER :: BB = 17.625          !  and Eskridge (1996)
      REAL, PARAMETER :: CC = 243.04
      REAL, PARAMETER :: MWWAT  =  18.0153    ! molwt of H2O vapor [ g/mol ]
      REAL, PARAMETER :: MWAIR  =  28.9628    ! molwt of dry air [ g/mol ]
      REAL, PARAMETER :: EPSWATER = MWWAT / MWAIR
      REAL  ESAT                              ! saturation vapor pressure
      REAL  BLKRH                             ! fractional relative humidity

C     Variables for calculating solute mass fraction & size-correction factors
      REAL RH                 ! BLKRH constrained to 45-99% range
      REAL*8 DRH, DRH2, DRH3  ! powers of RH stored for polynomial
      REAL*8 DRH4, DRH5       !  calculations
      REAL*8 C0_COEFF( 6 )    ! polycoeffs from Eq8 of Zhang et al (2005)
      REAL C0                 ! (Diam @ formation) / (Diam @ ambient RH)
      REAL C80                ! (Diam @ 80% RH) / (Diam @ ambient RH)    
      REAL*8 X_COEFF( 6 )     ! polycoeffs from Eq1 of Zhang et al (2005)
      REAL XSW                ! fractional solute weight fraction

C     Sea-salt speciation factors based on dry mass
      REAL, PARAMETER ::  SPCFAC_NA  = 0.3856
      REAL, PARAMETER ::  SPCFAC_CL  = 0.5389
      REAL, PARAMETER ::  SPCFAC_SO4 = 0.0755
      
C     Density of dry sea salt [kg/m3] using volume additivity
      REAL, PARAMETER ::  RHOSSDRY = 1.0  / ( (SPCFAC_NA+SPCFAC_CL) / RHOSEAS
     &                                    +    SPCFAC_SO4           / RHOSO4 )
      REAL RHOSW              ! sea-salt particle density at ambient RH [g/m3]

C     Maximum and minimum diameters for emission-flux integrations [m]
      REAL, PARAMETER :: DPMAX = 20.0E-6 ! upper limit used for the
                                         ! numerical integrations
      REAL DPMINO             ! minimum diameter for open-ocean flux
      REAL DPMAXO             ! maximum diameter for open-ocean flux
      REAL DPMINS             ! minimum diameter for surf-zone flux
      REAL DPMAXS             ! maximum diameter for surf-zone flux
   
C     Wind-speed-dependent scaling factors for emissions flux functions
      REAL WDSCALO            ! Gong-Monahan open-ocean flux
      REAL WDSCALS            ! deLeeuw surf-zone flux

C     Fraction of whitecap coverage
      REAL WCAP

C     Lognormal parameters and numerically-integrated fluxes at ambient RH
      INTEGER IPOS            ! position of ambient RH in data arrays
      REAL DGN1RH             ! geometric mean of accum-mode flux [m]
      REAL DGN2RH             ! geometric mean of coarse-mode flux [m]
      REAL SIG1RH             ! geom std deviation of accum-mode flux
      REAL SIG2RH             ! geom std deviation of coarse-mode flux
      REAL FN1RH              ! magnitude of accum-mode flux [1/m2/s]
      REAL FN2RH              ! magnitude of coarse-mode flux [1/m2/s]
      REAL VFLXRH             ! total particle volume [m3/m2/s]

C     Alternate forms of geometric std devs stored for efficiency
      REAL LNSG1, LNSG2       ! natural log of std dev
      REAL RDIV1, RDIV2       ! reciprocal of sqrt(2)*lnsg

C     Geometric mean diameters by volume/mass [m]
      REAL DGV1RH             ! accumulation mode
      REAL DGV2RH             ! coarse mode

C     Variables for converting 3rd moment emission rates to M0 and M2 rates
      REAL FACTNUMAC          ! conversion for accumulation mode M0
      REAL FACTNUMC           ! conversion for coarse mode M0
      REAL FACTM2AC           ! conversion for accumulation mode M2
      REAL FACTM2CO           ! conversion for coarse mode M2 

C     Variables for calculating fraction of mass emissions assigned to each mode
      REAL FFIN1              ! finite integral for accumulation mode
      REAL FFIN2              ! finite integral for coarse mode
      REAL FRACM1             ! mass fraction assigned to accumulation mode
      REAL FRACM2             ! mass fraction assigned to coarse mode

C     Mode-specific sea-salt mass emission rates [g/m3/s]
      REAL ESEASACC           ! accumulation mode
      REAL ESEASCOR           ! coarse mode

C     3rd moment emission rates [m3/m3/s]
      REAL M3ACC              ! accumulation mode, 3rd moment
      REAL M3COR              ! coarse mode, 3rd moment
      
C     Variables for calculating dry surface-area emission rate for accum mode
      REAL WET_M3AC, DRY_M3AC ! wet & dry 3rd moment emission rates [m3/m3/s] 
      REAL WET_M2AC, DRY_M2AC ! wet & dry 2nd moment emission rates [m3/m3/s] 
      REAL WET_M3CO, DRY_M3CO ! wet & dry 3rd moment emission rates [m3/m3/s] 
      REAL WET_M2CO, DRY_M2CO ! wet & dry 2nd moment emission rates [m3/m3/s] 


C     Mathematical constants; statement function for ERF
      REAL( 8 ), PARAMETER :: PI = 3.14159265358979324
      REAL, PARAMETER :: F6DPI = 6.0 / PI
      REAL, PARAMETER :: ONE3 = 1.0 / 3.0
      REAL, PARAMETER :: TWO3 = 2.0 / 3.0
      REAL ERF                ! complementary error function
      REAL XX                 ! dummy argument for ERF statement function
      ERF(XX) = SIGN( 1.0, XX) * SQRT(1.0 - EXP( -4.0 * XX * XX / PI ))

C -------------------------- BEGIN DATA STATEMENTS --------------------------

C *** RH-dependent values calculated using MATLAB

C     Geometric mean diameter of accumulation mode [um]
      DATA (DGN1(I),I=1,136)/
     &  0.2651, 0.2677, 0.2696, 0.2713, 0.2736, 0.2752, 0.2768, 0.2784, 
     &  0.2806, 0.2823, 0.2840, 0.2858, 0.2883, 0.2903, 0.2924, 0.2951, 
     &  0.2974, 0.2998, 0.3028, 0.3054, 0.3080, 0.3113, 0.3140, 0.3174, 
     &  0.3203, 0.3232, 0.3267, 0.3297, 0.3332, 0.3363, 0.3394, 0.3430, 
     &  0.3462, 0.3499, 0.3533, 0.3573, 0.3611, 0.3652, 0.3700, 0.3752, 
     &  0.3806, 0.3871, 0.3945, 0.4030, 0.4129, 0.4246, 0.4259, 0.4271, 
     &  0.4284, 0.4298, 0.4311, 0.4328, 0.4342, 0.4356, 0.4371, 0.4385, 
     &  0.4401, 0.4416, 0.4434, 0.4451, 0.4467, 0.4484, 0.4501, 0.4518, 
     &  0.4538, 0.4556, 0.4575, 0.4594, 0.4613, 0.4633, 0.4655, 0.4676, 
     &  0.4697, 0.4718, 0.4740, 0.4764, 0.4787, 0.4811, 0.4835, 0.4859, 
     &  0.4886, 0.4912, 0.4938, 0.4965, 0.4994, 0.5022, 0.5051, 0.5081, 
     &  0.5112, 0.5143, 0.5175, 0.5208, 0.5242, 0.5276, 0.5312, 0.5348, 
     &  0.5386, 0.5424, 0.5463, 0.5503, 0.5544, 0.5587, 0.5631, 0.5675, 
     &  0.5722, 0.5769, 0.5818, 0.5867, 0.5919, 0.5972, 0.6026, 0.6083, 
     &  0.6140, 0.6201, 0.6263, 0.6325, 0.6392, 0.6458, 0.6529, 0.6600, 
     &  0.6676, 0.6751, 0.6833, 0.6914, 0.7001, 0.7088, 0.7177, 0.7276, 
     &  0.7372, 0.7478, 0.7582, 0.7690, 0.7810, 0.7928, 0.8050, 0.8187/ 

C     Geometric mean diameter of coarse mode [um]
      DATA (DGN2(I),I=1,136)/
     &  1.7618, 1.8130, 1.8257, 1.8375, 1.8861, 1.8973, 1.9084, 1.9194, 
     &  1.9678, 1.9796, 1.9919, 2.0046, 2.0548, 2.0690, 2.0840, 2.1359, 
     &  2.1525, 2.1698, 2.2236, 2.2423, 2.2617, 2.3171, 2.3378, 2.3936, 
     &  2.4154, 2.4376, 2.4940, 2.5169, 2.5732, 2.5969, 2.6209, 2.6776, 
     &  2.7027, 2.7599, 2.7868, 2.8453, 2.8754, 2.9076, 2.9723, 3.0398, 
     &  3.0835, 3.1604, 3.2439, 3.3357, 3.4386, 3.5555, 3.5660, 3.5767, 
     &  3.5876, 3.5986, 3.6099, 3.6421, 3.6539, 3.6659, 3.6782, 3.6906, 
     &  3.7034, 3.7163, 3.7487, 3.7623, 3.7761, 3.7903, 3.8047, 3.8194, 
     &  3.8519, 3.8673, 3.8830, 3.8990, 3.9154, 3.9321, 3.9651, 3.9826, 
     &  4.0005, 4.0188, 4.0374, 4.0710, 4.0906, 4.1106, 4.1311, 4.1520, 
     &  4.1860, 4.2080, 4.2305, 4.2535, 4.2878, 4.3120, 4.3367, 4.3621, 
     &  4.3969, 4.4236, 4.4510, 4.4790, 4.5146, 4.5442, 4.5745, 4.6104, 
     &  4.6424, 4.6752, 4.7089, 4.7463, 4.7816, 4.8182, 4.8565, 4.8952, 
     &  4.9350, 4.9746, 5.0169, 5.0563, 5.1010, 5.1470, 5.1884, 5.2369, 
     &  5.2789, 5.3308, 5.3844, 5.4289, 5.4856, 5.5313, 5.5922, 5.6394, 
     &  5.7042, 5.7532, 5.8228, 5.8741, 5.9486, 6.0022, 6.0564, 6.1392, 
     &  6.1968, 6.2860, 6.3476, 6.4099, 6.5098, 6.5772, 6.6466, 6.7590/

C     Geometric standard deviation of accumulation mode
      DATA (SIG1(I),I=1,136)/
     &  1.7572, 1.7612, 1.7612, 1.7612, 1.7651, 1.7651, 1.7651, 1.7651, 
     &  1.7689, 1.7689, 1.7689, 1.7689, 1.7727, 1.7727, 1.7727, 1.7764, 
     &  1.7764, 1.7764, 1.7801, 1.7801, 1.7801, 1.7837, 1.7837, 1.7873, 
     &  1.7873, 1.7873, 1.7907, 1.7907, 1.7941, 1.7941, 1.7941, 1.7974, 
     &  1.7974, 1.8006, 1.8006, 1.8037, 1.8037, 1.8037, 1.8066, 1.8095, 
     &  1.8095, 1.8122, 1.8147, 1.8172, 1.8195, 1.8217, 1.8217, 1.8217, 
     &  1.8217, 1.8217, 1.8217, 1.8237, 1.8237, 1.8237, 1.8237, 1.8237, 
     &  1.8237, 1.8237, 1.8255, 1.8255, 1.8255, 1.8255, 1.8255, 1.8255, 
     &  1.8272, 1.8272, 1.8272, 1.8272, 1.8272, 1.8272, 1.8287, 1.8287, 
     &  1.8287, 1.8287, 1.8287, 1.8300, 1.8300, 1.8300, 1.8300, 1.8300, 
     &  1.8312, 1.8312, 1.8312, 1.8312, 1.8322, 1.8322, 1.8322, 1.8322, 
     &  1.8330, 1.8330, 1.8330, 1.8330, 1.8336, 1.8336, 1.8336, 1.8340, 
     &  1.8340, 1.8340, 1.8340, 1.8342, 1.8342, 1.8342, 1.8343, 1.8343, 
     &  1.8343, 1.8342, 1.8342, 1.8339, 1.8339, 1.8339, 1.8334, 1.8334, 
     &  1.8327, 1.8327, 1.8327, 1.8318, 1.8318, 1.8307, 1.8307, 1.8295, 
     &  1.8295, 1.8281, 1.8281, 1.8266, 1.8266, 1.8248, 1.8229, 1.8230, 
     &  1.8209, 1.8209, 1.8188, 1.8164, 1.8164, 1.8140, 1.8115, 1.8115/
           
C     Geometric standard deviation of coarse mode
      DATA (SIG2(I),I=1,136)/
     &  2.1039, 2.0834, 2.0834, 2.0834, 2.0639, 2.0639, 2.0639, 2.0639, 
     &  2.0452, 2.0452, 2.0452, 2.0452, 2.0273, 2.0273, 2.0273, 2.0101, 
     &  2.0101, 2.0101, 1.9937, 1.9937, 1.9937, 1.9780, 1.9780, 1.9631, 
     &  1.9631, 1.9631, 1.9488, 1.9488, 1.9351, 1.9351, 1.9351, 1.9222, 
     &  1.9222, 1.9098, 1.9098, 1.8981, 1.8981, 1.8981, 1.8870, 1.8765, 
     &  1.8765, 1.8666, 1.8573, 1.8486, 1.8405, 1.8329, 1.8329, 1.8329, 
     &  1.8329, 1.8329, 1.8329, 1.8259, 1.8259, 1.8259, 1.8259, 1.8259, 
     &  1.8259, 1.8259, 1.8195, 1.8195, 1.8195, 1.8195, 1.8195, 1.8195, 
     &  1.8137, 1.8137, 1.8137, 1.8137, 1.8137, 1.8137, 1.8085, 1.8085, 
     &  1.8085, 1.8085, 1.8085, 1.8038, 1.8038, 1.8038, 1.8038, 1.8038, 
     &  1.7997, 1.7997, 1.7997, 1.7997, 1.7963, 1.7963, 1.7963, 1.7963, 
     &  1.7934, 1.7934, 1.7934, 1.7934, 1.7913, 1.7913, 1.7913, 1.7897, 
     &  1.7897, 1.7897, 1.7897, 1.7888, 1.7889, 1.7889, 1.7886, 1.7886, 
     &  1.7886, 1.7891, 1.7891, 1.7904, 1.7903, 1.7903, 1.7923, 1.7924, 
     &  1.7952, 1.7952, 1.7952, 1.7988, 1.7989, 1.8034, 1.8034, 1.8089, 
     &  1.8090, 1.8155, 1.8154, 1.8230, 1.8230, 1.8318, 1.8417, 1.8417, 
     &  1.8529, 1.8529, 1.8656, 1.8799, 1.8799, 1.8960, 1.9138, 1.9138/
           
C     Number flux in accumulation mode at 1 m/s [1/m2/s]
      DATA (FN1(I),I=1,136)/
     &  89.326, 89.524, 89.524, 89.524, 89.709, 89.709, 89.709, 89.709, 
     &  89.881, 89.881, 89.881, 89.881, 90.040, 90.040, 90.040, 90.188, 
     &  90.188, 90.188, 90.325, 90.324, 90.324, 90.451, 90.451, 90.569, 
     &  90.569, 90.569, 90.678, 90.678, 90.779, 90.779, 90.779, 90.873, 
     &  90.873, 90.959, 90.959, 91.039, 91.039, 91.039, 91.113, 91.181, 
     &  91.181, 91.243, 91.301, 91.353, 91.401, 91.445, 91.445, 91.445, 
     &  91.445, 91.445, 91.445, 91.484, 91.484, 91.484, 91.484, 91.484, 
     &  91.484, 91.484, 91.520, 91.520, 91.520, 91.520, 91.520, 91.520, 
     &  91.551, 91.551, 91.551, 91.551, 91.551, 91.551, 91.579, 91.579, 
     &  91.579, 91.579, 91.579, 91.604, 91.604, 91.604, 91.604, 91.604, 
     &  91.625, 91.625, 91.625, 91.625, 91.642, 91.642, 91.642, 91.642, 
     &  91.656, 91.656, 91.656, 91.656, 91.667, 91.667, 91.667, 91.675, 
     &  91.675, 91.675, 91.675, 91.679, 91.679, 91.679, 91.680, 91.680, 
     &  91.680, 91.678, 91.678, 91.672, 91.672, 91.672, 91.663, 91.663, 
     &  91.650, 91.650, 91.650, 91.633, 91.633, 91.612, 91.612, 91.587, 
     &  91.587, 91.558, 91.558, 91.524, 91.524, 91.485, 91.441, 91.441, 
     &  91.392, 91.392, 91.337, 91.274, 91.274, 91.205, 91.129, 91.129/
           
C     Number flux in coarse mode at 1 m/s [1/m2/s]
      DATA (FN2(I),I=1,136)/
     &  11.315, 11.110, 11.110, 11.109, 10.919, 10.919, 10.919, 10.919, 
     &  10.740, 10.740, 10.740, 10.740, 10.573, 10.573, 10.573, 10.417, 
     &  10.417, 10.417, 10.270, 10.271, 10.271, 10.134, 10.134, 10.006, 
     &  10.006, 10.006, 9.886, 9.886, 9.774, 9.774, 9.774, 9.669, 
     &  9.669, 9.571, 9.571, 9.480, 9.480, 9.480, 9.395, 9.315, 
     &  9.315, 9.243, 9.174, 9.111, 9.053, 8.999, 8.999, 8.999, 
     &  9.000, 9.000, 9.000, 8.951, 8.951, 8.951, 8.951, 8.951, 
     &  8.951, 8.951, 8.907, 8.907, 8.907, 8.907, 8.907, 8.907, 
     &  8.867, 8.867, 8.867, 8.867, 8.867, 8.867, 8.832, 8.832, 
     &  8.832, 8.832, 8.832, 8.800, 8.800, 8.800, 8.800, 8.800, 
     &  8.773, 8.773, 8.773, 8.773, 8.751, 8.751, 8.751, 8.751, 
     &  8.732, 8.732, 8.732, 8.732, 8.718, 8.718, 8.718, 8.708, 
     &  8.708, 8.708, 8.708, 8.702, 8.703, 8.703, 8.701, 8.701, 
     &  8.701, 8.704, 8.704, 8.712, 8.712, 8.711, 8.724, 8.725, 
     &  8.742, 8.742, 8.742, 8.765, 8.765, 8.794, 8.794, 8.828, 
     &  8.828, 8.869, 8.869, 8.915, 8.915, 8.970, 9.031, 9.031, 
     &  9.101, 9.101, 9.179, 9.268, 9.268, 9.366, 9.476, 9.476/
           
C     Volume flux at 1 m/s numerically integrated over full size range [m3/m2/s]
      DATA (VFLUX(I),I=1,136)/
     &  3.3799e-16, 3.4012e-16, 3.4728e-16, 3.5407e-16, 3.5489e-16,  
     &  3.6124e-16, 3.6758e-16, 3.7402e-16, 3.7473e-16, 3.8152e-16,  
     &  3.8864e-16, 3.9617e-16, 3.9799e-16, 4.0632e-16, 4.1518e-16,  
     &  4.1825e-16, 4.2807e-16, 4.3846e-16, 4.4284e-16, 4.5420e-16,  
     &  4.6612e-16, 4.7169e-16, 4.8449e-16, 4.9073e-16, 5.0429e-16,  
     &  5.1829e-16, 5.2529e-16, 5.3992e-16, 5.4734e-16, 5.6260e-16,  
     &  5.7831e-16, 5.8648e-16, 6.0311e-16, 6.1215e-16, 6.3024e-16,  
     &  6.4087e-16, 6.6140e-16, 6.8386e-16, 6.9948e-16, 7.1775e-16,  
     &  7.4926e-16, 7.7574e-16, 8.0807e-16, 8.4803e-16, 8.9805e-16,  
     &  9.6150e-16, 9.7003e-16, 9.7878e-16, 9.8774e-16, 9.9692e-16,  
     &  1.0063e-15, 1.0026e-15, 1.0124e-15, 1.0224e-15, 1.0327e-15,  
     &  1.0432e-15, 1.0541e-15, 1.0652e-15, 1.0624e-15, 1.0739e-15,  
     &  1.0858e-15, 1.0980e-15, 1.1106e-15, 1.1235e-15, 1.1216e-15,  
     &  1.1351e-15, 1.1490e-15, 1.1633e-15, 1.1780e-15, 1.1931e-15,  
     &  1.1924e-15, 1.2083e-15, 1.2246e-15, 1.2415e-15, 1.2588e-15,  
     &  1.2593e-15, 1.2776e-15, 1.2964e-15, 1.3159e-15, 1.3360e-15,  
     &  1.3380e-15, 1.3592e-15, 1.3811e-15, 1.4038e-15, 1.4071e-15,  
     &  1.4311e-15, 1.4559e-15, 1.4816e-15, 1.4867e-15, 1.5139e-15,  
     &  1.5422e-15, 1.5715e-15, 1.5786e-15, 1.6098e-15, 1.6422e-15,  
     &  1.6509e-15, 1.6855e-15, 1.7215e-15, 1.7590e-15, 1.7706e-15,  
     &  1.8107e-15, 1.8526e-15, 1.8666e-15, 1.9115e-15, 1.9585e-15,  
     &  1.9754e-15, 2.0260e-15, 2.0447e-15, 2.0993e-15, 2.1565e-15,  
     &  2.1790e-15, 2.2411e-15, 2.2660e-15, 2.3335e-15, 2.4045e-15, 
     &  2.4346e-15, 2.5121e-15, 2.5456e-15, 2.6304e-15, 2.6678e-15, 
     &  2.7610e-15, 2.8028e-15, 2.9056e-15, 2.9525e-15, 3.0662e-15, 
     &  3.1191e-15, 3.1734e-15, 3.3053e-15, 3.3669e-15, 3.5144e-15, 
     &  3.5847e-15, 3.6572e-15, 3.8310e-15, 3.9146e-15, 4.0012e-15, 
     &  4.2082e-15/                                               

C *** Polynomial coefficients from Zhang et al. (2005)      
      DATA C0_COEFF  / 28.376, -205.44, 653.37, -1031.7, 803.18, -247.08 /
      DATA X_COEFF   / 3.1657, -19.079,  55.72, -83.998, 63.436, -19.248 /

C ------------------------------- BEGIN CODE -------------------------------
           
C *** Calculate fractional relative humidity in the grid cell
           
      ESAT  = AA * EXP( BB * ( BLKTA - 273.15 ) / ( BLKTA - 273.15 + CC ) )
      BLKRH = BLKPRS * BLKQV / ( ( EPSWATER + BLKQV ) *  ESAT )       
      BLKRH = MIN( 0.99, MAX( 0.005,  BLKRH ) )

C *** Calculate RH-dependent size-correction factors and solute mass fraction
C     using equations from Zhang et al. (2005), which are valid over the
C     45-99% RH range.  Store powers of RH for polynomial calculations.

      RH = BLKRH
      IF ( RH .GT. 0.99 ) THEN
          RH = .99
      ELSE IF ( RH .LT. 0.45 ) THEN
          RH = .45
      ENDIF
      DRH  = RH
      DRH2 = DRH * DRH
      DRH3 = DRH * DRH2
      DRH4 = DRH * DRH3
      DRH5 = DRH * DRH4

C *** Calculate size-correction factors from Equations 8 & 9 of Zhang et al.

      C0 = C0_COEFF( 1 ) + C0_COEFF( 2 ) * DRH  + C0_COEFF( 3 ) * DRH2
     &                   + C0_COEFF( 4 ) * DRH3 + C0_COEFF( 5 ) * DRH4
     &                   + C0_COEFF( 6 ) * DRH5
      C80 = C0 / 1.97

C *** Calculate solute mass fraction using Equation 1 of Zhang et al.

      XSW = X_COEFF( 1 ) + X_COEFF( 2 ) * DRH   + X_COEFF( 3 ) * DRH2
     &                   + X_COEFF( 4 ) * DRH3  + X_COEFF( 5 ) * DRH4
     &                   + X_COEFF( 6 ) * DRH5

C *** Calculate sea-salt-particle density [g/m3] at ambient RH, assuming
C     volume additivity of dry salt plus water

      RHOSW = 1.0E+3 / ( XSW/RHOSSDRY + (1.0-XSW)/RHOH2O )

C *** Set minimum and maximum diameters for integration using the size-
C     correction factors of Zhang et al. (2005)

      DPMINO = 0.12E-6 / C80          ! Gong-Monahan flux function is valid 
      DPMAXO = MIN(DPMAX, 20.E-6/C80) !  from 0.03-5.0 um dry radius; multiply
                                      !  by 4 to get wet diameter at 80% RH
      DPMINS = 1.6E-6 / C0            ! deLeeuw flux function is valid from 
      DPMAXS = MIN(DPMAX, 20.E-6/C0)  !  1.6-20 um diameter at formation
         
C *** Calculate scaling factors to correct the tabulated fluxes for the 10m
C     wind speed in this cell.  Note: tabulated values are based on 1 m/s.

      WDSCALO = MIN(U10, 20.0) ** 3.41        ! Gong-Monahan flux function
      WDSCALS = EXP( 0.23 * MIN(U10, 9.0) )   ! deLeeuw flux function

C *** Calculate whitecap coverage fraction 
      WCAP = 3.84e-6 * WDSCALO ! Eq. 5, Clarke et al. (2006) JGR

C *** Find position in the data arrays that corresponds to ambient RH in this
C     grid cell and fetch the appropriate data values

      IF ( BLKRH .LE. 0.45 ) THEN
         IPOS = 1
      ELSEIF ( BLKRH .LE. 0.90 ) THEN
         IPOS = NINT( 100. * BLKRH - 44.)
      ELSEIF ( BLKRH .LE. 0.99 ) THEN
         IPOS = NINT( 1000. * BLKRH - 854.)
      ELSE
         IPOS = 136
      ENDIF

      DGN1RH = DGN1  ( IPOS ) * 1.0E-6   ! convert to [m]
      DGN2RH = DGN2  ( IPOS ) * 1.0E-6   ! convert to [m]
      SIG1RH = SIG1  ( IPOS )
      SIG2RH = SIG2  ( IPOS ) 
      FN1RH  = FN1   ( IPOS )
      FN2RH  = FN2   ( IPOS )
      VFLXRH = VFLUX ( IPOS ) * (OFRAC + SFRAC / WCAP) 

C *** Save certain functions of the geometric standard deviations for later use

      LNSG1  = LOG(SIG1RH)
      LNSG2  = LOG(SIG2RH)
      RDIV1  = 1.0 / ( SQRT( 2.0 ) * LNSG1 )
      RDIV2  = 1.0 / ( SQRT( 2.0 ) * LNSG2 )

C *** Calculate geometric-mean diameters by volume using Equation 7.52 of
C     Seinfeld & Pandis (1998).

      DGV1RH = DGN1RH * EXP( 3.0 * LNSG1 * LNSG1 )
      DGV2RH = DGN2RH * EXP( 3.0 * LNSG2 * LNSG2 )

C *** Calculate factors for converting 3rd moment emission rates into number
C     and 2nd moment emission rates.  See Equations 7b and 7c of Binkowski 
C     & Roselle (2003)

      FACTNUMAC = EXP( 4.5 * LNSG1 * LNSG1 ) / DGV1RH ** 3
      FACTNUMC  = EXP( 4.5 * LNSG2 * LNSG2 ) / DGV2RH ** 3
      FACTM2AC  = EXP( 0.5 * LNSG1 * LNSG1 ) / DGV1RH
      FACTM2CO  = EXP( 0.5 * LNSG2 * LNSG2 ) / DGV2RH 

C *** Calculate modal volume fluxes [m3/m2/s] by evaluating finite integrals 
C     from DPMIN to DPMAX over each lognormal distribution.  Use resulting
C     values to calculate the fraction of the total number emissions to 
C     assign to each mode.  See Equations 19 and 20 of Uma Shankar's
C     "Integration of Sea-Salt Fluxes" memo.

      FFIN1  = 0.5 * FN1RH * DGV1RH ** 3.0 * EXP( -4.5 * LNSG1 * LNSG1 ) *
     &                       ( ERF( LOG( DPMAXO/DGV1RH ) * RDIV1 )  
     &                       - ERF( LOG( DPMINO/DGV1RH ) * RDIV1 ) )
      FFIN2  = 0.5 * FN2RH * DGV2RH ** 3.0 * EXP( -4.5 * LNSG2 * LNSG2 ) *
     &                       ( ERF( LOG( DPMAXO/DGV2RH ) * RDIV2 )  
     &                       - ERF( LOG( DPMINO/DGV2RH ) * RDIV2 ) )
      FRACM1 = FFIN1 / ( FFIN1 + FFIN2 )
      FRACM2 = 1.0 - FRACM1

C *** Calculate mode-specific mass emission rates [g/m3/s], by multiplying
C     numerically-integrated volume fluxes by the modal volume fractions,
C     scaling for wind speed, dividing by grid-cell height, and multiplying
C     by particle density.  Multiply by chemical speciation factors to
C     obtain speciated mass emissions.

      ESEASACC = VFLXRH * FRACM1 * WDSCALO * RLAY1HT * RHOSW
      ESEASCOR = VFLXRH * FRACM2 * WDSCALO * RLAY1HT * RHOSW

      SSOUTM(KNA,1)  = 0.0
      SSOUTM(KCL,1)  = 0.0
      SSOUTM(KSO4,1) = 0.0
      SSOUTM(KH2O,1) = 0.0
      SSOUTM(KNA,2)  = ESEASACC * XSW * SPCFAC_NA  ! accumulation mode sodium
      SSOUTM(KCL,2)  = ESEASACC * XSW * SPCFAC_CL  ! accumulation mode chloride
      SSOUTM(KSO4,2) = ESEASACC * XSW * SPCFAC_SO4 ! accumulation mode sulfate
      SSOUTM(KH2O,2) = ESEASACC * (1.0 - XSW)      ! accumulation mode water
      SSOUTM(KNA,3)  = ESEASCOR * XSW * SPCFAC_NA  ! coarse mode sodium  
      SSOUTM(KCL,3)  = ESEASCOR * XSW * SPCFAC_CL  ! coarse mode chloride
      SSOUTM(KSO4,3) = ESEASCOR * XSW * SPCFAC_SO4 ! coarse mode sulfate 
      SSOUTM(KH2O,3) = ESEASCOR * (1.0 - XSW)      ! coarse mode water 

C *** Calculate mode-specific 3rd moment emission rates [m3/m3/s]

      M3ACC = ESEASACC * F6DPI / RHOSW
      M3COR = ESEASCOR * F6DPI / RHOSW

C *** Calculate mode-specific emission rates of particle number [1/s]

      SSOUTN(1) = 0.0
      SSOUTN(2) = M3ACC * FACTNUMAC
      SSOUTN(3) = M3COR * FACTNUMC

C *** Calculate mode-specific dry surface area emission rates [m2/m3/s].  
C     Subtract water from 3rd moment to obtain dry 3rd moment emission rate.  
C     Calculate dry 2nd moment while holding the standard deviation constant.
C     Multiply dry 2nd moment by PI to obtain dry surface area emission rate.

      WET_M3AC    = M3ACC 
      WET_M2AC    = M3ACC * FACTM2AC
      DRY_M3AC    = WET_M3AC - F6DPI * SSOUTM(KH2O,2) / 1.0E+06
      DRY_M2AC    = WET_M2AC * ( DRY_M3AC / WET_M3AC ) ** (2.0/3.0)

      WET_M3CO    = M3COR
      WET_M2CO    = M3COR * FACTM2CO
      DRY_M3CO    = WET_M3CO - F6DPI * SSOUTM(KH2O,3) / 1.0E+06
      DRY_M2CO    = WET_M2CO * ( DRY_M3CO / WET_M3CO ) ** (2.0/3.0)

      SSOUTS(1) = 0.0
      SSOUTS(2) = PI * DRY_M2AC
      SSOUTS(3) = PI * DRY_M2CO

      RETURN
         
      END SUBROUTINE SSEMIS
