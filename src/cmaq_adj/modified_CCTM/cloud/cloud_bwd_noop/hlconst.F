
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/cloud/cloud_bwd_noop/hlconst.F,v 1.1 2011/02/21 20:57:23 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      REAL FUNCTION HLCONST ( CNAME, TEMP, EFFECTIVE, HPLUS )

C-----------------------------------------------------------------------
C
C  FUNCTION: return the Henry's law constant for the specified substance
C            at the given temperature
C
C  revision history:
C    who        when           what
C  ---------  --------  -------------------------------------
C  S.Roselle  08/15/97  code written for Models-3
C  J.Gipson   06/18/01  added Henry's Law constants 50-55 for saprc99
C  W.Hutzell  07/03/01  added Henry's Law constants 56-57 for Atrazine
C                       and the daughter products from Atrazine and OH
C                       reactions.
C  J.Gipson.  09/06/02  added Henry's Law constants 59-73   for toxics
C  S.Roselle  11/07/02  added capability for calculating the effective
C                       Henry's law constant and updated coefficients
C                       in Henry's law constant table
C  J.Gipson   08/06/03  added Henry's Law constants 77-79
C  G.Sarwar   11/21/04  added constants for chlorine chemistry (Henry's
C                       law constants 80-85 and dissociation constants
C                       14-16
C  R.Bullock  07/05/05  added Henry's Law constants 86-87 for mercury 
C                       with enthalpy calculated from cited laboratory 
C                       data fit to an Arrhenius equation
C  W.Hutzell  02/14/06  added HLC 88 to 116, dissociation constant for
C                       17 (hydrazine)
C  A.Carlton  09/20/06  updated Henry's Law constants for 1,7,19,20,21,30 
C                       O3, NO3, hexane, octane, nonane, methanol and
C                       isoprene reference
C  S.Roselle  10/10/07  changed pointers to parameters; reformatted
C                       variable declarations
C-----------------------------------------------------------------------

      IMPLICIT NONE

C...........INCLUDES and their descriptions

      INCLUDE SUBST_IODECL              ! I/O definitions and declarations
      INCLUDE SUBST_IOPARMS             ! I/O parameters definitions

C...........PARAMETERS and their descriptions:

      INTEGER, PARAMETER :: MXSPCS = 116     ! Number of substances
      INTEGER, PARAMETER :: MXDSPCS = 17     ! Number of dissociating species

C...pointers for the dissociation constants (array B and D)

      INTEGER, PARAMETER :: LSO2       =  1  ! SO2
      INTEGER, PARAMETER :: LHSO3      =  2  ! HSO3
      INTEGER, PARAMETER :: LHNO2      =  3  ! HNO3
      INTEGER, PARAMETER :: LHNO3      =  4  ! HNO3
      INTEGER, PARAMETER :: LCO2       =  5  ! CO2
      INTEGER, PARAMETER :: LHCO3      =  6  ! HCO3
      INTEGER, PARAMETER :: LH2O2      =  7  ! H2O2
      INTEGER, PARAMETER :: LHCHO      =  8  ! HCHO
      INTEGER, PARAMETER :: LHCOOH     =  9  ! HCOOH
      INTEGER, PARAMETER :: LHO2       = 10  ! HO2
      INTEGER, PARAMETER :: LNH4OH     = 11  ! NH4OH
      INTEGER, PARAMETER :: LH2O       = 12  ! H2O
      INTEGER, PARAMETER :: LATRA      = 13  ! Atrazine
      INTEGER, PARAMETER :: LCL2       = 14  ! CL2
      INTEGER, PARAMETER :: LHOCL      = 15  ! HOCL
      INTEGER, PARAMETER :: LHCL       = 16  ! HCL
      INTEGER, PARAMETER :: LHYDRAZINE = 17  ! Hydrazine

C...........ARGUMENTS and their descriptions

      CHARACTER*(*) CNAME               ! name of substance
      REAL          TEMP                ! temperature (K)
      LOGICAL       EFFECTIVE           ! true=compute the effective henry's law constant
      REAL          HPLUS               ! hydrogen ion concentration (mol/l)

C...........SCRATCH LOCAL VARIABLES and their descriptions:

      CHARACTER( 16 ), SAVE :: PNAME = 'HLCONST'  ! program name
      CHARACTER( 16 ), SAVE :: SUBNAME( MXSPCS )  ! list of substance names

      CHARACTER( 120 ) :: XMSG = ' '    ! exit status message string

      INTEGER       SPC                 ! species index

      REAL          HPLUSI              ! 1 / HPLUS
      REAL          HPLUS2I             ! 1 / HPLUS**2
      REAL          CLMINUS             ! chlorine ion conc [CL-]
      REAL          CLMINUSI            ! 1 / CLMINUS
      REAL          TFAC                ! (298-T)/(T*298)
      REAL          AKEQ1               ! temp var for dissociation constant
      REAL          AKEQ2               ! temp var for dissociation constant
      REAL          OHION               ! OH ion concentration
      REAL          KH                  ! temp var for henry's law constant

C...Henry's law constant data taken mostly from Rolf Sanders' Compilation of
C...  Henry's Law Constants for Inorganic and Organic Species of Potential
C...  Importance in Environment Chemistry 1999

      REAL, SAVE :: A( MXSPCS )         ! Henry's law constants at 298.15K (M/atm)
      REAL, SAVE :: E( MXSPCS )         ! enthalpy (like activation energy) (K)

C...dissociation constant data taken mostly from 6.A.1 of Seinfeld and Pandis
C...  Atmospheric Chemistry and Physics, 1997

      REAL, SAVE :: B( MXDSPCS )        ! dissociation constant at 298.15K (M or M2)
      REAL, SAVE :: D( MXDSPCS )        ! -dH/R (K)

      DATA SUBNAME(  1), A(  1), E(  1) / 'O3              ', 1.14E-02, 2.3E+03 / ! Kosak 1983
      DATA SUBNAME(  2), A(  2), E(  2) / 'HO2             ', 4.0E+03, 5.9E+03 /  ! Hanson et al. 1992
      DATA SUBNAME(  3), A(  3), E(  3) / 'H2O2            ', 8.3E+04, 7.4E+03 /  ! O'Sullivan et al. 1996
      DATA SUBNAME(  4), A(  4), E(  4) / 'NH3             ', 6.1E+01, 4.2E+03 /  ! Clegg and Brimblecombe 1989
      DATA SUBNAME(  5), A(  5), E(  5) / 'NO              ', 1.9E-03, 1.4E+03 /  ! Lide and Frederikse 1995
      DATA SUBNAME(  6), A(  6), E(  6) / 'NO2             ', 1.2E-02, 2.5E+03 /  ! Chameides 1984
      DATA SUBNAME(  7), A(  7), E(  7) / 'NO3             ', 0.6E+00, 0.0E+00 /  ! Rudich, Talukdar et al.1996
      DATA SUBNAME(  8), A(  8), E(  8) / 'N2O5            ', 1.0E+30, 0.0E+00 /  ! "inf" Sander and Crutzen 1996
      DATA SUBNAME(  9), A(  9), E(  9) / 'HNO2            ', 5.0E+01, 4.9E+03 /  ! Becker et al. 1996
      DATA SUBNAME( 10), A( 10), E( 10) / 'HNO3            ', 2.1E+05, 8.7E+03 /  ! Leieveld and Crutzen 1991
      DATA SUBNAME( 11), A( 11), E( 11) / 'HNO4            ', 1.2E+04, 6.9E+03 /  ! Regimbal and Mozurkewich 1997
      DATA SUBNAME( 12), A( 12), E( 12) / 'SO2             ', 1.4E+00, 2.9E+03 /  ! Linde and Frederikse 1995
      DATA SUBNAME( 13), A( 13), E( 13) / 'H2SO4           ', 1.0E+30, 0.0E+00 /  ! infinity
      DATA SUBNAME( 14), A( 14), E( 14) / 'METHANE         ', 1.4E-03, 1.6E+03 /  ! Linde and Frederikse 1995
      DATA SUBNAME( 15), A( 15), E( 15) / 'ETHANE          ', 1.9E-03, 2.3E+03 /  ! Linde and Frederikse 1995
      DATA SUBNAME( 16), A( 16), E( 16) / 'PROPANE         ', 1.5E-03, 2.7E+03 /  ! Linde and Frederikse 1995
      DATA SUBNAME( 17), A( 17), E( 17) / 'BUTANE          ', 1.1E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 18), A( 18), E( 18) / 'PENTANE         ', 8.1E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 19), A( 19), E( 19) / 'HEXANE          ', 0.1E-03, 7.5E+03 /  ! Ashworth, Howe et al 1988
      DATA SUBNAME( 20), A( 20), E( 20) / 'OCTANE          ', 2.9E-03, 7.8E+03 /  ! Hansen et al. 1993
      DATA SUBNAME( 21), A( 21), E( 21) / 'NONANE          ', 2.4E-03, 2.1E+02 /  ! Ashworth, Howe et al 1988
      DATA SUBNAME( 22), A( 22), E( 22) / 'DECANE          ', 1.4E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 23), A( 23), E( 23) / 'ETHENE          ', 4.7E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 24), A( 24), E( 24) / 'PROPENE         ', 4.8E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 25), A( 25), E( 25) / 'ISOPRENE        ', 2.8E-02, 0.0E+00 /  ! Karl, Lindinger et al 2003
      DATA SUBNAME( 26), A( 26), E( 26) / 'ACETYLENE       ', 4.1E-02, 1.8E+03 /  ! Wilhelm et al. 1977
      DATA SUBNAME( 27), A( 27), E( 27) / 'BENZENE         ', 1.6E-01, 4.1E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 28), A( 28), E( 28) / 'TOLUENE         ', 1.5E-01, 4.0E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 29), A( 29), E( 29) / 'O-XYLENE        ', 1.9E-01, 4.0E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 30), A( 30), E( 30) / 'METHANOL        ', 2.2E+02, 5.2E+03 /  ! Snider and Dawson 1985
      DATA SUBNAME( 31), A( 31), E( 31) / 'ETHANOL         ', 1.9E+02, 6.6E+03 /  ! Snider and Dawson 1985
      DATA SUBNAME( 32), A( 32), E( 32) / '2-CRESOL        ', 8.2E+02, 0.0E+00 /  ! Betterton 1992
      DATA SUBNAME( 33), A( 33), E( 33) / '4-CRESOL        ', 1.3E+02, 0.0E+00 /  ! Betterton 1992
      DATA SUBNAME( 34), A( 34), E( 34) / 'METHYLHYDROPEROX', 3.1E+02, 5.2E+03 /  ! O'Sullivan et al. 1996
      DATA SUBNAME( 35), A( 35), E( 35) / 'FORMALDEHYDE    ', 3.2E+03, 6.8E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 36), A( 36), E( 36) / 'ACETALDEHYDE    ', 1.4E+01, 5.6E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 37), A( 37), E( 37) / 'GENERIC_ALDEHYDE', 4.2E+03, 0.0E+00 /  ! Graedel and Goldberg 1983
      DATA SUBNAME( 38), A( 38), E( 38) / 'GLYOXAL         ', 3.6E+05, 0.0E+00 /  ! Zhou and Mopper 1990
      DATA SUBNAME( 39), A( 39), E( 39) / 'ACETONE         ', 3.0E+01, 4.6E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 40), A( 40), E( 40) / 'FORMIC_ACID     ', 8.9E+03, 6.1E+03 /  ! Johnson et al. 1996
      DATA SUBNAME( 41), A( 41), E( 41) / 'ACETIC_ACID     ', 4.1E+03, 6.3E+03 /  ! Johnson et al. 1996
      DATA SUBNAME( 42), A( 42), E( 42) / 'METHYL_GLYOXAL  ', 3.2E+04, 0.0E+00 /  ! Zhou and Mopper 1990
      DATA SUBNAME( 43), A( 43), E( 43) / 'CO              ', 9.9E-04, 1.3E+03 /  ! Linde and Frederikse 1995
      DATA SUBNAME( 44), A( 44), E( 44) / 'CO2             ', 3.6E-02, 2.2E+03 /  ! Zheng et al. 1997
      DATA SUBNAME( 45), A( 45), E( 45) / 'PAN             ', 2.8E+00, 6.5E+03 /  ! Kames et al. 1991
      DATA SUBNAME( 46), A( 46), E( 46) / 'MPAN            ', 1.7E+00, 0.0E+00 /  ! Kames and Schurath 1995
      DATA SUBNAME( 47), A( 47), E( 47) / 'OH              ', 3.0E+01, 4.5E+03 /  ! Hanson et al. 1992
      DATA SUBNAME( 48), A( 48), E( 48) / 'METHYLPEROXY_RAD', 2.0E+03, 6.6E+03 /  ! Lelieveld and Crutzen 1991
      DATA SUBNAME( 49), A( 49), E( 49) / 'PEROXYACETIC_ACI', 8.4E+02, 5.3E+03 /  ! O'Sullivan et al. 1996
      DATA SUBNAME( 50), A( 50), E( 50) / 'PROPANOIC_ACID  ', 5.7E+03, 0.0E+00 /  ! Kahn et al. 1995
      DATA SUBNAME( 51), A( 51), E( 51) / '2-NITROPHENOL   ', 7.0E+01, 4.6E+03 /  ! USEPA 1982
      DATA SUBNAME( 52), A( 52), E( 52) / 'PHENOL          ', 1.9E+03, 7.3E+03 /  ! USEPA 1982
      DATA SUBNAME( 53), A( 53), E( 53) / 'BIACETYL        ', 7.4E+01, 5.7E+03 /  ! Betteron 1991
      DATA SUBNAME( 54), A( 54), E( 54) / 'BENZALDEHYDE    ', 3.9E+01, 4.8E+03 /  ! Staudinger and Roberts 1996
      DATA SUBNAME( 55), A( 55), E( 55) / 'PINENE          ', 4.9E-02, 0.0E+00 /  ! Karl and Lindinger 1997
      DATA SUBNAME( 56), A( 56), E( 56) / 'ATRA            ', 4.1E+05, 6.0E+03 /  ! CIBA Corp (1989) and Scholtz (1999)
      DATA SUBNAME( 57), A( 57), E( 57) / 'DATRA           ', 4.1E+05, 6.0E+03 /  ! assumed same as Atrazine
      DATA SUBNAME( 58), A( 58), E( 58) / 'ADIPIC_ACID     ', 2.0E+08, 0.0E+00 /  ! Saxena and Hildemann (1996)
      DATA SUBNAME( 59), A( 59), E( 59) / 'ACROLEIN        ', 8.2E+00, 0.0E+00 /  ! Meylan and Howard (1991)
      DATA SUBNAME( 60), A( 60), E( 60) / '1,3-BUTADIENE   ', 1.4E-02, 0.0E+00 /  ! Mackay and Shiu (1981)
      DATA SUBNAME( 61), A( 61), E( 61) / 'ACRYLONITRILE   ', 7.3E+00, 0.0E+00 /  ! Meylan and Howard (1991)
      DATA SUBNAME( 62), A( 62), E( 62) / 'CARBONTETRACHLOR', 3.4E-02, 4.2E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 63), A( 63), E( 63) / 'PROPYLENE_DICHLO', 3.4E-01, 4.3E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 64), A( 64), E( 64) / '1,3DICHLORPROPEN', 6.5E-01, 4.2E+03 /  ! Wright et al (1992b)
      DATA SUBNAME( 65), A( 65), E( 65) / '1,1,2,2-CL4ETHAN', 2.4E+00, 3.2E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 66), A( 66), E( 66) / 'CHLOROFORM      ', 2.5E-01, 4.5E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 67), A( 67), E( 67) / '1,2DIBROMOETHANE', 1.5E+00, 3.9E+03 /  ! Ashworth et al (1988)
      DATA SUBNAME( 68), A( 68), E( 68) / '1,2DICHLOROETHAN', 7.3E-01, 4.2E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 69), A( 69), E( 69) / 'METHYLENE_CHLORI', 3.6E-01, 4.1E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 70), A( 70), E( 70) / 'PERCHLOROETHYLEN', 5.9E-02, 4.8E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 71), A( 71), E( 71) / 'TRICHLOROETHENE ', 1.0E-01, 4.6E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 72), A( 72), E( 72) / 'VINYL_CHLORIDE  ', 3.9E-02, 3.1E+03 /  ! Staudinger and Roberts (1996)
      DATA SUBNAME( 73), A( 73), E( 73) / 'ETHYLENE_OXIDE  ', 8.4E+00, 0.0E+00 /  ! CRC
      DATA SUBNAME( 74), A( 74), E( 74) / 'PPN             ', 2.9E+00, 0.0E+00 /  ! Kames and Schurath (1995)
      DATA SUBNAME( 75), A( 75), E( 75) / 'NAPHTHALENE     ', 2.0E+00, 3.6E+03 /  ! USEPA 1982
      DATA SUBNAME( 76), A( 76), E( 76) / 'QUINOLINE       ', 3.7E+03, 5.4E+03 /  ! USEPA 1982
      DATA SUBNAME( 77), A( 77), E( 77) / 'MEK             ', 2.0E+01, 5.0E+03 /  ! Zhou and Mopper 1990
      DATA SUBNAME( 78), A( 78), E( 78) / 'MVK             ', 4.1E+01, 0.0E+00 /  ! Iraci et al. 1998
      DATA SUBNAME( 79), A( 79), E( 79) / 'METHACROLEIN    ', 6.5E+00, 0.0E+00 /  ! Iraci et al. 1998
      DATA SUBNAME( 80), A( 80), E( 80) / 'CL2             ', 8.6E-02, 2.0E+03 /  ! ROLF SANDERS COMPILATION (1999)/KAVANAUGH AND TRUSSELL (1980)
      DATA SUBNAME( 81), A( 81), E( 81) / 'HOCL            ', 6.6E+02, 5.9E+03 /  ! ROLF SANDERS COMPILATION (1999)/HUTHWELKER ET AL (1995)
      DATA SUBNAME( 82), A( 82), E( 82) / 'HCL             ', 1.9E+01, 6.0E+02 /  ! ROLF SANDERS COMPILATION (1999)/DEAN (1992)
      DATA SUBNAME( 83), A( 83), E( 83) / 'FMCL            ', 1.1E+00, 0.0E+00 /  ! EPA SUITE PROGRAM/UNIT CONVERTED TO MATCH THE DEFINITION BY ROLF SANDERS.
      DATA SUBNAME( 84), A( 84), E( 84) / 'ICL1            ', 6.9E+01, 0.0E+00 /  ! EPA SUITE PROGRAM/UNIT CONVERTED TO MATCH THE DEFINITION BY ROLF SANDERS.
      DATA SUBNAME( 85), A( 85), E( 85) / 'ICL2            ', 6.9E+01, 0.0E+00 /  ! EPA SUITE PROGRAM/ASSUMED EQUAL TO THAT OF ICL1
      DATA SUBNAME( 86), A( 86), E( 86) / 'HG              ', 1.11E-01, 4.97E+03 /! Elemental Mercury from Clever et al. (1985)
      DATA SUBNAME( 87), A( 87), E( 87) / 'HGIIGAS         ', 1.41E+06, 5.26E+03 /! Hg(II) gas as mercuric chloride from Lindqvist and Rodhe (1985)
      DATA SUBNAME( 88), A( 88), E( 88) / 'TECDD_2378      ', 5.1E+00, 3.6E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 89), A( 89), E( 89) / 'PECDD_12378     ', 4.6E+00, 3.2E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 90), A( 90), E( 90) / 'HXCDD_123478    ', 8.1E+00, 2.9E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 91), A( 91), E( 91) / 'HXCDD_123678    ', 2.9E+00, 2.8E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 92), A( 92), E( 92) / 'HXCDD_123789    ', 6.5E+00, 2.7E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 93), A( 93), E( 93) / 'HPCDD_1234678   ', 1.2E+01, 2.4E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 94), A( 94), E( 94) / 'OTCDD           ', 9.8E+00, 2.3E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 95), A( 95), E( 95) / 'TECDF_2378      ', 8.5E+01, 3.7E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 96), A( 96), E( 96) / 'PECDF_12378     ', 5.2E+01, 2.9E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 97), A( 97), E( 97) / 'PECDF_23478     ', 1.8E+02, 3.0E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 98), A( 98), E( 98) / 'HXCDF_123478    ', 3.8E+01, 2.4E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME( 99), A( 99), E( 99) / 'HXCDF_123678    ', 9.0E+01, 2.9E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME(100), A(100), E(100) / 'HXCDF_234678    ', 1.0E+02, 2.6E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME(101), A(101), E(101) / 'HXCDF_123789    ', 5.6E+01, 2.6E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME(102), A(102), E(102) / 'HPCDF_1234678   ', 2.8E+01, 1.6E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME(103), A(103), E(103) / 'HPCDF_1234789   ', 8.0E+01, 2.1E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME(104), A(104), E(104) / 'OTCDF           ', 7.6E+01, 2.4E+03 /  ! Paasivirta et al. (1999)
      DATA SUBNAME(105), A(105), E(105) / 'NAPHTHOL        ', 3.60E+03, 0.0E+00 / ! Eabraham et al. (1994)
      DATA SUBNAME(106), A(106), E(106) / '1NITRONAPHTHALEN', 5.68E+02, 0.0E+00 / ! Altschuh et al. (1999)
      DATA SUBNAME(107), A(107), E(107) / '2NITRONAPHTHALEN', 6.42E+02, 0.0E+00 / ! HENRYWIN v3.10 (Meylan and Howard, 1991)
      DATA SUBNAME(108), A(108), E(108) / '14NAPHTHOQUINONE', 5.08E+05, 0.0E+00 / ! HENRYWIN v3.10 (Meylan and Howard, 1991)
      DATA SUBNAME(109), A(109), E(109) / '2,4-TOLUENE_DIIS', 7.25E+00, 0.0E+00 / ! HENRYWIN v3.10 (Meylan and Howard, 1991)
      DATA SUBNAME(110), A(110), E(110) / 'HEXAMETHYLE_DIIS', 2.08E+01, 0.0E+00 / ! HENRYWIN v3.10 (Meylan and Howard, 1991)
      DATA SUBNAME(111), A(111), E(111) / 'HYDRAZINE       ', 1.14E+03, 0.0E+00 / ! Daubert and Danner (1989), and Amoore and Hautala (1983)
      DATA SUBNAME(112), A(112), E(112) / 'MALEIC_ANHYDRIDE', 2.54E+02, 0.0E+00 / ! HENRYWIN v3.10 (Meylan and Howard, 1991)
      DATA SUBNAME(113), A(113), E(113) / 'TRIETHYLAMINE   ', 6.71E+00, 0.0E+00 / ! Yalkowsky and Dannenfelser (1992), and Riddick et al. (1986)
      DATA SUBNAME(114), A(114), E(114) / 'P_DICHLOROBENZEN', 2.38E+00, 0.0E+00 / ! MacKay and Shiu (1981), measured
      DATA SUBNAME(115), A(115), E(115) / 'M-XYLENE        ', 1.43E-01, 3.9E+03 / ! Staudinger and Roberts (2001)
      DATA SUBNAME(116), A(116), E(116) / 'P-XYLENE        ', 1.35E-01, 3.7E+03 / ! Staudinger and Roberts (2001)

      DATA B( LSO2   ), D( LSO2   ) / 1.30E-02,  1.96E+03 /  ! SO2*H2O<=>HSO3+H     : Smith and Martell (1976)
      DATA B( LHSO3  ), D( LHSO3  ) / 6.60E-08,  1.50E+03 /  ! HSO3<=>SO3+H         : Smith and Martell (1976)
      DATA B( LHNO2  ), D( LHNO2  ) / 5.10E-04, -1.26E+03 /  ! HNO2(aq)<=>NO2+H     : Schwartz and White (1981)
      DATA B( LHNO3  ), D( LHNO3  ) / 1.54E+01,  8.70E+03 /  ! HNO3(aq)<=>NO3+H     : Schwartz (1984)
      DATA B( LCO2   ), D( LCO2   ) / 4.30E-07, -1.00E+03 /  ! CO2*H2O<=>HCO3+H     : Smith and Martell (1976)
      DATA B( LHCO3  ), D( LHCO3  ) / 4.68E-11, -1.76E+03 /  ! HCO3<=>CO3+H         : Smith and Martell (1976)
      DATA B( LH2O2  ), D( LH2O2  ) / 2.20E-12, -3.73E+03 /  ! H2O2(aq)<=>HO2+H     : Smith and Martell (1976)
      DATA B( LHCHO  ), D( LHCHO  ) / 2.53E+03,  4.02E+03 /  ! HCHO(aq)<=>H2C(OH)2  : Le Hanaf (1968)
      DATA B( LHCOOH ), D( LHCOOH ) / 1.80E-04, -2.00E+01 /  ! HCOOH(aq)<=>HCOO+H   : Martell and Smith (1977)
      DATA B( LHO2   ), D( LHO2   ) / 3.50E-05,  0.00E+00 /  ! HO2(aq)<=>H+O2       : Perrin (1982)
      DATA B( LNH4OH ), D( LNH4OH ) / 1.70E-05, -4.50E+02 /  ! NH4*OH<=>NH4+OH      : Smith and Martell (1976)
      DATA B( LH2O   ), D( LH2O   ) / 1.00E-14, -6.71E+03 /  ! H2O<=>H+OH           : Smith and Martell (1976)
      DATA B( LATRA  ), D( LATRA  ) / 2.09E-02,  0.00E+00 /  ! C8H14ClN5<=>C8H13ClN5+H  : Weber (1970)
      DATA B( LCL2   ), D( LCL2   ) / 5.01E-04,  0.00E+00 /  ! CL2*H2O <=> HOCL + H + CL : LIN AND PEHKONEN, JGR, 103, D21, 28093-28102, NOVEMBER 20, 1998. ALSO SEE NOTE BELOW
      DATA B( LHOCL  ), D( LHOCL  ) / 3.16E-08,  0.00E+00 /  ! HOCL <=>H + OCL      : LIN AND PEHKONEN, JGR, 103, D21, 28093-28102, NOVEMBER 20, 1998
      DATA B( LHCL   ), D( LHCL   ) / 1.74E+06,  6.90E+03 /  ! HCL <=> H + CL       : Marsh and McElroy (1985)
      DATA B( LHYDRAZINE), D( LHYDRAZINE) / 1.11E-08,  0.00E+00 /  ! HYDRAZINE <=> HYDRAZINE+ + OH-  : Moliner and Street (1989)
!-------------------------------------------------------------------------------
! Note for dissociation constant for equation 14: CL2*H2O <=> HOCL + H + CL
! Need aqueous [CL-] concentration to calculate effective henry's law coefficient
! Used a value of 2.0 mM following Lin and Pehkonen, JGR, 103, D21, 28093-28102, November 20, 1998
!-------------------------------------------------------------------------------

C...........EXTERNAL FUNCTIONS and their descriptions:

      INTEGER, EXTERNAL :: INDEX1   ! array position for string matching
      INTEGER, EXTERNAL :: TRIMLEN  ! string length, excl. trailing blanks

C-----------------------------------------------------------------------
C  begin body of subroutine HLCONST

      SPC = INDEX1( CNAME, MXSPCS, SUBNAME )

C...error if species not found in table

      IF ( SPC .LE. 0 ) THEN
        XMSG = CNAME( 1:TRIMLEN( CNAME ) ) // ' not found in Henry''s '//
     &         ' Law Constant table in routine HLCONST.'
        CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT2 )
      END IF

C...compute the Henry's Law Constant

      TFAC = ( 298.0 - TEMP) / ( 298.0 * TEMP )
      KH = A( SPC ) * EXP( E( SPC ) * TFAC )
      HLCONST = KH

C...compute the effective Henry's law constants

      IF ( EFFECTIVE ) THEN

        IF ( HPLUS .LE. 0.0 ) THEN
          XMSG = 'Negative or Zero [H+] concentration specified ' //
     &           'in HLCONST '
          CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT2 )
        END IF

        HPLUSI = 1.0 / HPLUS
        HPLUS2I = HPLUSI * HPLUSI

C...assign a value for clminus.  use 2.0 mM based on Lin and Pehkonene, 1998, JGR

        CLMINUS   = 2.0E-03                ! chlorine ion conc [CL-]
        CLMINUSI  = 1.0 / CLMINUS          ! 1 / CLMINUS

        CHECK_NAME: SELECT CASE ( CNAME( 1:TRIMLEN( CNAME ) ) )

        CASE ('SO2')            !   SO2H2O <=> HSO3- + H+
                                ! & HSO3- <=> SO3= + H+

          AKEQ1 = B( LSO2 )  * EXP( D( LSO2 )  * TFAC )
          AKEQ2 = B( LHSO3 ) * EXP( D( LHSO3 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI + AKEQ1 * AKEQ2 * HPLUS2I )

        CASE ('HNO2')           ! HNO2(aq) <=> NO2- + H+

          AKEQ1 = B( LHNO2 ) * EXP( D( LHNO2 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('HNO3')           ! HNO3(aq) <=> NO3- + H+

          AKEQ1 = B( LHNO3 ) * EXP( D( LHNO3 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('CO2')            !   CO2H2O <=> HCO3- + H+
                                ! & HCO3- <=> CO3= + H+

          AKEQ1 = B( LCO2 )  * EXP( D( LCO2 )  * TFAC )
          AKEQ2 = B( LHCO3 ) * EXP( D( LHCO3 ) * TFAC )
          HLCONST = KH
     &            * ( 1.0 + AKEQ1 * HPLUSI + AKEQ1 * AKEQ2 * HPLUS2I )

        CASE ('H2O2')           ! H2O2(aq) <=> HO2- + H+

          AKEQ1 = B( LH2O2 ) * EXP( D( LH2O2 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('FORMALDEHYDE')   ! HCHO(aq) <=> H2C(OH)2(aq)

          AKEQ1 = B( LHCHO ) * EXP( D( LHCHO ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 )

        CASE ('FORMIC_ACID')    ! HCOOH(aq) <=> HCOO- + H+

          AKEQ1 = B( LHCOOH ) * EXP( D( LHCOOH ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('HO2')            ! HO2(aq) <=> H+ + O2-

          AKEQ1 = B( LHO2 ) * EXP( D( LHO2 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('NH3')            ! NH4OH <=> NH4+ + OH-

          AKEQ1 = B( LNH4OH ) * EXP( D( LNH4OH ) * TFAC )
          AKEQ2 = B( LH2O ) * EXP( D( LH2O ) * TFAC )
          OHION = AKEQ2 * HPLUSI
          HLCONST = KH * ( 1.0 + AKEQ1 / OHION )

        CASE ('HYDRAZINE')      ! HYDRAZINE <=> HYDRAZINE+ + OH-

          AKEQ1 = B( LHYDRAZINE ) * EXP( D( LHYDRAZINE ) * TFAC )
          AKEQ2 = B( LH2O ) * EXP( D( LH2O ) * TFAC )
          OHION = AKEQ2 * HPLUSI
          HLCONST = KH * ( 1.0 + AKEQ1 / OHION )

        CASE ('ATRA', 'DATRA')  !     ATRA(aq)  <=>  ATRA- + H
                                !  or DATRA(aq) <=> DATRA- + H

          AKEQ1   = B( LATRA ) * EXP( D( LATRA ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ( 'CL2' )          ! CL2*H2O <=> HOCL + H + CL
                                ! HOCL <=>H + OCL

          AKEQ1   = B( LCL2 )  * EXP( D( LCL2 ) * TFAC )
          AKEQ2   = B( LHOCL ) * EXP( D( LHOCL ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI * CLMINUSI
     &            + AKEQ1 * AKEQ2 * HPLUS2I * CLMINUSI )

        CASE ( 'HCL' )          ! HCL <=> H+ + CL-

          AKEQ1   = B( LHCL ) * EXP( D( LHCL ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ( 'HOCL' )         ! HOCL <=> H+ + OCL-

          AKEQ1   = B( LHOCL ) * EXP( D( LHOCL ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        END SELECT CHECK_NAME

      END IF

      RETURN
      END
