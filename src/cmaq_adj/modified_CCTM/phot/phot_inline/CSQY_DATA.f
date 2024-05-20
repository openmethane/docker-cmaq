
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/phot/phot_inline/CSQY_DATA.f,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE CSQY_DATA

C-----------------------------------------------------------------------
C
C  FUNCTION:  MODULE CSQY_DATA contains the absorption cross section
C     and quantum yield data for reference photolysis reactions.
C     Additionally, mapping pointers are included to relate chemical
C     mechanism photolysis reactions to the reference data.
C
C  PRECONDITIONS REQUIRED:
C     None
C
C  REVISION  HISTORY:
C       Date   Who             What
C     -------- ------------    -----------------------------------------
C     12/07    S.Roselle       Adapted data statements from module
C                              PHOTMOD12 and developed mapping/pointers
C                              to generalize code for other mechanisms.
C     05/2008  S.Roselle       added cs/qy data from SAPRC99 and CB05
C
C  NOTES:
C
C FSB The following DATA Statements contain the absorption cross
C     sections, quantum yields, and reference temperatures.  If there is
C     no variation with temperature, the reference temperature is set
C     to 350 [K] as a signal to not interpolate, but to use the value
C     in location (3).
C
C FSB Many values of absorption cross sections, quantum yields, and
C     reference temperatures are taken from Fast-JX at:
C
C        www.ess.uci/~prather/fastj.html
C
C     These values are taken from Version 5.0 (December 2004).  In some
C     cases the cross sections are equivalent cross sections that is the
C     product of cross section and quantum yields, Cs * Qy, then the
C     QY DATA are set to 1.0. In other cases where the quantum yields
C     have been calculated for this module, the quantum yields are
C     assigned specific values. These may be 1.0.
C
C *** Other values have been computed for this module by
C     Dr. Joseph Pinto of EPA using the latest IUPAC data sheets where
C     available. Certain species used in SAPRC99 and not available
C     elsewhere have been assigned to the 7 wavelength bins using the
C     flux weighted procedures of Fast-J
C
C *** The cross sections and quantum yields for ozone are taken from
C     Wild,O. X. Zhu, and M.J. Prather, Fast-J:  Accurate simulation of
C         In-and Below-Cloud Photolysis in tropospheric Chemical Models,
C         J. of Atmos. Chem., 37, 245-282, 2000.
C     They have been updated with the IUPAC and JPL recommendations in
C     version 5.0.  Independent calculations have been made to assure
C     their quality.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C.....PARAMETERS and their descriptions:

      INTEGER, PARAMETER :: NPHOT_REF = 68   ! # ref phot reactions
      INTEGER, PARAMETER :: NTEMP_REF = 3    ! # ref temperatures
      INTEGER, PARAMETER :: NWL_REF = 7      ! # ref wavelengths

      INTEGER, PARAMETER :: IACETONE          = 01   ! pointer to ACETONE
      INTEGER, PARAMETER :: IACETONE_SAPRC99  = 02   ! pointer to ACETONE_SAPRC99
      INTEGER, PARAMETER :: IACROLEIN         = 03   ! pointer to ACROLEIN
      INTEGER, PARAMETER :: IACROLEIN_SAPRC99 = 04   ! pointer to ACROLEIN_SAPRC99
      INTEGER, PARAMETER :: IBACL_ADJ         = 05   ! pointer to BACL_ADJ
      INTEGER, PARAMETER :: IBACL_ADJ_SAPRC99 = 06   ! pointer to BACL_ADJ_SAPRC99
      INTEGER, PARAMETER :: IBZCHO            = 07   ! pointer to BZCHO
      INTEGER, PARAMETER :: IBZCHO_SAPRC99    = 08   ! pointer to BZCHO_SAPRC99
      INTEGER, PARAMETER :: IC2CHO            = 09   ! pointer to C2CHO
      INTEGER, PARAMETER :: IC2CHO_SAPRC99    = 10   ! pointer to C2CHO_SAPRC99
      INTEGER, PARAMETER :: ICCHO_R           = 11   ! pointer to CCHO_R
      INTEGER, PARAMETER :: ICCHO_R_SAPRC99   = 12   ! pointer to CCHO_R_SAPRC99
      INTEGER, PARAMETER :: ICL2_FASTJX       = 13   ! pointer to CL2_FASTJX
      INTEGER, PARAMETER :: ICL2_IUPAC04      = 14   ! pointer to CL2_IUPAC04
      INTEGER, PARAMETER :: ICOOH             = 15   ! pointer to COOH
      INTEGER, PARAMETER :: ICOOH_SAPRC99     = 16   ! pointer to COOH_SAPRC99
      INTEGER, PARAMETER :: IFMCL_IUPAC04     = 17   ! pointer to FMCL_IUPAC04
      INTEGER, PARAMETER :: IGLY_ABS          = 18   ! pointer to GLY_ABS
      INTEGER, PARAMETER :: IGLY_ABS_SAPRC99  = 19   ! pointer to GLY_ABS_SAPRC99
      INTEGER, PARAMETER :: IGLY_R            = 20   ! pointer to GLY_R
      INTEGER, PARAMETER :: IGLY_R_SAPRC99    = 21   ! pointer to GLY_R_SAPRC99
      INTEGER, PARAMETER :: IH2O2             = 22   ! pointer to H2O2
      INTEGER, PARAMETER :: IH2O2_SAPRC99     = 23   ! pointer to H2O2_SAPRC99
      INTEGER, PARAMETER :: IHCHO_M           = 24   ! pointer to HCHO_M
      INTEGER, PARAMETER :: IHCHO_M_SAPRC99   = 25   ! pointer to HCHO_M_SAPRC99
      INTEGER, PARAMETER :: IHCHO_R           = 26   ! pointer to HCHO_R
      INTEGER, PARAMETER :: IHCHO_R_SAPRC99   = 27   ! pointer to HCHO_R_SAPRC99
      INTEGER, PARAMETER :: IHNO3             = 28   ! pointer to HNO3
      INTEGER, PARAMETER :: IHNO3_IUPAC04     = 29   ! pointer to HNO3_IUPAC04
      INTEGER, PARAMETER :: IHNO3_SAPRC99     = 30   ! pointer to HNO3_SAPRC99
      INTEGER, PARAMETER :: IHO2NO2           = 31   ! pointer to HO2NO2
      INTEGER, PARAMETER :: IHO2NO2_IUPAC04   = 32   ! pointer to HO2NO2_IUPAC04
      INTEGER, PARAMETER :: IHO2NO2_SAPRC99   = 33   ! pointer to HO2NO2_SAPRC99
      INTEGER, PARAMETER :: IHOCL_FASTJX      = 34   ! pointer to HOCL_FASTJX
      INTEGER, PARAMETER :: IHOCL_IUPAC04     = 35   ! pointer to HOCL_IUPAC04
      INTEGER, PARAMETER :: IHONO             = 36   ! pointer to HONO
      INTEGER, PARAMETER :: IHONO_IUPAC04     = 37   ! pointer to HONO_IUPAC04
      INTEGER, PARAMETER :: IHONO_NO          = 38   ! pointer to HONO_NO
      INTEGER, PARAMETER :: IHONO_NO2         = 39   ! pointer to HONO_NO2
      INTEGER, PARAMETER :: IHONO_NO2_SAPRC99 = 40   ! pointer to HONO_NO2_SAPRC99
      INTEGER, PARAMETER :: IHONO_NO_SAPRC99  = 41   ! pointer to HONO_NO_SAPRC99
      INTEGER, PARAMETER :: IIC3ONO2          = 42   ! pointer to IC3ONO2
      INTEGER, PARAMETER :: IIC3ONO2_SAPRC99  = 43   ! pointer to IC3ONO2_SAPRC99
      INTEGER, PARAMETER :: IKETONE           = 44   ! pointer to KETONE
      INTEGER, PARAMETER :: IKETONE_SAPRC99   = 45   ! pointer to KETONE_SAPRC99
      INTEGER, PARAMETER :: IMGLY_ABS         = 46   ! pointer to MGLY_ABS
      INTEGER, PARAMETER :: IMGLY_ABS_SAPRC99 = 47   ! pointer to MGLY_ABS_SAPRC99
      INTEGER, PARAMETER :: IMGLY_ADJ         = 48   ! pointer to MGLY_ADJ
      INTEGER, PARAMETER :: IMGLY_ADJ_SAPRC99 = 49   ! pointer to MGLY_ADJ_SAPRC99
      INTEGER, PARAMETER :: IMGLY_IUPAC04     = 50   ! pointer to MGLY_IUPAC04
      INTEGER, PARAMETER :: IN2O5_FASTJX      = 51   ! pointer to N2O5_FASTJX
      INTEGER, PARAMETER :: IN2O5_IUPAC04     = 52   ! pointer to N2O5_IUPAC04
      INTEGER, PARAMETER :: INO2              = 53   ! pointer to NO2
      INTEGER, PARAMETER :: INO2_SAPRC99      = 54   ! pointer to NO2_SAPRC99
      INTEGER, PARAMETER :: INO3NO            = 55   ! pointer to NO3NO
      INTEGER, PARAMETER :: INO3NO_SAPRC99    = 56   ! pointer to NO3NO_SAPRC99
      INTEGER, PARAMETER :: INO3NO2           = 57   ! pointer to NO3NO2
      INTEGER, PARAMETER :: INO3NO2_SAPRC99   = 58   ! pointer to NO3NO2_SAPRC99
      INTEGER, PARAMETER :: INTR_IUPAC04      = 59   ! pointer to NTR_IUPAC04
      INTEGER, PARAMETER :: IO3O1D            = 60   ! pointer to O3O1D
      INTEGER, PARAMETER :: IO3O1D_SAPRC99    = 61   ! pointer to O3O1D_SAPRC99
      INTEGER, PARAMETER :: IO3O3P            = 62   ! pointer to O3O3P
      INTEGER, PARAMETER :: IO3O3P_SAPRC99    = 63   ! pointer to O3O3P_SAPRC99
      INTEGER, PARAMETER :: IO3_O1D_IUPAC04   = 64   ! pointer to O3_O1D_IUPAC04
      INTEGER, PARAMETER :: IO3_O3P_IUPAC04   = 65   ! pointer to O3_O3P_IUPAC04
      INTEGER, PARAMETER :: IPACD_CB05        = 66   ! pointer to PACD_CB05
      INTEGER, PARAMETER :: IPAN_FASTJX       = 67   ! pointer to PAN_FASTJX
      INTEGER, PARAMETER :: IPAN_IUPAC04      = 68   ! pointer to PAN_IUPAC04

      INTEGER            :: IWLR             ! wavelength loop variable
      INTEGER            :: ITTR             ! temperature loop variable

C.....Names of the reference photolysis reactions

      CHARACTER*16, SAVE :: PNAME_REF( NPHOT_REF )

      DATA PNAME_REF( IACETONE          ) / 'ACETONE         ' /
      DATA PNAME_REF( IACETONE_SAPRC99  ) / 'ACETONE_SAPRC99 ' /
      DATA PNAME_REF( IACROLEIN         ) / 'ACROLEIN        ' /
      DATA PNAME_REF( IACROLEIN_SAPRC99 ) / 'ACROLEIN_SAPRC99' /
      DATA PNAME_REF( IBACL_ADJ         ) / 'BACL_ADJ        ' /
      DATA PNAME_REF( IBACL_ADJ_SAPRC99 ) / 'BACL_ADJ_SAPRC99' /
      DATA PNAME_REF( IBZCHO            ) / 'BZCHO           ' /
      DATA PNAME_REF( IBZCHO_SAPRC99    ) / 'BZCHO_SAPRC99   ' /
      DATA PNAME_REF( IC2CHO            ) / 'C2CHO           ' /
      DATA PNAME_REF( IC2CHO_SAPRC99    ) / 'C2CHO_SAPRC99   ' /
      DATA PNAME_REF( ICCHO_R           ) / 'CCHO_R          ' /
      DATA PNAME_REF( ICCHO_R_SAPRC99   ) / 'CCHO_R_SAPRC99  ' /
      DATA PNAME_REF( ICL2_FASTJX       ) / 'CL2_FASTJX      ' /
      DATA PNAME_REF( ICL2_IUPAC04      ) / 'CL2_IUPAC04     ' /
      DATA PNAME_REF( ICOOH             ) / 'COOH            ' /
      DATA PNAME_REF( ICOOH_SAPRC99     ) / 'COOH_SAPRC99    ' /
      DATA PNAME_REF( IFMCL_IUPAC04     ) / 'FMCL_IUPAC04    ' /
      DATA PNAME_REF( IGLY_ABS          ) / 'GLY_ABS         ' /
      DATA PNAME_REF( IGLY_ABS_SAPRC99  ) / 'GLY_ABS_SAPRC99 ' /
      DATA PNAME_REF( IGLY_R            ) / 'GLY_R           ' /
      DATA PNAME_REF( IGLY_R_SAPRC99    ) / 'GLY_R_SAPRC99   ' /
      DATA PNAME_REF( IH2O2             ) / 'H2O2            ' /
      DATA PNAME_REF( IH2O2_SAPRC99     ) / 'H2O2_SAPRC99    ' /
      DATA PNAME_REF( IHCHO_M           ) / 'HCHO_M          ' /
      DATA PNAME_REF( IHCHO_M_SAPRC99   ) / 'HCHO_M_SAPRC99  ' /
      DATA PNAME_REF( IHCHO_R           ) / 'HCHO_R          ' /
      DATA PNAME_REF( IHCHO_R_SAPRC99   ) / 'HCHO_R_SAPRC99  ' /
      DATA PNAME_REF( IHNO3             ) / 'HNO3            ' /
      DATA PNAME_REF( IHNO3_IUPAC04     ) / 'HNO3_IUPAC04    ' /
      DATA PNAME_REF( IHNO3_SAPRC99     ) / 'HNO3_SAPRC99    ' /
      DATA PNAME_REF( IHO2NO2           ) / 'HO2NO2          ' /
      DATA PNAME_REF( IHO2NO2_IUPAC04   ) / 'HO2NO2_IUPAC04  ' /
      DATA PNAME_REF( IHO2NO2_SAPRC99   ) / 'HO2NO2_SAPRC99  ' /
      DATA PNAME_REF( IHOCL_FASTJX      ) / 'HOCL_FASTJX     ' /
      DATA PNAME_REF( IHOCL_IUPAC04     ) / 'HOCL_IUPAC04    ' /
      DATA PNAME_REF( IHONO             ) / 'HONO            ' /
      DATA PNAME_REF( IHONO_IUPAC04     ) / 'HONO_IUPAC04    ' /
      DATA PNAME_REF( IHONO_NO          ) / 'HONO_NO         ' /
      DATA PNAME_REF( IHONO_NO2         ) / 'HONO_NO2        ' /
      DATA PNAME_REF( IHONO_NO2_SAPRC99 ) / 'HONO_NO2_SAPRC99' /
      DATA PNAME_REF( IHONO_NO_SAPRC99  ) / 'HONO_NO_SAPRC99 ' /
      DATA PNAME_REF( IIC3ONO2          ) / 'IC3ONO2         ' /
      DATA PNAME_REF( IIC3ONO2_SAPRC99  ) / 'IC3ONO2_SAPRC99 ' /
      DATA PNAME_REF( IKETONE           ) / 'KETONE          ' /
      DATA PNAME_REF( IKETONE_SAPRC99   ) / 'KETONE_SAPRC99  ' /
      DATA PNAME_REF( IMGLY_ABS         ) / 'MGLY_ABS        ' /
      DATA PNAME_REF( IMGLY_ABS_SAPRC99 ) / 'MGLY_ABS_SAPRC99' /
      DATA PNAME_REF( IMGLY_ADJ         ) / 'MGLY_ADJ        ' /
      DATA PNAME_REF( IMGLY_ADJ_SAPRC99 ) / 'MGLY_ADJ_SAPRC99' /
      DATA PNAME_REF( IMGLY_IUPAC04     ) / 'MGLY_IUPAC04    ' /
      DATA PNAME_REF( IN2O5_FASTJX      ) / 'N2O5_FASTJX     ' /
      DATA PNAME_REF( IN2O5_IUPAC04     ) / 'N2O5_IUPAC04    ' /
      DATA PNAME_REF( INO2              ) / 'NO2             ' /
      DATA PNAME_REF( INO2_SAPRC99      ) / 'NO2_SAPRC99     ' /
      DATA PNAME_REF( INO3NO            ) / 'NO3NO           ' /
      DATA PNAME_REF( INO3NO_SAPRC99    ) / 'NO3NO_SAPRC99   ' /
      DATA PNAME_REF( INO3NO2           ) / 'NO3NO2          ' /
      DATA PNAME_REF( INO3NO2_SAPRC99   ) / 'NO3NO2_SAPRC99  ' /
      DATA PNAME_REF( INTR_IUPAC04      ) / 'NTR_IUPAC04     ' /
      DATA PNAME_REF( IO3O1D            ) / 'O3O1D           ' /
      DATA PNAME_REF( IO3O1D_SAPRC99    ) / 'O3O1D_SAPRC99   ' /
      DATA PNAME_REF( IO3O3P            ) / 'O3O3P           ' /
      DATA PNAME_REF( IO3O3P_SAPRC99    ) / 'O3O3P_SAPRC99   ' /
      DATA PNAME_REF( IO3_O1D_IUPAC04   ) / 'O3_O1D_IUPAC04  ' /
      DATA PNAME_REF( IO3_O3P_IUPAC04   ) / 'O3_O3P_IUPAC04  ' /
      DATA PNAME_REF( IPACD_CB05        ) / 'PACD_CB05       ' /
      DATA PNAME_REF( IPAN_FASTJX       ) / 'PAN_FASTJX      ' /
      DATA PNAME_REF( IPAN_IUPAC04      ) / 'PAN_IUPAC04     ' /

C.....Temperatures, Absorption cross sections, and quantum yields
C...  for the reference photolysis reactions

      REAL, SAVE :: TEMP_REF( NTEMP_REF, NPHOT_REF )    ! temperatures

      REAL, SAVE :: CS_REF( NPHOT_REF, NTEMP_REF, NWL_REF ) ! cross sections
      REAL, SAVE :: QY_REF( NPHOT_REF, NTEMP_REF, NWL_REF ) ! quantum yields

C...Load data for the reference photolysis reactions

C...(NO2) :: NO2 + hv -> NO + O
C...  Cross sections and quantum yields are based upon the
C...  data in JPL-14 ( 02/01/03 ) and IUPAC dat sheet PNOx4
C...  dated 07/16 /2001.

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, INO2 ), ITTR=1,3 ) / 248.0, 273.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( INO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.984226E-19, 0.142188E-18, 0.184895E-18, 0.218586E-18,
     &  0.316859E-18, 0.560187E-18, 0.104973E-18 /
      DATA ( CS_REF( INO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.982100E-19, 0.142782E-18, 0.185154E-18, 0.218100E-18,
     &  0.323124E-18, 0.563060E-18, 0.105718E-18 /
      DATA ( CS_REF( INO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.979800E-19, 0.143400E-18, 0.185400E-18, 0.217600E-18,
     &  0.329800E-18, 0.566100E-18, 0.106500E-18 /

C...  quantum yields

      DATA ( QY_REF( INO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.7470, 0.0070 /
      DATA ( QY_REF( INO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.7550, 0.0080 /
      DATA ( QY_REF( INO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.7640, 0.0090 /

C...(NO2_SAPRC99) :: NO2 + hv = NO + O
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07
C...  for wl>424nm, cross sections are taken from Evaluation No. 15,
C...    JPL Publication 06-2, Table 4-12, T=294K
C...    http://jpldataeval.jpl.nasa.gov/pdf/Jpl15_Sectn4_PhotoChemData.pdf
C...    and quantum yields are set to 0.0 for wl>424nm

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, INO2_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( INO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.083E-19, 1.481E-19, 1.866E-19, 2.252E-19, 
     & 3.334E-19, 5.497E-19, 1.061E-19 /
      DATA ( CS_REF( INO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.083E-19, 1.481E-19, 1.866E-19, 2.252E-19, 
     & 3.334E-19, 5.497E-19, 1.061E-19 /
      DATA ( CS_REF( INO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.083E-19, 1.481E-19, 1.866E-19, 2.252E-19, 
     & 3.334E-19, 5.497E-19, 1.061E-19 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (product calculated on finer
C...    wavelength grid).  The quantum yield values below were then
C...    calculated for the 7 wavelength intervals by dividing the
C...    effective cross sections by the interval average cross sections
C...    (eQY=eCS/CS).

      DATA ( QY_REF( INO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.0000, 1.0000, 0.9989, 0.9907, 0.9901, 0.7713, 0.0038 /
      DATA ( QY_REF( INO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.0000, 1.0000, 0.9989, 0.9906, 0.9901, 0.7713, 0.0038 /
      DATA ( QY_REF( INO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.0000, 1.0000, 0.9989, 0.9906, 0.9901, 0.7713, 0.0038 /

C...(NO3NO) :: NO3 + hv -> NO + O2
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, INO3NO ), ITTR=1,3 ) / 230.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( INO3NO, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.835E-18 /
      DATA ( CS_REF( INO3NO, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.835E-18 /
      DATA ( CS_REF( INO3NO, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.835E-18 /

C...  quantum yields

      DATA ( QY_REF( INO3NO, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1140 /
      DATA ( QY_REF( INO3NO, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1140 /
      DATA ( QY_REF( INO3NO, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1140 /

C...(NO3NO_SAPRC99) :: NO3 + hv = NO + O2 (T=298)
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, INO3NO_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( INO3NO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 4.501E-21, 1.442E-18 /
      DATA ( CS_REF( INO3NO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 4.501E-21, 1.442E-18 /
      DATA ( CS_REF( INO3NO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 4.501E-21, 1.442E-18 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( INO3NO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 0.000E+00, 0.0571 /
      DATA ( QY_REF( INO3NO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 0.000E+00, 0.0571 /
      DATA ( QY_REF( INO3NO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 0.000E+00, 0.0571 /

C...(NO3NO2) :: NO3 + hv -> NO2 + O
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,INO3NO2 ),ITTR=1,3 ) / 230.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF(  INO3NO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.835E-18 /
      DATA ( CS_REF(  INO3NO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.835E-18 /
      DATA ( CS_REF(  INO3NO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.835E-18 /

C...  quantum yields

      DATA ( QY_REF(  INO3NO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.8860 /
      DATA ( QY_REF(  INO3NO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.8860 /
      DATA ( QY_REF(  INO3NO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.8860 /

C...(NO3NO2_SAPRC99) :: NO3 + hv = NO2 + O (T=298)
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, INO3NO2_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( INO3NO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 4.501E-21, 1.442E-18 /
      DATA ( CS_REF( INO3NO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000000E+00, 0.000E+00,
     & 0.000E+00, 4.501E-21, 1.442E-18 /
      DATA ( CS_REF( INO3NO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     & 0.000E+00, 4.501E-21, 1.442E-18 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( INO3NO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.5267 /
      DATA ( QY_REF( INO3NO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.5267 /
      DATA ( QY_REF( INO3NO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.0000, 0.5267 /

C...(O3O3P) :: O3 + hv -> O(3P) + O2
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IO3O3P ), ITTR=1,3 ) / 180.0, 260.0, 300.0 /

C...  absorption cross sections

      DATA ( CS_REF( IO3O3P, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.748000E-18, 0.236500E-18, 0.872200E-19, 0.369400E-19,
     &  0.429500E-20, 0.180400E-22, 0.163000E-20 /
      DATA ( CS_REF( IO3O3P, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.793100E-18, 0.257100E-18, 0.967300E-19, 0.414100E-19,
     &  0.545700E-20, 0.277500E-22, 0.163000E-20 /
      DATA ( CS_REF( IO3O3P, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.830500E-18, 0.277700E-18, 0.107500E-18, 0.472500E-19,
     &  0.678200E-20, 0.482400E-22, 0.163000E-20 /

C...  quantum yields

      DATA ( QY_REF( IO3O3P, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.10000, 0.10740, 0.55370, 0.91379, 0.92261, 1.00000, 1.00000 /
      DATA ( QY_REF( IO3O3P, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.10000, 0.10660, 0.49550, 0.85270, 0.91357, 1.00000, 1.00000 /
      DATA ( QY_REF( IO3O3P, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.10000, 0.10600, 0.43180, 0.76370, 0.90035, 1.00000, 1.00000 /

C...(O3O3P_SAPRC99) :: O3 + HV = O1D + O2
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07
C...  Absorption cross sections from NASA (1999)
C...  Quantum yields derived from O3->O1D quantum yields assuming total
C...  quantum yield is 1, though this is not adequately discussed in the
C...  evaluations.

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IO3O3P_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IO3O3P_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 8.158000E-19, 2.713000E-19, 1.034000E-19, 4.448000E-20,
     & 6.142000E-21, 1.770000E-23, 1.654000E-21 /
      DATA ( CS_REF( IO3O3P_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 8.158000E-19, 2.713000E-19, 1.034000E-19, 4.448000E-20,
     & 6.142000E-21, 1.770000E-23, 1.654000E-21 /
      DATA ( CS_REF( IO3O3P_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 8.158000E-19, 2.713000E-19, 1.034000E-19, 4.448000E-20,
     & 6.142000E-21, 1.770000E-23, 1.654000E-21 /

C...  quantum yields

      DATA ( QY_REF( IO3O3P_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 6.208000E-02, 4.479000E-02, 4.696000E-01, 7.889000E-01,
     & 9.572000E-01, 1.000000E+00, 1.000000E+00 /
      DATA ( QY_REF( IO3O3P_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 6.208000E-02, 4.479000E-02, 4.696000E-01, 7.889000E-01,
     & 9.572000E-01, 1.000000E+00, 1.000000E+00 /
      DATA ( QY_REF( IO3O3P_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 6.208000E-02, 4.479000E-02, 4.696000E-01, 7.889000E-01,
     & 9.572000E-01, 1.000000E+00, 1.000000E+00 /

C...(O3_O3P_IUPAC04) :: O3 + HV = O(3P) + O2
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet POx2, updated 2nd October 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IO3_O3P_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IO3_O3P_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 8.156E-19, 2.597E-19, 1.019E-19, 4.131E-20,
     & 6.153E-21, 4.215E-23, 1.554E-21 /
      DATA ( CS_REF( IO3_O3P_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 8.156E-19, 2.597E-19, 1.019E-19, 4.131E-20,
     & 6.153E-21, 4.215E-23, 1.554E-21 /
      DATA ( CS_REF( IO3_O3P_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 8.156E-19, 2.597E-19, 1.019E-19, 4.131E-20,
     & 6.153E-21, 4.215E-23, 1.554E-21 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( IO3_O3P_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.1000, 0.1051, 0.4567, 0.7695, 0.8984, 0.9243, 1.0000 /
      DATA ( QY_REF( IO3_O3P_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.1000, 0.1051, 0.4567, 0.7695, 0.8984, 0.9243, 1.0000 /
      DATA ( QY_REF( IO3_O3P_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.1000, 0.1051, 0.4567, 0.7695, 0.8984, 0.9243, 1.0000 /

C...(O3O1D) :: O3 + hv -> O(1D) + O2
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IO3O1D ), ITTR=1,3 ) / 180.0, 260.0, 300.0 /

C...  absorption cross sections

      DATA ( CS_REF( IO3O1D, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.7480E-18, 0.2365E-18, 0.8722E-19, 0.3694E-19,
     &  0.4295E-20, 0.1804E-22, 0.1630E-20 /
      DATA ( CS_REF( IO3O1D, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.7931E-18, 0.2571E-18, 0.9673E-19, 0.4141E-19,
     &  0.5457E-20, 0.2775E-22, 0.1630E-20 /
      DATA ( CS_REF( IO3O1D, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.8305E-18, 0.2777E-18, 0.1075E-18, 0.4725E-19,
     &  0.6782E-20, 0.4824E-22, 0.1630E-20 /

C...  quantum yields

      DATA ( QY_REF( IO3O1D, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.90000, 0.89260, 0.44630, 0.08621, 0.07739, 0.00000, 0.00000 /
      DATA ( QY_REF( IO3O1D, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.90000, 0.89340, 0.50450, 0.14730, 0.08643, 0.00000, 0.00000 /
      DATA ( QY_REF( IO3O1D, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.90000, 0.89540, 0.56820, 0.23630, 0.09965, 0.00000, 0.00000 /

C...(O3O1D_SAPRC99) :: O3 + HV = O1D + O2
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07
C...  Absorption cross sections from NASA (1999)
C...  Quantum yields from IUPAC, Supplement VI (1997)
C...  No quantum yield recommendation is given for wl>335.
C...  Assume they decrease linearly to zero at 340 nm.

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IO3O1D_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IO3O1D_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 8.158000E-19, 2.713000E-19, 1.034000E-19, 4.448000E-20,
     & 6.142000E-21, 1.770000E-23, 1.654000E-21 /
      DATA ( CS_REF( IO3O1D_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 8.158000E-19, 2.713000E-19, 1.034000E-19, 4.448000E-20,
     & 6.142000E-21, 1.770000E-23, 1.654000E-21 /
      DATA ( CS_REF( IO3O1D_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 8.158000E-19, 2.713000E-19, 1.034000E-19, 4.448000E-20,
     & 6.142000E-21, 1.770000E-23, 1.654000E-21 /

C...  quantum yields

      DATA ( QY_REF( IO3O1D_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 9.379000E-01, 9.552000E-01, 5.304000E-01, 2.111000E-01,
     & 4.281000E-02, 0.000000E+00, 0.000000E+00 /
      DATA ( QY_REF( IO3O1D_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 9.379000E-01, 9.552000E-01, 5.304000E-01, 2.111000E-01,
     & 4.281000E-02, 0.000000E+00, 0.000000E+00 /
      DATA ( QY_REF( IO3O1D_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 9.379000E-01, 9.552000E-01, 5.304000E-01, 2.111000E-01,
     & 4.281000E-02, 0.000000E+00, 0.000000E+00 /

C...(O3_O1D_IUPAC04) :: O3 + HV = O(1D) + O2
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet POx2, updated 2nd October 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IO3_O1D_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IO3_O1D_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 8.156E-19, 2.597E-19, 1.019E-19, 4.131E-20,
     & 6.153E-21, 4.215E-23, 1.554E-21 /
      DATA ( CS_REF( IO3_O1D_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 8.156E-19, 2.597E-19, 1.019E-19, 4.131E-20,
     & 6.153E-21, 4.215E-23, 1.554E-21 /
      DATA ( CS_REF( IO3_O1D_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 8.156E-19, 2.597E-19, 1.019E-19, 4.131E-20,
     & 6.153E-21, 4.215E-23, 1.554E-21 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( IO3_O1D_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.9000, 0.8949, 0.5430, 0.2306, 0.1017, 0.0758, 0.0000 /
      DATA ( QY_REF( IO3_O1D_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.9000, 0.8949, 0.5430, 0.2306, 0.1017, 0.0758, 0.0000 /
      DATA ( QY_REF( IO3_O1D_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.9000, 0.8949, 0.5430, 0.2306, 0.1017, 0.0758, 0.0000 /

C...(HONO) :: HONO + hv -> HO + NO
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHONO ), ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHONO, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.0000E+00, 0.0000E+00, 0.1200E-19, 0.3469E-19,
     &  0.1090E-18, 0.8645E-19, 0.0000E+00 /
      DATA ( CS_REF( IHONO, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.0000E+00, 0.0000E+00, 0.1200E-19, 0.3469E-19,
     &  0.1090E-18, 0.8645E-19, 0.0000E+00 /
      DATA ( CS_REF( IHONO, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.0000E+00, 0.0000E+00, 0.1200E-19, 0.3469E-19,
     &  0.1090E-18, 0.8645E-19, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHONO, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHONO, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHONO, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HONO_IUPAC04) :: HONO + hv = HO + NO
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet PNOx1_HONO, updated 16th July 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHONO_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHONO_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 4.861E-21, 1.620E-20, 3.144E-20,
     & 9.261E-20, 7.239E-20, 0.000E+00 /
      DATA ( CS_REF( IHONO_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 4.861E-21, 1.620E-20, 3.144E-20,
     & 9.261E-20, 7.239E-20, 0.000E+00 /
      DATA ( CS_REF( IHONO_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 4.861E-21, 1.620E-20, 3.144E-20,
     & 9.261E-20, 7.239E-20, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHONO_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHONO_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHONO_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HONO_NO) :: HONO + hv -> HO + NO
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IHONO_NO),ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHONO_NO, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.0000E+00, 0.0000E+00, 0.1200E-19, 0.3469E-19,
     &  0.1090E-18, 0.8645E-19, 0.0000E+00 /
      DATA ( CS_REF( IHONO_NO, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.0000E+00, 0.0000E+00, 0.1200E-19, 0.3469E-19,
     &  0.1090E-18, 0.8645E-19, 0.0000E+00 /
      DATA ( CS_REF( IHONO_NO, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.0000E+00, 0.0000E+00, 0.1200E-19, 0.3469E-19,
     &  0.1090E-18, 0.8645E-19, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHONO_NO, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHONO_NO, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHONO_NO, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HONO_NO_SAPRC99) :: HONO + hv = HO. + NO
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHONO_NO_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHONO_NO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 1.273000E-20, 3.499000E-20,
     & 1.090000E-19, 8.644000E-20, 0.000000E+00 /
      DATA ( CS_REF( IHONO_NO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 1.273000E-20, 3.499000E-20,
     & 1.090000E-19, 8.644000E-20, 0.000000E+00 /
      DATA ( CS_REF( IHONO_NO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 1.273000E-20, 3.499000E-20,
     & 1.090000E-19, 8.644000E-20, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHONO_NO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 2.883000E-01, 4.696000E-01,
     & 6.486000E-01, 6.842000E-01, 1.000000E+00 /
      DATA ( QY_REF( IHONO_NO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 2.883000E-01, 4.696000E-01,
     & 6.486000E-01, 6.842000E-01, 1.000000E+00 /
      DATA ( QY_REF( IHONO_NO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 2.883000E-01, 4.696000E-01,
     & 6.486000E-01, 6.842000E-01, 1.000000E+00 /

C...(HONO_NO2) :: HONO + hv -> H + NO2
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IHONO_NO2),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHONO_NO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     &  0.000000E+00, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IHONO_NO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     &  0.000000E+00, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IHONO_NO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     &  0.000000E+00, 0.000000E+00, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHONO_NO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /
      DATA ( QY_REF( IHONO_NO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /
      DATA ( QY_REF( IHONO_NO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /

C...(HONO_NO2_SAPRC99) :: HONO + hv = H. + NO2
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHONO_NO2_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHONO_NO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 1.273000E-20, 3.499000E-20,
     & 1.090000E-19, 3.903000E-20, 0.000000E+00 /
      DATA ( CS_REF( IHONO_NO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 1.273000E-20, 3.499000E-20,
     & 1.090000E-19, 3.903000E-20, 0.000000E+00 /
      DATA ( CS_REF( IHONO_NO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 1.273000E-20, 3.499000E-20,
     & 1.090000E-19, 3.903000E-20, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHONO_NO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 4.076000E-01, 5.304000E-01,
     & 3.514000E-01, 2.636000E-02, 0.000000E+00 /
      DATA ( QY_REF( IHONO_NO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 4.076000E-01, 5.304000E-01,
     & 3.514000E-01, 2.636000E-02, 0.000000E+00 /
      DATA ( QY_REF( IHONO_NO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000000E+00, 0.000000E+00, 4.076000E-01, 5.304000E-01,
     & 3.514000E-01, 2.636000E-02, 0.000000E+00 /

C...(HNO3) :: HNO3 + hv -> OH + NO2
C...  fzb (could be updated to latest fast-j)

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHNO3 ), ITTR=1,3) / 200.0, 300.0, 300.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHNO3, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.3371E-20, 0.1377E-20, 0.5451E-21, 0.2102E-21,
     &  0.2154E-22, 0.8768E-25, 0.0000E+00 /
      DATA ( CS_REF( IHNO3, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.4354E-20, 0.1923E-20, 0.8314E-21, 0.3589E-21,
     &  0.4764E-22, 0.2667E-24, 0.0000E+00 /
      DATA ( CS_REF( IHNO3, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.4354E-20, 0.1923E-20, 0.8314E-21, 0.3589E-21,
     &  0.4764E-22, 0.2667E-24, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHNO3, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHNO3, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHNO3, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HNO3_IUPAC04) :: HONO2 + hv = OH + NO2
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet PNOx2_HONO2, updated 16th July 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IHNO3_IUPAC04 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHNO3_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 4.328E-21, 1.943E-21, 8.356E-22, 3.627E-22,
     & 4.793E-23, 3.829E-25, 0.000E+00 /
      DATA ( CS_REF( IHNO3_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 4.328E-21, 1.943E-21, 8.356E-22, 3.627E-22,
     & 4.793E-23, 3.829E-25, 0.000E+00 /
      DATA ( CS_REF( IHNO3_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 4.328E-21, 1.943E-21, 8.356E-22, 3.627E-22,
     & 4.793E-23, 3.829E-25, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHNO3_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHNO3_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHNO3_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HNO3_SAPRC99) :: HNO3 + hv = products
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHNO3_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHNO3_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 4.328E-21, 1.943E-21, 8.356E-22, 3.627E-22,
     & 4.793E-23, 3.829E-25, 0.000E+00 /
      DATA ( CS_REF( IHNO3_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 4.328E-21, 1.943E-21, 8.356E-22, 3.627E-22,
     & 4.793E-23, 3.829E-25, 0.000E+00 /
      DATA ( CS_REF( IHNO3_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 4.328E-21, 1.943E-21, 8.356E-22, 3.627E-22,
     & 4.793E-23, 3.829E-25, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHNO3_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHNO3_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHNO3_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HO2NO2) :: HO2NO2 + hv -> products
C...  fzb (could be updated to latest fast-j)

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHO2NO2 ),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHO2NO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.2580E-19, 0.1102E-19, 0.5222E-20, 0.2794E-20,
     &  0.2301E-21, 0.0000E+00, 0.4739E-22 /
      DATA ( CS_REF( IHO2NO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.2580E-19, 0.1102E-19, 0.5222E-20, 0.2794E-20,
     &  0.2301E-21, 0.0000E+00, 0.4739E-22 /
      DATA ( CS_REF( IHO2NO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.2580E-19, 0.1102E-19, 0.5222E-20, 0.2794E-20,
     &  0.2301E-21, 0.0000E+00, 0.4739E-22 /

C...  quantum yields

      DATA ( QY_REF( IHO2NO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHO2NO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHO2NO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HO2NO2_IUPAC04) :: HOONO2 + hv = products
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet PNOx3_HO2NO2, updated 16th July 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHO2NO2_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHO2NO2_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 2.569E-20, 1.074E-20, 5.484E-21, 3.446E-21,
     & 6.359E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IHO2NO2_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 2.569E-20, 1.074E-20, 5.484E-21, 3.446E-21,
     & 6.359E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IHO2NO2_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 2.569E-20, 1.074E-20, 5.484E-21, 3.446E-21,
     & 6.359E-22, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHO2NO2_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHO2NO2_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHO2NO2_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HO2NO2_SAPRC99) :: HO2NO2 + hv = products
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHO2NO2_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHO2NO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 2.569E-20, 1.074E-20, 5.484E-21, 3.446E-21,
     & 6.359E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IHO2NO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 2.569E-20, 1.074E-20, 5.484E-21, 3.446E-21,
     & 6.359E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IHO2NO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 2.569E-20, 1.074E-20, 5.484E-21, 3.446E-21,
     & 6.359E-22, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHO2NO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHO2NO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHO2NO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(H2O2) :: H2O2 + hv -> 2OH
C...  fzb (could be updated to latest fast-j)

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IH2O2 ), ITTR=1,3 ) / 200.0, 300.0, 300.0 /

C...  absorption cross sections

      DATA ( CS_REF( IH2O2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.832000E-20, 0.500800E-20, 0.321500E-20, 0.211600E-20,
     &  0.801000E-21, 0.208800E-22, 0.000000E+00 /
      DATA ( CS_REF( IH2O2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.925800E-20, 0.573500E-20, 0.379700E-20, 0.258500E-20,
     &  0.104900E-20, 0.269900E-22, 0.000000E+00 /
      DATA ( CS_REF( IH2O2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.925800E-20, 0.573500E-20, 0.379700E-20, 0.258500E-20,
     &  0.104900E-20, 0.269900E-22, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IH2O2, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IH2O2, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IH2O2, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(H2O2_SAPRC99) :: H2O2 + hv = 2 OH
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IH2O2_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IH2O2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 9.309E-21, 5.753E-21, 3.902E-21, 2.715E-21,
     & 1.139E-21, 3.565E-23, 0.000E+00 /
      DATA ( CS_REF( IH2O2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 9.309E-21, 5.753E-21, 3.902E-21, 2.715E-21,
     & 1.139E-21, 3.565E-23, 0.000E+00 /
      DATA ( CS_REF( IH2O2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 9.309E-21, 5.753E-21, 3.902E-21, 2.715E-21,
     & 1.139E-21, 3.565E-23, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IH2O2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IH2O2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IH2O2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HCHO_R) :: HCHO + hv -> HCO + H
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IHCHO_R ),ITTR=1,3 ) / 223.0, 273.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHCHO_R, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.208974E-19, 0.235101E-19, 0.110142E-19, 0.178280E-19,
     &  0.612820E-20, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IHCHO_R, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.208731E-19, 0.235059E-19, 0.109941E-19, 0.178646E-19,
     &  0.613506E-20, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IHCHO_R, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.208609E-19, 0.235038E-19, 0.109841E-19, 0.178829E-19,
     &  0.613849E-20, 0.000000E+00, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHCHO_R, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.6700, 0.6990, 0.7180, 0.7110, 0.5650, 0.0010, 0.0000 /
      DATA ( QY_REF( IHCHO_R, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.6700, 0.6990, 0.7180, 0.7110, 0.5650, 0.0010, 0.0000 /
      DATA ( QY_REF( IHCHO_R, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.6700, 0.6990, 0.7180, 0.7110, 0.5650, 0.0010, 0.0000 /

C...(HCHO_R_SAPRC99) :: HCHO + hv = HCO + H
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IHCHO_R_SAPRC99 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHCHO_R_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  3.138E-20, 3.313E-20, 1.595E-20, 3.117E-20,
     &  1.660E-20, 7.129E-22, 0.000E+00 /
      DATA ( CS_REF( IHCHO_R_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  3.138E-20, 3.313E-20, 1.595E-20, 3.117E-20,
     &  1.660E-20, 7.129E-22, 0.000E+00 /
      DATA ( CS_REF( IHCHO_R_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  3.138E-20, 3.313E-20, 1.595E-20, 3.117E-20,
     &  1.660E-20, 7.129E-22, 0.000E+00 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( IHCHO_R_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.7543, 0.7794, 0.7737, 0.6830, 0.2322, 0.0000, 0.0000 /
      DATA ( QY_REF( IHCHO_R_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.7543, 0.7794, 0.7737, 0.6830, 0.2322, 0.0000, 0.0000 /
      DATA ( QY_REF( IHCHO_R_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.7543, 0.7794, 0.7737, 0.6830, 0.2322, 0.0000, 0.0000 /

C...(HCHO_M) :: HCHO + hv -> H2 + CO
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IHCHO_M ),ITTR=1,3 ) / 223.0, 273.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHCHO_M, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.997909E-20, 0.892819E-20, 0.424041E-20, 0.137440E-19,
     &  0.833517E-20, 0.134402E-21, 0.000000E+00 /
      DATA ( CS_REF( IHCHO_M, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.997123E-20, 0.892935E-20, 0.423299E-20, 0.137719E-19,
     &  0.835050E-20, 0.134213E-21, 0.000000E+00 /
      DATA ( CS_REF( IHCHO_M, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.996730E-20, 0.892993E-20, 0.422929E-20, 0.137859E-19,
     &  0.835817E-20, 0.134118E-21, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IHCHO_M, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.3130, 0.2980, 0.2820, 0.2890, 0.4070, 0.0990, 0.0000 /
      DATA ( QY_REF( IHCHO_M, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.3130, 0.2980, 0.2820, 0.2890, 0.4070, 0.0990, 0.0000 /
      DATA ( QY_REF( IHCHO_M, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.3130, 0.2980, 0.2820, 0.2890, 0.4070, 0.0990, 0.0000 /

C...(HCHO_M_SAPRC99) :: HCHO + hv = H2 + CO
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHCHO_M_SAPRC99 ), ITTR=1,3 ) /
     &  298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHCHO_M_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  3.138E-20, 3.313E-20, 1.595E-20, 3.117E-20,
     &  1.660E-20, 7.129E-22, 0.000E+00 /
      DATA ( CS_REF( IHCHO_M_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  3.138E-20, 3.313E-20, 1.595E-20, 3.117E-20,
     &  1.660E-20, 7.129E-22, 0.000E+00 /
      DATA ( CS_REF( IHCHO_M_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  3.138E-20, 3.313E-20, 1.595E-20, 3.117E-20,
     &  1.660E-20, 7.129E-22, 0.000E+00 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( IHCHO_M_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.2255, 0.2143, 0.2257, 0.3170, 0.5608, 0.1585, 0.0000 /
      DATA ( QY_REF( IHCHO_M_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.2255, 0.2143, 0.2257, 0.3170, 0.5608, 0.1585, 0.0000 /
      DATA ( QY_REF( IHCHO_M_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.2255, 0.2143, 0.2257, 0.3170, 0.5608, 0.1585, 0.0000 /

C...(CCHO_R) :: CH3CHO + hv -> CH3 + HCO
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,ICCHO_R ),ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( ICCHO_R, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.2155E-19, 0.1460E-19, 0.8378E-20, 0.3339E-20,
     &  0.1786E-21, 0.0000E+00, 0.0000E+00 /
      DATA ( CS_REF( ICCHO_R, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.2155E-19, 0.1460E-19, 0.8378E-20, 0.3339E-20,
     &  0.1786E-21, 0.0000E+00, 0.0000E+00 /
      DATA ( CS_REF( ICCHO_R, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.2155E-19, 0.1460E-19, 0.8378E-20, 0.3339E-20,
     &  0.1786E-21, 0.0000E+00, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( ICCHO_R, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICCHO_R, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICCHO_R, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(CCHO_R_SAPRC99) :: CCHO + hv -> CH3 + CHO
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07
C...  Cross sections below represent the effective values.
C...    The multiplication of absorption cross sections and quantum
C..     yields was performed on the finer wavelength data.

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, ICCHO_R_SAPRC99 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  effective cross sections (=cs*qy)

      DATA ( CS_REF( ICCHO_R_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  2.141E-20, 1.457E-20, 8.328E-21, 3.311E-21, 
     &  1.755E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( ICCHO_R_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  2.141E-20, 1.457E-20, 8.328E-21, 3.311E-21, 
     &  1.755E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( ICCHO_R_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  2.141E-20, 1.457E-20, 8.328E-21, 3.311E-21, 
     &  1.755E-22, 0.000E+00, 0.000E+00 /

C...  quantum yields (set to 1.0; included in CS values)

      DATA ( QY_REF( ICCHO_R_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICCHO_R_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICCHO_R_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(C2CHO) :: C2H5CHO + hv -> C2H5 + HCO
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IC2CHO ), ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IC2CHO, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.5538E-19, 0.4632E-19, 0.3577E-19, 0.2444E-19,
     &  0.5902E-20, 0.1250E-22, 0.0000E+00 /
      DATA ( CS_REF( IC2CHO, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.5538E-19, 0.4632E-19, 0.3577E-19, 0.2444E-19,
     &  0.5902E-20, 0.1250E-22, 0.0000E+00 /
      DATA ( CS_REF( IC2CHO, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.5538E-19, 0.4632E-19, 0.3577E-19, 0.2444E-19,
     &  0.5902E-20, 0.1250E-22, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( IC2CHO, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IC2CHO, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IC2CHO, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(C2CHO_SAPRC99) :: C2CHO + hv -> C2H5 + CHO
C...  (C2H5CHO)
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07
C...  Cross sections below represent the effective values.
C...    The multiplication of absorption cross sections and quantum
C..     yields was performed on the finer wavelength data.

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IC2CHO_SAPRC99 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  effective cross sections (=cs*qy)

      DATA ( CS_REF( IC2CHO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 2.814E-20, 3.702E-20, 2.119E-20, 1.072E-20, 
     & 1.395E-21, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IC2CHO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 2.814E-20, 3.702E-20, 2.119E-20, 1.072E-20, 
     & 1.395E-21, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IC2CHO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 2.814E-20, 3.702E-20, 2.119E-20, 1.072E-20, 
     & 1.395E-21, 0.000E+00, 0.000E+00 /

C...  quantum yields (set to 1.0; included in CS values)

      DATA ( QY_REF( IC2CHO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IC2CHO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IC2CHO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(ACETONE) :: ACETONE + hv -> CH3CO + CH3
C...  treated as a special case within PHOTMOD---data are loaded there
C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IACETONE ),ITTR=1,3) / 235.0, 273.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IACETONE, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     &  0.000000E+00, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IACETONE, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     &  0.000000E+00, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IACETONE, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
     &  0.000000E+00, 0.000000E+00, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IACETONE, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IACETONE, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IACETONE, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(ACETONE_SAPRC99) :: ACETONE + HV -> CH3CO. + CH3.
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IACETONE_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IACETONE_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 3.553000E-20, 2.339000E-20, 1.399000E-20, 7.494000E-21,
     & 8.422000E-22, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IACETONE_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 3.553000E-20, 2.339000E-20, 1.399000E-20, 7.494000E-21,
     & 8.422000E-22, 0.000000E+00, 0.000000E+00 /
      DATA ( CS_REF( IACETONE_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 3.553000E-20, 2.339000E-20, 1.399000E-20, 7.494000E-21,
     & 8.422000E-22, 0.000000E+00, 0.000000E+00 /

C...  quantum yields

      DATA ( QY_REF( IACETONE_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 2.229000E-01, 1.137000E-01, 5.744000E-02, 2.855000E-02,
     & 4.456000E-03, 0.000000E+00, 0.000000E+00 /
      DATA ( QY_REF( IACETONE_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 2.229000E-01, 1.137000E-01, 5.744000E-02, 2.855000E-02,
     & 4.456000E-03, 0.000000E+00, 0.000000E+00 /
      DATA ( QY_REF( IACETONE_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 2.229000E-01, 1.137000E-01, 5.744000E-02, 2.855000E-02,
     & 4.456000E-03, 0.000000E+00, 0.000000E+00 /

C...(KETONE) :: CH3COC2H5 + hv -> ACO3 + ETH
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IKETONE ),ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IKETONE, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.4157E-19, 0.2685E-19, 0.1565E-19, 0.7729E-20,
     &  0.8244E-21, 0.3711E-24, 0.0000E+00 /
      DATA ( CS_REF( IKETONE, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.4157E-19, 0.2685E-19, 0.1565E-19, 0.7729E-20,
     &  0.8244E-21, 0.3711E-24, 0.0000E+00 /
      DATA ( CS_REF( IKETONE, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.4157E-19, 0.2685E-19, 0.1565E-19, 0.7729E-20,
     &  0.8244E-21, 0.3711E-24, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( IKETONE, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IKETONE, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IKETONE, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(KETONE_SAPRC99) :: CH3COC2H5 + hv -> ACO3 + ETH
C...  Methyl Ethyl Ketone Absorption Cross Sections
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IKETONE_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IKETONE_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 4.241E-20, 2.705E-20, 1.552E-20, 7.629E-21,
     & 7.531E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IKETONE_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 4.241E-20, 2.705E-20, 1.552E-20, 7.629E-21,
     & 7.531E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IKETONE_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 4.241E-20, 2.705E-20, 1.552E-20, 7.629E-21,
     & 7.531E-22, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IKETONE_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IKETONE_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IKETONE_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(COOH) :: COOH + hv -> products
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, ICOOH ), ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( ICOOH, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.5598E-20, 0.3573E-20, 0.2437E-20, 0.1756E-20,
     &  0.7428E-21, 0.4236E-22, 0.0000E+00 /
      DATA ( CS_REF( ICOOH, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.5598E-20, 0.3573E-20, 0.2437E-20, 0.1756E-20,
     &  0.7428E-21, 0.4236E-22, 0.0000E+00 /
      DATA ( CS_REF( ICOOH, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.5598E-20, 0.3573E-20, 0.2437E-20, 0.1756E-20,
     &  0.7428E-21, 0.4236E-22, 0.0000E+00 /

C...  quantum yields

      DATA ( QY_REF( ICOOH, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICOOH, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICOOH, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(COOH_SAPRC99) :: CH3OOH + hv -> products
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, ICOOH_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( ICOOH_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 5.599E-21, 3.514E-21, 2.397E-21, 1.696E-21,
     & 7.240E-22, 5.345E-23, 0.000E+00 /
      DATA ( CS_REF( ICOOH_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 5.599E-21, 3.514E-21, 2.397E-21, 1.696E-21,
     & 7.240E-22, 5.345E-23, 0.000E+00 /
      DATA ( CS_REF( ICOOH_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 5.599E-21, 3.514E-21, 2.397E-21, 1.696E-21,
     & 7.240E-22, 5.345E-23, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( ICOOH_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICOOH_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICOOH_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(CLY_R) :: GLY_R + hv -> 2HCO
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IGLY_R ), ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IGLY_R, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.326037E-19, 0.314742E-19, 0.273529E-19, 0.234199E-19,
     &  0.830616E-20, 0.185293E-19, 0.706128E-20 /
      DATA ( CS_REF( IGLY_R, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.326037E-19, 0.314742E-19, 0.273529E-19, 0.234199E-19,
     &  0.830616E-20, 0.185293E-19, 0.706128E-20 /
      DATA ( CS_REF( IGLY_R, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.326037E-19, 0.314742E-19, 0.273529E-19, 0.234199E-19,
     &  0.830616E-20, 0.185293E-19, 0.706128E-20 /

C...  quantum yields

      DATA ( QY_REF( IGLY_R, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.4000, 0.4000, 0.4000, 0.4000, 0.4000, 0.0590, 0.0290 /
      DATA ( QY_REF( IGLY_R, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.4000, 0.4000, 0.4000, 0.4000, 0.4000, 0.0590, 0.0290 /
      DATA ( QY_REF( IGLY_R, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.4000, 0.4000, 0.4000, 0.4000, 0.4000, 0.0590, 0.0290 /

C...(GLY_R_SAPRC99) :: Glyoxal + hv = 2 HCO
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IGLY_R_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IGLY_R_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 3.311E-20, 3.061E-20, 2.763E-20, 2.058E-20,
     & 6.475E-21, 2.070E-20, 6.091E-22 /
      DATA ( CS_REF( IGLY_R_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 3.311E-20, 3.061E-20, 2.763E-20, 2.058E-20,
     & 6.475E-21, 2.070E-20, 6.091E-22 /
      DATA ( CS_REF( IGLY_R_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 3.311E-20, 3.061E-20, 2.763E-20, 2.058E-20,
     & 6.475E-21, 2.070E-20, 6.091E-22 /

C...  quantum yields

      DATA ( QY_REF( IGLY_R_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 1.000000E+00, 5.375000E-01, 0.000000E+00 /
      DATA ( QY_REF( IGLY_R_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 1.000000E+00, 5.375000E-01, 0.000000E+00 /
      DATA ( QY_REF( IGLY_R_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 1.000000E+00, 5.375000E-01, 0.000000E+00 /

C...(GLY_ABS) :: GLY_ABS + hv ->
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IGLY_ABS ),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IGLY_ABS, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.326037E-19, 0.314742E-19, 0.273529E-19, 0.234199E-19,
     &  0.830616E-20, 0.185293E-19, 0.706128E-20 /
      DATA ( CS_REF( IGLY_ABS, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.326037E-19, 0.314742E-19, 0.273529E-19, 0.234199E-19,
     &  0.830616E-20, 0.185293E-19, 0.706128E-20 /
      DATA ( CS_REF( IGLY_ABS, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.326037E-19, 0.314742E-19, 0.273529E-19, 0.234199E-19,
     &  0.830616E-20, 0.185293E-19, 0.706128E-20 /

C...  quantum yields

      DATA ( QY_REF( IGLY_ABS, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IGLY_ABS, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IGLY_ABS, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(GLY_ABS_SAPRC99) :: GLY_ABS + hv ->
C...  Glyoxal Absorption Cross Sections
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IGLY_ABS_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IGLY_ABS_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 3.311E-20, 3.061E-20, 2.763E-20, 2.058E-20,
     & 6.475E-21, 2.070E-20, 7.734E-21 /
      DATA ( CS_REF( IGLY_ABS_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 3.311E-20, 3.061E-20, 2.763E-20, 2.058E-20,
     & 6.475E-21, 2.070E-20, 7.734E-21 /
      DATA ( CS_REF( IGLY_ABS_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 3.311E-20, 3.061E-20, 2.763E-20, 2.058E-20,
     & 6.475E-21, 2.070E-20, 7.734E-21 /

C...  quantum yields

      DATA ( QY_REF( IGLY_ABS_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IGLY_ABS_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IGLY_ABS_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(MGLY_ABS) :: CH3COCHO + hv -> products (abs.CS only, qy=1)
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IMGLY_ABS),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IMGLY_ABS, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.4376E-19, 0.3464E-19, 0.2424E-19, 0.1789E-19,
     &  0.6219E-20, 0.2241E-19, 0.9673E-20 /
      DATA ( CS_REF( IMGLY_ABS, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.4376E-19, 0.3464E-19, 0.2424E-19, 0.1789E-19,
     &  0.6219E-20, 0.2241E-19, 0.9673E-20 /
      DATA ( CS_REF( IMGLY_ABS, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.4376E-19, 0.3464E-19, 0.2424E-19, 0.1789E-19,
     &  0.6219E-20, 0.2241E-19, 0.9673E-20 /

C...  quantum yields

      DATA ( QY_REF( IMGLY_ABS, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IMGLY_ABS, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IMGLY_ABS, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(MGLY_ABS_SAPRC99) :: CH3COCHO + hv -> products (abs.CS only, qy=1)
C...  Methyl Glyoxal Absorption Cross Sections
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IMGLY_ABS_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IMGLY_ABS_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 4.401E-20, 3.501E-20, 2.349E-20, 1.812E-20,
     & 6.010E-21, 3.752E-20, 7.985E-21 /
      DATA ( CS_REF( IMGLY_ABS_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 4.401E-20, 3.501E-20, 2.349E-20, 1.812E-20,
     & 6.010E-21, 3.752E-20, 7.985E-21 /
      DATA ( CS_REF( IMGLY_ABS_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 4.401E-20, 3.501E-20, 2.349E-20, 1.812E-20,
     & 6.010E-21, 3.752E-20, 7.985E-21 /

C...  quantum yields

      DATA ( QY_REF( IMGLY_ABS_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IMGLY_ABS_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IMGLY_ABS_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(MGLY_ADJ) :: CH3COCHO + hv -> products
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IMGLY_ADJ),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IMGLY_ADJ, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.4376E-19, 0.3464E-19, 0.2424E-19, 0.1789E-19,
     &  0.6219E-20, 0.2241E-19, 0.9673E-20 /
      DATA ( CS_REF( IMGLY_ADJ, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.4376E-19, 0.3464E-19, 0.2424E-19, 0.1789E-19,
     &  0.6219E-20, 0.2241E-19, 0.9673E-20 /
      DATA ( CS_REF( IMGLY_ADJ, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.4376E-19, 0.3464E-19, 0.2424E-19, 0.1789E-19,
     &  0.6219E-20, 0.2241E-19, 0.9673E-20 /

C...  quantum yields

      DATA ( QY_REF( IMGLY_ADJ, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.1500 /
      DATA ( QY_REF( IMGLY_ADJ, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.1500 /
      DATA ( QY_REF( IMGLY_ADJ, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.1500 /

C...(MGLY_ADJ_SAPRC99) :: MGLY + HV = PRODUCTS
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IMGLY_ADJ_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IMGLY_ADJ_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 4.401000E-20, 3.501000E-20, 2.349000E-20, 1.812000E-20,
     & 6.010000E-21, 3.752000E-20, 1.468000E-21 /
      DATA ( CS_REF( IMGLY_ADJ_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 4.401000E-20, 3.501000E-20, 2.349000E-20, 1.812000E-20,
     & 6.010000E-21, 3.752000E-20, 1.468000E-21 /
      DATA ( CS_REF( IMGLY_ADJ_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 4.401000E-20, 3.501000E-20, 2.349000E-20, 1.812000E-20,
     & 6.010000E-21, 3.752000E-20, 1.468000E-21 /

C...  quantum yields

      DATA ( QY_REF( IMGLY_ADJ_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 9.993000E-01, 3.758000E-01, 0.000000E+00 /
      DATA ( QY_REF( IMGLY_ADJ_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 9.993000E-01, 3.758000E-01, 0.000000E+00 /
      DATA ( QY_REF( IMGLY_ADJ_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 9.993000E-01, 3.758000E-01, 0.000000E+00 /

C...(MGLY_IUPAC04) :: CH3COCHO + hv ---> CH3CO + HCO
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation;
C...  IUPAC Stern-Volmer expression
C...  Data Sheet P6_CH3COCHO+hv.pdf, updated 16th Jan, 2003
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IMGLY_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IMGLY_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 4.371E-20, 3.459E-20, 2.417E-20, 1.785E-20,
     & 6.196E-21, 3.754E-20, 7.957E-21 /
      DATA ( CS_REF( IMGLY_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 4.371E-20, 3.459E-20, 2.417E-20, 1.785E-20,
     & 6.196E-21, 3.754E-20, 7.957E-21 /
      DATA ( CS_REF( IMGLY_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 4.371E-20, 3.459E-20, 2.417E-20, 1.785E-20,
     & 6.196E-21, 3.754E-20, 7.957E-21 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( IMGLY_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.2147, 0.0268 /
      DATA ( QY_REF( IMGLY_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.2147, 0.0268 /
      DATA ( QY_REF( IMGLY_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.2147, 0.0268 /

C...(BACL_ADJ) :: BACL_ADJ + hv -> products
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IBACL_ADJ),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IBACL_ADJ, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.255E-19, 0.159E-19, 0.903E-20, 0.601E-20,
     &  0.470E-20, 0.324E-20, 0.471E-20 /
      DATA ( CS_REF( IBACL_ADJ, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.255E-19, 0.159E-19, 0.903E-20, 0.601E-20,
     &  0.470E-20, 0.324E-20, 0.471E-20 /
      DATA ( CS_REF( IBACL_ADJ, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.255E-19, 0.159E-19, 0.903E-20, 0.601E-20,
     &  0.470E-20, 0.324E-20, 0.471E-20 /

C...  quantum yields

      DATA ( QY_REF( IBACL_ADJ, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IBACL_ADJ, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IBACL_ADJ, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(BACL_ADJ_SAPRC99) :: BACL + hv = products
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IBACL_ADJ_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IBACL_ADJ_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 2.605000E-20, 1.582000E-20, 8.965000E-21, 5.994000E-21,
     & 4.681000E-21, 3.258000E-20, 4.712000E-21 /
      DATA ( CS_REF( IBACL_ADJ_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 2.605000E-20, 1.582000E-20, 8.965000E-21, 5.994000E-21,
     & 4.681000E-21, 3.258000E-20, 4.712000E-21 /
      DATA ( CS_REF( IBACL_ADJ_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 2.605000E-20, 1.582000E-20, 8.965000E-21, 5.994000E-21,
     & 4.681000E-21, 3.258000E-20, 4.712000E-21 /

C...  quantum yields

      DATA ( QY_REF( IBACL_ADJ_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 1.000000E+00, 5.260000E-01, 9.086000E-04 /
      DATA ( QY_REF( IBACL_ADJ_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 1.000000E+00, 5.260000E-01, 9.086000E-04 /
      DATA ( QY_REF( IBACL_ADJ_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.000000E+00, 1.000000E+00, 1.000000E+00, 1.000000E+00,
     & 1.000000E+00, 5.260000E-01, 9.086000E-04 /

C...(BZCHO) :: BZCHO + hv -> products
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IBZCHO ),ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IBZCHO, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.919E-19, 0.661E-19, 0.673E-19,
     &  0.825E-19, 0.279E-19, 0.000E+00 /
      DATA ( CS_REF( IBZCHO, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.919E-19, 0.661E-19, 0.673E-19,
     &  0.825E-19, 0.279E-19, 0.000E+00 /
      DATA ( CS_REF( IBZCHO, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.000E+00, 0.919E-19, 0.661E-19, 0.673E-19,
     &  0.825E-19, 0.279E-19, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IBZCHO, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IBZCHO, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IBZCHO, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(BZCHO_SAPRC99) :: BZCHO + hv -> products
C...  Benzaldehyde absorption coefs in n-Hexane
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IBZCHO_SAPRC99 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IBZCHO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 9.261E-20, 6.614E-20, 6.732E-20,
     & 8.250E-20, 2.787E-20, 0.000E+00 /
      DATA ( CS_REF( IBZCHO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 9.261E-20, 6.614E-20, 6.732E-20,
     & 8.250E-20, 2.787E-20, 0.000E+00 /
      DATA ( CS_REF( IBZCHO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.000E+00, 9.261E-20, 6.614E-20, 6.732E-20,
     & 8.250E-20, 2.787E-20, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IBZCHO_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IBZCHO_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IBZCHO_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(ACROLEIN) :: C3H4O + hv -> products
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IACROLEIN),ITTR=1,3) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IACROLEIN, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.204E-19, 0.315E-19, 0.408E-19, 0.483E-19,
     &  0.575E-19, 0.118E-19, 0.000E+00 /
      DATA ( CS_REF( IACROLEIN, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.204E-19, 0.315E-19, 0.408E-19, 0.483E-19,
     &  0.575E-19, 0.118E-19, 0.000E+00 /
      DATA ( CS_REF( IACROLEIN, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.204E-19, 0.315E-19, 0.408E-19, 0.483E-19,
     &  0.575E-19, 0.118E-19, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IACROLEIN, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IACROLEIN, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IACROLEIN, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(ACROLEIN_SAPRC99) :: C3H4O + hv -> products
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IACROLEIN_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IACROLEIN_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 2.128E-20, 3.153E-20, 4.096E-20, 4.842E-20,
     & 5.748E-20, 1.178E-20, 0.000E+00 /
      DATA ( CS_REF( IACROLEIN_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 2.128E-20, 3.153E-20, 4.096E-20, 4.842E-20,
     & 5.748E-20, 1.178E-20, 0.000E+00 /
      DATA ( CS_REF( IACROLEIN_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 2.128E-20, 3.153E-20, 4.096E-20, 4.842E-20,
     & 5.748E-20, 1.178E-20, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IACROLEIN_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IACROLEIN_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IACROLEIN_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(IC3ONO2) :: i-C3H7ONO2 + hv -> products
C...  fzb

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR,IIC3ONO2),ITTR=1,3 ) / 350.0, 350.0, 350.0 /

C...  absorption cross sections

      DATA ( CS_REF( IIC3ONO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  0.121E-19, 0.635E-20, 0.328E-20, 0.171E-20,
     &  0.270E-21, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IIC3ONO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  0.121E-19, 0.635E-20, 0.328E-20, 0.171E-20,
     &  0.270E-21, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IIC3ONO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  0.121E-19, 0.635E-20, 0.328E-20, 0.171E-20,
     &  0.270E-21, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IIC3ONO2, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IIC3ONO2, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IIC3ONO2, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(IC3ONO2_SAPRC99) i-C3H7ONO2 + hv = products
C...  SAPRC-99 Photolysis data.  Supplied by William P. L. Carter.
C...  Created from PhotDat.xls on 29-Jan-2000 10:07

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IIC3ONO2_SAPRC99 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IIC3ONO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     & 1.248E-20, 6.322E-21, 3.254E-21, 1.702E-21,
     & 2.680E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IIC3ONO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     & 1.248E-20, 6.322E-21, 3.254E-21, 1.702E-21,
     & 2.680E-22, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IIC3ONO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     & 1.248E-20, 6.322E-21, 3.254E-21, 1.702E-21,
     & 2.680E-22, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IIC3ONO2_SAPRC99, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IIC3ONO2_SAPRC99, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IIC3ONO2_SAPRC99, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(CL2_FASTJX) :: CL2 + hv -> 2*CL
C...  from Fast-JXv5.7r
C...  (fort10_06b) UCI J-v8.4 (4/06) JPL02+irHNO4+IUPAC/NO2/VOC+Blitz +flux2006

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, ICL2_FASTJX ), ITTR=1,3 ) / 200.0, 300.0, 300.0 /

C...  absorption cross sections

      DATA ( CS_REF( ICL2_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  8.147E-20, 1.387E-19, 1.880E-19, 2.284E-19,
     &  2.549E-19, 6.256E-20, 6.399E-22 /
      DATA ( CS_REF( ICL2_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  8.674E-20, 1.403E-19, 1.849E-19, 2.205E-19,
     &  2.436E-19, 6.516E-20, 6.748E-22 /
      DATA ( CS_REF( ICL2_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  8.674E-20, 1.403E-19, 1.849E-19, 2.205E-19,
     &  2.436E-19, 6.516E-20, 6.748E-22 /

C...  quantum yields

      DATA ( QY_REF( ICL2_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICL2_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICL2_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(CL2_IUPAC04) :: CL2 + hv = 2*CL

C...IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet PCl11 Website: 15th December 2000
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, ICL2_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( ICL2_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 8.875E-20, 1.404E-19, 1.848E-19, 2.187E-19,
     & 2.411E-19, 6.480E-20, 6.167E-22 /
      DATA ( CS_REF( ICL2_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 8.875E-20, 1.404E-19, 1.848E-19, 2.187E-19,
     & 2.411E-19, 6.480E-20, 6.167E-22 /
      DATA ( CS_REF( ICL2_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 8.875E-20, 1.404E-19, 1.848E-19, 2.187E-19,
     & 2.411E-19, 6.480E-20, 6.167E-22 /

C...  quantum yields

      DATA ( QY_REF( ICL2_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICL2_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( ICL2_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HOCL_FASTJX) :: HOCL + hv -> HO + CL
C...  from Fast-JXv5.7r
C...  (fort10_06b) UCI J-v8.4 (4/06) JPL02+irHNO4+IUPAC/NO2/VOC+Blitz +flux2006

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHOCL_FASTJX ), ITTR=1,3 ) / 297.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHOCL_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  5.567E-20, 6.067E-20, 5.957E-20, 5.378E-20,
     &  3.123E-20, 6.548E-21, 1.195E-23 /
      DATA ( CS_REF( IHOCL_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  5.567E-20, 6.067E-20, 5.957E-20, 5.378E-20,
     &  3.123E-20, 6.548E-21, 1.195E-23 /
      DATA ( CS_REF( IHOCL_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  5.567E-20, 6.067E-20, 5.957E-20, 5.378E-20,
     &  3.123E-20, 6.548E-21, 1.195E-23 /

C...  quantum yields

      DATA ( QY_REF( IHOCL_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHOCL_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHOCL_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(HOCL_IUPAC04) :: HOCL + hv = HO + CL
C...  IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet PCl2 Website: 15th December 2000
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IHOCL_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IHOCL_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 5.578E-20, 6.068E-20, 5.953E-20, 5.373E-20,
     & 3.124E-20, 6.487E-21, 1.169E-23 /
      DATA ( CS_REF( IHOCL_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 5.578E-20, 6.068E-20, 5.953E-20, 5.373E-20,
     & 3.124E-20, 6.487E-21, 1.169E-23 /
      DATA ( CS_REF( IHOCL_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 5.578E-20, 6.068E-20, 5.953E-20, 5.373E-20,
     & 3.124E-20, 6.487E-21, 1.169E-23 /

C...  quantum yields

      DATA ( QY_REF( IHOCL_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHOCL_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IHOCL_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(PAN_FASTJX) :: CH3C(O)OONO2 + hv -> CH3C(O)OO + NO2
C...  from Fast-JXv5.7r
C...  (fort10_06b) UCI J-v8.4 (4/06) JPL02+irHNO4+IUPAC/NO2/VOC+Blitz +flux2006

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IPAN_FASTJX ), ITTR=1,3 ) / 250.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IPAN_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  2.438E-21, 9.252E-22, 4.355E-22, 2.288E-22,
     &  5.480E-23, 6.941E-25, 0.000E+00 /
      DATA ( CS_REF( IPAN_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  3.555E-21, 1.399E-21, 6.750E-22, 3.627E-22,
     &  9.257E-23, 1.265E-24, 0.000E+00 /
      DATA ( CS_REF( IPAN_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  3.555E-21, 1.399E-21, 6.750E-22, 3.627E-22,
     &  9.257E-23, 1.265E-24, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IPAN_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IPAN_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IPAN_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(PAN_IUPAC04) :: CH3C(O)OONO2 + HV = CH3C(O)OO + NO2
C...  Derived from IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet P21_CH3C(O)OONO2+hv, updated 16th July 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IPAN_IUPAC04 ), ITTR=1,3 ) /
     & 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IPAN_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 3.500E-21, 1.410E-21, 6.665E-22, 3.626E-22,
     & 9.164E-23, 1.295E-24, 0.000E+00 /
      DATA ( CS_REF( IPAN_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 3.500E-21, 1.410E-21, 6.665E-22, 3.626E-22,
     & 9.164E-23, 1.295E-24, 0.000E+00 /
      DATA ( CS_REF( IPAN_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 3.500E-21, 1.410E-21, 6.665E-22, 3.626E-22,
     & 9.164E-23, 1.295E-24, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IPAN_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IPAN_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IPAN_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(N2O5_FASTJX) :: N2O5 + hv -> NO2 + NO3
C...  from Fast-JXv5.7r
C...  (fort10_06b) UCI J-v8.4 (4/06) JPL02+irHNO4+IUPAC/NO2/VOC+Blitz +flux2006

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IN2O5_FASTJX ), ITTR=1,3 ) / 225.0, 300.0, 300.0 /

C...  absorption cross sections

      DATA ( CS_REF( IN2O5_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  3.823E-20, 1.998E-20, 1.170E-20, 7.246E-21,
     &  2.286E-21, 1.173E-22, 0.000E+00 /
      DATA ( CS_REF( IN2O5_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  5.404E-20, 3.317E-20, 2.226E-20, 1.551E-20,
     &  6.389E-21, 5.481E-22, 0.000E+00 /
      DATA ( CS_REF( IN2O5_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  5.404E-20, 3.317E-20, 2.226E-20, 1.551E-20,
     &  6.389E-21, 5.481E-22, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IN2O5_FASTJX, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IN2O5_FASTJX, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IN2O5_FASTJX, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(N2O5_IUPAC04) :: N2O5 + hv = NO2 + NO3
C...  From IUPAC Subcommittee on Gas Kinetic Data Evaluation
C...  Data Sheet PNOx7_N2O5, updated 16th July 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF(ITTR, IN2O5_IUPAC04 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IN2O5_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 5.264E-20, 3.274E-20, 2.180E-20, 1.509E-20,
     & 6.077E-21, 6.134E-22, 8.598E-26 /
      DATA ( CS_REF( IN2O5_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 5.264E-20, 3.274E-20, 2.180E-20, 1.509E-20,
     & 6.077E-21, 6.134E-22, 8.598E-26 /
      DATA ( CS_REF( IN2O5_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 5.264E-20, 3.274E-20, 2.180E-20, 1.509E-20,
     & 6.077E-21, 6.134E-22, 8.598E-26 /

C...  quantum yields
C...    effective quantum yields were computed by performing separate
C...    interval integrations for the cross sections and for the 
C...    effective cross sections (cs*qy) (calculated on the finer
C...    wavelength grid.  The effective quantum yield values below
C...    were then calculated for the 7 wavelength intervals by 
C...    dividing the effective cross sections by the interval average
C...    cross sections (eQY=eCS/CS).

      DATA ( QY_REF( IN2O5_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     & 0.9158, 0.9973, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IN2O5_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     & 0.9158, 0.9973, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IN2O5_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     & 0.9158, 0.9973, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(FMCL_IUPAC04) :: HC(O)Cl + hv -> HCO + CL

C...Derived from IUPAC Subcommittree on Gas Kinetic Data Evaluation
C...  Data Sheet PCl28 Website: 15th December 2000.
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk
C...  H. G. Libuda, F. Zabel, E. H. Fink, and K. H. Becker, J. Phys.
C...     Chem. 94, 5860 (1990).

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IFMCL_IUPAC04 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IFMCL_IUPAC04, 1, IWLR ), IWLR = 1,7 ) /
     &  5.234E-21, 1.410E-21, 2.207E-22, 8.239E-23,
     &  0.000E+00, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IFMCL_IUPAC04, 2, IWLR ), IWLR = 1,7 ) /
     &  5.234E-21, 1.410E-21, 2.207E-22, 8.239E-23,
     &  0.000E+00, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IFMCL_IUPAC04, 3, IWLR ), IWLR = 1,7 ) /
     &  5.234E-21, 1.410E-21, 2.207E-22, 8.239E-23,
     &  0.000E+00, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IFMCL_IUPAC04, 1, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IFMCL_IUPAC04, 2, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IFMCL_IUPAC04, 3, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(NTR_IUPAC04) :: i-C3H7ONO2 + hv -> iC3H7O + NO2

C...Derived from IUPAC Subcommittree on Gas Kinetic Data Evaluation
C...  Data Sheet P17_i-C3H7ONO2+hv, updated 16th July 2001
C...  Website: http://www.iupac-kinetic.ch.cam.ac.uk/

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, INTR_IUPAC04 ), ITTR =1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( INTR_IUPAC04, 1, IWLR ), IWLR = 1, 7 ) /
     &  1.222E-20, 6.081E-21, 2.981E-21, 1.475E-21,
     &  2.233E-22, 1.193E-24, 0.000E+00 /
      DATA ( CS_REF( INTR_IUPAC04, 2, IWLR ), IWLR = 1, 7 ) /
     &  1.222E-20, 6.081E-21, 2.981E-21, 1.475E-21,
     &  2.233E-22, 1.193E-24, 0.000E+00 /
      DATA ( CS_REF( INTR_IUPAC04, 3, IWLR ), IWLR = 1, 7 ) /
     &  1.222E-20, 6.081E-21, 2.981E-21, 1.475E-21,
     &  2.233E-22, 1.193E-24, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( INTR_IUPAC04, 1, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( INTR_IUPAC04, 2, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( INTR_IUPAC04, 3, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...(PACD_CB05) :: CH3COOOH + hv -> MEO2 + OH

C...Derived from data supplied by Greg Yarwood for CB05, 11/16/2007
C...  Gigu�re, P. A. and A. W. Olmos. Sur le spectre ultraviolet
C...  de l'acide perac�tique et l'hydrolyse des perac�tates.
C...  Can. J. Chem., 34, 689-691, 1956.

C...  reference temperatures (K)

      DATA ( TEMP_REF( ITTR, IPACD_CB05 ), ITTR=1,3 ) / 298.0, 298.0, 298.0 /

C...  absorption cross sections

      DATA ( CS_REF( IPACD_CB05, 1, IWLR ), IWLR = 1, 7 ) /
     &  7.713E-22, 4.112E-22, 2.411E-22, 1.520E-22,
     &  2.687E-23, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IPACD_CB05, 2, IWLR ), IWLR = 1, 7 ) /
     &  7.713E-22, 4.112E-22, 2.411E-22, 1.520E-22,
     &  2.687E-23, 0.000E+00, 0.000E+00 /
      DATA ( CS_REF( IPACD_CB05, 3, IWLR ), IWLR = 1, 7 ) /
     &  7.713E-22, 4.112E-22, 2.411E-22, 1.520E-22,
     &  2.687E-23, 0.000E+00, 0.000E+00 /

C...  quantum yields

      DATA ( QY_REF( IPACD_CB05, 1, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IPACD_CB05, 2, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      DATA ( QY_REF( IPACD_CB05, 3, IWLR ), IWLR = 1,7 ) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /

C...Setup the Mapping from CMAQ chemical reactions to the reference data

      INTEGER, PARAMETER :: NPHOT_MAP = 68  ! # mapped phot reactions

C...Names of the mapped photolysis reactions (available to chemical
C...  mechanisms) and their pointers to the reference photolysis rxn

      CHARACTER*16, SAVE :: PNAME_MAP( NPHOT_MAP )
      INTEGER, SAVE      :: PHOT_MAP ( NPHOT_MAP )

      DATA PNAME_MAP( 01 ), PHOT_MAP( 01 ) /'ACETONE         ', IACETONE          /
      DATA PNAME_MAP( 02 ), PHOT_MAP( 02 ) /'ACETONE_SAPRC99 ', IACETONE_SAPRC99  /
      DATA PNAME_MAP( 03 ), PHOT_MAP( 03 ) /'ACROLEIN        ', IACROLEIN         /
      DATA PNAME_MAP( 04 ), PHOT_MAP( 04 ) /'ACROLEIN_SAPRC99', IACROLEIN_SAPRC99 /
      DATA PNAME_MAP( 05 ), PHOT_MAP( 05 ) /'BACL_ADJ        ', IBACL_ADJ         /
      DATA PNAME_MAP( 06 ), PHOT_MAP( 06 ) /'BACL_ADJ_SAPRC99', IBACL_ADJ_SAPRC99 /
      DATA PNAME_MAP( 07 ), PHOT_MAP( 07 ) /'BZCHO           ', IBZCHO            /
      DATA PNAME_MAP( 08 ), PHOT_MAP( 08 ) /'BZCHO_SAPRC99   ', IBZCHO_SAPRC99    /
      DATA PNAME_MAP( 09 ), PHOT_MAP( 09 ) /'C2CHO           ', IC2CHO            /
      DATA PNAME_MAP( 10 ), PHOT_MAP( 10 ) /'C2CHO_SAPRC99   ', IC2CHO_SAPRC99    /
      DATA PNAME_MAP( 11 ), PHOT_MAP( 11 ) /'CCHO_R          ', ICCHO_R           /
      DATA PNAME_MAP( 12 ), PHOT_MAP( 12 ) /'CCHO_R_SAPRC99  ', ICCHO_R_SAPRC99   /
      DATA PNAME_MAP( 13 ), PHOT_MAP( 13 ) /'CL2_FASTJX      ', ICL2_FASTJX       /
      DATA PNAME_MAP( 14 ), PHOT_MAP( 14 ) /'CL2_IUPAC04     ', ICL2_IUPAC04      /
      DATA PNAME_MAP( 15 ), PHOT_MAP( 15 ) /'COOH            ', ICOOH             /
      DATA PNAME_MAP( 16 ), PHOT_MAP( 16 ) /'COOH_SAPRC99    ', ICOOH_SAPRC99     /
      DATA PNAME_MAP( 17 ), PHOT_MAP( 17 ) /'FMCL_IUPAC04    ', IFMCL_IUPAC04     /
      DATA PNAME_MAP( 18 ), PHOT_MAP( 18 ) /'GLY_ABS         ', IGLY_ABS          /
      DATA PNAME_MAP( 19 ), PHOT_MAP( 19 ) /'GLY_ABS_SAPRC99 ', IGLY_ABS_SAPRC99  /
      DATA PNAME_MAP( 20 ), PHOT_MAP( 20 ) /'GLY_R           ', IGLY_R            /
      DATA PNAME_MAP( 21 ), PHOT_MAP( 21 ) /'GLY_R_SAPRC99   ', IGLY_R_SAPRC99    /
      DATA PNAME_MAP( 22 ), PHOT_MAP( 22 ) /'H2O2            ', IH2O2             /
      DATA PNAME_MAP( 23 ), PHOT_MAP( 23 ) /'H2O2_SAPRC99    ', IH2O2_SAPRC99     /
      DATA PNAME_MAP( 24 ), PHOT_MAP( 24 ) /'HCHO_M          ', IHCHO_M           /
      DATA PNAME_MAP( 25 ), PHOT_MAP( 25 ) /'HCHO_M_SAPRC99  ', IHCHO_M_SAPRC99   /
      DATA PNAME_MAP( 26 ), PHOT_MAP( 26 ) /'HCHO_R          ', IHCHO_R           /
      DATA PNAME_MAP( 27 ), PHOT_MAP( 27 ) /'HCHO_R_SAPRC99  ', IHCHO_R_SAPRC99   /
      DATA PNAME_MAP( 28 ), PHOT_MAP( 28 ) /'HNO3            ', IHNO3             /
      DATA PNAME_MAP( 29 ), PHOT_MAP( 29 ) /'HNO3_IUPAC04    ', IHNO3_IUPAC04     /
      DATA PNAME_MAP( 30 ), PHOT_MAP( 30 ) /'HNO3_SAPRC99    ', IHNO3_SAPRC99     /
      DATA PNAME_MAP( 31 ), PHOT_MAP( 31 ) /'HO2NO2          ', IHO2NO2           /
      DATA PNAME_MAP( 32 ), PHOT_MAP( 32 ) /'HO2NO2_IUPAC04  ', IHO2NO2_IUPAC04   /
      DATA PNAME_MAP( 33 ), PHOT_MAP( 33 ) /'HO2NO2_SAPRC99  ', IHO2NO2_SAPRC99   /
      DATA PNAME_MAP( 34 ), PHOT_MAP( 34 ) /'HOCL_FASTJX     ', IHOCL_FASTJX      /
      DATA PNAME_MAP( 35 ), PHOT_MAP( 35 ) /'HOCL_IUPAC04    ', IHOCL_IUPAC04     /
      DATA PNAME_MAP( 36 ), PHOT_MAP( 36 ) /'HONO            ', IHONO             /
      DATA PNAME_MAP( 37 ), PHOT_MAP( 37 ) /'HONO_IUPAC04    ', IHONO_IUPAC04     /
      DATA PNAME_MAP( 38 ), PHOT_MAP( 38 ) /'HONO_NO         ', IHONO_NO          /
      DATA PNAME_MAP( 39 ), PHOT_MAP( 39 ) /'HONO_NO2        ', IHONO_NO2         /
      DATA PNAME_MAP( 40 ), PHOT_MAP( 40 ) /'HONO_NO2_SAPRC99', IHONO_NO2_SAPRC99 /
      DATA PNAME_MAP( 41 ), PHOT_MAP( 41 ) /'HONO_NO_SAPRC99 ', IHONO_NO_SAPRC99  /
      DATA PNAME_MAP( 42 ), PHOT_MAP( 42 ) /'IC3ONO2         ', IIC3ONO2          /
      DATA PNAME_MAP( 43 ), PHOT_MAP( 43 ) /'IC3ONO2_SAPRC99 ', IIC3ONO2_SAPRC99  /
      DATA PNAME_MAP( 44 ), PHOT_MAP( 44 ) /'KETONE          ', IKETONE           /
      DATA PNAME_MAP( 45 ), PHOT_MAP( 45 ) /'KETONE_SAPRC99  ', IKETONE_SAPRC99   /
      DATA PNAME_MAP( 46 ), PHOT_MAP( 46 ) /'MGLY_ABS        ', IMGLY_ABS         /
      DATA PNAME_MAP( 47 ), PHOT_MAP( 47 ) /'MGLY_ABS_SAPRC99', IMGLY_ABS_SAPRC99 /
      DATA PNAME_MAP( 48 ), PHOT_MAP( 48 ) /'MGLY_ADJ        ', IMGLY_ADJ         /
      DATA PNAME_MAP( 49 ), PHOT_MAP( 49 ) /'MGLY_ADJ_SAPRC99', IMGLY_ADJ_SAPRC99 /
      DATA PNAME_MAP( 50 ), PHOT_MAP( 50 ) /'MGLY_IUPAC04    ', IMGLY_IUPAC04     /
      DATA PNAME_MAP( 51 ), PHOT_MAP( 51 ) /'N2O5_FASTJX     ', IN2O5_FASTJX      /
      DATA PNAME_MAP( 52 ), PHOT_MAP( 52 ) /'N2O5_IUPAC04    ', IN2O5_IUPAC04     /
      DATA PNAME_MAP( 53 ), PHOT_MAP( 53 ) /'NO2             ', INO2              /
      DATA PNAME_MAP( 54 ), PHOT_MAP( 54 ) /'NO2_SAPRC99     ', INO2_SAPRC99      /
      DATA PNAME_MAP( 55 ), PHOT_MAP( 55 ) /'NO3NO           ', INO3NO            /
      DATA PNAME_MAP( 56 ), PHOT_MAP( 56 ) /'NO3NO_SAPRC99   ', INO3NO_SAPRC99    /
      DATA PNAME_MAP( 57 ), PHOT_MAP( 57 ) /'NO3NO2          ', INO3NO2           /
      DATA PNAME_MAP( 58 ), PHOT_MAP( 58 ) /'NO3NO2_SAPRC99  ', INO3NO2_SAPRC99   /
      DATA PNAME_MAP( 59 ), PHOT_MAP( 59 ) /'NTR_IUPAC04     ', INTR_IUPAC04      /
      DATA PNAME_MAP( 60 ), PHOT_MAP( 60 ) /'O3O1D           ', IO3O1D            /
      DATA PNAME_MAP( 61 ), PHOT_MAP( 61 ) /'O3O1D_SAPRC99   ', IO3O1D_SAPRC99    /
      DATA PNAME_MAP( 62 ), PHOT_MAP( 62 ) /'O3O3P           ', IO3O3P            /
      DATA PNAME_MAP( 63 ), PHOT_MAP( 63 ) /'O3O3P_SAPRC99   ', IO3O3P_SAPRC99    /
      DATA PNAME_MAP( 64 ), PHOT_MAP( 64 ) /'O3_O1D_IUPAC04  ', IO3_O1D_IUPAC04   /
      DATA PNAME_MAP( 65 ), PHOT_MAP( 65 ) /'O3_O3P_IUPAC04  ', IO3_O3P_IUPAC04   /
      DATA PNAME_MAP( 66 ), PHOT_MAP( 66 ) /'PACD_CB05       ', IPACD_CB05        /
      DATA PNAME_MAP( 67 ), PHOT_MAP( 67 ) /'PAN_FASTJX      ', IPAN_FASTJX       /
      DATA PNAME_MAP( 68 ), PHOT_MAP( 68 ) /'PAN_IUPAC04     ', IPAN_IUPAC04      /

      END MODULE CSQY_DATA
