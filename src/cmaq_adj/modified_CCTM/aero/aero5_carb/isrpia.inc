
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_carb/isrpia.inc,v 1.1.1.1 2010/06/14 16:02:59 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C This is ISORROPIA Version 1.7 (3/15/06).  More information is available
C at the ISORROPIA homepage: http://nenes.eas.gatech.edu/ISORROPIA/
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** INCLUDE FILE 'ISRPIA.INC'
C *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
C     AND VARIABLES. 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NCOMP=5,NIONS=7,NGASAQ=3,NSLDS=9,NPAIR=13,NZSR=100,
     &           NERRMX=25)
C
C *** INPUT VARIABLES **************************************************
C
      INTEGER METSTBL
      COMMON /INPT/ W(NCOMP), WAER(NCOMP), TEMP, RH, IPROB, METSTBL,
     &              NADJ
C
C *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************
C
      COMMON /ZSR / AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR),
     &              AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR),
     &              AWSA(NZSR), AWLC(NZSR)
C
C *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************
C
      INTEGER WFTYP
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3, 
     &              DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      COMMON /MDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2,
     &              DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2,
     &              DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3,
     &              DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8,
     &              DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,
     &              WFTYP
C
C *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
C
      DOUBLE PRECISION MOLAL, MOLALR, M0
      REAL IONIC
      LOGICAL CALAOU, CALAIN, FRST, DRYF
      COMMON /IONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR), ZZ(NPAIR),
     &              Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR),
     &              GASAQ(NGASAQ),
     &              EPSACT,       COH,           CHNO3,       CHCL,         
     &              WATER,        IONIC,         IACALC,      
     &              FRST,         CALAIN,        CALAOU,      DRYF
C
C *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
C
      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4, 
     &              CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC
C
C *** VARIABLES FOR GAS PHASE ******************************************
C
      COMMON /GAS / GNH3, GHNO3, GHCL 
C
C *** EQUILIBRIUM CONSTANTS ********************************************
C
      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,
     &              XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,
     &              XK42
C
C *** MOLECULAR WEIGHTS ************************************************
C
      DOUBLE PRECISION IMW
      COMMON /OTHR/ R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
C
C *** SOLUTION/INFO VARIABLES ******************************************
C
      CHARACTER SCASE*15
      COMMON /CASE/ SULRATW, SULRAT, SODRAT, SCASE
C
      COMMON /SOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT
C
C *** ERROR SYSTEM *****************************************************
C
      CHARACTER ERRMSG*40
      INTEGER   ERRSTK, NOFER   
      LOGICAL   STKOFL   
      COMMON /EROR/ STKOFL, NOFER, ERRSTK(NERRMX), ERRMSG(NERRMX)
C
C *** GENERIC VARIABLES ************************************************
C
      CHARACTER VERSION*15
      COMMON /CGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION
C
C *** END OF INCLUDE FILE **********************************************
C
