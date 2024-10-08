
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/cloud/cloud_acm_ae5_txhg/calczen.F,v 1.1.1.1 2010/06/14 16:03:03 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CALCZEN ( LAT, LONG, MDATE, MTIME, CSZ )
        
C*********************************************************************
C
C  this subroutine calculates cosine of solar zenith angle for a
C     particulat time and location.  Must specify:
C  INPUT:
c       LAT - latitude in decimal degrees
C       LONG - longitude in decimal degrees
C       MDATE - Julian Date at Greenwich - specify year (yyyy), day (ddd),
C               five-digit integer:  yyyyddd
C       GMT  - Greenwich mean time - decimal military eg.
C               22.75 = 45 min after ten pm gmt
C  OUTPUT
C       Zenith
C
C*********************************************************************

      IMPLICIT NONE

      INCLUDE SUBST_CONST        ! commonly used constants

C...........ARGUMENTS and their descriptions:


      REAL         LAT                ! northern latitude (decimal degrees)
      REAL         LONG               ! western longitude (decimal degrees)
      INTEGER      MDATE              ! date (yyyyddd)
      INTEGER      MTIME              ! Current model time , coded HHMMSS
      REAL         CSZ                ! cosine zenith angle
            
C...........LOCAL VARIABLES and their descriptions:

      INTEGER      IIYEAR             ! year (yyyy)
      INTEGER      LEAP               ! # leap years since 1974
      INTEGER      NYEARS             ! years since 1974
      INTEGER      NOLEAP             ! # nonleap years since 1974
      INTEGER      IN                 ! month pointer (mm-1)
      INTEGER      IMTH               ! month (mm)
      INTEGER      I                  ! month index
      INTEGER      IDAY               ! day (dd)
      INTEGER      JD                 ! julian days since 1974 ref

      REAL         GMT                ! Greenwich mean time (hours, dec.milt)
      REAL         LBGMT              ! 
      REAL         LZGMT              ! 
      REAL         ML                 ! geometric mean longitude (deg)
      REAL         RLT                ! latitude (radians)
      REAL         YREF               ! number of days to 1974
      REAL         YR                 ! number of days to yyyy
      REAL         D                  ! jul. days since 1974 + hh frac
      REAL         RML                ! geometric mean long (rad)
      REAL         W                  ! mean longitude of perigee (deg)
      REAL         EC                 ! eccentricity 
      REAL         EPSI               ! mean obliquity of ecliptic (deg)
      REAL         YT                 ! tan**2 (mean obl. of eclipt.)
      REAL         PEPSI              ! mean obliquity of ecliptic (rad)
      REAL         CW                 ! cosine mean long. perigee
      REAL         WR                 ! mean longitude of perigee (rad)
      REAL         SW                 ! sine mean long. perigee
      REAL         SSW                ! sine 2*mean long. perigee
      REAL         EYT                !
      REAL         FEQT               !
      REAL         FEQT1              !
      REAL         FEQT2              !
      REAL         FEQT3              !
      REAL         FEQT4              !
      REAL         FEQT5              !
      REAL         FEQT6              !
      REAL         FEQT7              !
      REAL         REQT               !
      REAL         EQT                !
      REAL         RA                 ! right ascension (degrees)
      REAL         TAB                !
      REAL         RRA                ! right ascension (radians)
      REAL         RDECL              ! declination angle (rad)
      REAL         ZPT                ! 
c     REAL         ZR                 ! zenith angle (radians)
 
      INTEGER, EXTERNAL :: TIME2SEC   ! converts MTIME into seconds

C*********************************************************************
C     begin body of subroutine CALZEN2        

C...convert to radians

      RLT = LAT * PI180

C...parse date

      IIYEAR = MDATE / 1000
      IDAY   = MOD( MDATE, 1000)

C...count days from Dec.31,1973 to Jan 1, YEAR, then add to 2,442,047.5

      YREF =  2442047.5
      NYEARS = IIYEAR - 1974
      LEAP = ( NYEARS + 1 ) / 4
      IF ( NYEARS .LE. -1 ) LEAP = ( NYEARS - 2 ) / 4
      NOLEAP = NYEARS - LEAP
      YR = YREF + 365.0 * NOLEAP + 366.0 * LEAP


C...calculate julian days since 1974

      JD  = IDAY + ( YR - YREF ) 
      GMT = REAL( TIME2SEC( MTIME ) )/3600.0

      D = REAL( JD ) + GMT / 24.0

C...calc geom mean longitude

      ML = 279.2801988 + 0.9856473354 * D + 2.267E-13 * D * D
      RML = ML * PI180

C...calc equation of time in sec
C... w = mean long of perigee
C... e = eccentricity
C... epsi = mean obliquity of ecliptic

      W = 282.4932328 + 4.70684E-5 * D + 3.39E-13 * D * D
      WR = W * PI180
      EC = 1.6720041E-2 - 1.1444E-9 * D - 9.4E-17 * D * D
      EPSI = 23.44266511 - 3.5626E-7 * D - 1.23E-15 * D * D
      PEPSI = EPSI * PI180
      YT = ( TAN( PEPSI / 2.0 ) )**2
      CW = COS( WR )
      SW = SIN( WR )
      SSW = SIN( 2.0 * WR )
      EYT = 2.0 * EC * YT
      FEQT1 = SIN( RML ) * ( -EYT * CW - 2.0 * EC * CW )
      FEQT2 = COS( RML ) * ( 2.0 * EC * SW - EYT * SW )
      FEQT3 = SIN( 2.0 * RML ) * ( YT - ( 5.0 * EC**2 / 4.0 ) 
     &      * ( CW**2 - SW**2 ) )
      FEQT4 = COS( 2.0 * RML ) * ( 5.0 * EC**2 * SSW / 4.0 )
      FEQT5 = SIN( 3.0 * RML ) * ( EYT * CW )
      FEQT6 = COS( 3.0 * RML ) * ( -EYT * SW )
      FEQT7 = -SIN( 4.0 * RML ) * ( 0.5 * YT**2 )
      FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7
      EQT = FEQT * 13751.0

C...convert eq of time from sec to deg

      REQT = EQT / 240.0

C...calc right ascension in rads

      RA = ML - REQT
      RRA = RA * PI180

C...calc declination in rads, deg

      TAB = 0.43360 * SIN( RRA )
      RDECL = ATAN( TAB )

C...calc local hour angle

      LBGMT = 12.0 - EQT / 3600.0 + LONG * 24.0 / 360.0
      LZGMT = 15.0 * ( GMT - LBGMT )
      ZPT = LZGMT * PI180
      CSZ = SIN( RLT ) * SIN( RDECL ) + COS( RLT ) * COS( RDECL )
     &    * COS( ZPT )
c      ZENITH = ZR / PI180
     
      RETURN
      END
