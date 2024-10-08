
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_ddm3d/s_zfdbc.f,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REAL FUNCTION S_ZFDBC (C1, C2, V1, V2)

c Zero Flux Divergence Boundary Condition (See Jon Pleim's JGR (1991) paper)
c For sensitivity, negative values are possible, so MAX statement removed.
c To eliminate reflections and other boundary anomolies
C Problem if V1 is outflow, but V2 is inflow

      IMPLICIT NONE
      REAL SMALL
      PARAMETER (SMALL = 1.0E-03 )   ! for small wind speed (m/s)
      REAL C1, C2, V1, V2
 
      IF ( ABS( V1 ) .GE. SMALL ) THEN
         IF ( V1 * V2 .GT. 0.0 ) THEN
            S_ZFDBC =  C1 - V2 / V1 * (C2 - C1) 
            ELSE
            S_ZFDBC = C1         ! nothing changes for wind divergence at edge
            END IF
         ELSE
         S_ZFDBC = C1            ! nothing changes for small wind speed
         END IF

      RETURN
      END
