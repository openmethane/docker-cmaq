
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline/czangle.F,v 1.1.1.1 2010/06/14 16:03:07 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE  CZANGLE( JDATE, JTIME, NX, NY, LAT, LON, COSZEN )

C----------------------------------------------------------------------
C Description:
C   Computes cosine of zenith angle for routine HRBIO()
C   The zenith angle is the angle measured from the zenith to the line of
C   sight to the sun
 
C Preconditions:
C   JDATE:JTIME represented in GMT
 
C Subroutines and Functions called:  none
 
C Revision History:
C   11/99: by Jeff Vukovich taken from v4.2 SMOKE prototype
C   10/06: yoj
C----------------------------------------------------------------------
C Modified from:

C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling System
C File: @(#)$Id: czangle.F,v 1.1.1.1 2010/06/14 16:03:07 sjr Exp $
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C smoke@unc.edu
C Pathname: $Source: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline/czangle.F,v $
C Last updated: $Date: 2010/06/14 16:03:07 $ 
C----------------------------------------------------------------------

      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
!include      SUBST_IODECL    # I/O definitions and declarations
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_CONST     ! constants

C Arguments:

      INTEGER, INTENT( IN )  :: JDATE   ! current simulation date (YYYYDDD)
      INTEGER, INTENT( IN )  :: JTIME   ! current simulation time (HHMMSS)
      INTEGER, INTENT( IN )  :: NX      ! no. columns
      INTEGER, INTENT( IN )  :: NY      ! no. rows

      REAL,    INTENT( IN )  :: LAT( NX,NY )  ! lat (deg) -90 <= LAT <= 90
      REAL,    INTENT( IN )  :: LON( NX,NY )  ! lon (deg) -180 <= LON <= 180

      REAL,    INTENT( OUT ) :: COSZEN( NX,NY ) !  cos of zenith angle

C Parameters:
 
      REAL, PARAMETER ::
     &            AA =   0.15,
     &            BB =   3.885,
     &            CC = - 1.253,
     &          SIGA = 279.9348,
     &          SDEC =   0.39784984,  ! SIN (23^26'37.8") the declination angle
     &           D60 = 1.0 / 60.0,
     &           D15 = 1.0 / 15.0,
     &           D24 = 1.0 / 24.0,
     &        ROTDAY = 360.0 / 365.242 ! fraction of a complete rotation per day

C Local variables:

      INTEGER    IOS, R, C
      REAL       SLA, GMT,  TK, DAD, DF,
     &           DESIN, DECOS, DESIN2, DECOS2, SIG, DECSIN, DECCOS,
     &           EQT, TST, HRANGL
                     
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      REAL, ALLOCATABLE, SAVE :: SINLAT( :,: )
      REAL, ALLOCATABLE, SAVE :: COSLAT( :,: )  

      CHARACTER( 16 ) :: PNAME = 'CZANGLE'   !  procedure name

C----------------------------------------------------------------------

C compute sine of lat and lon first time through

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         ALLOCATE( SINLAT( NX,NY  ), STAT=IOS )
         CALL CHECKMEM( IOS, 'SINLAT', PNAME )

         ALLOCATE( COSLAT( NX,NY ), STAT=IOS )
         CALL CHECKMEM( IOS, 'COSLAT', PNAME )

         DO R = 1, NY
            DO C = 1, NX
               SLA = PI180 * LAT( C,R )
               SINLAT( C,R ) = SIN( SLA )
               COSLAT( C,R ) = COS( SLA )
            END DO
         END DO

      END IF   ! if firstime

C Convert time to hours and add time-zone offset
      
      GMT    = FLOAT( JTIME / 10000 )                        !  hr part
     &       + D60 * ( FLOAT( MOD( JTIME / 100 , 100 ) )     !  min part
     &       + D60 *   FLOAT( MOD( JTIME, 100 ) ) )          !  sec part
      DAD    = GMT * D24 + MOD( JDATE, 1000 )
      DF     = ROTDAY * PI180 * DAD      !  The terrestrial-rotation angle
           
      DESIN  = SIN( DF )          !  SINE   of this angle
      DECOS  = COS( DF )          !  COSINE of this angle
           
      DESIN2 = SIN( DF + DF )     !  SINE   of twice the angle
      DECOS2 = COS( DF + DF )     !  COSINE of twice the angle
           
      SIG  =  DF
     &     +  PI180 * ( SIGA
     &                + 1.914827 * DESIN  - 0.079525 * DECOS
     &                + 0.019938 * DESIN2 - 0.00162  * DECOS2 )
           
C The sine and cosine of the declination
      
      DECSIN = SDEC * SIN( SIG )
      DECCOS = SQRT( 1.0 - DECSIN * DECSIN )
           
C The equation of time adjustment
      
      EQT = 0.123470 * DESIN  - 0.004289 * DECOS
     &    + 0.153809 * DESIN2 + 0.060783 * DECOS2

      DO R = 1, NY
         DO C = 1, NX

            TK     =  GMT + LON( C,R ) * D15   !  Distance in hours from LON=0
            TST    =  TK - EQT                 !  true solar time
            HRANGL =  PI180 * 15.0 * ABS( TST - 12.0 ) !  hour angle
                   
C Compute the sine of the solar elevation (cosine of zenith angle)
          
            COSZEN( C,R ) = DECSIN * SINLAT( C,R )
     &                    + DECCOS * COSLAT( C,R ) * COS( HRANGL )
                   
         END DO
      END DO

      RETURN
      END
