
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_ddm3d/AVG_CONC.F,v 1.1.1.1 2010/06/14 16:03:04 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE AVG_CONC

      IMPLICIT NONE

C Function: species, layer pointers and definitions for integral average
C CONC calculations

      INTEGER, SAVE :: N_ASPCS    ! Number of species saved to avg conc file
      CHARACTER( 16 ), SAVE :: AVG_SPCS( 900 ) ! avg conc file species list
      INTEGER, SAVE :: ACONC_BLEV ! Beginning level saved to avg conc file
      INTEGER, SAVE :: ACONC_ELEV ! Ending level saved to avg conc file

      INTEGER, SAVE, ALLOCATABLE :: AVG_CONC_MAP( : ) ! pointer into CGRID
      CHARACTER( 16 ), SAVE, ALLOCATABLE :: A_GC_SPC( : ) ! pointer into GC_SPC
      CHARACTER( 16 ), SAVE, ALLOCATABLE :: A_AE_SPC( : ) ! pointer into AE_SPC
      CHARACTER( 16 ), SAVE, ALLOCATABLE :: A_NR_SPC( : ) ! pointer into NR_SPC
      CHARACTER( 16 ), SAVE, ALLOCATABLE :: A_TR_SPC( : ) ! pointer into TR_SPC

C species classes configuration for average CONC 
      INTEGER, SAVE :: A_GC_STRT
      INTEGER, SAVE :: N_A_GC_SPC
      INTEGER, SAVE :: A_AE_STRT
      INTEGER, SAVE :: N_A_AE_SPC
      INTEGER, SAVE :: A_NR_STRT
      INTEGER, SAVE :: N_A_NR_SPC
      INTEGER, SAVE :: A_TR_STRT
      INTEGER, SAVE :: N_A_TR_SPC

      INTEGER, SAVE :: N_A_AE_SPCD

      CONTAINS

         SUBROUTINE A_CONC_DEFN ()

         USE VGRD_DEFN             ! vertical layer specifications

!        INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
         INCLUDE SUBST_GC_SPC      ! gas chemistry species table
         INCLUDE SUBST_AE_SPC      ! aerosol species table
         INCLUDE SUBST_NR_SPC      ! non-reactive species table
         INCLUDE SUBST_TR_SPC      ! tracer species table
         INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

         CHARACTER( 16 ) :: PNAME = 'A_CONC_DEFN'
         CHARACTER( 96 ) :: XMSG = ' '

         INTEGER OFF, VAR, V, NV, LVL
         INTEGER ALLOCSTAT, STATUS
         INTEGER :: JDATE = 0
         INTEGER :: JTIME = 0

         INTEGER, ALLOCATABLE :: GC_MAP( : )
         INTEGER, ALLOCATABLE :: AE_MAP( : )
         INTEGER, ALLOCATABLE :: NR_MAP( : )
         INTEGER, ALLOCATABLE :: TR_MAP( : )

         INTEGER, EXTERNAL :: INDEX1

         CHARACTER( 32 ) :: AVG_CONC_SPCS   = 'AVG_CONC_SPCS'
         CHARACTER( 32 ) :: ACONC_BLEV_ELEV = 'ACONC_BLEV_ELEV'
         CHARACTER( 16 ) :: V_LIST( 2 )

!        integer logdev
!        integer, external :: init3

         INTERFACE
            SUBROUTINE GET_ENVLIST ( ENV_VAR, NVARS, VAL_LIST )
               IMPLICIT NONE
               CHARACTER( * ),  INTENT ( IN )  :: ENV_VAR
               INTEGER,         INTENT ( OUT ) :: NVARS
               CHARACTER( 16 ), INTENT ( OUT ) :: VAL_LIST( : )
            END SUBROUTINE GET_ENVLIST
         END INTERFACE

C-----------------------------------------------------------------------

!        LOGDEV = INIT3()

C RETRIEVE THE SPECIES SAVED TO INTEGRAL AVERAGE CONCENTRATION FILE

         CALL GET_ENVLIST ( AVG_CONC_SPCS, N_ASPCS, AVG_SPCS )

         IF ( AVG_SPCS( 1 ) .EQ. 'ALL' ) THEN
            N_ASPCS = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC
            DO VAR = 1, N_GC_SPC
               V = VAR
               AVG_SPCS( V ) = GC_SPC( VAR )
            END DO
            DO VAR = 1, N_AE_SPC
               V = N_GC_SPC + VAR
               AVG_SPCS( V ) = AE_SPC( VAR )
            END DO
            DO VAR = 1, N_NR_SPC
               V = N_GC_SPC + N_AE_SPC + VAR
               AVG_SPCS( V ) = NR_SPC( VAR )
            END DO
            DO VAR = 1, N_TR_SPC
               V = N_GC_SPC + N_AE_SPC + N_NR_SPC + VAR
               AVG_SPCS( V ) = TR_SPC( VAR )
            END DO
         END IF

C Retrieve the layer range used in integral average concentration file

         CALL GET_ENVLIST ( ACONC_BLEV_ELEV, NV, V_LIST )
         IF ( NV .NE. 2 ) THEN
            XMSG = 'Environment variable error for ' // ACONC_BLEV_ELEV
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF
         READ( V_LIST( 1 ), '( I4 )' ) ACONC_BLEV
         READ( V_LIST( 2 ), '( I4 )' ) ACONC_ELEV

         IF ( ACONC_BLEV .LE. 0 .OR. ACONC_ELEV .GT. NLAYS ) THEN
            WRITE( XMSG,'( "Layer range", 2I4, " invalid for this model" )' )
     &      ACONC_BLEV, ACONC_ELEV
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         IF ( ACONC_BLEV .NE. 1 ) THEN
            WRITE( XMSG,'( "Layer", I3, " Not 1st layer in CGRID" )' )
     &      ACONC_BLEV
            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
            END IF

C Check species names against include files and create AVG_CONC_MAP, and
C get starting index in AGRID and total count for each species class

         ALLOCATE ( AVG_CONC_MAP( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'AVG_CONC_MAP memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( A_GC_SPC( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'A_GC_SPC memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( A_AE_SPC( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'A_AE_SPC memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( A_NR_SPC( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'A_NR_SPC memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( A_TR_SPC( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'A_TR_SPC memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( GC_MAP( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'GC_MAP memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( AE_MAP( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'AE_MAP memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( NR_MAP( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'NR_MAP memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

         ALLOCATE ( TR_MAP( N_ASPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'TR_MAP memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

C The selected species MUST be contiguous with the classes,
C and the classes must be in the sequence: GC, AE, NR, TR

         N_A_GC_SPC = 0
         N_A_AE_SPC = 0
         N_A_NR_SPC = 0
         N_A_TR_SPC = 0

         DO VAR = 1, N_ASPCS
            V = INDEX1 ( AVG_SPCS( VAR ), N_GC_SPC, GC_SPC )
            IF ( V .GT. 0 ) THEN
               N_A_GC_SPC = N_A_GC_SPC + 1
               A_GC_SPC( N_A_GC_SPC ) = AVG_SPCS( VAR )
               OFF = 0
               GC_MAP( N_A_GC_SPC ) = V + OFF
               ELSE
               V = INDEX1 ( AVG_SPCS( VAR ), N_AE_SPC, AE_SPC )
               IF ( V .GT. 0 ) THEN
                  N_A_AE_SPC = N_A_AE_SPC + 1
                  OFF = N_GC_SPC + 1   ! accounts for advected density
                  A_AE_SPC( N_A_AE_SPC ) = AVG_SPCS( VAR )
                  AE_MAP( N_A_AE_SPC ) = V + OFF
                  ELSE
                  V = INDEX1 ( AVG_SPCS( VAR ), N_NR_SPC, NR_SPC )
                  IF ( V .GT. 0 ) THEN
                     N_A_NR_SPC = N_A_NR_SPC + 1
                     OFF = N_GC_SPC + 1 + N_AE_SPC
                     A_NR_SPC( N_A_NR_SPC ) = AVG_SPCS( VAR )
                     NR_MAP( N_A_NR_SPC ) = V + OFF
                     ELSE
                     V = INDEX1 ( AVG_SPCS( VAR ), N_TR_SPC, TR_SPC )
                     IF ( V .GT. 0 ) THEN
                        N_A_TR_SPC = N_A_TR_SPC + 1
                        OFF = N_GC_SPC + 1 + N_AE_SPC + N_NR_SPC
                        A_TR_SPC( N_A_TR_SPC ) = AVG_SPCS( VAR )
                        TR_MAP( N_A_TR_SPC ) = V + OFF
                        ELSE
                        XMSG = 'Variable ' // AVG_SPCS( VAR )
     &                  // ' incorrect for this model'
                        CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
                        END IF
                     END IF
                  END IF
               END IF
            END DO

         A_GC_STRT = 1
         A_AE_STRT = A_GC_STRT + N_A_GC_SPC
         A_NR_STRT = A_AE_STRT + N_A_AE_SPC
         A_TR_STRT = A_NR_STRT + N_A_NR_SPC

         IF ( N_A_AE_SPC .GT. 0 ) THEN
            N_A_AE_SPCD = N_A_AE_SPC
            ELSE
            N_A_AE_SPCD = 1
            END IF

         VAR = 0
         DO V = 1, N_A_GC_SPC
            VAR = VAR + 1
            AVG_CONC_MAP( VAR ) = GC_MAP( V )
!        write( logdev,* ) '<<>>var, AVG_CONC_MAP: ', VAR, AVG_CONC_MAP( VAR )
            END DO

         DO V = 1, N_A_AE_SPC
            VAR = VAR + 1
            AVG_CONC_MAP( VAR ) = AE_MAP( V )
!        write( logdev,* ) '<<>>var, AVG_CONC_MAP: ', VAR, AVG_CONC_MAP( VAR )
            END DO

         DO V = 1, N_A_NR_SPC
            VAR = VAR + 1
            AVG_CONC_MAP( VAR ) = NR_MAP( V )
!        write( logdev,* ) '<<>>var, AVG_CONC_MAP: ', VAR, AVG_CONC_MAP( VAR )
            END DO

         DO V = 1, N_A_TR_SPC
            VAR = VAR + 1
            AVG_CONC_MAP( VAR ) = TR_MAP( V )
!        write( logdev,* ) '<<>>var, AVG_CONC_MAP: ', VAR, AVG_CONC_MAP( VAR )
            END DO

         IF ( N_A_GC_SPC .EQ. 0 ) THEN
            DEALLOCATE ( A_GC_SPC )
            DEALLOCATE ( GC_MAP )
            END IF
         IF ( N_A_AE_SPC .EQ. 0 ) THEN
            DEALLOCATE ( A_AE_SPC )
            DEALLOCATE ( AE_MAP )
            END IF
         IF ( N_A_NR_SPC .EQ. 0 ) THEN
            DEALLOCATE ( A_NR_SPC )
            DEALLOCATE ( NR_MAP )
            END IF
         IF ( N_A_TR_SPC .EQ. 0 ) THEN
            DEALLOCATE ( A_TR_SPC )
            DEALLOCATE ( TR_MAP )
            END IF

         END SUBROUTINE A_CONC_DEFN

      END MODULE AVG_CONC
