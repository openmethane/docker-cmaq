
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/couple/gencoor/couple.F,v 1.3 2011/07/26 15:32:24 ppercell Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)couple.F	1.1 /project/mod3/CMAQ/src/convert/couple/SCCS/s.couple.F 03 Jun 1997 11:41:33

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COUPLE ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Convert units and couple concentration values in CGRID for transport
 
C Preconditions:
 
C Subroutines and functions called:
C   INTERPX, M3EXIT
 
C Revision History:
C    Jeff Sep 97 - leave gas chem, non-reactive and tracer species in
C                  standard (ppmV) units for transport
C    2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   Jeff - Dec 00 - move CGRID_MAP into f90 module
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX
C        - Jun 01 - update units conversion calls and comments
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   Feb 2011 by Peter Percell at UH-IMAQS:
C        - duplicated first time block in DECOUPLE to support adjoint/4dvar
C   Jul 2011 by Peter Percell at UH-IMAQS:
C        - for backward runs of CMAQ/ADJOINT, the function INTERPX needs to
C          be told to use a backward interpolation interval for a first time
C          call at the last date & time in the MET_CRO_3D file
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines

      IMPLICIT NONE   

C Include files:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

!     REAL      :: CGRID( :,:,:,: )  ! concentrations
      REAL, POINTER :: CGRID( :,:,:,: )  ! concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
C Parameters:

C External Functions (not already declared by IODECL3.EXT):

C File Variables:
 
      REAL        JACOBM( NCOLS,NROWS,NLAYS )  !"total" Jacobian
      REAL        RHOJ  ( NCOLS,NROWS,NLAYS )  !"total" Jacobian * air density
      REAL :: RHOJ_CHK( NCOLS,NROWS,NLAYS, 1 ) ! RhoJ (Kg/m**2) from ADJ_VA_RHOJ_CHK file

C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'COUPLE'
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 96 ) :: XMSG = ' '

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      INTEGER, SAVE :: NQAE              ! number of micro-grams / m**3 species
      INTEGER, SAVE :: QAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      INTEGER, SAVE :: NNAE              ! number of  # / m**3 species
      INTEGER, SAVE :: NAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      INTEGER, SAVE :: NSAE              ! number of  m**2 / m**3 species
      INTEGER, SAVE :: SAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""

!     REAL        SGRID( NCOLS,NROWS,NLAYS,N_AE_SPCD ) ! aerosol subsets
      REAL, ALLOCATABLE :: SGRID( :,:,:,: ) ! aerosol subsets

      INTEGER     ALLOCSTAT
      INTEGER     GXOFF, GYOFF          ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      INTEGER, SAVE :: ASPC     ! pointer in CGRID to transported RHOJ
      INTEGER     NSPCS            ! number of species for subroutine calls
      INTEGER     OFF              ! loop offset to CGRID species
      INTEGER     C, R, L, SPC     ! loop counters

      INTEGER     LDATE, LTIME

      INTERFACE
         SUBROUTINE MGPM3_KGPM3 ( NSPCS, NLAYS, SGRID )
            IMPLICIT NONE
            INTEGER, INTENT( IN )  :: NSPCS, NLAYS
            REAL, INTENT( INOUT )  :: SGRID( :,:,:,: )
         END SUBROUTINE MGPM3_KGPM3
         SUBROUTINE KGPM3_MGPM3 ( NSPCS, NLAYS, SGRID )
            IMPLICIT NONE
            INTEGER, INTENT( IN )  :: NSPCS, NLAYS
            REAL, INTENT( INOUT )  :: SGRID( :,:,:,: )
         END SUBROUTINE KGPM3_MGPM3
      END INTERFACE
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
         ASPC = GC_STRT - 1 + N_GC_SPCD ! SPT

C create aerosol species pointers to distinguish micro-grams / m**3,
C # / m**3 (number density), and m**2 / m**3 (surface area) species
 
         NQAE = 0       ! no. of micro-grams / m**3 species
         NNAE = 0       ! no. of  # / m**3 species
         NSAE = 0       ! no. of  m**2 / m**3 species
         OFF = AE_STRT - 1
         DO SPC = 1, N_AE_SPC
            IF ( AE_SPC( SPC )( 1:3 ) .EQ. 'NUM' ) THEN
               NNAE = NNAE + 1
               NAE( NNAE ) = OFF + SPC
               ELSE IF ( AE_SPC( SPC )( 1:3 ) .EQ. 'SRF' ) THEN
               NSAE = NSAE + 1
               SAE( NSAE ) = OFF + SPC
               ELSE
               NQAE = NQAE + 1
               QAE( NQAE ) = OFF + SPC
               END IF
            END DO

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

C Open checkpoint file as read-only
         CALL CHK_FILE_OPEN ( ADJ_CPL_CHK, JDATE, JTIME, .FALSE. )

         END IF       ! if firstime

C Read Jacobian X Air Density (Jacobian =  sq. root det. metric tensor)

      VNAME = 'DENSA_J'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME,
     &                    RHOJ ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
 
      VNAME = 'JACOBM'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME,
     &                    JACOBM ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      DO L = 1, NLAYS
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               !RHOJ( C,R,L ) = CGRID( C,R,L,ASPC )
               CGRID( C,R,L,ASPC ) = RHOJ( C,R,L )
               END DO
            END DO
         END DO

C couple for advection - use density times the square root of the determinant
C of the metric tensor = RHOJ
 
C (air density X "total" jacobian) X mixing ratio [ppmV]
      NSPCS = N_GC_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = GC_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+SPC ) = CGRID( C,R,L,OFF+SPC )
     &                                      * RHOJ( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
C micro-grams/m**3 -> Kg/m**3
      NSPCS = NQAE
      IF ( NSPCS .GT. 0 ) THEN
         ALLOCATE ( SGRID( MY_NCOLS,MY_NROWS,NLAYS,NSPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating SGRID'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     SGRID( C,R,L,SPC ) = CGRID( C,R,L,QAE( SPC ) )
                     END DO
                  END DO
               END DO
            END DO
         CALL MGPM3_KGPM3 ( NSPCS, NLAYS, SGRID )
C ("total" jacobian) X [Kg/m**3]
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,QAE( SPC ) ) = SGRID( C,R,L,SPC )
     &                                         * JACOBM( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         DEALLOCATE ( SGRID )
         END IF

C CGRID in #/m**3 -> ("total" jacobian) X [#/m**3]
      NSPCS = NNAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,NAE( SPC ) ) = CGRID( C,R,L,NAE( SPC ) )
     &                                         * JACOBM( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
C CGRID in m**2/m**3 -> ("total" jacobian) X [m**2/m**3]
      NSPCS = NSAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,SAE( SPC ) ) = CGRID( C,R,L,SAE( SPC ) )
     &                                         * JACOBM( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C (air density X "total" jacobian) X mixing ratio [ppmV]
      NSPCS = N_NR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = NR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+SPC ) = CGRID( C,R,L,OFF+SPC )
     &                                      * RHOJ( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
C (air density X "total" jacobian) X mixing ratio [ppmV]
      NSPCS = N_TR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = TR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+SPC ) = CGRID( C,R,L,OFF+SPC )
     &                                      * RHOJ( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      RETURN
 
C............................................................................
C entry DECOUPLE
 
      ENTRY DECOUPLE_BWD ( CGRID, JDATE, JTIME, TSTEP )

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
         ASPC = GC_STRT - 1 + N_GC_SPCD ! SPT

C create aerosol species pointers to distinguish micro-grams / m**3,
C # / m**3 (number density), and m**2 / m**3 (surface area) species
 
         NQAE = 0       ! no. of micro-grams / m**3 species
         NNAE = 0       ! no. of  # / m**3 species
         NSAE = 0       ! no. of  m**2 / m**3 species
         OFF = AE_STRT - 1
         DO SPC = 1, N_AE_SPC
            IF ( AE_SPC( SPC )( 1:3 ) .EQ. 'NUM' ) THEN
               NNAE = NNAE + 1
               NAE( NNAE ) = OFF + SPC
               ELSE IF ( AE_SPC( SPC )( 1:3 ) .EQ. 'SRF' ) THEN
               NSAE = NSAE + 1
               SAE( NSAE ) = OFF + SPC
               ELSE
               NQAE = NQAE + 1
               QAE( NQAE ) = OFF + SPC
               END IF
            END DO

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

C For backward runs of CMAQ/ADJOINT, the function INTERPX needs to be told
C to use a backward interpolation interval for a first time call at the last
C date & time in the MET_CRO_3D file
         CALL LSTEPF( MET_CRO_3D, LDATE, LTIME )
         IF ( JDATE .EQ. LDATE .AND. JTIME .EQ. LTIME ) THEN
            CALL NEXTIME(LDATE, LTIME, -TSTEP(2))

            VNAME = 'DENSA_J'
            IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &           STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &           LDATE, LTIME,
     &           RHOJ ) ) THEN
               XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
               CALL M3EXIT( PNAME, LDATE, LTIME, XMSG, XSTAT1 )
               END IF
         
            VNAME = 'JACOBM'
            IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &           STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &           LDATE, LTIME,
     &           JACOBM ) ) THEN
               XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
               CALL M3EXIT( PNAME, LDATE, LTIME, XMSG, XSTAT1 )
               END IF
            END IF
         
C Open checkpoint file as read-only
         CALL CHK_FILE_OPEN ( ADJ_CPL_CHK, JDATE, JTIME, .FALSE. )

         END IF       ! if firstime

C Read advected RHOJ checkpoint at sync step
      CALL CHK_FILE_READ ( ADJ_CPL_CHK, JDATE, JTIME, RHOJ_CHK )

      VNAME = 'DENSA_J'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME,
     &                    RHOJ ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
 
      VNAME = 'JACOBM'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME,
     &                    JACOBM ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      DO L = 1, NLAYS
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               RHOJ( C,R,L ) = RHOJ_CHK( C,R,L,1 )
               END DO
            END DO
         END DO
 
C decouple for chemistry and diffusion
 
C mixing ratio [ppmV] / (air density X "total" jacobian) 
      NSPCS = N_GC_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = GC_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+SPC ) = CGRID( C,R,L,OFF+SPC )
     &                                      / RHOJ( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
C {J}Kg/m**3 -> {J}micro-grams/m**3
      NSPCS = NQAE
      IF ( NSPCS .GT. 0 ) THEN
         ALLOCATE ( SGRID( MY_NCOLS,MY_NROWS,NLAYS,NSPCS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating SGRID'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     SGRID( C,R,L,SPC ) = CGRID( C,R,L,QAE( SPC ) )
                     END DO
                  END DO
               END DO
            END DO
         CALL KGPM3_MGPM3 ( NSPCS, NLAYS, SGRID )
C [micro-grams/m**3] / ("total" jacobian)
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,QAE( SPC ) ) = SGRID( C,R,L,SPC )
     &                                         / JACOBM( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         DEALLOCATE ( SGRID )
         END IF

C CGRID in {J}[#/m**3] -> #/m**3
      NSPCS = NNAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,NAE( SPC ) ) = CGRID( C,R,L,NAE( SPC ) )
     &                                         / JACOBM( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
C CGRID in {J}[m**2/m**3] -> m**2/m**3
      NSPCS = NSAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,SAE( SPC ) ) = CGRID( C,R,L,SAE( SPC ) )
     &                                         / JACOBM( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C mixing ratio [ppmV] / (air density X "total" jacobian) 
      NSPCS = N_NR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = NR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+SPC ) = CGRID( C,R,L,OFF+SPC )
     &                                      / RHOJ( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
C mixing ratio [ppmV] / (air density X "total" jacobian) 
      NSPCS = N_TR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = TR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+SPC ) = CGRID( C,R,L,OFF+SPC )
     &                                      / RHOJ( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      RETURN
      END
