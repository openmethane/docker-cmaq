
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_txhg/rddepv.F,v 1.1.1.1 2010/06/14 16:03:10 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE RDDEPV ( MDATE, MTIME, MSTEP, CGRID, DEPV )

C-----------------------------------------------------------------------
C Function:
C   Get deposition velocities for vertical diffusion
 
C Preconditions:
 
C Subroutines and functions called:
C INTERP3, M3EXIT, AERO_DEPV
 
C Revision History:
C   Jeff - Sep 97 - targetted, call subroutine for aerosol species
C   Jeff - Dec 97 - add MDEPV
C   2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   Jeff - May 2000 - correct DEPV to "real world" values
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX
C   15 Mar 03 J.Young: swap depv and mdepv dimensions
C   21 Jul 03 J.Gipson: revised to compute avg DEPV for aero3 NR CG/SOA species
C   22 Aug 03 J.Young: aero in ppm units, eliminate MDEPV
C   23 Oct 03 J.Young: change depv from flux to m/s
C   14 Sep 04 J.Gipson: generalized avg DEPV calcs for aero3 SGTOT species
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   07 Jun 05 P.Bhave: added logic for 'AE4' mechanism
C   03 Apr 07 D.Wong: added logic to avoid divide by zero in calculating
C                     the fraction of SOA in j-mode
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE

!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_DEPV   ! gas chem dep vel surrogate names and map table
      INCLUDE SUBST_AE_DEPV   ! aerosol dep vel surrogate names and map table
      INCLUDE SUBST_NR_SPC
      INCLUDE SUBST_NR_DEPV   ! non-react dep vel surrogate names and map table
      INCLUDE SUBST_NR_N2AE
      INCLUDE SUBST_TR_DEPV   ! tracer dep vel surrogate names and map table
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_GC_SPC
      INCLUDE SUBST_AE_SPC
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block

C Arguments:

      INTEGER      MDATE, MTIME, MSTEP ! internal simulation date&time, tstep
      REAL, POINTER :: CGRID( :,:,:,: ) !  concentrations
      REAL      :: DEPV ( :,:,: )      ! deposition velocities

C Paramters:

      INTEGER, PARAMETER :: N_SGTOT = 10   ! No. of SOA+gas species in aero3

      REAL, PARAMETER :: RGAS1 = 1.0 / 8.314510 ! reciprocal of universal
                                                ! gas constant

      REAL, PARAMETER :: H_VAP = 156.0E+03   ! enthalpy of vaporization (J/mol)
      REAL, PARAMETER :: TREF = 298.0        ! Reference temperature for Cstar
      REAL, PARAMETER :: TREF1 = 1.0 / 298.0 ! Reciprocal ref temperature


C External Functions:
      INTEGER INDEX1                       !  Looks up name in a list
      INTEGER FINDEX                       !  Looks up number in a list

C File variables:

!     REAL         ADEPV( NCOLS,NROWS,N_AE_DEPVD ) ! Aerosol dep vel from
                                                   ! squared
      REAL, ALLOCATABLE, SAVE :: TEMP( :,: )       ! air temp (degK)
      REAL, ALLOCATABLE, SAVE :: PRES( :,: )       ! air pressure (Pa)
      REAL, ALLOCATABLE, SAVE :: ADEPV( :,:,: )    ! Aerosol dep vel from
                                                   ! subroutine AERO_DEPV
      REAL, ALLOCATABLE, SAVE :: BUFF( :,: )       ! read buffer
      REAL, ALLOCATABLE, SAVE :: BUFF2( :,:,: )    ! read buffer

C Local variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 16 ), SAVE :: PNAME = 'RDDEPV'
      CHARACTER( 96 ) :: XMSG = ' '

!     Names of aero3 SOA + gas species
      CHARACTER( 16 ) :: SGNAME( N_SGTOT ) = ( / 'SGTOT_ALK       ',
     &                                           'SGTOT_OLI_1     ', 
     &                                           'SGTOT_OLI_2     ', 
     &                                           'SGTOT_XYL_1     ', 
     &                                           'SGTOT_XYL_2     ', 
     &                                           'SGTOT_CSL       ', 
     &                                           'SGTOT_TOL_1     ',
     &                                           'SGTOT_TOL_2     ', 
     &                                           'SGTOT_TRP_1     ', 
     &                                           'SGTOT_TRP_2     '  / )

      INTEGER      ALLOCSTAT
      INTEGER      C, R, V, N                      ! loop induction variables
      INTEGER      STRT, FINI                      ! loop induction variables

      INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      LOGICAL, SAVE :: LAERO3 = .FALSE.     ! Flag for aero3
      LOGICAL, SAVE :: LAERO4 = .FALSE.     ! Flag for aero4
      INTEGER, SAVE :: PRJ, PRI             ! Primary organic AE indices
      INTEGER, SAVE :: ASJ, ASI             ! Anthropogenic SOA indices 
      INTEGER, SAVE :: BSJ, BSI             ! Biogenic SOA indices 
      INTEGER          SP                   ! Cgrid species index
      INTEGER          AE_OFFSET            ! Offset to AE start in CGRID
      INTEGER          IND                  ! Temp index
      INTEGER          IND2                 ! Temp index

      REAL   POA            ! Total primary organic AE
      REAL   ASOA           ! Total anthropogenic SOA in AE phase
      REAL   BSOA           ! Total biogenic SOA in AE phase
      REAL   JFRAC          ! Fraction of total SOA in j-mode
      REAL   AFRAC          ! Fraction of NR species in AE phase
      REAL   TORG           ! Total moles of organic absorbant in AE phase
      REAL   CTOT           ! Species sum of SOA and condesable gas (ug/m3)
      REAL   C0             ! Species saturation concentration (ug/m3) at T,P
      REAL   CAER           ! Species aerosol concentration (ug/m3)
      REAL   TEMP1          ! Reciprocal of temperature
      REAL   VP_TADJ        ! Temp adjustment factor for Cstar


      INTEGER, SAVE :: NR_OFFSET            ! Offset to NR start in CGRID

      REAL, SAVE  ::   MW_POA  = 220.0      ! MW of primary organic AE
      REAL, SAVE  ::   MW_ASOA = 150.0      ! MW of anthropogenic SOA
      REAL, SAVE  ::   MW_BSOA = 177.0      ! MW of biogenic SOA


      CHARACTER( 16 ), SAVE :: GASDV4SG( N_SGTOT ) ! list of DV surrogates 
                                                   ! for SGTOT

      INTEGER, SAVE :: N_GASDV4SG          ! No. of DV surrogates for SGTOT spec
      INTEGER, SAVE :: SG2NRSPC( N_SGTOT ) ! Index of SGTOT spec in NR_SPC list
      INTEGER, SAVE :: SG2NRDV(  N_SGTOT ) ! Index of SGTOT spec in NR_DEPV list
      INTEGER, SAVE :: SG2GASDV( N_SGTOT ) ! Index of SGTOT spec in DV surr list

      LOGICAL, SAVE :: SGSPEC( N_NR_DEPV ) ! Flag for SGTOT spec in NR_DEPV list

      REAL, SAVE    ::  PPM2UG(  N_SGTOT ) ! species ppm to ug/m3 conv facs
      REAL, SAVE    ::  CSTAR(   N_SGTOT ) ! species saturation concs at 298K

      

      INTERFACE
         SUBROUTINE AERO_DEPV ( CGRID, MDATE, MTIME, MSTEP, ADEPV )
            IMPLICIT NONE
!           REAL, INTENT( IN )     :: CGRID( :,:,:,: )
            REAL, POINTER          :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )  :: MDATE, MTIME, MSTEP
            REAL                   :: ADEPV( :,:,: )
         END SUBROUTINE AERO_DEPV
      END INTERFACE

C-----------------------------------------------------------------------
      DATA CSTAR   / 0.3103, 111.11, 1000.0, 2.165, 64.946,
     &               0.2611,  1.716, 47.855, 0.865, 11.804 /


      IF ( FIRSTIME )  THEN

         FIRSTIME  =  .FALSE.

         IF( INDEX( MECHNAME, 'AE3' ) .NE. 0 ) LAERO3 = .TRUE.
         IF( INDEX( MECHNAME, 'AE4' ) .NE. 0 ) LAERO4 = .TRUE.

         ALLOCATE ( BUFF( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating BUFF'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         IF ( N_AE_DEPV .GT. 0 ) THEN
            ALLOCATE ( ADEPV( N_AE_DEPV,MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating ADEPV'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
            END IF

         CALL SUBHFILE ( MET_CRO_2D, GXOFF, GYOFF,
     &                   STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2 )
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )


         SGSPEC( : ) = .FALSE.

         IF( LAERO3 .OR. LAERO4 ) THEN

            AE_OFFSET = N_GC_SPC + 1            ! Offset to 1st AE spec in cgrid
            NR_OFFSET = N_GC_SPC + 1 + N_AE_SPC ! Offset to 1st NR spec in cgrid

c..check that all SOA+gas totals are in list of N2AE species &
c..get NR_SPC indices
            SG2NRSPC( : ) = 0      ! Index of STOT species in NR_SPC list
            DO N = 1, N_SGTOT

               IND = INDEX1( SGNAME( N ), N_NR_N2AE, NR_N2AE )

               IF( IND .NE. 0 ) THEN
                  SG2NRSPC( N ) = NR_N2AE_MAP( IND )
                  ELSE
                  XMSG = 'Could not find ' // SGNAME( N ) // 
     &               'in non-reactives table'
                  CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT3 ) 
                  END IF

               END DO

c..map the SGTOT species to the NR_DEPV species & get a list of gas DVs 
c..that are used by them
            N_GASDV4SG    = 0     ! No. of DV surrogates needed for SGTOT spec
            SG2NRDV(  : ) = 0     ! Index of SGTOT species in NR_DEPV list
            SG2GASDV( : ) = 0     ! Index of DV surrogate to use for SGTOT spec
            GASDV4SG( : ) = ' '   ! List of DV surrogates needed for SGTOT spec
            DO N = 1, N_SGTOT

               IND = FINDEX( SG2NRSPC( N ), N_NR_DEPV, NR_DEPV_MAP ) 

               IF( IND .NE. 0 ) THEN

                  SG2NRDV( N ) = IND
                  IND2 = INDEX1( NR_DEPV( IND ), N_GASDV4SG, GASDV4SG )

                  IF( IND2 .EQ. 0 ) THEN
                     N_GASDV4SG = N_GASDV4SG + 1
                     GASDV4SG( N_GASDV4SG ) = NR_DEPV( IND )
                     SG2GASDV( N ) = N_GASDV4SG
                     ELSE
                     SG2GASDV( N ) = IND2
                     END IF

                  END IF
  
               END DO

c..Flag NR_DEPV species that are also SGTOT species
            DO N = 1, N_SGTOT
               IF( SG2NRDV( N ) .NE. 0 ) SGSPEC( SG2NRDV( N ) ) = .TRUE.
            ENDDO

                    
c..get indices of aero3 AE species in cgrid & set ppm to ug/m3 conversion facs
            PRJ = AE_OFFSET + INDEX1( 'AORGPAJ', N_AE_SPC, AE_SPC )
            PRI = AE_OFFSET + INDEX1( 'AORGPAI', N_AE_SPC, AE_SPC )
            ASJ = AE_OFFSET + INDEX1( 'AORGAJ',  N_AE_SPC, AE_SPC )
            ASI = AE_OFFSET + INDEX1( 'AORGAI',  N_AE_SPC, AE_SPC )
            BSJ = AE_OFFSET + INDEX1( 'AORGBJ',  N_AE_SPC, AE_SPC )
            BSI = AE_OFFSET + INDEX1( 'AORGBI',  N_AE_SPC, AE_SPC )

            DO N = 1, N_SGTOT
               PPM2UG( N ) = NR_MOLWT( SG2NRSPC( N ) ) * RGAS1
               END DO

c..Allocate array for temp, press, & DV surrogates used for SGTOT species
            ALLOCATE ( TEMP( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TEMP'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
            
            ALLOCATE ( PRES( MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating TEMP'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF

            ALLOCATE ( BUFF2( MY_NCOLS,MY_NROWS,N_GASDV4SG ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating BUFF2'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF

            END IF       ! if (laero3 .or. laero4)

         END IF          !  if firstime

      N = 0
      STRT = 1
      FINI = N_GC_DEPV
      DO V = STRT, FINI
         N = N + 1

C                                        vvvvvvv - dep vel surrogate name
         IF ( .NOT. INTERPX( MET_CRO_2D, GC_DEPV( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
!    &                       MDATE, MTIME, DEPV( 1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
            XMSG = 'Could not read  ' // GC_DEPV( N ) // ' from ' // MET_CRO_2D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DEPV( V,C,R ) = GC_DEPV_FAC( N ) * BUFF( C,R ) 
               END DO
            END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         WRITE( *,* ) '  GC Depv: ', V, '  ', GC_DEPV( N )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO

      IF ( N_AE_DEPV .GT. 0 ) THEN
         CALL AERO_DEPV( CGRID, MDATE, MTIME, MSTEP, ADEPV )

         STRT = N_GC_DEPV + 1
         FINI = N_GC_DEPV + N_AE_DEPV

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               N = 0
               DO V = STRT, FINI
                  N = N + 1
                  DEPV( V,C,R ) = AE_DEPV_FAC( N ) * ADEPV( N,C,R ) 
                  END DO
               END DO
            END DO

         END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         WRITE( *,* ) '  AE Depv: ', V, '  ', AE_DEPV( N )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      N = 0
      STRT = N_GC_DEPV + N_AE_DEPV + 1
      FINI = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV
      DO V = STRT, FINI
         N = N + 1

         IF( (LAERO3 .OR. LAERO4) .AND. SGSPEC( N ) ) CYCLE

C                                        vvvvvvv - dep vel surrogate name
         IF ( .NOT. INTERPX( MET_CRO_2D, NR_DEPV( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
!    &                       MDATE, MTIME, DEPV( 1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
            XMSG = 'Could not read  ' // NR_DEPV( N ) // ' from ' // MET_CRO_2D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DEPV( V,C,R ) = NR_DEPV_FAC( N ) * BUFF( C,R )
               END DO
            END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         WRITE( *,* ) '  NR Depv: ', V, '  ', NR_DEPV( N )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO

      N = 0
      STRT = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV + 1
      FINI = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV + N_TR_DEPV 
      DO V = STRT, FINI
         N = N + 1

C                                         vvvvvvv - dep vel surrogate name
         IF ( .NOT. INTERPX( DEPV_TRAC_1, TR_DEPV( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
!    &                       MDATE, MTIME, DEPV( 1,1,V ) ) ) THEN
     &                       MDATE, MTIME, BUFF ) ) THEN
            XMSG = 'Could not read  ' // TR_DEPV( N ) // ' from ' // DEPV_TRAC_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DEPV( V,C,R ) = TR_DEPV_FAC( N ) * BUFF( C,R ) 
               END DO
            END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         WRITE( *,* ) '  TR Depv: ', V, '  ', TR_DEPV( N )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO


       IF( .NOT. ( LAERO3 .OR. LAERO4 ) ) RETURN

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This section computes deposition velocities for the NR secondary organic
c  aerosol species included in aero3. These model species are the sum of the 
c  SOA precursor condensable gas (CG) and SOA species itself.
c  A weighted average deposition velocity is calculated since these model
c  species are the sum of species in two phases.
c  For the gas phase, the deposition velocity for organic acids is used,
c  and for SOA, the DVs for the i- and j-mode VMASS are used.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c.. get the temp and pressure

      IF ( .NOT. INTERPX( MET_CRO_3D, 'TA', PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1, 1,
     &                    MDATE, MTIME, TEMP ) ) THEN
         XMSG = 'Could not interpolate TA from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

      IF ( .NOT. INTERPX( MET_CRO_3D, 'PRES', PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1, 1,
     &                    MDATE, MTIME, PRES ) ) THEN
         XMSG = 'Could not interpolate PRES from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

c..read the gas-phase depvs to use for SGTOT
      DO N = 1, N_GASDV4SG
         IF ( .NOT. INTERPX( MET_CRO_2D, GASDV4SG( N ), PNAME,
     &                       STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                       MDATE, MTIME, BUFF2( 1,1,N ) ) ) THEN
            XMSG = 'Could not read  ' // GASDV4SG( N ) // ' from ' // MET_CRO_2D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
      END DO

c..set the loop indices to do all NR species except first (NH3)

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS

            POA  = CGRID( C,R,1,PRJ )  + CGRID( C,R,1,PRI )
            ASOA = CGRID( C,R,1,ASJ )  + CGRID( C,R,1,ASI )
            BSOA = CGRID( C,R,1,BSJ )  + CGRID( C,R,1,BSI )

c..fraction of SOA in j-mode
            IF ( ( ASOA .EQ. 0.0 ) .AND. ( BSOA .EQ. 0.0 ) ) THEN
               JFRAC = 0.0
               ELSE
               JFRAC = ( CGRID( C, R, 1, ASJ ) + CGRID( C, R, 1, BSJ ) )
     &               / ( ASOA + BSOA )
               END IF

c..total amount of organic absorbant
            TORG = POA / MW_POA + ASOA / MW_ASOA + BSOA / MW_BSOA  

            TEMP1 = 1.0 / TEMP( C, R )
            VP_TADJ = TREF * TEMP1 * EXP( H_VAP * RGAS1 * ( TREF1 - TEMP1 ) )

            DO N = 1, N_SGTOT

               IF ( SG2NRDV( N ) .EQ. 0 ) CYCLE   ! Skip if not in NR_DEPV list


! Indices of sgtot pos in cgrid & sgtot pos in DEPV list
                  SP = NR_OFFSET + SG2NRSPC( N )  
                  V  = N_GC_DEPV + N_AE_DEPV + SG2NRDV( N )

c..get SGTOT concentration in ug/m3
                  CTOT = CGRID( C, R, 1, SP ) * PPM2UG( N ) * 
     &                   PRES( C, R ) * TEMP1

c..adjust saturation concentration for temp
                  C0 = CSTAR( N ) * VP_TADJ

c..estimate SOA fraction of SOA+CG (afrac)
c     caer = ctot * torg / ( torg + c0 / mw_cg( n ) )
c     afrac = caer / ctot
c     therefore afrac = torg / ( torg + c0 / mw_cg( n ) )

                  AFRAC = TORG / ( TORG + C0 / NR_MOLWT( SG2NRSPC( N ) ) )

c..compute weighted avg DV (for SOA protion, DV weighted by modal fraction)

                  DEPV( V,C,R ) = ( 1.0 - AFRAC ) * 
     &                  BUFF2( C,R, SG2GASDV( N ) ) + AFRAC *
     &                  ( JFRAC * ADEPV( 1,C,R ) + 
     &                  ( 1.0 - JFRAC ) * ADEPV( 2,C,R ) )

               END DO

            END DO

         END DO
   
      RETURN
      END
