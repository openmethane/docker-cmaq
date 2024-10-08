
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhgsim/PTMAP.F,v 1.1.1.1 2010/06/14 16:03:09 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE PTMAP

C-----------------------------------------------------------------------
C Function: species mappings from input point source file groups to
C           the intermediate VDEMIS_PT and the model VDEMIS arrays

C Revision History:
C     29 Nov  2007 J.Young: initial implementation
C     17 Sept 2009 B.Hutzell: corrected mapping algorithm to account
C                             for multiple usage of an emission species
C                             from an individual stack      
C-----------------------------------------------------------------------

      USE UDTYPES, ONLY: IARRY1, CARRY1 ! implementation for multiple
                                        ! groups of stack input files
      IMPLICIT NONE

      INTEGER, SAVE :: N_SPC_PTEM = 0   ! merged no. of unique species for non_PM

C Species names from input file used for point source non-PM emissions mapping
      INTEGER,          ALLOCATABLE, SAVE :: PTEM_MAP( : )
      INTEGER,          ALLOCATABLE, SAVE :: SPC_PTEM_MAP( : )
      CHARACTER( 16 ),  ALLOCATABLE, SAVE :: STKEM_NAME( : ) ! diagnostic file

C Species names from input file used for point source PM emissions mapping
      INTEGER,          ALLOCATABLE, SAVE :: PTPM_MAP( : )
      CHARACTER( 16 ),  ALLOCATABLE, SAVE :: STKPM_NAME( : ) ! diagnostic file

C Mapping for point source non-PM emissions (maps only GC, NR and TR)
      TYPE ( IARRY1 ),  ALLOCATABLE, SAVE :: STKEM_MAP( : ) ! STKEM index
      TYPE ( IARRY1 ),  ALLOCATABLE, SAVE :: GAS_MAP  ( : ) ! VDEMIS_PT index

      INTEGER, SAVE :: N_SPC_PTPM = 0   ! merged no. of unique species for PM

C Mapping for point source PM emissions (maps pt. src. PM spcs to PMEMIS_PT)
      TYPE ( IARRY1 ),  ALLOCATABLE, SAVE :: STKPM_MAP( : ) ! STKPM index
      TYPE ( IARRY1 ),  ALLOCATABLE, SAVE :: PM_MAP   ( : ) ! PMEMIS_PT index

C-----------------------------------------------------------------------

      CONTAINS
 
         FUNCTION PTMAP_INIT ( NPTGRPS ) RESULT ( SUCCESS )

         USE STK_EMIS  ! hourly point source emissions
         USE PMEM_DEFN ! defines emissions file(s) PM mapping
         USE UDTYPES, ONLY: CARRY1 ! implementation for multiple
                                   ! groups of stack input files
         IMPLICIT NONE

C Includes:
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names & map table
         INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names & map table
         INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names & map table
         INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names & map table

         INCLUDE SUBST_GC_SPC   ! gas names & map table
         INCLUDE SUBST_AE_SPC   ! aerosol names & map table
         INCLUDE SUBST_NR_SPC   ! non-react names & map table
         INCLUDE SUBST_TR_SPC   ! tracer names & map table

         INTEGER, INTENT( IN )  :: NPTGRPS
         LOGICAL SUCCESS

         INTEGER, EXTERNAL :: INDEX1, SETUP_LOGDEV

         TYPE :: STRNDX
            CHARACTER( 16 ) :: STR
            INTEGER         :: NDX
         END TYPE STRNDX

         TYPE ( STRNDX ), ALLOCATABLE :: V( : )
         TYPE ( STRNDX ), ALLOCATABLE :: VC( : )
         TYPE ( STRNDX )              :: TV

         TYPE ( CARRY1 ), ALLOCATABLE :: STKEM_VAR( : )   ! just to record on log
         TYPE ( CARRY1 ), ALLOCATABLE :: STKPM_VAR( : )   ! just to record on log

         INTEGER, PARAMETER :: N_GAS_EMIS = N_GC_EMIS + N_NR_EMIS + N_TR_EMIS

         INTEGER          :: GMAP( N_GAS_EMIS*NPTGRPS )
         INTEGER          :: GLST( N_GAS_EMIS,NPTGRPS )
         CHARACTER( 16 )  :: GVNM( N_GAS_EMIS*NPTGRPS )
         CHARACTER( 16 )  :: GSPC( N_GAS_EMIS*NPTGRPS )
         INTEGER          :: PMAP( N_EMIS_PM*NPTGRPS )
         INTEGER          :: PLST( N_EMIS_PM,NPTGRPS )
         CHARACTER( 16 )  :: PVNM( N_EMIS_PM*NPTGRPS )
         CHARACTER( 16 )  :: PSPC( N_EMIS_PM*NPTGRPS )
         CHARACTER( 16 )  :: VNAME

         CHARACTER( 16 )       :: PNAME = 'PTMAP_INIT'  ! procedure name
         INTEGER I, I2, J, K, L, N, NM, NCULL, IOS, LOGDEV
         LOGICAL          :: QUIT = .FALSE.

!        LOGICAL :: NO_PMEMIS = .FALSE.   ! set to .TRUE. if N_AE_EMIS = 0

C-----------------------------------------------------------------------

         LOGDEV = SETUP_LOGDEV()

         SUCCESS = .TRUE.

C point source emissions to gas-phase species map

         ALLOCATE ( STKEM_MAP( NPTGRPS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'STKEM_MAP', PNAME )

         ALLOCATE ( GAS_MAP( NPTGRPS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'GAS_MAP', PNAME )

         ALLOCATE ( STKEM_VAR( NPTGRPS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'STKEM_VAR', PNAME )

!        WRITE( LOGDEV,1005 )

         NM = 0; L = 0


         DO N = 1, NPTGRPS
            
            I  = 0
            LOOP_GC: DO J = 1, N_GC_EMIS
               VNAME = GC_EMIS( J )
               K = INDEX1( VNAME, NVARS_EMIS( N ), STKSPC( N )%ARRY )
               IF ( K .GT. 0 ) THEN

!WTH skip if emission species used multiple times for current point group
                   DO I2 = NM+1, L 
                      IF( VNAME .EQ. GVNM( I2 ))CYCLE LOOP_GC
                   ENDDO
D                  WRITE(LOGDEV,*)N,J,' ',VNAME,' added to GVNM array '
!WTH
                  I = I + 1
                  L = L + 1
                  GMAP( L ) = J
                  GLST( I,N ) = K
                  GVNM( L ) = VNAME
                  GSPC( L ) = GC_SPC( GC_EMIS_MAP( J ) )
D                 WRITE( LOGDEV,1009 ) N, J, K, I, L, GVNM( L ), GLST( I,N ), GMAP( L )
               END IF
            END DO LOOP_GC

            LOOP_NR: DO J = 1, N_NR_EMIS
               VNAME = NR_EMIS( J )
               K = INDEX1( VNAME, NVARS_EMIS( N ), STKSPC( N )%ARRY )
               IF ( K .GT. 0 ) THEN

!WTH skip if emission species used multiple times for current point group
                   DO I2 = NM+1, L
                        IF( VNAME .EQ. GVNM( I2 ))CYCLE LOOP_NR
                   ENDDO
D                  WRITE(LOGDEV,*)N,J,' ',VNAME,' added to GVNM array '
!WTH
                  I = I + 1
                  L = L + 1
                  GMAP( L ) = J + N_GC_EMIS + N_AE_EMIS
                  GLST( I,N ) = K
                  GVNM( L ) = VNAME
                  GSPC( L ) = NR_SPC( NR_EMIS_MAP( J ) )
D                 WRITE( LOGDEV,1009 ) N, J, K, I, L, GVNM( L ), GLST( I,N ), GMAP( L )
               END IF
            END DO LOOP_NR

            LOOP_TR: DO J = 1, N_TR_EMIS
               VNAME = TR_EMIS( J )
               K = INDEX1( VNAME, NVARS_EMIS( N ), STKSPC( N )%ARRY )
               IF ( K .GT. 0 ) THEN
!WTH skip if emission species used multiple times for current point group
                   DO I2 = NM+1, L
                     IF( VNAME .EQ. GVNM( I2 ))CYCLE LOOP_TR
                   ENDDO
D                  WRITE(LOGDEV,*)N,J,' ',VNAME,' added to GVNM array '
!WTH
                  I = I + 1
                  L = L + 1
                  GMAP( L ) = J + N_GC_EMIS + N_AE_EMIS + N_NR_EMIS
                  GLST( I,N ) = K
                  GVNM( L ) = VNAME
                  GSPC( L ) = TR_SPC( TR_EMIS_MAP( J ) )
D                 WRITE( LOGDEV,1009 ) N, J, K, I, L, GVNM( L ), GLST( I,N ), GMAP( L )
               END IF
            END DO LOOP_TR

            NM = NM + I
            
            STKEM_MAP( N )%LEN = I
            ALLOCATE ( STKEM_MAP( N )%ARRY( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKEM_MAP', PNAME )
            STKEM_MAP( N )%ARRY = 0   ! array

            STKEM_VAR( N )%LEN = I
            ALLOCATE ( STKEM_VAR( N )%ARRY( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKEM_VAR', PNAME )
            STKEM_VAR( N )%ARRY = ' '   ! array

            GAS_MAP( N )%LEN = I
            ALLOCATE ( GAS_MAP( N )%ARRY( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GAS_MAP', PNAME )
            GAS_MAP( N )%ARRY = 0   ! array

         END DO   ! NPTGRPS

D        WRITE( LOGDEV,'(5X, A, I3)' ) 'NM =', NM
D         WRITE( LOGDEV,* ) ' '
         L = 0
         DO N = 1, NPTGRPS
D            WRITE( LOGDEV,1011 ) N, STKEM_MAP( N )%LEN
            DO I = 1, STKEM_MAP( N )%LEN
               L = L + 1
               STKEM_MAP( N )%ARRY( I ) = GLST( I,N )
               STKEM_VAR( N )%ARRY( I ) = GVNM( L )
            END DO
            DO I = 1, STKEM_MAP( N )%LEN, 8
               I2 = MIN( I + 7, STKEM_MAP( N )%LEN )
D               WRITE( LOGDEV,1015 ) N, ( STKEM_MAP( N )%ARRY( J ), J = I, I2 )
D               WRITE( LOGDEV,1019 ) N, ( STKEM_VAR( N )%ARRY( J ), J = I, I2 )
            END DO
         END DO

         ALLOCATE( V( NM ), STAT=IOS )
         ALLOCATE( VC( NM ), STAT=IOS )

         DO I = 1, NM
            V( I )%STR = GVNM( I )
            V( I )%NDX = GMAP( I )
         END DO

C Sort ...
D         WRITE( LOGDEV,'(5X, A)' ) 'Concatenated Species List:'
D         WRITE( LOGDEV,'(5X, I3, 1X, A, I3)' ) 1, V( 1 )
         DO I = 2, NM
D            WRITE( LOGDEV,'(5X, I3, 1X, A, I3)' ) I, V( I )
            TV = V( I )
            DO J = I - 1, 1, -1
               IF ( V( J )%NDX .LE. TV%NDX ) GO TO 101
               V( J+1 ) = V( J )
            END DO
            J = 0
101         CONTINUE
            V( J+1 ) = TV
         END DO

D         WRITE( LOGDEV,'(5X, A)' ) 'Sorted Species List:'
D         DO I = 1, NM
D            WRITE( LOGDEV,'(I3, 5X, A, I3)' )I, V( I )
D         END DO

C Cull ...
         VC( 1 ) = V( 1 )
         NCULL = 1
         DO I = 2, NM
            IF ( V( I )%STR .NE. VC( NCULL )%STR ) THEN
               NCULL = NCULL + 1
               VC( NCULL ) = V( I )
            END IF
         END DO
D         WRITE( LOGDEV,'(5X, A, I3)' ) 'Unique Gas =', NCULL

         N_SPC_PTEM = NCULL

         ALLOCATE ( PTEM_MAP( (N_GC_EMIS + N_NR_EMIS + N_TR_EMIS) ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PTEM_MAP', PNAME )

         ALLOCATE ( SPC_PTEM_MAP( (N_GC_EMIS + N_NR_EMIS + N_TR_EMIS) ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PTEM_MAP', PNAME )

         ALLOCATE ( STKEM_NAME( N_SPC_PTEM ), STAT=IOS )
         CALL CHECKMEM( IOS, 'STKEM_NAME', PNAME )

         GSPC = 'UNDEFINED'
         DO I = 1, N_SPC_PTEM
            PTEM_MAP( I ) = VC( I )%NDX

            IF( PTEM_MAP( I ) .LE. N_GC_EMIS )THEN
                GSPC(I)=GC_SPC( GC_EMIS_MAP( PTEM_MAP( I ) )) 
            ENDIF
            IF(PTEM_MAP( I ) .GT. N_GC_EMIS + N_AE_EMIS 
     &               .AND. PTEM_MAP( I ) .LE. N_GC_EMIS + N_AE_EMIS + N_NR_EMIS )THEN
                GSPC(I) = NR_SPC( NR_EMIS_MAP( PTEM_MAP( I ) - N_GC_EMIS - N_AE_EMIS) )
            ENDIF
            IF(PTEM_MAP( I ) .GE. N_GC_EMIS + N_AE_EMIS + N_NR_EMIS  + 1)THEN
                GSPC(I) = TR_SPC( TR_EMIS_MAP( PTEM_MAP( I ) - N_GC_EMIS - N_AE_EMIS - N_NR_EMIS ))
            ENDIF

            STKEM_NAME( I ) = VC( I )%STR
D            WRITE( LOGDEV,1041 )I, STKEM_NAME( I ), PTEM_MAP( I ),GSPC(I)
         END DO
         
         
CCC WTH
         GSPC = 'UNDEFINED'
D        WRITE( LOGDEV, *)'To match the number: N_GC_EMIS + N_NR_EMIS + N_TR_EMIS:',
D    &      (N_GC_EMIS + N_NR_EMIS + N_TR_EMIS)
D        WRITE( LOGDEV, *)'Redefining GAS EMISSIONS to MODEL SPECIES MAP '
         L = 0
         DO J = 1, N_GC_EMIS
            K = INDEX1( GC_EMIS( J ), N_SPC_PTEM, STKEM_NAME )
            IF( K .GT. 0)THEN
                L = L + 1
                GSPC( J ) = GC_SPC( GC_EMIS_MAP( J ) )
                SPC_PTEM_MAP( L ) = J
                PTEM_MAP( L ) =  K
D               WRITE( LOGDEV,1045 )L, PTEM_MAP( L ), STKEM_NAME( K ),  
D    &          SPC_PTEM_MAP( L ), GSPC( J )
            ELSE
                L = L + 1
                SPC_PTEM_MAP( L ) = -1
                PTEM_MAP( L ) =  -1
                WRITE( LOGDEV,1046 )J, GC_EMIS( J ),GC_SPC( GC_EMIS_MAP( J ) )
            ENDIF
         ENDDO
         DO J = 1, N_NR_EMIS
            K = INDEX1( NR_EMIS( J ), N_SPC_PTEM, STKEM_NAME )
            IF( K .GT. 0)THEN
                L = L + 1
                GSPC( J ) = NR_SPC( NR_EMIS_MAP( J ) )
                SPC_PTEM_MAP( L ) = J + N_GC_EMIS + N_AE_EMIS
                PTEM_MAP( L ) =  K
D               WRITE( LOGDEV,1045 )L, PTEM_MAP( L ), STKEM_NAME( K ),  
D    &          SPC_PTEM_MAP( L ), GSPC( J )
            ELSE
                L = L + 1
                SPC_PTEM_MAP( L ) = -1
                PTEM_MAP( L ) =  -1
                WRITE( LOGDEV,1046 )J, NR_EMIS( J ),NR_SPC( GC_EMIS_MAP( J ) )
            ENDIF
         ENDDO
         DO J = 1, N_TR_EMIS
            K = INDEX1( TR_EMIS( J ), N_SPC_PTEM, STKEM_NAME )
            IF( K .GT. 0)THEN
                L = L + 1
                GSPC( J ) = TR_SPC( TR_EMIS_MAP( J ) )
                SPC_PTEM_MAP( L ) = J + N_GC_EMIS + N_AE_EMIS + N_NR_EMIS
                PTEM_MAP( L ) =  K
D               WRITE( LOGDEV,1045 )L, PTEM_MAP( L ), STKEM_NAME( K ),  
D    &          SPC_PTEM_MAP( L ), GSPC( J )
            ELSE
                L = L + 1
                SPC_PTEM_MAP( L ) = -1
                PTEM_MAP( L ) =  -1
                WRITE( LOGDEV,1046 )J, TR_EMIS( J ),TR_SPC( TR_EMIS_MAP( J ) )
            ENDIF
         ENDDO
                
CCC WTH         

         DO N = 1, NPTGRPS
            J = 0
            DO I = 1, N_SPC_PTEM
               K = INDEX1( STKEM_NAME( I ), STKSPC( N )%LEN, STKSPC( N )%ARRY )
               IF ( K .NE. 0 ) THEN
                  J = J + 1
                  GAS_MAP( N )%ARRY( J ) = I
D                  WRITE( LOGDEV,1031 ) N, J, GAS_MAP( N )%ARRY( J )
               END IF
            END DO
!           WRITE( LOGDEV,'(5X, A, I4)' ) 'J =', J
         END DO

         DO N = 1, NPTGRPS
            DO I = 1, STKEM_MAP( N )%LEN
D                  WRITE( LOGDEV,1031 ) N, I, GAS_MAP( N )%ARRY( I )
                  IF( GAS_MAP( N )%ARRY( I ) .LT. 1)THEN
                    WRITE(LOGDEV,*)'WARNING: GAS_MAP( ',N,' )%ARRY( ',I,' ) = 0'
                  ENDIF
            END DO
         END DO

C point source emissions to PM species map

         ALLOCATE ( STKPM_MAP( NPTGRPS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'STKPM_MAP', PNAME )

         ALLOCATE ( PM_MAP( NPTGRPS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PM_MAP', PNAME )

         ALLOCATE ( STKPM_VAR( NPTGRPS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'STKPM_VAR', PNAME )

!        IF ( NO_PMEMIS ) THEN
         IF ( N_AE_EMIS .LE. 0 ) THEN
            DO N = 1, NPTGRPS
               STKPM_MAP( N )%LEN = 0
               ALLOCATE ( STKPM_MAP( N )%ARRY( 1 ), STAT=IOS )
               CALL CHECKMEM( IOS, 'STKPM_MAP', PNAME )
               STKPM_MAP( N )%ARRY = 0   ! array
            END DO
            WRITE( LOGDEV,'(/ 5X, A/ )' ) 'Not processing point source PM emissions!'
            RETURN
         END IF
            
!        WRITE( LOGDEV,1007 )

         NM = 0; L = 0
         DO N = 1, NPTGRPS
           
            I = 0
            LOOP_PM: DO J = 1, N_EMIS_PM
               VNAME = PM_EMIS( J )
               K = INDEX1( VNAME, NVARS_EMIS( N ), STKSPC( N )%ARRY )
               IF ( K .GT. 0 ) THEN
!WTH skip if emission species used multiple times for current point group
                   DO I2 = NM+1, L
                     IF( VNAME .EQ. PVNM( I2 ))CYCLE LOOP_PM
                   ENDDO
D                  WRITE(LOGDEV,*)N,J,' ',VNAME,' added to PVNM array '
!WTH
                  I = I + 1
                  L = L + 1
                  PMAP( L ) = J
                  PLST( I,N ) = K
                  PVNM( L ) = VNAME
D                 WRITE( LOGDEV,1009 ) N, J, K, I, L, PVNM( L ), PLST( I,N ), PMAP( L )
               END IF
            END DO LOOP_PM

            NM = NM + I

            STKPM_MAP( N )%LEN = I
            ALLOCATE ( STKPM_MAP( N )%ARRY( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKPM_MAP', PNAME )
            STKPM_MAP( N )%ARRY = 0   ! array

            STKPM_VAR( N )%LEN = I
            ALLOCATE ( STKPM_VAR( N )%ARRY( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKPM_VAR', PNAME )
            STKPM_VAR( N )%ARRY = ' '   ! array

            PM_MAP( N )%LEN = I
            ALLOCATE ( PM_MAP( N )%ARRY( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PM_MAP', PNAME )
            PM_MAP( N )%ARRY = 0   ! array

         END DO   ! NPTGRPS

D        WRITE( LOGDEV,'(5X, A, I3)' ) 'NM =', NM
D         WRITE( LOGDEV,* ) ' '
         L = 0
         DO N = 1, NPTGRPS
D            WRITE( LOGDEV,1021 ) N, STKPM_MAP( N )%LEN
            DO I = 1, STKPM_MAP( N )%LEN
               L = L + 1
               STKPM_MAP( N )%ARRY( I ) = PLST( I,N )
               STKPM_VAR( N )%ARRY( I ) = PVNM( L )
            END DO
D            WRITE( LOGDEV,1025 ) N, STKPM_MAP( N )%ARRY
D            WRITE( LOGDEV,1029 ) N, STKPM_VAR( N )%ARRY

            IF ( N_AE_EMIS .GT. 0 ) THEN   ! assume PM species required
C Check that the PM species names match EMIS_1
               DO I = 1, N_EMIS_PM
                  J = INDEX1( PM_EMIS( I ), STKPM_VAR( N )%LEN, STKPM_VAR( N )%ARRY )
                  IF ( J .LE. 0 ) THEN
                     WRITE( LOGDEV,'(5X, A, I3)' ) TRIM( PM_EMIS( I ) )
     &                                        // ' Not in Stack Emissions file', N
!                    QUIT = .TRUE.
                  END IF
               END DO
            END IF

         END DO

         IF ( QUIT ) THEN
            SUCCESS = .FALSE.; RETURN
         END IF

         ALLOCATE( V( NM ), STAT=IOS )
         ALLOCATE( VC( NM ), STAT=IOS )

         DO I = 1, NM
            V( I )%STR = PVNM( I )
            V( I )%NDX = PMAP( I )
         END DO

C Sort ...
D         WRITE( LOGDEV,'(5X, A)' ) 'Concatenated Species List:'
D         WRITE( LOGDEV,'(5X, I3, 1X, A, I3)' ) 1, V( 1 )
         DO I = 2, NM
D            WRITE( LOGDEV,'(5X, I3, 1X, A, I3)' ) I, V( I )
            TV = V( I )
            DO J = I - 1, 1, -1
               IF ( V( J )%NDX .LE. TV%NDX ) GO TO 201
               V( J+1 ) = V( J )
            END DO
            J = 0
201         CONTINUE
            V( J+1 ) = TV
         END DO

D         WRITE( LOGDEV,'(5X, A)' ) 'Sorted Species List:'
D         DO I = 1, NM
D            WRITE( LOGDEV,'(5X, A, I3)' ) V( I )
D         END DO

C Cull ...
         VC( 1 ) = V( 1 )
         NCULL = 1
         DO I = 2, NM
            IF ( V( I )%STR .NE. VC( NCULL )%STR ) THEN
               NCULL = NCULL + 1
               VC( NCULL ) = V( I )
            END IF
         END DO
D         WRITE( LOGDEV,'(5X, A, I3)' ) 'Unique PM =', NCULL

         N_SPC_PTPM = NCULL

         ALLOCATE ( PTPM_MAP( N_SPC_PTPM ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PTPM_MAP', PNAME )

         ALLOCATE ( STKPM_NAME( N_SPC_PTPM ), STAT=IOS )
         CALL CHECKMEM( IOS, 'STKPM_NAME', PNAME )

         DO I = 1, N_SPC_PTPM
            PTPM_MAP( I ) = VC( I )%NDX
            STKPM_NAME( I ) = VC( I )%STR
D            WRITE( LOGDEV,1043 ) I, STKPM_NAME( I ), PTPM_MAP( I )
         END DO

         DO N = 1, NPTGRPS
            J = 0
            DO I = 1, N_SPC_PTPM
               K = INDEX1( STKPM_NAME( I ), STKSPC( N )%LEN, STKSPC( N )%ARRY )
               IF ( K .NE. 0 ) THEN
                  J = J + 1
                  PM_MAP( N )%ARRY( J ) = I
D                  WRITE( LOGDEV,1033 ) N, J, PM_MAP( N )%ARRY( J )
               END IF
            END DO
D           WRITE( LOGDEV,'(5X, A, I4)' ) 'J =', J
         END DO

         RETURN

1005     FORMAT( / 5X, 'N  J  K  I  L ', 4X, 'GVNM', 15X, 'GLST', 2X, 'GMAP' )
1007     FORMAT( 5X, 'N  J  K  I  L ', 4X, 'PVNM', 15X, 'PLST', 2X, 'PMAP' )
1009     FORMAT( 3X, 5I3, 5X, A16, 2I6 )
1011     FORMAT( 5X, 'STKEM_MAP(', I2, ' )%LEN:', 1X, I3 )
1015     FORMAT( 5X, 'STKEM_MAP(', I2, ' )%ARRY:', 13X, 8I4 )
1019     FORMAT( 5X, 'STKEM_VAR(', I2, ' )%ARRY:', 1X,  8A16 )
1021     FORMAT( 5X, 'STKPM_MAP(', I2, ' )%LEN:', 1X, I3 )
1025     FORMAT( 5X, 'STKPM_MAP(', I2, ' )%ARRY:', 13X, 8(:I4) )
1029     FORMAT( 5X, 'STKPM_VAR(', I2, ' )%ARRY:', 1X, 8(:A16) )
1031     FORMAT( 5X, 'GAS_MAP(', I2, ' )%ARRY(', I3, ' ) =', I3 )
1033     FORMAT( 5X, 'PM_MAP(', I2, ' )%ARRY(', I3, ' ) =', I3 )
1041     FORMAT( 5X, I3, 1X,'STKEM_NAME, PTEM_MAP = ', A16, I3,' Model Species ', A16 )
1043     FORMAT( 5X, I3, 1X,'STKPM_NAME, PTPM_MAP = ', A16, I3 )
1045     FORMAT( I3, 2X, 'PTEM_MAP = ',I3,' STKEM_NAME ',A16,1X,' SPC_PTEM_MAP = ',I3,
     &                  ' into Model Species ', A16 )
1046     FORMAT('WARNING: ',I3,' No Point Emissions of ',A16,' for Species ',A16)
         END FUNCTION PTMAP_INIT

      END MODULE PTMAP

