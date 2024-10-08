
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

C RCS file, release, date & time of last delta, author, state, [locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_txhg/map_ae_trace.F,v 1.1.1.1 2010/06/14 16:03:00 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


         SUBROUTINE MAP_AE_TRACE( N_AERO_MAP, AERO_MAP, JDATE, JTIME, LOGDEV )
C**********************************************************************
C
C  Function:  Find and map trace aerosol species 
C  CALLED BY: AERO_DRIVER
C
C  REVISION HISTORY: 01/01/06 : B.Hutzell - Initial version
C
C**********************************************************************

        USE AE_TRACE

         IMPLICIT NONE
C.....INCLUDES:

         INCLUDE SUBST_CONST         ! constants
         INCLUDE SUBST_AE_SPC        ! aerosol species table
         INCLUDE SUBST_GC_SPC        ! gas chemistry species table
         INCLUDE SUBST_NR_SPC        ! non-reactive species table

         INCLUDE SUBST_IOPARMS       ! I/O parameters definitions
         INCLUDE SUBST_IOFDESC       ! file header data structure
         INCLUDE SUBST_IODECL        ! I/O definitions and declarations


C.....ARGUMENTS:

         INTEGER :: N_AERO_MAP, AERO_MAP( N_AE_SPCD )
         
         INTEGER :: JDATE        ! current model date , coded YYYYDDD
         INTEGER :: JTIME        ! current model time , coded HHMMSS
         INTEGER :: LOGDEV       ! unit number for the log file


C.....LOCAL VARIABLES:

         CHARACTER(16)  :: PNAME =  'MAP_AE_TRACE  '     ! name of routine
         CHARACTER(16)  :: EMPTY
         CHARACTER(16)  :: XNAME                         ! SCRATCH variable
         CHARACTER(16)  :: WNAME                         ! SCRATCH variable
         CHARACTER(16)  :: VNAME                         ! SCRATCH variable
         CHARACTER(128) :: XMSG = 'FATAL ERROR in MAP_AE_TRACE'

         INTEGER :: TEMP_MAP( 3, N_AE_SPCD)  ! stores locations in GRID
         INTEGER :: LOCA_MAP( 3, N_AE_SPCD)  ! stores locations in AE_TRACE_CONC

         INTEGER :: MARKER, N, M       ! indexes
         INTEGER :: I, J, K, L         ! loop counters
         INTEGER :: LEN_NAME           ! number of nonblank characters in species name


C.....EXTERNAL FUNCTIONS:

         INTEGER, EXTERNAL :: INDEX1   ! IOAPI function that return an array index

C**********************************************************************

        

C..Get CGRID offsets and allocate array for saved and change concentrations

         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C..Initialize variables
 

         ALLOCATE(AE_TRACE_CONC( N_AE_SPCD ))
         ALLOCATE(  AE_TRACE_MAP(N_AE_SPCD )) 
         ALLOCATE(AE_TRACE_NAME( N_AE_SPCD ))

         N_TRACE = 0
         N_TRACE_I = 0
         N_TRACE_J = 0
         N_TRACE_K = 0

         TEMP_MAP      = -1
         LOCA_MAP      = -1
         AE_TRACE_MAP  = -1
         AE_TRACE_NAME = " "
         
         N_TRIO = 0
         N_PAIR_IJ = 0
         N_PAIR_IK = 0
         N_PAIR_JK = 0


         AERO_LOOP: DO I = 1, N_AE_SPCD

            VNAME = AE_SPC( I )

            IF(VNAME(1:1) .NE. 'A')CYCLE AERO_LOOP    ! skip non mass conc. variable

            MAP_LOOP: DO J = 1, N_AERO_MAP

               MARKER = I - 1 + AE_STRT

               IF(MARKER .EQ. AERO_MAP( J ))CYCLE AERO_LOOP

            ENDDO MAP_LOOP

            LEN_NAME = LEN_TRIM( VNAME )

            N_TRACE = N_TRACE + 1
            AE_TRACE_MAP(  N_TRACE ) = MARKER
            AE_TRACE_NAME( N_TRACE ) = VNAME

C        Identify the mode of the mass concentration variable

            IF(VNAME(LEN_NAME:LEN_NAME) .EQ. 'I')THEN

               N_TRACE_I = N_TRACE_I + 1
               TEMP_MAP( 1, N_TRACE_I )  = MARKER
               LOCA_MAP( 1, N_TRACE_I )  = N_TRACE

            ENDIF

            IF(VNAME(LEN_NAME:LEN_NAME) .EQ. 'J')THEN

               N_TRACE_J = N_TRACE_J + 1 
               TEMP_MAP( 2, N_TRACE_J )  = MARKER
               LOCA_MAP( 2, N_TRACE_J )  = N_TRACE

            ENDIF

            IF(VNAME(LEN_NAME:LEN_NAME) .NE. 'I' .AND.
     &                VNAME(LEN_NAME:LEN_NAME) .NE. 'J')THEN

               N_TRACE_K = N_TRACE_K + 1 
               TEMP_MAP( 3, N_TRACE_K )  = MARKER
               LOCA_MAP( 3, N_TRACE_K )  = N_TRACE

            ENDIF
            

          ENDDO AERO_LOOP

          IF(N_TRACE .LT. 1)THEN
            WRITE(LOGDEV,*)'NO TRACE SPECIES FOUND IN AEROSOL'
            RETURN
          ENDIF

          ALLOCATE(TRACE_TRIO( N_TRACE ))
          ALLOCATE(AE_TRACE_TRIO(  N_TRACE , 3))
          ALLOCATE(TRACE_PAIR_IJ(  N_TRACE , 2))
          ALLOCATE(TRACE_PAIR_IK(  N_TRACE , 2))
          ALLOCATE(TRACE_PAIR_JK(  N_TRACE , 2))


          TRACE_TRIO    = " "
          AE_TRACE_TRIO = -1
          TRACE_PAIR_IJ = -1
          TRACE_PAIR_IK = -1
          TRACE_PAIR_JK = -1


          WRITE(LOGDEV, *)' REPORT ON AE_TRACE MAPS '
          WRITE(LOGDEV,1000)'     SPECIES    ','  MODE   ','CGRID INDEX '

          IF(N_TRACE_I .GT. 0)THEN
             ALLOCATE(AE_TRACE_I(N_TRACE_I))
C
C   Search for an individual substance across modes. The search assumes that the names
C   are found between the first and last characters of the aerosol species names.
C
             LOOP_TRACE_I: DO I = 1, N_TRACE_I      ! start with a substance in Aitken mode

                AE_TRACE_I( I )   = TEMP_MAP( 1, I)
                J = AE_TRACE_I( I ) + 1 - AE_STRT 
                LEN_NAME = LEN_TRIM( AE_SPC( J ) )  - 1
                VNAME = AE_SPC( J )(2:LEN_NAME) ! assuming sub

                LOOP_IJ: DO K = 1, N_TRACE_J         !  search across accumulation mode
                   J = TEMP_MAP( 2, K) + 1 - AE_STRT 
                   LEN_NAME = LEN_TRIM( AE_SPC( J ) )  - 1
                   WNAME    = AE_SPC( J )(2:LEN_NAME)

                   IF(VNAME .EQ. WNAME)THEN
                      N_PAIR_IJ = N_PAIR_IJ + 1
                      TRACE_PAIR_IJ(N_PAIR_IJ , 1) = LOCA_MAP( 1, I)
                      TRACE_PAIR_IJ(N_PAIR_IJ , 2) = LOCA_MAP( 2, K)
C
C                    Determine if substance also exists in coarse mode.
C
                      LOOP_TRIO: DO L = 1, N_TRACE_K  
                         J = TEMP_MAP( 3, L) + 1 - AE_STRT 
                         LEN_NAME = LEN_TRIM( AE_SPC( J ) )  - 1
                         XNAME    = AE_SPC( J )(2:LEN_NAME)                         

                         IF(VNAME .EQ. XNAME)THEN
                            N_TRIO = N_TRIO + 1
                            TRACE_TRIO(   N_TRIO  )  = VNAME
                            AE_TRACE_TRIO(N_TRIO, 1) = LOCA_MAP( 1, I)
                            AE_TRACE_TRIO(N_TRIO, 2) = LOCA_MAP( 2, K)
                            AE_TRACE_TRIO(N_TRIO, 3) = LOCA_MAP( 3, L)
                         ENDIF
                      ENDDO LOOP_TRIO
                    ENDIF
                ENDDO LOOP_IJ         
C
C    Use the above method and identify substances that exist in only in Aitken
C    and coarse modes.
C
                LOOP_IK: DO K = 1, N_TRACE_K
                   J = TEMP_MAP( 3, K) + 1 - AE_STRT 
                   LEN_NAME = LEN_TRIM( AE_SPC( J ) )  - 1
                   WNAME    = AE_SPC( J )(2:LEN_NAME)

                   IF(VNAME .EQ. WNAME)THEN
                      N_PAIR_IK = N_PAIR_IK + 1
                      TRACE_PAIR_IK(N_PAIR_IK , 1) = LOCA_MAP( 1, I)
                      TRACE_PAIR_IK(N_PAIR_IK , 2) = LOCA_MAP( 3, K)
                   ENDIF
                ENDDO LOOP_IK

                WNAME = ' '
                IF(TRACE_PAIR_IJ(I,2) .GT. 0)THEN
                   WNAME = AE_TRACE_NAME(TRACE_PAIR_IJ(I,2))
                ENDIF

                XNAME = ' '
                IF(TRACE_PAIR_IK(I,2) .GT. 0)THEN
                   XNAME = AE_TRACE_NAME(TRACE_PAIR_IK(I,2))
                ENDIF

                WRITE(LOGDEV, 1100)VNAME,' AITKEN ',AE_TRACE_I( I ),
     &                             WNAME, XNAME
             ENDDO LOOP_TRACE_I
          ENDIF
C
C    Use the above method and identify substances that exist in only in the 
C    accumulation and coarse modes.
C
          IF(N_TRACE_J .GT. 0)THEN

             ALLOCATE(AE_TRACE_J(N_TRACE_J))
             LOOP_TRACE_J: DO I = 1, N_TRACE_J

                AE_TRACE_J( I ) = TEMP_MAP( 2, I)
                J = AE_TRACE_J( I ) + 1 - AE_STRT 
                LEN_NAME = LEN_TRIM( AE_SPC( J ) )  - 1
                VNAME = AE_SPC( J )(2:LEN_NAME)

                LOOP_JK: DO K = 1, N_TRACE_K

                   J = TEMP_MAP( 3, K) + 1 - AE_STRT 
                   LEN_NAME = LEN_TRIM( AE_SPC( J ) )  - 1
                   WNAME    = AE_SPC( J )(2:LEN_NAME)
                   IF(VNAME .EQ. WNAME)THEN
                      N_PAIR_JK = N_PAIR_JK + 1
                      TRACE_PAIR_JK(N_PAIR_JK , 1) = LOCA_MAP( 2, I)
                      TRACE_PAIR_JK(N_PAIR_JK , 2) = LOCA_MAP( 3, K)
                   ENDIF

               ENDDO LOOP_JK

                WNAME = ' '
                IF(TRACE_PAIR_JK(N_PAIR_JK,2) .GT. 0)THEN
                   WNAME = AE_TRACE_NAME(TRACE_PAIR_JK(I,2))
                ENDIF
                J = AE_TRACE_J( I ) + 1 - AE_STRT
                WRITE(LOGDEV, 1100)AE_SPC( J ),' ACCUM  ',AE_TRACE_J( I ),
     &                             WNAME

             ENDDO LOOP_TRACE_J
          ENDIF

          IF(N_TRACE_K .GT. 0)THEN

             ALLOCATE(AE_TRACE_C(N_TRACE_K))
             LOOP_TRACE_K: DO I = 1, N_TRACE_K

                AE_TRACE_C( I ) = TEMP_MAP( 3, I)
                J = AE_TRACE_C( I ) + 1 - AE_STRT 
                WRITE(LOGDEV, 1100)AE_SPC( J ),' COARSE ',AE_TRACE_C( I )

             ENDDO LOOP_TRACE_K

          ENDIF

          WRITE(LOGDEV,*)' TRIO FOUND OVER MODES '

          DO I = 1, N_TRIO

            VNAME = AE_TRACE_NAME(AE_TRACE_TRIO(I, 1))
            WNAME = AE_TRACE_NAME(AE_TRACE_TRIO(I, 2))
            XNAME = AE_TRACE_NAME(AE_TRACE_TRIO(I, 3))
            WRITE(LOGDEV,1200)TRACE_TRIO(I),VNAME,WNAME,XNAME

          ENDDO

1000     FORMAT(A16,2X,A9,4X,A12,2X,A34)
1100     FORMAT(A16,2X,A9,12X,I4,2(2X,A16))
1200     FORMAT(4(A16,2X))

         RETURN
         END SUBROUTINE MAP_AE_TRACE

