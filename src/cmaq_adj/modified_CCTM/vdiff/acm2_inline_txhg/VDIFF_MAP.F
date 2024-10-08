
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_txhg/VDIFF_MAP.F,v 1.1.1.1 2010/06/14 16:03:08 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE VDIFF_MAP

C-----------------------------------------------------------------------
C Function: index mapping between CGRID, diffusion, dry dep, and emissions
C           species

C Module depends on prior call to DEPV_INIT to get NDEPV_SPC

C Revision History:
C     10 Apr 2008 J.Young: initial implementation using subroutine vdiff_map.F

C-----------------------------------------------------------------------

      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table

      INCLUDE SUBST_GC_DDEP   ! gas chem dry dep species and map table
      INCLUDE SUBST_AE_DDEP   ! aerosol dry dep species and map table
      INCLUDE SUBST_NR_DDEP   ! non-react dry dep species and map table
      INCLUDE SUBST_TR_DDEP   ! tracer dry dep species and map table

C global diffusion species
      INTEGER, PARAMETER :: N_SPC_DIFF = N_GC_DIFF
     &                                 + N_AE_DIFF
     &                                 + N_NR_DIFF
     &                                 + N_TR_DIFF

C global dry dep species
      INTEGER, PARAMETER :: N_SPC_DDEP = N_GC_DDEP
     &                                 + N_AE_DDEP
     &                                 + N_NR_DDEP
     &                                 + N_TR_DDEP

      INTEGER, SAVE :: DIFF_MAP( N_SPC_DIFF+1 ) ! global diff map to CGRID
      INTEGER, SAVE :: DF2DV   ( N_SPC_DIFF+1 ) ! map from diff spc to depv spc
      INTEGER, SAVE :: DF2EM   ( N_SPC_DIFF+1 ) ! map from diff spc to emis spc
      INTEGER, SAVE :: DD2DV   ( N_SPC_DDEP+1 ) ! map from ddep spc to depv spc
      INTEGER, ALLOCATABLE, SAVE :: DEPV_MAP( : ) ! global depv map to CGRID
      INTEGER, ALLOCATABLE, SAVE :: DV2DF( : )  ! map from depv spc to diff spc

      CHARACTER( 16 ),              SAVE :: DDEP_SPC( N_SPC_DDEP+1 )
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: DEPV_SPC( : )

      CHARACTER( 16 ),              SAVE :: DIFF_SPC( N_SPC_DIFF )
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: DV2DF_SPC( : )

      CONTAINS

C-----------------------------------------------------------------------
         FUNCTION VDIFF_MAP_INIT( NDEPV_SPC ) RESULT ( SUCCESS )

         USE CGRID_SPCS   ! CGRID species number and offsets
         USE AERO_EMIS    ! aerosol emissions defns and fncts

         IMPLICIT NONE

         INCLUDE SUBST_GC_SPC    ! gas species table
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
         INCLUDE SUBST_GC_DEPV   ! gas chem dep vel surrogate names and map table

         INCLUDE SUBST_AE_SPC    ! aerosol species table
         INCLUDE SUBST_AE_DEPV   ! aerosol dep vel surrogate names and map table

         INCLUDE SUBST_NR_SPC    ! nr species table
         INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
         INCLUDE SUBST_NR_DEPV   ! non-react dep vel surrogate names and map table

         INCLUDE SUBST_TR_SPC    ! tracer species table
         INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names and map table
         INCLUDE SUBST_TR_DEPV   ! tracer dep vel surrogate names and map table

         INCLUDE SUBST_IOPARMS
         INCLUDE SUBST_IODECL    ! I/O definitions and declarations ! debug
#ifdef emis_chem
         INCLUDE SUBST_EMPR_CH   ! emissions processing in chem
#else
         INCLUDE SUBST_EMPR_VD   ! emissions processing in vdif
#endif

C Arguments:
         INTEGER, INTENT( IN ) :: NDEPV_SPC ! no. of depv species
         LOGICAL SUCCESS

C global dep vel species
!        INTEGER, PARAMETER :: N_DEPV_SPC = N_GC_DEPV
!    &                                    + N_AE_DEPV
!    &                                    + N_NR_DEPV
!    &                                    + N_TR_DEPV

C global emis species
         INTEGER N_SPC_EMIS
!        INTEGER, PARAMETER :: N_SPC_EMIS = NEMIS
!!   &                                    + N_AE_EMIS
!    &                                    + N_NR_EMIS
!    &                                    + N_TR_EMIS

C External Functions not previously declared in IODECL3.EXT:

         INTEGER, EXTERNAL :: FINDEX

C Local variables:
         
         CHARACTER( 80 )  :: XMSG       ! Message text
         CHARACTER( 16 )  :: PNAME = 'VDIFF_MAP'

         INTEGER NEMIS_AE
         INTEGER STRT, FINI
         INTEGER S, V, N, I, INDX
         INTEGER ALLOCSTAT

         LOGICAL LERROR                    ! Error flag
D        INTEGER LOGDEV   ! debug

C-----------------------------------------------------------------------

D        LOGDEV = INIT3()   ! debug
         SUCCESS = .TRUE.

         IF ( N_AE_SPC .GT. 0 ) THEN
            NEMIS_AE = N_AE_EMIS  ! from AERO_EMIS module
         ELSE
            NEMIS_AE = 0
         END IF

         N_SPC_EMIS = NEMIS
     &              + NEMIS_AE
     &              + N_NR_EMIS
     &              + N_TR_EMIS

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

         ALLOCATE ( DEPV_MAP( NDEPV_SPC+1 ),
     &              DEPV_SPC( NDEPV_SPC ),
!    &              DV2DF( NDEPV_SPC ), STAT = ALLOCSTAT )
     &              DV2DF( NDEPV_SPC ),
     &              dv2df_spc( ndepv_spc ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
!           XMSG = 'Failure allocating DEPV_MAP, DEPV_SPC or DV2DF'
            XMSG = 'Failure allocating DEPV_MAP, DEPV_SPC, DV2DF or dv2df_spc'
            CALL M3WARN( PNAME, 0, 0, XMSG )
            SUCCESS = .FALSE.; RETURN
         END IF

C create global maps to CGRID

         N = 0; S = 0; I = 0
         IF ( EMISVD ) THEN
            DO V = 1, N_GC_DIFF ! assumes emis species a subset of vdif species
               N = N + 1
               INDX = FINDEX( GC_DIFF_MAP( V ), N_GC_EMIS, GC_EMIS_MAP )
               IF ( INDX .GT. 0 ) THEN
                  I = I + 1
                  DF2EM( N ) = S + INDX
               ELSE
                  DF2EM( N ) = N_SPC_EMIS + 1
               END IF
            END DO
         ELSE            ! emissions processing done in chem
            DO V = 1, N_GC_DIFF 
               N = N + 1
               DF2EM( N ) = N_SPC_EMIS + 1
            END DO
         END IF

         S = S + I     ! S should now be = N_GC_EMIS
         I = 0
         DO V = 1, N_AE_DIFF
            N = N + 1
            INDX = FINDEX( AE_DIFF_MAP( V ), N_AE_EMIS, AE_EMIS_MAP )
            IF ( N_AE_SPC .GT. 0 .AND. INDX .GT. 0 ) THEN
               I = I + 1
               DF2EM( N ) = S + INDX
            ELSE
               DF2EM( N ) = N_SPC_EMIS + 1
            END IF
         END DO

         S = S + I     ! S should now be = N_GC_EMIS + N_AE_EMIS
         I = 0
         DO V = 1, N_NR_DIFF
            N = N + 1
            INDX = FINDEX( NR_DIFF_MAP( V ), N_NR_EMIS, NR_EMIS_MAP )
            IF ( INDX .GT. 0 ) THEN
               I = I + 1
               DF2EM( N ) = S + INDX
            ELSE
               DF2EM( N ) = N_SPC_EMIS + 1
            END IF
         END DO

         S = S + I     ! S should now be = N_GC_EMIS + N_AE_EMIS + N_NR_EMIS
         DO V = 1, N_TR_DIFF
            N = N + 1
            INDX = FINDEX( TR_DIFF_MAP( V ), N_TR_EMIS, TR_EMIS_MAP )
            IF ( INDX .GT. 0 ) THEN
               DF2EM( N ) = S + INDX
            ELSE
               DF2EM( N ) = N_SPC_EMIS + 1
            END IF
         END DO

         I = 0; S = 0; N = 0
         DO V = 1, N_GC_DIFF ! assumes dep.vel species subset of vdif species
            N = N + 1
            INDX = FINDEX( GC_DIFF_MAP( V ), N_GC_DEPV, GC_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               I = I + 1
               DF2DV( N ) = S + INDX
            ELSE
               DF2DV( N ) = NDEPV_SPC + 1
            END IF
         END DO
         S = S + I     ! S should now be = N_GC_DEPV

         I = 0
         DO V = 1, N_AE_DIFF
            N = N + 1
            INDX = FINDEX( AE_DIFF_MAP( V ), N_AE_DEPV, AE_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               I = I + 1
               DF2DV( N ) = S + INDX
            ELSE
               DF2DV( N ) = NDEPV_SPC + 1
            END IF
         END DO
         S = S + I     ! S should now be = N_GC_DEPV + N_AE_DEPV

         I = 0
         DO V = 1, N_NR_DIFF
            N = N + 1
            INDX = FINDEX( NR_DIFF_MAP( V ), N_NR_DEPV, NR_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               I = I + 1
               DF2DV( N ) = S + INDX
            ELSE
               DF2DV( N ) = NDEPV_SPC + 1
            END IF
         END DO
         S = S + I     ! S should now be = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV

         DO V = 1, N_TR_DIFF
            N = N + 1
            INDX = FINDEX( TR_DIFF_MAP( V ), N_TR_DEPV, TR_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               DF2DV( N ) = S + INDX
            ELSE
               DF2DV( N ) = NDEPV_SPC + 1
            END IF
         END DO

D        DO N = 1, N_SPC_DIFF   ! debug
D           WRITE( LOGDEV,'(A,3I4)' ) '  N, DF2EM(N), DF2DV(N): ',
D    &                                   N, DF2EM( N ), DF2DV( N )
D        END DO

         LERROR = .FALSE.
         S = 0
         V = 0
         STRT = 1
         FINI = N_GC_DDEP
         DO N = STRT, FINI
            V = V + 1
            INDX = FINDEX ( GC_DDEP_MAP( N ), N_GC_DEPV, GC_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               DD2DV ( V ) = S + INDX
               DDEP_SPC( V ) = GC_DDEP( N )
            ELSE
               XMSG = 'ERROR: No deposition velocity for DDEP species ' //
     &               GC_DDEP( N )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

         S = N_GC_DEPV
         STRT = 1
         FINI = N_AE_DDEP
         DO N = STRT, FINI
            V = V + 1
            INDX = FINDEX ( AE_DDEP_MAP( N ), N_AE_DEPV, AE_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               DD2DV ( V ) = S + INDX
               DDEP_SPC( V ) = AE_DDEP( N )
            ELSE
               XMSG = 'ERROR: No deposition velocity for DDEP species ' //
     &               AE_DDEP( N )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

         S = N_GC_DEPV + N_AE_DEPV
         STRT = 1
         FINI = N_NR_DDEP
         DO N = STRT, FINI
            V = V + 1
            INDX = FINDEX ( NR_DDEP_MAP( N ), N_NR_DEPV, NR_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               DD2DV ( V ) = S + INDX
               DDEP_SPC( V ) = NR_DDEP( N )
            ELSE
               XMSG = 'ERROR: No deposition velocity for DDEP species ' //
     &               NR_DDEP( N )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

         S = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV
         STRT = 1
         FINI = N_TR_DEPV
         DO N = STRT, FINI
            V = V + 1
            INDX = FINDEX ( TR_DDEP_MAP( N ), N_TR_DEPV, TR_DEPV_MAP )
            IF ( INDX .GT. 0 ) THEN
               DD2DV ( V ) = S + INDX
               DDEP_SPC( V ) = TR_DDEP( N )
            ELSE
               XMSG = 'ERROR: No deposition velocity for DDEP species ' //
     &               TR_DDEP( N )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

D        DO N = 1, N_SPC_DDEP   ! debug
D           WRITE( LOGDEV,'(A,2I4,3X,A)' ) ' N, DD2DV(N), DDEP_SPC(N): ',
D    &                                       N, DD2DV( N ), DDEP_SPC( N )
D        END DO

         S = 0
         DO V = 1, N_GC_DEPV
            S = S + 1
            DEPV_MAP( S ) = GC_STRT - 1 + GC_DEPV_MAP( V )
!           DEPV_SPC( S ) = GC_DEPV( V )
            DEPV_SPC( S ) = GC_SPC( GC_DEPV_MAP( V ) )
         END DO
         DO V = 1, N_AE_DEPV
            S = S + 1
            DEPV_MAP( S ) = AE_STRT - 1 + AE_DEPV_MAP( V )
!           DEPV_SPC( S ) = AE_DEPV( V )
            DEPV_SPC( S ) = AE_SPC( AE_DEPV_MAP( V ) )
         END DO
         DO V = 1, N_NR_DEPV
            S = S + 1
            DEPV_MAP( S ) = NR_STRT - 1 + NR_DEPV_MAP( V )
!           DEPV_SPC( S ) = NR_DEPV( V )
            DEPV_SPC( S ) = NR_SPC( NR_DEPV_MAP( V ) )
         END DO
         DO V = 1, N_TR_DEPV
            S = S + 1
            DEPV_MAP( S ) = TR_STRT - 1 + TR_DEPV_MAP( V )
!           DEPV_SPC( S ) = TR_DEPV( V )
            DEPV_SPC( S ) = TR_SPC( TR_DEPV_MAP( V ) )
         END DO

D        DO N = 1, NDEPV_SPC   ! debug
D           WRITE( LOGDEV,'(A,2I4,3X,A)' ) ' N, DEPV_MAP(N), DEPV_SPC(N): ',
D    &                                       N, DEPV_MAP( N ), DEPV_SPC( N )
D        END DO

         S = 0
         DO V = 1, N_GC_DIFF
            S = S + 1
            DIFF_MAP( S ) = GC_STRT - 1 + GC_DIFF_MAP( V )
            DIFF_SPC( S ) = GC_SPC( GC_DIFF_MAP( V ) )
         END DO
         DO V = 1, N_AE_DIFF
            S = S + 1
            DIFF_MAP( S ) = AE_STRT - 1 + AE_DIFF_MAP( V )
            DIFF_SPC( S ) = AE_SPC( AE_DIFF_MAP( V ) )
         END DO
         DO V = 1, N_NR_DIFF
            S = S + 1
            DIFF_MAP( S ) = NR_STRT - 1 + NR_DIFF_MAP( V )
            DIFF_SPC( S ) = NR_SPC( NR_DIFF_MAP( V ) )
         END DO
         DO V = 1, N_TR_DIFF
            S = S + 1
            DIFF_MAP( S ) = TR_STRT - 1 + TR_DIFF_MAP( V )
            DIFF_SPC( S ) = TR_SPC( TR_DIFF_MAP( V ) )
         END DO

D        DO N = 1, N_SPC_DIFF   ! debug
D           WRITE( LOGDEV,'(A,2I4,2X,A)' ) ' N,      DIFF_MAP(N), DIFF_SPC(N): ',
D    &                                  N,      DIFF_MAP( N ), DIFF_SPC( N )
D        END DO


         S = 0
         V = 0
         DO N = 1, N_GC_DEPV
            V = V + 1
            INDX = FINDEX ( GC_DEPV_MAP( N ), N_GC_DIFF, GC_DIFF_MAP )
            IF ( INDX .GT. 0 ) THEN
               DV2DF ( V ) = S + INDX
               DV2DF_SPC( V ) = GC_DIFF( INDX )
            ELSE
               XMSG = 'ERROR: deposition velocity specified for non-diff species '
     &               // GC_SPC( GC_DEPV_MAP( N ) )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

         S = N_GC_DIFF
         DO N = 1, N_AE_DEPV
            V = V + 1
            INDX = FINDEX ( AE_DEPV_MAP( N ), N_AE_DIFF, AE_DIFF_MAP )
            IF ( INDX .GT. 0 ) THEN
               DV2DF ( V ) = S + INDX
               DV2DF_SPC( V ) = AE_DIFF( INDX )
            ELSE
               XMSG = 'ERROR: deposition velocity specified for non-diff species '
     &               // AE_SPC( AE_DEPV_MAP( N ) )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

         S = N_GC_DIFF + N_AE_DIFF
         DO N = 1, N_NR_DEPV
            V = V + 1
            INDX = FINDEX ( NR_DEPV_MAP( N ), N_NR_DIFF, NR_DIFF_MAP )
            IF ( INDX .GT. 0 ) THEN
               DV2DF ( V ) = S + INDX
               DV2DF_SPC( V ) = NR_DIFF( INDX )
            ELSE
               XMSG = 'ERROR: deposition velocity specified for non-diff species '
     &               // NR_SPC( NR_DEPV_MAP( N ) )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

         S = N_GC_DIFF + N_AE_DIFF + N_NR_DIFF
         DO N = 1, N_TR_DEPV
            V = V + 1
            INDX = FINDEX ( TR_DEPV_MAP( N ), N_TR_DIFF, TR_DIFF_MAP )
            IF ( INDX .GT. 0 ) THEN
               DV2DF ( V ) = S + INDX
               DV2DF_SPC( V ) = TR_DIFF( INDX )
            ELSE
               XMSG = 'ERROR: deposition velocity specified for non-diff species '
     &               // TR_SPC( TR_DEPV_MAP( N ) )
               CALL M3MESG( XMSG )
               LERROR = .TRUE.
            END IF
         END DO

D        DO N = 1, NDEPV_SPC   ! debug
D           WRITE( LOGDEV,'(A,2I4,2X,A)' ) ' N, DV2DF(N), DV2DF_SPC(N): ',
D    &                                  N, DV2DF( N ), DV2DF_SPC( N )
D        END DO

         IF ( LERROR ) THEN
            XMSG = 'Stopping because of species mapping error(s) in VDIFF'
            CALL M3WARN( PNAME, 0, 0, XMSG )
            SUCCESS = .FALSE.
         END IF

         RETURN

         END FUNCTION VDIFF_MAP_INIT

      END MODULE VDIFF_MAP
