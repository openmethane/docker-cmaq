
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_txhg/vdiff_map.F,v 1.1.1.1 2010/06/14 16:03:10 sjr Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VDIFF_MAP ( DF2EM, DF2DV, DD2DV, DEPV_MAP, DIFF_MAP, DDEP_SPC,
     &                       DV2DF )

!   Jeff - Dec 00 - move CGRID_MAP into f90 module
!   FSB  - Apr 01 - Turn on AE emissions
!   Jeff - May 01 - Use AERO_EMIS; generalize
!   Jeff - Aug 01 - dyn alloc
!   Jeff - Mar 03 - bug fix for no aero case
!   Jerry - Oct 03 - add dv2df mapping array
!   Jerry - Oct 03 - add error checking for dvep/ddep species

      USE CGRID_SPCS   ! CGRID species number and offsets
      USE AERO_EMIS    ! aerosol emissions defns and fncts

      IMPLICIT NONE

!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC    ! gas species table
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_GC_DEPV   ! gas chem dep vel surrogate names and map table
      INCLUDE SUBST_GC_DDEP   ! gas chem dry dep species and map table
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_SPC    ! aerosol species table
!     INCLUDE SUBST_AE_EMIS   ! aerosol emis surrogate names and map table
      INCLUDE SUBST_AE_DEPV   ! aerosol dep vel surrogate names and map table
      INCLUDE SUBST_AE_DDEP   ! aerosol dry dep species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_SPC    ! nr species table
      INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
      INCLUDE SUBST_NR_DEPV   ! non-react dep vel surrogate names and map table
      INCLUDE SUBST_NR_DDEP   ! non-react dry dep species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table
      INCLUDE SUBST_TR_SPC    ! tracer species table
      INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names and map table
      INCLUDE SUBST_TR_DEPV   ! tracer dep vel surrogate names and map table
      INCLUDE SUBST_TR_DDEP   ! tracer dry dep species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table
!     INCLUDE SUBST_EMLYRS_ID ! emissions layers parameter
      INCLUDE SUBST_IOPARMS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef emis_chem
      INCLUDE SUBST_EMPR_CH   ! emissions processing in chem
#else
      INCLUDE SUBST_EMPR_VD   ! emissions processing in vdif
#endif

! global dry dep species
      INTEGER, PARAMETER :: N_SPC_DDEP = N_GC_DDEP
     &                                 + N_AE_DDEP
     &                                 + N_NR_DDEP
     &                                 + N_TR_DDEP

      INTEGER DD2DV( N_SPC_DDEP+1 )  ! map from ddep spec to depv spec

      CHARACTER( 16 ) :: DDEP_SPC( N_SPC_DDEP+1 )  ! ddep species names

! global dep vel species
      INTEGER, PARAMETER :: N_SPC_DEPV = N_GC_DEPV
     &                                 + N_AE_DEPV
     &                                 + N_NR_DEPV
     &                                 + N_TR_DEPV
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHARACTER( 16 ) :: DEPV_SPC( N_SPC_DEPV+1 )  ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      INTEGER      DEPV_MAP( N_SPC_DEPV+1 )  ! global ddep map to CGRID

! global emis species
      INTEGER N_SPC_EMIS
!     INTEGER, PARAMETER :: N_SPC_EMIS = NEMIS
c    &                                 + N_AE_EMIS
!    &                                 + N_NR_EMIS
!    &                                 + N_TR_EMIS

! global diffusion species
      INTEGER, PARAMETER :: N_SPC_DIFF = N_GC_DIFF
     &                                 + N_AE_DIFF
     &                                 + N_NR_DIFF
     &                                 + N_TR_DIFF

      INTEGER DIFF_MAP( N_SPC_DIFF+1 )  ! global diff map to CGRID
      INTEGER DF2DV   ( N_SPC_DIFF+1 )  ! map from diff spec to depv spec
      INTEGER DF2EM   ( N_SPC_DIFF+1 )  ! map from diff spec to emis spec
      INTEGER DV2DF   ( N_SPC_DEPV )    ! map from depv spec to diff spec

C External Functions not previously declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: FINDEX

C Local variables:
      
      CHARACTER( 80 )  ::   MSG       ! Message text
      CHARACTER( 16 )  :: PNAME = 'VDIFF_MAP'

      INTEGER NEMIS_AE
      INTEGER STRT, FINI
      INTEGER S, V, N, I, INDX

      LOGICAL LERROR                    ! Error flag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INTEGER      LOGDEV   ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOGDEV = INIT3()   ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF ( N_AE_SPC .GT. 0 ) THEN
         NEMIS_AE = N_AE_EMIS  ! from AERO_EMIS module
         ELSE
         NEMIS_AE = 0
         END IF

      N_SPC_EMIS = NEMIS
     &           + NEMIS_AE
     &           + N_NR_EMIS
     &           + N_TR_EMIS

C Get CGRID offsets
 
      CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C create global maps to CGRID

      N = 0
      S = 0
      I = 0
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

      I = 0
      S = 0
      N = 0
      DO V = 1, N_GC_DIFF ! assumes dep.vel species subset of vdif species
         N = N + 1
         INDX = FINDEX( GC_DIFF_MAP( V ), N_GC_DEPV, GC_DEPV_MAP )
         IF ( INDX .GT. 0 ) THEN
            I = I + 1
            DF2DV( N ) = S + INDX
            ELSE
            DF2DV( N ) = N_SPC_DEPV + 1
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
            DF2DV( N ) = N_SPC_DEPV + 1
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
            DF2DV( N ) = N_SPC_DEPV + 1
            END IF
         END DO
      S = S + I     ! S should now be = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV

      DO V = 1, N_TR_DIFF
         N = N + 1
         INDX = FINDEX( TR_DIFF_MAP( V ), N_TR_DEPV, TR_DEPV_MAP )
         IF ( INDX .GT. 0 ) THEN
            DF2DV( N ) = S + INDX
            ELSE
            DF2DV( N ) = N_SPC_DEPV + 1
            END IF
         END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DO N = 1, N_SPC_DIFF   ! debug
c        WRITE( LOGDEV,'(A,3I4)' ) '  N, DF2EM(N), DF2DV(N): ',
c    &                                N, DF2EM( N ), DF2DV( N )
c        END DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
            MSG = 'ERROR: No deposition velocity for DDEP species ' //
     &            GC_DDEP( N )
            CALL M3MESG( MSG )
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
            MSG = 'ERROR: No deposition velocity for DDEP species ' //
     &            AE_DDEP( N )
            CALL M3MESG( MSG )
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
            MSG = 'ERROR: No deposition velocity for DDEP species ' //
     &            NR_DDEP( N )
            CALL M3MESG( MSG )
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
            MSG = 'ERROR: No deposition velocity for DDEP species ' //
     &            TR_DDEP( N )
            CALL M3MESG( MSG )
            LERROR = .TRUE.
            END IF
         END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DO N = 1, N_SPC_DDEP   ! debug
c        WRITE( LOGDEV,'(A,2I4,3X,A)' ) ' N, DD2DV(N), DDEP_SPC(N): ',
c    &                                    N, DD2DV( N ), DDEP_SPC( N )
c        END DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      S = 0
      DO V = 1, N_GC_DEPV
         S = S + 1
         DEPV_MAP( S ) = GC_STRT - 1 + GC_DEPV_MAP( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        DEPV_SPC( S ) = GC_DEPV( V )   ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO
      DO V = 1, N_AE_DEPV
         S = S + 1
         DEPV_MAP( S ) = AE_STRT - 1 + AE_DEPV_MAP( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        DEPV_SPC( S ) = AE_DEPV( V )   ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO
      DO V = 1, N_NR_DEPV
         S = S + 1
         DEPV_MAP( S ) = NR_STRT - 1 + NR_DEPV_MAP( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        DEPV_SPC( S ) = NR_DEPV( V )   ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO
      DO V = 1, N_TR_DEPV
         S = S + 1
         DEPV_MAP( S ) = TR_STRT - 1 + TR_DEPV_MAP( V )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        DEPV_SPC( S ) = TR_DEPV( V )   ! debug
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DO N = 1, N_SPC_DEPV   ! debug
c        WRITE( LOGDEV,'(A,2I4,3X,A)' ) ' N, DEPV_MAP(N), DEPV_SPC(N): ',
c    &                                    N, DEPV_MAP( N ), DEPV_SPC( N )
c        END DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      S = 0
      DO V = 1, N_GC_DIFF
         S = S + 1
         DIFF_MAP( S ) = GC_STRT - 1 + GC_DIFF_MAP( V )
         END DO
      DO V = 1, N_AE_DIFF
         S = S + 1
         DIFF_MAP( S ) = AE_STRT - 1 + AE_DIFF_MAP( V )
         END DO
      DO V = 1, N_NR_DIFF
         S = S + 1
         DIFF_MAP( S ) = NR_STRT - 1 + NR_DIFF_MAP( V )
         END DO
      DO V = 1, N_TR_DIFF
         S = S + 1
         DIFF_MAP( S ) = TR_STRT - 1 + TR_DIFF_MAP( V )
         END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DO N = 1, N_SPC_DIFF   ! debug
c        WRITE( LOGDEV,'(A,2I4)' ) ' N,      DIFF_MAP(N): ',
c    &                               N,      DIFF_MAP( N )
c        END DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      S = 0
      V = 0
      STRT = 1
      FINI = N_GC_DEPV
      DO N = STRT, FINI
         V = V + 1
         INDX = FINDEX ( GC_DEPV_MAP( N ), N_GC_DIFF, GC_DIFF_MAP )
         IF ( INDX .GT. 0 ) THEN
            DV2DF ( V ) = S + INDX
            ELSE
            MSG = 'ERROR: Deposition velocity specified for non-DIFF species '
     &            // GC_SPC( GC_DEPV_MAP( N ) )
            CALL M3MESG( MSG )
            LERROR = .TRUE.
            END IF
         END DO

      S = N_GC_DIFF
      STRT = 1
      FINI = N_AE_DEPV
      DO N = STRT, FINI
         V = V + 1
         INDX = FINDEX ( AE_DEPV_MAP( N ), N_AE_DIFF, AE_DIFF_MAP )
         IF ( INDX .GT. 0 ) THEN
            DV2DF ( V ) = S + INDX
            ELSE
            MSG = 'ERROR: Deposition velocity specified for non-DIFF species '
     &            // AE_SPC( AE_DEPV_MAP( N ) )
            CALL M3MESG( MSG )
            LERROR = .TRUE.
            END IF
         END DO

      S = N_GC_DIFF + N_AE_DIFF
      STRT = 1
      FINI = N_NR_DEPV
      DO N = STRT, FINI
         V = V + 1
         INDX = FINDEX ( NR_DEPV_MAP( N ), N_NR_DIFF, NR_DIFF_MAP )
         IF ( INDX .GT. 0 ) THEN
            DV2DF ( V ) = S + INDX
            ELSE
            MSG = 'ERROR: Deposition velocity specified for non-DIFF species '
     &            // NR_SPC( NR_DEPV_MAP( N ) )
            CALL M3MESG( MSG )
            LERROR = .TRUE.
            END IF
         END DO

      S = N_GC_DIFF + N_AE_DIFF + N_NR_DIFF
      STRT = 1
      FINI = N_TR_DEPV
      DO N = STRT, FINI
         V = V + 1
         INDX = FINDEX ( TR_DEPV_MAP( N ), N_TR_DIFF, TR_DIFF_MAP )
         IF ( INDX .GT. 0 ) THEN
            DV2DF ( V ) = S + INDX
            ELSE
            MSG = 'ERROR: Deposition velocity specified for non-DIFF species '
     &            // TR_SPC( TR_DEPV_MAP( N ) )
            CALL M3MESG( MSG )
            LERROR = .TRUE.
            END IF
         END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DO N = 1, N_SPC_DEPV   ! debug
c        WRITE( LOGDEV,'(A,2I4)' ) ' N, DV2DF(N): ',
c    &                               N, DV2DF( N )
c        END DO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      IF( LERROR ) THEN
        MSG = 'Stopping because of species mapping error(s) in VDIFF'
        CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT3 )
      ENDIF

      RETURN
      END
