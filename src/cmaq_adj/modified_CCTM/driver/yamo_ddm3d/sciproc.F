
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_ddm3d/sciproc.F,v 1.1.1.1 2010/06/14 16:03:05 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef sens
      SUBROUTINE SCIPROC ( CGRID, SENGRID, JDATE, JTIME, TSTEP, ASTEP,
     &                     NPMAX, IPT, IPARM, ISPLIT,
     &                     IAMOUNT, ILAYER, STDATE, IDATE,
     &                     ITIME, IREGION, IRXN, SEN_PAR )
#else
      SUBROUTINE SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
#endif

C-----------------------------------------------------------------------
C Function:
C    Controls all of the physical and chemical processes for a grid
C    Operator splitting symmetric around chemistry
 
C Preconditions:
C    Dates and times represented YYYYDDD:HHMMSS.
C    No "skipped" dates and times.  All boundary input variables (layered or
C    non-layered) have the same perimeter structure with a thickness of NTHIK
 
C Subroutines and functions called:
C    All physical and chemical subroutines, 
C    DECOUPLE, COUPLE
 
C Revision History:
C    Oct. 24, 1995 by M. Talat Odman and Clint L. Ingram at NCSC: created
C    13 Dec 97 - Jeff - uncouple diffusion processes
C    27 Jun 98 - Jeff - sync step = chem step
C     7 Jul 01 - Shawn - mv cloud processing before chem
C       Jan 02 - Jeff - dyn alloc; remove PCGRID argument to ping
C    23 Jun 03 J.Young: for layer dependent advection tstep
C    18 Aug 03 J. Pleim - move vdiff before advection
C    29 Aug 03 J.Young: eliminate symmetric processing option and unused
C                       SUBST_GRID_ID string
C    30 May 05 J.Young: mass-conserving advection (yamo)
C     7 Jul 07 J.Young: No cksummer call after HADV - delay possible truncation of
C                       HADV negative conc's until after ZADV
C    25 Jul 08 S.Roselle: removed PING from modeled processes
C    15 Nov 08 S.Napelenok: DDM-3D implementation for CMAQv4.7
C-----------------------------------------------------------------------

!     USE CGRID_SPCS   ! CGRID species number and offsets

      IMPLICIT NONE   

C Include files:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
!     INCLUDE SUBST_AE_SPC      ! aerosol species table
!     INCLUDE SUBST_NR_SPC      ! non-reactive species table
!     INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PACTL_ID    ! PA control parameters

C Arguments:
      REAL, POINTER :: CGRID( :,:,:,: )

      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

!     INTEGER      ASTEP( NLAYS )  ! layer advection time step
      INTEGER      ASTEP( : )      ! layer advection time step

#ifdef sens
      REAL, POINTER           :: SENGRID( :,:,:,:,: )   !
      INTEGER NPMAX                                     ! 
      INTEGER, POINTER        :: IPT( : )               !
      INTEGER, POINTER        :: IPARM( :,:)            !
      INTEGER, POINTER        :: ISPLIT ( :,: )         !
      REAL, POINTER           :: IAMOUNT ( :,:,:)       !
      INTEGER, POINTER        :: ILAYER( :,: )          !
      INTEGER STDATE                                    !  
      INTEGER, POINTER        :: IDATE ( :,: )          !
      INTEGER, POINTER        :: ITIME( :,: )           !
      REAL, POINTER           :: IREGION( :,:,:,: )     !
      INTEGER, POINTER        :: IRXN( :,: )            !
      CHARACTER( 8 ), POINTER :: SEN_PAR( : )           !
#endif

      INTEGER LAY


C Parameters:

C External Functions (not already declared by IODECL3.EXT):

      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC
      LOGICAL, EXTERNAL :: ENVYN

C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'SCIPROC'

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 120 ) :: XMSG = ' '
      CHARACTER( 36 ) :: NMSG = 'After NEXTIME: returned JDATE, JTIME'

      INTEGER, SAVE :: LOGDEV

      INTEGER      SDATE        ! current science process date, coded YYYYDDD
      INTEGER      STIME        ! current science process time, coded HHMMSS
 
      INTEGER      ALLOCSTAT
      INTEGER      STATUS          ! ENV... status
      CHARACTER( 80 ) :: VARDESC   ! environment variable description
      CHARACTER( 80 ) :: MSG = ' '
      CHARACTER( 16 ) :: CTM_CKSUM = 'CTM_CKSUM'     ! env var for cksum on
      LOGICAL, SAVE   :: CKSUM     ! flag for cksum on, default = [T]



      INTERFACE
         SUBROUTINE CKSUMMER ( PNAME, CGRID, JDATE, JTIME )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
         END SUBROUTINE CKSUMMER
         SUBROUTINE PA_UPDATE ( PNAME, CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE PA_UPDATE
#ifdef sens
         SUBROUTINE HADV ( CGRID, SENGRID, JDATE, JTIME, TSTEP, ASTEP,
     &                     NPMAX, IPT, IPARM, STDATE, IDATE,
     &                     ITIME, IREGION, SEN_PAR )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: ASTEP( : )
            INTEGER, INTENT( IN )     :: NPMAX
            INTEGER, POINTER          :: IPT( : )
            INTEGER, POINTER          :: IPARM( :,: )
            INTEGER, INTENT( IN )     :: STDATE
            INTEGER, POINTER          :: IDATE( :,: )
            INTEGER, POINTER          :: ITIME( :,: )
            REAL, POINTER             :: IREGION( :,:,:,: )
            CHARACTER( * ), POINTER   :: SEN_PAR( : )
         END SUBROUTINE HADV
#else
         SUBROUTINE HADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: ASTEP( : )
         END SUBROUTINE HADV
#endif
#ifdef sens
         SUBROUTINE ZADV ( CGRID, SENGRID,JDATE, JTIME, TSTEP, NPMAX )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
         END SUBROUTINE ZADV
#else
         SUBROUTINE ZADV ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE ZADV
#endif
#ifdef sens
         SUBROUTINE COUPLE ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
         END SUBROUTINE COUPLE
         SUBROUTINE DECOUPLE ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
         END SUBROUTINE DECOUPLE
#else
         SUBROUTINE COUPLE ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE COUPLE
         SUBROUTINE DECOUPLE ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE DECOUPLE
#endif

#ifdef sens
         SUBROUTINE HDIFF ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
         END SUBROUTINE HDIFF
#else
         SUBROUTINE HDIFF ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE HDIFF
#endif

#ifdef sens
         SUBROUTINE VDIFF ( CGRID, SENGRID, JDATE, JTIME, TSTEP,
     &                      NPMAX, IPT, IPARM, ISPLIT, IAMOUNT,
     &                      ILAYER, STDATE, IDATE, ITIME, IREGION,
     &                      SEN_PAR )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
            INTEGER, POINTER          :: IPT( : )
            INTEGER, POINTER          :: IPARM( :,: )
            INTEGER, POINTER          :: ISPLIT( :,: )
            REAL, POINTER             :: IAMOUNT( :,:,: )
            INTEGER, POINTER          :: ILAYER( :,: )
            INTEGER, INTENT( IN )     :: STDATE 
            INTEGER, POINTER          :: IDATE( :,: )
            INTEGER, POINTER          :: ITIME( :,: )
            REAL, POINTER             :: IREGION( :,:,:,: )
            CHARACTER( * ), POINTER   :: SEN_PAR( : )
         END SUBROUTINE VDIFF
#else
         SUBROUTINE VDIFF ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE VDIFF
#endif

#ifdef sens
         SUBROUTINE CLDPROC ( CGRID, SENGRID, JDATE, JTIME, 
     &                        TSTEP, NPMAX, SEN_PAR)
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 ) 
            INTEGER, INTENT( IN )     :: NPMAX
            CHARACTER( * ), POINTER   :: SEN_PAR( : ) 
         END SUBROUTINE CLDPROC
#else
         SUBROUTINE CLDPROC ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE CLDPROC
#endif
#ifdef sens
         SUBROUTINE CHEM ( CGRID, SENGRID, JDATE, JTIME, TSTEP,
     &                     NPMAX, IPT, IPARM, STDATE, IDATE, ITIME,
     &                     IREGION, IRXN )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX
            INTEGER, POINTER          :: IPT( : )
            INTEGER, POINTER          :: IPARM( :,: )
            INTEGER, INTENT( IN )     :: STDATE
            INTEGER, POINTER          :: IDATE( :,: )
            INTEGER, POINTER          :: ITIME( :,: )
            REAL, POINTER             :: IREGION( :,:,:,: )
            INTEGER, POINTER          :: IRXN( :,: )
         END SUBROUTINE CHEM
#else
         SUBROUTINE CHEM ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE CHEM
#endif
#ifdef sens
         SUBROUTINE AERO ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            REAL, POINTER             :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
            INTEGER, INTENT( IN )     :: NPMAX

         END SUBROUTINE AERO
#else
         SUBROUTINE AERO ( CGRID, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 2 )
         END SUBROUTINE AERO
#endif

      END INTERFACE
  
C-----------------------------------------------------------------------

C If ISPCA .ne. 0, then air is advected and concs. are adjusted

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

         CKSUM = .TRUE.         ! default
         VARDESC = 'Cksum on flag'
         CKSUM = ENVYN( CTM_CKSUM, VARDESC, CKSUM, STATUS )
            IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
            IF ( STATUS .EQ. 1 ) THEN
              XMSG = 'Environment variable improperly formatted'
              CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
              ELSE IF ( STATUS .EQ. -1 ) THEN
              MSG='Environment variable set, but empty...Using default:'
              WRITE( LOGDEV, '(5X, A, I9)' ) MSG, JTIME
              ELSE IF ( STATUS .EQ. -2 ) THEN
              MSG = 'Environment variable not set ... Using default:'
              WRITE( LOGDEV, '(5X, A, I9)' ) MSG, JTIME
              END IF

         END IF       ! if firstime

C Physical Processes for any Grid

#ifdef sens
      CALL VDIFF ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX, IPT, IPARM,
     &             ISPLIT, IAMOUNT, ILAYER, STDATE, IDATE, ITIME, IREGION,
     &             SEN_PAR )
#else
      CALL VDIFF ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'VDIFF', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'VDIF', CGRID, JDATE, JTIME, TSTEP )

C couple CGRID for advection

#ifdef sens
      CALL COUPLE ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
#else
      CALL COUPLE ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'COUPLE', CGRID, JDATE, JTIME )

#ifdef sens
      CALL HADV ( CGRID, SENGRID, JDATE, JTIME, TSTEP, ASTEP, NPMAX,
     &            IPT, IPARM, STDATE, IDATE, ITIME, IREGION, SEN_PAR )
#else
      CALL HADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
#endif
!     IF ( CKSUM ) CALL CKSUMMER ( 'HADV', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'HADV', CGRID, JDATE, JTIME, TSTEP )

#ifdef sens
      CALL ZADV ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
#else
      CALL ZADV ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'ADV', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'ZADV', CGRID, JDATE, JTIME, TSTEP )

#ifdef sens
      CALL HDIFF ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
#else
      CALL HDIFF ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'HDIFF', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'HDIF', CGRID, JDATE, JTIME, TSTEP )

C decouple CGRID for cloud and chemistry

      SDATE = JDATE
      STIME = JTIME
      CALL NEXTIME ( SDATE, STIME, TSTEP( 2 ) )

#ifdef sens
      CALL DECOUPLE ( CGRID, SENGRID, SDATE, STIME,TSTEP, NPMAX )
#else
      CALL DECOUPLE ( CGRID, SDATE, STIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'DECOUPLE', CGRID, JDATE, JTIME )

#ifdef sens
      CALL CLDPROC ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX, SEN_PAR )
#else
      CALL CLDPROC ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'CLDPROC', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'CLDS', CGRID, JDATE, JTIME, TSTEP )

#ifdef sens
      CALL CHEM (CGRID,SENGRID, JDATE, JTIME, TSTEP, NPMAX, IPT, IPARM,
     &           STDATE, IDATE, ITIME, IREGION, IRXN )
#else
      CALL CHEM ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'CHEM', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'CHEM', CGRID, JDATE, JTIME, TSTEP ) 

#ifdef sens
      CALL AERO ( CGRID, SENGRID, JDATE, JTIME, TSTEP, NPMAX )
#else
      CALL AERO ( CGRID, JDATE, JTIME, TSTEP )
#endif
      IF ( CKSUM ) CALL CKSUMMER ( 'AERO', CGRID, JDATE, JTIME )
      IF ( LIPR ) CALL PA_UPDATE ( 'AERO', CGRID, JDATE, JTIME, TSTEP )

      CALL NEXTIME ( JDATE, JTIME, TSTEP( 2 ) )
      WRITE( LOGDEV,'(/ 5X, A, I8, I7.6)' ) NMSG, JDATE, JTIME

      RETURN
      END
