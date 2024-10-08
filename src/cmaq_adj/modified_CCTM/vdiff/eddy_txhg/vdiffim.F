
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/eddy_txhg/vdiffim.F,v 1.1.1.1 2010/06/14 16:03:11 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)vdiffim.F 1.8 /project/mod3/CMAQ/src/vdiff/eddy/SCCS/s.vdiffim.F 25 Jul 1997 12:57:45

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VDIFF ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   calculates and writes dry deposition.
C   calculates vertical diffusion
C   controlled by flag THETA, using Crank-Nicolson difference scheme
C   THETA : Crank-Nicolson index  [ 1, fully implicit | 0, fully explicit ]

C The associated tri-diagonal system is stored in 3 arrays
C   DI: diagonal
C   LI: sub-diagonal
C   UI: super-diagonal
C   BI: right hand side function
C   XI: return solution from tridiagonal solver

C     [ DI(1) UI(1) 0     0     0 ...       0     ]
C     [ LI(2) DI(2) UI(2) 0     0 ...       .     ]
C     [ 0     LI(3) DI(3) UI(3) 0 ...       .     ]
C     [ .           .     .     .           .     ] XI(i) = BI(i)
C     [ .                 .     .     .     0     ]
C     [ .                       .     .     .     ]
C     [ 0                             LI(n) DI(n) ]

C   where n = NLAYS

C Preconditions:
C   Z, EDDYV calculated by preprocessor.

C Subroutines and Functions Called:
C   INIT3, SEC2TIME, TIME2SEC, WRITE3, NEXTIME,
C   M3EXIT, EDYINTB, TRIDIAG, PA_UPDATE_EMIS, PA_UPDATE_DDEP

C Revision History:
C   Written by   :    Hwang & Kiran & Carlie
C   Adapted 7/1/93 by CJC to use INTERP3(); interpolates only to the
C   middle of the model time step, for compatibility with RADM LCM
C   2/17/94 by Dongming Hwang at NCSC: 
C          1. Add dry deposition species mapping arrays DD2C, DF2C, DF2DD, DF2EM
C          2. Add emission process.
C   8/16/94 by Dongming Hwang at NCSC: CM template
C   12/08/95 Jeff Young for upgraded Models-3 I/O API, 24/04/96 other mods
C   24 Aug 96 Jeff - emis in vdif/chem option, split out EMLYRS include
C   12 Dec 96 Jerry - added proc. analysis capability
C   17 Apr 97 Jeff - gas chem, aersols, non-reactive, and tracer species
C    3 May 97 Jeff - aerosol emissions never done here
C   19 Jul 97 Jeff - improve N_SPC_* for dimensioning, loops
C   21 Jul 97 Jeff - split out into subroutines
C    9 Sep 97 Jeff - generalized coord.
C    9 Oct 97 Jeff - produce combined factor for aerosol species (in DD_CONV)
C   13 Dec 97 Jeff - now expects CGRID in uncoupled mixing ratios (not nec. unif
C   30 Mar 98 Jeff - now expects emis in mole/s, not g/s (VDEMIS)
C    4 Jun 98 Jeff - write aerosol number deposition as #/Ha
C    2 October, 1998 by Al Bourgeois at LM: parallel implementation
C   12 Jan, 1999 by David Wong at LM: removed NCOLS and NROWS in the argument
C                                     list
C   14 Apr 00 Jerry - now computes ddep for all depv species rather than
C                     ddep species only ( for PA )
C   22 Nov 00 J.Young: GLOBAL_RSUM -> Dave Wong's f90 stenex GLOBAL_SUM
C   05 Apr 01 F. Binkowski: Initial modification for Aerosol Emissions.
C   26 Jun 01 J.Young: aerosol emissions - VDEMIS passed to pa_update with
C                      all emissions rates converted to ppmV/sec
C                      For parallel applicataion - NOTE: VDEMIS must have
C                      the layer dimension the same as the full domain because
C                      the arguments to SE_DATA_COPY require it.
C   30 Aug 01 J.Young: dyn alloc - Use HGRD_DEFN
C   03 Sep 01 David Wong: new pario: allow only PE 0 to call OPDDEP
C   25 Mar 03 J.Young: fix bug for no aero; improve data locality for speedup
C   20 Aug 03 J.Young: fix aerosol treatment - ppm units, not mass density
C   23 Oct 03 J.Young: change depv from flux to m/s, fix integration units;
C                      update conv_grid for MDATE/MTIME, not JDATE/JTIME
C   31 Oct 03 J.Gipson: fix ddep_pa species array dimension ( for PA ),
C                          must be n_spc_depv, not n_spc_ddep
C   31 Oct 03 J.Gipson: use concs at begin & end of time step for ddep calc
C   14 Nov 03 J.Young: fix total to vertical Jacobian product for eddyv
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module (GRID_CONF)
C    7 Jun 05 P.Bhave: added call to OPSSEMIS if MECHNAME='AE4';
C                      added TSTEP to RDEMIS_AE call vector
C-----------------------------------------------------------------------

      USE AERO_EMIS           ! inherits GRID_CONF

      USE SUBST_MODULES              ! stenex
!     USE SUBST_GLOBAL_SUM_MODULE    ! stenex

      IMPLICIT NONE

!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_RXCMMN    ! model mechanism name

      INCLUDE SUBST_GC_SPC    ! gas chemistry species table
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_GC_DEPV   ! gas chem dep vel surrogate names and map table
      INCLUDE SUBST_GC_DDEP   ! gas chem dry dep species and map table
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table

      INCLUDE SUBST_AE_SPC    ! aerosol species table
!     INCLUDE SUBST_AE_EMIS   ! aerosol emis surrogate names and map table
      INCLUDE SUBST_AE_DEPV   ! aerosol dep vel surrogate names and map table
      INCLUDE SUBST_AE_DDEP   ! aerosol dry dep species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table

      INCLUDE SUBST_NR_SPC    ! non-reactive species table
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
#ifdef emis_chem
      INCLUDE SUBST_EMPR_CH   ! emissions processing in chem
#else
      INCLUDE SUBST_EMPR_VD   ! emissions processing in vdif 
#endif
      INCLUDE SUBST_PACTL_ID  ! PA control parameters
      INCLUDE SUBST_CONST     ! constants
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include      SUBST_IODECL    # I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID  ! coordinate and domain definitions (req IOPARMS)

      CHARACTER( 120 ) :: XMSG = ' '

C Arguments:

!     REAL         CGRID( NCOLS,NROWS,NLAYS,* )  !  concentrations
!     REAL      :: CGRID( :,:,:,: )              !  concentrations
      REAL, POINTER :: CGRID( :,:,:,: )              !  concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C Parameters:

C explicit, THETA = 0, implicit, THETA = 1
      REAL, PARAMETER :: THETA = 0.5,  ! Semi-implicit (Crank-Nicolson)
     &                   THBAR = 1.0 - THETA
      REAL THRAT  ! THBAR/THETA

      INTEGER, PARAMETER :: N_SPC_DDEP = N_GC_DDEP
     &                                 + N_AE_DDEP
     &                                 + N_NR_DDEP
     &                                 + N_TR_DDEP

C global dep vel species
      INTEGER, PARAMETER :: N_SPC_DEPV = N_GC_DEPV
     &                                 + N_AE_DEPV
     &                                 + N_NR_DEPV
     &                                 + N_TR_DEPV

C global diffusion species
      INTEGER, PARAMETER :: N_SPC_DIFF = N_GC_DIFF
     &                                 + N_AE_DIFF
     &                                 + N_NR_DIFF
     &                                 + N_TR_DIFF

C global emissions species
      INTEGER, SAVE      :: N_SPC_EMIS
!     INTEGER, PARAMETER :: N_SPC_EMIS = NEMIS
!!   &                                 + N_AE_EMIS
!    &                                 + N_NR_EMIS
!    &                                 + N_TR_EMIS

C number of CGRID species

      INTEGER, SAVE :: N_SPC_CGRID

C number of species on the PM emissions input file. Set in OPEMIS
C the value changes with the type of emissions file.
      INTEGER, SAVE :: NAESPCEMIS

      REAL, PARAMETER :: M2PHA = 1.0E+04       ! 1 hectare = 1.0e4 m**2
      REAL, PARAMETER :: CMLMR = 1.0E+06       ! ppmV/Molar Mixing Ratio
      REAL, PARAMETER :: CNVTD = M2PHA / CMLMR / MWAIR ! combined ddep
                                                       ! conversion factor
      REAL, PARAMETER :: GPKG = 1.0E+03        ! g/Kg
      REAL, PARAMETER :: MGPG = 1.0E+06        ! micro-g/g

C External Functions not previously declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: SECSDIFF, SEC2TIME, TIME2SEC
      LOGICAL, EXTERNAL :: ENVYN

C File variables:

      REAL         RDEPVHT( NCOLS,NROWS )        ! air dens / dep vel height
      REAL         RJACM  ( NCOLS,NROWS,NLAYS )  ! reciprocal mid-layer Jacobian
      REAL         RVJACMF( NCOLS,NROWS,NLAYS )  ! 1/ mid-full layer vert Jac
      REAL         RRHOJ  ( NCOLS,NROWS,NLAYS )  ! reciprocal density X Jacobian
      REAL         DENS1  ( NCOLS,NROWS )        ! layer 1 air density
!     REAL         DEPV   ( NCOLS,NROWS,N_SPC_DEPV+1 ) ! deposition velocities
!     REAL         MDEPV  ( NCOLS,NROWS,N_SPC_DEPV+1 ) ! deposition velocities
                                                       ! X air density for all
                                                       ! but aerosol species
      REAL         DEPV  ( N_SPC_DEPV+1,NCOLS,NROWS ) ! dep vel X dens/msfx2

C Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'VDIFFIM'
      CHARACTER( 16 ), SAVE :: DDEP_SPC( N_SPC_DDEP + 1 )
      CHARACTER( 16 ), SAVE :: CTM_SSEMDIAG = 'CTM_SSEMDIAG'  ! env var for SSEMDIAG file
      CHARACTER( 80 ) :: VARDESC                ! environment variable description

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: WSTEP  = 0               ! local write counter
      INTEGER  STATUS                           ! ENV... status

      REAL          DX1, DX2                    ! CX x1- and x2-cell widths
      REAL          DX3F ( NLAYS )              ! layer thickness at face
!     REAL, SAVE :: RDX3F( NLAYS )              ! reciprocal layer thickness
      REAL, ALLOCATABLE, SAVE :: RDX3F( : )     ! reciprocal layer thickness
      REAL          X3M  ( NLAYS )              ! middle layer heigth
      REAL          DX3M ( NLAYS )              ! layer thickness at middle
!     REAL, SAVE :: RDX3M( NLAYS )              ! reciprocal layer thickness
      REAL, ALLOCATABLE, SAVE :: RDX3M( : )     ! reciprocal layer thickness
      REAL          CONVPA           ! conversion factor to pressure in Pascals
      REAL          CONVEM           ! conversion for emissions rates to Kg/s
      REAL, SAVE :: CNVTE                       ! combined conversion factor
      REAL       :: CNVTR                       ! combined conversion factor

      REAL, ALLOCATABLE, SAVE :: CNGRD( :,:,:,: )  ! cgrid replacement

      REAL, ALLOCATABLE, SAVE :: DDEP( :,:,: )  ! ddep accumulator
      REAL, ALLOCATABLE, SAVE :: DDEP_PA( :,:,: )! ddep for process analysis
      REAL          WRDD( NCOLS,NROWS )         ! ddep write buffer

      INTEGER, SAVE :: DEPV_MAP( N_SPC_DEPV+1 ) ! global depv map to CGRID
      INTEGER, SAVE :: DIFF_MAP( N_SPC_DIFF+1 ) ! global diff map to CGRID
      INTEGER, SAVE :: DF2DV   ( N_SPC_DIFF+1 ) ! map from diff spc to depv spc
      INTEGER, SAVE :: DF2EM   ( N_SPC_DIFF+1 ) ! map from diff spc to emis spc
      INTEGER, SAVE :: DD2DV   ( N_SPC_DDEP+1 ) ! map from ddep spc to depv spc
      INTEGER, SAVE :: DV2DF   ( N_SPC_DEPV )   ! map from depv spc to diff spc

      INTEGER, SAVE :: ELAYS                    ! no. of emis integration layers
                                                ! ELAYS must be .LT. NLAYS
      INTEGER          EMISLYRS                 ! no. of file emissions layers

!     REAL         VDEMIS ( NCOLS,NROWS,NLAYS,N_SPC_EMIS+1 )    ! file input
      REAL, ALLOCATABLE, SAVE :: VDEMIS( :,:,:,: ) ! total emissions array
      REAL, ALLOCATABLE, SAVE :: EMIS_PA( :,:,:,: ) ! emis for process analysis
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE( :,:,:,: ) ! aerosol emissions
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR( :,:,:,: ) ! nonreactive gas emis
      REAL, ALLOCATABLE, SAVE :: VDEMIS_TR( :,:,:,: ) ! tracer emissions

      LOGICAL, SAVE :: EM_TRAC = .FALSE.        ! do tracer emissions?
      LOGICAL, SAVE :: SSEMDIAG                 ! flag for creating SSEMIS output file
      INTEGER, SAVE :: NEMIS_AE                 ! no. of aero emis species

      REAL, SAVE :: DD_CONV( N_SPC_DEPV+1 )     ! ddep spc conversion factors
      
      REAL         DD_FAC( N_SPC_DEPV)          ! combined subexpression
      REAL           DDBF( N_SPC_DEPV)          ! secondary DDEP
      REAL           CONC( N_SPC_DIFF,NLAYS )   ! secondary CGRID expression
      REAL            LI ( NLAYS )              ! see under Function, above
      REAL            DI ( N_SPC_DIFF,NLAYS )   ! see under Function, above
      REAL            BI ( N_SPC_DIFF,NLAYS )   ! see under Function, above
      REAL            UI ( NLAYS )              ! see under Function, above
      REAL            XI ( N_SPC_DIFF,NLAYS )   ! see under Function, above
      REAL            PK ( NLAYS )              ! UI subfactor
      REAL            QK ( NLAYS )              ! LI subfactor
      REAL           EMIS( N_SPC_DIFF,NLAYS )   ! emissions subexpression
      REAL           DVD1( N_SPC_DIFF )         ! dep vel subexpression
      REAL           DVB1( N_SPC_DIFF )         ! dep vel subexpression
      REAL         DVFACB( NCOLS,NROWS )        ! dep vel factor in BI (r.h.s.)
      REAL         DVFACD( NCOLS,NROWS )        ! dep vel factor in DI (diag.)
      REAL         EDDYV ( NCOLS,NROWS,NLAYS )  ! from EDYINTB
      REAL         SEDDY ( NLAYS,NCOLS,NROWS )  ! flipped EDDYV
      INTEGER      NSTEPS( NCOLS,NROWS )        ! diffusion time steps
      REAL         DT    ( NCOLS,NROWS )        ! eddy diff. delta T
      REAL         DELT                         ! DT
      REAL         DTDENS1                      ! DT * layer 1 air density
      REAL         DTSEC                        ! model time step in seconds

      INTEGER, SAVE :: LOGDEV
 
      INTEGER      ALLOCSTAT
      INTEGER      C, R, L, S, V, N             ! loop induction variables
      INTEGER      STRT, FINI                   ! loop induction variables
      INTEGER      MDATE, MTIME, MSTEP          ! internal simulation date&time
                                                ! array in vert. mixing
      LOGICAL, SAVE :: EDDY_STATS = .FALSE.
      REAL DT_AVG                               ! avg eddy delta T 
      REAL NSTP_AVG                             ! avg no. of integration steps

      CHARACTER( 16 ) :: VNAME

      INTERFACE
         SUBROUTINE RDMET( MDATE, MTIME, RDEPVHT, RJACM, RVJACMF, RRHOJ, DENS1 )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME
            REAL, INTENT( OUT )         :: RDEPVHT( :,: )
            REAL, INTENT( OUT )         :: RJACM  ( :,:,: )
            REAL, INTENT( OUT )         :: RVJACMF( :,:,: )
            REAL, INTENT( OUT )         :: RRHOJ  ( :,:,: )
            REAL, INTENT( OUT )         :: DENS1  ( :,: )
         END SUBROUTINE RDMET
         SUBROUTINE RDDEPV ( MDATE, MTIME, MSTEP, CGRID, DEPV )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, MSTEP
            REAL, POINTER               :: CGRID( :,:,:,: )
            REAL, INTENT( OUT )         :: DEPV( :,:,: )
         END SUBROUTINE RDDEPV
         SUBROUTINE RDEMIS_GC ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, VDEMIS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, EMISLYRS, NSPC_EMIS
            REAL, INTENT( OUT )         :: VDEMIS( :,:,:,: )
         END SUBROUTINE RDEMIS_GC
         SUBROUTINE RDEMIS_NR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, VDEMIS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, EMISLYRS, NSPC_EMIS
            REAL, INTENT( OUT )         :: VDEMIS( :,:,:,: )
         END SUBROUTINE RDEMIS_NR
         SUBROUTINE RDEMIS_TR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, VDEMIS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, EMISLYRS, NSPC_EMIS
            REAL, INTENT( OUT )         :: VDEMIS( :,:,:,: )
         END SUBROUTINE RDEMIS_TR
         SUBROUTINE PA_UPDATE_EMIS ( PNAME, VDEMIS, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, INTENT( IN )           :: VDEMIS( :,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            INTEGER, INTENT( IN )        :: TSTEP( 2 )
         END SUBROUTINE PA_UPDATE_EMIS
         SUBROUTINE PA_UPDATE_DDEP ( PNAME, DDEP, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, INTENT( IN )           :: DDEP( :,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            INTEGER, INTENT( IN )        :: TSTEP( 2 )
         END SUBROUTINE PA_UPDATE_DDEP
         SUBROUTINE CONV_CGRID ( CGRID, JDATE, JTIME, CNGRD )
            IMPLICIT NONE
            REAL, POINTER :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            REAL, INTENT( OUT ) :: CNGRD( :,:,:,: )
         END SUBROUTINE CONV_CGRID
         SUBROUTINE REV_CGRID ( CNGRD, JDATE, JTIME, CGRID )
            IMPLICIT NONE
            REAL, INTENT( IN ) :: CNGRD( :,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            REAL, POINTER :: CGRID( :,:,:,: )
         END SUBROUTINE REV_CGRID
         SUBROUTINE EDYINTB ( EDDYV, DT, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            REAL, INTENT( OUT )         :: EDDYV( :,:,: )
            REAL, INTENT( OUT )         :: DT( :,: )
            INTEGER, INTENT( IN )       :: JDATE, JTIME, TSTEP
         END SUBROUTINE EDYINTB
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

C for emissions (form COORD.EXT) .......................................

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD ! in m.
            DX2 = DG2M * YCELL_GD
     &          * COS( PI180*( YORIG_GD + YCELL_GD * FLOAT( GL_NROWS/2 ))) ! in m.
            ELSE
            DX1 = XCELL_GD        ! in m.
            DX2 = YCELL_GD        ! in m.
            END IF

C create global maps

         CALL VDIFF_MAP ( DF2EM, DF2DV, DD2DV, DEPV_MAP, DIFF_MAP, DDEP_SPC,
     &                    DV2DF )

C set vertical layer definitions from COORD.EXT

         ALLOCATE ( RDX3F( NLAYS ), STAT = ALLOCSTAT )
         ALLOCATE ( RDX3M( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating RDX3F or RDX3M'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         DO L = 1, NLAYS
            DX3F( L )  = X3FACE_GD( L ) - X3FACE_GD( L-1 )
            RDX3F( L ) = 1.0 / DX3F( L )
            X3M( L ) = 0.5 * ( X3FACE_GD( L ) + X3FACE_GD( L-1 ) )
            END DO

         DO L = 1, NLAYS - 1
            RDX3M( L ) = 1.0 / ( X3M( L+1 ) - X3M( L ) )
            END DO
         RDX3M( NLAYS ) = 0.0

C set molecular weights

         S = 0
         DO V = 1, N_GC_DEPV
            S = S + 1
            DD_CONV( S ) = CNVTD * GC_MOLWT( GC_DEPV_MAP( V ) )
            END DO 
         DO V = 1, N_AE_DEPV
            S = S + 1
            IF ( AE_SPC( AE_DEPV_MAP( V ) )( 1:3 ) .EQ. 'NUM' ) THEN
!              DD_CONV( S ) = M2PHA ! irrelevant, since not deposited
               DD_CONV( S ) = CNVTD * AVO * 1.0E+03    ! --> #/Ha
               ELSE IF ( AE_SPC( AE_DEPV_MAP( V ) )( 1:3 ) .EQ. 'SRF' ) THEN
!              DD_CONV( S ) = M2PHA ! irrelevant, since not deposited
               DD_CONV( S ) = M2PHA * 1.0E+03 / MWAIR  ! --> M**2/Ha
               ELSE
!              DD_CONV( S ) = M2PHA / GPKG / MGPG
               DD_CONV( S ) = CNVTD * AE_MOLWT( AE_DEPV_MAP( V ) )
               END IF
            END DO
         DO V = 1, N_NR_DEPV
            S = S + 1
            DD_CONV( S ) = CNVTD * NR_MOLWT( NR_DEPV_MAP( V ) )
            END DO
         DO V = 1, N_TR_DEPV
            S = S + 1
            DD_CONV( S ) = CNVTD * TR_MOLWT( TR_DEPV_MAP( V ) )
            END DO

C Open the met files

         CALL OPMET ( JDATE, JTIME, CONVPA )
 
C Open Emissions files

         CALL OPEMIS ( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM, EMISLYRS )


         ELAYS = MIN ( EMISLYRS, NLAYS - 1 )

C Set output file characteristics based on COORD.EXT and open the dry dep file

         IF ( MYPE .EQ. 0 ) CALL OPDDEP ( JDATE, JTIME, TSTEP( 1 ), N_SPC_DDEP )

C Get sea-salt-emission diagnostic file flag.

         SSEMDIAG = .FALSE.         ! default
         VARDESC = 'Flag for writing the sea-salt-emission diagnostic file'
         SSEMDIAG = ENVYN( CTM_SSEMDIAG, VARDESC, SSEMDIAG, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG =
     &          'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         END IF

C Open sea-salt emission file if running the AE4 or AE5 aerosol mechanism

         IF ( ( INDEX ( MECHNAME, 'AE3' ) .EQ. 0 ) .AND. SSEMDIAG ) THEN
            IF ( MYPE .EQ. 0 ) CALL OPSSEMIS ( JDATE, JTIME, TSTEP( 1 ) )
         END IF

C Allocate and initialize dry deposition array

!        ALLOCATE ( DDEP( MY_NCOLS,MY_NROWS,N_SPC_DEPV ), STAT = ALLOCSTAT )
         ALLOCATE ( DDEP( N_SPC_DEPV,MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DDEP'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         DDEP = 0.0

         IF ( N_AE_SPC .GT. 0 ) THEN
            WRITE( LOGDEV,'( /5X, A )' ) 'Aerosol Emissions Processing in '
     &                                // 'Vertical diffusion ...'
            NEMIS_AE = N_AE_EMIS  ! from AE_EMIS.EXT
            ELSE
            NEMIS_AE = 0
            END IF

         N_SPC_EMIS = NEMIS
     &              + NEMIS_AE
     &              + N_NR_EMIS
     &              + N_TR_EMIS

         ALLOCATE ( VDEMIS ( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         IF ( N_AE_SPC .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_AE( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'VDEMIS_AE memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF

         IF ( N_NR_EMIS .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_NR( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'VDEMIS_NR memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF

         IF ( EM_TRAC .AND. N_TR_EMIS .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_TR( N_TR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'VDEMIS_TR memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF

         IF ( LIPR ) THEN
            ALLOCATE ( EMIS_PA( MY_NCOLS,MY_NROWS,ELAYS,N_SPC_EMIS+1 ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'EMIS_PA memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
!           ALLOCATE ( DDEP_PA( MY_NCOLS,MY_NROWS,N_SPC_DDEP ),
            ALLOCATE ( DDEP_PA( MY_NCOLS,MY_NROWS,N_SPC_DEPV ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'DDEP_PA memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF

C combined gas emssions conversion factor

         CNVTE = CMLMR * CONVPA * CONVEM * MWAIR / ( DX1 * DX2 )

         N_SPC_CGRID = SIZE ( CGRID,4 )

         ALLOCATE ( CNGRD( N_SPC_CGRID,NLAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CNGRD'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         CNGRD = 0.0

         END IF          !  if Firstime

      MDATE = JDATE
      MTIME = JTIME
      MSTEP = TIME2SEC( TSTEP( 2 ) )
      DTSEC = FLOAT( MSTEP )
      CALL NEXTIME ( MDATE, MTIME, SEC2TIME( MSTEP / 2 ) )

C read & interpolate met data

      CALL RDMET ( MDATE, MTIME, RDEPVHT, RJACM, RVJACMF, RRHOJ, DENS1 )

C read & interpolate deposition velocities

      CALL RDDEPV ( MDATE, MTIME, TSTEP( 2 ), CGRID, DEPV )

C Initialize deposition velocities for nondeposited species to zero

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DEPV( N_SPC_DEPV+1,C,R ) = 0.0 ! accounts for dry dep. species
            END DO                         ! names as a subset of the
         END DO                            ! vert. diffused species list

C read & interpolate emissions (create VDEMIS in the species class order)

      VDEMIS = 0.0

      CALL RDEMIS_GC ( MDATE, MTIME, ELAYS, NEMIS, VDEMIS )

C reactive gases (conversion to ppmv/s) VDEMIS in this order from RDEMIS

      STRT = 1
      FINI = NEMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R ) * CNVTR
                  END DO
               END DO
            END DO
         END DO

C aerosol emissions - all units conversions done in RDEMIS_AE for aerosols

      IF ( N_AE_SPC .GT. 0 ) THEN
C        RDEMIS_AE in f90 module AERO_EMIS
         CALL RDEMIS_AE ( MDATE, MTIME, TSTEP, ELAYS, RJACM, VDEMIS, VDEMIS_AE )
         END IF

      STRT = NEMIS + 1
      FINI = NEMIS + NEMIS_AE
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS( V,L,C,R ) = VDEMIS_AE( S,L,C,R )
                  END DO
               END DO
            END DO
         END DO

C non-reactive gases (conversion to ppmv/s) VDEMIS in this order from RDEMIS

      IF ( N_NR_EMIS .GT. 0 ) THEN
         CALL RDEMIS_NR ( MDATE, MTIME, ELAYS, N_NR_EMIS, VDEMIS_NR )
         END IF

      STRT = NEMIS + NEMIS_AE + 1
      FINI = NEMIS + NEMIS_AE + N_NR_EMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS( V,L,C,R ) = VDEMIS_NR( S,L,C,R ) * CNVTR
                  END DO
               END DO
            END DO
         END DO

C tracer gases (conversion to ppmv/s)

      IF ( EM_TRAC .AND. N_TR_EMIS .GT. 0 ) THEN
         CALL RDEMIS_TR ( MDATE, MTIME, ELAYS, N_TR_EMIS, VDEMIS_TR )
         END IF

      STRT = NEMIS + NEMIS_AE + N_NR_EMIS + 1
      FINI = NEMIS + NEMIS_AE + N_NR_EMIS + N_TR_EMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS( V,L,C,R ) = VDEMIS_TR( S,L,C,R ) * CNVTR
                  END DO
               END DO
            END DO
         END DO

C zero out emissions values for species not included in diffused list

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               VDEMIS( N_SPC_EMIS+1,L,C,R ) = 0.0  ! accounts for emissions
               END DO                              ! species names as a subset
            END DO                                 ! of the vert. diffused
         END DO                                    ! species list

      IF ( LIPR ) THEN
         DO S = 1, N_SPC_EMIS+1
            DO L = 1, ELAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     EMIS_PA( C,R,L,S ) = VDEMIS( S,L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         CALL PA_UPDATE_EMIS ( 'VDIF', EMIS_PA, JDATE, JTIME, TSTEP )
         END IF

      CALL EDYINTB ( EDDYV, DT, JDATE, JTIME, TSTEP( 2 ) )

C EDDYV returned = Kz, where Kz is in m**2/sec

C get number of steps based on eddy time

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            NSTEPS( C,R ) = INT ( DTSEC / DT( C,R ) ) + 1
            DT( C,R ) = DTSEC / FLOAT( NSTEPS( C,R ) )
            END DO
         END DO

      IF ( EDDY_STATS ) THEN
         NSTP_AVG = 0.0
         DT_AVG = 0.0
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               NSTP_AVG = NSTP_AVG + FLOAT( NSTEPS( C,R ) )
               DT_AVG = DT_AVG + DT( C,R )
               END DO
            END DO
         WRITE( LOGDEV,1003 ) SUBST_GLOBAL_SUM( DT_AVG ) /
     &                        FLOAT( GL_NCOLS*GL_NROWS ),
     &                        SUBST_GLOBAL_SUM( NSTP_AVG ) /
     &                        FLOAT( GL_NCOLS*GL_NROWS )
1003     FORMAT( 5X, 'Grid-average eddy DT (sec) & no. of integration steps:',
     &           2(1PE10.2) )
         END IF

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DELT = DT( C,R )
            DO L = 1, NLAYS
               SEDDY( L,C,R ) = EDDYV( C,R,L )
     &                        * RVJACMF( C,R,L ) * RDX3M( L ) * DELT
               END DO
            END DO
         END DO

C Convert non-molar mixing ratio species and re-order CGRID

      CALL CONV_CGRID ( CGRID, MDATE, MTIME, CNGRD )

C compute species-independent deposition velocity combination conversion
C factors DVFACB, DVFACD, and ddep accumulation factor

      IF ( THETA .GT. 0.0 ) THEN
         THRAT = THBAR / THETA
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DVFACD( C,R ) = THETA * DT( C,R ) * RDEPVHT( C,R )
               DVFACB( C,R ) = THRAT * DVFACD( C,R )
               END DO
            END DO
         ELSE
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DVFACD( C,R ) = 0.0
               DVFACB( C,R ) = DT( C,R ) * RDEPVHT( C,R )
               END DO
            END DO
         END IF

      IF ( WSTEP .EQ. 0 ) DDEP = 0.0

      DO 345 R = 1, MY_NROWS
         DO 344 C = 1, MY_NCOLS
            DELT = DT( C,R )
            DTDENS1 = DELT * DENS1( C,R )

            DO L = 1, NLAYS
               DO V = 1, N_SPC_DIFF
                  CONC( V,L ) = CNGRD( DIFF_MAP( V ),L,C,R )
!                 CONC( V,L ) = CGRID( C,R,L,DIFF_MAP( V ) )
                  END DO
               END DO

            DO V = 1, N_SPC_DEPV
!              DDBF( V ) = DDEP( C,R,V )
               DDBF( V ) = DDEP( V,C,R )
               DD_FAC( V ) = DTDENS1 * DD_CONV( V ) * DEPV( V,C,R )
               END DO

            EMIS = 0.0
            DO L = 1, ELAYS
               DO V = 1, N_SPC_DIFF
                  EMIS( V,L ) = VDEMIS( DF2EM( V ),L,C,R ) * DELT
                  END DO
               END DO

            DO V = 1, N_SPC_DIFF
               DVD1( V ) = DEPV( DF2DV( V ),C,R ) * DVFACD( C,R )
               DVB1( V ) = DEPV( DF2DV( V ),C,R ) * DVFACB( C,R )
               END DO

            PK( 1 ) = SEDDY( 1,C,R ) * RDX3F( 1 )
            UI( 1 ) = -THETA * PK( 1 )
            DO L = 2, NLAYS
               PK( L )   = SEDDY( L,  C,R ) * RDX3F( L )
               QK( L )   = SEDDY( L-1,C,R ) * RDX3F( L )
               UI( L )   = -THETA * PK( L )
               LI( L )   = -THETA * QK( L )
               DI( 1,L ) = 1.0 - UI( L ) - LI( L )  ! UI(NLAYS) = 0
               END DO

            DO V = 1, N_SPC_DIFF
               DI( V,1 ) = 1.0 - UI( 1 ) + DVD1( V )
               DI( V,NLAYS ) = DI( 1,NLAYS )
               END DO

            DO L = 2, NLAYS - 1
               DO V = 1, N_SPC_DIFF
                  DI( V,L ) = DI( 1,L )
                  END DO
               END DO

            DO 301 N = 1, NSTEPS( C,R )

               DO V = 1, N_SPC_DEPV
                  DDBF( V ) = DDBF( V ) + THBAR * DD_FAC( V )
     &                      * CONC( DV2DF( V ),1 )
                  END DO

               L = 1
               DO V = 1, N_SPC_DIFF
                  BI( V,L ) = CONC( V,L )
     &                      + THBAR * PK( L )
     &                      * ( CONC( V,L+1 ) - CONC( V,L ) )
     &                      - DVB1( V ) * CONC( V,L )
     &                      + EMIS( V,L )
                  END DO

               L = NLAYS
               DO V = 1, N_SPC_DIFF
                  BI( V,L ) = CONC( V,L )
     &                      - THBAR * QK( L )
     &                      * ( CONC( V,L ) - CONC( V,L-1 ) )
                  END DO

               DO L = 2, NLAYS - 1
                  DO V = 1, N_SPC_DIFF
                     BI( V,L ) = CONC( V,L )
     &                         + THBAR
     &                         * ( PK( L )
     &                         * ( CONC( V,L+1 ) - CONC( V,L ) )
     &                         -   QK( L )
     &                         * ( CONC( V,L   ) - CONC( V,L-1 ) ) )
     &                         + EMIS( V,L )
                     END DO
                  END DO

               CALL TRIDIAG ( LI, DI, UI, BI, XI )

               DO L = 1, NLAYS
                  DO V = 1, N_SPC_DIFF
!                    S = DIFF_MAP( V )
!                    IF ( XI( V,L ) .LT. 0.0 ) THEN
!                       WRITE( *,* ) C,R,L,S, XI( V,L )
!                       END IF
                     CONC( V,L ) = XI( V,L )
                     END DO
                  END DO

               DO V = 1, N_SPC_DEPV
                  DDBF( V ) = DDBF( V ) + THETA * DD_FAC( V )
     &                      * CONC( DV2DF( V ),1 )
                  END DO

301            CONTINUE      !  end time steps loop

            DO L = 1, NLAYS
               DO V = 1, N_SPC_DIFF
                  CNGRD( DIFF_MAP( V ),L,C,R ) = CONC( V,L )
!                 CGRID( C,R,L,DIFF_MAP( V ) ) = CONC( V,L )
                  END DO
               END DO
            DO V = 1, N_SPC_DEPV
!              DDEP( C,R,V ) = DDBF( V )
               DDEP( V,C,R ) = DDBF( V )
               END DO

344         CONTINUE         !  end loop on col C
345      CONTINUE            !  end loop on row R

C Revert non-molar mixing ratio species and re-order CGRID

      CALL REV_CGRID ( CNGRD, MDATE, MTIME, CGRID )

C If last call this hour:  write accumulated depositions:

      WSTEP = WSTEP + TIME2SEC( TSTEP( 2 ) )
      IF ( WSTEP .GE. TIME2SEC( TSTEP( 1 ) ) ) THEN
         MDATE = JDATE
         MTIME = JTIME
         CALL NEXTIME( MDATE, MTIME, TSTEP( 2 ) )
         WSTEP = 0

         DO V = 1, N_SPC_DDEP
            S = DD2DV( V )
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  WRDD( C,R ) = DDEP( S,C,R )
                  END DO
               END DO

            IF ( .NOT. WRITE3( CTM_DRY_DEP_1, DDEP_SPC( V ),
     &                 MDATE, MTIME, WRDD ) ) THEN
               XMSG = 'Could not write ' // CTM_DRY_DEP_1 // ' file'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF

            END DO

         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', CTM_DRY_DEP_1,
     &         'for date and time', MDATE, MTIME

         IF ( LIPR ) THEN
!           DO V = 1, N_SPC_DDEP
            DO V = 1, N_SPC_DEPV
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     DDEP_PA( C,R,V ) = DDEP( V,C,R )
                     END DO
                  END DO
               END DO
            CALL PA_UPDATE_DDEP ( 'VDIF', DDEP_PA, JDATE, JTIME, TSTEP )
            END IF

C re-set dry deposition array to zero

         DDEP = 0.0

         END IF

      RETURN
      END
