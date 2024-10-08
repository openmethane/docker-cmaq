
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/plsprd.f,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      SUBROUTINE PLSPRD( DTHDZ, ZF, KZ, CEFSTK, HTMIX, PLTOP, PLBOT )
       SUBROUTINE PLSPRD( DTHDZ, ZF, KZ, CEFSTK,        PLTOP, PLBOT )
     
C-----------------------------------------------------------------------
C Description:  
C   Calculates the initial vertical spread of a plume; modified
C   from Gillani's model.
 
C Preconditions:
 
C Subroutines and Functions Called:
 
C Revision History:
C   Initially provided by J. Godowitch ( EPA, 9/03 )
 
C-----------------------------------------------------------------------
C Modified from:

C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id: plsprd.f,v 1.1 2011/02/22 01:10:45 ppercell Exp $
C COPYRIGHT (C) 2003, MCNC Environmental Modeling Center
C All Rights Reserved
C See file COPYRIGHT for conditions of use.
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C smoke@emc.mcnc.org
C Pathname: $Source: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/plsprd.f,v $
C Last updated: $Date: 2011/02/22 01:10:45 $ 
C-----------------------------------------------------------------------

       IMPLICIT NONE
       
C Arguments:

       REAL,    INTENT ( IN ) :: DTHDZ( KZ ) ! potential temperature lapse rate (K/m)
       REAL,    INTENT ( IN ) :: ZF( 0:KZ )  ! full-layer heights (m)
       INTEGER, INTENT ( IN ) :: KZ          ! number of emissions layers
       REAL,    INTENT ( IN ) :: CEFSTK      ! effective stack height (m)
!      REAL,    INTENT ( IN ) :: HTMIX       ! mixing height (m)
       REAL,    INTENT( OUT ) :: PLTOP       ! plume top (m)
       REAL,    INTENT( OUT ) :: PLBOT       ! plume bottom (m)
       
C Parameters:
       REAL, PARAMETER :: SZ0FAC = 3.545     ! factor used to derive plume depth
       REAL, PARAMETER :: SPRFAC = 15.       ! empirical coefficient for vertical spread
       REAL, PARAMETER :: GAMMA  = -0.0098   ! adiabatic lapse rate (K/m)
       
C Local Variables:
       INTEGER K
       
       REAL    SIGZ0
       REAL    DTDZ
       REAL    DPTH
       
C-----------------------------------------------------------------------

C Get ambient temperature above plume rise height (effective stack height)
      K = 0
      DO 
         K = K + 1
         IF ( K .EQ. KZ .OR. CEFSTK .LE. ZF( K ) ) EXIT
      END DO
      DTDZ  = DTHDZ( K ) + GAMMA
       
C Compute initial vertical spread
      SIGZ0 = MAX( 10.0, SPRFAC * EXP( -117.0 * DTDZ ) )
      DPTH  = SZ0FAC * SIGZ0
       
C Compute plume top and bottom heights; plume is either completely
C within or outside mixing layer
      PLTOP = CEFSTK + DPTH / 2.0
      PLBOT = CEFSTK - DPTH / 2.0
       
C Make sure plume bottom is at least zero
      PLBOT = MAX( 0.0, PLBOT )
       
C Make sure that plume top and bottom height < ZF (KZ) .       
      PLTOP = MIN( ZF( KZ ), PLTOP )
      PLBOT = MIN( ZF( KZ )-1.0, PLBOT )

      RETURN
       
      END SUBROUTINE PLSPRD
