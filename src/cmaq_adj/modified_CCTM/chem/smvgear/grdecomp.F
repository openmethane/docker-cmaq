
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/smvgear/grdecomp.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)grdecomp.F        1.1 /project/mod3/CMAQ/src/chem/smvgear/SCCS/s.grdecomp.  F 07 Jul 1997 12:45:23

       SUBROUTINE DECOMP

C***********************************************************************
C
C  FUNCTION:  To decompose the matrix [P] into lower- and upper
C             triangular form to facilitate solution of the set of 
C             linear simultaneous equations of the form [A]{x}={b}.
C
C  PRECONDITIONS: None
C                                                                     
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, June, 1995.
C                      Based on  the code originally developed by 
C                      M. Jacobson, (Atm. Env., Vol 28, No 2, 1994).
C                    Revised 3/14/96 by Jerry Gipson to conform to
C                      the Models-3 minimum IOV configuration
C                    Revised December 1996 by Jerry Gipson to conform
C                      to the Models-3 interim CTM that includes emissions
C                      in chemistry.
C                    Modified June, 1997 by Jerry Gipson to be consistent
C                      with beta CTM
C                    Modified September, 1997 by Jerry Gipson to be
C                      consistent with the targetted CTM
C                    16 Aug 01 J.Young: Use HGRD_DEFN
C                    31 Jan 05 J.Young: get BLKSIZE from dyn alloc horizontal
C                    & vertical domain specifications module (GRID_CONF)
C***********************************************************************

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE
      
C..INCLUDES:
!     INCLUDE SUBST_HGRD_ID
!     INCLUDE SUBST_VGRD_ID 
!     INCLUDE SUBST_BLKPRM
      INCLUDE SUBST_GC_SPC
      INCLUDE SUBST_RXCMMN 
      INCLUDE 'GRPARMS.EXT'
      INCLUDE 'GRVARS1.EXT'
      INCLUDE 'GRVARS2.EXT'
      INCLUDE 'GRVARS3.EXT'
      
C..ARGUMENTS: None

C..PARAMETERS: None

C..EXTERNAL FUNCTIONS: None

C..SAVED LOCAL VARIABLES: None

C..SCRATCH LOCAL VARIABLES:
      INTEGER IAR             ! Pointer to diagonal terms
      INTEGER IC              ! Loop index for ops in decomp loop 1
      INTEGER IDLO            ! Start index for decomp loop 1
      INTEGER IDHI            ! End index for decomp loop 1
      INTEGER IJ0             ! Pointer to ij term 1 in decomp loop 1
      INTEGER IJ1             ! Pointer to ij term 2 in decomp loop 1
      INTEGER IJA             ! Pointer to ij term 1 in decomp loop 2
      INTEGER IJB             ! Pointer to ij term 2 in decomp loop 2
      INTEGER IK0             ! Pointer to ik term 1 in decomp loop 1
      INTEGER IK1             ! Pointer to ik term 2 in decomp loop 1
      INTEGER J               ! Loop index for number of species
      INTEGER JC              ! Loop index for ops in decomp loop 2
      INTEGER JHI1            ! End index for 2-term decomp loop 2
      INTEGER JHI2            ! End index for 1-term decomp loop 2
      INTEGER JLO1            ! Start index for 2-term decomp loop 2
      INTEGER JLO2            ! Start index for 1-term decomp loop 2
      INTEGER KJ0             ! Pointer to kj term 1 in decomp loop 1
      INTEGER KJ1             ! Pointer to kj term 2 in decomp loop 1
      INTEGER NCELL           ! Loop index for number of cells      
c***********************************************************************      
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  First loop of L-U decomposition 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      JHI2 = JZLO( NCSP )
      DO 100 J = 1, ISCHAN
         IDLO = IDEC1LO( J, NCSP )
         IDHI = IDEC1HI( J, NCSP )
         DO 20 IC = IDLO, IDHI
            IJ0 = IJDECA( IC )
            IJ1 = IJDECB( IC )
            IK0 = IKDECA( IC )
            IK1 = IKDECB( IC )
            KJ0 = KJDECA( IC )
            KJ1 = KJDECB( IC )
            DO 20  NCELL = 1, NUMCELLS
               CC2( NCELL, IJ0 ) = CC2( NCELL, IJ0 ) - 
     &                             CC2( NCELL, IK0 ) * CC2( NCELL, KJ0 )
               CC2( NCELL, IJ1 ) = CC2( NCELL, IJ1 ) - 
     &                             CC2( NCELL, IK1 ) * CC2( NCELL, KJ1 )
   20    CONTINUE
    
c...vdiag = 1 / current diagonal term of the decomposed matrix
         IAR = JARRAYPT( J, J, NCSP )
         DO 40 NCELL = 1, NUMCELLS
            VDIAG( NCELL, J )  = 1.0 / CC2( NCELL, IAR )
   40    CONTINUE
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Second loop of decompostion. The elements of L are divided by the
c  diagonal element, and the preocess is divided into parts to improve
c  vectorization.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         JLO1 = JHI2 + 1
         JHI1 = JHI2 + JHIZ1( J, NCSP ) 
         JLO2 = JHI1 + 1 
         JHI2 = JHI1 + JHIZ2( J, NCSP )
         
c...determine 2 terms at a time
         DO 60 JC = JLO1, JHI1
            IJA = JZEROA( JC )
            IJB = JZEROB( JC )
            DO 60 NCELL = 1, NUMCELLS
               CC2( NCELL, IJA ) = CC2( NCELL, IJA ) * VDIAG( NCELL, J )  
               CC2( NCELL, IJB ) = CC2( NCELL, IJB ) * VDIAG( NCELL, J )  
  60     CONTINUE
 
c...determine 1 term at a time 
         DO 80 JC = JLO2, JHI2 
            IJA = JZEROA( JC )
            DO 80 NCELL = 1, NUMCELLS
               CC2( NCELL, IJA ) = CC2( NCELL, IJA ) * VDIAG( NCELL, J )  
   80    CONTINUE
  100 CONTINUE
      RETURN
      END
