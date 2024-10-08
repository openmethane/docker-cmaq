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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ros3/rbfeval.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE RBFEVAL( NCSP, YIN, YDOT )

C***********************************************************************
C
C  FUNCTION:  To compute YDOT = dc/dt for each species. YDOT is the
C             net rate of change in species concentrations resulting
C             from chemical production minus chemical loss.
C
C  PRECONDITIONS: None
C                                                                     
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, August, 2004
C                    Based on the SMVGEAR code originally developed by 
C                    M. Jacobson, (Atm. Env., Vol 28, No 2, 1994).
C                    31 Jan 05 J.Young: get BLKSIZE from dyn alloc horizontal
C                    & vertical domain specifications module (GRID_CONF)
C
C***********************************************************************

      USE GRID_CONF            ! horizontal & vertical domain specifications
      USE RBDATA               ! ROS3 solver data

      IMPLICIT NONE

C..INCLUDES:
      INCLUDE SUBST_GC_SPC     ! Gas chem species names and MWs
      INCLUDE SUBST_RXCMMN     ! Mechanism reaction common block

      
C..ARGUMENTS:
      INTEGER NCSP                         ! Index of chem mech to use
                                           ! 1=gas/day, 2=gas/night

      REAL( 8 ) ::  YIN(  BLKSIZE, * )     ! Species concs, ppm
      REAL( 8 ) ::  YDOT( BLKSIZE, * )     ! Species rates of change, ppm/min


C..PARAMETERS: None


C..EXTERNAL FUNCTIONS: None


C..SAVED LOCAL VARIABLES: None


C..SCRATCH LOCAL VARIABLES:
      INTEGER ISP            ! Loop index for species
      INTEGER ISP1,ISP2,ISP3 ! Pointers to species numbers
      INTEGER NCELL          ! Loop index for number of cells
      INTEGER NP             ! Loop index for number of products
      INTEGER NR             ! Loop index for number of reactants
      INTEGER NRK            ! Reaction number
      INTEGER NRX            ! Loop index for number of reactions

      REAL( 8 )    ::   RXRAT( BLKSIZE, NRXNS )   ! Computed rxn rates

C***********************************************************************      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Initialize dc/dt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO ISP = 1, ISCHAN
         DO NCELL = 1, NUMCELLS
            YDOT( NCELL, ISP ) = 0.0
         ENDDO
      ENDDO
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Loop over reactions to calculate dc/dt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO 100 NRX = 1, NUSERAT( NCSP )

        NRK = NKUSERAT( NRX, NCSP )

c..Calculate reaction rate
        IF( NREACT( NRK ) .EQ. 1 ) THEN
            ISP1 = IRM2( NRK, 1, NCS )
            DO NCELL = 1, NUMCELLS
               RXRAT( NCELL, NRK ) = RKI( NCELL, NRK ) * 
     &                               YIN( NCELL, ISP1 )
            ENDDO
         ELSEIF( NREACT( NRK ) .EQ. 2 ) THEN
            ISP1 = IRM2( NRK, 1, NCS )
            ISP2 = IRM2( NRK, 2, NCS )
            DO NCELL = 1, NUMCELLS
               RXRAT( NCELL, NRK ) = RKI( NCELL, NRK ) * 
     &                               YIN( NCELL, ISP1 ) *
     &                               YIN( NCELL, ISP2 )
            ENDDO
         ELSEIF( NREACT( NRK ) .EQ. 3 ) THEN
            ISP1 = IRM2( NRK, 1, NCS )
            ISP2 = IRM2( NRK, 2, NCS )
            ISP3 = IRM2( NRK, 3, NCS )
            DO NCELL = 1, NUMCELLS
               RXRAT( NCELL, NRK ) = RKI( NCELL, NRK ) * 
     &                               YIN( NCELL, ISP1 ) *
     &                               YIN( NCELL, ISP2 ) *
     &                               YIN( NCELL, ISP3 )
            ENDDO 
         ELSEIF( NREACT( NRK ) .EQ. 0 ) THEN
            DO NCELL = 1, NUMCELLS
               RXRAT( NCELL, NRK ) = RKI( NCELL, NRK )
            ENDDO
         ENDIF
         
c..Subtract loss terms from dc/dt for this reaction 
         DO NR = 1, NREACT( NRK )
            ISP1 = IRM2( NRK, NR, NCS )
            DO NCELL = 1, NUMCELLS
               YDOT( NCELL, ISP1 ) = YDOT( NCELL, ISP1 ) - 
     &                                RXRAT( NCELL, NRK )
            ENDDO
         ENDDO
  
c..Add production terms to dc/dt for this reaction
          DO NP = 1, NPRDCT( NRK )
            ISP1 = IRM2( NRK, NP + 3, NCS )
            DO NCELL = 1, NUMCELLS
             YDOT( NCELL, ISP1 ) = YDOT( NCELL, ISP1 ) + 
     &                              SC( NRK, NP ) * RXRAT( NCELL, NRK )
            ENDDO
         ENDDO

  100 CONTINUE               ! END LOOP OVER REACTIONS
  

      RETURN
      END

