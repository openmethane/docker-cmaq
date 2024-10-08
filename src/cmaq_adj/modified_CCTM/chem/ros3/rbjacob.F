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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/ros3/rbjacob.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE RBJACOB( NCSP, YIN )

C***********************************************************************
C
C  FUNCTION: To compute [J] where J is the Jacobian matrix,
C            (i.e., Jij = d[dCi/dt]/dCj )
C
C  PRECONDITIONS: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, August, 2004
C                    Based on the SMVGEAR code originally developed by 
C                    M. Jacobson, (Atm. Env., Vol 28, No 2, 1994).
C
C                    31 Jan 05 J.Young: get BLKSIZE from dyn alloc horizontal
C                    & vertical domain specifications module (GRID_CONF)
C
C***********************************************************************

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE RBDATA                ! ROS3 solver data

      IMPLICIT NONE
      

C..INCLUDES:
      INCLUDE SUBST_GC_SPC      ! Gas chem species names and MWs
      INCLUDE SUBST_RXCMMN      ! Mechanism reaction common block

      
C..ARGUMENTS:
      INTEGER   :: NCSP                    ! Index of chem mech to use
                                           ! 1=gas/day, 2=gas/night

      REAL( 8 )    :: YIN( BLKSIZE, * )    ! Species concs, ppm


C..PARAMETERS: None


C..EXTERNAL FUNCTIONS: None


C..SAVED LOCAL VARIABLES: None


C..SCRATCH LOCAL VARIABLES:
      INTEGER IALP           ! Pointer to location of PD term in EXPLIC
      INTEGER IAR            ! Loop index for non-zero entries in [P]
      INTEGER IARP           ! Pointer to location of PD term in [P]
      INTEGER IARRY          ! Pointer to end of [P] entries
      INTEGER ISCP           ! Pointer to stoichiometric coefficient
      INTEGER ISPC           ! Loop index for species
      INTEGER JR1, JR2, JR3  ! Pointer to reactant species conc.
      INTEGER NCELL          ! Loop index for number of cells
      INTEGER NL             ! Loop index for loss PD terms
      INTEGER NLD            ! Number of loss PD terms for each rxn.
      INTEGER NP             ! Loop index for prod PD terms
      INTEGER NPD            ! Number of prod PD terms for each rxn.
      INTEGER NRK            ! Reaction number
      INTEGER NRX            ! Loop index for number of reactions
      INTEGER NONDIAG        ! Pointer to end of off-diagonal entries
      INTEGER NONDIAG1       ! Pointer to start of diagonal entries
      
      REAL( 8 ) :: CR2                   ! Temporary product for 3 reactant rxns
      REAL( 8 ) :: FRACN                 ! Stoichiometric coeff. times b*h
      REAL( 8 ) :: EXPLIC( BLKSIZE, 3 )  ! Reaction partial derivative terms

C***********************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Zero out Jacobian ( stored in sparse matrix array cc2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IARRY = IARRAY( NCSP ) 
      NONDIAG = IARRY - ISCHAN  
      NONDIAG1 = NONDIAG + 1
      DO 20 IAR = 1, NONDIAG
         DO 20 NCELL = 1, NUMCELLS
            CC2( NCELL, IAR ) = 0.
  20  CONTINUE
      DO 40 IAR = NONDIAG1, IARRY
      DO 40 NCELL = 1, NUMCELLS
         CC2( NCELL, IAR ) = 0.
  40  CONTINUE
  
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Loop over reaction rates adding partial derivative terms; EXPLIC
c  holds the PD terms according to number of reactants
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 240 NRX = 1, NUSERAT( NCSP )

         NRK = NKUSERAT( NRX, NCSP )
         
c...partial derivative term for reactions with 1 reactant
         IF( NREACT( NRK ) .EQ. 1 ) THEN
            DO 100 NCELL = 1, NUMCELLS
               EXPLIC( NCELL, 1 )  = RKI( NCELL, NRK ) 
  100       CONTINUE
  
c...partial derivative terms for reactions with 2 reactants
         ELSEIF( NREACT( NRK ) .EQ. 2 ) THEN
            JR1 = IRM2( NRK, 1, NCS )
            JR2 = IRM2( NRK, 2, NCS )
            DO 120 NCELL  = 1, NUMCELLS
               EXPLIC( NCELL, 1 )  = RKI(   NCELL, NRK ) * 
     &                               YIN( NCELL, JR2 )
               EXPLIC( NCELL, 2 )  = RKI(   NCELL, NRK ) * 
     &                               YIN( NCELL, JR1 )
 120        CONTINUE
 
c.....partial derivative terms for reactions with 3 reactants
         ELSEIF( NREACT( NRK ).EQ.3 ) THEN
            JR1 = IRM2( NRK, 1, NCS )
            JR2 = IRM2( NRK, 2, NCS )
            JR3 = IRM2( NRK, 3, NCS )
            DO 140 NCELL = 1, NUMCELLS
               CR2 = RKI( NCELL, NRK ) * YIN( NCELL, JR2 )
               EXPLIC( NCELL, 1 ) = CR2 * YIN( NCELL, JR3 )
               EXPLIC( NCELL, 2 ) = RKI(   NCELL, NRK ) * 
     &                              YIN( NCELL, JR1 ) * 
     &                              YIN( NCELL, JR3 ) 
               EXPLIC( NCELL, 3 ) = CR2 * YIN( NCELL, JR1 )
  140       CONTINUE
         ENDIF
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Add PD terms to [J] for this reaction
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...loss terms
         NLD = NDERIVL( NRK, NCSP )         
         DO 160 NL = 1, NLD

            IARP = JARRL( NRK, NL, NCSP )
            IALP = JLIAL( NRK, NL, NCSP )
            DO 160 NCELL = 1, NUMCELLS
               CC2( NCELL, IARP ) = CC2( NCELL, IARP ) -  
     &                              EXPLIC( NCELL, IALP ) 

  160       CONTINUE    ! End loop over loss terms
  
c...production terms with stoichiomteric coeff EQ 1.0 and NE 1.0
         NPD = NDERIVP( NRK, NCSP )
         DO 220 NP = 1, NPD

            IARP = JARRP( NRK, NP, NCSP )
            IALP = JPIAL( NRK, NP, NCSP )

c..production terms with unit stoichiometry
            IF( ICOEFF( NRK, NP, NCSP ) .EQ. 0 ) THEN
               DO 180 NCELL = 1, NUMCELLS
                  CC2( NCELL, IARP ) = CC2( NCELL, IARP ) +  
     &                                 EXPLIC( NCELL, IALP ) 
  180          CONTINUE

c..production terms with non-unit stoichiometry
            ELSE

               ISCP = ICOEFF( NRK, NP, NCSP )
               FRACN = SC( NRK, ISCP ) 
               DO 200 NCELL = 1, NUMCELLS
                  CC2( NCELL, IARP ) = CC2( NCELL, IARP ) + FRACN * 
     &                                 EXPLIC( NCELL, IALP ) 
  200          CONTINUE
 
            ENDIF

  220    CONTINUE      ! End loop over production terms

  240 CONTINUE         ! End loop over reactions


      RETURN 
      END
