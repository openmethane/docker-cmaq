
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/smvgear/grsprse.F,v 1.1.1.1 2010/06/14 16:03:02 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)grsprse.F 1.1 /project/mod3/CMAQ/src/chem/smvgear/SCCS/s.grsprse.F 07 Jul 1997 12:45:30

      SUBROUTINE JSPARSE
 
C***********************************************************************      
C
C  FUNCTION: To define array pointers for sparse matrix storage by
C            doing symbolic LU decomposition
C
C  PRECONDITIONS: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, June, 1995. 
C                      Based on  the code originally developed by 
C                      M. Jacobson, (Atm. Env., Vol 28, No 2, 1994).
C                    Revised 3/24/96 By Jerry Gipson to conform to 
C                      the Models-3 minimum IOV configuration
C                    Revised December 1996 by Jerry Gipson to conform
C                      to the Models-3 interim CTM that includes emissions
C                      in chemistry.
C                    Revised April 1997 to distinguish NSPCS from NSPCSD
C                    Revised April 1997 to conform to Models-3 framework
C                    Modified June, 1997 by Jerry Gipson to be consistent
C                      with beta CTM
C                    Modified September, 1997 by Jerry Gipson to be 
C                      consistent with the targeted CTM
C                    16 Aug 01 J.Young: Use HGRD_DEFN
C***********************************************************************

      USE HGRD_DEFN             ! horizontal domain specifications

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
      
C..ARGUMENTS: None

C..PARAMETERS: None

C..EXTERNAL FUNCTIONS: None

C..SAVED LOCAL VARIABLES: 
      LOGICAL, SAVE :: LFIRST = .TRUE. ! Flag for first call to this subroutine

      INTEGER, SAVE :: IFNEVER = 0  ! Flag for counter initialization
      INTEGER, SAVE :: NDLMAX  = 0  ! Max # of PD loss terms in any reaction
      INTEGER, SAVE :: NDPMAX  = 0  ! Max # of PD prod terms in any reaction

C..SCRATCH LOCAL VARIABLES:
      CHARACTER( 16 ) :: PNAME = 'GRSPRSE'     ! Program name
      CHARACTER( 80 ) :: MSG       ! Mesaage text for output to log

      INTEGER EXITSTAT       ! Exit status code
      INTEGER I,J,K,I1,J1,I2 ! Matrix loop indices
      INTEGER IA, IB         ! I,J index holders for decomp loop 2
      INTEGER INEW, JNEW     ! Index for sorted species number
      INTEGER IOLD, JOLD     ! Index for old species number
      INTEGER IPA, KPA       ! I,K index holders for decomp loop 1
      INTEGER IPB, KPB       ! I,K index holders for decomp loop 1
      INTEGER IPROD, JP      ! Species number of a product
      INTEGER IREACT, IR, JR ! Species number of a reactant
      INTEGER ISP, ISP2      ! Species loop indices
      INTEGER JRE, JPR, IRE  ! Indices for nonzero Jacobian entries 
      INTEGER JZ3, JZ4       ! Counter for calcs in backsub groupings
      INTEGER NP, IAP        ! Product loop indices
      INTEGER NR, IAL, JAL   ! Reactant loop indices
      INTEGER IAR            ! Pointer to location of PD term
      INTEGER IARRAY2        ! Final # of matrix entries w/ Sp. Mat
      INTEGER ICB            ! Counter for # of terms in decomp loop 1
      INTEGER ICBSUM         ! Running count of calcs for j index 
                             ! in decomp loop 1
      INTEGER ICCOUNT        ! Two term op count for decomp loop 1
      INTEGER ICNT           ! Total op counter for decomp loop 1
      INTEGER ICNTA          ! op. counter for decomp loop 1 w/ Sp Mat 
      INTEGER ICNTB          ! op. counter for decomp loop 1 w/ Sp Mat
      INTEGER IFSUN          ! Day/night loop index
      INTEGER IJSTEP         ! Number of terms to calc in decomp loops
      INTEGER IMINNEW        ! Index holder for sort routine
      INTEGER IMINOLD        ! Index holder for sort routine
      INTEGER IPORR          ! Species number of a product or reactant      
      INTEGER, SAVE :: IZERO = 0   ! Integer zero
      INTEGER JCB            ! Counter for # of terms in decomp loop 2
      INTEGER JCCOUNT        ! Two term op count for decomp loop 2
      INTEGER JCNT           ! Total op counter for decomp loop 2 
      INTEGER JCNTA          ! op. counter for decomp loop 2 w/o Sp Mat
      INTEGER JCNTB          ! op. counter for decomp loop 2 w/ Sp Mat
      INTEGER JZ             ! Loop index for backsub loops
      INTEGER KA             ! Loop index for decomposition loops
      INTEGER KCNT           ! op. counter for bksub loop 1 w/ Sp. Mat.
      INTEGER KCNTA          ! op. counter for bksub loop 1 w/o Sp Mat
      INTEGER KNTARRAY       ! Final # of matrix entries w/o Sp. Mat
      INTEGER KOUNT0         ! Initial # of matrix entries w/ Sp. Mat
      INTEGER KOUNT0A        ! Initial # of matrix entries w/o Sp. Mat
      INTEGER KZ             ! # of nonzero calcs in backsub loop 1
      INTEGER NCSP           ! Mechanism number NCS+1=day NCS+2=night
      INTEGER NK             ! Reaction number
      INTEGER NLS            ! Number of loss PD terms
      INTEGER NOCHANG        ! Count of number of species not reacting
      INTEGER NPR            ! Number of prod PD terms
      INTEGER NQQ            ! Loop index for Gear order      
      INTEGER NRPP           ! Reactant plus product loop index
      INTEGER NRX            ! Reaction loop index
      INTEGER NU             ! Active reaction count holder
      INTEGER MCNT           ! op. counter for bksub loop 2 w/ Sp. Mat.
      INTEGER MCNTA          ! op. counter for bksub loop 2 w/o Sp. Mat.
      INTEGER MINVALU        ! Current number of PD terms in sort
      INTEGER MZ             ! # of nonzero calcs in backsub loop 2
      INTEGER ICLO(        NCS2 ) ! Pointer to # of ops in decomp loop 1
      INTEGER ISAPORL( N_GC_SPC ) ! Count of PD terms for each species
      INTEGER IZEROK(  MXCOUNT2 ) ! Pointer to decomp loop 1 k index
      INTEGER IZEROI(  MXCOUNT1 ) ! Pointer to decomp loop 1 i index
      INTEGER JCLO(        NCS2 ) ! Pointer to # of ops in decomp loop 2
      INTEGER JZERO(   MXCOUNT1 ) ! Pointer to decomp loop 2 i index
      INTEGER NSPECT(       NCS ) ! Number of species in mechanism ncs
      INTEGER ISPARDER( N_GC_SPC, N_GC_SPC ) ! Indicator of a PD term in the 
                                             ! Jacobian matrix
      INTEGER IZILCH(   N_GC_SPC,     NCS2 ) ! # of nonzero calcs in decomp
                                             ! loop 1
      INTEGER JZILCH(   N_GC_SPC,     NCS2 ) ! # of nonzero calcs in decomp
                                             ! loop 2
      INTEGER LZERO(    N_GC_SPC, N_GC_SPC ) ! Symbolic Jacobian matrix
C***********************************************************************                                             
      DATA PERTST / 
     &  2.0,  4.5,  7.333, 10.42,  13.7,     17.15,     1.0,
     &  3.0,  6.0,  9.167, 12.5,   15.98,     1.0,      1.0,
     &  1.0,  1.0,  0.5,    0.1667, 0.04133,  0.008267, 1.0    /

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Initialize some variables on first call
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF( LFIRST ) THEN

         DO NCSP = 1, NCS2
            DO ISP = 1, N_GC_SPC
               IZILCH( ISP, NCSP ) = 0
               JZILCH( ISP, NCSP ) = 0
               JHIZ1 ( ISP, NCSP ) = 0
               JHIZ2 ( ISP, NCSP ) = 0
               KZILCH( ISP, NCSP ) = 0
               MZILCH( ISP, NCSP ) = 0
            ENDDO
         ENDDO

         DO NCSP = 1, NCS2
            NUSERAT( NCSP ) = 0
            DO NK = 1, NRXNS
               NDERIVL( NK, NCSP ) = 0
               NDERIVP( NK, NCSP ) = 0
            ENDDO
         ENDDO
         
         DO NCSP = 1, NCS
            ISCHANG( NCSP ) = 0
         ENDDO
         
      ENDIF
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Initialize Prod/loss and PD tabulator arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NCSP = NCS

      NSPECT( NCS ) = N_GC_SPC
      DO ISP = 1, NSPECT( NCS )
         ISAPORL( ISP ) = 0
      ENDDO

      DO ISP = 1, NSPECT( NCS )
         DO ISP2 = 1, NSPECT( NCS )
            ISPARDER( ISP, ISP2 ) = 0
         ENDDO
      ENDDO
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set the number of Partial derivative terms in the Jacobian and
c  count the number of terms for each species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO NRX = 1, NRXNS
        DO NR = 1, 3
            IREACT = IRR( NRX, NR )
            IF( IREACT .NE. 0 ) THEN
               DO NRPP = 1, 3 + MXPRD
                  IPORR = IRR( NRX, NRPP )
                  IF( IPORR .NE. 0 ) ISPARDER( IPORR, IREACT ) = 1
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      DO IREACT = 1, NSPECT( NCS ) 
         DO IPORR = 1, NSPECT( NCS )
            IF( ISPARDER( IPORR, IREACT ) .EQ. 1 ) 
     &            ISAPORL( IPORR ) = ISAPORL( IPORR ) + 1
         ENDDO
      ENDDO
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Sort the species, putting all with zero partial derivative 
c  terms at the bottom and those with fewest PD terms at top.
c  Set arrays for species with zero PD terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NOCHANG = NSPECT( NCS )
      DO JOLD = 1, NSPECT( NCS )
         IF( ISAPORL( JOLD ) .GT. 0 ) THEN
            ISCHANG( NCS ) = ISCHANG( NCS ) + 1
            JNEW = ISCHANG( NCS )
            INEW2OLD( JNEW, NCS ) = JOLD
            IOLD2NEW( JOLD, NCS ) = JNEW
         ELSE
            INEW2OLD( NOCHANG, NCS ) = JOLD
            IOLD2NEW( JOLD, NCS ) = NOCHANG
            NOCHANG = NOCHANG - 1
         ENDIF
      ENDDO
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Now sort by number of PD terms, fewest at position 1, most at
c  the end position. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO JNEW = 1, ISCHANG( NCS )
c  uncomment the following three lines to turn off ordering
c         INEW2OLD( JNEW, NCS ) = JNEW
c         IOLD2NEW( JNEW, NCS ) = JNEW
c         IF( JNEW .NE. 0 ) GO TO 180
c  uncomment the above three lines to turn off ordering
         JOLD = INEW2OLD( JNEW, NCS )
         MINVALU = ISAPORL( JOLD )
         IMINOLD = JOLD
         IMINNEW = JNEW

         DO INEW = JNEW + 1, ISCHANG( NCS )
            IOLD = INEW2OLD( INEW, NCS )
            IF( ISAPORL( IOLD ) .LT. MINVALU ) THEN
               MINVALU = ISAPORL( IOLD )
               IMINOLD = IOLD
               IMINNEW = INEW
            ENDIF
         ENDDO

         INEW2OLD( IMINNEW, NCS ) = JOLD
         INEW2OLD( JNEW, NCS )    = IMINOLD
         IOLD2NEW( JOLD, NCS )    = IMINNEW
         IOLD2NEW( IMINOLD, NCS ) = JNEW
      ENDDO
               
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Fill the irm2 array using the new species order developed above.
c  Also determine active reactions for day and then night (i.e., photo
c  reactions determined by BTEST=.TRUE. are not included for nighttime)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO NRX = 1, NRXNS
         DO NR = 1, NREACT( NRX )
            IREACT = IRR( NRX, NR )
            IRM2( NRX, NR, NCS ) = IOLD2NEW( IREACT,NCS ) 
         ENDDO

         DO NP = 1, NPRDCT( NRX )
            IPROD = IRR( NRX, NP + 3 )
            IRM2( NRX, NP+3, NCS ) = IOLD2NEW( IPROD, NCS )
         ENDDO
         
         IF( NREACT( NRX ) .GT. 0 ) THEN
            NUSERAT( NCS ) = NUSERAT( NCS ) + 1
            NU = NUSERAT( NCS )
            NKUSERAT( NU, NCS ) = NRX
            IF( .NOT. ( BTEST ( IRXBITS( NRX ),1 ) ) ) THEN
               NUSERAT( NCS + 1 ) = NUSERAT( NCS + 1 ) + 1
               NU = NUSERAT( NCS + 1 )
               NKUSERAT( NU, NCS + 1 ) = NRX
            ENDIF
         ENDIF
      ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do symbolic LU decomposition to determine sparse storage array
c  structures. Done twice, first for day and then for night. An entry
c  of 1 in lzero means a non-negative entry in the Jacobian. First
c  put ones on the diagonal and zeroes everywhere else.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 700 IFSUN = 1, 2
         NCSP = IFSUN
         DO I = 1, N_GC_SPC
            DO J = 1, N_GC_SPC
               LZERO( J, I ) = 0
            ENDDO
            LZERO( I, I ) = 1
         ENDDO
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fill in the rest of the entries in the Jacobian
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO NRX = 1, NUSERAT( NCSP )
            NK = NKUSERAT( NRX, NCSP )
            DO NR = 1, NREACT( NK )
               IRE = IRM2( NK, NR, NCS )
               DO JAL = 1, NREACT( NK )
                  JRE = IRM2( NK, JAL, NCS )
                  LZERO( JRE, IRE ) = 1
               ENDDO
               DO IAP = 1, NPRDCT( NK )
                  JPR = IRM2( NK, 3+IAP, NCS )
                  LZERO( JPR, IRE ) = 1 
               ENDDO
           ENDDO
         ENDDO
   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set up arrays for decomposition / back-substitution of sparse     
c  matrices by removing all calculations involving a zero.          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF ( IFNEVER.EQ.0 ) THEN
            IFNEVER = 1
            ICNT    = 0 
            JCNT    = 0 
            ICCOUNT = 0
            JCCOUNT = 0
         ENDIF
         KOUNT0A = 0
         KOUNT0  = 0
         ICNTA   = 0
         ICNTB   = 0
         JCNTA   = 0
         JCNTB   = 0
         KCNTA   = 0
         MCNTA   = 0
         KCNT    = 0
         MCNT    = 0
         IARRAY2 = 0
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Count number of entries w/ and w/o sparse matrix storage
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
         DO J = 1, ISCHANG( NCS )
            DO K = 1, ISCHANG( NCS )
               KOUNT0A = KOUNT0A + 1
               IF( LZERO( J, K ) .EQ. 1 ) KOUNT0 = KOUNT0 + 1
            ENDDO
         ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do the symbolic decomposition (ludcmp) converting [A] to [L][U] 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ICLO( NCSP ) = ICNT + 1
         JCLO( NCSP ) = JCNT + 1
         DO J = 1, ISCHANG( NCS )
            J1 = J - 1
            
c...  First loop of decomposition
            DO I = 2, ISCHANG( NCS ) 
               I1 = J1 
               IF( I .LE. J1 ) I1 = I - 1
               DO K = 1, I1
                  ICNTA = ICNTA + 1
                  IF( LZERO( I, K ) .EQ. 1 .AND. LZERO( K, J ) .EQ. 1 )
     &                  THEN
                     IZILCH( J, NCSP ) = IZILCH( J, NCSP ) + 1
                     ICNT           = ICNT + 1
                     ICNTB          = ICNTB + 1
                     IZEROK( ICNT ) = K   
                     IZEROI( ICNT ) = I
                     LZERO( I, J )  = 1 
                  ENDIF
               ENDDO
            ENDDO
c... Second loop of decomposition 
            DO I = J + 1, ISCHANG( NCS ) 
               JCNTA = JCNTA + 1
               IF( LZERO( I, J ) .EQ. 1 ) THEN
                  JZILCH( J, NCSP ) = JZILCH( J, NCSP ) + 1
                  JCNT              = JCNT  + 1
                  JCNTB             = JCNTB + 1
                  JZERO( JCNT )     = I  
               ENDIF
            ENDDO 
         ENDDO
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do symbolic back-substition for solving [L][U]{x}={b}. Store data
c  in sparse matrix pointer jarraypt.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... First loop of back-substitution
         DO I = 2, ISCHANG( NCS )
            I1 = I - 1
            DO J = 1, I1    
               KCNTA = KCNTA + 1
               IF( LZERO( I, J ) .EQ. 1 ) THEN 
                  KZILCH( I, NCSP ) = KZILCH( I, NCSP ) + 1
                  KCNT = KCNT + 1
                  IARRAY2 = IARRAY2 + 1
                  KZERO( IARRAY2, NCSP ) = J
                  JARRAYPT( I, J, NCSP ) = IARRAY2 
               ENDIF
            ENDDO
         ENDDO 

c... Second loop of back-substitution 
         DO I = ISCHANG( NCS ) - 1, 1, -1
            I2 = I + 1
            DO J = I + 1, ISCHANG( NCS )
               MCNTA = MCNTA + 1
               IF( LZERO( I, J ) .EQ. 1 ) THEN 
                  MZILCH( I, NCSP )      = MZILCH( I, NCSP ) + 1
                  MCNT                   = MCNT + 1
                  IARRAY2                = IARRAY2 + 1
                  KZERO( IARRAY2, NCSP ) = J
                  JARRAYPT( I, J, NCSP ) = IARRAY2 
               ENDIF
            ENDDO
         ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fill jarraypt with remaining diagonal array points and save counts
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO I = 1, ISCHANG( NCS ) 
            IARRAY2 = IARRAY2 + 1
            JARRAYPT( I, I, NCSP ) = IARRAY2 
         ENDDO
         IARRAY( NCSP ) = IARRAY2 
         KNTARRAY = KCNTA + MCNTA + ISCHANG( NCS )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Do decomposition again to change arrays to use jarraypt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         JCB = JCLO( NCSP ) 
         JZLO( NCSP ) = JCCOUNT
         ICBSUM = ICLO( NCSP ) - 1 
         IJSTEP = 2   
         DO J = 1, ISCHANG( NCS )

c...First loop of decomposition
            IDEC1LO( J, NCSP ) = ICCOUNT + 1
            ICB = ICBSUM  + 1
            ICBSUM = ICBSUM + IZILCH( J, NCSP ) 

            DO KA = 1, IZILCH( J, NCSP ), IJSTEP
               ICCOUNT = ICCOUNT + 1
               IPA = IZEROI( ICB ) 
               KPA = IZEROK( ICB ) 
               IJDECA( ICCOUNT ) = JARRAYPT( IPA,   J, NCSP ) 
               IKDECA( ICCOUNT ) = JARRAYPT( IPA, KPA, NCSP )
               KJDECA( ICCOUNT ) = JARRAYPT( KPA,   J, NCSP )
               IF( ICB + 1 .LE. ICBSUM ) THEN
                  IPB = IZEROI( ICB + 1 ) 
                  KPB = IZEROK( ICB + 1 ) 
                  IJDECB( ICCOUNT ) = JARRAYPT( IPB,   J, NCSP ) 
                  IKDECB( ICCOUNT ) = JARRAYPT( IPB, KPB, NCSP )
                  KJDECB( ICCOUNT ) = JARRAYPT( KPB,   J, NCSP )
               ENDIF
               ICB = ICB + IJSTEP   
            ENDDO

            IDEC1HI( J, NCSP ) = ICCOUNT  
            
c...Second loop of decomposition
            JZ = JZILCH( J, NCSP )

            DO I = 1, JZ - 1, 2
               JCCOUNT           = JCCOUNT + 1
               JHIZ1( J, NCSP )  = JHIZ1( J, NCSP ) + 1
               IA                = JZERO( JCB )
               IB                = JZERO( JCB + 1 )
               JZEROA( JCCOUNT ) = JARRAYPT( IA, J, NCSP )
               JZEROB( JCCOUNT ) = JARRAYPT( IB, J, NCSP )
               JCB = JCB + 2
            ENDDO

            IF( MOD( JZ, 2 ) .EQ. 1 ) THEN 
               JCCOUNT           = JCCOUNT + 1
               JHIZ2( J, NCSP )  = JHIZ2( J, NCSP ) + 1
               IA                = JZERO( JCB )
               JZEROA( JCCOUNT ) = JARRAYPT( IA, J, NCSP )
               JCB               = JCB + 1 
            ENDIF
         ENDDO
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Group terms to increase efficiency in back-substition
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... First back-substitution loop
         DO I = 1, ISCHANG( NCS ) 
            KZ               = KZILCH( I, NCSP )
            KZHI0( I, NCSP ) = KZ - 4 
            JZ3              = 0

            DO JZ = 1, KZHI0( I, NCSP ), 5     
               JZ3 = JZ + 4
            ENDDO  

            KZLO1( I, NCSP )  = JZ3 + 1
            KZHI1( I, NCSP )  = KZ  - 1 
            JZ4 = JZ3 

            DO JZ = JZ3 + 1, KZ - 1, 2    
               JZ4 = JZ + 1
            ENDDO

            KZLO2( I, NCSP )   = JZ4 + 1
         ENDDO
 
c... Second loop of back-substitution
         DO I = ISCHANG( NCS ), 1, -1
            MZ = MZILCH( I, NCSP ) 
            MZHI0( I, NCSP ) = MZ - 4  
            JZ3 = 0 

            DO JZ = 1, MZHI0( I, NCSP ), 5  
               JZ3 = JZ + 4 
            ENDDO

            MZLO1( I, NCSP ) = JZ3 + 1
            MZHI1( I, NCSP ) = MZ  - 1
            JZ4 = JZ3 

            DO JZ = JZ3+1, MZ-1, 2 
               JZ4 = JZ + 1 
            ENDDO

            MZLO2( I, NCSP )   = JZ4 + 1
         ENDDO
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check dimensions and print out array savings if ldebug on 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF( ICNT.GT.MXCOUNT2 .OR. JCNT .GT. MXCOUNT1 .OR. 
     &         IARRAY2 .GT. MXARRAY .OR. ICCOUNT .GT. MXCOUNT2 .OR.
     &         JCCOUNT .GT. MXARRAY ) THEN
            WRITE( MSG, 94000 ) 
            CALL M3MESG( MSG )
            WRITE( MSG, 94020 ) MXCOUNT2, ICNT 
            CALL M3MESG( MSG )
            WRITE( MSG, 94040 ) MXCOUNT1, JCNT 
            CALL M3MESG( MSG )
            WRITE( MSG, 94060 ) MXARRAY, IARRAY2 
            CALL M3MESG( MSG )
            WRITE( MSG, 94080 ) MXARRAY, ICCOUNT 
            CALL M3MESG( MSG )
            WRITE( MSG, 94100 ) MXARRAY, JCCOUNT 
            CALL M3MESG( MSG )
            EXITSTAT = 2
            CALL M3EXIT( PNAME, IZERO, IZERO, ' ', EXITSTAT )
         ENDIF           
         IF( LDEBUG ) THEN
            WRITE( IUNDBG, 93000 ) NCSP, KOUNT0A, KOUNT0, KNTARRAY, 
     &                             IARRAY2, ICNTA, ICNTB, JCNTA, JCNTB,
     &                              KCNTA, KCNT, MCNTA, MCNT
         ENDIF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set final arrays for partial derivative calculations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO NRX = 1, NUSERAT( NCSP )
            NK = NKUSERAT( NRX, NCSP )
            DO IAL = 1, NREACT( NK )
               IR = IRM2( NK, IAL, NCS )

               DO JAL = 1, NREACT( NK )
                  JR = IRM2( NK, JAL, NCS )
                  IAR = JARRAYPT( JR, IR, NCSP )
                  NDERIVL( NK, NCSP ) = NDERIVL( NK, NCSP ) + 1
                  NLS = NDERIVL( NK, NCSP )
                  JARRL( NK, NLS, NCSP ) = IAR
                  JLIAL( NK, NLS, NCSP ) = IAL
                  NDLMAX = MAX( NLS, NDLMAX )
               ENDDO
               
               DO IAP = 1, NPRDCT( NK )
                  JP = IRM2( NK, IAP + 3, NCS )
                  IAR = JARRAYPT( JP, IR, NCSP )
                  NDERIVP( NK, NCSP ) = NDERIVP( NK, NCSP ) + 1
                  NPR = NDERIVP( NK, NCSP )
                  JARRP(  NK, NPR, NCSP ) = IAR
                  JPIAL(  NK, NPR, NCSP ) = IAL
                  ICOEFF( NK, NPR, NCSP ) = 0
                  IF( ABS( SC( NK, IAP ) - 1.0 ) .GT. 1.0E-06 ) THEN
                     ICOEFF( NK, NPR, NCSP ) = IAP
                  ENDIF
                  NDPMAX = MAX( NPR, NDPMAX )
               ENDDO
            ENDDO     
         ENDDO
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check dimensions of PD arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF( NDPMAX .GT. MXRP .OR. NDLMAX .GT. MXRR ) THEN
            WRITE( MSG, 94000 ) 
            CALL M3MESG( MSG )
            WRITE( MSG, 94200 ) MXRP, NDPMAX 
            CALL M3MESG( MSG )
            WRITE( MSG, 94220 ) MXRR, NDLMAX 
            CALL M3MESG( MSG )
            EXITSTAT = 2
            CALL M3EXIT( PNAME, IZERO, IZERO, ' ', EXITSTAT ) 
         ENDIF
  700 CONTINUE
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set some parameters for the Gear integration method 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF( LFIRST ) THEN 
         LFIRST   = .FALSE. 
         MSTEP    = 3
         MAXORD   = 5
         MBETWEEN = 50

         DO NQQ = 1, 7
            ENQQ1(  NQQ )    = 0.5E+00 / FLOAT( NQQ     )
            ENQQ2(  NQQ )    = 0.5E+00 / FLOAT( NQQ + 1 ) 
            ENQQ3(  NQQ )    = 0.5E+00 / FLOAT( NQQ + 2 )
            CONPST( NQQ )    = 1.0E+00 / ( PERTST( NQQ, 1 ) *
     &                         ENQQ3( NQQ ) ) 
            CONP15( NQQ )    = 1.5E+00 * CONPST( NQQ )
            PERTST( NQQ, 1 ) = PERTST( NQQ, 1 ) * PERTST( NQQ, 1 )
            PERTST( NQQ, 2 ) = PERTST( NQQ, 2 ) * PERTST( NQQ, 2 )
            PERTST( NQQ, 3 ) = PERTST( NQQ, 3 ) * PERTST( NQQ, 3 )
         ENDDO

         DO I2 = 1, 6 
            ASET( I2, 2 ) = 1.0E+00
            ASET( I2, 8 ) = 0.
         ENDDO
         ASET( 1, 1 )   = 1.0                                                            
         ASET( 2, 1 )   = 2.0        /    3.0      
         ASET( 2, 3 )   = 1.0        /    3.0      
         ASET( 3, 1 )   = 6.0        /   11.0    
         ASET( 3, 3 )   = 6.0        /   11.0       
         ASET( 3, 4 )   = 1.0        /   11.0                                                
         ASET( 4, 1 )   = 12.0       /   25.0                                                
         ASET( 4, 3 )   =   .70        
         ASET( 4, 4 )   =   .20       
         ASET( 4, 5 )   =   .020      
         ASET( 5, 1 )   =   60.0     /  137.0    
         ASET( 5, 3 )   =  225.0     /  274.0      
         ASET( 5, 4 )   =   85.0     /  274.0      
         ASET( 5, 5 )   =   15.0     /  274.0      
         ASET( 5, 6 )   =    1.0     /  274.0                                                
         ASET( 6, 1 )   =  180.0     /  441.0     
         ASET( 6, 3 )   =  406.0     /  441.0       
         ASET( 6, 4 )   =  735.0     / 1764.0      
         ASET( 6, 5 )   =  175.0     / 1764.0    
         ASET( 6, 6 )   =   21.0     / 1764.0      
         ASET( 6, 7 )   =    1.0     / 1764.0      
      ENDIF
      RETURN
      
C********************** FORMAT STATEMENTS ******************************      
93000 FORMAT( 1X,/'PARAM    POSS MATRIX POINTS -- NONZEROS -- NCSP=',I4/
     &        1X, 'INITMAT  ',4X,I8,9X,I8/                               
     &        1X, 'FINMAT   ',4X,I8,9X,I8/   
     &        1X, 'DECOMP1  ',4X,I8,9X,I8/ 
     &        1X, 'DECOMP2  ',4X,I8,9X,I8/ 
     &        1X, 'BACKSB1  ',4X,I8,9X,I8/ 
     &        1X, 'BACKSB2  ',4X,I8,9X,I8/)  
94000 FORMAT( 1X,'One of the dimensions below is too small:')
94020 FORMAT( 1X,'DIMENSION: MXCOUNT2 = ',I6,' VARIABLE: ICNT    = ',I6)  
94040 FORMAT( 1X,'DIMENSION: MXCOUNT1 = ',I6,' VARIABLE: JCNT    = ',I6)  
94060 FORMAT( 1X,'DIMENSION: MXARRAY  = ',I6,' VARIABLE: IARRAY2 = ',I6)  
94080 FORMAT( 1X,'DIMENSION: MXARRAY  = ',I6,' VARIABLE: ICCOUNT = ',I6)  
94100 FORMAT( 1X,'DIMENSION: MXARRAY  = ',I6,' VARIABLE: JCCOUNT = ',I6)
94200 FORMAT( 1X,'DIMENSION: MXRP     = ',I6,' VARIABLE: NDPMAX  = ',I6)
94220 FORMAT( 1X,'DIMENSION: MXRR     = ',I6,' VARIABLE: NDLMAX  = ',I6)
      END
