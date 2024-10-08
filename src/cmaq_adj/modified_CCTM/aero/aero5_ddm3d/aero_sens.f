
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_ddm3d/aero_sens.f,v 1.1.1.1 2010/06/14 16:02:59 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE AERO_SENS( NSPCSDA, CBLK, GAS, AERLIQ, AERSLD, SBLK, 
     &                      NPMAX )
c
c
c  Jim Boylan, Yeuh-Jin Yang, Ted Russell, and others at Georgia Tech - Aug 99
c     -- derived the aerosol equilibrium sensitivity equations for ISORROPIA
c     -- implemented and tested for URM-1ATM model
c
c  Sergey L. Napelenok - Apr 03
c     -- adapted for the CMAQ ver4.3
c
c  Sergey L. Napelenok - Jul 06
c     -- updated for CMAQ ver4.5
c  
c  Sergey L. Napelenok - May 08
c     -- updated for CMAQ ver4.7
c     -- changed the call structure (now called from AEROPROC)
c     -- now accepts the SBLK array assigned in AERO
c     -- now uses AERO_INFO to get index assignments
c     -- cleared a bunch of unused variables
c
c
c  Variables GAS, AERLIQ, AERSLD are passed from ISORROPIA
c
c  GAS		GAS concentration vector (mole/m3-air)
c		GAS(1)		NH3(g)
c		GAS(2)		HNO3(g)
c		GAS(3)		HCl(g)
c
c  AERLIQ  	aqueous concentration vector (mole/m3-air)
c		AERLIQ(1)	H+(aq)
c		AERLIQ(2)	Na+(aq)
c		AERLIQ(3)	NH4+(aq)
c		AERLIQ(4)	Cl-(aq)
c		AERLIQ(5)	SO4=(aq)
c		AERLIQ(6)	HSO4-(aq)
c		AERLIQ(7)	NO3-(aq)
c		AERLIQ(8)	H2O(aq)
c		AERLIQ(9)	NH3(aq) undissociated
c		AERLIQ(10)	HCl(aq) undissociated
c		AERLIQ(11)	HNO3(aq) undissociated
c		AERLIQ(12) 	OH-(aq)
c
c  AERSLD  	solid concentration vector (mole/m3-air)
c		AERSLD(1)	NaNO3(s)
c		AERSLD(2)	NH4NO3(s)
c		AERSLD(3)	NaCl(s)
c		AERSLD(4)	NH4Cl(s)
c		AERSLD(5)	Na2SO4(s)
c		AERSLD(6)	(NH4)2SO4(s)
c		AERSLD(7)	NaHSo4(s)
c		AERSLD(8)	NH4HSO4(s)
c		AERSLD(9)	(NH4)4H(SO4)2(s)
c
c  asen		speciated aerosol sensitivities (total of I & J modes)
c  		asen(1)	aerosol sulfate
c  		asen(2)	aerosol nitrate
c  		asen(3)	aerosol ammonia
c  		asen(4)	aerosol sodium
c  		asen(5)	aerosol chlorine
c  		asen(6)	aerosol hydrogen
c
c  nsize	size of a_coef
c  NPMAX	size of sensitivity parameter
c  naq		number of aqueous aerosol species
c  nsolid	number of solid aerosol species
c  dt		time step (min)
c
      USE AERO_INFO           ! aero include files

      IMPLICIT NONE

C sln =-=-=-=-=-=-=  DDM-3D sensitivity variables

      INTEGER NSPCSDA       ! number of species in CBLK
      REAL, DIMENSION( NSPCSDA )        :: CBLK     ! main array of concentrations 
      REAL, DIMENSION( NSPCSDA, NPMAX ) :: SBLK     ! sensitivity equivalent of CBLK

      REAL(KIND=8), DIMENSION(3)        :: GAS      ! gas-phase ISOROPIA output
      REAL(KIND=8), DIMENSION(12)       :: AERLIQ   ! liq aerosol ISOROPIA output
      REAL(KIND=8), DIMENSION(9)        :: AERSLD   ! solid ISOROPIA output

      INTEGER NPMAX                                 ! defined number of sensivity parameters
      INTEGER NP                                    ! counter for NPMAX

      REAL, PARAMETER                   :: cmin = 1.0E-30     ! minimum concentration

      REAL(KIND=8), PARAMETER           :: f0 = 1.005d0
      REAL(KIND=8)                      :: f  = 1.000d0

      INTEGER, PARAMETER                      :: nsize = 26             ! full matrix size
      REAL(KIND=8), DIMENSION(nsize)          :: acon
      REAL(KIND=8), DIMENSION(nsize,nsize)    :: a_coef
      REAL(KIND=8), DIMENSION(nsize*nsize)    :: atemp
      REAL(KIND=8), DIMENSION(nsize)          :: rhs

      INTEGER, PARAMETER                   :: gas1 = 1,  gas2 = 3      
      INTEGER, PARAMETER                   :: aq1  = 4,   aq2 = 10     
      INTEGER, PARAMETER                   :: sld1 = 11, sld2 = 20
      INTEGER, DIMENSION(nsize)            :: row_flag, col_flag, ipvt
      INTEGER nzero, nrow, ncol, nrowxncol

      INTEGER i, j, k
      INTEGER info, ind

      REAL(KIND=8), DIMENSION(6)        :: asen

      REAL fji                                      ! inorganic Aitken  mode fraction

C sln =-=-=-=-=-=-=  end DDM-3D sensitivity variables

c
c calculate modal concentration based modal fractions
c

      fji = ( CBLK( VSO4AI ) + CBLK( VNH4AI ) + CBLK( VNO3AI ) ) 
     &   / MAX ( ( CBLK( VSO4AJ ) + CBLK( VNH4AJ ) + CBLK( VNO3AJ )
     &   +   CBLK( VSO4AI ) + CBLK( VNH4AI ) + CBLK( VNO3AI ) ), 
     &   cmin )
      fji = MIN( MAX( 0.0, fji ), 1.0 )  

  90  continue

c
c  get coeff in the matrix
c
      call coeff(gas, aerliq, aersld, acon, a_coef, f,nsize)

c
c  initialize row_flag and col_flagc
c
      row_flag = 1
      col_flag = 1
      ipvt = 0

 100  continue
c
c  initialize the number of zero-sensitivities
c
      nzero = 0
c
c  check if the matrix size needs adjustment
c
      do i = 1, 20
         if(acon(i) .eq. 0.0d0) goto 180
      end do
      goto 500
c
c (1) gas phase, 1st ~ 3rd species
c     if acon(j)=0, both corresponding column and row
c     can be eliminated from the a_coef matrix
c
  180 continue

      do 200 j = gas1, gas2
         if(acon(j) .eq. 0.0) then
            col_flag(j) = 0
            row_flag(j) = 0
            do 205 i = 1, nsize
               a_coef(i, j) = 0.0d0
               a_coef(j, i) = 0.0d0
  205       continue
        endif

  200 continue
c
c (2) solid phase, 11th ~ 20th species
c     if acon(j)=0, both corresponding column and row
c     can be eliminated from the a_coef matrix
c
      do 210 j = sld1, sld2
         if(acon(j) .eq. 0.0d0) then
            col_flag(j) = 0
            row_flag(j) = 0
            do 215 i = 1, nsize
               a_coef(i, j) = 0.0d0
               a_coef(j, i) = 0.0d0
  215       continue
         endif
  210 continue
c
c (3) aqueous phase, 4th ~ 10th species
c     if acon(j)=0, only the corresponding column
c     can be eliminated from the a_coef matrix
c
      do 220 j = aq1, aq2
         if(acon(j) .eq. 0.0d0) then
            col_flag(j) = 0
            do 225 i = 1, nsize
               a_coef(i, j) = 0.0d0
  225       continue
         endif
  220 continue
c
c sln - special case with no Na+, CL-
      if ( col_flag(2) .eq. 0 .and. col_flag(4) .eq. 0 ) then
        row_flag(7) = 0
        row_flag(8) = 0 
        a_coef(7, 24) = 0.0d0
        a_coef(8, 25) = 0.0d0
      endif 
c
c  all columns have been marked for elimination,
c  check rows with all zero elements
c
      do 240 i = 1, nsize
            row_flag(i) = 0
            do 250 j = 1, nsize
               if (a_coef(i, j) .ne. 0.0d0) row_flag(i) = 1
  250       continue
  240 continue
c
c  check the total number of nonzero rows and columns,
c  if not equal, then stop
c
  500 continue

      nrow = 0
      ncol = 0

      do 260 i = 1, nsize
         if(row_flag(i) .eq. 1) nrow = nrow + 1
         if(col_flag(i) .eq. 1) ncol = ncol + 1
  260 continue

       if(nrow .ne. ncol) goto 9900

c  start matrix size reduction
c  "1" and "0" for row and column reduction respectively
c

      call coeffresize(a_coef, row_flag, col_flag, nrow, ncol,
     &                 nsize)


c
c  note: nrow = ncol, store the size-reduced matrix as a temp array
c
            k = 0
            do j = 1, ncol
               do i = 1, nrow
                  k = k + 1
                  atemp(k) = a_coef(i, j)
               end do
            end do
c
        nrowxncol = nrow*ncol
c
c  solve {A}[x] = [R] using the Harwell LU-decomposition
c
        info = 0
        call dgefa(atemp, nrow, ncol, ipvt, info)
        if(info .ne. 0) then
           if(f .gt. 1.22d0) then
              write(14,*)'INFO =',info 
              ind = info
              do i=1,nrow
                write(14,1313) (a_coef(i, j),j=1,ncol)
              end do
1313          format (26(e12.4,1x))
              write(14,*),"GAS"
              do i = 1, 3
                write(14,*),GAS(i)
              end do
              write(14,*),"AERLIQ"
              do i = 1, 12
                write(14,*),AERLIQ(i)
              end do
              write(14,*),"AERSLD"
              do i = 1, 9
                write(14,*),AERSLD(i)
              end do
              stop
           endif
           f = f0*f
           goto 90
        end if

7777    FORMAT(1x,26(1x,E13.5))

 
c  loop over each sensitivity
 
        do 1100 NP = 1, NPMAX

c  get all rhs terms, vary for each p(i)
 
            call get_rhs(rhs, NSPCSDA, SBLK(:, NP), f, nsize )

            call rhs_resize(row_flag, rhs, nrow, nsize)

            call dgesl(atemp, nrow, ncol, ipvt, rhs, 0)
 
c  put 3 gas sensitivities back to the SBLK() vector, update asen(6) vector

            call asupdt(NSPCSDA, SBLK(1,NP), row_flag, rhs, nrow, 
     &                  asen, f, nsize)

c --- convert to ug/m^3 

           asen(1) =asen(1) * dble(MWSO4)
           asen(2) =asen(2) * dble(MWNO3)
           asen(3) =asen(3) * dble(MWNH4)
           asen(4) =asen(4) * dble(MWNA)
           asen(5) =asen(5) * dble(MWCL)
c          asen(6) =asen(6) * 1.d0

1099  CONTINUE
c --- Distribute mass between the 2 sizes  !!!! ip = NP

           SBLK( VSO4AI,NP ) = sngl( asen( 1 ) ) * fji
           SBLK( VSO4AJ,NP ) = sngl( asen( 1 ) ) - SBLK( VSO4AI,NP )

           SBLK( VNO3AI,NP ) = sngl( asen( 2 ) ) * fji
           SBLK( VNO3AJ,NP ) = sngl( asen( 2 ) ) - SBLK( VNO3AI,NP )

           SBLK( VNH4AI,NP ) = sngl( asen( 3 ) ) * fji
           SBLK( VNH4AJ,NP ) = sngl( asen( 3 ) ) - SBLK( VNH4AI,NP )

c          SBLK( VNAJ,NP   ) = sngl( asen( 4 ) ) * fji 
c          SBLK( VNAI,NP   ) = sngl( asen( 4 ) ) - SBLK( VNAJ,NP   )

c          SBLK( VCLJ,NP   ) = sngl( asen( 5 ) ) * fji 
c          SBLK( VCLI,NP   ) = sngl( asen( 5 ) ) - SBLK( VCLJ,NP   )

 1100  continue                                 !close of NP
 
      goto 9999
 9900 continue

      print*,'the no of nonzero rows and cols are inconsistent'
      print*,'in aerosol sensitivity module'
      print*,'sensitivity parameter (NP) = ', NP
      print*,'nrow= ', nrow,'ncol= ', ncol
      do i = 1, nsize
         print*,'i=',i,'row_flag=',row_flag(i),'col_flag=',col_flag(i)
      end do

             write(14,*),"GAS"
             do i = 1, 3
               write(14,*),GAS(i)
             end do

             write(14,*),"AERLIQ"
             do i = 1, 12
               write(14,*),AERLIQ(i)
             end do

             write(14,*),"AERSLD"
             do i = 1, 9
               write(14,*),AERSLD(i)
             end do
          
             do i = 1, nsize
               write(14,9901) (a_coef(i,j),j=1,nsize)
             end do
9901         format(26(1x,e12.4)) 

      stop
 9999 continue
      return
c
      end
c
c
c .................................................................
c
      subroutine coeff(gas, aerliq, aersld, acon, a_coef, f, nsize)

      USE AERO_INFO           ! aero include files

      implicit none

      integer nsize
      integer i, j
      real*8 acon(nsize), a_coef(nsize, nsize)
      real*8 HION, f
      real*8 gas(3)
      real*8 aerliq(12)
      real*8 aersld(9)
      real, parameter :: cmin=1.0e-30

c
c  acon         concentration array
c  a_coef	constant matrix coefficients
c  nsize	size of constant matrix
c
c
c  initialize coefficient matrix
c
      a_coef = 0.0d0

  100 continue
c
      do 105 i = 1, nsize
         acon(i) = 0.0d0
  105 continue
c
c  get gas concentrations (nanomoles/m3-air)
c
c*** H+ and OH-  ***
c limit pH to [0.1,13]
      HION = aerliq(1)/(aerliq(8)*dble(MWWAT)/1000.d+00)
      HION = max( 1.d-13, min( 0.79433d+00, HION ) )
      acon( 6) = HION           * (aerliq(8)*dble(MWWAT)/1000.d+00 ) * 1.d+9
      acon( 7) = (1.0d-14/HION) * (aerliq(8)*dble(MWWAT)/1000.d+00 ) * 1.d+9

c*** HNO3 ***
c
      acon( 1) = DBLE(gas(2)*1.0e+09)
c
c*** HCL  ***
c
      acon( 2) = DBLE(gas(3)*1.0e+09)
c
c*** NH3  ***
c
      acon( 3) = DBLE(gas(1)*1.0e+09)
c
c*** NA+  ***
c
      acon( 4) = DBLE(aerliq(2)*1.0e+09)
c
c*** NH4+ ***
c
      acon( 5) = DBLE(aerliq(3)*1.0e+09)
c
c*** SO4= ***
c
      acon( 8) = DBLE(aerliq(5)*1.0e+09)
c
c*** NO3- ***
c
      acon( 9) = DBLE(aerliq(7)*1.0e+09)
c
c*** CL-  ***
c
      acon(10) = DBLE(aerliq(4)*1.0e+09)
c
c*** HSO4- ***
c
      acon(11) = DBLE(aerliq(6)*1.0e+09)
c
c   get solid phase concentrations (nanomoles/M3-air)
c
c*** NaHSO4 ***
c
      acon(12) = DBLE(aersld(7)*1.0e+09)
c
c*** NH4HSO4 ***
c
c
      acon(13) = DBLE(aersld(8)*1.0e+09)
c
c*** NASO4 ***
c
      acon(14) = DBLE(aersld(5)*1.0e+09)
c
c*** (NH4)2SO4 ***
c
      acon(15) = DBLE(aersld(6)*1.0e+09)
c
c*** (NH4)3H(SO4)2 ***
c
      acon(16) = DBLE(aersld(9)*1.0e+09)
c
c*** NH4NO3 ***
c
      acon(17) = DBLE(aersld(2)*1.0e+09)
c
c*** NH4CL ***
c
      acon(18) = DBLE(aersld(4)*1.0e+09)
c
c*** NANO3 ***
c
      acon(19) = DBLE(aersld(1)*1.0e+09)
c
c*** NACL ***
c
      acon(20) = DBLE(aersld(3)*1.0e+09)

c
c  set up coefficients
c
      a_coef( 1,  1) = -acon( 9)*acon( 6)*f*f
      a_coef( 1,  6) =  acon( 1)*acon( 9)*f*f
      a_coef( 1,  9) =  acon( 1)*acon( 6)*f*f
      a_coef( 2,  2) = -acon( 6)*acon(10)*f*f
      a_coef( 2,  6) =  acon( 2)*acon(10)*f*f
      a_coef( 2, 10) =  acon( 2)*acon( 6)*f*f
      a_coef( 3,  3) = -acon( 5)*acon( 7)*f*f
      a_coef( 3,  5) =  acon( 3)*acon( 7)*f*f
      a_coef( 3,  7) =  acon( 3)*acon( 5)*f*f
      a_coef( 4, 21) =  1.0d0
      a_coef( 5,  1) =  1.0d0
      a_coef( 5, 22) =  1.0d0
      a_coef( 6,  3) =  1.0d0
      a_coef( 6, 23) =  1.0d0
      a_coef( 7, 24) =  1.0d0
      a_coef( 8,  2) =  1.0d0
      a_coef( 8, 25) =  1.0d0
      a_coef( 9,  4) =  1.0d0
      a_coef( 9,  5) =  1.0d0
      a_coef( 9,  6) =  1.0d0
      a_coef( 9,  7) = -1.0d0
      a_coef( 9,  8) = -2.0d0
      a_coef( 9,  9) = -1.0d0
      a_coef( 9, 10) = -1.0d0
      a_coef( 9, 11) = -1.0d0
      a_coef(10,  6) =  acon( 7)*f
      a_coef(10,  7) =  acon( 6)*f
      a_coef(11,  6) =  acon( 8)*acon(11)*f*f
      a_coef(11,  8) =  acon( 6)*acon(11)*f*f
      a_coef(11, 11) = -acon( 6)*acon( 8)*f*f
      a_coef(12,  4) =  acon(11)*f
      a_coef(12, 11) =  acon( 4)*f
      a_coef(13,  5) =  acon(11)*f
      a_coef(13, 11) =  acon( 5)*f
      a_coef(14,  4) =  2.0d0*acon( 8)*f
      a_coef(14,  8) =  acon( 4)*f
      a_coef(15,  5) =  2.0d0*acon( 8)*f
      a_coef(15,  8) =  acon( 5)*f
      a_coef(16,  5) =  3.0d0*acon( 8)*acon(11)*f*f
      a_coef(16,  8) =  acon( 5)*acon(11)*f*f
      a_coef(16, 11) =  acon( 5)*acon( 8)*f*f
      a_coef(17,  1) =  acon( 3)*f
      a_coef(17,  3) =  acon( 1)*f
      a_coef(18,  2) =  acon( 3)*f
      a_coef(18,  3) =  acon( 2)*f
      a_coef(19,  4) =  acon( 9)*f
      a_coef(19,  9) =  acon( 4)*f
      a_coef(20,  4) =  acon(10)*f
      a_coef(20, 10) =  acon( 4)*f
      a_coef(21,  8) =  1.0d0
      a_coef(21, 11) =  1.0d0
      a_coef(21, 12) =  1.0d0
      a_coef(21, 13) =  1.0d0
      a_coef(21, 14) =  1.0d0
      a_coef(21, 15) =  1.0d0
      a_coef(21, 16) =  2.0d0
      a_coef(21, 21) = -1.0d0
      a_coef(22,  9) =  1.0d0
      a_coef(22, 17) =  1.0d0
      a_coef(22, 19) =  1.0d0
      a_coef(22, 22) = -1.0d0
      a_coef(23,  5) =  1.0d0
      a_coef(23, 13) =  1.0d0
      a_coef(23, 15) =  2.0d0
      a_coef(23, 16) =  3.0d0
      a_coef(23, 17) =  1.0d0
      a_coef(23, 18) =  1.0d0
      a_coef(23, 23) = -1.0d0
      a_coef(24,  4) =  1.0d0
      a_coef(24, 12) =  1.0d0
      a_coef(24, 14) =  2.0d0
      a_coef(24, 19) =  1.0d0
      a_coef(24, 20) =  1.0d0
      a_coef(24, 24) = -1.0d0
      a_coef(25, 10) =  1.0d0
      a_coef(25, 18) =  1.0d0
      a_coef(25, 20) =  1.0d0
      a_coef(25, 25) = -1.0d0
      a_coef(26, 21) = -2.0d0
      a_coef(26, 22) = -1.0d0
      a_coef(26, 23) =  1.0d0
      a_coef(26, 24) =  1.0d0
      a_coef(26, 25) = -1.0d0
      a_coef(26, 26) =  1.0d0
 
      return
      end

c .................................................................

      subroutine get_rhs(rhs, NSPCSDA, snt,f,nsize ) 

      USE AERO_INFO

      implicit none

      INTEGER nsize                                      ! full matrix size
      REAL(KIND=8), DIMENSION(nsize)     :: rhs          ! "right-hand" side
      INTEGER NSPCSDA       ! number of species in CBLK
      REAL, DIMENSION( NSPCSDA )         :: snt          ! SBLK for 1 parameter (ug/m3)
      REAL(KIND=8)                       :: f

      REAL(KIND=8)                       :: rhs4, rhs5, rhs6, rhs7, rhs8
      
c  initialize rhs terms
 
      rhs = 0.0d+00
 
c  note: d(total sulfate)/dp is same before and after aerosol module,
c  so are other aerosol species
 
c  total sulfate: AS  (sulfuric acid added here)
 
      rhs4 = DBLE( snt( VSULF ) * 1.0e+3 / MWH2SO4                      ! H2SO4
     &           + ( snt( VSO4AJ ) + snt( VSO4AI ) ) * 1.0e+3 / MWSO4 ) ! ASO4J + ASO4I
 
c  total nitrate: HNO3 + AN
 
      rhs5 = DBLE( snt( VHNO3 ) * 1.0e+3 / MWHNO3                       ! HNO3
     &           + ( snt( VNO3AJ ) + snt( VNO3AI ) ) * 1.0e+3 / MWNO3 ) ! ANO3J + ANO3I

 
c  total amomnia: NH3 + AA
 
      rhs6 = DBLE( snt( VNH3 ) * 1.0e+3 / MWNH3                         ! NH3
     &           + ( snt( VNH4AJ ) + snt( VNH4AI ) ) * 1.0e+3 / MWNH4  )! ANO3J + ANH4I

 
c  total sodium: ANa

c     rhs7 = DBLE( ( snt( VNAJ ) + snt( VNAI ) ) * 1.0e+3 / MWNA )      ! ANAJ + ANAI
 
c  total cloride: HCl + ACl
 
c     rhs8 = DBLE( snt( VHCL ) * 1.0e+3 / MWHCL                         ! HCL
c    &           + ( snt( VCLJ ) + snt( VCLI ) ) * 1.0e+3 / MWCL  )     ! ACLJ + ACLI

c populate rhs

      rhs(4) = rhs4 * f
      rhs(5) = rhs5 * f
      rhs(6) = rhs6 * f
      rhs(7) = rhs7 * f
      rhs(8) = rhs8 * f
 
      return
      end

c .................................................................

      subroutine coeffresize(a_coef, row_flag, col_flag, nrow, ncol,
     &                       nsize)

      IMPLICIT NONE

c
      integer nsize
      integer row_flag(nsize), col_flag(nsize)
      integer i, j, irow, icol, nrow, ncol
      double precision a_coef(nsize, nsize), temp_coef(nsize, nsize)
c
c
c
      i = 0
      do 100 irow = 1, nsize
         j = 0
         if(row_flag(irow) .eq. 0) goto 100
            i = i + 1
            if(i .gt. nrow) write(6,900) i, nrow
            do 200 icol = 1, nsize
               if(col_flag(icol) .eq. 0) goto 200
               j = j + 1
               if(j .gt. ncol) write(6,901) j, ncol
               temp_coef(i, j) = a_coef(irow, icol)
 200        continue
 100  continue
c
c put temp_coef back into a_coef
c
      do 300 i = 1, nrow
      do 300 j = 1, ncol
         a_coef(i, j) = temp_coef(i, j)
 300  continue
 400  continue
c
  900 format('FAILURE OF COEFF RESIZE',i3,'.GT. NROW=',i3)
  901 format('FAILURE OF COEFF RESIZE',i3,'.GT. NCOL=',i3)
c
      return
      end

c .................................................................

      subroutine rhs_resize(row_flag, rhs, nrow,nsize)

      IMPLICIT NONE

c  rworcl = 1 for row, and 0 for column
c

      integer nsize
      integer row_flag(nsize), temp_flag(nsize)
      integer i, is, icount, nrow, nzero, nsize1
      double precision rhs(nsize)
c      real rhs(*)
c
c  no of rows with all zeros
c
      nzero = nsize - nrow
      nsize1 = nsize - 1
c
c  copy row_flag as a temp array
c
      do 100 i = 1, nsize
         temp_flag(i) = row_flag(i)
  100 continue
c
      i = 0
      icount = 0

  300 continue
c
      i = i + 1
c
      if(icount .ge. nzero) go to 500

      if(temp_flag(i) .eq. 0) then
c
c  adjust rhs vector based on temp_flag()
c
         do is = i, nsize1
            temp_flag(is) = temp_flag(is + 1)
            rhs(is) = rhs(is + 1)
         end do

         icount = icount + 1

         i = i - 1

      endif

      go to 300

  500 continue
c
      return
      end

c .................................................................

      subroutine asupdt(NSPCSDA, snt, row_flag, rhs, nrow,
     &                  asen, f, nsize)

      USE AERO_INFO

      implicit none

      INTEGER nsize                                      ! full matrix size
      REAL(KIND=8), DIMENSION(nsize)     :: rhs, rhs2    ! "right-hand" side
      INTEGER NSPCSDA       ! number of species in CBLK
      REAL, DIMENSION( NSPCSDA )         :: snt          ! SBLK for 1 parameter (ug/m3)
      REAL(KIND=8)                       :: f

      INTEGER, DIMENSION(nsize)          :: row_flag
      REAL(KIND=8), DIMENSION(6)         :: asen

      INTEGER nrow
      INTEGER is, i, k, n

      REAL, PARAMETER                   :: cmin = 1.0E-30     ! minimum concentration

      is = 0
      n = nrow

      asen = 0.0d0

c convert from nanomol/m3 to umol/m3
      do k = 1, nrow
         rhs2(k) = rhs(k)*1.d-03
      end do
c
      do 400 i = 1, nsize
         if(row_flag(i) .eq. 1) then
              is = is + 1
              if(is .gt. n) goto 9999

              if(i .eq. 1) then
                 snt( VHNO3 ) = sngl( rhs2(is) * MWHNO3 / f )  ! (HNO3)
              elseif(i .eq. 2) then
c                snt( VHCL  ) = sngl( rhs2(is) * MWHCL  / f )  ! (HCL)
              elseif(i .eq. 3) then
                 snt( VNH3 )  = sngl( rhs2(is) * MWNH3  / f )  ! (NH3)
              elseif(i .eq. 21) then
                 asen(1) = rhs2(is)/f
              elseif(i .eq. 22) then
                 asen(2) = rhs2(is)/f
              elseif(i .eq. 23) then
                 asen(3) = rhs2(is)/f
              elseif(i .eq. 24) then
                 asen(4) = rhs2(is)/f
              elseif(i .eq. 25) then
                 asen(5) = rhs2(is)/f
              elseif(i .eq. 26) then
                 asen(6) = rhs2(is)/f
              endif
         endif
  400 continue
      return
c
 9999 continue
      print*,'*** the row_flag() and rhs() are inconsistant ***'
      print*,'*** stop at function asupdt() ***'
      stop
      end

c .................................................................

      subroutine trid(a, b, c, u, r, n)

      implicit none

      integer n, j, i
      double precision a(n), b(n), c(n), r(n), u(n)
      double precision bet, gam(4)
c
c	|b1	c1	0	0 | |u1|	|r1|
c	|a2	b2	c2	0 | |u2|	|r2|
c	|0	a3	b2	c2| |u3|   =	|r3|
c	|0	0	a4	b4| |u4|	|r4|
c
c
c      if(b(1).eq.0.)pause 'tridiag: rewrite equations'
c
      u(1) = 0.0d0
      u(2) = 0.0d0
      u(3) = 0.0d0
      u(4) = 0.0d0
c

      bet = b(1)
      u(1) = r(1)/bet
c
      do 11 j = 2, n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j) * gam(j)
        if(bet .eq. 0.d0) goto 999
        u(j) = (r(j) - a(j) * u(j-1))/bet
11    continue
c
c  back substitution
c
      do 12 j = n-1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
12    continue
      return
c
 999  continue
      print*,'BET IS ZERO AT TRIDIAG()'
      write(6, 99) (a(i), i = 1, n)
      write(6, 99) (b(i), i = 1, n)
      write(6, 99) (c(i), i = 1, n)
      print*,'gam(', j,')= ', gam(j), 'c(', j-1,')= ', c(j-1)
      print*,'a(', j,')= ', a(j), 'b(', j,')= ', b(j)
      stop 333
  99  format(4(1pe14.5))
      end

c .................................................................


      subroutine dgefa(a,lda,n,ipvt,info)

      implicit none

      integer lda,n,ipvt(n),info
      double precision a(lda,n)

c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

c .................................................................

      subroutine dgesl(a,lda,n,ipvt,b,job)
      
      implicit none

      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

c .................................................................

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none

      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

c .................................................................

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

c .................................................................

      subroutine dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none

      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

c .................................................................

      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
