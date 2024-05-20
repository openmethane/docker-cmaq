C //////////////////////////////////////////////////////////////////
C  Subroutine GETPAR calculates the 3rd moments (M3), masses, aerosol
C  densities, and geometric mean diameters (Dg) of all 3 modes, and the 
C  natural logs of geometric standard deviations (Sg) of the 
C  Aitken and accumulation modes.
 
C  The input logical variable, M3_WET_FLAG, dictates whether the 
C  calculations in GETPAR are to assume that the aerosol is "wet" or
C  "dry."  In the present context, a "wet" aerosol consists of all
C  chemical components of the aerosol.  A "dry" aerosol excludes 
C  particle-bound water and also excludes secondary organic aerosol.
 
C  NOTE! 2nd moment concentrations (M2) are passed into GETPAR in the 
C  CBLK array and are modified within GETPAR only in the event that
C  the Sg value of a given mode has gone outside of the acceptable
C  range (1.05 to 2.50).  The GETPAR calculations implicitly assume
C  that the input value of M2 is consistent with the input value of
C  M3_WET_FLAG.  If, for example, the input M2 value was calculated
C  for a "dry" aerosol and the M3_WET_FLAG is .TRUE., GETPAR would
C  incorrectly adjust the M2 concentrations!
 
C-----------------------------------------------------------------------
      Subroutine getpar( m3_wet_flag, limit_sg  )

      Use aero_data
      Use met_data

      Implicit None

C Arguments:
      Logical, Intent( In ) :: m3_wet_flag ! true = include H2O and SOA in 3rd moment
                                           ! false = exclude H2O and SOA from 3rd moment

      Logical, Intent( In ) :: limit_sg  ! fix coarse and accum Sg's to the input value?

C Output variables:
C  updates arrays in aero_data module
C  moment3_conc   3rd moment concentration [ ug /m**3 ]
C  aeromode_mass  mass concentration: [ ug / m**3 ]
C  aeromode_dens  avg particle density [ kg / m**3 ]
C  aeromode_diam  geometric mean diameter [ m ]
C  aeromode_sdev  log of geometric standard deviation

C Local Variables:
      Real( 8 ) :: xxm0        ! temporary storage of moment 0 conc's
      Real( 8 ) :: xxm2        ! temporary storage of moment 2 conc's
      Real( 8 ) :: xxm3        ! temporary storage of moment 3 conc's
      Real( 8 ) :: xfsum       ! (ln(M0)+2ln(M3))/3; used in Sg calcs
      Real( 8 ) :: lxfm2       ! ln(M2); used in Sg calcs
      Real( 8 ) :: l2sg        ! square of ln(Sg); used in diameter calcs
      Real      :: es36        ! exp(4.5*l2sg); used in diameter calcs
      Real      :: m3augm      ! temp variable for wet 3rd moment calcs

      Real( 8 ), Parameter :: one3d = 1.0D0 / 3.0D0
      Real( 8 ), Parameter :: two3d = 2.0D0 / 3.0D0

      Real,      Parameter :: one3  = 1.0 / 3.0
      Real,      Parameter :: dgmin = 1.0E-09   ! minimum particle diameter [ m ]
      Real,      Parameter :: densmin = 1.0E03  ! minimum particle density [ kg/m**3 ]

      Real( 8 ) :: minl2sg( n_mode )   ! min value of ln(sg)**2 for each mode 
      Real( 8 ) :: maxl2sg( n_mode )   ! max value of ln(sg)**2 for each mode

      Real      :: factor
      Integer   :: n, spc   ! loop counters

C-----------------------------------------------------------------------

C *** Set bounds for ln(Sg)**2

      Do n = 1 , n_mode
        If ( limit_sg .eqv. .True. ) Then
          minl2sg( n ) = aeromode_sdev( n ) ** 2
          maxl2sg( n ) = aeromode_sdev( n ) ** 2
        Else
          minl2sg( n ) = Log( min_sigma_g ) ** 2
          maxl2sg( n ) = Log( max_sigma_g ) ** 2
        EndIf
      End Do

C *** Calculate aerosol 3rd moment concentrations [ m**3 / m**3 ]

      Do n = 1, n_mode
         moment3_conc( n )  = 0.0
         aeromode_mass( n ) = 0.0

         Do spc = 1, n_aerospc
            If ( (aerospc( spc )%name( n ) .Ne. ' ' ) .And.
     &         ( .Not. aerospc( spc )%iswet .Or. m3_wet_flag) ) Then
               factor = 1.0E-9 * f6dpi / aerospc( spc )%density
               moment3_conc( n )  = moment3_conc( n ) + ( factor * aerospc_conc( spc,n ) )
               aeromode_mass( n ) = aeromode_mass( n ) + aerospc_conc( spc,n )
            End If
         End Do

      End Do

C *** Calculate modal average particle densities [ kg m**-3 ]

      Do n = 1, n_mode    
        aeromode_dens( n ) = Max( Real( densmin,8 ),
     &                      1.0E-9 * f6dpi * aeromode_mass( n ) /moment3_conc( n )  )
      End Do

C *** Calculate geometric standard deviations as follows:
c        ln^2(Sg) = 1/3*ln(M0) + 2/3*ln(M3) - ln(M2)
c     NOTES: 
c      1. Equation 10-5a of [Binkowski:1999] and Equation 5a of 
c         Binkowski&Roselle(2003) contain typographical errors.
c      2. If the square of the logarithm of the geometric standard
c         deviation is out of an acceptable range, reset this value and
c         adjust the second moments to be consistent with this value.
c         In this manner, M2 is artificially increased when Sg exceeds
c         the maximum limit.  M2 is artificially decreased when Sg falls
c         below the minimum limit.

C *** Aitken Mode:

      Do n = 1, n_mode
         xxm0 = moment0_conc( n ) 
         xxm2 = moment2_conc( n ) 
         xxm3 = moment3_conc( n ) 

         xfsum = one3d * Log( xxm0 ) + two3d * Log( xxm3 )

         lxfm2 = Log( xxm2 )
         l2sg =  xfsum - lxfm2

         l2sg = Max( l2sg, minl2sg( n ) )
         l2sg = Min( l2sg, maxl2sg( n ) )

         lxfm2 = xfsum - l2sg
         moment2_conc( n )  = Exp ( lxfm2 )
         aeromode_sdev( n ) = Sqrt( l2sg )

         ES36 = Exp( 4.5 * l2sg )
         aeromode_diam( n ) = Max( dgmin, ( moment3_conc( n )
     &                      / ( moment0_conc( n ) * es36 ) ) ** one3 )

      End Do


      Return
      End Subroutine getpar

      
