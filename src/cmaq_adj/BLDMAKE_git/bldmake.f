!-----------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in    !
!  continuous development by various groups and is based on information !
!  from these groups: Federal Government employees, contractors working !
!  within a United States Government contract, and non-Federal sources  !
!  including research institutions.  These groups give the Government   !
!  permission to use, prepare derivative works of, and distribute copies!
!  of their work in the CMAQ system to the public and to permit others  !
!  to Do so.  The United States Environmental Protection Agency         !
!  therefore grants similar permission to use the CMAQ system software, !
!  but users are requested to provide copies of derivative works or     !
!  products designed to operate in the CMAQ system to the United States !
!  Government without restrictions as to use by others.  Software       !
!  that is used with the CMAQ system but distributed under the GNU      !
!  General Public License or the GNU Lesser General Public License is   !
!  subject to their copyright restrictions.                             !
!-----------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!     PROGRAM bldmake
!     Generate a Makefile for source files extracted from a CVS repository
!     Additional option for building Makefile in git repository
!     originally written in C by Steve Thorpe
!     rewritten in Fortran by Steve Howard
!     redone to meet CMAQ coding standards by Jeff Young (Nov 2012)
!-------------------------------------------------------------------------------

      Program bldmake

      Use ModelCFG

      Implicit None

      Character( EXT_LEN ) :: cfgFile
      Integer :: lfn = 11
      Integer status
      Integer n

      ! call Setup routine to process command line arguments
      Call setup( cfgFile )

      ! open cfgFile
      Open (unit=lfn, file=cfgFile, status='old', iostat=status)
      If ( status .Ne. 0 ) Then
        Write( *,'(" Open error number:",i5)') status
        Call error_msg( 'Cannot open FILE [' // Trim(cfgFile) //']' )
      End If

      ! read CFG file
      Call readCFG( lfn )
      Close (unit=lfn)

      ! CVS repository
      If ( repo .Eq. 'CVS' ) Then
        ! extract files for each module
        If ( checkout ) Then
          Call cvs_checkout( status )
        Else
          Call cvs_export( status )
        End If
      End If   !! CVS repository

      ! GIT repository
      If ( repo .Eq. 'GIT' ) Then
        Call git_export( status )
      End If

      ! create Makefile
      Open (unit=lfn, file='Makefile', iostat=status)
      If ( status .ne. 0 ) Call error_msg( 'Cannot create FILE [Makefile]' )

      Call makefile( lfn, cfgFile )

      Close (unit=lfn)

      ! If not makefo, Then run the make command to compile
      If ( .not. makefo ) Then
        Call RunMake( status )
      End If

      Stop
      End Program bldmake

!-------------------------------------------------------------------------------
!     Setup routine:  gets input file and run options from command line
!-------------------------------------------------------------------------------
      Subroutine setup( cfgFile )

      Use ModelCfg

      Implicit None

      ! arguments
      Character(*) cfgFile

      ! functions
      Integer :: IARGC

      ! local variables
      Integer :: status
      Integer :: nargs
      Integer :: n
      Character(32) :: argv

      ! date and time variables
      Character(8)  :: cdate
      Character(10) :: ctime
      Character(5)  :: czone
      Integer       :: dateValues(8)

      ! set defaults
      verbose   = .False.
      debug     = .False.
      checkout  = .False.
      makefo    = .False.
      git_local = .False.
      repo      = 'GIT'
      reporoot  = ' '
      cvsroot   = ' '

      ! check number of arguments on command line
      nargs = IARGC()  ! non-standard compiler extension returns the number of arguments
      If ( nargs .Eq. 0 ) Then
        Call help_msg('No arguments on command line')
        Stop
      End If

      ! get last argument (the bldit-created config file) and write to cfgFile
      Call GETARG( nargs, cfgFile )  ! non-standard compiler extension returns argument specified - the last argumet in this case
      If ( cfgFile(1:1) .Eq. '-' ) Then   ! config file not the last argument
        Call ucase( cfgFile )

        If ( cfgFile .Eq. '-HELP' ) Then
          Call help_msg(' ')
          Call cfgHelp( ); Stop
        End If

        Call help_msg('Invalid configuration file argument:' // Trim(cfgFile) )
        Stop
      End If

      ! check for run options
      Do n = 1, nargs-1
        Call GETARG( n, argv )

        If ( argv(1:1) .Ne. '-' ) Then
          Call help_msg('Invalid arguments on command line:' // Trim(argv) )
          Stop
        End If

        Call ucase( argv )

        If ( argv .Eq. '-HELP' ) Then
          Call help_msg('Help option:' // Trim(argv) )
          Call cfgHelp( ); Stop
        End If

        If ( argv .Eq. '-MAKEFO' ) Then     ! Make file only
          makefo = .True.; Cycle
        End If

        If ( argv .Eq. '-GIT_LOCAL' ) Then  ! do not copy source files to BLD directory
          git_local = .True.; Cycle
        End If

        If ( argv .Eq. '-CVS' ) Then
          repo = 'CVS'; Cycle
        End If

        If ( argv .Eq. '-CO' ) Then         ! use CVS "checkout" instead of "export"
          checkout = .True.; Cycle
        End If

        If ( argv .Eq. '-DEBUG' ) Then
          debug = .True.; Cycle
        End If

        If ( argv .Eq. '-VERBOSE' ) Then
          verbose = .True.; Cycle
        End If

        Call help_msg('Invalid arguments [' // Trim(argv) // '] on command line')

      End Do

      ! get CVSROOT
      Call GETENV('CVSROOT', cvsroot)

      ! If CVS Then CVSROOT must be defined 
      If ( repo .Eq. 'CVS' .and. Len_Trim(cvsroot) .Eq. 0 ) Then
        Write(*,'(/,"**ERROR** System variable [CVSROOT] not defined.",/)')
        Write(*,'(" ( System variable CVSROOT must be set to your cvs repository. )",/)')
        Stop
      End If

      ! If CVSROOT is set, Then set repo to CVS
      If ( Len_Trim(cvsroot) .Gt. 0 ) Then
        repo = 'CVS'
        If ( debug .or. verbose ) Then
          Write(*,'("CVSROOT set to:",a,/)') Trim( cvsroot )
        End If
      End If

      ! If REPOROOT is defined, use it, Else set to current directory
      If ( repo .Eq. 'GIT' ) Then
        Call GETENV('REPOROOT', reporoot)
        If ( Len_Trim(reporoot) .Eq. 0 ) Then
          Call PWD( reporoot, status )
          If ( status.ne.0 ) reporoot = './'
        End If
        If ( debug .or. verbose ) Then
          Write(*,'(''REPOROOT set to:'',a,/)') Trim( reporoot )
        End If
      End If  !! GIT Repository

      ! Get system date and time
      Call date_and_time( cdate, ctime, czone, dateValues )
      Write(currentDate, '(i2.2,"/",i2.2,"/",2i4.2,":",i2.2,":",i2.2)')
     &   dateValues(2), dateValues(3), dateValues(1),
     &   dateValues(5), dateValues(6), dateValues(7)

      Return
      End Subroutine setup

!-------------------------------------------------------------------------------
!     Help message:  Prints command line format and options and stops run
!-------------------------------------------------------------------------------
      Subroutine help_msg( msg )

      Implicit None

      ! arguments
      Character(*) msg

      Write(*,'(/,a)') Trim( msg )

      Write(*,'(/"Usage: bldmake [-<option>...] filename")')

      Write(*,'(/"where <option> is one of the following:")')
      Write(*,'("  -verbose   Echo actions")')
      Write(*,'("  -debug     Echo all actions")')
      Write(*,'("  -makefo    Creates Makefile without building")')
      Write(*,'("  -git_local Does NOT copy source files to BLD directory")')
      Write(*,'("  -cvs       Uses CVS repository")')
      Write(*,'("  -co        Uses CVS checkout instead of export")')
      Write(*,'("  -help      Displays help screen")')
      Write(*,'(//)')

      End Subroutine help_msg

!-------------------------------------------------------------------------------
!     Error:  Prints error string and stops run
!-------------------------------------------------------------------------------
      Subroutine error_msg( msg )

      Implicit None

      ! arguments
      Character(*) msg

      Write(*,'(/"*** Program terminated on Error ***"/)')
      Write(*,'(5x,a//)') Trim( msg )

      Stop
      End Subroutine error_msg

!-------------------------------------------------------------------------------
!     Makefile routine: creates Makefile from CFG data
!-------------------------------------------------------------------------------
      Subroutine makefile( lfn, cfgFile )

      Use ModelCfg

      Implicit None

      ! arguments
      Integer lfn
      Character(*) cfgFile

      If ( verbose ) Then
        Write(*,'(/"Generating Makefile"/)')
      End If

      ! print header lines
      Write(lfn,'("#   Makefile generated using program [bldmake]")')
      Write(lfn,'("#")')
      Write(lfn,'("#   Generation date [",a,"]")') Trim( currentDate )
      Write(lfn,'("#   Configuration file [",a,"]")') Trim( cfgFile )

      If ( repo .Eq. 'CVS' ) Then
        Write(lfn,'("#   Using CVS repository [",a,"]")') Trim( cvsroot )
      End If

      If ( repo .Eq. 'GIT' ) Then
        Write(lfn,'("#   Using GIT repository")')
      End If

      Write(lfn,'("#")')

      Write(lfn,'(/"MODEL = ",a)') Trim( model )

      Write(lfn,'(/"FC = ",a)') Trim( f_compiler )
      Write(lfn,'( "CC = ",a)') Trim( c_compiler )

      Write(lfn,'(/"f_FLAGS    = ",a)') Trim( f_flags )
      Write(lfn,'( "F_FLAGS    = ",a)') Trim( Fflags )
      Write(lfn,'( "f90_FLAGS  = ",a)') Trim( f90_flags )
      Write(lfn,'( "F90_FLAGS  = ",a)') Trim( F90flags )
      Write(lfn,'( "C_FLAGS    = ",a)') Trim( c_flags )

      If ( verbose ) Then
        Write(*,'("  Compilers defined")')
      End If

      Write(lfn,'(/"LINKER     = ",a)') Trim( linker )
      Write(lfn,'( "LINK_FLAGS = ",a)') Trim( link_flags )

      If ( Len_Trim(reporoot) .Gt. 0 ) Then
        Write(lfn,'(/"REPOROOT = ",a)') Trim( reporoot )
      End If

      If ( Len_Trim(VPATH) .Gt. 0 ) Then
        Call writeVPATH( lfn )
      End If

      Call writeCPP( lfn )
      If ( verbose ) Then
        Write(*,'("  CPP Flags defined")')
      End If

      Call writeLIB( lfn )
      If ( verbose ) Then
        Write(*,'("  Libraries defined")')
      End If

      Call writeINC( lfn )
      If ( verbose ) Then
        Write(*,'("  Includes defined")')
      End If

      Call writeOBJS( lfn )
      If ( verbose ) Then
        Write(*,'("  Objects defined")')
      End If

      Call writeRules( lfn )
      If ( verbose ) Then
        Write(*,'("  Make rules defined")')
      End If

      Call writeDEP( lfn )
      If ( verbose ) Then
        Write(*,'("  USE/MODULE dependencies defined")')
      End If

      Write(*,'(/"Makefile generated")')

      Return
      End Subroutine makefile

!-------------------------------------------------------------------------------
!     WriteVPATH routine:  Writes each directory on its own line
!-------------------------------------------------------------------------------
      Subroutine writeVPATH( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer lfn

      ! functions
      Integer getFieldCount

      ! local variables
      Integer nfields
      Integer n
      Character( EXT_LEN ) field

      Write(lfn,'(/,"#   Search PATH for source files")')

      nfields = getFieldCount( VPATH, ':' )

      Call getField( VPATH, ':', 1, field )
      Write(lfn,'("VPATH = ",a,":",$)') Trim(field)

      ! print each field at a time
      Do n = 2, nfields
        Call getField( VPATH, ':', n, field )
        If ( Len_Trim(field) .Gt. 0 )
     &     Write(lfn,'(1x,a,/,2x,a,":",$)') backslash, Trim( field )
      End Do

      Write(lfn,'(1x)')

      Return
      End Subroutine writeVPATH

!-------------------------------------------------------------------------------
!     WriteCPP routine:  Writes each '-D' on its own line
!-------------------------------------------------------------------------------
      Subroutine writeCPP( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer lfn

      ! functions
      Integer getFieldCount

      ! local variables
      Integer nfields
      Integer n
      Character( EXT_LEN ) :: field

      Write(lfn,'(/,''CPP = '',a)') Trim( cpp )

      nfields = getFieldCount( cpp_flags, ' ' )

      If ( nfields .Le. 1 ) Then
        Write(lfn,'(''CPP_FLAGS = '',a)') Trim( cpp_flags )
        Return
      End If

      Write(lfn,'(''CPP_FLAGS ='',$)')

      ! print each field at a time
      Do n = 1, nfields
        Call getField( cpp_flags, ' ', n, field )
        Write(lfn,'(1x,a,/,2x,a,$)') backslash, Trim( field )
      End Do

      Write(lfn,'(1x)')

      Return
      End Subroutine writeCPP

!-------------------------------------------------------------------------------
!     WriteLIB routine:  Writes libraries line to Makefile
!-------------------------------------------------------------------------------
      Subroutine writeLIB( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer lfn

      ! functions
      Integer getFieldCount

      ! local variables
      Integer nfields
      Integer n
      Integer pos
      Character( EXT_LEN ) :: field
      Character( EXT_LEN ) :: librec
      Character( EXT_LEN ) :: libname
      Character( EXT_LEN ) :: libs


      !!Write(lfn,'(/,"#   Library paths")')
      Write(lfn,'(1x)')

      nfields = getFieldCount( libraries, ' ' )

      If ( nfields .Le. 1 ) Then
        Write(lfn,'("LIBRARIES  = ",a)') Trim( libraries )
        Return
      End If

      libs = 'LIBRARIES ='

      ! parse library fields
      librec = ' '
      libname = ' '
      Do n = 1, nfields
        Call getField( libraries, ' ', n, field )

        If ( n .Gt. 1 .And. field(1:2) .Eq. '-L' ) Then
          If ( libname .Eq. ' ' ) Write(libname,'(''LIB'',i2.2)') n/2

          Write(lfn,'(a)') Trim(libname) // ' = ' // Trim(librec)
          libs = Trim(libs) // ' $(' // Trim(libname) // ')'
          librec = ' '
          libname = ' '
        End If

        librec = Trim( librec ) // ' ' // Trim( field )

        pos = Index( librec, '-l' )
        If ( libname .Eq. ' ' .And. pos .Gt. 0 ) Then
           libname = librec(pos+2:)
           Call ucase( libname )
        End If

      End Do

      If ( libname .Eq. ' ' ) Write(libname,'("LIB",i2.2)') n/2
      Write(lfn,'(a)') Trim(libname) // ' = ' // Trim(librec)
      libs = Trim(libs) // ' $(' // Trim(libname) // ')'

      !!Write(lfn,'(/,''#   Libraries'')')
      Write(lfn,'(/,a)') Trim(libs)

      Return
      End Subroutine writeLIB

!-------------------------------------------------------------------------------
!     WriteINC routine:  Writes include lines to Makefile
!-------------------------------------------------------------------------------
      Subroutine writeINC( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer lfn

      ! functions
      Integer getNumberOfFields

      ! local variables
      Integer n
      Integer i
      Integer j
      Integer k
      Integer pos
      Integer pos2
      Character( EXT_LEN ) :: base_inc
      Character( EXT_LEN ) :: mech_inc
      Character( EXT_LEN ) :: pa_inc

      Character( EXT_LEN ) :: pathName(3) = (/'BASE_INC', 'MECH_INC', 'PA_INC  '/)
      Character( EXT_LEN ) :: pathChk(3) =
     &                     (/'SUBST_CONST   ', 'SUBST_RXCMMN  ', 'SUBST_PACTL_ID'/)
      Integer              :: pathInd(3) = (/ 2, 1, 1 /)
      Character( EXT_LEN ) :: pathStr(3)
      Logical              :: hasPaths

      Character( EXT_LEN ) :: path

      If ( n_includes .Eq. 0 ) Return

      pathStr = ' '
      hasPaths = .False.

      ! find path strings
      Do n = 1, n_includes
        Do i = 1, Size( pathName )
          If ( include(n)%name .Eq. pathChk(i) ) Then
            path = include(n)%path
            Do j = 1, pathInd(i)
              pos = Index( path, '/', .True.)
              If ( pos .Le. 0 ) Exit
              path = path(1:pos-1)
            End Do    ! pathInd loop
            pathStr(i) = path
            hasPaths = .True.
          End If
        End Do  ! i loop
      End Do  ! n loop

      ! If paths found, write them
      If ( hasPaths ) Then
        Write(lfn,'(1x)')

        Do i = 1, Size( pathName )
          If ( pathStr(i) .Ne. ' ' ) Then
            Write(lfn,'(a," = ",a)') Trim(pathName(i)), Trim(pathStr(i))
          End If
        End Do
      End If     ! has paths

      ! write include lines
      Write(lfn,'(/"INCLUDES = ",$)')

      Do n = 1, n_includes
        path = include(n)%path

        Do i = 1, Size( pathName )
          If ( pathStr(i) .Ne. ' ' .And. Index( path, Trim(pathStr(i)) ) .Gt. 0 ) Then
            pos = Index( path, Trim(pathStr(i)) )
            pos2 = pos + Len_Trim( pathStr(i) )
            If ( pos .Gt. 1 ) path = path(1:pos-1) // '$('
     &                             // Trim(pathName(i)) // ')' // path(pos2:)
            If ( pos .Le. 1 ) path = '$('
     &                             // Trim(pathName(i)) // ')' // path(pos2:)
          End If
        End Do

        Write(lfn,'(1x,a,/,''  -D'',a,''='',a,''"'',a,a,''"'',$)') backslash,
     &       Trim(include(n)%name), backslash, Trim(path), backslash

        !! define include(n)%path2 to path with macro subsitutions
        include(n)%path2 = path

      End Do

      Write(lfn,'(1x)')

      Return
      End Subroutine writeINC

!-------------------------------------------------------------------------------
!     WriteOBJS routine:  Writes objects files by modules to Makefile
!-------------------------------------------------------------------------------
      Subroutine writeOBJS( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer lfn

      ! functions
      Integer getNumberOfFields

      ! local variables
      Integer :: nfiles
      Integer :: nfields
      Integer :: n
      Integer :: i
      Integer :: pos
      Character( FLD_LEN ) :: filename( MAX_FILES )
      Character( FLD_LEN ) :: modname
      Character( FLD_LEN ) :: obj
      Character( FLD_LEN ) :: objStr = ' '

      ! get list of all global modules
      Call orderfiles( module(1), .True., nfiles, filename )
      modname = 'GLOBAL_MODULES'

      If ( nfiles .Gt. 0 ) Then
        Write(lfn,'(/,a," =",$)') Trim( modname )
        objStr = '$(' // Trim(modname) // ')'
        Do i = 1, nfiles
          pos = Index(filename(i), '.')
          If ( pos .Le. 0 ) Cycle
          obj = filename(i)(1:pos) // 'o'
          Write(lfn,'(1x,a,/,2x,a,$)') backslash, Trim(obj)
        End Do
      End If

      ! loop thru each module and build list of objects
      Do n = 1, n_modules
        Call orderfiles( module(n), .False., nfiles, filename )

        If ( nfiles .Gt. 0 ) Then
          modname = module(n)%name
          Call ucase( modname )
          Write(lfn,'(//,a," =",$)') Trim(modname)
          objStr = Trim(objStr) // ' $(' // Trim(modname) // ')'
          Do i = 1, nfiles
            pos = Index(filename(i), '.')
            If ( pos .Le. 0 ) Cycle
            obj = filename(i)(1:pos) // 'o'
            Write(lfn,'(1x,a,/,2x,a,$)') backslash, Trim(obj)
          End Do
        End If
      End Do

      Write(lfn,'(//"OBJS =",$)')

      nfields = getNumberOfFields( objStr, ' ')
      Do n = 1, nfields
        Call getField( objStr, ' ', n, obj )
        Write(lfn,'(1x,a,/,2x,a,$)') backslash, Trim(obj)
      End Do

      Write(lfn,'(1x)')

      Return
      End Subroutine writeOBJS

!-------------------------------------------------------------------------------
!     WriteDEP routine:  Writes USE/MODULE and INCLUDE dependencies lines to Makefile
!-------------------------------------------------------------------------------
      Subroutine writeDEP( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer lfn

      ! functions
      Integer getNumberOfFields

      ! local variables
      Character( FLD_LEN ) :: objname
      Character( FLD_LEN ) :: reqStr 
      Character( 32 )      :: modName
      Character( 32 )      :: modFile
      Integer :: nfields
      Integer :: n
      Integer :: i
      Integer :: j
      Integer :: pos
      Character( 1 )       :: tab = char(9)

      Logical :: incdep( n_includes )
      Integer :: ndep
      Character( FLD_LEN ) :: depfile(MAX_FILES+n_includes)

      Write(lfn,'("# dependencies"/)')

      ! loop thru each archive module and write list of source file dependencies
      Do n = 1, n_modules
        Do i = 1, module(n)%nfiles
          ndep = 0
          depfile = ''

          ! build object filename
          objname = module(n)%file(i)%name
          pos = Index(objname, '.')
          If ( pos .Gt. 0 ) objname = objname(1:pos) // 'o'

          ! parse USEs string to get module names 
          nfields = getNumberOfFields( module(n)%file(i)%uses, ':' )
      write(*,*) ' '
          If ( nfields .Gt. 2 ) Then 
      write(*,*) 'file,nfields-1 ', Trim( module(n)%file(i)%name ), ' ', nfields-1
            Do j = 2, nfields-1
              Call getField( module(n)%file(i)%uses, ':', j, modName )
              Call getModFile( modName, modFile ) 
              If ( Len_Trim(modFile) .Gt. 0 ) Then 
      write(*,*) 'modName,modFile ', j, trim( modName ), '  ', trim( modFile )
                ndep = ndep+ 1
                depfile( ndep ) = modFile
      else
      write(*,*) 'modName         ', j, trim( modName ), '  -------------'
              End If
            End Do
          End If   ! has USEs files

          ! check for include file dependencies
          Call getIncDep( module(n)%file(i)%path, incdep )
          Do j = 1, n_includes
            If ( incdep(j) ) Then
              ndep = ndep+ 1
              depfile( ndep ) = include(j)%path2
            End If
          End Do
              
          ! write dependencies string
          reqStr = ''
          Do j = 1, ndep
             If ( Len_Trim(reqStr) .Eq. 0 ) Then
               reqStr = Trim(objName) // ':' // tab // depFile(j)
             Else If ( Len_Trim(reqStr) .Le. 60 ) Then
               reqStr = Trim(reqStr) // ' ' // depFile(j)
             Else
               Write(lfn,'(a,1x,a)') Trim(reqStr), backslash
               reqStr = tab // tab // depFile(j)
             End If
          End Do

          If ( Len_Trim(reqStr) .Gt. 0 ) Write(lfn,'(a)') Trim(reqStr) 

        End Do    ! file loop
      End Do    ! module loop

      Return
      End Subroutine writeDEP 

!-------------------------------------------------------------------------------
!     WriteRules routine:  Writes rules to Makefile
!-------------------------------------------------------------------------------
      Subroutine writeRules( lfn )

      Use ModelCFG

      Implicit None

      ! arguments
      Integer            :: lfn

      ! local variables
      Integer            :: n
      Character(FLD_LEN) :: record
      Character(1)       :: tab = char(9)


      ! build SUFFIXES record
      record = '.SUFFIXES:'

      Do n = 1, Size(extension)
        record = Trim(record) // ' ' // extension(n)
      End Do

      Write(lfn,'(/,a)') Trim( record )


      Write(lfn,'(/"$(MODEL): $(OBJS)")')
      Write(lfn,'(a,"$(LINKER) $(LINK_FLAGS) $(OBJS) $(LIBRARIES) -o $@"/)') tab

      Write(lfn,'(".F.o:")')
      Write(lfn,'(a,"$(FC) -c $(F_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<"/)') tab

      Write(lfn,'(".f.o:")')
      Write(lfn,'(a,"$(FC) -c $(f_FLAGS) $<"/)') tab

      Write(lfn,'(".F90.o:")')
      Write(lfn,'(a,"$(FC) -c $(F90_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<"/)') tab

      Write(lfn,'(".f90.o:")')
      Write(lfn,'(a,"$(FC) -c $(f90_FLAGS) $<"/)') tab

      Write(lfn,'(".c.o:")')
      Write(lfn,'(a,"$(CC) -c $(C_FLAGS) $<"/)') tab

      Write(lfn,'("clean:")')
      !!Write(lfn,'(a,"rm -f $(OBJS) $(MODEL) *.mod"/)') tab
      Write(lfn,'(a,"rm -f $(OBJS) *.mod"/)') tab

      Write(lfn,'(1x)')

      Return
      End Subroutine writeRules

