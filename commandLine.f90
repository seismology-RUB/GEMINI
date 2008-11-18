!--------------------------------------------------------------------
!  Command line module
!
!  maxopt:  total number of possible options
!  maxmanarg:  total number of mandatory arguments
!  optid:   array of option names
!  optarg:  array of arguments to options
!  hasarg:  array of 1 and 0 indicating whether an option has an argument
!  optset:  array of 1 and 0 indicating which option was set on the command line
!  lastarg: index of last argument that belongs to the options
!           next argument is a mandatory one
!
!  for initialization use strings "opts" and "mask" to enter
!  a list of options and hasarg-flags:
!  opts:    string filled with blank separated options without minus signs
!  mask:    string with blank separated hasarg-flags
!--------------------------------------------------------------------
 module commandLine
	private :: clInit
	interface new
		module procedure clInit
	end interface
	interface dealloc
		module procedure clDealloc
	end interface
	type cmdLine
		private
		integer :: maxopt,maxmanarg,lastarg
		integer, dimension(:), pointer :: hasarg => null()
		integer, dimension(:), pointer :: optset => null()
		character (len=6), dimension(:), pointer :: optid => null()
		character (len=132), dimension(:), pointer :: optarg => null()
		character (len=132), dimension(:), pointer :: manarg => null()
	end type cmdLine
!
	contains
!-------------------------------------------------------------------------------
	subroutine clInit(this,maxopt,maxmanarg,opts,mask,printhelp)
	implicit none
	type (cmdLine) :: this
	integer :: maxopt, maxmanarg
	character (len=*) :: opts, mask
	external printhelp
	integer :: i,jarg,lastarg,jopt
	integer :: iargc
	character (len=6) argopt
!
	allocate(this%optid(maxopt),this%hasarg(maxopt),this%optarg(maxopt),&
	       & this%optset(maxopt))
	if (maxmanarg > 0) allocate(this%manarg(maxmanarg))
!
	read(opts,*,end=11) (this%optid(i),i=1,maxopt)
	read(mask,*,end=12) (this%hasarg(i),i=1,maxopt)
!
!  process command line
!
	do i=1,maxopt
		this%optset(i)=0
	enddo
	jarg=1
	lastarg=0
	do while (jarg <= iargc())           ! loop over optional command line arguments
		call getarg(jarg,argopt)
		jopt=1
		do while (jopt <= maxopt)     ! check which option has been set
			if(argopt == '-'//this%optid(jopt)) then
				this%optset(jopt)=1
				if(this%hasarg(jopt) == 1) then
					if(jarg+1 > iargc()) then
						print *,'commandLIne: Argument for option ',this%optid(jopt),' is missing!'
						call printhelp
						stop
					endif
					call getarg(jarg+1,this%optarg(jopt))
					jarg=jarg+1
				endif
				lastarg=jarg
				exit
			endif
			jopt=jopt+1
		enddo
		jarg=jarg+1
	enddo
	this%lastarg=lastarg
!
!  get mandatory arguments only if -h is not set and iargc() > 0
!
	if(this%optset(1) == 0 .and. iargc() > 0) then
		do i=1,maxmanarg
			if(i+lastarg > iargc()) then
				print *,'commandLine: not enough mandatory arguments'
				call printhelp
				stop
			endif
			call getarg(i+lastarg,this%manarg(i))
		enddo
	endif
!
!  if -h is set, call printhelp and stop
!
	if(this%optset(1) == 1) then
		call printhelp
		stop
	endif
!
!  if iargc() == 0, trigger printhelp and stop
!
	if(iargc() == 0 .and. maxmanarg > 0) then
		call printhelp
		stop
	endif
	return
!
 11	print *,'commandLine: maxopt and option string are inconsistent'
 	stop
 12	print *,'commandLine: maxopt and hasarg string are inconsistent'
	stop
	end subroutine clInit
!-------------------------------------------------------------------------------
!  delete allocated pointer
!
	subroutine clDealloc(this)
	type (cmdLine) :: this
!
	if(associated(this%optid)) deallocate(this%optid)
	if(associated(this%hasarg)) deallocate(this%hasarg)
	if(associated(this%optarg)) deallocate(this%optarg)
	if(associated(this%optset)) deallocate(this%optset)
	if(associated(this%manarg)) deallocate(this%manarg)
	end subroutine clDealloc
!-------------------------------------------------------------------------------
!  return optset(i)
!
	logical function clOptset(this,n)
	type (cmdLine) :: this
	integer :: n
	if(this%optset(n) == 1) then
		clOptset = .true.
	else
		clOptset = .false.
	endif
	end function clOptset
!-------------------------------------------------------------------------------
!  return optarg(i)
!
	character (len=132) function clOptarg(this,n)
	type (cmdLine) :: this
	integer :: n
	clOptarg = this%optarg(n)
	end function clOptarg
!-------------------------------------------------------------------------------
!  return manarg(i)
!
	character (len=132) function clManarg(this,n)
	type (cmdLine) :: this
	integer :: n
	clManarg = this%manarg(n)
	end function clManarg
!-------------------------------------------------------------------------------
  end module commandLine
