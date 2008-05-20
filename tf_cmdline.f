c this is <tf_cmdline.f>
cS
c------------------------------------------------------------------------------
c
c 24/06/97 by Thomas Forbriger (IfG Stuttgart)
c
c evaluate commandline arguments
c
c REVISIONS and CHANGES
c    24/06/97   V1.0   Thomas Forbriger
c    25/06/97   V1.1   now also gets correct last option on command line
c    07/07/99   V1.2   whooo - found major bug! optset was not cleared
c                      correctly - just optset(maxopt) was set to .false.
c    05/08/99    WF    declaration of ident to *(*) and argument to *6
c                      -> max length of recognized option string is 6
c
c==============================================================================
c
c  tf_cmdline
c
c  evaluate commandline input
c
c  optstart          first argument to be used
c  lastarg           last command-line argument used
c  ident(maxopt)     option identifiers
c  argu(maxopt)      returnvalues
c  optset(maxopt)    option is set
c  hasarg(maxopt)    are there arguments to read
c
      subroutine tf_cmdline(optstart, lastarg,
     &    maxopt, ident, argu, optset, hasarg)
c
c  declare parameters
c
      integer maxopt, optstart, lastarg
      character*(*) ident(maxopt)
      character*(*) argu(maxopt)
      logical optset(maxopt), hasarg(maxopt)
c
cE
c  declare variables
c
      integer opt,arg, iargc
      character*6 argument
c
      do 3 opt=1,maxopt
        optset(opt)=.FALSE.
    3 continue
      arg=optstart
      lastarg=optstart-1
    1 continue
        call getarg(arg, argument)
        opt=1
    2   continue
          if (argument.eq.ident(opt)) then
            optset(opt)=.TRUE.
            if (hasarg(opt)) then
              if ((arg+1).gt.iargc()) then
                print *,'WARNING (tf_cmdline): missing argument for option ',
     &                  ident(opt),' (defaults to ''none'')'
                argu(opt)='none'
              else
                call getarg(arg+1, argu(opt))
              endif
              arg=arg+1
            endif
            lastarg=arg
            opt=maxopt
          endif
          opt=opt+1
        if (opt.le.maxopt) goto 2
        arg=arg+1
      if (arg.le.iargc()) goto 1
      return
      end
c
c ----- END OF tf_cmdline.f -----
