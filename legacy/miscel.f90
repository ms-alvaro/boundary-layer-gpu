      !> \file miscel.f
      !!
      !!  File with miscellaneous functions and subroutines
      !!
      !! \date 04-10-2017 by M.Moriche \n
      !!       Updated to use with gfortran
        
      !! NOTE: this file was writen originally by M. Moriche for IBcode
      !! and it has been adapted by G. Arranz to be used in PCcode
      !!
      !! All credits to M. Moriche

      !> subroutine to set blanks in character strings
      !!
      !! tipically to initialize strings
      subroutine blank(string,n)
      implicit none
      character(1), parameter:: val=' '
      integer n
      character(n) string
      integer i
      do i=1,n
         string(i:i) = val
      enddo
      end subroutine

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! IO routines
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !!! sub:getvar beg !!!
      subroutine getvar(envvarname,string)

      use params, only : paramsfilename, IO_TMP

      implicit none
    
      character (LEN=* ) :: envvarname
      character (LEN=256) :: string
      character (LEN=256) :: line_str
      integer            :: eq_pos=1
      logical            :: exists
      
      string=' '

      !.....Check parameter file for variable

      inquire(FILE=trim(adjustl(paramsfilename)),EXIST=exists)

      if (exists) then
      open(IO_TMP,FILE=trim(adjustl(paramsfilename)),STATUS='OLD',ERR=1)

      do while (.true.)
        read(IO_TMP,'(A)',ERR=1,END=1) line_str
        eq_pos=SCAN(line_str,'=')

        if (TRIM(ADJUSTL(line_str(1:eq_pos-1))).eq.envvarname) then
          string=TRIM(line_str(eq_pos+1:))
          eq_pos=SCAN(string,'!')
          if (eq_pos.ne.0) then
              string = trim(string(1:eq_pos-1))
          endif
          exit
        endif

      enddo
1     close(IO_TMP)

      endif
      endsubroutine
      !!! sub:getvar end !!!
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !!! sub:get_dbl beg !!!
      subroutine get_dbl(envvarname,var,found)

      implicit none

      character (LEN= *) :: envvarname
      real(8)            :: var
      logical            :: found
      character (LEN=256) :: string

      call getvar(envvarname,string)
      if (string.ne.' ') then
        read  (string,*) var
        found=.true.
      else
        found=.false.
      endif

      endsubroutine
      !!! sub:get_dbl end !!!
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !!! sub:get_dbl_arr beg !!!
      subroutine get_dbl_arr(envvarname,var,n,found)

      implicit none

      character (LEN= *) :: envvarname
      integer            :: n
      real(8)            :: var(n)
      logical            :: found
      character (LEN=256) :: string
      ! aux
      integer i
      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,*) var(1:n)
        found=.true.
      else
        found=.false.
      endif

      endsubroutine
      !!! sub:get_dbl_arr end !!!
      
      !!! sub:get_int_arr2 beg !!!
      subroutine get_dbl_arr2(envvarname,var,m,n,found)

      implicit none

      character (LEN= *) :: envvarname
      integer            :: m,n
      real(8)            :: var(m,n)
      logical            :: found
      character (LEN=256) :: string

      integer i, j, k
      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,*) var
        found=.true.
      else
        found =.false.
      endif

      endsubroutine
      !!! sub:get_dbl_arr2 end !!!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !!! sub:get_int beg !!!
      subroutine get_int(envvarname,var,found)


      implicit none

      character (LEN= *) :: envvarname
      integer            :: var
      logical            :: found
      character (LEN=256) :: string

      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,*) var
        found=.true.
      else
        found =.false.
      endif

      endsubroutine
      !!! sub:get_int end !!!

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !!! sub:get_int_arr beg !!!
      subroutine get_int_arr(envvarname,var,n,found)

      implicit none

      character (LEN= *) :: envvarname
      integer            :: n
      integer            :: var(n)
      logical            :: found
      character (LEN=256) :: string

      integer i
      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,*) var
        found=.true.
      else
        found =.false.
      endif


      endsubroutine
      !!! sub:get_int_arr end !!!

      !!! sub:get_int_arr2 beg !!!
      subroutine get_int_arr2(envvarname,var,m,n,found)

      implicit none

      character (LEN= *) :: envvarname
      integer            :: m,n
      integer            :: var(m,n)
      logical            :: found
      character (LEN=256) :: string

      integer i, j
      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,*) var
        found=.true.
      else
        found =.false.
      endif

      endsubroutine
      !!! sub:get_int end !!!
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !!! sub:get_str beg !!!
      subroutine get_str(envvarname,var,n,found)


      implicit none

      character (LEN= *) :: envvarname
      integer n
      character (LEN= *) :: var
      logical            :: found
      character (LEN=256) :: string

      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,'(a)') var
        found=.true.
      else
        found=.false.
      endif

      endsubroutine
      !!! sub:get_str end !!!

      !!! sub:get_str0 beg !!!
      subroutine get_str0(envvarname,var,n,found)


      implicit none

      character (LEN= *) :: envvarname
      integer n
      character (LEN= *) :: var
      logical            :: found
      character (LEN=256) :: string

      call getvar(envvarname,string)

      if (string.ne.' ') then
        read  (string,'(a)') var
        found=.true.
      else
        found=.false.
      endif

      endsubroutine
      !!! sub:get_str end !!!
