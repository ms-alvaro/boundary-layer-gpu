      !> 
      !!
      !! NOTE: this file was writen originally by M. Moriche for IBcode
      !! and it has been adapted by G. Arranz to be used in PCcode
      !!
      !! All credits to M. Moriche
      !!
      !!
      !! \details
      !! 
      !! contents of the file:
      !!  - module params
      !!  - subroutine initparams()

      !> \author M.Moriche
      !!
      !! \date 27-02-2014 by M.Moriche \n
      !!       Created
      !! \date 04-09-2015 by M.Moriche \n
      !!       Modified to consolidate input
      !! \date 09-09-2015 by A.Gonzalo \n
      !!       Set value of cflmax
      !! \date 17-03-2016 by A.Gonzalo \n
      !!       Added Dimlessnums array (Equivalent of Re number for ps)
      !!       Removed uref
      !! \date 29-10-2019 by A.Gonzalo \n
      !!       Added cfl2max and cfldiffmax that are used when
      !!       FPP_NONNEWTONIAN > 0
      !! \date 08-11-2019 by A.Gonzalo \n
      !!       Removed cfldiffmax that was used when
      !!       FPP_NONNEWTONIAN > 0
      !! \date 21-04-2020 by J.M. Catalan
      !!       Added FPP_TSIG
      !! \date 10-06-2020 by A.Gonzalo \n
      !!       Added foufilename and qmitfilename for tucanH simulations
      !! \date 15-06-2020 by M.Guerrero \n
      !!       Added ndhneumann for Neumann Lag bc implementation
      !! \date 23-06-2020 by M.Guerrero \n
      !!       Added mLoop for multi-forcing algorithm implementation on
      !!       IB points
      !! \brief Module with code parameters (INPUT)
      !!
      !! \details
      !!
      !! Input has been consolidated with general IO routines:
      !! - get_double
      !! - get_double_array
      !! - get_int
      !! - get_string
      !!
      !! Contains every variable that must be an input for the code.
      !!
      !!
      module params

      use global, only : Int32, Int64

      integer, parameter:: IO_TMP = 10 !< temporal file unit 

      integer(Int32), parameter :: ndim = 3 !< number of dimensions

      integer, parameter:: baselength=256   !< character strings length
      integer, parameter:: filelength=32768 !< character file length
                                            !! (400 lines of 80 chars)

      character(baselength) paramsfilename !< File with the name of the parameters file

      integer nxyz(ndim)            !< number of grid points in x, y, z
      real(int64)  boxsize(ndim+1)  !< Length of domain and stretching in y-dir
                                    !< [Lx Ly Lz alpha]
      real(int64) alphas(ndim+3)    !< alpha_mean_[xyz] + std and freq_mult
                                    !< slip lenghts for robin
      endmodule params

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      subroutine initparams()

      use params
      use global
      !use ifport, only : MAKEDIRQQ  ! removed for gfortran compatibility

      implicit none

      logical :: f, f1
      integer :: bl, fl
      integer :: i, Nout
      character(len=32) :: arg

      bl = baselength
      fl = filelength

      ! read name of parameters file 
      f = .FALSE. ! init flag
      i = 0
      do
        call get_command_argument(i, arg)
        if (len_trim(arg) == 0) exit

        select case (arg)
            case ('-i')
                call get_command_argument(i+1,paramsfilename)
                f = .true. ! input_file found 
            end select
        i = i+1
      end do
      if (f.eqv..FALSE.) then
          stop 'No input file provided. Use *.exe -i <input_file>'
      endif

      INQUIRE(FILE= paramsfilename, EXIST= f )
      if ( f.eqv..false. ) then
          write(*,*) 'input file: ', &
             adjustl(trim(paramsfilename)), ' does not exist'
          stop
      endif

      call get_dbl('nu'           , nu              , f)
      call get_dbl('CFL'          , CFL             , f)

      call get_dbl('Vmax'         , Vbs_max         , f)
      call get_dbl('x0'           , x_bs            , f)
      call get_dbl('sigma'        , sigma_bs        , f)
      call get_dbl('phi'          , phi_bs          , f)

      call get_int('nsteps'       , nsteps          , f)
      call get_int('nsave'        , nsave           , f)
      call get_int('nstats'       , nstats          , f)
      call get_int('nmonitor'     , nmonitor        , f)

      call get_int_arr('nxyz'     , nxyz , ndim     , f)
      call get_dbl_arr('boxsize'  , boxsize, ndim+1 , f)

      call get_int('inflow_flag'  , inflow_boundary_flag, f)
      call get_int('top_flag'     , top_boundary_flag, f)

      call get_str('inflow_file'     , file_inflow         , 200, f)
      if (f.eqv..FALSE.) then
          stop ' ERROR: you must specify a BL profile in inflow_file '
      endif


      call get_int('Lund_ix'      , i_rescale       , f)
      call get_dbl('Lund_deltai'  , delta_inlet     , f)
      call get_dbl('Lund_T'       , T_resc          , f)


      call get_dbl('dPdx'         , dPdx            , f)
      call get_dbl('dPdz'         , dPdz            , f)

      !call get_str('LES'          , LES_model_str   ,f)
      call get_int('LES'          , LES_model       , f)
      call get_int('WM'           , iwall_model     , f)
      call get_int('TauwModel'    , istress_model   , f)
      call get_int('nutBC'        , Dirichlet_nu_t  , f)

      ! wall model
      call get_int('WMnutflag'    , iwall_model_nut,  f)
      call get_int('WMnut'        , frac_vis_wall_model,  f)

      call get_dbl('Amplitude'    , Amplitude_perturbations, f)

      call get_int('init_step'    , nstep_init_input, f)
      call get_int('init_rand'    , random_init     , f)
      call get_str('filein'       , filein,      200, f1)
      if ( (random_init.eq.1) .and. (f1.eqv..TRUE.) ) then
          stop ' ERROR: init_rand = 1, but filein exists'
      endif
      ! Assign input number if given
      if (nstep_init_input.ne.-45) Then
         nstep_init = nstep_init_input
      endif
  
      call get_str('fileout'      , fileout,     200, f)
        
      ! Create a folder for the parent directory of the file if it does not exist
      Nout = len(trim(adjustl(fileout))) ! len of fileout string

      i = Nout
      f = .false.
      do while (f.eqv..false.)
         if (fileout(i:i).eq.'/') then
             f = .true.
         endif
         i = i - 1
      end do
        
      write(*,*) fileout, Nout, fileout(1:i)

      inquire( FILE=trim( adjustl( fileout(1:i) ) ), EXIST=f1 )
      if (f1.eqv..FALSE.) then
          call system( 'mkdir -p '//trim( adjustl( fileout(1:i) ) ) )
      endif
    
      call get_str('timeinflow_file' , file_temporal_inlet, 200, f)
    
      call get_int('RKscheme'        , itime_step   , f)        

      call get_dbl_arr('alphas'   , alphas, ndim+2 , f)
      if (f.eqv..TRUE.) then
          alpha_mean_x = alphas(1)
          alpha_mean_y = alphas(2)
          alpha_mean_z = alphas(3)
          alpha_std    = alphas(4)
          freq_mult    = alphas(5)
      endif

      nx_global = nxyz(1)
      ny_global = nxyz(2)
      nz_global = nxyz(3)

      Lx_rand    = boxsize(1)
      Ly_rand    = boxsize(2)
      Lz_rand    = boxsize(3)
      alpha_rand = boxsize(4)
    

      endsubroutine 

