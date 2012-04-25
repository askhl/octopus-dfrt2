# 1 "xc_ks_inversion.F90"
# 1 "/Users/umbe/Work/Projects/DFT/octopus-dfrt2/src/system//"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "xc_ks_inversion.F90"
!! Copyright (C) 2010 H. Appel, N. Helbig
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: $

# 1 "../../src/include/global.h" 1
!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: global.h 8345 2011-10-03 22:16:37Z dstrubbe $

# 1 "../../src/include/config_F90.h" 1
# 21 "../../src/include/global.h" 2


! If the compiler accepts long Fortran lines, it is better to use that
! capacity, and build all the preprocessor definitions in one line. In
! this way, the debuggers will provide the right line numbers.
# 35 "../../src/include/global.h"
! If the compiler accepts line number markers, then "CARDINAL" will
! put them. Otherwise, just a new line. Note that the "cardinal" and
! "newline" words are substituted by the program preprocess.pl by the
! ampersand and by a real new line just before compilation.
# 50 "../../src/include/global.h"
! The assertions are ignored if the code is compiled in not-debug mode (NDEBUG
! is defined). Otherwise it is merely a logical assertion that, when fails,
! prints out the assertion string, the file, and the line. The subroutine
! aassert_die is in the global_m module.
# 65 "../../src/include/global.h"
! Some compilers will not have the sizeof intrinsic.







! In octopus, one should normally use the SAFE_(DE)ALLOCATE macros below, which emit
! a helpful error if the the allocation or deallocation fails. The "MY_DEALLOCATE" macro
! is only used in this file; in the code, one should use SAFE_DEALLOCATE_P for pointers
! and SAFE_DEALLOCATE_A for arrays.
# 114 "../../src/include/global.h"
! This was used in the past and should not be used any more.



! The following macros facilitate the use of real or complex variables,
! and the possibility of compiling the code in single or double precision.
# 152 "../../src/include/global.h"
! The code directories should be defined here, and not hard coded in the Fortran files.
# 164 "../../src/include/global.h"
! The MPI1 and MPI2 standards are different regarding the MPI_IN_PLACE constant. In
! the code, just use the MPI_IN_PLACE_OR defined here.
# 175 "../../src/include/global.h"
! the TOSTRING macro converts a macro into a string
! do not use the STRINGIFY macro





! Whenever a procedure is not called too many times, one should start it
! and finish it with the PUSH_SUB and POP_SUB macros, which are these
! pieces of code that call the push_sub and pop_sub routines defined
! in the messages_m module.
# 195 "../../src/include/global.h"
! the leading dimension of the array


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
# 21 "xc_ks_inversion.F90" 2

module xc_ks_inversion_m
  use datasets_m
  use density_m
  use derivatives_m
  use eigensolver_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use multicomm_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use xc_f90_lib_m
  use xc_m
  use xc_functl_m

  implicit none

  private
  public :: &
    xc_ks_inversion_t, &
    xc_ks_inversion_init, &
    xc_ks_inversion_end, &
    xc_ks_inversion_write_info, &
    xc_ks_inversion_calc, &
    invertks_2part, &
    invertks_iter

  ! KS inversion methods/algorithms
  integer, public, parameter :: &
    XC_INV_METHOD_VS_ITER = 1, &
    XC_INV_METHOD_TWO_PARTICLE = 2

  ! the KS inversion levels
  integer, public, parameter :: &
    XC_KS_INVERSION_NONE = 1, &
    XC_KS_INVERSION_ADIABATIC = 2, &
    XC_KS_INVERSION_TD_EXACT = 3

  type xc_ks_inversion_t
     integer :: method
     integer :: level
     real(8), pointer :: vxc_previous_step(:,:)
     type(states_t) :: aux_st
     type(hamiltonian_t) :: aux_hm
     type(eigensolver_t) :: eigensolver
  end type xc_ks_inversion_t


contains

  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_init(ks_inv, family, gr, geo, d, mc)
    type(xc_ks_inversion_t), intent(out) :: ks_inv
    integer, intent(in) :: family
    type(grid_t), intent(inout) :: gr
    type(states_dim_t), intent(in) :: d
    type(geometry_t), intent(inout) :: geo
    type(multicomm_t), intent(in) :: mc

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_init"); 
 endif

    if(iand(family, XC_FAMILY_KS_INVERSION) .eq. 0) then
      ks_inv%level = XC_KS_INVERSION_NONE
      if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_init"); 
 endif
      return
    end if






    call messages_experimental("Kohn-Sham inversion")

    !%Variable InvertKSmethod
    !%Type integer
    !%Default iterative
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects whether the exact two-particle method or the iterative scheme
    !% is used to invert the density to get the KS potential.
    !%Option iterative 1
    !% Iterative scheme for v_s.
    !%Option two_particle 2
    !% Exact two-particle scheme.
    !%Option iterativevxc 3
    !% Iterative scheme for v_xc.
    !%End
    call parse_integer(datasets_check('InvertKSmethod'), &
            XC_INV_METHOD_VS_ITER, ks_inv%method)

    if(ks_inv%method < XC_INV_METHOD_VS_ITER &
      .or. ks_inv%method > XC_INV_METHOD_TWO_PARTICLE) then
      call input_error('InvertKSmethod')
      call messages_fatal(1)
    endif

    !%Variable KSInversionLevel
    !%Type integer
    !%Default ks_inversion_adiabatic
    !%Section Hamiltonian::XC
    !%Description
    !% At what level shall <tt>Octopus</tt> handle the KS inversion
    !%Option ks_inversion_none 1
    !% Do not compute KS inversion
    !%Option ks_inversion_adiabatic 2
    !% Compute exact adiabatic vxc
    !%End
    call messages_obsolete_variable('KS_Inversion_Level', 'KSInversionLevel')
    call parse_integer(datasets_check('KSInversionLevel'), XC_KS_INVERSION_ADIABATIC, ks_inv%level)
    if(.not.varinfo_valid_option('KSInversionLevel', ks_inv%level)) call input_error('KSInversionLevel')

    if(ks_inv%level.ne.XC_KS_INVERSION_NONE) then
      ! initialize auxilary random wavefunctions
      call states_null(ks_inv%aux_st)
      call states_init(ks_inv%aux_st, gr, geo)
      call states_allocate_wfns(ks_inv%aux_st, gr%mesh)
      call states_generate_random(ks_inv%aux_st, gr%mesh)
      ! initialize densities, hamiltonian and eigensolver
      call states_densities_init(ks_inv%aux_st, gr, geo, mc)
      call states_exec_init(ks_inv%aux_st, mc)
      call hamiltonian_init(ks_inv%aux_hm, gr, geo, ks_inv%aux_st, INDEPENDENT_PARTICLES, XC_FAMILY_NONE)
      call eigensolver_init(ks_inv%eigensolver, gr, ks_inv%aux_st)
    end if

    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_init"); 
 endif
  end subroutine xc_ks_inversion_init


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_end(ks_inv, gr, geo)
    type(xc_ks_inversion_t), intent(inout) :: ks_inv
    type(grid_t), intent(inout) :: gr
    type(geometry_t), intent(inout) :: geo

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_end"); 
 endif

    if(ks_inv%level .ne. XC_KS_INVERSION_NONE) then
      ! cleanup
      call eigensolver_end(ks_inv%eigensolver)
      call hamiltonian_end(ks_inv%aux_hm, gr, geo)
      call states_end(ks_inv%aux_st)
    end if

    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_end"); 
 endif
  end subroutine xc_ks_inversion_end


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_write_info(ks_inversion, iunit)
    type(xc_ks_inversion_t), intent(in) :: ks_inversion
    integer, intent(in) :: iunit

    if(ks_inversion%level.eq.XC_KS_INVERSION_NONE) return

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_write_info"); 
 endif
    call messages_print_var_option(iunit, 'KSInversionLevel', ks_inversion%level)

    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"xc_ks_inversion_write_info"); 
 endif
  end subroutine xc_ks_inversion_write_info


  ! ---------------------------------------------------------
  subroutine invertks_2part(target_rho, nspin, aux_hm, gr, st, eigensolver)

    type(grid_t), intent(in) :: gr
    type(states_t), intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: aux_hm
    type(eigensolver_t), intent(inout) :: eigensolver
    integer, intent(in) :: nspin
    real(8), intent(in) :: target_rho(1:gr%mesh%np, 1:nspin)

    integer :: ii, jj
    integer :: ndim, np
    real(8) :: spacing(1:3), stabilizer
    real(8), allocatable :: sqrtrho(:,:), laplace(:,:), vks(:,:)

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"invertks_2part"); 
 endif

    call parse_float(datasets_check('InvertKSStabilizer'), M_HALF, stabilizer)

    ndim = gr%sb%dim
    spacing = gr%mesh%spacing
    np = gr%mesh%np

    allocate(sqrtrho(1:gr%der%mesh%np_part, 1:nspin), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(sqrtrho(1:gr%der%mesh%np_part, 1:nspin)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "sqrtrho(1:gr%der%mesh%np_part, 1:nspin)", & 
 "xc_ks_inversion.F90", & 
 221, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 221); 
# 221 "xc_ks_inversion.F90"
    allocate(vks(1:np, 1:nspin), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(vks(1:np, 1:nspin)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "vks(1:np, 1:nspin)", & 
 "xc_ks_inversion.F90", & 
 222, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 222); 
# 222 "xc_ks_inversion.F90"
    allocate(laplace(1:gr%der%mesh%np, 1:nspin), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(laplace(1:gr%der%mesh%np, 1:nspin)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "laplace(1:gr%der%mesh%np, 1:nspin)", & 
 "xc_ks_inversion.F90", & 
 223, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 223); 
# 223 "xc_ks_inversion.F90"

    sqrtrho = 0.0_8

    do jj = 1, nspin
      do ii = 1, gr%der%mesh%np
        sqrtrho(ii, jj) = sqrt(target_rho(ii, jj))
      enddo
    enddo

    do jj = 1, nspin
      call dderivatives_lapl(gr%der, sqrtrho(:,jj), laplace(:,jj))
    enddo

    do jj = 1, nspin
      do ii = 1, np
          vks(ii, jj) = laplace(ii, jj)/(M_TWO*sqrtrho(ii, jj)) + st%eigenval(1,jj)
      enddo
    enddo

    do jj = 1, nspin
      aux_hm%vxc(:,jj) = vks(:,jj) - aux_hm%ep%vpsl(:) - aux_hm%vhartree(:)
      aux_hm%vhxc(:,jj) = aux_hm%vxc(:,jj) + aux_hm%vhartree(:)
    enddo

    call hamiltonian_update(aux_hm, gr%mesh)

    call eigensolver_run(eigensolver, gr, st, aux_hm, 1)

    call density_calc(st, gr, st%rho)

    if(allocated(sqrtrho)) then; 
 global_sizeof = sizeof(sqrtrho); 
 deallocate(sqrtrho, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("sqrtrho", & 
 "xc_ks_inversion.F90", & 
 254, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 254); 
# 254 "xc_ks_inversion.F90"; 
 end if
    if(allocated(laplace)) then; 
 global_sizeof = sizeof(laplace); 
 deallocate(laplace, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("laplace", & 
 "xc_ks_inversion.F90", & 
 255, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 255); 
# 255 "xc_ks_inversion.F90"; 
 end if
    if(allocated(vks)) then; 
 global_sizeof = sizeof(vks); 
 deallocate(vks, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("vks", & 
 "xc_ks_inversion.F90", & 
 256, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 256); 
# 256 "xc_ks_inversion.F90"; 
 end if

    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"invertks_2part"); 
 endif
  end subroutine invertks_2part


  ! ---------------------------------------------------------
  subroutine invertks_iter(target_rho, nspin, aux_hm, gr, st, eigensolver)
    type(grid_t), intent(in) :: gr
    type(states_t), intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: aux_hm
    type(eigensolver_t), intent(inout) :: eigensolver
    integer, intent(in) :: nspin
    real(8), intent(in) :: target_rho(1:gr%mesh%np, 1:nspin)

    integer :: ii, jj, ierr, idiffmax
    integer :: iunit, iunit2, verbosity, counter, np
    real(8) :: rr
    real(8) :: alpha, beta, convergence, alphascale, betascale
    real(8) :: stabilizer, convdensity, diffdensity, aa
    real(8), allocatable :: vhxc(:,:)

    character(len=20) :: alpha_str
    character(len=20) :: beta_str
    character(len=256) :: fname
    character(len=256) :: filename2

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"invertks_iter"); 
 endif

    np = gr%mesh%np

    !%Variable InvertKSConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Absolute difference between the calculated and the target density in the KS
    !% inversion. Has to be larger than the convergence of the density in the SCF run.
    !%End
    call parse_float(datasets_check('InvertKSConvAbsDens'), 1e-5_8, convdensity)

    !%Variable InvertKSStabilizer
    !%Type float
    !%Default 0.5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Additive constant <i>c</i> in the iterative calculation of the KS potential
    !% (v(alpha+1)=rho(alpha)+c)/(rho_target+c)*v(alpha)
    !% ensures that very small densities do not cause numerical problems.
    !%End
    call parse_float(datasets_check('InvertKSStabilizer'), M_HALF, stabilizer)

    !%Variable InvertKSVerbosity
    !%Type integer
    !%Default 0
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects what is output during the calculation of the KS potential.
    !%Option 0
    !% Only outputs the converged density and KS potential.
    !%Option 1
    !% Same as 0 but outputs the maximum difference to the target density in each
    !% iteration in addition.
    !%Option 2
    !% Same as 1 but outputs the density and the KS potential in each iteration in
    !% addition.
    !%End
    call parse_integer(datasets_check('InvertKSVerbosity'), 0, verbosity)
    if(verbosity < 0 .or. verbosity > 2) then
      call input_error('InvertKSVerbosity')
      call messages_fatal(1)
    endif

    allocate(vhxc(1:np, 1:nspin), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(vhxc(1:np, 1:nspin)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "vhxc(1:np, 1:nspin)", & 
 "xc_ks_inversion.F90", & 
 329, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 329); 
# 329 "xc_ks_inversion.F90"

    vhxc(1:np,1:nspin) = aux_hm%vhxc(1:np,1:nspin)

    if(verbosity == 1 .or. verbosity == 2) then
      iunit = io_open('InvertKSconvergence', action = 'write')
    endif

    diffdensity = 1.0_8
    counter = 0
    alpha = 0.1_8
    beta = 0.1_8
    convergence = 0.1_8
    alphascale = 10.0_8
    betascale = 100.0_8

    if(verbosity == 2) then
      write(alpha_str, '(f8.4)') alpha
      write(beta_str, '(f8.4)') beta
      filename2 = "diffdens_" // trim(adjustl(alpha_str)) // "_" // trim(adjustl(beta_str)) // ".dat"
      iunit2 = io_open(trim(filename2), action='write')
    end if

    do while(diffdensity > convdensity)

      counter = counter + 1

      if(verbosity == 2) then
        write(fname,'(i6.6)') counter
        call dio_function_output(io_function_fill_how("AxisX"), &
             ".", "vhxc"//fname, gr%mesh, aux_hm%vhxc(:,1), units_out%energy, ierr)
        call dio_function_output(io_function_fill_how("AxisX"), &
             ".", "rho"//fname, gr%mesh, st%rho(:,1), units_out%length**(-gr%sb%dim), ierr)
      endif

      call hamiltonian_update(aux_hm, gr%mesh)
      call eigensolver_run(eigensolver, gr, st, aux_hm, 1)
      call density_calc(st, gr, st%rho)

      ! Inversion according to Phys. Rev. Lett. 100, 153004 (2008), Eq. (6)
      do jj = 1, np
        call mesh_r(gr%mesh, jj, rr)
        if (abs(rr).lt. 1e-1_8) rr = 1e-1_8
        vhxc(jj, 1:nspin) = vhxc(jj, 1:nspin) + (st%rho(jj,1:nspin) - target_rho(jj,1:nspin))*alpha*rr**beta
      end do

      aux_hm%vhxc(1:np,1:nspin) = vhxc(1:np, 1:nspin)

      do jj = 1, nspin
        aux_hm%vxc(:,jj) = vhxc(:,jj) - aux_hm%vhartree(:)
      enddo

      diffdensity = 0.0_8
      do jj = 1, nspin
        do ii = 1, np
          if (abs(st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
            diffdensity = abs(st%rho(ii,jj)-target_rho(ii,jj))
            idiffmax=ii
          endif
        enddo
      enddo

      if (diffdensity .lt. convergence) then
        if(verbosity == 2) then
          write(iunit2,*) counter, diffdensity
        end if
        convergence = convergence / 10.0_8
      endif
      beta = min(1.0_8, diffdensity*betascale)
      alpha = 1.0_8 - diffdensity

      if(verbosity == 1 .or. verbosity == 2) then
        write(iunit,'(i6.6)', ADVANCE = 'no') counter
        write(iunit,'(es18.10)') diffdensity

        call flush(iunit)

      endif

    end do

    !calculate final density

    call hamiltonian_update(aux_hm, gr%mesh)
    call eigensolver_run(eigensolver, gr, st, aux_hm, 1)
    call density_calc(st, gr, st%rho)

    write(message(1),'(a,I8)') "Iterative KS inversion, iterations needed:", counter
    call messages_info(1)

    call io_close(iunit)
    if(verbosity == 2) then
      call io_close(iunit2)
    end if

    if(allocated(vhxc)) then; 
 global_sizeof = sizeof(vhxc); 
 deallocate(vhxc, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("vhxc", & 
 "xc_ks_inversion.F90", & 
 424, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 424); 
# 424 "xc_ks_inversion.F90"; 
 end if

    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"invertks_iter"); 
 endif

  end subroutine invertks_iter

  ! ---------------------------------------------------------
  subroutine precond_kiks(mesh, np, nspin, st, target_rho, vhxc_out)
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: np, nspin
    type(states_t), intent(in) :: st
    real(8), intent(in) :: target_rho(1:np, 1:nspin)
    real(8), intent(out) :: vhxc_out(1:np, 1:nspin,1:1)

    integer :: ip, iprime, ii, jj, ivec, jdim
    real(8) :: numerator, diffrho, epsij, occij, inverse
    real(8) :: vol_element
    real(8) :: ki(1:np, 1:np)
    real(8) :: eigenvals(1:np), inverseki(1:np,1:np)
    real(8), allocatable :: matrixmul(:,:), kired(:,:)
    real(8), allocatable :: psii(:, :), psij(:, :)

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"precond_kiks"); 
 endif

    numerator = 0.0_8
    vhxc_out = 0.0_8

    !do ip = 1, np
    ! diffrho = st%rho(ip, 1) - target_rho(ip, 1)
    ! numerator = numerator + diffrho**2
    !end do

    ki = 0.0_8

    allocate(psii(1:mesh%np, 1:st%d%dim), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(psii(1:mesh%np, 1:st%d%dim)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "psii(1:mesh%np, 1:st%d%dim)", & 
 "xc_ks_inversion.F90", & 
 458, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 458); 
# 458 "xc_ks_inversion.F90"
    allocate(psij(1:mesh%np, 1:st%d%dim), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(psij(1:mesh%np, 1:st%d%dim)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "psij(1:mesh%np, 1:st%d%dim)", & 
 "xc_ks_inversion.F90", & 
 459, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 459); 
# 459 "xc_ks_inversion.F90"

    do jj = 1, st%nst

      call states_get_state(st, mesh, jj, 1, psij)

      do ii = jj + 1, st%nst

        call states_get_state(st, mesh, ii, 1, psii)

        epsij = 1.0_8 / (st%eigenval(jj, 1) - st%eigenval(ii, 1))
 occij = st%occ(jj, 1) - st%occ(ii, 1)
        do iprime = 1, np
          do ip = 1, np
            ki(ip, iprime) = ki(ip, iprime) + occij*epsij*(psii(ip, 1)*psij(ip, 1))*(psii(iprime, 1)*psij(iprime, 1))
   end do
        end do

      end do
    end do

    if(allocated(psii)) then; 
 global_sizeof = sizeof(psii); 
 deallocate(psii, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("psii", & 
 "xc_ks_inversion.F90", & 
 480, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 480); 
# 480 "xc_ks_inversion.F90"; 
 end if
    if(allocated(psij)) then; 
 global_sizeof = sizeof(psij); 
 deallocate(psij, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("psij", & 
 "xc_ks_inversion.F90", & 
 481, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 481); 
# 481 "xc_ks_inversion.F90"; 
 end if

    call lalg_eigensolve(np, ki, eigenvals)

    !do ip = 1, np
    ! if(abs(eigenvals(ip))>1d-10) then
    ! neigenval = ip
    ! endif
    !enddo

    allocate(matrixmul(1:np, 1:st%nst), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(matrixmul(1:np, 1:st%nst)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "matrixmul(1:np, 1:st%nst)", & 
 "xc_ks_inversion.F90", & 
 491, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 491); 
# 491 "xc_ks_inversion.F90"
    allocate(kired(1:np, 1:st%nst), stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & 
 global_sizeof = sizeof(kired(1:np, 1:st%nst)); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_allocate(& 
 "kired(1:np, 1:st%nst)", & 
 "xc_ks_inversion.F90", & 
 492, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call alloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 492); 
# 492 "xc_ks_inversion.F90"

    do ivec = 1, st%nst
      inverse = 1.0_8/eigenvals(ivec)
      do ip = 1, np
        matrixmul(ip, ivec) = ki(ip,ivec)*inverse
 kired(ip, ivec) = ki(ip, ivec)
      enddo
    enddo

    inverseki = matmul(matrixmul, transpose(kired))

    vhxc_out = 0.0_8

    vol_element = 1.0_8
    do jdim = 1, 3
      if (mesh%spacing(jdim) > 1.e-10_8) vol_element = vol_element*mesh%spacing(jdim)
    end do

    do iprime = 1, np
      diffrho = target_rho(iprime, 1) - st%rho(iprime, 1)
      do ip = 1, np
 vhxc_out(ip, 1, 1) = vhxc_out(ip, 1, 1) + inverseki(ip, iprime)*diffrho
 write(200,*) ip, iprime, inverseki(ip, iprime)
      enddo
    enddo

    do ip = 1, np
      write(100,*) ip, vhxc_out(ip, 1, 1), target_rho(ip, 1) - st%rho(ip, 1)
    enddo


    if(allocated(matrixmul)) then; 
 global_sizeof = sizeof(matrixmul); 
 deallocate(matrixmul, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("matrixmul", & 
 "xc_ks_inversion.F90", & 
 524, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 524); 
# 524 "xc_ks_inversion.F90"; 
 end if
    if(allocated(kired)) then; 
 global_sizeof = sizeof(kired); 
 deallocate(kired, stat=global_alloc_err); 
 if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & 
 call profiling_memory_deallocate("kired", & 
 "xc_ks_inversion.F90", & 
 525, & 
 global_sizeof); 
 if(global_alloc_err.ne.0) & 
 call dealloc_error(global_sizeof, & 
 "xc_ks_inversion.F90", & 
 525); 
# 525 "xc_ks_inversion.F90"; 
 end if


    !call flush(200)


    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"precond_kiks"); 
 endif

  end subroutine precond_kiks


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_calc(ks_inversion, gr, hm, st, ex, ec, vxc, time)
    type(xc_ks_inversion_t), intent(inout) :: ks_inversion
    type(grid_t), intent(inout) :: gr
    type(hamiltonian_t), intent(in) :: hm
    type(states_t), intent(inout) :: st
    real(8), intent(inout) :: ex, ec
    real(8), intent(inout) :: vxc(:,:) ! vxc(gr%mesh%np, st%d%nspin)
    real(8), optional, intent(in) :: time

    integer :: ii

    if(ks_inversion%level == XC_KS_INVERSION_NONE) return

    if(in_debug_mode) then; 
 call push_sub("xc_ks_inversion.F90"//"." & 
 //"X(xc_ks_inversion_calc)"); 
 endif

    call density_calc(st, gr, st%rho)

    if(present(time)) then
      write(message(1),'(A,F18.12)') 'xc_ks_inversion_calc - time:', time
      call messages_info(1)
    end if

    ks_inversion%aux_hm%energy%intnvxc = 0.0_8
    ks_inversion%aux_hm%energy%hartree = 0.0_8
    ks_inversion%aux_hm%energy%exchange = 0.0_8
    ks_inversion%aux_hm%energy%correlation = 0.0_8

    ks_inversion%aux_hm%vhartree = hm%vhartree

    if (present(time) .and. time > 0.0_8) then
! write(*,*) 'debug 1'
      do ii = 1, st%d%nspin
        ks_inversion%aux_hm%vhxc(:,ii) = ks_inversion%vxc_previous_step(:,ii)
      enddo
    else
      do ii = 1, st%d%nspin
        ks_inversion%aux_hm%vxc(:,ii) = 0.0_8 !hm%ep%vpsl(:)
        ks_inversion%aux_hm%vhxc(:,ii) = ks_inversion%aux_hm%vhartree(:) + ks_inversion%aux_hm%vxc(:,ii)
      enddo
    end if
    !ks_inversion%aux_hm%ep%vpsl(:) = 0.0_8 ! hm%ep%vpsl(:)
    ks_inversion%aux_hm%ep%vpsl(:) = hm%ep%vpsl(:)

    ! compute ks inversion, vhxc contains total KS potential

    select case (ks_inversion%method)
    ! adiabatic ks inversion
    case(XC_INV_METHOD_TWO_PARTICLE)
      call invertks_2part(ks_inversion%aux_st%rho, st%d%nspin, ks_inversion%aux_hm, gr, &
                         ks_inversion%aux_st, ks_inversion%eigensolver)
    case(XC_INV_METHOD_VS_ITER)
      call invertks_iter(st%rho, st%d%nspin, ks_inversion%aux_hm, gr, &
                         ks_inversion%aux_st, ks_inversion%eigensolver)
    end select

    !subtract external and Hartree potentials, ATTENTION: subtracts true external potential not adiabatic one

    do ii = 1, st%d%nspin
      ks_inversion%aux_hm%vxc(:,ii) = ks_inversion%aux_hm%vhxc(:,ii) - hm%vhartree(:)
    enddo

    vxc = ks_inversion%aux_hm%vxc

    if (present(time)) then
! write(*,*) 'debug 2'
      ks_inversion%vxc_previous_step = ks_inversion%aux_hm%vhxc
! write(*,*) 'debug 3'
    end if

    if(in_debug_mode) then; 
 call pop_sub("xc_ks_inversion.F90"//"." & 
 //"X(xc_ks_inversion_calc)"); 
 endif

  end subroutine xc_ks_inversion_calc


end module xc_ks_inversion_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
