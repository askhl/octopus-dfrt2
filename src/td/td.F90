!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: td.F90 9915 2013-01-31 23:03:15Z dstrubbe $

#include "global.h"

module td_m
  use batch_m
  use calc_mode_m
  use cpmd_m
  use datasets_m
  use density_m
  use energy_calc_m
  use epot_m
  use forces_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use ground_state_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use ion_dynamics_m
  use kick_m
  use lasers_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_m
  use messages_m
  use mpi_m
  use multicomm_m
  use parser_m
  use PES_m
  use profiling_m
  use projector_m
  use restart_m
  use scf_m
  use simul_box_m
  use species_m
  use spectrum_m
  use states_m
  use states_calc_m
  use states_dim_m
  use states_io_m
  use system_m
  use propagator_m
  use td_write_m
  use types_m
  use unit_m
  use unit_system_m
  use v_ks_m
  use varinfo_m

  implicit none

  private
  public ::               &
    td_t,                 &
    td_run,               &
    td_run_init,          &
    td_init,              &
    td_end

  ! Parameters.
  integer, parameter :: &
       EHRENFEST = 1,   &
       BO        = 2,   &
       CP        = 3

  type td_t
    type(propagator_t)   :: tr             !< contains the details of the time-evolution
    type(scf_t)          :: scf
    type(ion_dynamics_t) :: ions
    type(cpmd_t)         :: cp_propagator
    FLOAT                :: dt             !< time step
    integer              :: max_iter       !< maximum number of iterations to perform
    integer              :: iter           !< the actual iteration
    logical              :: recalculate_gs !< Recalculate ground-state along the evolution.

    type(PES_t)          :: PESv

    FLOAT                :: mu
    integer              :: dynamics
    integer              :: energy_update_iter
  end type td_t


contains

  subroutine td_run_init()

    PUSH_SUB(td_run_init)

    call calc_mode_set_parallelization(P_STRATEGY_STATES, default = .true.)

    POP_SUB(td_run_init)
  end subroutine td_run_init

  ! ---------------------------------------------------------
  subroutine td_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(td_t)                :: td
    type(td_write_t)          :: write_handler
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo
    logical                   :: stopping, update_energy, cmplxscl
    integer                   :: iter, ii, ierr, scsteps, ispin
    real(8)                   :: etime
    logical                   :: generate
    type(gauge_force_t)       :: gauge_force
    type(profile_t),     save :: prof
    FLOAT, allocatable :: vold(:, :), Imvold(:, :)

    PUSH_SUB(td_run)

    cmplxscl = hm%cmplxscl%space

    ! some shortcuts
    gr  => sys%gr
    geo => sys%geo
    st  => sys%st

    if(simul_box_is_periodic(gr%mesh%sb)) call messages_experimental('Time propagation for periodic systems')

    call td_init(td, sys, hm)

    ! Allocate wavefunctions during time-propagation
    if(td%dynamics == EHRENFEST) then
      !complex wfs are required for Ehrenfest
      call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX, alloc_Left = cmplxscl)
      if(st%open_boundaries) then
        ASSERT(associated(gr%ob_grid%lead))
        call states_allocate_intf_wfns(st, gr%ob_grid%lead(:)%mesh)
      end if
    else
      call states_allocate_wfns(st, gr%mesh, alloc_Left = cmplxscl)
    end if

    ! CP has to be initialized after wavefunction type is set
    if(td%dynamics == CP) call cpmd_init(td%cp_propagator, sys%gr, sys%st)

    call init_wfs()

    ! Calculate initial forces and kinetic energy
    if(ion_dynamics_ions_move(td%ions)) then
      if(td%iter > 0) then
        call td_read_coordinates()
        call hamiltonian_epot_generate(hm, gr, geo, st, time = td%iter*td%dt)
      end if

      call forces_calculate(gr, geo, hm%ep, st, td%iter*td%dt)

      geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
    end if

    ! Calculate initial value of the gauge vector field
    call gauge_field_init(hm%ep%gfield, gr%sb)

    if (gauge_field_is_applied(hm%ep%gfield)) then

      if(td%iter > 0) then
        call td_read_gauge_field()
      else
        call gauge_field_init_vec_pot(hm%ep%gfield, gr%sb, st)
      end if

      call hamiltonian_update(hm, gr%mesh, time = td%dt*td%iter)

      call gauge_field_get_force(gr, geo, hm%ep%proj, hm%phase, st, gauge_force)
    end if

    call td_write_init(write_handler, gr, st, hm, geo, &
         ion_dynamics_ions_move(td%ions), gauge_field_is_applied(hm%ep%gfield), hm%ep%kick, td%iter, td%max_iter, td%dt)

    if(td%iter == 0) call td_run_zero_iter()

    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call messages_print_stress(stdout, "Time-Dependent Simulation")
    call print_header()

    if(td%PESv%calc_rc .or. td%PESv%calc_mask) then
       if (fromScratch) then
          call PES_init_write(td%PESv,gr%mesh,st)
       else
          call PES_restart_read(td%PESv, st)
       endif
    endif

    if(st%d%pack_states .and. hamiltonian_apply_packed(hm, gr%mesh)) call states_pack(st)

    ii = 1
    etime = loct_clock()
    ! This is the time-propagation loop. It starts at t=0 and finishes
    ! at td%max_iter*dt. The index i runs from 1 to td%max_iter, and
    ! step "iter" means propagation from (iter-1)*dt to iter*dt.
    propagation: do iter = td%iter, td%max_iter

      stopping = clean_stop(sys%mc%master_comm)
      call profiling_in(prof, "TIME_STEP")

      if(iter > 1) then
        if( ((iter-1)*td%dt <= hm%ep%kick%time) .and. (iter*td%dt > hm%ep%kick%time) ) then
          if(.not. cmplxscl) then
            call kick_apply(gr, st, td%ions, geo, hm%ep%kick)
          else
            call kick_apply(gr, st, td%ions, geo, hm%ep%kick, hm%cmplxscl%theta)
          end if
          call td_write_kick(gr, hm, sys%outp, geo, iter)
        end if
      end if

      !Apply mask absorbing boundaries
      if(hm%ab == MASK_ABSORBING) call zvmask(gr, hm, st) 

      ! time iterate wavefunctions
      select case(td%dynamics)
      case(EHRENFEST)
        if(ion_dynamics_ions_move(td%ions)) then
          call propagator_dt(sys%ks, hm, gr, st, td%tr, iter*td%dt, td%dt, td%mu, td%max_iter, iter, gauge_force, &
            ions = td%ions, geo = sys%geo, scsteps = scsteps)
        else
          call propagator_dt(sys%ks, hm, gr, st, td%tr, iter*td%dt, td%dt, td%mu, td%max_iter, iter, gauge_force, &
            geo=sys%geo, scsteps = scsteps)
        end if
      case(BO)
        ! move the hamiltonian to time t
        call ion_dynamics_propagate(td%ions, sys%gr%sb, sys%geo, iter*td%dt, td%dt)
        call hamiltonian_epot_generate(hm, gr, sys%geo, st, time = iter*td%dt)
        ! now calculate the eigenfunctions
        call scf_run(td%scf, sys%mc, sys%gr, geo, st, sys%ks, hm, sys%outp, &
          gs_run = .false., verbosity = VERB_COMPACT, iters_done = scsteps)
      case(CP)
        if(states_are_real(st)) then
          call dcpmd_propagate(td%cp_propagator, sys%gr, hm, st, iter, td%dt)
        else
          call zcpmd_propagate(td%cp_propagator, sys%gr, hm, st, iter, td%dt)
        end if
        scsteps = 1
      end select

      ! update density
      if(.not. propagator_dens_is_propagated(td%tr)) then 
        if(.not. cmplxscl) then
          call density_calc(st, gr, st%rho)
        else
          call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
        end if  
      end if
      
      generate = .false.

      if(ion_dynamics_ions_move(td%ions)) then
        if(td%dynamics == CP .or. .not. propagator_ions_are_propagated(td%tr)) then
          call ion_dynamics_propagate(td%ions, sys%gr%sb, sys%geo, iter*td%dt, td%dt)
          generate = .true.
        end if
      end if

      if(gauge_field_is_applied(hm%ep%gfield) .and. .not. propagator_ions_are_propagated(td%tr)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, td%dt)
      end if

      if(generate .or. geometry_species_time_dependent(geo)) then
        call hamiltonian_epot_generate(hm, gr, sys%geo, st, time = iter*td%dt)
      end if

      update_energy = (td%dynamics == BO) .or. (mod(iter, td%energy_update_iter) == 0) .or. (iter == td%max_iter)

      if(update_energy .or. propagator_requires_vks(td%tr)) then

        ! save the vhxc potential for later
        if(.not. propagator_requires_vks(td%tr)) then
          SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
          do ispin = 1, st%d%nspin
            call lalg_copy(gr%mesh%np, hm%vhxc(:, ispin), vold(:, ispin))
          end do
          if(cmplxscl) then
            SAFE_ALLOCATE(Imvold(1:gr%mesh%np, 1:st%d%nspin))
            do ispin = 1, st%d%nspin
              call lalg_copy(gr%mesh%np, hm%Imvhxc(:, ispin), Imvold(:, ispin))
            end do
          end if
        end if
   


        ! update Hamiltonian and eigenvalues (fermi is *not* called)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval = update_energy, time = iter*td%dt, calc_energy = update_energy)

        ! Get the energies.
        if(update_energy) call energy_calc_total(hm, sys%gr, st, iunit = -1)

        if (td%dynamics == CP) then
          if(states_are_real(st)) then
            call dcpmd_propagate_vel(td%cp_propagator, sys%gr, hm, st, td%dt)
          else
            call zcpmd_propagate_vel(td%cp_propagator, sys%gr, hm, st, td%dt)
          end if
        end if

        ! restore the vhxc
        if(.not. propagator_requires_vks(td%tr)) then
          do ispin = 1, st%d%nspin
            call lalg_copy(gr%mesh%np, vold(:, ispin), hm%vhxc(:, ispin))
          end do
          SAFE_DEALLOCATE_A(vold)
          if(cmplxscl) then
            do ispin = 1, st%d%nspin
              call lalg_copy(gr%mesh%np, Imvold(:, ispin), hm%Imvhxc(:, ispin))
            end do
            SAFE_DEALLOCATE_A(Imvold)
          end if
        end if
      end if

      ! Recalculate forces, update velocities...
      if(ion_dynamics_ions_move(td%ions)) then
        if(td%dynamics /= BO) call forces_calculate(gr, sys%geo, hm%ep, st, iter*td%dt)

        call ion_dynamics_propagate_vel(td%ions, sys%geo, atoms_moved = generate)

        if(generate) call hamiltonian_epot_generate(hm, gr, sys%geo, st, time = iter*td%dt)

        geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        if(td%dynamics /= BO) call gauge_field_get_force(gr, geo, hm%ep%proj, hm%phase, st, gauge_force)
        call gauge_field_propagate_vel(hm%ep%gfield, gauge_force, td%dt)
      end if

      !Photoelectron stuff 
      if(td%PESv%calc_rc .or. td%PESv%calc_mask ) &
           call PES_calc(td%PESv, gr%mesh, st, ii, td%dt, iter)

      call td_write_iter(write_handler, gr, st, hm, geo, hm%ep%kick, td%dt, iter)

      ! write down data
      call check_point()

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(prof)
      if (stopping) exit

    end do propagation

    if(st%d%pack_states .and. hamiltonian_apply_packed(hm, gr%mesh)) call states_unpack(st)

    call td_write_end(write_handler)
    call end_()

#ifdef HAVE_MPI
    ! wait for all processors to finish
    if(st%parallel_in_states) then
      call MPI_Barrier(st%mpi_grp%comm, mpi_err)
    end if
#endif

    POP_SUB(td_run)

  contains

    subroutine print_header

      if(td%dynamics /= CP) then
        if(.not.cmplxscl) then
          write(message(1), '(a7,1x,a14,a14,a10,a17)') 'Iter ', 'Time ', 'Energy ', 'SC Steps', 'Elapsed Time '
        else
          write(message(1), '(a7,1x,a14,a14,a14,a10,a17)') &
                                      'Iter ', 'Time ', 'Re(Energy) ','Im(Energy) ', 'SC Steps', 'Elapsed Time '
        end if
      else
        write(message(1), '(a7,1x,a14,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'CP Energy ', 'Elapsed Time '
      end if
      call messages_info(1)
      call messages_print_stress(stdout)

    end subroutine print_header

    ! ---------------------------------------------------------
    subroutine check_point
      PUSH_SUB(td_run.check_point)

      ! write info
      if(td%dynamics /= CP) then
        if(.not. cmplxscl) then
          write(message(1), '(i7,1x,2f14.6,i10,f14.3)') iter, &
               units_from_atomic(units_out%time, iter*td%dt), &
               units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
               scsteps, loct_clock() - etime
        else
          write(message(1), '(i7,1x,3f14.6,i10,f14.3)') iter, &
               units_from_atomic(units_out%time, iter*td%dt), &
               units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
               units_from_atomic(units_out%energy, hm%energy%Imtotal), &
               scsteps, loct_clock() - etime
        end if     
      else
        write(message(1), '(i7,1x,3f14.6,f14.3, i10)') iter, &
             units_from_atomic(units_out%time, iter*td%dt), &
             units_from_atomic(units_out%energy, hm%energy%total + geo%kinetic_energy), &
             units_from_atomic(units_out%energy, &
                hm%energy%total + geo%kinetic_energy + cpmd_electronic_energy(td%cp_propagator)), &
             loct_clock() - etime
      end if
      call messages_info(1)
      etime = loct_clock()
      ii = ii + 1
      if(ii==sys%outp%iter+1 .or. iter == td%max_iter .or. stopping) then ! output
        if(iter == td%max_iter) sys%outp%iter = ii - 1
        ii = 1
        call td_save_restart(iter)
        call td_write_data(write_handler, gr, st, hm, sys%ks%xc, sys%outp, geo, iter, td%dt)
	!Photoelectron output and restart dump
        call PES_output(td%PESv, gr%mesh, st, iter, sys%outp, td%dt,gr,geo)
        call PES_restart_write(td%PESv, st)
        if( (ion_dynamics_ions_move(td%ions)) .and. td%recalculate_gs) then
          call messages_print_stress(stdout, 'Recalculating the ground state.')
          fromScratch = .false.
          call ground_state_run(sys, hm, fromScratch)
          call restart_read(trim(restart_dir)//'td', st, gr, ierr, iter=iter)
          call messages_print_stress(stdout, "Time-dependent simulation proceeds")
          call print_header()
        end if
      end if

      POP_SUB(td_run.check_point)
    end subroutine check_point

   ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(td_run.end_)

      ! free memory
      if(td%dynamics == CP) call cpmd_end(td%cp_propagator)
      call states_deallocate_wfns(st)
      call ion_dynamics_end(td%ions)
      call td_end(td)

      POP_SUB(td_run.end_)
    end subroutine end_

    ! ---------------------------------------------------------
    subroutine init_wfs()

      integer :: i, is, ierr, ist, ik, jst, freeze_orbitals
      character(len=50) :: filename
      FLOAT :: x
      type(block_t) :: blk
      type(states_t) :: stin
      CMPLX, allocatable :: rotation_matrix(:, :), zv_old(:), zpsi(:,:)

      PUSH_SUB(td_run.init_wfs)

      if(.not.fromscratch) then
        call restart_read(trim(tmpdir)//'td', st, gr, ierr, iter=td%iter, read_left = st%have_left_states)
        if(ierr.ne.0) then
          message(1) = "Could not load "//trim(tmpdir)//"td: Starting from scratch"
          call messages_warning(1)

          fromScratch = .true.
          td%iter = 0
        end if
        ! extract the interface wave function
        if(st%open_boundaries) call restart_get_ob_intf(st, gr)
      end if

      if(td%iter >= td%max_iter) then
        message(1) = "All requested iterations have already been done. Use FromScratch = yes if you want to redo them."
        call messages_info(1)
      endif

      if(.not. fromscratch .and. td%dynamics == CP) then
        call cpmd_restart_read(td%cp_propagator, gr, st, ierr)

        if(ierr.ne.0) then
          message(1) = "Could not load "//trim(restart_dir)//"td/cpmd: Starting from scratch"
          call messages_warning(1)

          fromScratch = .true.
          td%iter = 0
        end if
      end if

      if(.not. fromscratch) then
        ! read potential from previous interactions
        
        if(cmplxscl) then
          SAFE_ALLOCATE(zv_old(1:gr%mesh%np))
        end if
        
        do i = 1, 2
          do is = 1, st%d%nspin
            write(filename,'(a,i2.2,i3.3)') trim(tmpdir)//'td/vprev_', i, is
            if(cmplxscl) then
              call zio_function_input(trim(filename)//'.obf', gr%mesh, zv_old(1:gr%mesh%np), ierr)
              td%tr%v_old(1:gr%mesh%np, is, i)   =  real(zv_old(1:gr%mesh%np))
              td%tr%Imv_old(1:gr%mesh%np, is, i) = aimag(zv_old(1:gr%mesh%np))
            else
              call dio_function_input(trim(filename)//'.obf', gr%mesh, td%tr%v_old(1:gr%mesh%np, is, i), ierr)
            end if
            ! If we do not succeed, try netcdf
            if(ierr > 0) call dio_function_input(trim(filename)//'.ncdf', gr%mesh, td%tr%v_old(1:gr%mesh%np, is, i), ierr)
            if(ierr > 0) then
              write(message(1), '(3a)') 'Unsuccessful read of "', trim(filename), '"'
              call messages_fatal(1)
            end if
          end do
        end do
        if(cmplxscl) then
          SAFE_DEALLOCATE_A(zv_old)
        end if

      end if

      if(fromScratch) then
        if(.not. st%only_userdef_istates) then
          call restart_read(trim(restart_dir)//GS_DIR, st, gr, ierr, exact = .true.)
          if(ierr.ne.0) then
            write(message(1), '(3a)') 'Unsuccessful read of states.'
            call messages_fatal(1)
          end if
          ! extract the interface wave function
          if(st%open_boundaries) call restart_get_ob_intf(st, gr)
        end if

        ! check if we should deploy user-defined wavefunctions.
        ! according to the settings in the input file the routine
        ! overwrites orbitals that were read from restart/gs
        if(parse_isdef(datasets_check('UserDefinedStates')).ne.0) then
          call restart_read_user_def_orbitals(gr%mesh, st)
        end if
        
!         if(st%have_left_states) then
!           ! At the beginning of the time evolution left and right states are the same 
!           ! Rational: At the moment cmplxscl left and right states are symmetric 
!           ! since the original (unscaled) Hamiltonian is real.
!           ! In the future we may consider to extend the ground state scf
!           ! to deal with the general case (this would be needed, for instance, to include
!           ! magnetic fields) 
!           cSAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:st%d%dim))
!           do ik = st%d%kpt%start, st%d%kpt%end
!             do ist = st%st_start, st%st_end
!               call states_get_state(st, gr%mesh, ist,  ik, zpsi)
!               call states_set_state(st, gr%mesh, ist,  ik, zpsi, left = .true.)
!             end do
!           end do    
!           cSAFE_DEALLOCATE_A(zpsi)       
!         end if
        

        !%Variable TransformStates
        !%Type block
        !%Default no
        !%Section States
        !%Description
        !% Before starting the <tt>td</tt> calculation, the initial states (that are
        !% read from the <tt>restart/gs</tt> directory, which should have been
        !% generated in a previous ground-state calculation) can be "transformed"
        !% among themselves. The block <tt>TransformStates</tt> gives the transformation matrix
        !% to be used. The number of rows of the matrix should equal the number
        !% of the states present in the time-dependent calculation (the independent
        !% spin and <i>k</i>-point subspaces are all transformed equally); the number of
        !% columns should be equal to the number of states present in the
        !% <tt>restart/gs</tt> directory. This number may be different: for example,
        !% one could have run previously in <tt>unocc</tt> mode in order to obtain unoccupied
        !% Kohn-Sham states, and therefore <tt>restart/gs</tt> will contain more states.
        !% These states can be used in the transformation.
        !%
        !% Note that the code will not check the orthonormality of the new states!
        !%
        !% Each line provides the coefficients of the new states, in terms of
        !% the old ones.
        !%End
        if(parse_isdef(datasets_check('TransformStates')).ne.0) then
          if(parse_block(datasets_check('TransformStates'), blk) == 0) then
            call states_copy(stin, st)
            SAFE_DEALLOCATE_P(stin%zpsi)
            call restart_look_and_read(stin, gr)
            SAFE_ALLOCATE(rotation_matrix(1:st%nst, 1:stin%nst))
            rotation_matrix = M_z0
            do ist = 1, st%nst
              do jst = 1, parse_block_cols(blk, ist-1)
                call parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            call states_rotate(gr%mesh, st, stin, rotation_matrix)
            SAFE_DEALLOCATE_A(rotation_matrix)
            call states_end(stin)
          else
            message(1) = '"TransformStates" has to be specified as block.'
            call messages_info(1)
            call input_error('TransformStates')
          end if
        end if

      end if



      !%Variable TDFreezeOrbitals
      !%Type integer
      !%Default 0
      !%Section Time-Dependent
      !%Description
      !% You have the possibility of "freezing" a number of orbitals during a time-propagation.
      !% The Hartree and exchange-correlation potential due to these orbitals (which
      !% will be the lowest-energy ones) will be added during the propagation, but the orbitals
      !% will not be propagated.
      !%
      !% <b>WARNING: NOT TESTED YET.</b>
      !%Option sae -1
      !% Single-active-electron approximation. This option is only valid for time-dependent
      !% calculations (<tt>CalculationMode = td</tt>). Also, the nuclei should not move.
      !% The idea is that all orbitals except the last one are frozen. The orbitals are to
      !% be read from a previous ground-state calculation. The active orbital is then treated
      !% as independent (whether it contains one electron or two) -- although it will
      !% feel the Hartree and exchange-correlation potentials from the ground-state electronic
      !% configuration.
      !%
      !% It is almost equivalent to setting <tt>TDFreezeOrbitals = N-1</tt>, where <tt>N</tt> is the number
      !% of orbitals, but not completely.
      !%End
      call parse_integer(datasets_check('TDFreezeOrbitals'), 0, freeze_orbitals)

      if(freeze_orbitals /= 0) call messages_experimental('TDFreezeOrbitals')

      if(.not. cmplxscl) then
        call density_calc(st, gr, st%rho)
      else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
      end if

      if(freeze_orbitals > 0) then
        ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
        call states_freeze_orbitals(st, gr, sys%mc, freeze_orbitals)
        write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
          ' orbitals have been frozen.', st%nst, ' will be propagated.'
        call messages_info(1)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      elseif(freeze_orbitals < 0) then
        ! This means SAE approximation. We calculate the Hxc first, then freeze all
        ! orbitals minus one.
        write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
        call messages_info(1)
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
        call states_freeze_orbitals(st, gr, sys%mc, n = st%nst-1)
        call v_ks_freeze_hxc(sys%ks)
        call density_calc(st, gr, st%rho)
      else
        ! Normal run.
        call v_ks_calc(sys%ks, hm, st, sys%geo, calc_eigenval=.true., time = td%iter*td%dt)
      end if

      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_Bcast(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, mpi_err)
      end if
#endif
      call hamiltonian_span(hm, minval(gr%mesh%spacing(1:gr%mesh%sb%dim)), x)
      call energy_calc_total(hm, gr, st)

      POP_SUB(td_run.init_wfs)
    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      PUSH_SUB(td_run.td_run_zero_iter)

      call td_write_iter(write_handler, gr, st, hm, geo, hm%ep%kick, td%dt, 0)

      ! I apply the delta electric field *after* td_write_iter, otherwise the
      ! dipole matrix elements in write_proj are wrong
      if(hm%ep%kick%time .eq. M_ZERO) then
        if(.not. cmplxscl) then
          call kick_apply(gr, st, td%ions, geo, hm%ep%kick)
        else
          call kick_apply(gr, st, td%ions, geo, hm%ep%kick, hm%cmplxscl%theta)
        end if
        call td_write_kick(gr, hm, sys%outp, geo, 0)
      end if
      call propagator_run_zero_iter(hm, gr, td%tr)
      call td_write_data(write_handler, gr, st, hm, sys%ks%xc, sys%outp, geo, 0)

      POP_SUB(td_run.td_run_zero_iter)
    end subroutine td_run_zero_iter


    ! ---------------------------------------------------------
    !> reads the pos and vel from coordinates file
    subroutine td_read_coordinates() 
      integer :: iatom, iter, iunit
      PUSH_SUB(td_run.td_read_coordinates)

      call io_assign(iunit)
      open(unit = iunit, file = io_workpath('td.general/coordinates'), action='read', status='old')

      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/coordinates'))//"'."
        message(2) = "Starting simulation from initial geometry."
        call messages_warning(2)
        POP_SUB(td_run.td_read_coordinates)
        return
      end if

      call io_skip_header(iunit)
      do iter = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no portable seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%x(1:gr%mesh%sb%dim)
        geo%atom(iatom)%x(:) = units_to_atomic(units_inp%length, geo%atom(iatom)%x(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%v(1:gr%mesh%sb%dim)
        geo%atom(iatom)%v(:) = units_to_atomic(units_inp%velocity, geo%atom(iatom)%v(:))
      end do
      do iatom = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(iatom)%f(1:gr%mesh%sb%dim)
        geo%atom(iatom)%f(:) = units_to_atomic(units_inp%force, geo%atom(iatom)%f(:))
      end do

      call io_close(iunit)
      POP_SUB(td_run.td_read_coordinates)
    end subroutine td_read_coordinates

    ! ---------------------------------------------------------
    subroutine td_read_gauge_field()

      integer :: iter, iunit
      FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM), dummy(1:MAX_DIM)

      PUSH_SUB(td_run.td_read_gauge_field)

      call io_assign(iunit)
      open(unit = iunit, file = io_workpath('td.general/gauge_field'), &
        action='read', status='old')
      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/gauge_field'))//"'."
        message(2) = "Starting simulation from initial values."
        call messages_warning(2)
        POP_SUB(td_run.td_read_gauge_field)
        return
      end if

      call io_skip_header(iunit)
      do iter = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

      vecpot = M_ZERO
      vecpot_vel = M_ZERO

      read(iunit, '(3es20.12)', advance='no') vecpot(1:gr%mesh%sb%dim)
      vecpot(1:gr%mesh%sb%dim) = units_to_atomic(units_inp%energy, vecpot(1:gr%mesh%sb%dim))
      ! A ~ B * length ~ E * length ~ force * length ~ energy (all for e = 1)
      read(iunit, '(3es20.12)', advance='no') vecpot_vel(1:gr%mesh%sb%dim)
      vecpot_vel(1:gr%mesh%sb%dim) = units_to_atomic(units_inp%energy / units_inp%time, vecpot_vel(1:gr%mesh%sb%dim))
      read(iunit, '(3es20.12)', advance='no') dummy(1:gr%mesh%sb%dim) ! skip the accel field.

      call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
      call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
      call hamiltonian_update(hm, gr%mesh, time = td%iter*td%dt)

      call io_close(iunit)
      POP_SUB(td_run.td_read_gauge_field)
    end subroutine td_read_gauge_field

    ! ---------------------------------------------------------
    subroutine td_save_restart(iter)
      integer, intent(in) :: iter

      integer :: ii, is, ierr
      character(len=256) :: filename
      CMPLX, allocatable :: zv_old(:)

      if(.not. write_restart()) return

      PUSH_SUB(td_run.td_save_restart)

      ! first write resume file
      call restart_write(trim(tmpdir)//'td', st, gr, geo, ierr, iter)
      if(ierr.ne.0) then
        message(1) = 'Unsuccessful write of "'//trim(tmpdir)//'td"'
        call messages_fatal(1)
      end if
      
      if(cmplxscl) then
        SAFE_ALLOCATE(zv_old(1:gr%mesh%np))
      end if

      ! write potential from previous interactions
      do ii = 1, 2
        do is = 1, st%d%nspin
          write(filename,'(a6,i2.2,i3.3)') 'vprev_', ii, is
          if(cmplxscl) then
            zv_old = td%tr%v_old(1:gr%mesh%np, is, ii) + M_zI * td%tr%Imv_old(1:gr%mesh%np, is, ii)
            call zio_function_output(restart_format, trim(tmpdir)//"td", &
              filename, gr%mesh, zv_old, unit_one, ierr, &
              is_tmp = .true., grp = st%dom_st_kpt_mpi_grp)
          else
            call dio_function_output(restart_format, trim(tmpdir)//"td", &
              filename, gr%mesh, td%tr%v_old(1:gr%mesh%np, is, ii), unit_one, ierr, &
              is_tmp = .true., grp = st%dom_st_kpt_mpi_grp)
          end if
          ! the unit is energy actually, but this only for restart, and can be kept in atomic units
          ! for simplicity
          if(ierr.ne.0) then
            write(message(1), '(3a)') 'Unsuccessful write of "', trim(filename), '"'
            call messages_fatal(1)
          end if
        end do
      end do

      if(td%dynamics == CP) call cpmd_restart_write(td%cp_propagator, gr, st)

      if(cmplxscl) then
        SAFE_DEALLOCATE_A(zv_old)
      end if
      

      POP_SUB(td_run.td_save_restart)
    end subroutine td_save_restart

  end subroutine td_run

#include "td_init_inc.F90"

end module td_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
