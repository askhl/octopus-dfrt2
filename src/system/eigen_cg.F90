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
!! $Id: eigen_cg.F90 5954 2009-10-17 20:53:52Z xavier $

#include "global.h"

module eigen_cg_m
  use comm_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use lalg_adv_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use preconditioners_m
  use profiling_m
  use states_m
  use states_calc_m

  implicit none

  private
  public ::                 &
    deigensolver_cg2,       &
    zeigensolver_cg2,       &
    deigensolver_cg2_new,   &
    zeigensolver_cg2_new,   &
    eigensolver_bicg,       &
    eigensolver_direct
contains

  subroutine eigensolver_direct(gr, st, hm, pre, tol, niter, converged, ik, diff)
    type(grid_t),           intent(in)    :: gr
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(preconditioner_t), intent(in)    :: pre
    FLOAT,                  intent(in)    :: tol
    integer,                intent(inout) :: niter
    integer,                intent(inout) :: converged
    integer,                intent(in)    :: ik
    FLOAT,        optional, intent(out)   :: diff(1:st%nst)
        
    CMPLX, allocatable :: psi(:, :), h_psi(:,:), h_rr(:,:), cL_rr(:,:), cR_rr(:,:), zeigenval(:), manyzeigenval(:)
    FLOAT, allocatable :: eigenval(:), manyeigenval(:), sortkey(:)
    integer, allocatable :: sortindices(:)
    FLOAT :: spacingsquared
    CMPLX :: kinetic_phase, tmp, tmp2
    logical :: fdkinetic

    integer       :: ib, jb, p, errcode
    
    print*, 'hello world direct'

    SAFE_ALLOCATE(  psi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(h_psi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(h_rr(1:gr%mesh%np, 1:gr%mesh%np))
    SAFE_ALLOCATE(cL_rr(1:gr%mesh%np, 1:gr%mesh%np))
    SAFE_ALLOCATE(cR_rr(1:gr%mesh%np, 1:gr%mesh%np))
    SAFE_ALLOCATE(zeigenval(1:st%nst))
    SAFE_ALLOCATE(manyzeigenval(1:gr%mesh%np))
    SAFE_ALLOCATE(sortkey(1:gr%mesh%np))
    SAFE_ALLOCATE(sortindices(1:gr%mesh%np))
    h_psi=M_z0 !R_TOTYPE(M_ZERO) ! presumably this sets h_psi = 0
    h_rr=M_z0
    psi=M_z0

    kinetic_phase = exp(-2.0 * M_zI * hm%cmplxscl_th)
    
    spacingsquared = (gr%mesh%spacing(1) * gr%mesh%spacing(1))

    fdkinetic = .false.
    !fdkinetic = .true.

    if (.not.fdkinetic) then
       do ib = 1, gr%mesh%np
          do jb = 1, gr%mesh%np
             tmp = 0
             do p = 1, (gr%mesh%np - 1) / 2
                tmp2 = p * M_PI / (gr%mesh%np * gr%mesh%spacing(1))
                tmp = tmp + cos((p*2*M_PI*(ib-jb))/gr%mesh%np)*2*tmp2*tmp2
             end do
             h_rr(jb, ib) = tmp * 2 * kinetic_phase / gr%mesh%np
          end do
       end do
    end if

    do ib = 1, gr%mesh%np
       ! kinetic fd stencil
       if (fdkinetic) then
          h_rr(ib, ib) = (1.0, 0.0) / spacingsquared * kinetic_phase
          if (ib > 1) then
             h_rr(ib, ib - 1) = (-0.5, 0.0) / spacingsquared * kinetic_phase
             h_rr(ib - 1, ib) = (-0.5, 0.0) / spacingsquared * kinetic_phase
          end if
       end if
       h_rr(ib, ib) = h_rr(ib, ib) + hm%hm_base%potential(ib, 1) + (0.0, 1.0) * hm%hm_base%Impotential(ib, 1)
    end do

    cL_rr = h_rr
    cR_rr = h_rr
    
    call lalg_eigensolve_nonh(gr%mesh%np, cR_rr, manyzeigenval, errcode, 'R')
    if (errcode.ne.0) then
       print*, 'something went wrong, errcode', errcode
    end if
    call lalg_eigensolve_nonh(gr%mesh%np, cL_rr, manyzeigenval, errcode, 'L')
    if (errcode.ne.0) then
       print*, 'something went wrong, errcode', errcode
    end if
    
    !sortkey(:) = -imag(manyzeigenval(:))
    !sortkey(:) = real(manyzeigenval(:)) - imag(manyzeigenval(:))
    !sortkey(:) = real(manyzeigenval(:))
    sortkey(:) = abs(manyzeigenval(:))
    call sort(sortkey, sortindices)
    
    do p = 1, st%nst
       zeigenval(p) = manyzeigenval(sortindices(p))
    end do
    
    do p = 1, st%nst
       do ib = 1, gr%mesh%np
          st%psi%zL(ib, 1, p, 1) = cL_rr(ib, p)
          st%psi%zR(ib, 1, p, 1) = cR_rr(ib, p)
       end do
       st%zeigenval%Re(p, ik) = real(zeigenval(p))
       st%zeigenval%Im(p, ik) = aimag(zeigenval(p))
    end do

    if (present(diff)) diff = M_ZERO

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(h_psi)
    SAFE_DEALLOCATE_A(h_rr)
    SAFE_DEALLOCATE_A(zeigenval)
    SAFE_DEALLOCATE_A(manyeigenval)
    SAFE_DEALLOCATE_A(manyzeigenval)
    SAFE_DEALLOCATE_A(sortkey)
    SAFE_DEALLOCATE_A(sortindices)
    SAFE_DEALLOCATE_A(cL_rr)
    SAFE_DEALLOCATE_A(cR_rr)
  end subroutine eigensolver_direct


subroutine eigensolver_bicg (gr, st, hm, pre, tol, niter, converged, ik, diff)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  FLOAT,        optional, intent(out)   :: diff(1:st%nst)

  CMPLX, allocatable :: h_psi(:,:), g(:,:), g0(:,:),  cg(:,:), ppsi(:,:), psi(:, :)
  CMPLX   :: es(2), a0, b0, gg, gg0, gg1, gamma, theta, norma, cg0
  !real(8)  :: cg0, e0, res
  real(8)  :: res
  integer  :: p, iter, maxter, idim, ip
  CMPLX   :: sb(3), tmpzeigenval, e0

  PUSH_SUB(eigensolver_bicg)

  maxter = niter
  niter = 0

  SAFE_ALLOCATE(  psi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(h_psi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(   cg(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(    g(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(   g0(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( ppsi(1:gr%mesh%np, 1:st%d%dim))
  h_psi = M_z0
  cg    = M_z0
  g     = M_z0
  g0    = M_z0
  ppsi  = M_z0

  do idim = 1, st%d%dim
    cg(1:gr%mesh%np_part, idim) = M_z0
  end do

  ! Set the diff to zero, since it is intent(out)
  if(present(diff)) diff(1:st%nst) = M_z0

  ! Start of main loop, which runs over all the eigenvectors searched
  ASSERT(converged >= 0)

  eigenfunction_loop : do p = converged + 1, st%nst

    call states_get_state(st, gr%mesh, p, ik, psi)

    ! Orthogonalize starting eigenfunctions to those already calculated...
    if(p > 1) call zstates_orthogonalize_single(st, gr%mesh, p - 1, ik, psi, normalize = .true.)

    ! Calculate starting gradient: |hpsi> = H|psi>
    call zhamiltonian_apply(hm, gr%der, psi, h_psi, p, ik)

    ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
    tmpzeigenval = zmf_dotp (gr%mesh, st%d%dim, psi, h_psi, dotu = .true.)
    st%zeigenval%Re(p, ik) = real(tmpzeigenval)
    st%zeigenval%Im(p, ik) = aimag(tmpzeigenval)
    !st%eigenval(p, ik) = REAL(tmpzeigenval)

    ! Starts iteration for this band
    iter_loop: do iter = 1, maxter

      ! inverse preconditioner....
      call  zpreconditioner_apply(pre, gr, hm, ik, h_psi, g)
      call  zpreconditioner_apply(pre, gr, hm, ik, psi, ppsi)

      es(1) = zmf_dotp (gr%mesh, st%d%dim, psi, g, reduce = .false., dotu = .true.)
      es(2) = zmf_dotp (gr%mesh, st%d%dim, psi, ppsi, reduce = .false., dotu = .true.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%vp%comm, es, dim = 2)

      es(1) = es(1)/es(2)

      do idim = 1, st%d%dim
         ! the REAL below was R_TOPREC
         ! apparently the only effect of the cmlx(real/aimag) is to reconstruct
         ! the number with a specific precision
        !call lalg_axpy(gr%mesh%np, cmplx(real(-es(1)), aimag(-es(1)), REAL_PRECISION), ppsi(:, idim), g(:, idim))
         call lalg_axpy(gr%mesh%np, -es(1), ppsi(:, idim), g(:, idim))
      end do

      ! Orthogonalize to lowest eigenvalues (already calculated)
      if(p > 1) call zstates_orthogonalize_single(st, gr%mesh, p - 1, ik, g, normalize = .false.)

      if(iter .ne. 1) then
        gg1 = zmf_dotp (gr%mesh, st%d%dim, g, g0, reduce = .false., dotu = .true.)
      else
        gg1 = M_ZERO
      end if

      ! Approximate inverse preconditioner...
      call  zpreconditioner_apply(pre, gr, hm, ik, g(:,:), g0(:,:))

      gg = zmf_dotp (gr%mesh, st%d%dim, g, g0, reduce = .false., dotu = .true.)

      if(gr%mesh%parallel_in_domains) then
        sb(1) = gg1
        sb(2) = gg
        call comm_allreduce(gr%mesh%vp%comm, sb, dim = 2)
        gg1 = sb(1)
        gg  = sb(2)
      end if

      if( abs(gg) < M_EPSILON ) then
        if(converged == p - 1) converged = p ! only consider the first converged eigenvectors
        !st%eigenval(p, ik) = es(1)
        st%zeigenval%Re(p, ik) = real(es(1))
        st%zeigenval%Im(p, ik) = aimag(es(1))
        res = sqrt(abs(gg)) ! Why???  res is assigned later.  Maybe not always?
        exit
      end if

      ! Starting or following iterations...
      if(iter .eq. 1) then
        gg0 = gg

        do idim = 1, st%d%dim
          call lalg_copy(gr%mesh%np, g(:,idim), cg(:, idim))
        end do
      else
        !gamma = gg/gg0        ! (Fletcher-Reeves)
        gamma = (gg - gg1)/gg0   ! (Polack-Ribiere)
        gg0 = gg
        
        norma = gamma*cg0*sin(theta)
        
        forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
          cg(ip, idim) = gamma*cg(ip, idim) + g(ip, idim) - norma*psi(ip, idim)
        end forall
        ! R_ADD and R_MUL
        call profiling_count_operations(st%d%dim*gr%mesh%np*(2*2 + 2*6))

      end if

      ! cg contains now the conjugate gradient
      call zhamiltonian_apply(hm, gr%der, cg, ppsi, p, ik)

      ! Line minimization.
      a0 = zmf_dotp (gr%mesh, st%d%dim, psi, ppsi, reduce = .false., dotu = .true.)
      b0 = zmf_dotp (gr%mesh, st%d%dim, cg, ppsi, reduce = .false., dotu = .true.)
      ! XXXXXXXXX choose the right root
      cg0 = sqrt(zmf_dotp (gr%mesh, st%d%dim, cg, cg, reduce = .false., dotu = .true.))

      if(gr%mesh%parallel_in_domains) then
        sb(1) = a0
        sb(2) = b0
        sb(3) = cg0**2
        call comm_allreduce(gr%mesh%vp%comm, sb, dim = 3)
        a0 = sb(1)
        b0 = sb(2)
        cg0 = sqrt(sb(3))
      end if

      a0 = M_TWO * a0 / cg0
      b0 = b0/cg0**2
      e0 = st%zeigenval%Re(p, ik) + M_zI * st%zeigenval%Im(p, ik)

      ! what the heck should we do here?  Probably the angle should be
      ! complex or something
      theta = atan(REAL(a0/(e0 - b0)))/M_TWO
      es(1) = M_HALF*((e0-b0)*cos(M_TWO*theta) + a0*sin(M_TWO*theta) + e0 + b0)
      es(2) = -M_HALF*((e0-b0)*cos(M_TWO*theta) + a0*sin(M_TWO*theta) - (e0 + b0))

      ! Choose the minimum solutions.
      ! What about imaginary part?
      if (REAL(es(2)) < REAL(es(1))) then
         theta = theta + M_PI/M_TWO
         st%zeigenval%Re(p, ik) = real(es(2))
         st%zeigenval%Im(p, ik) = aimag(es(2))
      else
         st%zeigenval%Re(p, ik) = real(es(1))
         st%zeigenval%Im(p, ik) = aimag(es(1))
      end if
      !st%eigenval(p, ik) = min(REAL(es(1)), REAL(es(2)))

      ! Upgrade psi...
      a0 = cos(theta)
      b0 = sin(theta)/cg0

      forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
        psi(ip, idim) = a0*psi(ip, idim) + b0*cg(ip, idim)
        h_psi(ip, idim) = a0*h_psi(ip, idim) + b0*ppsi(ip, idim)
      end forall
      
      ! R_ADD and R_MUL
      call profiling_count_operations(st%d%dim*gr%mesh%np*(2*2 + 4*6))

      !res = zstates_residue(gr%mesh, st%d%dim, h_psi, st%eigenval(p, ik), psi)

      if(in_debug_mode) then
        write(message(1), '(a,i4,a,i4,a,i4,a,f12.6)') 'Debug: BiCG Eigensolver - ik', ik, ' ist ', p, ' iter ', iter!, ' res ', res
        call messages_info(1)
      end if

      ! Test convergence.
      !if(res < tol) then
      if(iter > 190) then
        if(converged == p - 1) converged = p ! only consider the first converged eigenvectors
        exit iter_loop
      end if

    end do iter_loop

    do ip = 1, gr%mesh%np_part
       st%psi%zL(ip, 1, p, 1) = psi(ip, 1)
       st%psi%zR(ip, 1, p, 1) = psi(ip, 1)
    end do
    
    call states_set_state(st, gr%mesh, p, ik, psi)
    
    niter = niter + iter + 1

    if(present(diff)) then
      diff(p) = 0.0 !res
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik - 1) +  p, st%nst*st%d%nik)
    end if

  end do eigenfunction_loop

  print*, "eigs"
  do p = 1, st%nst
     print*, st%eigenval(p, ik)
  end do

  ! Deallocation of variables
  SAFE_DEALLOCATE_A(h_psi)
  SAFE_DEALLOCATE_A(g)
  SAFE_DEALLOCATE_A(g0)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(ppsi)

  POP_SUB(eigensolver_bicg)
end subroutine eigensolver_bicg
 



#include "real.F90"
#include "eigen_cg_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "eigen_cg_inc.F90"
#include "undef.F90"

end module eigen_cg_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
