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
    CMPLX :: cmplxscl_phase, cmplxscl_phase2

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

    cmplxscl_phase = exp((0.0, 1.0) * hm%cmplxscl_th)
    cmplxscl_phase2 = exp((0.0, -2.0) * hm%cmplxscl_th)
    
    spacingsquared = (gr%mesh%spacing(1) * gr%mesh%spacing(1))

    do ib = 1, gr%mesh%np
       h_rr(ib, ib) = (1.0, 0.0) / spacingsquared * cmplxscl_phase2
       h_rr(ib, ib) = h_rr(ib, ib) + hm%hm_base%potential(ib, 1) + (0.0, 1.0) * hm%hm_base%Impotential(ib, 1)
        if (ib > 1) then
          h_rr(ib, ib - 1) = (-0.5, 0.0) / spacingsquared * cmplxscl_phase2
          h_rr(ib - 1, ib) = (-0.5, 0.0) / spacingsquared * cmplxscl_phase2
       end if
    end do

    cL_rr = h_rr
    cR_rr = h_rr
    
    call lalg_eigensolve_nonh(gr%mesh%np, h_rr, manyzeigenval, errcode, 'R')
    if (errcode.ne.0) then
       print*, 'something went wrong, errcode', errcode
    end if
!     call lalg_eigensolve_nonh(gr%mesh%np, cL_rr, manyzeigenval, errcode, 'L')
!     if (errcode.ne.0) then
!        print*, 'something went wrong, errcode', errcode
!     end if


    
!     sortkey(:) = real(manyzeigenval(:))
!      sortkey(:) = -aimag(manyzeigenval(:))
    sortkey(:) = abs(manyzeigenval(:))
!     sortkey(:) = real(manyzeigenval(:)) -aimag(manyzeigenval(:)) 
    call sort(sortkey, sortindices)
     
     do p = 1, gr%mesh%np
       write (*,*) p, manyzeigenval(sortindices(p))
     end do
    
    do p = 1, st%nst
       zeigenval(p) = manyzeigenval(sortindices(p))
    end do
    
!     sortkey1(:) = real(zeigenval(:))
!     call sort(sortkey, sortindices)
!     do p = 1, st%nst
!        zeigenval(p) = zeigenval(sortindices(p))
!     end do


    do p = 1, st%nst
       do ib = 1, gr%mesh%np
!           st%psi%zL(ib, 1, p, 1) = cL_rr(ib, p)
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
  

  subroutine eigensolver_bicg(gr, st, hm, pre, tol, niter, converged, ik, diff)
    type(grid_t),           intent(in)    :: gr
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(preconditioner_t), intent(in)    :: pre
    FLOAT,                  intent(in)    :: tol
    integer,                intent(inout) :: niter
    integer,                intent(inout) :: converged
    integer,                intent(in)    :: ik
    FLOAT,        optional, intent(out)   :: diff(1:st%nst)

    print*, 'hello world bicg'
    
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
