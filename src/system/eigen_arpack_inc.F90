!!!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
	
	
subroutine X(eigen_solver_arpack)(gr, st, hm, tol_, niter, ncv, converged, ik, diff)
  type(grid_t),        intent(in)    :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  FLOAT,               intent(in)    :: tol_
  integer,             intent(inout) :: niter
  integer,             intent(in)    :: ncv
  integer,             intent(inout) :: converged
  integer,             intent(in)    :: ik
  FLOAT,     optional, intent(out)   :: diff(1:st%nst)
	
  logical, allocatable :: select(:)
  FLOAT, allocatable :: ax(:), d(:, :), resid(:), v(:, :),   &
    workd(:), workev(:), workl(:)
  integer :: ldv, nev, iparam(11), ipntr(14), ido, n, lworkl, info, ierr, &
    i, j, ishfts, maxitr, mode1
  FLOAT :: tol, sigmar, sigmai
	
	!!!!WARNING: No support for spinors, yet. No support for complex wavefunctions.
  PUSH_SUB('eigen_arpack.eigen_solver_arpack')
	

	
#if defined(HAVE_MPI)
  if(gr%m%parallel_in_domains) then
    message(1) = 'Error: Arpack-Solver not parallelized for domain decomposition.'
    call messages_fatal(1)
    !  FIXME: Need to adjust m%x and m%vol_pp occurences in the code below
    !         appropriately for domain decomposition. Also parallelization
    !         of the vectors has to be taken care of.
  end if
#endif
	
  ldv = gr%mesh%np
  n = gr%mesh%np
  nev = st%nst
  lworkl  = 3*ncv**2+6*ncv
  SAFE_ALLOCATE(ax(ldv))
  SAFE_ALLOCATE(d(ncv, 3))
  SAFE_ALLOCATE(resid(ldv))
  SAFE_ALLOCATE(v(ldv, ncv))
  SAFE_ALLOCATE(workd(3*ldv))
  SAFE_ALLOCATE(workev(3*ncv))
  SAFE_ALLOCATE(workl(lworkl))
  SAFE_ALLOCATE(select(ncv))
	
  select = .true.
  tol    = tol_
  ido    = 0
  info = 1
	
  do i = 1, gr%mesh%np
    resid(i) = sum(st%X(psi)(i, 1, 1:st%nst, ik))*sqrt(gr%mesh%vol_pp(i))
  end do
	
  ishfts = 1
  maxitr = niter
  mode1 = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
	
  do
    call dnaupd ( ido, 'I', n, 'SR', nev, tol, resid, ncv, &
      v, ldv, iparam, ipntr, workd, workl, lworkl, info )
    if( abs(ido).ne.1) exit
    call av (n, workd(ipntr(1)), workd(ipntr(2)))
  end do
  ! If info is larger than zero, it may not be an error (i.e., not all eigenvectors
  ! were converged)
  if(info .lt. 0) then
    write(message(1),'(a,i5)') 'Error with ARPACK _naupd, info = ', info
    write(message(2),'(a)')    'Check the documentation of _naupd.'
    call messages_fatal(2)
  end if
	
  call dneupd ( .true., 'A', select, d, d(1,2), v, ldv, &
    sigmar, sigmai, workev, 'I', n, 'SR', nev, tol, &
    resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
  lworkl, ierr )

  if(ierr .ne. 0) then
    write(message(1),'(a,i5)') 'Error with ARPACK _neupd, info = ', info
    write(message(2),'(a)')    'Check the documentation of _neupd.'
    call messages_fatal(2)
  end if

  ! This sets the number of converged eigenvectors.
  converged =  iparam(5)
  !call dmout(6, converged, 3, d, ncv, -6, 'Ritz values (Real, Imag) and residual residuals')
  ! This sets niter to the number of matrix-vector operations.
  niter = iparam(9)
  do j = 1, min(st%nst, converged)
    do i = 1, gr%mesh%np
      st%X(psi)(i, 1, j, ik) = v(i, j)/sqrt(gr%mesh%vol_pp(i))
    end do
    st%eigenval(j, ik) = d(j, 1)
    if(workl(ipntr(11)+j-1)< M_EPSILON) then
      diff(j) = M_ZERO
    else
      diff(j) = workl(ipntr(11)+j-1)
    end if
  end do

  SAFE_DEALLOCATE_A(ax)
  SAFE_DEALLOCATE_A(d)
  SAFE_DEALLOCATE_A(resid)
  SAFE_DEALLOCATE_A(v)
  SAFE_DEALLOCATE_A(workd)
  SAFE_DEALLOCATE_A(workev)
  SAFE_DEALLOCATE_A(workl)
  SAFE_DEALLOCATE_A(select)



   POP_SUB('eigen_arpack.eigen_solver_arpack')
contains

  ! ---------------------------------------------------------
  subroutine av (n, v, w)
    integer :: n
    FLOAT :: v(n), w(n)
    integer :: i, NP, NP_PART
    R_TYPE, allocatable :: psi(:, :), hpsi(:, :)
    
    NP = gr%mesh%np
    NP = gr%mesh%np_part

    SAFE_ALLOCATE(psi(NP_PART, 1))
    SAFE_ALLOCATE(hpsi(NP_PART, 1))

    do i = 1, NP
      psi(i, 1) = v(i)/sqrt(gr%mesh%vol_pp(i))
    end do
    do i = NP+1, NP_PART
      psi(i, 1) = M_ZERO
    end do
    
    call X(hamiltonian_apply) (hm, gr%der, psi, hpsi, st%nst, ik)
    
        
    do i = 1, NP
      w(i) = hpsi(i, 1)*sqrt(gr%mesh%vol_pp(i))
    end do


    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)

  end subroutine av

end subroutine X(eigen_solver_arpack)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: