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
  R_TYPE, allocatable  :: resid(:), v(:, :),   &
                          workd(:), workev(:), workl(:), zd(:), &
                          psi(:,:)
                     
  integer :: ldv, nev, iparam(11), ipntr(14), ido, n, lworkl, info, ierr, &
             i, j, ishfts, maxitr, mode1, ist
  FLOAT :: tol, sigmar, sigmai
  FLOAT, allocatable :: rwork(:), d(:, :) 
  CMPLX :: sigma 
  integer :: mpi_comm 
  	
	!!!!WARNING: No support for spinors, yet. 
 
  PUSH_SUB(eigen_arpack.eigen_solver_arpack)
	


  if(st%parallel_in_states) then
    message(1) = 'Arpack-Solver not parallelized for states decomposition.'
    message(2) = 'Change ParallelizationStrategy and rerun.'
    call messages_fatal(2)
  end if

	mpi_comm = mpi_world%comm
  if (gr%mesh%parallel_in_domains) mpi_comm = gr%mesh%mpi_grp%comm
  

  n = gr%mesh%np
!   n = gr%mesh%np_part
  ldv = n
  nev = st%nst
  lworkl  = 3*ncv**2+6*ncv

  SAFE_ALLOCATE(d(ncv+1, 3))
  SAFE_ALLOCATE(resid(ldv))       !residual vector 
  SAFE_ALLOCATE(v(ldv, ncv))      !Arnoldi basis vectors / Eigenstates
  SAFE_ALLOCATE(workd(3*ldv))
  SAFE_ALLOCATE(workev(3*ncv))
  SAFE_ALLOCATE(workl(lworkl))
  SAFE_ALLOCATE(select(ncv))
  
  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
  
#if defined(R_TCOMPLEX)
  SAFE_ALLOCATE(rwork(ncv))
  SAFE_ALLOCATE(zd(ncv+1))
#endif
	
  select(:) = .true.
  tol  = M_ZERO !tol_
  ido  = 0
!  info = 1 !resid contains the initial residual vector
  info = 0 
  
  print *,mpi_world%rank,  "tol", tol
  print *,mpi_world%rank, "Ncv", ncv, "nev", nev, "n", n

  
  do i = 1, ldv
!      resid(i) = sum(st%X(psi)(i, 1, 1:st%nst, ik))*sqrt(gr%mesh%vol_pp(1))
      resid(i) = R_TOTYPE(M_ONE)
  end do

	
  ishfts = 1
  maxitr = niter
  mode1 = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
	
  do
#if defined(R_TCOMPLEX)
 #if defined(HAVE_PARPACK) && defined(HAVE_MPI)

    call pznaupd  ( mpi_comm, &
          ido, 'I', n, 'SR', nev, tol, resid, ncv, &
          v, ldv, iparam, ipntr, workd, workl, lworkl, &
          rwork,info )

 #else
    call znaupd  ( & 
          ido, 'I', n, 'SR', nev, tol, resid, ncv, &
          v, ldv, iparam, ipntr, workd, workl, lworkl, &
          rwork,info )
 #endif

#else 
  #if defined(HAVE_PARPACK) && defined(HAVE_MPI)
    call pdnaupd  ( mpi_comm, &
          ido, 'I', n, 'SR', nev, tol, resid, ncv, &
          v, ldv, iparam, ipntr, workd, workl, lworkl, & 
          info )
  
  #else
    call dnaupd  ( & 
          ido, 'I', n, 'SR', nev, tol, resid, ncv, &
          v, ldv, iparam, ipntr, workd, workl, lworkl, & 
          info )
  #endif
#endif      
      
    if( abs(ido).ne.1) exit
    
    call av (n, workd(ipntr(1)), workd(ipntr(2))) ! calculate H * psi
    
  end do
  
  call arpack_check_error('naupd', info)
  

#if defined(R_TCOMPLEX) 
 #if defined(HAVE_PARPACK) && defined(HAVE_MPI)
  call pzneupd  (mpi_comm, .true., 'A', select, zd, v, ldv, sigma, &
        workev, 'I', n, 'SR', nev, tol, resid, ncv, & 
        v, ldv, iparam, ipntr, workd, workl, lworkl, &
        rwork, info)
        d(:,1)=real(zd(:))
        d(:,2)=aimag(zd(:))
        d(:,3)=M_ZERO
 #else
  call zneupd  (.true., 'A', select, zd, v, ldv, sigma, &
        workev, 'I', n, 'SR', nev, tol, resid, ncv, & 
        v, ldv, iparam, ipntr, workd, workl, lworkl, &
        rwork, info)
        d(:,1)=real(zd(:))
        d(:,2)=aimag(zd(:))
        d(:,3)=M_ZERO
 #endif
        
#else	
 #if defined(HAVE_PARPACK) && defined(HAVE_MPI)
  call pdneupd (mpi_comm, .true., 'A', select, d, d(1,2), v, ldv, &
       sigmar, sigmai, workev, 'I', n, 'SR', nev, tol, &
       resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
       lworkl, info )
 
 #else
  call dneupd ( .true., 'A', select, d, d(1,2), v, ldv, &
       sigmar, sigmai, workev, 'I', n, 'SR', nev, tol, &
       resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
       lworkl, info )
 #endif
#endif
       
  call arpack_check_error('neupd', info)


  ! This sets the number of converged eigenvectors.
  converged =  iparam(5)

!   call dmout(6, converged, 3, d, ncv+1, -6, 'Ritz values (Real, Imag) and residual residuals')

  ! This sets niter to the number of matrix-vector operations.
  niter = iparam(9)
  do j = 1, converged
    do i = 1, gr%mesh%np
      psi(i,1) = v(i, j)/sqrt(gr%mesh%volume_element) 
    end do
    do i = gr%mesh%np + 1, gr%mesh%np_part
      psi(i,1) = R_TOTYPE(M_ZERO) 
    end do
    
    call states_set_state(st, gr%mesh, j, ik, psi)
    
!     print *,"st", j, "norm", sqrt(X(mf_dotp)(gr%mesh, st%d%dim, psi, psi, dotu = .true.))

        
    st%eigenval(j, ik) = d(j, 1)
    if(associated(st%zeigenval%Im))then 
      st%zeigenval%Im(j, ik) = d(j, 2)
    end if

    if(abs(workl(ipntr(11)+j-1))< M_EPSILON) then
      diff(j) = M_ZERO
    else
      diff(j) = workl(ipntr(11)+j-1)
    end if
  end do

  !Fill unconverged states
  do j = converged + 1, st%nst
    do i = 1, gr%mesh%np
!       st%X(psi)(i, 1, j, ik) = R_TOTYPE(M_ONE)
      psi(i,1) = R_TOTYPE(M_ONE) 
    end do
    call states_set_state(st, gr%mesh, j, ik, psi)

    st%eigenval(j, ik) = M_HUGE
    if(associated(st%zeigenval%Im))then 
      st%zeigenval%Im(j, ik) = M_HUGE
    end if
    diff(j) = M_HUGE
  end do



  SAFE_DEALLOCATE_A(d)
  SAFE_DEALLOCATE_A(resid)
  SAFE_DEALLOCATE_A(v)
  SAFE_DEALLOCATE_A(workd)
  SAFE_DEALLOCATE_A(workev)
  SAFE_DEALLOCATE_A(workl)
  SAFE_DEALLOCATE_A(select)
  
  SAFE_DEALLOCATE_A(psi)
  
#if defined(R_TCOMPLEX)
  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(zd)  
#endif



   POP_SUB(eigen_arpack.eigen_solver_arpack)
contains

  ! ---------------------------------------------------------
  subroutine av (n, v, w)
    integer, intent(in) :: n
    R_TYPE,  intent(in) :: v(n)
    R_TYPE,  intent(out):: w(n)
    
    integer :: i, NP, NP_PART
    R_TYPE, allocatable :: psi(:, :), hpsi(:, :)
    
   
    
    NP = gr%mesh%np
    NP_PART = gr%mesh%np_part

    ASSERT(n == NP .or. n == NP_PART)

    SAFE_ALLOCATE(psi(NP_PART, hm%d%dim))
    SAFE_ALLOCATE(hpsi(NP_PART, hm%d%dim))

    do i = 1, NP
      psi(i, 1) = v(i)/sqrt(gr%mesh%volume_element)
    end do
    do i = NP+1, NP_PART
      psi(i, 1) = M_ZERO
    end do
    
    call X(hamiltonian_apply) (hm, gr%der, psi, hpsi, 1, ik)
    
    do i = 1, NP
      w(i) = hpsi(i, 1)*sqrt(gr%mesh%volume_element)
    end do


    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)

  end subroutine av
  
  !----------------------------------------------------
  subroutine arpack_check_error(sub, info)
    integer,           intent(in) :: info
    character(len= *), intent(in) :: sub
    
    integer :: msg_lines
    logical :: OK
    
    msg_lines = 1
    OK = .false.
  
    if (sub == 'neupd') then
      
      select case (info)
        case (0)
          OK = .true.
          
        case (1)
          write(message(2),'(a)'), 'The Schur form computed by LAPACK routine csheqr'
          write(message(3),'(a)'), 'could not be reordered by LAPACK routine ztrsen.'
          write(message(4),'(a)'), 'Re-enter subroutine pzneupd with IPARAM(5)=NCV and'
          write(message(5),'(a)'), 'increase the size of the array D to have'
          write(message(6),'(a)'), 'dimension at least dimension NCV and allocate at least NCV'
          write(message(7),'(a)'), 'columns for Z. NOTE: Not necessary if Z and V share'
          write(message(8),'(a)'), 'the same space. Please notify the authors if this error'
          write(message(9),'(a)'), 'occurs.'
          msg_lines = 9
         
        case(-1) 
          write(message(2),'(a)'), 'N must be positive.'
          msg_lines = 2
                         
        case(-2)
          write(message(2),'(a)'), 'NEV must be positive.'
          msg_lines = 2
          
        case(-3)
          write(message(2),'(a)'), 'NCV-NEV >= 2 and less than or equal to N.'
          msg_lines = 2
          
        case(-5)
          write(message(2),'(a)'), 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI"'
          msg_lines = 2
          
        case(-6)
          write(message(2),'(a)'), 'BMAT must be one of "I" or "G".'
          msg_lines = 2
          
        case(-7)
          write(message(2),'(a)'), 'Length of private work WORKL array is not sufficient.'
          msg_lines = 2
          
        case(-8)
          write(message(2),'(a)'), 'Error return from LAPACK eigenvalue calculation.'
          write(message(3),'(a)'), 'This should never happened.'
          msg_lines = 3
          
        case(-9)
          write(message(2),'(a)'), 'Error return from calculation of eigenvectors.'
          write(message(3),'(a)'), 'Informational error from LAPACK routine ztrevc.'
          msg_lines = 3
                               
        case(-10)
          write(message(2),'(a)'), 'IPARAM(7) must be 1,2,3'
          msg_lines = 2
                        
        case(-11)
          write(message(2),'(a)'), 'PARAM(7) = 1 and BMAT = "G" are incompatible.'
          msg_lines = 2       
                            
        case(-12)
          write(message(2),'(a)'), 'HOWMNY = "S" not yet implemented'
          msg_lines = 2
          
        case(-13)
          write(message(2),'(a)'), 'OWMNY must be one of "A" or "P" if RVEC = .true.'
          msg_lines = 2
          
        case(-14)
          write(message(2),'(a)'), 'PZNAUPD did not find any eigenvalues to sufficient'
          write(message(3),'(a)'), 'accuracy.'
          msg_lines = 3
          
        case(-15)
          write(message(2),'(a)'), 'ZNEUPD got a different count of the number of converged'
          write(message(3),'(a)'), 'Ritz values than ZNAUPD got.  This indicates the user'
          write(message(4),'(a)'), 'probably made an error in passing data from ZNAUPD to'
          write(message(5),'(a)'), 'ZNEUPD or that the data was modified before entering'
          write(message(6),'(a)'), 'ZNEUPD.'
          msg_lines = 6          
          
      end select
      
    else if( sub == 'naupd') then
  
      select case (info)
        case (0)
          OK = .true.
          
        case (1)
          write(message(2),'(a)'), 'Maximum number of iterations taken.'
          write(message(3),'(a)'), 'All possible eigenvalues of OP has been found. IPARAM(5)'
          write(message(4),'(a)'), 'returns the number of wanted converged Ritz values.'
          msg_lines = 4
          OK = .true.
          
        case (2)        
          write(message(2),'(a)'), 'No longer an informational error. Deprecated starting'
          write(message(3),'(a)'), 'with release 2 of ARPACK.'
          msg_lines = 3
          OK = .true.
          
        case (3)
          write(message(2),'(a)'), 'No shifts could be applied during a cycle of the'
          write(message(3),'(a)'), 'Implicitly restarted Arnoldi iteration. One possibility'
          write(message(4),'(a)'), 'is to increase the size of NCV relative to NEV.'
          write(message(5),'(a)'), 'See remark 4 below.'
          msg_lines = 5
                
        case (-1)
           write(message(2),'(a)'), 'N must be positive.'
           msg_lines = 2
           
        case (-2)       
           write(message(2),'(a)'), 'NEV must be positive.'
           msg_lines = 2
           
        case (-3)
           write(message(2),'(a)'), 'NCV-NEV >= 2 and less than or equal to N.'
           msg_lines = 2
           
        case (-4)
           write(message(2),'(a)'), 'The maximum number of Arnoldi update iteration'          
           write(message(3),'(a)'), 'must be greater than zero.'
           msg_lines = 3
                
        case (-5)
           write(message(2),'(a)'), 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI"'
           msg_lines = 2
           
        case (-6)
           write(message(2),'(a)'), 'BMAT must be one of "I" or "G".'
           msg_lines = 2
           
        case (-7)
           write(message(2),'(a)'), 'Length of private work array is not sufficient.'
           msg_lines = 2
           
        case (-8)
           write(message(2),'(a)'), 'Error return from LAPACK eigenvalue calculation;'
           msg_lines = 2
           
        case (-9)
           write(message(2),'(a)'), 'Starting vector is zero.'
           msg_lines = 2
           
        case (-10)
           write(message(2),'(a)'), 'IPARAM(7) must be 1,2,3.'
           msg_lines = 2
           
        case (-11)
           write(message(2),'(a)'), 'IPARAM(7) = 1 and BMAT = "G" are incompatable.'
           msg_lines = 2
           
        case (-12)
           write(message(2),'(a)'), 'IPARAM(1) must be equal to 0 or 1.'
           msg_lines = 2
           
        case (-9999)
           write(message(2),'(a)'), 'Could not build an Arnoldi factorization.'
           write(message(3),'(a)'), 'User input error highly likely.  Please'
           write(message(4),'(a)'), 'check actual array dimensions and layout.'
           write(message(5),'(a)'), 'IPARAM(5) returns the size of the current Arnoldi'
           write(message(6),'(a)'), 'factorization.'
           msg_lines = 6
          
      end select
        
        
    else
      write(message(1),'(a)') 'Unrecognized arpack subroutine '
      call messages_fatal(1)  
      
    end if
 
    if(.not. OK) then
      write(message(1),'(a,a,a,i5)') 'Error with P/ARPACK ', sub, ', info = ', info
      call messages_fatal(msg_lines)
    else if(msg_lines >= 2) then      
      write(message(1),'(a)') 'P/ARPACK eigensolver:'
      call messages_warning(msg_lines)
    end if
    
  end subroutine arpack_check_error
 
        
         

end subroutine X(eigen_solver_arpack)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
