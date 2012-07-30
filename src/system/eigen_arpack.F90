!! Copyright (C) 2008 X. Andrade
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
!! $Id: eigen_rmmdiis.F90 5954 2009-10-17 20:53:52Z xavier $

#include "global.h"


module eigen_arpack_m

#if defined(HAVE_ARPACK) || defined(HAVE_PARPACK)   
use batch_m
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
  use mesh_batch_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use profiling_m
  use states_m
  use states_calc_m

  implicit none

  private
  public ::                     &
    eigen_arpack_t,             &
    deigen_solver_arpack,       &
    zeigen_solver_arpack

    type eigen_arpack_t
      integer          :: arnoldi_vectors !< number of Arnoldi vectors
      character(len=2) :: sort            !< which eigenvalue sorting 
      integer          :: init_resid      !< inital residual strategy 
    end type eigen_arpack_t

  contains
    
  subroutine arpack_debug(debug_level)
    integer, intent(in) :: debug_level

   
! Modified from ARPACK debug.h 
! I don't think this is going to change too much.. well at least it didn't since 1997 :)
!
!\SCCS Information: @(#) 
! FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
!
!     %---------------------------------%
!     | See debug.doc for documentation |
!     %---------------------------------%
    integer  logfil, ndigit, mgetv0, &
             msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
             mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
             mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
    common /debug/ &
             logfil, ndigit, mgetv0, &
             msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
             mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
             mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

    PUSH_SUB(arpack_debug)
    
    ndigit = -3
    logfil = 6
    mnaitr = 0
    mnapps = 0
    mnaupd = 3
    mnaup2 = 3
    mneigh = 0
    mneupd = 3
    
    POP_SUB(arpack_debug)
  end subroutine arpack_debug
    
  !----------------------------------------------------
  subroutine arpack_check_error(sub, info)
    integer,           intent(in) :: info
    character(len= *), intent(in) :: sub
    
    integer :: msg_lines
    logical :: OK
    
    PUSH_SUB(arpack_check_error)
    
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
      write(message(1),'(a,a,a,i5)') 'P/ARPACK ',sub, ', info = ', info
      call messages_warning(msg_lines)
    end if
    
    POP_SUB(arpack_check_error)
  end subroutine arpack_check_error


#include "real.F90" 
#include "eigen_arpack_inc.F90" 
#include "undef.F90" 

#include "complex.F90" 
#include "eigen_arpack_inc.F90" 
#include "undef.F90" 

#else 
! this avoids compilers complaining about empty module 
  integer, public :: arpack_dummy 
#endif 

  end module eigen_arpack_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
