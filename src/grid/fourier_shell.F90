!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: fourier_shell.F90 6722 2010-06-13 12:44:43Z acastro $

#include "global.h"

module fourier_shell_m
  use cube_function_m
  use cube_m
  use fft_m
  use global_m
  use math_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m

  implicit none

  private
  public ::                           &
       fourier_shell_t,               &
       fourier_shell_init,            &
       fourier_shell_end

  type fourier_shell_t
    integer          :: ngvectors
    FLOAT            :: ekin_cutoff
    integer, pointer :: coords(:, :)
    integer, pointer :: red_gvec(:, :)   
  end type fourier_shell_t

contains

  !> convention = .true.  : G-vectors are in (-n/2, n/2]
  !> convention = .false. : G-vectors are in [-n/2, n/2)
  subroutine fourier_shell_init(this, cube, mesh, convention)
    type(fourier_shell_t), intent(out)   :: this
    type(cube_t),          intent(in)    :: cube
    type(mesh_t),          intent(in)    :: mesh
    logical,               intent(in)    :: convention

    integer :: ig, ix, iy, iz, ixx(1:3), imap
    FLOAT :: dg(1:3), gmax2, gvec(1:3)
    FLOAT, allocatable :: modg2(:)
    integer, allocatable :: map(:), ucoords(:, :), ured_gvec(:, :)

    PUSH_SUB(fourier_shell_init)

    dg(1:3) = M_PI/(cube%rs_n_global(1:3)/2*mesh%spacing(1:3))
    gmax2 = (dg(1)*(cube%rs_n_global(1)/2))**2
    this%ekin_cutoff = gmax2/M_TWO

    SAFE_ALLOCATE(modg2(1:product(cube%rs_n_global(1:3))))
    SAFE_ALLOCATE(ucoords(1:3, 1:product(cube%rs_n_global(1:3))))
    SAFE_ALLOCATE(ured_gvec(1:3, 1:product(cube%rs_n_global(1:3))))
    
    ig = 0
    do ix = 1, cube%rs_n_global(1)
      ixx(1) = pad_feq(ix, cube%rs_n_global(1), .true.)
      if(.not. convention .and. 2 * ixx(1) == cube%rs_n_global(1)) ixx(1) = -ixx(1)
      do iy = 1, cube%rs_n_global(2)
        ixx(2) = pad_feq(iy, cube%rs_n_global(2), .true.)
        if(.not. convention .and. 2 * ixx(2) == cube%rs_n_global(2)) ixx(2) = -ixx(2)
        do iz = 1, cube%rs_n_global(3)
          ixx(3) = pad_feq(iz, cube%rs_n_global(3), .true.)
          if(.not. convention .and. 2 * ixx(3) == cube%rs_n_global(3)) ixx(3) = -ixx(3)

          gvec(1:3) = dg(1:3)*ixx(1:3)
          if(sum(gvec(1:3)**2) <= gmax2 + CNST(1e-10)) then
            INCR(ig, 1)
            ucoords(1:3, ig) = (/ ix, iy, iz /)
            ured_gvec(1:3, ig) = ixx(1:3)
            modg2(ig) = sum(gvec(1:3)**2)
          end if

        end do
      end do
    end do

    this%ngvectors = ig

    SAFE_ALLOCATE(this%coords(1:3, 1:this%ngvectors))
    SAFE_ALLOCATE(this%red_gvec(1:3, 1:this%ngvectors))
    SAFE_ALLOCATE(map(1:this%ngvectors))

    do ig = 1, this%ngvectors
      map(ig) = ig
    end do
    
    call sort(modg2(1:this%ngvectors), map)

    do ig = 1, this%ngvectors
      imap = map(ig)
      this%coords(1:3, ig) = ucoords(1:3, imap)
      this%red_gvec(1:3, ig) = ured_gvec(1:3, imap)
    end do
    
    SAFE_DEALLOCATE_A(ucoords)
    SAFE_DEALLOCATE_A(ured_gvec)
    SAFE_DEALLOCATE_A(modg2)
    SAFE_DEALLOCATE_A(map)

    POP_SUB(fourier_shell_init)
  end subroutine fourier_shell_init

  ! -----------------------------------------------------
  
  subroutine fourier_shell_end(this)
    type(fourier_shell_t), intent(inout) :: this
    
    PUSH_SUB(fourier_shell_end)

    SAFE_DEALLOCATE_P(this%coords)
    SAFE_DEALLOCATE_P(this%red_gvec)

    POP_SUB(fourier_shell_end)
  end subroutine fourier_shell_end

end module fourier_shell_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
