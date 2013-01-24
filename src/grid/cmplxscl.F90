!! Copyright (C) 2013 U. De Giovannini
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
!! $Id: linear_response.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module cmplxscl_m
  use datasets_m
  use global_m
  use parser_m
  use messages_m
  use varinfo_m

  implicit none
  
  private


  public :: &
       cmplxscl_t, &
       cmplxscl_init, & 
       cmplxscl_end 
       
  !> Complex scaling module    
  type cmplxscl_t
    logical            :: space       !< scale spatial coordinates
    logical            :: time        !< scale time coordinate
  
     FLOAT             :: theta       !< spatial coordinates scaling angle
     FLOAT             :: alphaR     !< time coordinate scaling angle for right states
     FLOAT             :: alphaL     !< time coordinate scaling angle for left states

  end type cmplxscl_t

contains

  !-----------------------------------------------------------------
  subroutine cmplxscl_init(this)
    type(cmplxscl_t),    intent(out) :: this


    PUSH_SUB(cmplxscl_init)

    !%Variable ComplexScalingSpace
    !%Type logical
    !%Default false
    !%Section Hamiltonian
    !%Description
    !% (experimental) If set to yes, a complex scaled Hamiltonian will be used. 
    !% This is implemented through a global coordinate transformation r-> r e^{i \theta}.
    !% When <tt>TheoryLevel=DFT</tt> Density functional resonance theory DFRT is employed.  
    !% In order to reveal resonances <tt>ComplexScalingTheta</tt> bigger than zero should be set.
    !% D. L. Whitenack and A. Wasserman, Phys. Rev. Lett. 107, 163002 (2011).
    !%End
    call parse_logical(datasets_check('ComplexScalingSpace'), .false., this%space)

    !%Variable ComplexScalingTime
    !%Type logical
    !%Default false
    !%Section Hamiltonian
    !%Description
    !% (experimental) If set to yes, the time evolution will be performed with a 
    !% complex time variable. 
    !% This is implemented by a coordinate transformation t-> t e^{i \alpha_r} for 
    !% right states and t-> t e^{i \alpha_l} for left states.
    !% J. Bengtsson, E. Lindroth, and S. Selstø, Phys. Rev. A 85, 013419 (2012).
    !%End
    call parse_logical(datasets_check('ComplexScalingTime'), .false., this%time)

    
    !%Variable ComplexScalingTheta
    !%Type float 
    !%Default 0.3
    !%Section Hamiltonian
    !%Description
    !% The spatial coordinate complex scaling angle \theta.
    !% Allowed values must be in the range 0 <= \theta < \pi/4. 
    !%End
    call parse_float(datasets_check('ComplexScalingTheta'), CNST(0.3), this%theta)
    if(this%theta < M_ZERO .or. this%theta > M_PI/CNST(4.0)) call input_error('ComplexScalingTheta')



    !%Variable ComplexScalingAlphaRight
    !%Type float 
    !%Default 2*theta
    !%Section Hamiltonian
    !%Description
    !% The time coordinate complex scaling angle \alpha_r used to evolve 
    !% right states.  
    !%End
    call parse_float(datasets_check('ComplexScalingAlphaRight'), M_TWO*this%theta, this%alphaR)

    !%Variable ComplexScalingAlphaLeft
    !%Type float 
    !%Default -2*theta
    !%Section Hamiltonian
    !%Description
    !% The time coordinate complex scaling angle \alpha_l used to evolve 
    !% left states.  
    !%End
    call parse_float(datasets_check('ComplexScalingAlphaLeft'), M_TWO*this%theta, this%alphaL)



    if (this%space .or. this%time) then
      call messages_print_stress(stdout, "Complex Scaling")
    end if
    if (this%space) then
      write(message(1), '(a,f12.3)') 'Complex scaling angle theta = ', this%theta
      call messages_info(1)
    end if
      
    if (this%time) then
      write(message(1), '(a,f12.3)') 'Complex scaling angle alphaR = ', this%alphaR
      write(message(2), '(a,f12.3)') 'Complex scaling angle alphaL = ', this%alphaL
      call messages_info(2)
    end if

    if (this%space .or. this%time) then
      call messages_print_stress(stdout)
    end if


    POP_SUB(cmplxscl_init)

  end subroutine cmplxscl_init


  !-----------------------------------------------------------------
  subroutine cmplxscl_end(this)
    type(cmplxscl_t),     intent(inout) :: this

    PUSH_SUB(cmplxscl_end)


    POP_SUB(cmplxscl_end)
  end subroutine cmplxscl_end


  !-----------------------------------------------------------------
  subroutine cmplxscl_copy(in, out)
    type(cmplxscl_t),     intent(in) :: in
    type(cmplxscl_t),     intent(out) :: out

    PUSH_SUB(cmplxscl_copy)

    out%space  = in%space
    out%time   = in%time
    out%theta  = in%theta
    out%alphaL = in%alphaL
    out%alphaR = in%alphaR


    POP_SUB(cmplxscl_copy)
  end subroutine cmplxscl_copy


end module cmplxscl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
