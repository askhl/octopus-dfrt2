!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: libxc.f90 3550 2007-11-19 14:32:49Z marques $

#if SINGLE_PRECISION
#  define XC_F90(x) xc_s_f90_ ## x
#else
#  define XC_F90(x) xc_f90_ ## x
#endif


!-------------------------------------------------------------------
module XC_F90(types_m)
#if SINGLE_PRECISION
  integer, public, parameter :: xc_f90_kind = selected_real_kind(4)
#else
  integer, public, parameter :: xc_f90_kind = selected_real_kind(14)
#endif

  type XC_F90(pointer_t)
    private
    integer, pointer :: buffer
  end type XC_F90(pointer_t)

end module XC_F90(types_m)


!-------------------------------------------------------------------
module XC_F90(lib_m)
  
  use XC_F90(types_m)
  use libxc_funcs_m

  implicit none

  public

  ! Families of xc functionals
  integer, parameter ::     &
    XC_FAMILY_UNKNOWN       =  -1,  &
    XC_FAMILY_NONE          =   0,  &
    XC_FAMILY_LDA           =   1,  &
    XC_FAMILY_GGA           =   2,  &
    XC_FAMILY_MGGA          =   4,  &
    XC_FAMILY_LCA           =   8,  &
    XC_FAMILY_OEP           =  16,  &
    XC_FAMILY_HYB_GGA       =  32

  integer, parameter ::             &
    XC_UNPOLARIZED          =   1,  &  ! Spin unpolarized
    XC_POLARIZED            =   2      ! Spin polarized

  integer, parameter ::             &
    XC_NON_RELATIVISTIC     =   0,  &  ! Functional includes or not relativistic
    XC_RELATIVISTIC         =   1      ! corrections. Only available in some functionals.

  ! Kinds
  integer, parameter ::             &
    XC_EXCHANGE             =   0,  &
    XC_CORRELATION          =   1,  &
    XC_EXCHANGE_CORRELATION =   2,  &
    XC_KINETIC              =   3

  integer, parameter ::             &
    XC_FLAGS_HAVE_EXC       =    1,  &
    XC_FLAGS_HAVE_VXC       =    2,  &
    XC_FLAGS_HAVE_FXC       =    4,  &
    XC_FLAGS_HAVE_KXC       =    8,  &
    XC_FLAGS_HAVE_LXC       =   16,  &
    XC_FLAGS_1D             =   32,  &
    XC_FLAGS_2D             =   64,  &
    XC_FLAGS_3D             =  128,  &
    XC_FLAGS_STABLE         =  512,  &
    XC_FLAGS_DEVELOPMENT    = 1024
 
  ! This value was redefined as XC_GGA_X_LB, we define it here to keep
  ! compatibility.
  integer, parameter :: XC_GGA_XC_LB = 160

  !----------------------------------------------------------------
  interface
    subroutine XC_F90(version)(major, minor)
      integer, intent(out) :: major, minor
    end subroutine XC_F90(version)
  end interface

  !----------------------------------------------------------------
  interface
    integer function XC_F90(info_number)(info)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in) :: info
    end function XC_F90(info_number)
    
    integer function XC_F90(info_kind)(info)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in) :: info
    end function XC_F90(info_kind)

    subroutine XC_F90(info_name)(info, s)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)),  intent(in)  :: info
      character(len=*), intent(out) :: s
    end subroutine XC_F90(info_name)

    integer function XC_F90(info_family)(info)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: info
    end function XC_F90(info_family)

    integer function XC_F90(info_flags)(info)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: info
    end function XC_F90(info_flags)

    subroutine XC_F90(info_refs)(info, number, str, s)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)),  intent(in)    :: info
      integer,                  intent(inout) :: number ! number of the reference. Must be 0 in the first call
      type(XC_F90(pointer_t)),  intent(inout) :: str    ! this will hold a (char **) pointer
      character(len=*),         intent(out)   :: s      ! the string that is output
    end subroutine XC_F90(info_refs)

    integer function XC_F90(family_from_id)(id)
      use XC_F90(types_m)
      integer, intent(in) :: id
    end function XC_F90(family_from_id)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(func_init)(p, info, functional, nspin)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(out) :: p
      type(XC_F90(pointer_t)), intent(out) :: info
      integer,                 intent(in)  :: functional
      integer,                 intent(in)  :: nspin
    end subroutine XC_F90(func_init)

    subroutine XC_F90(func_end)(p)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout) :: p
    end subroutine XC_F90(func_end)
  end interface


  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(lda)(p, np, rho, zk, vrho, fxc, kxc)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: zk    ! the energy per unit particle
      real(xc_f90_kind),       intent(out) :: vrho  ! v(nspin) the potential
      real(xc_f90_kind),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
      real(xc_f90_kind),       intent(out) :: kxc   ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine XC_F90(lda)

    subroutine XC_F90(lda_exc)(p, np, rho, zk)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: zk    ! the energy per unit particle
    end subroutine XC_F90(lda_exc)

    subroutine XC_F90(lda_exc_vxc)(p, np, rho, e, v)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: e     ! the energy per unit particle
      real(xc_f90_kind),       intent(out) :: v     ! v(nspin) the potential
    end subroutine XC_F90(lda_exc_vxc)

    subroutine XC_F90(lda_vxc)(p, np, rho, v)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: v     ! v(nspin) the potential
    end subroutine XC_F90(lda_vxc)

    subroutine XC_F90(lda_fxc)(p, np, rho, fxc)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine XC_F90(lda_fxc)

    subroutine XC_F90(lda_kxc)(p, np, rho, kxc)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),       intent(out) :: kxc
    end subroutine XC_F90(lda_kxc)
  end interface
  

  interface
    subroutine XC_F90(lda_x_1d_set_par)(p, interaction, bb)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout)  :: p
      integer,                 intent(in)     :: interaction
      real(xc_f90_kind),       intent(in)     :: bb
    end subroutine XC_F90(lda_x_1d_set_par)

    subroutine XC_F90(lda_c_xalpha_set_par)(p, alpha)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout)  :: p
      real(xc_f90_kind),      intent(in)     :: alpha
    end subroutine XC_F90(lda_c_xalpha_set_par)

    subroutine XC_F90(lda_x_set_par)(p, alpha, relativistic, omega)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout)  :: p
      real(xc_f90_kind),       intent(in)     :: alpha  ! of Xalpha, set to 4/3 to obtain standard LDA
      integer,                 intent(in)     :: relativistic
      real(xc_f90_kind),       intent(in)     :: omega
    end subroutine XC_F90(lda_x_set_par)

    subroutine XC_F90(lda_c_1d_csc_set_par)(p, interaction, bb)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout)  :: p
      integer,                 intent(in)     :: interaction
      real(xc_f90_kind),       intent(in)     :: bb
    end subroutine XC_F90(lda_c_1d_csc_set_par)

    subroutine XC_F90(lda_c_2d_prm_set_par)(p, N)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout)  :: p
      real(xc_f90_kind),       intent(in)     :: N
    end subroutine XC_F90(lda_c_2d_prm_set_par)
  end interface

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga)(p, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: zk
      real(xc_f90_kind),       intent(out) :: vrho
      real(xc_f90_kind),       intent(out) :: vsigma
      real(xc_f90_kind),       intent(out) :: v2rho2
      real(xc_f90_kind),       intent(out) :: v2rhosigma
      real(xc_f90_kind),       intent(out) :: v2sigma2
    end subroutine XC_F90(gga)

    subroutine XC_F90(gga_exc)(p, np, rho, sigma, zk)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,                 intent(in)  :: np
      real(xc_f90_kind),       intent(in)  :: rho
      real(xc_f90_kind),       intent(in)  :: sigma
      real(xc_f90_kind),       intent(out) :: zk
    end subroutine XC_F90(gga_exc)

    subroutine XC_F90(gga_exc_vxc)(p, np, rho, sigma, zk, vrho, vsigma)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
    end subroutine XC_F90(gga_exc_vxc)

    subroutine XC_F90(gga_vxc)(p, np, rho, sigma, vrho, vsigma)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
    end subroutine XC_F90(gga_vxc)

    subroutine XC_F90(gga_fxc)(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2sigma2
    end subroutine XC_F90(gga_fxc)
  end interface

  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_lb_set_par)(p, modified, threshold, ip, qtot)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: modified   ! should we use the modified version
      real(xc_f90_kind),    intent(in)  :: threshold  ! if so, the threshold to use the asymtotic version
      real(xc_f90_kind),    intent(in)  :: ip         ! ionization potential
      real(xc_f90_kind),    intent(in)  :: qtot       ! total charge
    end subroutine XC_F90(gga_lb_set_par)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_lb_modified)(p, np, rho, grho, r, dedd)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   ! rho(nspin) the density
      real(xc_f90_kind),    intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind),    intent(in)  :: r     ! distance from center of finite system
      real(xc_f90_kind),    intent(out) :: dedd
    end subroutine XC_F90(gga_lb_modified)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_x_wpbeh_set_par)(p, omega)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      real(xc_f90_kind),       intent(in)  :: omega       ! range separation
    end subroutine XC_F90(gga_x_wpbeh_set_par)
  end interface
  
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(gga_x_hjs_set_par)(p, omega)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      real(xc_f90_kind),       intent(in)  :: omega       ! range separation
    end subroutine XC_F90(gga_x_hjs_set_par)
  end interface

  !----------------------------------------------------------------
  interface
    subroutine XC_F90(hyb_exx_coef)(p, coef)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      real(xc_f90_kind),       intent(out) :: coef
    end subroutine XC_F90(hyb_exx_coef)

    subroutine XC_F90(hyb_cam_coef)(p, omega, alpha, beta)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      real(xc_f90_kind),       intent(out) :: omega, alpha, beta
    end subroutine XC_F90(hyb_cam_coef)
  end interface


  !----------------------------------------------------------------
  interface
    subroutine XC_F90(hyb_gga_xc_hse_set_par)(p, beta, omega)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      real(xc_f90_kind),       intent(in)  :: beta        ! mixing
      real(xc_f90_kind),       intent(in)  :: omega       ! range separation
    end subroutine XC_F90(hyb_gga_xc_hse_set_par)

    subroutine XC_F90(hyb_gga_xc_pbeh_set_par)(p, alpha)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      real(xc_f90_kind),       intent(in)  :: alpha       ! mixing
    end subroutine XC_F90(hyb_gga_xc_pbeh_set_par)
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine XC_F90(mgga)(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   
      real(xc_f90_kind),    intent(in)  :: sigma 
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau   
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2lapl2
      real(xc_f90_kind),    intent(out) :: v2tau2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2rholapl
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2sigmalapl
      real(xc_f90_kind),    intent(out) :: v2sigmatau
      real(xc_f90_kind),    intent(out) :: v2lapltau
    end subroutine XC_F90(mgga)

    subroutine XC_F90(mgga_exc)(p, np, rho, sigma, lapl, tau, zk)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
    end subroutine XC_F90(mgga_exc)

    subroutine XC_F90(mgga_exc_vxc)(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   
      real(xc_f90_kind),    intent(in)  :: sigma 
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: zk
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
    end subroutine XC_F90(mgga_exc_vxc)

    subroutine XC_F90(mgga_vxc)(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho   
      real(xc_f90_kind),    intent(in)  :: sigma 
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: vrho
      real(xc_f90_kind),    intent(out) :: vsigma
      real(xc_f90_kind),    intent(out) :: vlapl
      real(xc_f90_kind),    intent(out) :: vtau
    end subroutine XC_F90(mgga_vxc)

    subroutine XC_F90(mgga_fxc)(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(in)  :: p
      integer,              intent(in)  :: np
      real(xc_f90_kind),    intent(in)  :: rho
      real(xc_f90_kind),    intent(in)  :: sigma
      real(xc_f90_kind),    intent(in)  :: lapl
      real(xc_f90_kind),    intent(in)  :: tau
      real(xc_f90_kind),    intent(out) :: v2rho2
      real(xc_f90_kind),    intent(out) :: v2sigma2
      real(xc_f90_kind),    intent(out) :: v2lapl2
      real(xc_f90_kind),    intent(out) :: v2tau2
      real(xc_f90_kind),    intent(out) :: v2rhosigma
      real(xc_f90_kind),    intent(out) :: v2rholapl
      real(xc_f90_kind),    intent(out) :: v2rhotau
      real(xc_f90_kind),    intent(out) :: v2sigmalapl
      real(xc_f90_kind),    intent(out) :: v2sigmatau
      real(xc_f90_kind),    intent(out) :: v2lapltau
    end subroutine XC_F90(mgga_fxc)
  end interface

  interface
    subroutine XC_F90(mgga_x_tb09_set_par)(p, cc)
      use XC_F90(types_m)
      type(XC_F90(pointer_t)), intent(inout)  :: p
      real(xc_f90_kind),    intent(in)     :: cc
    end subroutine XC_F90(mgga_x_tb09_set_par)
  end interface


end module XC_F90(lib_m)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
