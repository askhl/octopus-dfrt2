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
!! $Id: vxc_inc.F90 9152 2012-06-21 08:30:22Z umberto $

subroutine dxc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, ex, ec, deltaxc, vxc, vtau)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: deltaxc         !< The XC derivative descontinuity
  FLOAT, optional,      intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)

  integer :: n_block

  FLOAT, allocatable :: l_zk(:)        ! Local block of the energy functional (with the correct memory order for libxc)
  FLOAT, allocatable :: l_dens(:,:)    ! Local block for the density 
  FLOAT, allocatable :: l_dedd(:,:)    ! Local block of the exchange or correl. potential(with the correct memory order for libxc)
  FLOAT, allocatable :: l_sigma(:,:)   
  FLOAT, allocatable :: l_vsigma(:,:)  
  FLOAT, allocatable :: l_tau(:,:)
  FLOAT, allocatable :: l_ldens(:,:)
  FLOAT, allocatable :: l_dedtau(:,:)
  FLOAT, allocatable :: l_dedldens(:,:)

  FLOAT, allocatable :: dens(:,:)      ! Density
  FLOAT, allocatable :: dedd(:,:)      ! (Functional) Derivative of the exchange or correlation energy with
  ! respect to the density (vector used to store the exchange or the correlation potential)
  FLOAT, allocatable :: ex_per_vol(:)  ! Exchange energy per unit volume 
  FLOAT, allocatable :: ec_per_vol(:)  ! Correlation energy per unit volume 
  FLOAT, allocatable :: gdens(:,:,:)   ! Gradient of the density
  FLOAT, allocatable :: dedgd(:,:,:)   ! (Functional) Derivative of the exchange or correlation energy with
  !respect to the gradient of the density.
  FLOAT, allocatable :: current(:,:,:) ! Paramagnetic or total current
  FLOAT, allocatable :: ldens(:,:)     ! Laplacian of the density
  FLOAT, allocatable :: tau(:,:)       ! Kinetic energy density
  FLOAT, allocatable :: dedldens(:,:)  ! (Functional) Derivative of the exchange or correlation energy with
  !respect to the laplacian of the density.
  FLOAT, allocatable :: symmtmp(:, :)  ! Temporary vector for the symmetrizer
  FLOAT, allocatable :: vx(:, :)
  FLOAT, allocatable :: gf(:,:)
  integer :: ib, ib2, ip, isp, families, ixc, spin_channels, is
  integer, save :: xc_get_vxc_counter = 0
  FLOAT   :: rr,ipot_to_pass
  logical :: gga, mgga
  type(profile_t), save :: prof
  logical :: calc_energy
  type(xc_functl_t), pointer :: functl(:)
  type(symmetrizer_t) :: symmetrizer

  PUSH_SUB(dxc_get_vxc)
  call profiling_in(prof, "dXC_LOCAL")

  ASSERT(present(ex) .eqv. present(ec))
  calc_energy = present(ex)

  xc_get_vxc_counter = xc_get_vxc_counter + 1
  !xprint*, "xc_get_vxc call number ", xc_get_vxc_counter

  !Pointer-shortcut for xcs%functl
  !It helps to remember that for xcs%functl(:,:)
  ! (1,:) => exchange,    (2,:) => correlation
  ! (:,1) => unpolarized, (:,2) => polarized
  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if

  ! is there anything to do ?
  families = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA
  if(iand(xcs%family, families) == 0) then
    POP_SUB(dxc_get_vxc)
    call profiling_out(prof)
    return
  endif

  n_block = 1000

  ! initialize a couple of handy variables
  gga  = iand(xcs%family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA).ne.0
  mgga = iand(xcs%family, XC_FAMILY_MGGA).ne.0

  !Read the spin channels
  !Index 1 refers to "exchange mode" of type "xc_functl_t". For this purpose using 1 or 2 makes no difference. 
  spin_channels = functl(1)%spin_channels

  if(xcs%xc_density_correction == LR_X) then
    SAFE_ALLOCATE(vx(1:der%mesh%np, 1:spin_channels))
  end if

  call lda_init()
  if(gga .or. xcs%xc_density_correction == LR_X) call  gga_init()
  if(mgga) call mgga_init()

  ! Get the gradient and the Laplacian of the density and the kinetic-energy density
  ! We do it here instead of doing it in gga_init and mgga_init in order to 
  ! avoid calling the subroutine states_calc_quantities twice
  if((gga .and. (.not. mgga)) .or. xcs%xc_density_correction == LR_X) then
    ! get gradient of the density (this is faster than calling states_calc_quantities)
    do isp = 1, spin_channels 
      call dderivatives_grad(der, dens(:, isp), gdens(:, :, isp)) 
    end do
  else if(mgga) then
    if (xc_get_vxc_counter .le. 2) then 
      ipot_to_pass = M_ONE
    else
      ipot_to_pass = ioniz_pot
    end if
    !print *, "Ioniz potential to pass =", ipot_to_pass
    ! We calculate everything from the wavefunctions to benefit from
    ! the error cancellation between the gradient of the density and
    ! tau.
    !
    ! FIXME: Probably this is wrong for non-local corrections or other
    ! cases when the density is not directly generated by the
    ! orbitals.
    call states_calc_quantities(der, st, gi_kinetic_energy_density = tau, density_gradient = gdens, density_laplacian = ldens)

    ! We have to symmetrize everything as they are calculated from the
    ! wavefunctions.
    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmtmp(1:der%mesh%np, 1:3))
      call symmetrizer_init(symmetrizer, der%mesh)
      do isp = 1, spin_channels
        call dsymmetrizer_apply(symmetrizer, tau(:, isp), symmtmp(:, 1))
        tau(1:der%mesh%np, isp) = symmtmp(1:der%mesh%np, 1)
        call dsymmetrizer_apply(symmetrizer, ldens(:, isp), symmtmp(:, 1))
        ldens(1:der%mesh%np, isp) = symmtmp(1:der%mesh%np, 1)
        call dsymmetrizer_apply_vector(symmetrizer, gdens(:, :, isp), symmtmp)
        gdens(1:der%mesh%np, 1:3, isp) = symmtmp(1:der%mesh%np, 1:3)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmtmp)
    end if

    if(functl(1)%id == XC_MGGA_X_TB09 .and. der%mesh%sb%periodic_dim == 3) then
      call calc_tb09_c()
    end if

  end if

  space_loop: do ip = 1, der%mesh%np, n_block

    !Resize the dimension of the last block when the number of the mesh points
    !it is not a perfect divider of the dimension of the blocks.
    if(ip + n_block > der%mesh%np) n_block = der%mesh%np - ip + 1 

    ! make a local copy with the correct memory order for libxc
    ib2 = ip
    do ib = 1, n_block
      l_dens(1:spin_channels, ib) = dens(ib2, 1:spin_channels)
      ib2 = ib2 + 1
    end do

    if(gga) then
      ib2 = ip
      do ib = 1, n_block
        l_sigma(1, ib) = sum(gdens(ib2, 1:der%mesh%sb%dim, 1)*gdens(ib2, 1:der%mesh%sb%dim, 1))
        if(ispin /= UNPOLARIZED) then
          ! memo: please check the following indices
          l_sigma(2, ib) = sum(gdens(ib2, 1:der%mesh%sb%dim, 1)*gdens(ib2, 1:der%mesh%sb%dim, 2)) 
          l_sigma(3, ib) = sum(gdens(ib2, 1:der%mesh%sb%dim, 2)*gdens(ib2, 1:der%mesh%sb%dim, 2))
        end if
        ib2 = ib2 + 1
      end do
    end if

    if(mgga) then
      ib2 = ip
      do ib = 1, n_block
        l_tau  (1:spin_channels, ib) =   tau(ib2, 1:spin_channels)
        l_ldens(1:spin_channels, ib) = ldens(ib2, 1:spin_channels)
        ib2 = ib2 + 1
      end do
    end if

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = 1, 2

      if(.not. present(vxc)) then ! get only the xc energy

        if(iand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC).ne.0) then
          select case(functl(ixc)%family)

          case(XC_FAMILY_LDA)
            call XC_F90(lda_exc)(functl(ixc)%conf, n_block, l_dens(1,1), l_zk(1))

          case(XC_FAMILY_GGA)
            call XC_F90(gga_exc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_zk(1))

          case(XC_FAMILY_HYB_GGA)
            message(1) = 'Hyb-GGAs are currently disabled.'
            call messages_fatal(1)

          case(XC_FAMILY_MGGA)
            call XC_F90(mgga_exc)(functl(ixc)%conf, n_block, &
              l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), l_zk(1))

          case default
            cycle
          end select

        else ! Do not have an energy functional
          l_zk(:) = M_ZERO
        end if

      else if(calc_energy .and. iand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC).ne.0) then
        ! we get the xc energy and potential
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_zk(1), l_dedd(1,1))

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call XC_F90(gga_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
            l_zk(1), l_dedd(1,1), l_vsigma(1,1))

        case(XC_FAMILY_MGGA)
          call XC_F90(mgga_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), &
            l_zk(1), l_dedd(1,1), l_vsigma(1,1), l_dedldens(1,1), l_dedtau(1,1))

        case default
          cycle
        end select

      else ! we just get the potential
        l_zk(:) = M_ZERO

        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_dedd(1,1))

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          l_vsigma = M_ZERO

          if(functl(ixc)%id == XC_GGA_X_LB) then
            call mesh_r(der%mesh, ip, rr)
            call XC_F90(gga_lb_modified)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
              rr, l_dedd(1,1))
          else
            call XC_F90(gga_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
              l_dedd(1,1), l_vsigma(1,1))
          end if

        case(XC_FAMILY_MGGA)
          call XC_F90(mgga_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), &
            l_dedd(1,1), l_vsigma(1,1), l_dedldens(1,1), l_dedtau(1,1))

        case default
          cycle
        end select

      end if

      if(calc_energy) then
        ib2 = ip
        if(functl(ixc)%type == XC_EXCHANGE) then
          do ib = 1, n_block
            ex_per_vol(ib2) = ex_per_vol(ib2) + sum(l_dens(1:spin_channels, ib)) * l_zk(ib)
            ib2 = ib2 + 1
          end do
        else
          do ib = 1, n_block
            ec_per_vol(ib2) = ec_per_vol(ib2) + sum(l_dens(1:spin_channels, ib)) * l_zk(ib)
            ib2 = ib2 + 1
          end do
        end if
      end if

      ! store results
      if(present(vxc)) then

        if(xcs%xc_density_correction == LR_X .and. &
          (functl(ixc)%type == XC_EXCHANGE .or. functl(ixc)%type == XC_EXCHANGE_CORRELATION)) then

          ! in this case, vx goes to a different array
          ib2 = ip
          do ib = 1, n_block
            vx(ib2, 1:spin_channels) = l_dedd(1:spin_channels, ib)
            ib2 = ib2 + 1
          end do
          
        else
          ib2 = ip
          do ib = 1, n_block
            dedd(ib2, 1:spin_channels) = dedd(ib2, 1:spin_channels) + l_dedd(1:spin_channels, ib)
            ib2 = ib2 + 1
          end do

        end if

        if((functl(ixc)%family == XC_FAMILY_GGA).or.(functl(ixc)%family == XC_FAMILY_MGGA)) then
          ib2 = ip
          do ib = 1, n_block
            dedgd(ib2,:,1) = dedgd(ib2,:,1) + M_TWO*l_vsigma(1, ib)*gdens(ib2,:,1)
            if(ispin /= UNPOLARIZED) then
              dedgd(ib2,:,1) = dedgd(ib2,:,1) + l_vsigma(2, ib)*gdens(ib2,:,2)
              dedgd(ib2,:,2) = dedgd(ib2,:,2) +  &
                M_TWO*l_vsigma(3, ib)*gdens(ib2,:,2) + l_vsigma(2, ib)*gdens(ib2,:,1)
            end if
            ib2 = ib2 + 1
          end do
        end if

        if(functl(ixc)%family == XC_FAMILY_MGGA) then
          ib2 = ip
          do ib = 1, n_block
            dedldens(ib2, 1:spin_channels) = dedldens(ib2, 1:spin_channels) + l_dedldens(1:spin_channels, ib)
            vtau    (ib2, 1:spin_channels) = vtau    (ib2, 1:spin_channels) + l_dedtau  (1:spin_channels, ib)
            ib2 = ib2 + 1
          end do
        end if
      end if

    end do functl_loop
  end do space_loop

  if(present(deltaxc)) then
    deltaxc = M_ZERO
  end if

  if(xcs%xc_density_correction == LR_XC) then
    call xc_density_correction_calc(xcs, der, spin_channels, rho, dedd, deltaxc = deltaxc)
  end if

  if(xcs%xc_density_correction == LR_X) then
    call xc_density_correction_calc(xcs, der, spin_channels, rho, vx, deltaxc = deltaxc)
    dedd(1:der%mesh%np, 1:spin_channels) = dedd(1:der%mesh%np, 1:spin_channels) + vx(1:der%mesh%np, 1:spin_channels)

    if(calc_energy) then
      ! get the energy density from Levy-Perdew
      SAFE_ALLOCATE(gf(1:der%mesh%np, 1:3))
      
      do is = 1, spin_channels
        do ip = 1, der%mesh%np
          ex_per_vol(ip) = vx(ip, is)*(CNST(3.0)*rho(ip, is) &
            + sum(der%mesh%x(ip, 1:der%mesh%sb%dim)*gdens(ip, 1:der%mesh%sb%dim, is)))
        end do
      end do
      
      SAFE_DEALLOCATE_A(gf)
    end if
  end if
  
!  print*, deltaxc

  ! this has to be done in inverse order
  if(present(vxc)) then
    if(mgga) call mgga_process()
    if( gga) call  gga_process()
    call lda_process()
  end if

  if(calc_energy) then
    ! integrate energies per unit volume
    ex = ex + dmf_integrate(der%mesh, ex_per_vol)
    ec = ec + dmf_integrate(der%mesh, ec_per_vol)
  end if

  ! clean up allocated memory
  call lda_end()
  if(gga .or. xcs%xc_density_correction == LR_X) call  gga_end()
  if(mgga) call mgga_end()

  POP_SUB(dxc_get_vxc)
  call profiling_out(prof)

contains

  ! ---------------------------------------------------------
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates density and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: ii
    FLOAT   :: d(2), dtot, dpol

    PUSH_SUB(xc_get_vxc.lda_init)

    ! allocate some general arrays
    SAFE_ALLOCATE(l_dens(1:spin_channels, 1:n_block))
    SAFE_ALLOCATE(l_zk(1:n_block))

    SAFE_ALLOCATE(dens(1:der%mesh%np_part, 1:spin_channels))
    dens       = M_ZERO

    if(calc_energy) then
      SAFE_ALLOCATE(ex_per_vol(1:der%mesh%np))
      SAFE_ALLOCATE(ec_per_vol(1:der%mesh%np))
      ex_per_vol = M_ZERO
      ec_per_vol = M_ZERO
    end if

    if(present(vxc)) then
      SAFE_ALLOCATE(l_dedd(1:spin_channels, 1:n_block))
      SAFE_ALLOCATE(dedd(1:der%mesh%np_part, 1:spin_channels))
      dedd = M_ZERO
    end if

    do ii = 1, der%mesh%np
      d(1:spin_channels) = rho(ii, 1:spin_channels)

      select case(ispin)
      case(UNPOLARIZED)
        dens(ii, 1) = max(d(1), M_ZERO)
      case(SPIN_POLARIZED)
        dens(ii, 1) = max(d(1), M_ZERO)
        dens(ii, 2) = max(d(2), M_ZERO)
      case(SPINORS)
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(ii, 3)**2 + rho(ii, 4)**2))
        dens(ii, 1) = max(M_HALF*(dtot + dpol), M_ZERO)
        dens(ii, 2) = max(M_HALF*(dtot - dpol), M_ZERO)
      end select
    end do

    POP_SUB(xc_get_vxc.lda_init)
  end subroutine lda_init


  ! ---------------------------------------------------------
  ! deallocate variables allocated in lda_init
  subroutine lda_end()
    PUSH_SUB(xc_get_vxc.lda_end)

    SAFE_DEALLOCATE_A(l_dens)
    SAFE_DEALLOCATE_A(l_zk)

    SAFE_DEALLOCATE_A(dens)
    SAFE_DEALLOCATE_A(ex_per_vol)
    SAFE_DEALLOCATE_A(ec_per_vol)
    if(present(vxc)) then
      SAFE_DEALLOCATE_A(l_dedd)
      SAFE_DEALLOCATE_A(dedd)
    end if

    POP_SUB(xc_get_vxc.lda_end)
  end subroutine lda_end


  ! ---------------------------------------------------------
  ! calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: ip
    FLOAT :: d(2), dpol, vpol

    PUSH_SUB(xc_get_vxc.lda_process)

    if(ispin == SPINORS) then
      ! rotate back (do not need the rotation matrix for this).
      do ip = 1, der%mesh%np
        d(1:spin_channels) = rho(ip, 1:spin_channels)

        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(ip, 3)**2 + rho(ip, 4)**2))
        vpol = (dedd(ip, 1) - dedd(ip, 2))*(d(1) - d(2))/(dpol + tiny)

        vxc(ip, 1) = vxc(ip, 1) + M_HALF*(dedd(ip, 1) + dedd(ip, 2) + vpol)
        vxc(ip, 2) = vxc(ip, 2) + M_HALF*(dedd(ip, 1) + dedd(ip, 2) - vpol)
        vxc(ip, 3) = vxc(ip, 3) + (dedd(ip, 1) - dedd(ip, 2))*rho(ip, 3)/(dpol + tiny)
        vxc(ip, 4) = vxc(ip, 4) + (dedd(ip, 1) - dedd(ip, 2))*rho(ip, 4)/(dpol + tiny)
      end do
    elseif(ispin == SPIN_POLARIZED) then
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:, 1), vxc(:, 1))
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:, 2), vxc(:, 2))
    else
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:, 1), vxc(:, 1))
    end if

    POP_SUB(xc_get_vxc.lda_process)
  end subroutine lda_process


  ! ---------------------------------------------------------
  ! initialize GGAs
  !   *) allocates gradient of the density (gdens), dedgd, and its local variants

  subroutine gga_init()
    integer :: ii

    PUSH_SUB(xc_get_vxc.gga_init)

    ii = 1
    if(ispin /= UNPOLARIZED) ii = 3

    ! allocate variables
    SAFE_ALLOCATE(l_sigma(1:ii, 1:n_block))
    SAFE_ALLOCATE(gdens(1:der%mesh%np, 1:3, 1:spin_channels))
    gdens = M_ZERO

    if(present(vxc)) then
      SAFE_ALLOCATE(l_vsigma(1:ii, 1:n_block))
      SAFE_ALLOCATE(dedgd(1:der%mesh%np_part, 1:3, 1:spin_channels))
      dedgd = M_ZERO
    end if

    do ii = 1, 2
      if(functl(ii)%id == XC_GGA_X_LB) then
        call XC_F90(gga_lb_set_par)(functl(ii)%conf, &
          functl(ii)%LB94_modified, functl(ii)%LB94_threshold, ioniz_pot, qtot)
      end if
    end do

    POP_SUB(xc_get_vxc.gga_init)
  end subroutine gga_init


  ! ---------------------------------------------------------
  ! cleans up memory allocated in gga_init
  subroutine gga_end()
    PUSH_SUB(xc_get_vxc.gga_end)

    SAFE_DEALLOCATE_A(l_sigma)
    SAFE_DEALLOCATE_A(gdens)
    if(present(vxc)) then
      SAFE_DEALLOCATE_A(l_vsigma)
      SAFE_DEALLOCATE_A(dedgd)
    end if

    POP_SUB(xc_get_vxc.gga_end)
  end subroutine gga_end


  ! ---------------------------------------------------------
  ! calculates the GGA contribution to vxc
  subroutine gga_process()
    integer :: ip, is
    FLOAT, allocatable :: gf(:,:)

    PUSH_SUB(xc_get_vxc.gga_process)

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    SAFE_ALLOCATE(gf(1:der%mesh%np, 1:1))
    do is = 1, spin_channels
      call dderivatives_div(der, dedgd(:, :, is), gf(:, 1))
      call lalg_axpy(der%mesh%np, -M_ONE, gf(:,1), dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(gf)

    ! If LB94, we can calculate an approximation to the energy from
    ! Levy-Perdew relation PRA 32, 2010 (1985)
    if(calc_energy .and. functl(1)%id == XC_GGA_X_LB) then
      SAFE_ALLOCATE(gf(1:der%mesh%np, 1:3))

      do is = 1, spin_channels
        call dderivatives_grad(der, dedd(:, is), gf(:,:))
        do ip = 1, der%mesh%np
          ex_per_vol(ip) = ex_per_vol(ip) - dens(ip, is) * sum(der%mesh%x(ip,:)*gf(ip,:))
        end do
      end do

      SAFE_DEALLOCATE_A(gf)
    end if

    POP_SUB(xc_get_vxc.gga_process)
  end subroutine gga_process


  ! ---------------------------------------------------------
  ! initialize meta-GGAs
  !   *) allocate the kinetic-energy density, dedtau, and local variants

  subroutine mgga_init()
    PUSH_SUB(xc_get_vxc.mgga_init)

    ! allocate variables
    SAFE_ALLOCATE( tau(1:der%mesh%np, 1:spin_channels))
    SAFE_ALLOCATE( current (1:der%mesh%np, 1:der%mesh%sb%dim, 1:spin_channels) )
    SAFE_ALLOCATE(ldens(1:der%mesh%np, 1:spin_channels))

    SAFE_ALLOCATE(l_tau  (1:spin_channels, 1:n_block))
    SAFE_ALLOCATE(l_ldens(1:spin_channels, 1:n_block))

    if(present(vxc)) then
      SAFE_ALLOCATE(dedldens(1:der%mesh%np_part, 1:spin_channels))
      dedldens = M_ZERO

      SAFE_ALLOCATE(l_dedtau  (1:spin_channels, 1:n_block))
      SAFE_ALLOCATE(l_dedldens(1:spin_channels, 1:n_block))
    end if

    POP_SUB(xc_get_vxc.mgga_init)
  end subroutine mgga_init


  ! ---------------------------------------------------------

  subroutine calc_tb09_c()
    FLOAT, allocatable :: gnon(:)
    FLOAT gn(MAX_DIM), n, tb09_c
    integer :: ii

    PUSH_SUB(xc_get_vxc.calc_tb09_c)

    SAFE_ALLOCATE(gnon(1:der%mesh%np))

    do ii = 1, der%mesh%np
      if(ispin == UNPOLARIZED) then
        n = dens(ii, 1)
        gn(1:der%mesh%sb%dim) = gdens(ii, 1:der%mesh%sb%dim, 1)
      else
        n = dens(ii, 1) + dens(ii, 2)
        gn(1:der%mesh%sb%dim) = gdens(ii, 1:der%mesh%sb%dim, 1) + gdens(ii, 1:der%mesh%sb%dim, 2)
      end if

      if (n <= CNST(1e-7)) then 
        gnon(ii) = CNST(0.0)
        ! here you will have to print the true gnon(ii) with the correspondent mesh point ii
      else
        gnon(ii) = sqrt(sum((gn(1:der%mesh%sb%dim)/n)**2))
      end if
    end do

    tb09_c =  -CNST(0.012) + CNST(1.023)*sqrt(dmf_integrate(der%mesh, gnon)/der%mesh%sb%rcell_volume)

    write(message(1), '(a,f13.6)') "Info: In the functional TB09 c = ", tb09_c
    call messages_info(1)

    call  XC_F90(mgga_x_tb09_set_par)(functl(1)%conf, tb09_c)

    SAFE_DEALLOCATE_A(gnon)

    POP_SUB(xc_get_vxc.calc_tb09_c)
  end subroutine calc_tb09_c


  ! ---------------------------------------------------------
  ! clean up memory allocates in mgga_init
  subroutine mgga_end()
    PUSH_SUB(xc_get_vxc.mgga_end)

    SAFE_DEALLOCATE_A(tau)
    SAFE_DEALLOCATE_A(current)
    SAFE_DEALLOCATE_A(ldens)

    SAFE_DEALLOCATE_A(l_tau)
    SAFE_DEALLOCATE_A(l_ldens)

    if(present(vxc)) then
      SAFE_DEALLOCATE_A(dedldens)

      SAFE_DEALLOCATE_A(l_dedtau)
      SAFE_DEALLOCATE_A(l_dedldens)
    end if

    POP_SUB(xc_get_vxc.mgga_end)
  end subroutine mgga_end


  ! ---------------------------------------------------------
  ! calculate the mgga contribution to vxc
  subroutine mgga_process()
    integer :: is
    FLOAT, allocatable :: lf(:,:)

    PUSH_SUB(xc_get_vxc.mgga_process)

    ! add the Laplacian of the functional derivative of Exc with respect to
    ! the laplacian of the density.

    SAFE_ALLOCATE(lf(1:der%mesh%np, 1:1))
    do is = 1, spin_channels
      call dderivatives_lapl(der, dedldens(:, is), lf(:, 1))
      call lalg_axpy(der%mesh%np, M_ONE, lf(:, 1), dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(lf)

    POP_SUB(xc_get_vxc.mgga_process)
  end subroutine mgga_process

end subroutine dxc_get_vxc

! -----------------------------------------------------

subroutine xc_density_correction_calc(xcs, der, nspin, density, vxc, deltaxc)
  type(xc_t),          intent(in)    :: xcs
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: nspin
  FLOAT,               intent(in)    :: density(:, :)
  FLOAT,               intent(inout) :: vxc(:, :)
  FLOAT, optional,     intent(out)   :: deltaxc

  logical :: find_root, done
  integer :: ip, iunit, ierr
  integer, save :: iter = 0
  FLOAT,   save :: ncsave
  character(len=30) :: number
  FLOAT   :: qxc, ncutoff, qxc_old, ncutoff_old, deriv, deriv_old, qxcfin
  FLOAT   :: x1, x2, x3, f1, f2, f3, dd, vol, mindd, maxdd
  FLOAT, allocatable :: nxc(:), lrvxc(:)
  type(profile_t), save :: prof
  FLOAT, parameter :: thres = CNST(1e-6)

  PUSH_SUB('vxc_inc.xc_density_correction_calc')

  call profiling_in(prof, "XC_DENSITY_CORRECTION")

  SAFE_ALLOCATE(nxc(1:der%mesh%np))
  SAFE_ALLOCATE(lrvxc(1:der%mesh%np_part))

  if(nspin > 1) then
    call messages_not_implemented('XC density correction with spin polarization')
  end if

  forall(ip = 1:der%mesh%np) lrvxc(ip) = CNST(-1.0)/(CNST(4.0)*M_PI)*vxc(ip, 1)
  call dderivatives_lapl(der, lrvxc, nxc)

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "rho", der%mesh, density(:, 1), unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "vxcorig", der%mesh, vxc(:, 1), unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "nxc", der%mesh, nxc, unit_one, ierr)

  if(xcs%xcd_optimize_cutoff) then

    x1 = CNST(1.0e-8)
    qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)
    deriv = HUGE(deriv)
    done = .false.

    INCR(iter, 1)
    if(mpi_world%rank == 0) then
      write(number, '(i4)') iter
      iunit = io_open('qxc.'//trim(adjustl(number)), action='write')
    end if
    do
      if(.not. done) then
        ncutoff_old = x1
        qxc_old = qxc
        deriv_old = deriv
      end if

      x1 = x1*CNST(1.01)

      if(mpi_world%rank == 0) then
        write(iunit, *) x1, qxc
      end if

      if(x1 > CNST(1.0)) exit

      qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)

      if(qxc == qxc_old) cycle

      deriv = (qxc - qxc_old)/(x1 - ncutoff_old)

      if(.not. done .and. abs(qxc) >= 1.0_8) then
        find_root = .true.
        done = .true.
        ncutoff = x1
      end if

      if(xcs%xcd_minimum .and. .not. done .and. abs(qxc_old) - abs(qxc) > thres) then
        find_root = .false.
        done = .true.
        ncutoff = x1
        print*, x1, 0
      end if

    end do

    if(mpi_world%rank == 0) call io_close(iunit)

    if(iter > 1) x3 = ncsave

    if(find_root) then
      x1 = ncutoff
      x2 = ncutoff_old
      x3 = ncutoff
      f1 = 1.0_8 + get_qxc(der%mesh, nxc, density(:, 1), x1)
      f2 = 1.0_8 + get_qxc(der%mesh, nxc, density(:, 1), x2)

      do ip = 1, 20
        if(abs(f1 - f2) < 1e-16_8) exit
        x3 = x2 - f2*(x2 - x1)/(f2 - f1)
        f3 = 1.0_8 + get_qxc(der%mesh, nxc, density(:, 1), x3)
        if(abs(f3) < 1e-6_8) exit
        x1 = x2
        f1 = f2
        x2 = x3
        f2 = f3
      end do

      if(x3 <= ncutoff) ncutoff = x3
    end if

    ncsave = x3

  else
    ncutoff = xcs%xcd_ncutoff
  end if

  qxcfin = get_qxc(der%mesh, nxc, density(:, 1), ncutoff)

  do ip = 1, der%mesh%np
    if(density(ip, 1) < ncutoff) then
      nxc(ip) = -nxc(ip)
    else
      nxc(ip) = M_ZERO
    end if
  end do

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "nxcmod", der%mesh, nxc, unit_one, ierr)

  if(mpi_world%rank == 0) then
    print*, "Iter",    iter, ncutoff, qxcfin
  end if

  call dpoisson_solve(psolver, lrvxc, nxc)

  if(xcs%xcd_normalize .and. abs(qxcfin) > CNST(1e-10)) then
    do ip = 1, der%mesh%np
      lrvxc(ip) = lrvxc(ip)/abs(qxcfin)
    end do
  end if

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Y, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_PLANE_X, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_PLANE_Y, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_PLANE_Z, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)

  forall(ip = 1:der%mesh%np) vxc(ip, 1) = vxc(ip, 1) + lrvxc(ip)

  maxdd = -HUGE(maxdd)
  mindd =  HUGE(maxdd)
  vol = M_ZERO
  do ip = 1, der%mesh%np
    if(density(ip, 1) >= ncutoff) then
      vol = vol + der%mesh%volume_element
      maxdd = max(lrvxc(ip), maxdd)
      mindd = min(lrvxc(ip), mindd)
    else
      lrvxc(ip) = M_ZERO
    end if
  end do

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Y, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)

  dd = dmf_integrate(der%mesh, lrvxc)/vol

  if(mpi_world%rank == 0) then
    print*, "DD",  -CNST(2.0)*dd, -CNST(2.0)*mindd, -CNST(2.0)*maxdd
  end if

  if(present(deltaxc)) deltaxc = -CNST(2.0)*dd

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "fnxc", der%mesh, nxc, unit_one, ierr)

  call profiling_out(prof)

  SAFE_DEALLOCATE_A(lrvxc)
  SAFE_DEALLOCATE_A(nxc)

  POP_SUB('vxc_inc.xc_density_correction_calc')
end subroutine xc_density_correction_calc

! -----------------------------------------------------

FLOAT function get_qxc(mesh, nxc, density, ncutoff)  result(qxc)
  type(mesh_t), intent(in) :: mesh
  FLOAT,        intent(in) :: nxc(:)
  FLOAT,        intent(in) :: density(:)
  FLOAT,        intent(in) :: ncutoff

  integer :: ip
  FLOAT, allocatable :: nxc2(:)

  PUSH_SUB('vxc_inc.get_qxc')

  SAFE_ALLOCATE(nxc2(1:mesh%np))

  do ip = 1, mesh%np
    if(density(ip) < ncutoff) then
      nxc2(ip) = 0.0
    else
      nxc2(ip) = nxc(ip)
    end if
  end do

  qxc = dmf_integrate(mesh, nxc2)

  SAFE_DEALLOCATE_A(nxc2)

  POP_SUB('vxc_inc.get_qxc')
end function get_qxc


subroutine stitch_convex(mesh, get_branch, functionvalues, istart)

  interface 
     CMPLX function get_branch(x, branch)
       CMPLX :: x
       integer :: branch
     end function get_branch
  end interface

  type(mesh_t), intent(in) :: mesh
    
  CMPLX, intent(inout) :: functionvalues(:)
  integer, intent(in)  :: istart

  integer :: jstart(3), j, jj, nr, nl, N, jmin, jmax, tmp1, tmp2, i
  integer, allocatable :: idx(:), manyjmin(:), manyjmax(:)

  ! First stitch a single line along one direction
  call stitchline(mesh, get_branch, functionvalues, istart, 1, jmin, jmax)

  call index_to_coords(mesh%idx, mesh%sb%dim, istart, jstart)

  ! Now we stitch the next direction.  In case this is 3d, we will need
  ! to stitch along a third direction also; for this reason we will take
  ! care to remember the bounding indices
  SAFE_ALLOCATE(idx(1:1 + jmax - jmin))
  SAFE_ALLOCATE(manyjmin(jmin:jmax))
  SAFE_ALLOCATE(manyjmax(jmin:jmax))

  call mesh_subset_indices(mesh, (/jmin, jstart(2), 0/), (/jmax, jstart(2), 0/), idx)
  do j=jmin, jmax
     call stitchline(mesh, get_branch, functionvalues, idx(1 + j - jmin), 2, manyjmin(j), manyjmax(j))
  end do

  ! Now we have stitched one entire plane containing jstart.
  ! We can loop over all points in that plane and just stitch along the
  ! third direction now.
  if (mesh%sb%dim.eq.3) then
     do j=jmin, jmax
        do jj=manyjmin(j), manyjmax(j)
           i = index_from_coords(mesh%idx, mesh%sb%dim, (/j, jj, 0/))
           call stitchline(mesh, get_branch, functionvalues, i, 3, tmp1, tmp2)
        end do
     end do
  end if
  
  SAFE_DEALLOCATE_A(idx)
  SAFE_DEALLOCATE_A(manyjmin)
  SAFE_DEALLOCATE_A(manyjmax)

end subroutine stitch_convex

integer function countpoints(mesh, istart, displacement) result(npts)
  type(mesh_t), intent(in) :: mesh
  integer,      intent(in) :: istart
  integer,      intent(in) :: displacement(:)
  
  integer :: i
  
  i = istart
  npts = -1
  do while (i.ne.-1)
     npts = npts + 1
     i = translate_point(mesh, i, displacement)
  end do
end function countpoints
  
  subroutine stitchline(mesh, get_branch, functionvalues, istart, axis, jmin, jmax)
    type(mesh_t), intent(in) :: mesh

  interface 
     CMPLX function get_branch(x, branch)
       CMPLX :: x
       integer :: branch
     end function get_branch
  end interface

  CMPLX, intent(inout) :: functionvalues(:)
    integer, intent(in)      :: istart
    integer, intent(in)      :: axis ! x=1, y=2, z=3
    integer, intent(out)     :: jmin, jmax
    
    integer :: direction(3)
    integer :: otherdirection(3), nl, nr, vecjmax(3), vecjmin(3), vecjstart(3)
    integer, allocatable :: idx(:)
    
    vecjmin = 0
    vecjmax = 0
    vecjstart = 0
    
    direction = 0
    otherdirection = 0
    direction(axis) = 1
    otherdirection(axis) = -1

    nr = countpoints(mesh, istart, direction)
    nl = countpoints(mesh, istart, otherdirection)
    
    call index_to_coords(mesh%idx, mesh%sb%dim, istart, vecjstart)
    vecjmin = vecjstart
    vecjmax = vecjstart
    vecjmin(axis) = vecjmin(axis) - nl
    vecjmax(axis) = vecjmax(axis) + nr
    
    SAFE_ALLOCATE(idx(1:nr + nl + 1))
    call mesh_subset_indices(mesh, vecjmin, vecjmax, idx)
    call stitch(get_branch, functionvalues, nl + 1, idx)
    SAFE_DEALLOCATE_A(idx)
    jmin = vecjmin(axis)
    jmax = vecjmax(axis)
  end subroutine stitchline



! Subroutine to stitch discontinuous values of a multiple-valued function
! together to a single continuous, single-valued function by smoothly
! joining at the branch cuts.
subroutine stitch(get_branch, functionvalues, istart, idx)
  
  ! Function for getting values of multiple-valued functions.
  ! Each value of the parameter 'branch' corresponds to one such value.
  interface 
     CMPLX function get_branch(x, branch)
       CMPLX :: x
       integer :: branch
     end function get_branch
  end interface
  
  CMPLX, intent(inout) :: functionvalues(:)
  integer, intent(in) :: istart
  integer, intent(in) :: idx(:)
  
  integer :: currentbranch, npts, i
  CMPLX :: prev_value, err

  npts = size(idx, 1)

  ! First loop forwards from zero and stitch along the way
  currentbranch = 0
  prev_value = functionvalues(idx(istart))
  do i=istart + 1, npts
     call stitch_single_point()
  end do
  
  ! Now loop backwards
  currentbranch = 0
  prev_value = functionvalues(idx(istart))
  do i=istart - 1, 1, -1
     call stitch_single_point()
  end do
  
contains

  subroutine stitch_single_point()
    CMPLX :: v1, v2, v3, v
    integer :: j, adj
    v1 = get_branch(functionvalues(idx(i)), currentbranch)
    v2 = get_branch(functionvalues(idx(i)), currentbranch - 1)
    v3 = get_branch(functionvalues(idx(i)), currentbranch + 1)
    
    adj = 0
    v = v1
    if (abs(v2 - prev_value).lt.abs(v - prev_value)) then
       v = v2
       adj = -1
    end if
    if (abs(v3 - prev_value).lt.abs(v - prev_value)) then
       v = v3
       adj = +1
    end if
    currentbranch = currentbranch + adj
    functionvalues(idx(i)) = v
    prev_value = v
  end subroutine stitch_single_point
  
end subroutine stitch

! For evaluating values of multiple-valued functions when one value,
! e.g. the principal value, is known.  Used to stitch
CMPLX function get_root3_branch(x, branch) result(y)
  CMPLX, intent(in) :: x
  integer, intent(in) :: branch
  
  y = x * exp(branch * M_TWO * M_zI * M_PI / M_THREE)
end function get_root3_branch

CMPLX function get_root6_branch(x, branch) result(y)
  CMPLX, intent(in) :: x
  integer, intent(in) :: branch
  
  y = x * exp(branch * M_zI * M_PI / M_THREE)
end function get_root6_branch

CMPLX function get_logarithm_branch(x, branch) result(y)
  CMPLX, intent(in) :: x
  integer, intent(in) :: branch
  
  y = x + branch * M_TWO * M_zI * M_PI
end function get_logarithm_branch


subroutine zxc_complex_lda(mesh, rho, vxc, ex, ec, Imrho, Imvxc, Imex, Imec, cmplxscl_th)
  type(mesh_t), intent(in) :: mesh
  FLOAT, intent(in)        :: rho(:, :)
  FLOAT, intent(inout)     :: vxc(:, :)
  FLOAT, intent(inout)     :: ex
  FLOAT, intent(inout)     :: ec
  FLOAT, intent(in)        :: Imrho(:, :)
  FLOAT, intent(inout)     :: Imvxc(:, :)
  FLOAT, intent(inout)     :: Imex
  FLOAT, intent(inout)     :: Imec
  FLOAT, intent(in)        :: cmplxscl_th

  ! Exchange potential prefactor
  FLOAT, parameter :: Wx = -0.98474502184269641

  ! LDA correlation parameters
  FLOAT, parameter :: gamma = 0.031091, alpha1 = 0.21370, beta1 = 7.5957, beta2 = 3.5876, beta3 = 1.6382, beta4 = 0.49294

  CMPLX, allocatable   :: zrho_arr(:), zvx_arr(:), zvc_arr(:), rootrs(:), Q0(:), Q1(:), dQ1drs(:), epsc(:), depsdrs(:)
  integer, allocatable :: stitch_idx(:)
  CMPLX                :: dimphase, tmp, zex, zec
  FLOAT                :: scaling_origin(MAX_DIM)
  FLOAT                :: dmin_unused
  integer              :: N, izero, rankmin_unused, i

  
  N = size(rho, 1)
  SAFE_ALLOCATE(zrho_arr(1:N))
  SAFE_ALLOCATE(zvx_arr(1:N))
  SAFE_ALLOCATE(zvc_arr(1:N))
  SAFE_ALLOCATE(rootrs(1:N))
  SAFE_ALLOCATE(Q0(1:N))
  SAFE_ALLOCATE(Q1(1:N))
  SAFE_ALLOCATE(dQ1drs(1:N))
  SAFE_ALLOCATE(epsc(1:N))
  SAFE_ALLOCATE(depsdrs(1:N))
  SAFE_ALLOCATE(stitch_idx(1:N))

  do i=1, N
     stitch_idx(i) = i
  end do

  dimphase = exp(-mesh%sb%dim * M_zI * cmplxscl_th)
  scaling_origin = M_z0
  
  izero = mesh_nearest_point(mesh, scaling_origin, dmin_unused, rankmin_unused)
  
  zrho_arr(:) = rho(:, 1) + M_zI * Imrho(:, 1)

  ! Okay, now it starts.  First initialize the only point we know something
  ! definite about, namely the one where the scaling transformation does
  ! nothing.


  zvx_arr(:) = Wx * (zrho_arr(:) * dimphase)**(M_ONE / M_THREE)
  if(mesh%sb%dim.ne.1) then
     call stitch_convex(mesh, get_root3_branch, zvx_arr, izero)
  else
     call stitch(get_root3_branch, zvx_arr, izero, stitch_idx)
  end if

  zex = M_THREE / M_FOUR * sum(zvx_arr(:) * zrho_arr(:)) * mesh%volume_element

  ! Right.  Next is correlation which is much more complicated.  We
  ! use the PW91 parametrization.  We have to stitch the square root
  ! of the Wigner-Seitz radius and also the logarithm.
  
  rootrs(:) = (M_THREE / (M_FOUR * M_PI * zrho_arr(:) * dimphase))**(M_ONE / M_SIX)
  if(mesh%sb%dim.ne.1) then
     call stitch_convex(mesh, get_root6_branch, rootrs, izero)
  else
     call stitch(get_root6_branch, rootrs, izero, stitch_idx)
  end if
  Q0(:) = -M_TWO * gamma * (M_ONE + alpha1 * rootrs(:)**2)
  Q1(:) = M_TWO * gamma * rootrs(:) * (beta1 + rootrs(:) * (beta2 + rootrs(:) * (beta3 + rootrs(:) * beta4)))
  
  dQ1drs(:) = gamma * (beta1 / rootrs(:) + M_TWO * beta2 + rootrs(:) * (M_THREE * beta3 + M_FOUR * beta4 * rootrs(:)))

  epsc(:) = log(M_ONE + M_ONE / Q1)
  if(mesh%sb%dim.ne.1) then
     call stitch_convex(mesh, get_logarithm_branch, epsc, izero)
  else
     call stitch(get_logarithm_branch, epsc, izero, stitch_idx)
  end if
  epsc(:) = Q0(:) * epsc(:)

  depsdrs(:) = -M_TWO * gamma * alpha1 * epsc(:) / Q0(:) - Q0(:) * dQ1drs(:) / (Q1(:) * (Q1(:) + M_ONE))

  zvc_arr(:) = epsc(:) - rootrs(:)**2 * depsdrs(:) / M_THREE

  zec = sum(epsc(:) * zrho_arr(:)) * mesh%volume_element

  ex = real(zex)
  Imex = aimag(zex)
  ec = real(zec)
  Imec = aimag(zec)

  vxc(:, 1) = real(zvx_arr(:)) + real(zvc_arr(:))
  Imvxc(:, 1) = aimag(zvx_arr(:)) + aimag(zvc_arr(:))

  SAFE_DEALLOCATE_A(zrho_arr)
  SAFE_DEALLOCATE_A(zvx_arr)
  SAFE_DEALLOCATE_A(zvc_arr)
  SAFE_DEALLOCATE_A(rootrs)
  SAFE_DEALLOCATE_A(Q0)
  SAFE_DEALLOCATE_A(Q1)
  SAFE_DEALLOCATE_A(dQ1drs)
  SAFE_DEALLOCATE_A(epsc)
  SAFE_DEALLOCATE_A(depsdrs)
  SAFE_DEALLOCATE_A(stitch_idx)
end subroutine zxc_complex_lda


! ----------------------------------------------------------------------------- 
! This is the complex scaled interface for xc functionals.
! It will eventually be merged with the other one dxc_get_vxc after some test
! -----------------------------------------------------------------------------
subroutine zxc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, ex, ec, vxc, vtau, Imrho, Imex, Imec, Imvxc, Imvtau, cmplxscl_th)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)
  FLOAT,                intent(in)    :: Imrho(:, :)     !< cmplxscl: Electronic density 
  FLOAT, optional,      intent(inout) :: Imex            !< cmplxscl: Exchange energy.
  FLOAT, optional,      intent(inout) :: Imec            !< cmplxscl: Correlation energy
  FLOAT, optional,      intent(inout) :: Imvxc(:,:)      !< cmplxscl: XC potential
  FLOAT, optional,      intent(inout) :: Imvtau(:,:)     !< cmplxscl: Derivative wrt (two times kinetic energy density)
  FLOAT,                intent(in)    :: cmplxscl_th     !< complex scaling angle

  
  CMPLX, pointer :: zpot(:), zrho_tot(:)
  CMPLX          :: ztmp
  Integer        :: isp
  type(xc_functl_t), pointer :: functl(:)
  logical         :: calc_energy

  PUSH_SUB(zxc_get_vxc)

  print *, "LDA calc energy exc"
  ASSERT(present(ex) .eqv. present(ec))
  calc_energy = present(ex)

  !Pointer-shortcut for xcs%functl
  !It helps to remember that for xcs%functl(:,:)
  ! (1,:) => exchange,    (2,:) => correlation
  ! (:,1) => unpolarized, (:,2) => polarized
  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if


  
  if(functl(1)%id == XC_LDA_XC_CMPLX) then
    
    call zxc_complex_lda(der%mesh, rho, vxc, ex, ec, Imrho, Imvxc, Imex, Imec, cmplxscl_th)

    ! Exact exchange for 2 particles [vxc(r) = 1/2 * vh(r)]
    ! we keep it here for debug purposes
    if(.false.) then
      SAFE_ALLOCATE(zpot(1:size(vxc,1)))
      SAFE_ALLOCATE(zrho_tot(1:size(vxc,1)))

      zrho_tot = M_z0
      do isp = 1, ispin
        zrho_tot(:) = zrho_tot(:)+ rho(:,isp) +M_zI * Imrho(:,isp)
      end do

      call zpoisson_solve(psolver, zpot, zrho_tot, theta = cmplxscl_th)

      zpot = - zpot /CNST(2.0)
      vxc(:,1) = real(zpot(:)) 
      Imvxc(:,1) = aimag(zpot(:))
  
      if(calc_energy) then
        ztmp = M_HALF *zmf_dotp(der%mesh, zrho_tot, zpot, dotu = .true. )
        ex =  real(ztmp)
        Imex = aimag(ztmp)
        ec   = M_ZERO
        Imec = M_ZERO    
      end if
      SAFE_DEALLOCATE_P(zrho_tot)
      SAFE_DEALLOCATE_P(zpot)
    end if
    
  else if(functl(1)%family == XC_FAMILY_NONE) then

    vxc = M_ZERO
    if(calc_energy) then
      ex   = M_ZERO
      Imex = M_ZERO
      ec   = M_ZERO
      Imec = M_ZERO
    end if

  else  
    write(message(1), '(a)') 'The selected XCFunctional will not work with ComplexScaling = yes.'
    write(message(2), '(a)') 'Use XCFunctional = lda_xc_cmplx.'
    call messages_fatal(2)     
  end if

  POP_SUB(zxc_get_vxc)
end subroutine zxc_get_vxc


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
