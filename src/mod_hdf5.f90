module hdf5_output

   use sdata, only: dp
#ifdef HDF5_OUTPUT
   use hdf5
#endif

   implicit none

   save

   private

   public :: hdf5_configure, hdf5_is_active, hdf5_is_compiled
   public :: hdf5_write_step, hdf5_close

   logical :: hdf5_is_active = .false.
   logical :: h5_file_open = .false.
   logical :: h5_has_thermal = .false.
   logical :: h5_has_transient = .false.
   logical :: h5_has_control_rods = .false.
   logical :: h5_has_boron = .false.
   character(len = 200) :: h5_file_name = ''
   character(len = 200) :: h5_input_name = ''

#ifdef HDF5_OUTPUT
   integer(HID_T) :: h5_file_id
#endif

contains

   !******************************************************************************!

   logical function hdf5_is_compiled()

#ifdef HDF5_OUTPUT
      hdf5_is_compiled = .true.
#else
      hdf5_is_compiled = .false.
#endif

   end function hdf5_is_compiled

   !******************************************************************************!

   subroutine hdf5_configure(file_name, input_name, has_thermal, has_transient, &
      has_control_rods, has_boron)

      character(len = *), intent(in) :: file_name
      character(len = *), intent(in) :: input_name
      logical, intent(in) :: has_thermal
      logical, intent(in) :: has_transient
      logical, intent(in) :: has_control_rods
      logical, intent(in) :: has_boron

#ifdef HDF5_OUTPUT
      h5_file_name = trim(file_name)
      h5_input_name = trim(input_name)
      h5_has_thermal = has_thermal
      h5_has_transient = has_transient
      h5_has_control_rods = has_control_rods
      h5_has_boron = has_boron
      hdf5_is_active = .true.
#else
      hdf5_is_active = .false.
#endif

   end subroutine hdf5_configure

   !******************************************************************************!

   subroutine hdf5_write_step(step, time_s, reactivity, boron_concentration, power_w)

#ifdef HDF5_OUTPUT
      use sdata, only: ng, nnod, mode, ke, bcon, bpos, pos0, ssize, ctbeta, &
         f0, siga, sigf, nuf, vdel, pow, ppow, ftem, tfm, mtem, cden, nt
#endif

      integer, intent(in) :: step
      real(dp), intent(in) :: time_s
      real(dp), optional, intent(in) :: reactivity
      real(dp), optional, intent(in) :: boron_concentration
      real(dp), optional, intent(in) :: power_w

#ifdef HDF5_OUTPUT
      integer(HID_T) :: steps_id, step_id, global_id, node_id
      integer :: ierr
      character(len = 6) :: step_name
      character(len = 20) :: normalization
      real(dp) :: rho
      real(dp) :: current_power
      real(dp) :: raw_power
      real(dp) :: boron
      real(dp), dimension(1) :: scalar
      real(dp), dimension(:), allocatable :: relative_power
      real(dp), dimension(:), allocatable :: power_density
      real(dp), dimension(:), allocatable :: rod_position_cm
      real(dp), dimension(:), allocatable :: centerline
      real(dp), dimension(:), allocatable :: fuel_surface
      real(dp), dimension(:), allocatable :: clad_surface
      real(dp), dimension(:, :), allocatable :: flux
      real(dp), dimension(:, :), allocatable :: reaction_rate

      if (.not. hdf5_is_active) return

      if (.not. h5_file_open) call hdf5_open_file()

      write(step_name, '(I6.6)') step
      if (h5_has_thermal) then
         normalization = 'absolute'
      else
         normalization = 'relative'
      end if

      current_power = 0._dp
      if (h5_has_thermal) then
         if (present(power_w)) then
            current_power = power_w
         else
            current_power = pow * ppow * 0.01_dp
         end if
      end if

      allocate(relative_power(nnod))
      call get_relative_power(relative_power, raw_power)

      allocate(flux(nnod, ng))
      if (h5_has_thermal .and. raw_power > 0._dp) then
         flux = f0 * current_power / raw_power
      else
         flux = f0
      end if

      call h5gopen_f(h5_file_id, 'steps', steps_id, ierr)
      call h5gcreate_f(steps_id, trim(step_name), step_id, ierr)
      call h5gcreate_f(step_id, 'global', global_id, ierr)
      call h5gcreate_f(step_id, 'node', node_id, ierr)

      call write_real_scalar(global_id, 'time_s', time_s)
      call write_string(global_id, 'normalization', trim(normalization))
      if (h5_has_thermal) then
         call write_string(global_id, 'flux_units', 'n/cm2/s')
         call write_string(global_id, 'reaction_rate_units', '1/s')
         call write_string(global_id, 'power_density_units', 'kW/cm3')
      else
         call write_string(global_id, 'flux_units', 'relative')
         call write_string(global_id, 'reaction_rate_units', 'relative')
         call write_string(global_id, 'power_density_units', 'relative_per_cm3')
      end if

      if (present(reactivity)) then
         rho = reactivity
         call write_real_scalar(global_id, 'reactivity', rho)
         call write_real_scalar(global_id, 'reactivity_pcm', rho * 1.e5_dp)
         if (h5_has_transient .and. ctbeta > 0._dp) then
            call write_real_scalar(global_id, 'reactivity_dollar', rho / ctbeta)
         end if
      else if (mode /= 'FIXEDSRC' .and. ke /= 0._dp) then
         rho = (ke - 1._dp) / ke
         call write_real_scalar(global_id, 'reactivity', rho)
         call write_real_scalar(global_id, 'reactivity_pcm', rho * 1.e5_dp)
      end if

      if (h5_has_thermal) then
         call write_real_scalar(global_id, 'power_W', current_power)
      end if

      if (h5_has_boron) then
         if (present(boron_concentration)) then
            boron = boron_concentration
         else
            boron = bcon
         end if
         call write_real_scalar(global_id, 'boron_concentration_ppm', boron)
      end if

      if (h5_has_control_rods .and. allocated(bpos)) then
         call write_real_1d(global_id, 'control_rod_bank_position_step', bpos)
         allocate(rod_position_cm(size(bpos)))
         rod_position_cm = pos0 + bpos * ssize
         call write_real_1d(global_id, 'control_rod_bank_position_cm', rod_position_cm)
         deallocate(rod_position_cm)
      end if

      if (h5_has_thermal) then
         call write_real_scalar(global_id, 'max_fuel_temperature_K', maxval(ftem))
         call write_real_scalar(global_id, 'averaged_fuel_temperature_K', active_average(ftem))
         call write_real_scalar(global_id, 'max_fuel_centerline_temperature_K', maxval(tfm(:, 1)))
         call write_real_scalar(global_id, 'max_moderator_temperature_K', maxval(mtem))
         call write_real_scalar(global_id, 'averaged_moderator_temperature_K', active_average(mtem))
         call write_real_scalar(global_id, 'max_coolant_density_g_cm3', maxval(cden))
         call write_real_scalar(global_id, 'averaged_coolant_density_g_cm3', active_average(cden))
      end if

      call write_real_2d(node_id, 'flux', flux)

      allocate(power_density(nnod))
      if (h5_has_thermal) then
         power_density = relative_power * current_power * 1.e-3_dp / vdel
      else
         power_density = relative_power / vdel
      end if
      call write_real_1d(node_id, 'power_density', power_density)
      deallocate(power_density)

      allocate(reaction_rate(nnod, ng))
      reaction_rate = flux * nuf
      call write_real_2d(node_id, 'nu_fission_rate', reaction_rate)
      deallocate(reaction_rate)

      if (h5_has_thermal) then
         call write_real_1d(node_id, 'fuel_temperature_K', ftem)
         call write_real_1d(node_id, 'moderator_temperature_K', mtem)
         call write_real_1d(node_id, 'coolant_density_g_cm3', cden)

         allocate(centerline(nnod), fuel_surface(nnod), clad_surface(nnod))
         centerline = tfm(:, 1)
         fuel_surface = tfm(:, nt - 1)
         clad_surface = tfm(:, nt + 1)
         call write_real_1d(node_id, 'fuel_centerline_temperature_K', centerline)
         call write_real_1d(node_id, 'clad_surface_temperature_K', clad_surface)
         deallocate(centerline, fuel_surface, clad_surface)
      end if

      call h5gclose_f(node_id, ierr)
      call h5gclose_f(global_id, ierr)
      call h5gclose_f(step_id, ierr)
      call h5gclose_f(steps_id, ierr)

      deallocate(flux, relative_power)
#endif

   end subroutine hdf5_write_step

   !******************************************************************************!

   subroutine hdf5_close()

#ifdef HDF5_OUTPUT
      integer :: ierr

      if (h5_file_open) then
         call h5fclose_f(h5_file_id, ierr)
         call h5close_f(ierr)
         h5_file_open = .false.
      end if
#endif

   end subroutine hdf5_close

#ifdef HDF5_OUTPUT

   !******************************************************************************!

   subroutine hdf5_open_file()

      integer(HID_T) :: metadata_id, geometry_id, steps_id
      integer :: ierr

      call h5open_f(ierr)
      call h5fcreate_f(trim(h5_file_name), H5F_ACC_TRUNC_F, h5_file_id, ierr)
      h5_file_open = .true.

      call h5gcreate_f(h5_file_id, 'metadata', metadata_id, ierr)
      call write_metadata(metadata_id)
      call h5gclose_f(metadata_id, ierr)

      call h5gcreate_f(h5_file_id, 'geometry', geometry_id, ierr)
      call write_geometry(geometry_id)
      call h5gclose_f(geometry_id, ierr)

      call h5gcreate_f(h5_file_id, 'steps', steps_id, ierr)
      call h5gclose_f(steps_id, ierr)

   end subroutine hdf5_open_file

   !******************************************************************************!

   subroutine write_metadata(metadata_id)

      use sdata, only: mode, ng, nnod, nx, ny, nz, nxx, nyy, nzz, pow, ppow

      integer(HID_T), intent(in) :: metadata_id

      call write_int_scalar(metadata_id, 'schema_version', 1)
      call write_string(metadata_id, 'code_name', 'KOMODO')
      call write_string(metadata_id, 'input_file', trim(h5_input_name))
      call write_string(metadata_id, 'mode', trim(mode))
      call write_int_scalar(metadata_id, 'has_thermal_hydraulics', merge(1, 0, h5_has_thermal))
      call write_int_scalar(metadata_id, 'has_transient', merge(1, 0, h5_has_transient))
      call write_int_scalar(metadata_id, 'number_of_groups', ng)
      call write_int_scalar(metadata_id, 'number_of_nodes', nnod)
      call write_int_scalar(metadata_id, 'nx', nx)
      call write_int_scalar(metadata_id, 'ny', ny)
      call write_int_scalar(metadata_id, 'nz', nz)
      call write_int_scalar(metadata_id, 'nxx', nxx)
      call write_int_scalar(metadata_id, 'nyy', nyy)
      call write_int_scalar(metadata_id, 'nzz', nzz)

      if (h5_has_thermal) then
         call write_real_scalar(metadata_id, 'full_power_W', pow)
         call write_real_scalar(metadata_id, 'percent_power', ppow)
         call write_real_scalar(metadata_id, 'operating_power_W', pow * ppow * 0.01_dp)
      end if

   end subroutine write_metadata

   !******************************************************************************!

   subroutine write_geometry(geometry_id)

      use sdata, only: ix, iy, iz, mat, vdel, xdel, ydel, zdel, rpos, rdel, fbmap

      integer(HID_T), intent(in) :: geometry_id

      call write_int_1d(geometry_id, 'ix', ix)
      call write_int_1d(geometry_id, 'iy', iy)
      call write_int_1d(geometry_id, 'iz', iz)
      call write_int_1d(geometry_id, 'material', mat)
      call write_real_1d(geometry_id, 'node_volume_cm3', vdel)
      call write_real_1d(geometry_id, 'x_width_cm', xdel)
      call write_real_1d(geometry_id, 'y_width_cm', ydel)
      call write_real_1d(geometry_id, 'z_width_cm', zdel)

      if (h5_has_thermal) then
         if (allocated(rpos)) call write_real_1d(geometry_id, 'fuel_radial_mesh_position_m', rpos)
         if (allocated(rdel)) call write_real_1d(geometry_id, 'fuel_radial_mesh_delta_m', rdel)
      end if

      if (h5_has_control_rods) then
         if (allocated(fbmap)) call write_int_2d(geometry_id, 'control_rod_bank_map', fbmap)
      end if

   end subroutine write_geometry

   !******************************************************************************!

   subroutine get_relative_power(relative_power, raw_power)

      use sdata, only: ng, nnod, f0, sigf, vdel

      real(dp), dimension(:), intent(out) :: relative_power
      real(dp), intent(out) :: raw_power

      integer :: g, n
      real(dp) :: p

      relative_power = 0._dp
      do g = 1, ng
         do n = 1, nnod
            p = f0(n, g) * sigf(n, g) * vdel(n)
            if (p < 0._dp) p = 0._dp
            relative_power(n) = relative_power(n) + p
         end do
      end do

      raw_power = 0._dp
      do n = 1, nnod
         raw_power = raw_power + relative_power(n)
      end do

      if (raw_power > 0._dp) then
         do n = 1, nnod
            relative_power(n) = relative_power(n) / raw_power
         end do
      end if

   end subroutine get_relative_power

   !******************************************************************************!

   real(dp) function active_average(values)

      use sdata, only: nnod, ng, nuf, vdel

      real(dp), dimension(:), intent(in) :: values

      integer :: n
      real(dp) :: weighted_sum
      real(dp) :: volume_sum

      weighted_sum = 0._dp
      volume_sum = 0._dp
      do n = 1, nnod
         if (nuf(n, ng) > 0._dp) then
            weighted_sum = weighted_sum + values(n) * vdel(n)
            volume_sum = volume_sum + vdel(n)
         end if
      end do

      if (volume_sum > 0._dp) then
         active_average = weighted_sum / volume_sum
      else
         active_average = 0._dp
      end if

   end function active_average

   !******************************************************************************!

   subroutine write_real_scalar(loc_id, name, value)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      real(dp), intent(in) :: value

      real(dp), dimension(1) :: data

      data(1) = value
      call write_real_1d(loc_id, name, data)

   end subroutine write_real_scalar

   !******************************************************************************!

   subroutine write_int_scalar(loc_id, name, value)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      integer, intent(in) :: value

      integer, dimension(1) :: data

      data(1) = value
      call write_int_1d(loc_id, name, data)

   end subroutine write_int_scalar

   !******************************************************************************!

   subroutine write_real_1d(loc_id, name, data)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      real(dp), dimension(:), intent(in) :: data

      integer(HID_T) :: space_id, dataset_id
      integer(HSIZE_T), dimension(1) :: dims
      integer :: ierr

      dims(1) = size(data, 1)
      call h5screate_simple_f(1, dims, space_id, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dataset_id, ierr)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dims, ierr)
      call h5dclose_f(dataset_id, ierr)
      call h5sclose_f(space_id, ierr)

   end subroutine write_real_1d

   !******************************************************************************!

   subroutine write_real_2d(loc_id, name, data)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      real(dp), dimension(:, :), intent(in) :: data

      integer(HID_T) :: space_id, dataset_id
      integer(HSIZE_T), dimension(2) :: dims
      integer :: ierr
      real(dp), dimension(:, :), allocatable :: buffer

      allocate(buffer(size(data, 2), size(data, 1)))
      buffer = transpose(data)
      dims(1) = size(buffer, 1)
      dims(2) = size(buffer, 2)
      call h5screate_simple_f(2, dims, space_id, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dataset_id, ierr)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, buffer, dims, ierr)
      call h5dclose_f(dataset_id, ierr)
      call h5sclose_f(space_id, ierr)
      deallocate(buffer)

   end subroutine write_real_2d

   !******************************************************************************!

   subroutine write_int_1d(loc_id, name, data)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      integer, dimension(:), intent(in) :: data

      integer(HID_T) :: space_id, dataset_id
      integer(HSIZE_T), dimension(1) :: dims
      integer :: ierr

      dims(1) = size(data, 1)
      call h5screate_simple_f(1, dims, space_id, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, space_id, dataset_id, ierr)
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dims, ierr)
      call h5dclose_f(dataset_id, ierr)
      call h5sclose_f(space_id, ierr)

   end subroutine write_int_1d

   !******************************************************************************!

   subroutine write_int_2d(loc_id, name, data)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      integer, dimension(:, :), intent(in) :: data

      integer(HID_T) :: space_id, dataset_id
      integer(HSIZE_T), dimension(2) :: dims
      integer :: ierr
      integer, dimension(:, :), allocatable :: buffer

      allocate(buffer(size(data, 2), size(data, 1)))
      buffer = transpose(data)
      dims(1) = size(buffer, 1)
      dims(2) = size(buffer, 2)
      call h5screate_simple_f(2, dims, space_id, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, space_id, dataset_id, ierr)
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, buffer, dims, ierr)
      call h5dclose_f(dataset_id, ierr)
      call h5sclose_f(space_id, ierr)
      deallocate(buffer)

   end subroutine write_int_2d

   !******************************************************************************!

   subroutine write_string(loc_id, name, value)

      integer(HID_T), intent(in) :: loc_id
      character(len = *), intent(in) :: name
      character(len = *), intent(in) :: value

      integer(HID_T) :: type_id, space_id, dataset_id
      integer(HSIZE_T), dimension(1) :: dims
      integer(SIZE_T) :: string_size
      integer :: ierr
      character(len = 200) :: buffer

      buffer = ''
      buffer = trim(value)
      dims(1) = 1
      string_size = len(buffer)

      call h5tcopy_f(H5T_FORTRAN_S1, type_id, ierr)
      call h5tset_size_f(type_id, string_size, ierr)
      call h5screate_simple_f(1, dims, space_id, ierr)
      call h5dcreate_f(loc_id, trim(name), type_id, space_id, dataset_id, ierr)
      call h5dwrite_f(dataset_id, type_id, buffer, dims, ierr)
      call h5dclose_f(dataset_id, ierr)
      call h5sclose_f(space_id, ierr)
      call h5tclose_f(type_id, ierr)

   end subroutine write_string

#endif

end module hdf5_output
