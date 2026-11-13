module dlib

   use sdata, only: dp
   implicit none
   private

   public :: dlib_initialize, dlib_xs_updt, dlib_update_node_history
   public :: dlib_get_node_burnup, dlib_get_node_history, dlib_set_node_state
   public :: dlib_get_hm_density, dlib_get_material_name, dlib_is_fuel

   real(dp), parameter :: key_tolerance = 1.e-8_dp
   real(dp), parameter :: burnup_tolerance = 1.e-8_dp
   integer, parameter :: maximum_interaction_corners = 64
   character(len = 40), dimension(6), parameter :: family_names = (/ character(len = 40) :: &
      'CURRENT_BORON', 'CURRENT_MODERATOR_TEMPERATURE', 'CURRENT_FUEL_TEMPERATURE', &
      'HISTORY_BORON', 'HISTORY_MODERATOR_TEMPERATURE', 'HISTORY_FUEL_TEMPERATURE' /)

   type :: dlib_curve
      character(len = 40) :: kind = ''
      character(len = 64) :: label = ''
      real(dp), dimension(6) :: key = 0._dp
      integer :: nkey = 0
      real(dp), dimension(:), allocatable :: burnup
      real(dp), dimension(:, :), allocatable :: values
   end type dlib_curve

   type :: dlib_material
      integer :: composition_id = 0
      character(len = 64) :: name = ''
      character(len = 16) :: kind = ''
      real(dp) :: reference_fuel_temperature = 0._dp
      real(dp) :: reference_moderator_temperature = 0._dp
      real(dp) :: reference_boron = 0._dp
      real(dp) :: hm_smeared_density = 0._dp
      integer :: base_curve = 0
      logical :: has_reflector_interactions = .false.
      integer, dimension(6) :: family_point_count = 0
      integer, dimension(:, :), allocatable :: family_point_curve
      real(dp), dimension(:, :), allocatable :: family_point_key
      integer :: interaction_count = 0
      integer, dimension(6) :: interaction_axis_count = 0
      logical, dimension(6, 6) :: interaction_tied = .false.
      integer, dimension(:), allocatable :: interaction_curve
      integer, dimension(:, :), allocatable :: interaction_axis_index
      integer, dimension(:), allocatable :: interaction_order
      real(dp), dimension(:, :), allocatable :: interaction_axes
      real(dp), dimension(:), allocatable :: fission_spectrum
      character(len = 64), dimension(:), allocatable :: field_names
      integer :: curve_count = 0
      integer :: curve_capacity = 0
      type(dlib_curve), dimension(:), allocatable :: curves
   end type dlib_material

   type(dlib_material), dimension(:), allocatable, save :: library_materials
   integer, dimension(:), allocatable, save :: mat_to_library
   real(dp), dimension(:), allocatable, save :: node_burnup
   real(dp), dimension(:), allocatable, save :: node_history_boron
   real(dp), dimension(:), allocatable, save :: node_history_moderator_temperature
   real(dp), dimension(:), allocatable, save :: node_history_fuel_temperature
   integer, dimension(:), allocatable, save :: index_transport, index_siga
   integer, dimension(:), allocatable, save :: index_nusigf, index_kappa_fission
   integer, dimension(:), allocatable, save :: index_adf
   integer, dimension(:, :), allocatable, save :: index_scattering
   integer, save :: library_ng = 0, maximum_curve_count = 0
   integer, save :: ounit = 6
   logical, save :: library_loaded = .false.

contains

   subroutine dlib_initialize(output_unit)
      use sdata, only: ng, nmat, nnod, mat, chi, ccnuf, ccsigf, dlib_library, &
         dlib_composition_id, dlib_initial_burnup

      implicit none

      integer, intent(in) :: output_unit
      integer :: i, n, material_index
      real(dp) :: base_value

      if (library_loaded) return
      ounit = output_unit

      call read_library(trim(dlib_library))
      if (library_ng /= ng) call dlib_error('GROUPS changed between input parsing and library loading')
      if (.not. allocated(dlib_composition_id)) call dlib_error('DLIB material mapping is missing')
      if (.not. allocated(dlib_initial_burnup)) call dlib_error('DLIB initial burnup data are missing')
      if (size(dlib_composition_id) /= nmat .or. size(dlib_initial_burnup) /= nmat) then
         call dlib_error('material mapping size does not match NMAT')
      end if

      allocate(mat_to_library(nmat))
      do i = 1, nmat
         material_index = find_material(dlib_composition_id(i))
         if (material_index == 0) then
            call dlib_error('composition ID is not present in the depletion library')
         end if
         mat_to_library(i) = material_index
         if (library_materials(material_index)%kind == 'FUEL') then
            if (library_materials(material_index)%hm_smeared_density <= 0._dp) then
               call dlib_error('fuel material must have positive HM smeared density')
            end if
         end if
      end do

      if (allocated(chi)) then
         if (size(chi, 1) /= nmat .or. size(chi, 2) /= ng) then
            call dlib_error('fission spectrum array dimensions do not match DLIB')
         end if
      else
         allocate(chi(nmat, ng))
      end if
      ccnuf = .true.
      ccsigf = .true.
      do i = 1, nmat
         material_index = mat_to_library(i)
         chi(i, :) = library_materials(material_index)%fission_spectrum
         base_value = curve_value(material_index, library_materials(material_index)%base_curve, &
            index_nusigf(1), dlib_initial_burnup(i))
         if (any(library_materials(material_index)%fission_spectrum /= 0._dp) .and. &
            base_value > 0._dp) ccnuf = .false.
         base_value = curve_value(material_index, library_materials(material_index)%base_curve, &
            index_kappa_fission(1), dlib_initial_burnup(i))
         if (base_value > 0._dp) ccsigf = .false.
      end do
      if (ccnuf) call dlib_error('library contains no positive nu-fission cross section')
      if (ccsigf) call dlib_error('library contains no positive fission energy cross section')

      allocate(node_burnup(nnod), node_history_boron(nnod), &
         node_history_moderator_temperature(nnod), node_history_fuel_temperature(nnod))
      do n = 1, nnod
         if (mat(n) < 1 .or. mat(n) > nmat) call dlib_error('node material ID is outside the DLIB mapping')
         node_burnup(n) = dlib_initial_burnup(mat(n))
         material_index = mat_to_library(mat(n))
         ! Treat the specified initial inventory as having followed reference conditions.
         node_history_boron(n) = library_materials(material_index)%reference_boron
         node_history_moderator_temperature(n) = &
            library_materials(material_index)%reference_moderator_temperature
         node_history_fuel_temperature(n) = library_materials(material_index)%reference_fuel_temperature
      end do

      library_loaded = .true.
      write(ounit, *) ' FULL DLIB CURVES   : ', trim(dlib_library)
      write(ounit, *) ' LIBRARY MATERIALS  : ', size(library_materials)
      write(ounit, *) ' ...Depletion library curves are loaded...'

   end subroutine dlib_initialize

   !******************************************************************************!

   subroutine read_library(path)
      implicit none

      character(len = *), intent(in) :: path
      integer, parameter :: library_unit = 998
      integer :: ios, version, count_materials
      integer :: expected_materials
      logical :: found_signature, found_groups, found_end
      character(len = 4096) :: line
      character(len = 64) :: keyword
      type(dlib_material) :: material

      open(unit = library_unit, file = path, status = 'OLD', action = 'READ', iostat = ios)
      if (ios /= 0) call dlib_error('cannot open depletion library: ' // trim(path))

      found_signature = .false.
      found_groups = .false.
      found_end = .false.
      library_ng = 0
      maximum_curve_count = 0
      expected_materials = count_library_materials(path)
      if (expected_materials < 1) call dlib_error('library has no MATERIAL entries')
      allocate(library_materials(expected_materials))
      count_materials = 0

      do
         call next_library_line(library_unit, line, ios)
         if (ios < 0) exit
         if (ios /= 0) call dlib_error('error reading depletion library')

         if (.not. found_signature) then
            read(line, *, iostat = ios) keyword, version
            if (ios /= 0) call dlib_error('invalid depletion library signature')
            if (trim(keyword) /= 'KOMODO_DEPL' .or. version /= 1) then
               call dlib_error('unsupported library format; expected KOMODO_DEPL 1')
            end if
            found_signature = .true.
            cycle
         end if

         read(line, *, iostat = ios) keyword
         if (ios /= 0) cycle
         select case (trim(keyword))
          case ('GROUPS')
            if (found_groups) call dlib_error('duplicate GROUPS entry in library')
            read(line, *, iostat = ios) keyword, library_ng
            if (ios /= 0 .or. library_ng < 1) call dlib_error('invalid GROUPS entry in library')
            found_groups = .true.
          case ('MATERIAL')
            if (.not. found_groups) call dlib_error('GROUPS must precede MATERIAL entries')
            call read_material(library_unit, line, material)
            call append_material(material, count_materials)
            maximum_curve_count = max(maximum_curve_count, material%curve_count)
          case ('END_LIBRARY')
            found_end = .true.
            exit
          case default
            continue
         end select
      end do

      close(library_unit)
      if (.not. found_signature) call dlib_error('depletion library is empty')
      if (.not. found_groups) call dlib_error('GROUPS entry is missing from library')
      if (.not. found_end) call dlib_error('END_LIBRARY entry is missing from library')
      if (count_materials < 1) call dlib_error('library has no MATERIAL entries')
      if (count_materials /= expected_materials) then
         call dlib_error('MATERIAL count changed while reading the depletion library')
      end if
      call validate_library_fields()

   end subroutine read_library

   !******************************************************************************!

   integer function count_library_materials(path)
      implicit none

      character(len = *), intent(in) :: path
      integer, parameter :: count_unit = 995
      integer :: ios
      character(len = 4096) :: line
      character(len = 64) :: keyword

      count_library_materials = 0
      open(unit = count_unit, file = path, status = 'old', action = 'read', iostat = ios)
      if (ios /= 0) call dlib_error('cannot open depletion library: ' // trim(path))
      do
         call next_library_line(count_unit, line, ios)
         if (ios < 0) exit
         if (ios /= 0) call dlib_error('error reading depletion library')
         keyword = ''
         read(line, *, iostat = ios) keyword
         if (ios == 0 .and. trim(keyword) == 'MATERIAL') then
            count_library_materials = count_library_materials + 1
         end if
      end do
      close(count_unit)
   end function count_library_materials

   !******************************************************************************!

   subroutine read_material(unit, header, material)
      implicit none

      integer, intent(in) :: unit
      character(len = *), intent(in) :: header
      type(dlib_material), intent(out) :: material

      integer :: ios, nfield, i, j, k
      logical :: found_fields, found_spectrum, found_hm, found_ref_fuel, found_ref_mod, found_ref_boron
      character(len = 4096) :: line
      character(len = 64) :: keyword, token

      read(header, *, iostat = ios) keyword, material%composition_id, material%name, material%kind
      if (ios /= 0) call dlib_error('invalid MATERIAL entry in library')
      if (material%composition_id < 1) call dlib_error('library composition IDs must be positive')
      if (material%kind /= 'FUEL' .and. material%kind /= 'REFLECTOR') then
         call dlib_error('MATERIAL kind must be FUEL or REFLECTOR')
      end if

      found_fields = .false.
      found_spectrum = .false.
      found_hm = .false.
      found_ref_fuel = .false.
      found_ref_mod = .false.
      found_ref_boron = .false.
      nfield = 0

      do
         call next_library_line(unit, line, ios)
         if (ios < 0) call dlib_error('unexpected end of file inside MATERIAL')
         if (ios /= 0) call dlib_error('error reading MATERIAL data')
         read(line, *, iostat = ios) keyword
         if (ios /= 0) cycle

         select case (trim(keyword))
          case ('END_MATERIAL')
            exit
          case ('FIELDS')
            if (found_fields) call dlib_error('duplicate FIELDS entry in MATERIAL')
            read(line, *, iostat = ios) keyword, nfield
            if (ios /= 0 .or. nfield < 1) call dlib_error('invalid FIELDS entry')
            allocate(material%field_names(nfield))
            read(line, *, iostat = ios) keyword, nfield, (material%field_names(i), i = 1, nfield)
            if (ios /= 0) call dlib_error('field count does not match FIELDS entry')
            found_fields = .true.
          case ('FISSION_SPECTRUM')
            if (found_spectrum) call dlib_error('duplicate FISSION_SPECTRUM entry in MATERIAL')
            allocate(material%fission_spectrum(library_ng))
            read(line, *, iostat = ios) keyword, (material%fission_spectrum(i), i = 1, library_ng)
            if (ios /= 0) call dlib_error('invalid FISSION_SPECTRUM entry')
            found_spectrum = .true.
          case ('HM_SMEARED_DENSITY_KG_M3')
            if (found_hm) call dlib_error('duplicate HM_SMEARED_DENSITY_KG_M3 entry in MATERIAL')
            token = trim(adjustl(line(len_trim(keyword) + 1:)))
            if (trim(token) == 'N/A') then
               material%hm_smeared_density = 0._dp
            else
               read(token, *, iostat = ios) material%hm_smeared_density
               if (ios /= 0) call dlib_error('invalid HM_SMEARED_DENSITY_KG_M3 value')
            end if
            found_hm = .true.
          case ('REFERENCE_FUEL_TEMPERATURE_K')
            read(line, *, iostat = ios) keyword, material%reference_fuel_temperature
            if (ios /= 0) call dlib_error('invalid reference fuel temperature')
            found_ref_fuel = .true.
          case ('REFERENCE_MODERATOR_TEMPERATURE_K')
            read(line, *, iostat = ios) keyword, material%reference_moderator_temperature
            if (ios /= 0) call dlib_error('invalid reference moderator temperature')
            found_ref_mod = .true.
          case ('REFERENCE_BORON_PPM')
            read(line, *, iostat = ios) keyword, material%reference_boron
            if (ios /= 0) call dlib_error('invalid reference boron concentration')
            found_ref_boron = .true.
          case ('CURVE')
            if (.not. found_fields) call dlib_error('FIELDS must precede CURVE entries')
            call read_curve(unit, line, material)
          case ('SOURCE_CAX', 'BURNUP_UNIT', 'MACRO_XS_UNIT', 'KAPPA_FISSION_UNIT', &
             'ADF_UNIT', 'BETA_UNIT', 'LAMBDA_UNIT', 'INV_VELOCITY_UNIT', &
             'TEMPERATURE_UNIT', 'BORON_UNIT', 'HM_SMEARED_DENSITY_UNIT')
            continue
          case default
            call dlib_error('unrecognized MATERIAL entry: ' // trim(keyword))
         end select
      end do

      if (.not. found_fields) call dlib_error('FIELDS entry is missing from MATERIAL')
      if (.not. found_spectrum) call dlib_error('FISSION_SPECTRUM entry is missing from MATERIAL')
      if (.not. found_hm) call dlib_error('HM_SMEARED_DENSITY_KG_M3 entry is missing from MATERIAL')
      if (.not. found_ref_fuel .or. .not. found_ref_mod .or. .not. found_ref_boron) then
         call dlib_error('reference operating conditions are incomplete in MATERIAL')
      end if
      if (material%kind == 'FUEL' .and. material%hm_smeared_density <= 0._dp) then
         call dlib_error('FUEL MATERIAL requires positive HM smeared density')
      end if
      if (.not. ieee_finite(material%hm_smeared_density)) then
         call dlib_error('HM smeared density must be finite')
      end if
      if (.not. ieee_finite(material%reference_fuel_temperature) .or. &
         .not. ieee_finite(material%reference_moderator_temperature) .or. &
         .not. ieee_finite(material%reference_boron)) then
         call dlib_error('reference operating conditions must be finite')
      end if
      if (material%curve_count < 1) call dlib_error('MATERIAL has no curves')

      material%base_curve = 0
      do i = 1, material%curve_count
         if (material%curves(i)%kind == 'BASE') then
            if (material%base_curve /= 0) call dlib_error('duplicate BASE curve in MATERIAL')
            material%base_curve = i
         end if
      end do
      if (material%base_curve == 0) call dlib_error('BASE curve is missing from MATERIAL')

      do i = 1, material%curve_count
         do j = i + 1, material%curve_count
            if (material%curves(i)%kind /= material%curves(j)%kind) cycle
            if (material%curves(i)%kind == 'ROD') then
               if (material%curves(i)%label == material%curves(j)%label) then
                  call dlib_error('duplicate ROD curve in MATERIAL')
               end if
            else if (material%curves(i)%nkey == material%curves(j)%nkey) then
               if (all(abs(material%curves(i)%key(1:material%curves(i)%nkey) - &
                  material%curves(j)%key(1:material%curves(j)%nkey)) <= key_tolerance)) then
                  call dlib_error('duplicate branch curve key in MATERIAL')
               end if
            end if
         end do
      end do

      call build_material_lookup(material)

   end subroutine read_material

   !******************************************************************************!

   subroutine build_material_lookup(material)
      implicit none

      type(dlib_material), intent(inout) :: material
      integer :: family, i, j, count_points, curve_index
      integer :: curve_index_value
      real(dp) :: reference_key, key_value
      logical :: found

      allocate(material%family_point_curve(material%curve_count + 1, 6))
      allocate(material%family_point_key(material%curve_count + 1, 6))
      material%family_point_count = 0
      do family = 1, 6
         select case (family)
          case (1, 4)
            reference_key = material%reference_boron
          case (2, 5)
            reference_key = material%reference_moderator_temperature
          case (3, 6)
            reference_key = material%reference_fuel_temperature
         end select

         count_points = 1
         material%family_point_curve(1, family) = 0
         material%family_point_key(1, family) = reference_key
         do i = 1, material%curve_count
            if (material%curves(i)%kind /= family_names(family)) cycle
            found = .false.
            do j = 1, count_points
               if (abs(material%curves(i)%key(1) - material%family_point_key(j, family)) <= &
                  key_tolerance) then
                  found = .true.
                  exit
               end if
            end do
            if (found) then
               if (material%family_point_curve(j, family) /= 0) then
                  call dlib_error('duplicate branch curve key in MATERIAL')
               end if
               material%family_point_curve(j, family) = i
               material%family_point_key(j, family) = material%curves(i)%key(1)
            else
               count_points = count_points + 1
               material%family_point_curve(count_points, family) = i
               material%family_point_key(count_points, family) = material%curves(i)%key(1)
            end if
         end do

         material%family_point_count(family) = count_points
         do i = 2, count_points
            curve_index_value = material%family_point_curve(i, family)
            key_value = material%family_point_key(i, family)
            j = i - 1
            do while (j >= 1)
               if (material%family_point_key(j, family) <= key_value) exit
               material%family_point_curve(j + 1, family) = material%family_point_curve(j, family)
               material%family_point_key(j + 1, family) = material%family_point_key(j, family)
               j = j - 1
            end do
            material%family_point_curve(j + 1, family) = curve_index_value
            material%family_point_key(j + 1, family) = key_value
         end do
      end do

      call build_interaction_lookup(material)
      material%has_reflector_interactions = material%interaction_count > 0

   end subroutine build_material_lookup

   !******************************************************************************!

   subroutine build_interaction_lookup(material)
      implicit none

      type(dlib_material), intent(inout) :: material

      integer :: i, d, axis_index, count, new_index
      integer, dimension(6) :: temporary_count
      real(dp), dimension(:, :), allocatable :: temporary_axes
      real(dp) :: distance, best_distance

      material%interaction_count = 0
      do i = 1, material%curve_count
         if (material%curves(i)%kind == 'REFLECTOR_INTERACTION') then
            material%interaction_count = material%interaction_count + 1
         end if
      end do
      if (material%interaction_count == 0) return

      allocate(material%interaction_curve(material%interaction_count))
      allocate(material%interaction_axis_index(6, material%interaction_count))
      allocate(material%interaction_order(material%interaction_count))
      allocate(temporary_axes(material%interaction_count, 6))
      temporary_count = 0
      count = 0
      do i = 1, material%curve_count
         if (material%curves(i)%kind /= 'REFLECTOR_INTERACTION') cycle
         count = count + 1
         material%interaction_curve(count) = i
         do d = 1, 6
            axis_index = 0
            do new_index = 1, temporary_count(d)
               if (abs(temporary_axes(new_index, d) - material%curves(i)%key(d)) <= key_tolerance) then
                  axis_index = new_index
                  exit
               end if
            end do
            if (axis_index == 0) then
               temporary_count(d) = temporary_count(d) + 1
               axis_index = temporary_count(d)
               temporary_axes(axis_index, d) = material%curves(i)%key(d)
            end if
            material%interaction_axis_index(d, count) = axis_index
         end do
      end do

      do d = 1, 6
         call sort_axis(temporary_axes(:, d), temporary_count(d))
      end do
      material%interaction_axis_count = temporary_count
      allocate(material%interaction_axes(maxval(temporary_count), 6))
      material%interaction_axes = 0._dp
      do d = 1, 6
         material%interaction_axes(1:temporary_count(d), d) = temporary_axes(1:temporary_count(d), d)
      end do
      deallocate(temporary_axes)

      ! Map each curve's original key to the sorted axis index. Choose the first
      ! axis within tolerance, matching the runtime bracketing convention.
      do i = 1, material%interaction_count
         do d = 1, 6
            axis_index = 0
            best_distance = huge(best_distance)
            do new_index = 1, temporary_count(d)
               distance = abs(material%interaction_axes(new_index, d) - &
                  material%curves(material%interaction_curve(i))%key(d))
               if (distance <= key_tolerance .and. distance < best_distance) then
                  axis_index = new_index
                  best_distance = distance
               end if
            end do
            if (axis_index == 0) call dlib_error('could not index REFLECTOR_INTERACTION curve key')
            material%interaction_axis_index(d, i) = axis_index
         end do
      end do

      material%interaction_tied = .false.
      do d = 1, 6
         do new_index = 1, 6
            material%interaction_tied(d, new_index) = .true.
            do i = 1, material%interaction_count
               count = material%interaction_curve(i)
               if (abs(material%curves(count)%key(d) - material%curves(count)%key(new_index)) > &
                  key_tolerance) then
                  material%interaction_tied(d, new_index) = .false.
                  exit
               end if
            end do
         end do
      end do

      do i = 1, material%interaction_count
         material%interaction_order(i) = i
      end do
      call sort_interaction_order(material, material%interaction_order, 1, material%interaction_count)
   end subroutine build_interaction_lookup

   !******************************************************************************!

   subroutine read_curve(unit, header, material)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none

      integer, intent(in) :: unit
      character(len = *), intent(in) :: header
      type(dlib_material), intent(inout) :: material

      integer :: ios, npoints, i, f
      type(dlib_curve) :: curve
      character(len = 4096) :: line
      character(len = 64) :: keyword, mode

      mode = ''
      read(header, *, iostat = ios) keyword, curve%kind
      if (ios /= 0) call dlib_error('invalid CURVE entry')
      select case (trim(curve%kind))
       case ('BASE')
         read(header, *, iostat = ios) keyword, curve%kind, mode
         if (ios /= 0 .or. trim(mode) /= 'ABSOLUTE') call dlib_error('BASE curve must be ABSOLUTE')
       case ('CURRENT_BORON', 'CURRENT_MODERATOR_TEMPERATURE', 'CURRENT_FUEL_TEMPERATURE', &
          'HISTORY_BORON', 'HISTORY_MODERATOR_TEMPERATURE', 'HISTORY_FUEL_TEMPERATURE')
         curve%nkey = 1
         mode = ''
         read(header, *, iostat = ios) keyword, curve%kind, curve%key(1), mode
         if (ios /= 0 .or. trim(mode) /= 'DELTA') call dlib_error('branch curves must be DELTA curves')
       case ('ROD')
         mode = ''
         read(header, *, iostat = ios) keyword, curve%kind, curve%label, mode
         if (ios /= 0 .or. trim(mode) /= 'DELTA') call dlib_error('ROD curves must be DELTA curves')
       case ('REFLECTOR_INTERACTION')
         mode = ''
         read(header, *, iostat = ios) keyword, curve%kind, curve%label, &
            (curve%key(i), i = 1, 6), mode
         if (ios /= 0 .or. trim(mode) /= 'DELTA') then
            call dlib_error('REFLECTOR_INTERACTION curves must have six keys and be DELTA curves')
         end if
         if (trim(curve%label) /= 'TFU_HTFU_TMO_HTMO_BOR_HBOR') then
            call dlib_error('unsupported REFLECTOR_INTERACTION key descriptor')
         end if
         curve%nkey = 6
       case default
         call dlib_error('unsupported curve type: ' // trim(curve%kind))
      end select
      do i = 1, curve%nkey
         if (.not. ieee_is_finite(curve%key(i))) call dlib_error('curve keys must be finite')
      end do

      call next_library_line(unit, line, ios)
      if (ios /= 0) call dlib_error('missing POINTS entry after CURVE')
      npoints = 0
      keyword = ''
      read(line, *, iostat = ios) keyword, npoints
      if (ios /= 0 .or. trim(keyword) /= 'POINTS' .or. npoints < 1) then
         call dlib_error('invalid POINTS entry')
      end if

      allocate(curve%burnup(npoints), curve%values(size(material%field_names), npoints))
      do i = 1, npoints
         call next_library_line(unit, line, ios)
         if (ios /= 0) call dlib_error('unexpected end of file inside CURVE')
         read(line, *, iostat = ios) curve%burnup(i), &
            (curve%values(f, i), f = 1, size(material%field_names))
         if (ios /= 0) call dlib_error('invalid CURVE point or field count')
         if (.not. ieee_is_finite(curve%burnup(i))) call dlib_error('curve burnup must be finite')
         if (curve%burnup(i) < 0._dp) call dlib_error('curve burnup must be nonnegative')
         do f = 1, size(material%field_names)
            if (.not. ieee_is_finite(curve%values(f, i))) then
               call dlib_error('CURVE contains a non-finite field value')
            end if
         end do
         if (i > 1) then
            if (curve%burnup(i) <= curve%burnup(i - 1)) then
               call dlib_error('curve burnup points must be strictly increasing')
            end if
         end if
      end do

      call next_library_line(unit, line, ios)
      if (ios /= 0) call dlib_error('missing END_CURVE entry')
      read(line, *, iostat = ios) keyword
      if (ios /= 0 .or. trim(keyword) /= 'END_CURVE') call dlib_error('missing END_CURVE entry')
      call append_curve(material, curve)

   end subroutine read_curve

   !******************************************************************************!

   subroutine append_curve(material, curve)
      implicit none

      type(dlib_material), intent(inout) :: material
      type(dlib_curve), intent(in) :: curve
      type(dlib_curve), dimension(:), allocatable :: expanded
      integer :: new_capacity

      if (.not. allocated(material%curves)) then
         material%curve_capacity = 4
         allocate(material%curves(material%curve_capacity))
      else if (material%curve_count == material%curve_capacity) then
         new_capacity = max(4, 2 * material%curve_capacity)
         allocate(expanded(new_capacity))
         expanded(1:material%curve_count) = material%curves(1:material%curve_count)
         call move_alloc(expanded, material%curves)
         material%curve_capacity = new_capacity
      end if
      material%curve_count = material%curve_count + 1
      material%curves(material%curve_count) = curve

   end subroutine append_curve

   !******************************************************************************!

   subroutine append_material(material, count)
      implicit none

      type(dlib_material), intent(in) :: material
      integer, intent(inout) :: count
      integer :: i

      do i = 1, count
         if (library_materials(i)%composition_id == material%composition_id) then
            call dlib_error('duplicate MATERIAL composition ID in library')
         end if
      end do
      if (count >= size(library_materials)) call dlib_error('MATERIAL count exceeds first-pass count')
      count = count + 1
      library_materials(count) = material

   end subroutine append_material

   !******************************************************************************!

   subroutine next_library_line(unit, line, ios)
      implicit none

      integer, intent(in) :: unit
      character(len = *), intent(out) :: line
      integer, intent(out) :: ios
      integer :: comment_position

      do
         read(unit, '(A)', iostat = ios) line
         if (ios /= 0) return
         comment_position = index(line, '#')
         if (comment_position > 0) line = line(1:comment_position - 1)
         line = trim(adjustl(line))
         if (len_trim(line) > 0) return
      end do

   end subroutine next_library_line

   !******************************************************************************!

   subroutine validate_library_fields()
      implicit none

      integer :: m, g, h, i
      character(len = 64) :: field

      do m = 2, size(library_materials)
         if (size(library_materials(m)%field_names) /= size(library_materials(1)%field_names)) then
            call dlib_error('FIELD count differs between MATERIAL entries')
         end if
         do i = 1, size(library_materials(1)%field_names)
            if (library_materials(m)%field_names(i) /= library_materials(1)%field_names(i)) then
               call dlib_error('FIELD ordering differs between MATERIAL entries')
            end if
         end do
      end do

      do m = 1, size(library_materials)
         do i = 1, size(library_materials(m)%field_names)
            do h = i + 1, size(library_materials(m)%field_names)
               if (library_materials(m)%field_names(i) == library_materials(m)%field_names(h)) then
                  call dlib_error('duplicate FIELD name in MATERIAL')
               end if
            end do
         end do
      end do

      allocate(index_transport(library_ng), index_siga(library_ng))
      allocate(index_nusigf(library_ng), index_kappa_fission(library_ng), index_adf(library_ng))
      allocate(index_scattering(library_ng, library_ng))
      do g = 1, library_ng
         write(field, '(A,I0)') 'transport_', g
         index_transport(g) = field_location(1, trim(field))
         write(field, '(A,I0)') 'siga_', g
         index_siga(g) = field_location(1, trim(field))
         write(field, '(A,I0)') 'nusigf_', g
         index_nusigf(g) = field_location(1, trim(field))
         write(field, '(A,I0)') 'kappa_fission_', g
         index_kappa_fission(g) = field_location(1, trim(field))
         write(field, '(A,I0)') 'adf_', g
         index_adf(g) = field_location(1, trim(field))
         do h = 1, library_ng
            write(field, '(A,I0,I0)') 'sigs_', g, h
            index_scattering(g, h) = field_location(1, trim(field))
         end do
      end do

      do m = 1, size(library_materials)
         if (.not. allocated(library_materials(m)%fission_spectrum)) then
            call dlib_error('FISSION_SPECTRUM is missing from MATERIAL')
         end if
         if (size(library_materials(m)%fission_spectrum) /= library_ng) then
            call dlib_error('FISSION_SPECTRUM group count does not match GROUPS')
         end if
         do g = 1, library_ng
            if (.not. ieee_finite(library_materials(m)%fission_spectrum(g))) then
               call dlib_error('FISSION_SPECTRUM contains a non-finite value')
            end if
            if (library_materials(m)%fission_spectrum(g) < 0._dp) then
               call dlib_error('FISSION_SPECTRUM cannot contain negative values')
            end if
         end do
         if (library_materials(m)%kind == 'FUEL' .and. &
            sum(library_materials(m)%fission_spectrum) <= 0._dp) then
            call dlib_error('FUEL MATERIAL requires a positive fission spectrum')
         end if
      end do

   end subroutine validate_library_fields

   !******************************************************************************!

   integer function field_location(material_index, name)
      implicit none

      integer, intent(in) :: material_index
      character(len = *), intent(in) :: name
      integer :: i

      field_location = 0
      do i = 1, size(library_materials(material_index)%field_names)
         if (trim(library_materials(material_index)%field_names(i)) == trim(name)) then
            field_location = i
            return
         end if
      end do
      call dlib_error('required FIELD is missing: ' // trim(name))

   end function field_location

   !******************************************************************************!

   integer function find_material(composition_id)
      implicit none

      integer, intent(in) :: composition_id
      integer :: i

      find_material = 0
      do i = 1, size(library_materials)
         if (library_materials(i)%composition_id == composition_id) then
            find_material = i
            return
         end if
      end do

   end function find_material

   !******************************************************************************!

   subroutine dlib_xs_updt(xbcon, xftem, xmtem, xbpos)
      use sdata, only: nnod, ng, mat, sigtr, siga, nuf, sigf, sigs, dc, nb

      implicit none

      real(dp), intent(in) :: xbcon
      real(dp), dimension(:), intent(in) :: xftem, xmtem, xbpos

      integer :: n, g, h, family, material_index, selected_interactions, field_count
      integer, dimension(maximum_interaction_corners) :: interaction_indices
      integer, dimension(2, 6) :: family_indices
      real(dp), dimension(maximum_interaction_corners) :: interaction_weights
      real(dp), dimension(6) :: interaction_key
      real(dp), dimension(2, 6) :: family_weights
      real(dp), dimension(6) :: family_query
      real(dp), dimension(:), allocatable :: field_values
      real(dp), dimension(:, :), allocatable :: curve_values
      logical, dimension(:), allocatable :: curve_ready
      real(dp) :: current_fuel_temperature

      if (.not. library_loaded) call dlib_error('library data were not initialized before XS update')
      if (size(xftem) /= nnod .or. size(xmtem) /= nnod) then
         call dlib_error('temperature arrays do not match the node count')
      end if
      if (.not. ieee_finite(xbcon)) call dlib_error('current boron concentration is not finite')

      field_count = size(library_materials(1)%field_names)
      allocate(field_values(field_count), curve_values(field_count, maximum_curve_count))
      allocate(curve_ready(maximum_curve_count))

      do n = 1, nnod
         material_index = mat_to_library(mat(n))
         if (.not. ieee_finite(xftem(n)) .or. .not. ieee_finite(xmtem(n))) then
            call dlib_error('current node temperature is not finite')
         end if
         if (.not. ieee_finite(node_history_boron(n)) .or. &
            .not. ieee_finite(node_history_moderator_temperature(n)) .or. &
            .not. ieee_finite(node_history_fuel_temperature(n))) then
            call dlib_error('history condition is not finite')
         end if
         current_fuel_temperature = xftem(n)
         if (library_materials(material_index)%kind == 'REFLECTOR') then
            ! Thermal feedback supplies a coolant-derived ftem even in non-fuel nodes.
            current_fuel_temperature = library_materials(material_index)%reference_fuel_temperature
         end if
         interaction_key = (/current_fuel_temperature, node_history_fuel_temperature(n), &
            xmtem(n), node_history_moderator_temperature(n), xbcon, node_history_boron(n)/)
         selected_interactions = 0
         if (library_materials(material_index)%has_reflector_interactions) then
            call interaction_weights_for_node(material_index, interaction_key, interaction_indices, &
               interaction_weights, selected_interactions)
         end if
         family_query = (/xbcon, xmtem(n), current_fuel_temperature, node_history_boron(n), &
            node_history_moderator_temperature(n), node_history_fuel_temperature(n)/)
         do family = 1, 6
            call family_weights_for_node(material_index, family, family_query(family), &
               family_indices(:, family), family_weights(:, family))
         end do
         curve_ready = .false.
         call evaluate_fields(material_index, node_burnup(n), family_indices, family_weights, &
            interaction_indices, interaction_weights, selected_interactions, curve_values, &
            curve_ready, field_values)
         do g = 1, ng
            sigtr(n, g) = field_values(index_transport(g))
            siga(n, g) = field_values(index_siga(g))
            nuf(n, g) = field_values(index_nusigf(g))
            sigf(n, g) = field_values(index_kappa_fission(g))
            dc(n, g, :) = field_values(index_adf(g))
            do h = 1, ng
               sigs(n, g, h) = field_values(index_scattering(g, h))
            end do
         end do
      end do

      deallocate(field_values, curve_values, curve_ready)
      if (nb > 0) call apply_dlib_rods(xbpos)

   end subroutine dlib_xs_updt

   !******************************************************************************!

   subroutine evaluate_fields(material_index, burnup, family_indices, family_weights, &
      interaction_indices, interaction_weights, ninteraction, &
      curve_cache, curve_ready, field_values)
      implicit none

      integer, intent(in) :: material_index, ninteraction
      real(dp), intent(in) :: burnup
      integer, dimension(2, 6), intent(in) :: family_indices
      real(dp), dimension(2, 6), intent(in) :: family_weights
      integer, dimension(:), intent(in) :: interaction_indices
      real(dp), dimension(:), intent(in) :: interaction_weights
      real(dp), dimension(:, :), intent(inout) :: curve_cache
      logical, dimension(:), intent(inout) :: curve_ready
      real(dp), dimension(:), intent(out) :: field_values

      integer :: family, branch, i, base_curve, curve_index
      real(dp), dimension(size(field_values)) :: correction

      base_curve = library_materials(material_index)%base_curve
      call cached_curve_values(material_index, base_curve, burnup, curve_cache, curve_ready)
      field_values = curve_cache(:, base_curve)

      do family = 1, 3
         correction = 0._dp
         do branch = 1, 2
            curve_index = family_indices(branch, family)
            if (curve_index == 0) cycle
            call cached_curve_values(material_index, curve_index, burnup, curve_cache, curve_ready)
            correction = correction + family_weights(branch, family) * curve_cache(:, curve_index)
         end do
         field_values = field_values + correction
      end do

      if (ninteraction == 0) then
         do family = 4, 6
            correction = 0._dp
            do branch = 1, 2
               curve_index = family_indices(branch, family)
               if (curve_index == 0) cycle
               call cached_curve_values(material_index, curve_index, burnup, curve_cache, curve_ready)
               correction = correction + family_weights(branch, family) * curve_cache(:, curve_index)
            end do
            field_values = field_values + correction
         end do
      else
         do i = 1, ninteraction
            curve_index = interaction_indices(i)
            call cached_curve_values(material_index, curve_index, burnup, curve_cache, curve_ready)
            field_values = field_values + interaction_weights(i) * curve_cache(:, curve_index)
         end do
      end if

   end subroutine evaluate_fields

   !******************************************************************************!

   subroutine cached_curve_values(material_index, curve_index, burnup, curve_cache, curve_ready)
      implicit none

      integer, intent(in) :: material_index, curve_index
      real(dp), intent(in) :: burnup
      real(dp), dimension(:, :), intent(inout) :: curve_cache
      logical, dimension(:), intent(inout) :: curve_ready

      if (curve_ready(curve_index)) return
      call curve_values_at_burnup(material_index, curve_index, burnup, curve_cache(:, curve_index))
      curve_ready(curve_index) = .true.

   end subroutine cached_curve_values

   !******************************************************************************!

   subroutine family_weights_for_node(material_index, family, query_key, indices, weights)
      implicit none

      integer, intent(in) :: material_index, family
      real(dp), intent(in) :: query_key
      integer, dimension(2), intent(out) :: indices
      real(dp), dimension(2), intent(out) :: weights

      integer :: i, count_points, low, high
      real(dp) :: denominator

      count_points = library_materials(material_index)%family_point_count(family)
      if (count_points < 1) call dlib_error('branch condition axis is empty')
      indices = 0
      weights = 0._dp
      do i = 1, count_points
         if (abs(query_key - library_materials(material_index)%family_point_key(i, family)) <= &
            key_tolerance) then
            indices(1) = library_materials(material_index)%family_point_curve(i, family)
            weights(1) = 1._dp
            return
         end if
      end do

      if (count_points == 1) then
         indices(1) = library_materials(material_index)%family_point_curve(1, family)
         weights(1) = 1._dp
         return
      end if

      if (query_key < library_materials(material_index)%family_point_key(1, family)) then
         low = 1
         high = 2
      else if (query_key > &
         library_materials(material_index)%family_point_key(count_points, family)) then
         low = count_points - 1
         high = count_points
      else
         low = 1
         high = 2
         do i = 2, count_points
            if (query_key < library_materials(material_index)%family_point_key(i, family)) then
               low = i - 1
               high = i
               exit
            end if
         end do
      end if

      indices(1) = library_materials(material_index)%family_point_curve(low, family)
      indices(2) = library_materials(material_index)%family_point_curve(high, family)
      denominator = library_materials(material_index)%family_point_key(high, family) - &
         library_materials(material_index)%family_point_key(low, family)
      weights(2) = (query_key - library_materials(material_index)%family_point_key(low, family)) / denominator
      weights(1) = 1._dp - weights(2)

   end subroutine family_weights_for_node

   !******************************************************************************!

   real(dp) function curve_value(material_index, curve_index, field_index, burnup)
      implicit none

      integer, intent(in) :: material_index, curve_index, field_index
      real(dp), intent(in) :: burnup
      integer :: npoints, lower, upper
      real(dp) :: fraction

      npoints = size(library_materials(material_index)%curves(curve_index)%burnup)
      call burnup_interval(material_index, curve_index, burnup, lower, fraction)
      if (npoints == 1) then
         curve_value = library_materials(material_index)%curves(curve_index)%values(field_index, 1)
         return
      end if

      upper = lower + 1
      curve_value = library_materials(material_index)%curves(curve_index)%values(field_index, lower) + &
         fraction * (library_materials(material_index)%curves(curve_index)%values(field_index, upper) - &
         library_materials(material_index)%curves(curve_index)%values(field_index, lower))

   end function curve_value

   !******************************************************************************!

   subroutine curve_values_at_burnup(material_index, curve_index, burnup, values)
      implicit none

      integer, intent(in) :: material_index, curve_index
      real(dp), intent(in) :: burnup
      real(dp), dimension(:), intent(out) :: values
      integer :: npoints, lower
      real(dp) :: fraction

      npoints = size(library_materials(material_index)%curves(curve_index)%burnup)
      call burnup_interval(material_index, curve_index, burnup, lower, fraction)
      if (npoints == 1) then
         values = library_materials(material_index)%curves(curve_index)%values(:, 1)
      else
         values = library_materials(material_index)%curves(curve_index)%values(:, lower) + fraction * &
            (library_materials(material_index)%curves(curve_index)%values(:, lower + 1) - &
            library_materials(material_index)%curves(curve_index)%values(:, lower))
      end if

   end subroutine curve_values_at_burnup

   !******************************************************************************!

   subroutine burnup_interval(material_index, curve_index, burnup, lower, fraction)
      implicit none

      integer, intent(in) :: material_index, curve_index
      real(dp), intent(in) :: burnup
      integer, intent(out) :: lower
      real(dp), intent(out) :: fraction

      integer :: npoints, i
      real(dp) :: query_burnup, minimum_burnup, maximum_burnup

      npoints = size(library_materials(material_index)%curves(curve_index)%burnup)
      minimum_burnup = library_materials(material_index)%curves(curve_index)%burnup(1)
      maximum_burnup = library_materials(material_index)%curves(curve_index)%burnup(npoints)
      if (.not. ieee_finite(burnup) .or. &
         burnup < minimum_burnup - burnup_tolerance .or. &
         burnup > maximum_burnup + burnup_tolerance) then
         call dlib_error('requested burnup is outside a DLIB curve range')
      end if
      query_burnup = max(minimum_burnup, min(burnup, maximum_burnup))
      lower = 1
      fraction = 0._dp
      if (npoints == 1) return

      if (query_burnup >= maximum_burnup) then
         lower = npoints - 1
      else
         do i = 2, npoints
            if (query_burnup <= library_materials(material_index)%curves(curve_index)%burnup(i)) then
               lower = i - 1
               exit
            end if
         end do
      end if
      fraction = (query_burnup - library_materials(material_index)%curves(curve_index)%burnup(lower)) / &
         (library_materials(material_index)%curves(curve_index)%burnup(lower + 1) - &
         library_materials(material_index)%curves(curve_index)%burnup(lower))

   end subroutine burnup_interval

   !******************************************************************************!

   subroutine interaction_weights_for_node(material_index, query, indices, weights, nselected)
      implicit none

      integer, intent(in) :: material_index
      real(dp), dimension(6), intent(in) :: query
      integer, dimension(:), intent(out) :: indices
      real(dp), dimension(:), intent(out) :: weights
      integer, intent(out) :: nselected

      integer :: i, j, d, nactive, mask, bit, curve_index, tied_dimension, axis_matches
      integer, dimension(6) :: active, group_rep, lower_index, upper_index, target_index
      real(dp), dimension(6) :: lower_key, upper_key, fraction
      real(dp) :: weight, axis_limit
      logical :: found, unique_axis_tuple

      nselected = 0
      unique_axis_tuple = .true.
      do d = 1, 6
         axis_matches = 0
         do i = 1, library_materials(material_index)%interaction_axis_count(d)
            axis_limit = 2._dp * key_tolerance + epsilon(1._dp) * &
               max(1._dp, abs(library_materials(material_index)%interaction_axes(i, d)), abs(query(d)))
            if (abs(library_materials(material_index)%interaction_axes(i, d) - query(d)) <= &
               axis_limit) then
               axis_matches = axis_matches + 1
               target_index(d) = i
            end if
         end do
         if (axis_matches /= 1) unique_axis_tuple = .false.
      end do
      if (unique_axis_tuple) then
         curve_index = find_interaction_corner(material_index, target_index)
         if (curve_index > 0) then
            found = all(abs(library_materials(material_index)%curves(curve_index)%key - query) <= &
               key_tolerance)
            if (found) then
               nselected = 1
               indices(1) = curve_index
               weights(1) = 1._dp
               return
            end if
         end if
      end if

      ! Retain the original file-order check when tolerance neighborhoods overlap
      ! or a clustered axis does not identify a unique exact state.
      do i = 1, library_materials(material_index)%interaction_count
         curve_index = library_materials(material_index)%interaction_curve(i)
         found = all(abs(library_materials(material_index)%curves(curve_index)%key - query) <= key_tolerance)
         if (found) then
            nselected = 1
            indices(1) = curve_index
            weights(1) = 1._dp
            return
         end if
      end do

      do d = 1, 6
         call bracket_axis(library_materials(material_index)%interaction_axes(:, d), &
            library_materials(material_index)%interaction_axis_count(d), query(d), .true., &
            lower_key(d), upper_key(d), fraction(d))
         lower_index(d) = interaction_axis_index(material_index, d, lower_key(d))
         upper_index(d) = interaction_axis_index(material_index, d, upper_key(d))
      end do

      nactive = 0
      group_rep = 0
      do d = 1, 6
         if (upper_key(d) > lower_key(d) + key_tolerance) then
            tied_dimension = 0
            do bit = 1, nactive
               j = active(bit)
               if (abs(query(d) - query(j)) > key_tolerance) cycle
               if (library_materials(material_index)%interaction_tied(d, j)) then
                  tied_dimension = j
                  exit
               end if
            end do
            if (tied_dimension > 0) then
               group_rep(d) = tied_dimension
            else
               nactive = nactive + 1
               active(nactive) = d
               group_rep(d) = d
            end if
         end if
      end do

      do mask = 0, 2**nactive - 1
         target_index = lower_index
         weight = 1._dp
         do bit = 1, nactive
            d = active(bit)
            if (btest(mask, bit - 1)) then
               do j = 1, 6
                  if (group_rep(j) == d) target_index(j) = upper_index(j)
               end do
               weight = weight * fraction(d)
            else
               weight = weight * (1._dp - fraction(d))
            end if
         end do
         if (abs(weight) <= tiny(weight)) cycle

         curve_index = find_interaction_corner(material_index, target_index)
         if (curve_index == 0) then
            call dlib_error('REFLECTOR_INTERACTION data do not form a complete interpolation cell')
         end if
         nselected = nselected + 1
         if (nselected > size(indices) .or. nselected > size(weights)) then
            call dlib_error('REFLECTOR_INTERACTION interpolation stencil exceeds 64 corners')
         end if
         indices(nselected) = curve_index
         weights(nselected) = weight
      end do

   end subroutine interaction_weights_for_node

   !******************************************************************************!

   integer function interaction_axis_index(material_index, dimension, value)
      implicit none
      integer, intent(in) :: material_index, dimension
      real(dp), intent(in) :: value
      integer :: low, high, middle, count_points

      interaction_axis_index = 0
      count_points = library_materials(material_index)%interaction_axis_count(dimension)
      low = 1
      high = count_points
      do while (low < high)
         middle = (low + high) / 2
         if (library_materials(material_index)%interaction_axes(middle, dimension) < value - key_tolerance) then
            low = middle + 1
         else
            high = middle
         end if
      end do
      if (low <= count_points) then
         if (abs(library_materials(material_index)%interaction_axes(low, dimension) - value) <= key_tolerance) then
            interaction_axis_index = low
            return
         end if
      end if
      call dlib_error('could not locate a REFLECTOR_INTERACTION axis point')
   end function interaction_axis_index

   !******************************************************************************!

   integer function find_interaction_corner(material_index, target_index)
      implicit none
      integer, intent(in) :: material_index
      integer, dimension(6), intent(in) :: target_index
      integer :: low, high, middle, record, curve_index, dimension
      real(dp), dimension(6) :: target_key
      logical :: matches

      do dimension = 1, 6
         target_key(dimension) = library_materials(material_index)%interaction_axes( &
            target_index(dimension), dimension)
      end do

      low = 1
      high = library_materials(material_index)%interaction_count
      do while (low <= high)
         middle = (low + high) / 2
         record = library_materials(material_index)%interaction_order(middle)
         if (all(library_materials(material_index)%interaction_axis_index(:, record) == target_index)) then
            curve_index = library_materials(material_index)%interaction_curve(record)
            matches = all(abs(library_materials(material_index)%curves(curve_index)%key - &
               target_key) <= key_tolerance)
            if (matches .and. .not. has_duplicate_interaction_tuple(material_index, middle)) then
               find_interaction_corner = curve_index
               return
            end if
            exit
         else if (integer_tuple_less(library_materials(material_index)%interaction_axis_index(:, record), &
            target_index)) then
            low = middle + 1
         else
            high = middle - 1
         end if
      end do

      ! Tolerance-clustered axis values can place a curve in an adjacent tuple.
      ! Keep a correctness-preserving fallback for those uncommon boundary cases.
      do record = 1, library_materials(material_index)%interaction_count
         curve_index = library_materials(material_index)%interaction_curve(record)
         if (all(abs(library_materials(material_index)%curves(curve_index)%key - target_key) <= key_tolerance)) then
            find_interaction_corner = curve_index
            return
         end if
      end do
      find_interaction_corner = 0
   end function find_interaction_corner

   !******************************************************************************!

   logical function has_duplicate_interaction_tuple(material_index, order_position)
      implicit none
      integer, intent(in) :: material_index, order_position
      integer :: record, previous_record, next_record

      record = library_materials(material_index)%interaction_order(order_position)
      has_duplicate_interaction_tuple = .false.
      if (order_position > 1) then
         previous_record = library_materials(material_index)%interaction_order(order_position - 1)
         if (all(library_materials(material_index)%interaction_axis_index(:, record) == &
            library_materials(material_index)%interaction_axis_index(:, previous_record))) then
            has_duplicate_interaction_tuple = .true.
            return
         end if
      end if
      if (order_position < library_materials(material_index)%interaction_count) then
         next_record = library_materials(material_index)%interaction_order(order_position + 1)
         if (all(library_materials(material_index)%interaction_axis_index(:, record) == &
            library_materials(material_index)%interaction_axis_index(:, next_record))) then
            has_duplicate_interaction_tuple = .true.
         end if
      end if
   end function has_duplicate_interaction_tuple

   !******************************************************************************!

   logical function integer_tuple_less(left, right)
      implicit none
      integer, dimension(:), intent(in) :: left, right
      integer :: d

      integer_tuple_less = .false.
      do d = 1, size(left)
         if (left(d) < right(d)) then
            integer_tuple_less = .true.
            return
         else if (left(d) > right(d)) then
            return
         end if
      end do
   end function integer_tuple_less

   !******************************************************************************!

   recursive subroutine sort_interaction_order(material, order, left, right)
      implicit none
      type(dlib_material), intent(in) :: material
      integer, dimension(:), intent(inout) :: order
      integer, intent(in) :: left, right
      integer :: i, j, temporary
      integer, dimension(6) :: pivot

      if (left >= right) return
      pivot = material%interaction_axis_index(:, order((left + right) / 2))
      i = left
      j = right
      do
         do while (i <= right)
            if (.not. integer_tuple_less(material%interaction_axis_index(:, order(i)), pivot)) exit
            i = i + 1
         end do
         do while (j >= left)
            if (.not. integer_tuple_less(pivot, material%interaction_axis_index(:, order(j)))) exit
            j = j - 1
         end do
         if (i > j) exit
         temporary = order(i)
         order(i) = order(j)
         order(j) = temporary
         i = i + 1
         j = j - 1
      end do
      if (left < j) call sort_interaction_order(material, order, left, j)
      if (i < right) call sort_interaction_order(material, order, i, right)
   end subroutine sort_interaction_order

   !******************************************************************************!

   subroutine sort_axis(axis, n)
      implicit none
      real(dp), dimension(:), intent(inout) :: axis
      integer, intent(in) :: n
      integer :: i, j
      real(dp) :: value

      do i = 2, n
         value = axis(i)
         j = i - 1
         do while (j >= 1)
            if (axis(j) <= value) exit
            axis(j + 1) = axis(j)
            j = j - 1
         end do
         axis(j + 1) = value
      end do
   end subroutine sort_axis

   !******************************************************************************!

   subroutine bracket_axis(axis, n, query, allow_extrapolation, lower, upper, fraction)
      implicit none
      real(dp), dimension(:), intent(in) :: axis
      integer, intent(in) :: n
      real(dp), intent(in) :: query
      logical, intent(in) :: allow_extrapolation
      real(dp), intent(out) :: lower, upper, fraction
      integer :: i

      if (n < 1) call dlib_error('empty REFLECTOR_INTERACTION axis')
      if (.not. ieee_finite(query)) call dlib_error('REFLECTOR_INTERACTION condition is not finite')
      if (n == 1) then
         if (abs(query - axis(1)) > key_tolerance .and. .not. allow_extrapolation) then
            call dlib_error('requested condition is outside a REFLECTOR_INTERACTION range')
         end if
         lower = axis(1)
         upper = axis(1)
         fraction = 0._dp
         return
      end if
      do i = 1, n
         if (abs(query - axis(i)) <= key_tolerance) then
            lower = axis(i)
            upper = axis(i)
            fraction = 0._dp
            return
         end if
      end do
      if ((query < axis(1) - key_tolerance .or. query > axis(n) + key_tolerance) .and. &
         .not. allow_extrapolation) then
         call dlib_error('requested condition is outside a REFLECTOR_INTERACTION range')
      end if
      if (query < axis(1)) then
         if (allow_extrapolation .and. query < axis(1) - key_tolerance) then
            lower = axis(1)
            upper = axis(2)
            fraction = (query - lower) / (upper - lower)
         else
            lower = axis(1)
            upper = axis(1)
            fraction = 0._dp
         end if
         return
      else if (query > axis(n)) then
         if (allow_extrapolation .and. query > axis(n) + key_tolerance) then
            lower = axis(n - 1)
            upper = axis(n)
            fraction = (query - lower) / (upper - lower)
         else
            lower = axis(n)
            upper = axis(n)
            fraction = 0._dp
         end if
         return
      else
         do i = 2, n
            if (query < axis(i)) then
               lower = axis(i - 1)
               upper = axis(i)
               exit
            end if
         end do
      end if
      fraction = (query - lower) / (upper - lower)
   end subroutine bracket_axis

   !******************************************************************************!

   subroutine apply_dlib_rods(bpos)
      use sdata, only: nxx, nyy, nzz, xyz, zdel, fbmap, coreh, pos0, ssize, mat

      implicit none

      real(dp), dimension(:), intent(in) :: bpos
      integer :: i, j, k, bank, node, material_index, rod_curve
      real(dp) :: rod_height, depth, fraction
      character(len = 64) :: bank_label
      character(len = 16) :: bank_number

      do j = 1, nyy
         do i = 1, nxx
            bank = fbmap(i, j)
            if (bank <= 0) cycle
            if (bank > size(bpos)) call dlib_error('control rod map references an invalid bank')
            write(bank_number, '(I0)') bank
            bank_label = 'CR' // trim(bank_number)
            rod_height = coreh - pos0 - bpos(bank) * ssize
            depth = 0._dp
            do k = nzz, 1, -1
               node = xyz(i, j, k)
               if (rod_height >= depth .and. rod_height <= depth + zdel(k)) then
                  fraction = (rod_height - depth) / zdel(k)
                  if (fraction > 0._dp .and. dlib_is_fuel(mat(node))) then
                     material_index = mat_to_library(mat(node))
                     rod_curve = find_rod_curve(material_index, trim(bank_label))
                     if (rod_curve == 0) call dlib_error('rod curve ' // trim(bank_label) // &
                        ' is missing for a rodded material')
                     call add_rod_delta(node, material_index, rod_curve, fraction)
                  end if
                  exit
               else
                  if (dlib_is_fuel(mat(node))) then
                     material_index = mat_to_library(mat(node))
                     rod_curve = find_rod_curve(material_index, trim(bank_label))
                     if (rod_curve == 0) call dlib_error('rod curve ' // trim(bank_label) // &
                        ' is missing for a rodded material')
                     call add_rod_delta(node, material_index, rod_curve, 1._dp)
                  end if
               end if
               depth = depth + zdel(k)
            end do
         end do
      end do

   end subroutine apply_dlib_rods

   !******************************************************************************!

   integer function find_rod_curve(material_index, bank_label)
      implicit none
      integer, intent(in) :: material_index
      character(len = *), intent(in) :: bank_label
      integer :: i, matching_curves, only_curve

      find_rod_curve = 0
      matching_curves = 0
      only_curve = 0
      do i = 1, library_materials(material_index)%curve_count
         if (library_materials(material_index)%curves(i)%kind /= 'ROD') cycle
         if (trim(library_materials(material_index)%curves(i)%label) == trim(bank_label)) then
            find_rod_curve = i
            return
         end if
         matching_curves = matching_curves + 1
         only_curve = i
      end do
      ! CASMO can store one generic rod branch as CR1 while KOMODO has several core banks.
      if (matching_curves == 1) find_rod_curve = only_curve
   end function find_rod_curve

   !******************************************************************************!

   subroutine add_rod_delta(node, material_index, curve_index, fraction)
      use sdata, only: ng, sigtr, siga, nuf, sigf, sigs, dc
      implicit none

      integer, intent(in) :: node, material_index, curve_index
      real(dp), intent(in) :: fraction
      integer :: g, h
      real(dp) :: burnup

      burnup = node_burnup(node)
      do g = 1, ng
         sigtr(node, g) = sigtr(node, g) + fraction * curve_value(material_index, curve_index, index_transport(g), burnup)
         siga(node, g) = siga(node, g) + fraction * curve_value(material_index, curve_index, index_siga(g), burnup)
         nuf(node, g) = nuf(node, g) + fraction * curve_value(material_index, curve_index, index_nusigf(g), burnup)
         sigf(node, g) = sigf(node, g) + fraction * &
            curve_value(material_index, curve_index, index_kappa_fission(g), burnup)
         dc(node, g, :) = dc(node, g, :) + fraction * &
            curve_value(material_index, curve_index, index_adf(g), burnup)
         do h = 1, ng
            sigs(node, g, h) = sigs(node, g, h) + fraction * &
               curve_value(material_index, curve_index, index_scattering(g, h), burnup)
         end do
      end do
   end subroutine add_rod_delta

   !******************************************************************************!

   logical function dlib_is_fuel(material_id)
      integer, intent(in) :: material_id
      integer :: material_index

      dlib_is_fuel = .false.
      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (material_id < 1 .or. material_id > size(mat_to_library)) then
         call dlib_error('material ID is outside the DLIB mapping')
      end if
      material_index = mat_to_library(material_id)
      dlib_is_fuel = library_materials(material_index)%kind == 'FUEL'
   end function dlib_is_fuel

   !******************************************************************************!

   real(dp) function dlib_get_hm_density(material_id)
      integer, intent(in) :: material_id
      integer :: material_index

      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (material_id < 1 .or. material_id > size(mat_to_library)) then
         call dlib_error('material ID is outside the DLIB mapping')
      end if
      material_index = mat_to_library(material_id)
      dlib_get_hm_density = library_materials(material_index)%hm_smeared_density
   end function dlib_get_hm_density

   !******************************************************************************!

   subroutine dlib_get_node_burnup(burnup)
      implicit none
      real(dp), dimension(:), intent(out) :: burnup

      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (size(burnup) /= size(node_burnup)) call dlib_error('node burnup array has the wrong size')
      burnup = node_burnup
   end subroutine dlib_get_node_burnup

   !******************************************************************************!

   subroutine dlib_get_node_history(history_boron, history_moderator_temperature, history_fuel_temperature)
      implicit none
      real(dp), dimension(:), intent(out) :: history_boron
      real(dp), dimension(:), intent(out) :: history_moderator_temperature
      real(dp), dimension(:), intent(out) :: history_fuel_temperature

      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (size(history_boron) /= size(node_history_boron) .or. &
         size(history_moderator_temperature) /= size(node_history_moderator_temperature) .or. &
         size(history_fuel_temperature) /= size(node_history_fuel_temperature)) then
         call dlib_error('node history arrays have the wrong size')
      end if
      history_boron = node_history_boron
      history_moderator_temperature = node_history_moderator_temperature
      history_fuel_temperature = node_history_fuel_temperature
   end subroutine dlib_get_node_history

   !******************************************************************************!

   subroutine dlib_set_node_state(burnup, history_boron, history_moderator_temperature, &
      history_fuel_temperature)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      real(dp), dimension(:), intent(in) :: burnup, history_boron
      real(dp), dimension(:), intent(in) :: history_moderator_temperature, history_fuel_temperature

      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (size(burnup) /= size(node_burnup) .or. size(history_boron) /= size(node_burnup) .or. &
         size(history_moderator_temperature) /= size(node_burnup) .or. &
         size(history_fuel_temperature) /= size(node_burnup)) then
         call dlib_error('replacement node-state arrays have the wrong size')
      end if
      if (any(.not. ieee_is_finite(burnup)) .or. any(burnup < 0._dp) .or. &
         any(.not. ieee_is_finite(history_boron)) .or. &
         any(.not. ieee_is_finite(history_moderator_temperature)) .or. &
         any(.not. ieee_is_finite(history_fuel_temperature))) then
         call dlib_error('replacement node state contains invalid values')
      end if
      node_burnup = burnup
      node_history_boron = history_boron
      node_history_moderator_temperature = history_moderator_temperature
      node_history_fuel_temperature = history_fuel_temperature
   end subroutine dlib_set_node_state

   !******************************************************************************!

   subroutine dlib_get_material_name(material_id, name)
      implicit none
      integer, intent(in) :: material_id
      character(len = *), intent(out) :: name
      integer :: material_index

      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (material_id < 1 .or. material_id > size(mat_to_library)) then
         call dlib_error('material ID is outside the DLIB mapping')
      end if
      material_index = mat_to_library(material_id)
      name = library_materials(material_index)%name
   end subroutine dlib_get_material_name

   !******************************************************************************!

   subroutine dlib_update_node_history(burnup, current_boron, current_ftem, current_mtem)
      implicit none
      real(dp), dimension(:), intent(in) :: burnup
      real(dp), intent(in) :: current_boron
      real(dp), dimension(:), intent(in) :: current_ftem, current_mtem
      integer :: n
      real(dp) :: burnup_increment, previous_burnup

      if (.not. library_loaded) call dlib_error('library data are not initialized')
      if (size(burnup) /= size(node_burnup)) call dlib_error('node burnup array has the wrong size')
      if (size(current_ftem) /= size(node_burnup) .or. size(current_mtem) /= size(node_burnup)) then
         call dlib_error('history temperature arrays do not match the node count')
      end if
      if (.not. ieee_finite(current_boron)) call dlib_error('current boron concentration is not finite')
      do n = 1, size(burnup)
         if (.not. ieee_finite(burnup(n)) .or. burnup(n) < 0._dp .or. &
            burnup(n) < node_burnup(n)) then
            call dlib_error('node burnup must be finite, nonnegative, and nondecreasing')
         end if
         if (.not. ieee_finite(current_ftem(n)) .or. .not. ieee_finite(current_mtem(n))) then
            call dlib_error('current node temperature is not finite')
         end if

         previous_burnup = node_burnup(n)
         burnup_increment = burnup(n) - previous_burnup
         if (burnup_increment > 0._dp) then
            ! The initial burnup contributes its reference condition to this average.
            node_history_boron(n) = (previous_burnup * node_history_boron(n) + &
               burnup_increment * current_boron) / burnup(n)
            node_history_moderator_temperature(n) = &
               (previous_burnup * node_history_moderator_temperature(n) + &
               burnup_increment * current_mtem(n)) / burnup(n)
            node_history_fuel_temperature(n) = &
               (previous_burnup * node_history_fuel_temperature(n) + &
               burnup_increment * current_ftem(n)) / burnup(n)
            if (.not. ieee_finite(node_history_boron(n)) .or. &
               .not. ieee_finite(node_history_moderator_temperature(n)) .or. &
               .not. ieee_finite(node_history_fuel_temperature(n))) then
               call dlib_error('updated history condition is not finite')
            end if
         end if
         node_burnup(n) = burnup(n)
      end do
   end subroutine dlib_update_node_history

   !******************************************************************************!

   logical function ieee_finite(value)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      real(dp), intent(in) :: value
      ieee_finite = ieee_is_finite(value)
   end function ieee_finite

   !******************************************************************************!

   subroutine dlib_error(error_text)
      implicit none
      character(len = *), intent(in) :: error_text

      write(ounit, *) '   ERROR: %DLIB ', trim(error_text)
      write(*, *) '   ERROR: %DLIB ', trim(error_text)
      stop
   end subroutine dlib_error

end module dlib
