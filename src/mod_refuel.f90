module refuel

   use sdata, only: dp
   implicit none
   private

   public :: refuel_read_card, refuel_apply

   integer, parameter :: id_length = 32
   integer, parameter :: cycle_unit = 997
   character(len = id_length), dimension(:, :), allocatable, save :: refuel_assembly
   character(len = 256), save :: refuel_cycle_path = ''
   logical, save :: refuel_active = .false.

   type :: cycle_state
      integer :: nx = 0, ny = 0, nzz = 0, nnod = 0
      integer :: nlocalx = 0, nlocaly = 0
      character(len = id_length), dimension(:, :), allocatable :: assembly_id
      real(dp), dimension(:), allocatable :: x_width, y_width, z_width
   end type cycle_state

contains

   subroutine refuel_read_card(unit, nx, ny, output_unit, terminal_output)
      implicit none

      integer, intent(in) :: unit, nx, ny, output_unit
      logical, intent(in) :: terminal_output

      integer :: i, j, ii, jj, ios, line_number
      character(len = 2) :: indicator
      character(len = 200) :: line
      character(len = 200), dimension(nx) :: row_entry
      character(len = 200) :: token
      integer :: position, token_start, token_length, entry_count, line_length

      if (allocated(refuel_assembly)) deallocate(refuel_assembly)
      allocate(refuel_assembly(nx, ny))

      read(unit, '(A2,I5,A200)', iostat = ios) indicator, line_number, line
      if (ios /= 0) call refuel_error(output_unit, terminal_output, &
         '%REFU is missing the prior-cycle filename')
      refuel_cycle_path = trim(adjustl(line))
      if (len_trim(refuel_cycle_path) == 0) call refuel_error(output_unit, terminal_output, &
         '%REFU prior-cycle filename cannot be blank')

      do j = ny, 1, -1
         read(unit, '(A2,I5,A200)', iostat = ios) indicator, line_number, line
         if (ios /= 0) call refuel_error(output_unit, terminal_output, &
            'could not read an assembly-map row in %REFU')
         row_entry = ''
         position = 1
         entry_count = 0
         line_length = len_trim(line)
         do while (position <= line_length)
            if (line(position:position) == ' ' .or. line(position:position) == achar(9)) then
               position = position + 1
               cycle
            end if
            token_start = position
            do while (position <= line_length)
               if (line(position:position) == ' ' .or. line(position:position) == achar(9)) exit
               position = position + 1
            end do
            token_length = position - token_start
            entry_count = entry_count + 1
            if (entry_count > nx) call refuel_error(output_unit, terminal_output, &
               'each %REFU assembly-map row must contain exactly NX entries')
            if (token_length > id_length) call refuel_error(output_unit, terminal_output, &
               'assembly-map entry exceeds the 32-character ID limit')
            row_entry(entry_count) = line(token_start:position - 1)
         end do
         if (entry_count /= nx) call refuel_error(output_unit, terminal_output, &
            'each %REFU assembly-map row must contain exactly NX entries')
         do i = 1, nx
            token = trim(adjustl(row_entry(i)))
            if (len_trim(token) == 0) call refuel_error(output_unit, terminal_output, &
               'assembly-map entries cannot be blank')
            refuel_assembly(i, j) = token
            if (same_word(token, '0')) refuel_assembly(i, j) = '0'
            if (same_word(token, 'XXX')) refuel_assembly(i, j) = 'XXX'
         end do
      end do

      do j = 1, ny
         do i = 1, nx
            if (same_word(refuel_assembly(i, j), '0') .or. &
               same_word(refuel_assembly(i, j), 'XXX')) cycle
            do jj = 1, j
               do ii = 1, nx
                  if (jj == j .and. ii >= i) exit
                  if (same_word(refuel_assembly(ii, jj), refuel_assembly(i, j))) then
                     call refuel_error(output_unit, terminal_output, &
                        'each prior-cycle assembly ID in the %REFU map must be unique')
                  end if
               end do
            end do
         end do
      end do

      refuel_active = .true.
      write(output_unit, '(A)') ' %REFU prior-cycle state and assembly map are read.'
      if (terminal_output) write(*, '(A)') ' %REFU prior-cycle state and assembly map are read.'
   end subroutine refuel_read_card

   !******************************************************************************!

   subroutine refuel_apply(output_unit, terminal_output, fuel_material, material_name_by_material, &
      burnup, history_boron, history_moderator_temperature, history_fuel_temperature)
      use sdata, only: nx, ny, nxx, nyy, nzz, nnod, nmat, ix, iy, mat, xyz, &
         xdiv, ydiv, xdel, ydel, zdel
      implicit none

      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      logical, dimension(:), intent(in) :: fuel_material
      character(len = *), dimension(:), intent(in) :: material_name_by_material
      real(dp), dimension(:), intent(inout) :: burnup, history_boron
      real(dp), dimension(:), intent(inout) :: history_moderator_temperature
      real(dp), dimension(:), intent(inout) :: history_fuel_temperature

      type(cycle_state) :: old_cycle
      integer, dimension(:), allocatable :: x_owner, y_owner
      integer, dimension(:), allocatable :: x_offset, y_offset
      integer, dimension(:, :), allocatable :: source_to_dest_x, source_to_dest_y
      integer, dimension(:, :), allocatable :: transferred_count
      logical, dimension(:, :), allocatable :: has_fuel
      logical, dimension(:), allocatable :: seen_node, seen_destination
      integer :: node, i, j, axial, global_x, global_y
      integer :: source_x, source_y, destination_x, destination_y, new_node
      integer :: asm_i, asm_j, record_local_x, record_local_y, record_axial
      integer :: reused_count, fresh_count
      integer :: ios
      character(len = 32768) :: line
      character(len = 64) :: keyword, material
      character(len = id_length) :: record_assembly_id
      real(dp) :: node_burnup, node_history_boron
      real(dp) :: node_history_moderator_temperature, node_history_fuel_temperature

      if (.not. refuel_active) return
      if (size(fuel_material) /= nmat .or. size(material_name_by_material) /= nmat .or. &
         size(burnup) /= nnod .or. size(history_boron) /= nnod .or. &
         size(history_moderator_temperature) /= nnod .or. &
         size(history_fuel_temperature) /= nnod) then
         call refuel_error(output_unit, terminal_output, 'node-state arrays do not match the geometry')
      end if
      if (size(refuel_assembly, 1) /= nx .or. size(refuel_assembly, 2) /= ny) then
         call refuel_error(output_unit, terminal_output, 'assembly-map dimensions do not match %GEOM')
      end if

      call read_cycle_state(trim(refuel_cycle_path), output_unit, terminal_output, old_cycle)
      if (old_cycle%nx /= nx .or. old_cycle%ny /= ny .or. old_cycle%nzz /= nzz .or. &
         old_cycle%nnod /= nnod) then
         call refuel_error(output_unit, terminal_output, &
            'prior-cycle dimensions do not match the current %GEOM node layout')
      end if

      allocate(x_owner(nxx), y_owner(nyy))
      if (sum(xdiv) /= nxx .or. sum(ydiv) /= nyy) then
         call refuel_error(output_unit, terminal_output, 'geometry node subdivisions do not match %GEOM')
      end if
      call build_axis_map(nx, xdiv, x_owner)
      call build_axis_map(ny, ydiv, y_owner)
      allocate(x_offset(nx), y_offset(ny))
      x_offset(1) = 0
      do i = 2, nx
         x_offset(i) = x_offset(i - 1) + xdiv(i - 1)
      end do
      y_offset(1) = 0
      do j = 2, ny
         y_offset(j) = y_offset(j - 1) + ydiv(j - 1)
      end do
      allocate(has_fuel(nx, ny))
      has_fuel = .false.

      do node = 1, nnod
         if (mat(node) < 1 .or. mat(node) > nmat) then
            call refuel_error(output_unit, terminal_output, 'node material ID is outside NMAT')
         end if
         i = x_owner(ix(node))
         j = y_owner(iy(node))
         has_fuel(i, j) = has_fuel(i, j) .or. fuel_material(mat(node))
      end do

      ! Validate fuel/non-fuel positions and require the source layout and widths to match.
      do j = 1, ny
         do i = 1, nx
            if (same_word(refuel_assembly(i, j), '0')) then
               if (has_fuel(i, j)) call refuel_error(output_unit, terminal_output, &
                  'a fuel-containing position cannot be 0 in the %REFU map')
            else
               if (.not. has_fuel(i, j)) call refuel_error(output_unit, terminal_output, &
                  'XXX or a prior-cycle ID must identify a fuel-containing position')
            end if
            if (.not. has_fuel(i, j)) cycle
            if (xdiv(i) /= old_cycle%nlocalx .or. ydiv(j) /= old_cycle%nlocaly) then
               call refuel_error(output_unit, terminal_output, &
                  'assembly node subdivisions do not match the prior-cycle file')
            end if
            global_x = x_offset(i) + 1
            if (.not. all(widths_match(xdel(global_x:global_x + old_cycle%nlocalx - 1), &
               old_cycle%x_width))) then
               call refuel_error(output_unit, terminal_output, &
                  'assembly x-node widths do not match the prior-cycle file')
            end if
            global_y = y_offset(j) + 1
            if (.not. all(widths_match(ydel(global_y:global_y + old_cycle%nlocaly - 1), &
               old_cycle%y_width))) then
               call refuel_error(output_unit, terminal_output, &
                  'assembly y-node widths do not match the prior-cycle file')
            end if
         end do
      end do
      if (.not. all(widths_match(zdel, old_cycle%z_width))) then
         call refuel_error(output_unit, terminal_output, &
            'axial node widths do not match the prior-cycle file')
      end if

      reused_count = 0
      fresh_count = 0
      allocate(source_to_dest_x(old_cycle%nx, old_cycle%ny), &
         source_to_dest_y(old_cycle%nx, old_cycle%ny), transferred_count(old_cycle%nx, old_cycle%ny))
      source_to_dest_x = 0
      source_to_dest_y = 0
      transferred_count = 0
      do j = 1, ny
         do i = 1, nx
            if (same_word(refuel_assembly(i, j), '0')) cycle
            if (same_word(refuel_assembly(i, j), 'XXX')) then
               fresh_count = fresh_count + 1
               cycle
            end if
            call find_source_assembly(old_cycle, refuel_assembly(i, j), source_x, source_y, &
               output_unit, terminal_output)
            reused_count = reused_count + 1
            source_to_dest_x(source_x, source_y) = i
            source_to_dest_y(source_x, source_y) = j
         end do
      end do

      allocate(seen_node(nnod), seen_destination(nnod))
      seen_node = .false.
      seen_destination = .false.
      do node = 1, old_cycle%nnod
         call next_cycle_line(cycle_unit, line, ios)
         if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle node data are incomplete')
         keyword = ''
         read(line, *, iostat = ios) keyword, i, asm_i, asm_j, record_local_x, record_local_y, &
            record_axial, record_assembly_id, material, node_burnup, node_history_boron, &
            node_history_moderator_temperature, node_history_fuel_temperature
         if (ios /= 0 .or. trim(keyword) /= 'NODE') then
            call refuel_error(output_unit, terminal_output, 'invalid prior-cycle node-state record')
         end if
         if (i < 1 .or. i > old_cycle%nnod .or. asm_i < 1 .or. asm_i > old_cycle%nx .or. &
            asm_j < 1 .or. asm_j > old_cycle%ny .or. record_local_x < 1 .or. &
            record_local_y < 1 .or. record_axial < 1 .or. record_axial > old_cycle%nzz) then
            call refuel_error(output_unit, terminal_output, 'prior-cycle node coordinates are out of range')
         end if
         if (seen_node(i)) call refuel_error(output_unit, terminal_output, &
            'duplicate node ID in prior-cycle state')
         if (.not. same_word(record_assembly_id, old_cycle%assembly_id(asm_i, asm_j))) then
            call refuel_error(output_unit, terminal_output, 'prior-cycle node assembly ID does not match its map')
         end if
         if (len_trim(material) == 0) call refuel_error(output_unit, terminal_output, &
            'prior-cycle node material name is blank')
         if (.not. same_word(old_cycle%assembly_id(asm_i, asm_j), '0')) then
            if (record_local_x > old_cycle%nlocalx .or. record_local_y > old_cycle%nlocaly) then
               call refuel_error(output_unit, terminal_output, &
                  'prior-cycle fuel node coordinates exceed the stored assembly layout')
            end if
         end if
         if (.not. finite_value(node_burnup) .or. node_burnup < 0._dp .or. &
            .not. finite_value(node_history_boron) .or. &
            .not. finite_value(node_history_moderator_temperature) .or. &
            .not. finite_value(node_history_fuel_temperature)) then
            call refuel_error(output_unit, terminal_output, 'prior-cycle node state contains invalid values')
         end if
         seen_node(i) = .true.

         destination_x = source_to_dest_x(asm_i, asm_j)
         destination_y = source_to_dest_y(asm_i, asm_j)
         if (destination_x == 0 .or. destination_y == 0) cycle
         if (record_local_x > xdiv(destination_x) .or. record_local_y > ydiv(destination_y)) then
            call refuel_error(output_unit, terminal_output, &
               'prior-cycle node layout is incomplete for an assembly in %REFU')
         end if
         global_x = x_offset(destination_x) + record_local_x
         global_y = y_offset(destination_y) + record_local_y
         new_node = xyz(global_x, global_y, record_axial)
         if (new_node < 1 .or. new_node > nnod) call refuel_error(output_unit, terminal_output, &
            'current geometry has no node corresponding to a prior-cycle fuel node')
         if (seen_destination(new_node)) call refuel_error(output_unit, terminal_output, &
            'duplicate prior-cycle coordinates map to the same current node')
         seen_destination(new_node) = .true.
         transferred_count(asm_i, asm_j) = transferred_count(asm_i, asm_j) + 1
         if (mat(new_node) < 1 .or. mat(new_node) > nmat) then
            call refuel_error(output_unit, terminal_output, 'node material ID is outside NMAT')
         end if
         if (.not. same_word(material, material_name_by_material(mat(new_node)))) then
            call refuel_error(output_unit, terminal_output, &
               'material mismatch for prior-cycle assembly ' // trim(refuel_assembly(destination_x, destination_y)) // &
               ' at destination position; check %GEOM material assignments and assembly orientation')
         end if
         burnup(new_node) = node_burnup
         history_boron(new_node) = node_history_boron
         history_moderator_temperature(new_node) = node_history_moderator_temperature
         history_fuel_temperature(new_node) = node_history_fuel_temperature
      end do
      if (any(.not. seen_node)) call refuel_error(output_unit, terminal_output, &
         'prior-cycle state is missing node records')
      do source_y = 1, old_cycle%ny
         do source_x = 1, old_cycle%nx
            if (source_to_dest_x(source_x, source_y) == 0) cycle
            if (transferred_count(source_x, source_y) /= &
               old_cycle%nlocalx * old_cycle%nlocaly * old_cycle%nzz) then
               call refuel_error(output_unit, terminal_output, &
                  'prior-cycle node layout is incomplete for an assembly in %REFU')
            end if
         end do
      end do
      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle end marker is missing')
      keyword = ''
      read(line, *, iostat = ios) keyword
      if (ios /= 0 .or. trim(keyword) /= 'END_CYCLE') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle end marker')
      end if
      close(cycle_unit)

      write(output_unit, '(A,I0,A,I0,A)') ' %REFU state applied: ', reused_count, &
         ' reused assemblies, ', fresh_count, ' new assemblies.'
      if (terminal_output) write(*, '(A,I0,A,I0,A)') ' %REFU state applied: ', reused_count, &
         ' reused assemblies, ', fresh_count, ' new assemblies.'

      deallocate(x_owner, y_owner, x_offset, y_offset, has_fuel)
      deallocate(source_to_dest_x, source_to_dest_y, transferred_count, seen_node, seen_destination)
      call destroy_cycle_state(old_cycle)
   end subroutine refuel_apply

   !******************************************************************************!

   subroutine read_cycle_state(path, output_unit, terminal_output, state)
      implicit none

      character(len = *), intent(in) :: path
      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      type(cycle_state), intent(out) :: state

      integer :: ios, unit_version, i, j, ii, jj, asm_i, asm_j
      character(len = 32768) :: line
      character(len = 64) :: keyword
      character(len = 64) :: burnup_unit
      character(len = 64) :: unit_boron, unit_moderator_temperature, unit_fuel_temperature
      character(len = id_length) :: assembly_id
      real(dp) :: average_burnup

      open(unit = cycle_unit, file = path, status = 'old', action = 'read', iostat = ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, &
         'cannot open prior-cycle state file: ' // trim(path))

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle file is empty')
      keyword = ''
      read(line, *, iostat = ios) keyword, unit_version
      if (ios /= 0 .or. trim(keyword) /= 'KOMODO_CYCLE' .or. unit_version /= 1) then
         call refuel_error(output_unit, terminal_output, 'unsupported prior-cycle file format')
      end if

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle dimensions are missing')
      keyword = ''
      read(line, *, iostat = ios) keyword, state%nx, state%ny, state%nzz, state%nnod
      if (ios /= 0 .or. trim(keyword) /= 'DIMENSIONS' .or. min(state%nx, state%ny, &
         state%nzz, state%nnod) < 1) then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle dimensions')
      end if

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle burnup units are missing')
      keyword = ''
      burnup_unit = ''
      call read_text_token(line, 1, keyword)
      call read_text_token(line, 2, burnup_unit)
      if (trim(keyword) /= 'BURNUP_UNIT' .or. trim(burnup_unit) /= 'MWd/kgHM') then
         call refuel_error(output_unit, terminal_output, 'prior-cycle burnup units are not MWd/kgHM')
      end if
      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle history units are missing')
      keyword = ''
      unit_boron = ''
      unit_moderator_temperature = ''
      unit_fuel_temperature = ''
      read(line, *, iostat = ios) keyword, unit_boron, unit_moderator_temperature, unit_fuel_temperature
      if (ios /= 0 .or. trim(keyword) /= 'HISTORY_UNITS' .or. &
         trim(unit_boron) /= 'boron_ppm' .or. &
         trim(unit_moderator_temperature) /= 'moderator_temperature_K' .or. &
         trim(unit_fuel_temperature) /= 'fuel_temperature_K') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle history units')
      end if

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle assembly node layout is missing')
      read(line, *, iostat = ios) keyword, state%nlocalx, state%nlocaly
      if (ios /= 0 .or. trim(keyword) /= 'ASSEMBLY_NODE_LAYOUT' .or. &
         min(state%nlocalx, state%nlocaly) < 1) then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle assembly node layout')
      end if
      allocate(state%x_width(state%nlocalx), state%y_width(state%nlocaly), state%z_width(state%nzz))

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle x-node widths are missing')
      read(line, *, iostat = ios) keyword, state%x_width
      if (ios /= 0 .or. trim(keyword) /= 'X_NODE_WIDTHS') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle x-node widths')
      end if
      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle y-node widths are missing')
      read(line, *, iostat = ios) keyword, state%y_width
      if (ios /= 0 .or. trim(keyword) /= 'Y_NODE_WIDTHS') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle y-node widths')
      end if
      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle z-node widths are missing')
      read(line, *, iostat = ios) keyword, state%z_width
      if (ios /= 0 .or. trim(keyword) /= 'Z_NODE_WIDTHS') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle z-node widths')
      end if
      if (any(.not. finite_value(state%x_width)) .or. any(state%x_width <= 0._dp) .or. &
         any(.not. finite_value(state%y_width)) .or. any(state%y_width <= 0._dp) .or. &
         any(.not. finite_value(state%z_width)) .or. any(state%z_width <= 0._dp)) then
         call refuel_error(output_unit, terminal_output, 'prior-cycle node widths contain invalid values')
      end if

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle assembly map is missing')
      read(line, *, iostat = ios) keyword
      if (ios /= 0 .or. trim(keyword) /= 'ASSEMBLIES') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle assembly map header')
      end if
      allocate(state%assembly_id(state%nx, state%ny))
      state%assembly_id = ''
      do i = 1, state%nx * state%ny
         call next_cycle_line(cycle_unit, line, ios)
         if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle assembly map is incomplete')
         read(line, *, iostat = ios) keyword, asm_i, asm_j, assembly_id, average_burnup
         if (ios /= 0 .or. trim(keyword) /= 'ASSEMBLY') then
            call refuel_error(output_unit, terminal_output, 'invalid prior-cycle assembly map record')
         end if
         if (.not. finite_value(average_burnup) .or. average_burnup < 0._dp) then
            call refuel_error(output_unit, terminal_output, 'prior-cycle assembly burnup is invalid')
         end if
         if (asm_i < 1 .or. asm_i > state%nx .or. asm_j < 1 .or. asm_j > state%ny) then
            call refuel_error(output_unit, terminal_output, 'prior-cycle assembly coordinates are out of range')
         end if
         if (len_trim(state%assembly_id(asm_i, asm_j)) > 0) then
            call refuel_error(output_unit, terminal_output, 'duplicate position in prior-cycle assembly map')
         end if
         state%assembly_id(asm_i, asm_j) = assembly_id
      end do
      do j = 1, state%ny
         do i = 1, state%nx
            if (len_trim(state%assembly_id(i, j)) == 0) then
               call refuel_error(output_unit, terminal_output, 'prior-cycle assembly map has missing positions')
            end if
            if (same_word(state%assembly_id(i, j), '0')) cycle
            do jj = 1, j
               do ii = 1, state%nx
                  if (jj == j .and. ii >= i) exit
                  if (same_word(state%assembly_id(ii, jj), state%assembly_id(i, j))) then
                     call refuel_error(output_unit, terminal_output, 'prior-cycle assembly IDs are not unique')
                  end if
               end do
            end do
         end do
      end do

      call next_cycle_line(cycle_unit, line, ios)
      if (ios /= 0) call refuel_error(output_unit, terminal_output, 'prior-cycle node-state header is missing')
      read(line, *, iostat = ios) keyword
      if (ios /= 0 .or. trim(keyword) /= 'NODE_STATE') then
         call refuel_error(output_unit, terminal_output, 'invalid prior-cycle node-state header')
      end if
      ! Leave the file positioned at the first NODE record. The caller streams
      ! each record into the destination state instead of retaining a full copy.
   end subroutine read_cycle_state

   !******************************************************************************!

   subroutine find_source_assembly(state, wanted_id, source_x, source_y, output_unit, terminal_output)
      implicit none
      type(cycle_state), intent(in) :: state
      character(len = *), intent(in) :: wanted_id
      integer, intent(out) :: source_x, source_y
      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      integer :: i, j

      source_x = 0
      source_y = 0
      do j = 1, state%ny
         do i = 1, state%nx
            if (same_word(state%assembly_id(i, j), wanted_id)) then
               source_x = i
               source_y = j
               return
            end if
         end do
      end do
      call refuel_error(output_unit, terminal_output, &
         'assembly ID ' // trim(wanted_id) // ' was not found in the prior-cycle file')
   end subroutine find_source_assembly

   !******************************************************************************!

   subroutine next_cycle_line(unit, line, ios)
      implicit none
      integer, intent(in) :: unit
      character(len = *), intent(out) :: line
      integer, intent(out) :: ios

      do
         read(unit, '(A)', iostat = ios) line
         if (ios /= 0) return
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         return
      end do
   end subroutine next_cycle_line

   !******************************************************************************!

   subroutine read_text_token(line, token_number, token)
      implicit none
      character(len = *), intent(in) :: line
      integer, intent(in) :: token_number
      character(len = *), intent(out) :: token

      integer :: position, token_start, current_token

      token = ''
      position = 1
      current_token = 0
      do while (position <= len_trim(line))
         if (line(position:position) == ' ' .or. line(position:position) == achar(9)) then
            position = position + 1
            cycle
         end if
         token_start = position
         do while (position <= len_trim(line))
            if (line(position:position) == ' ' .or. line(position:position) == achar(9)) exit
            position = position + 1
         end do
         current_token = current_token + 1
         if (current_token == token_number) then
            if (position - token_start > len(token)) return
            token = line(token_start:position - 1)
            return
         end if
      end do
   end subroutine read_text_token

   !******************************************************************************!

   subroutine build_axis_map(nassembly, divisions, owner)
      implicit none
      integer, intent(in) :: nassembly
      integer, dimension(:), intent(in) :: divisions
      integer, dimension(:), intent(out) :: owner
      integer :: assembly, local, global

      global = 0
      do assembly = 1, nassembly
         do local = 1, divisions(assembly)
            global = global + 1
            owner(global) = assembly
         end do
      end do
   end subroutine build_axis_map

   !******************************************************************************!

   elemental logical function widths_match(left, right)
      implicit none
      real(dp), intent(in) :: left, right

      widths_match = abs(left - right) <= 1.e-10_dp * max(1._dp, abs(left), abs(right))
   end function widths_match

   !******************************************************************************!

   subroutine destroy_cycle_state(state)
      implicit none
      type(cycle_state), intent(inout) :: state
      if (allocated(state%assembly_id)) deallocate(state%assembly_id)
      if (allocated(state%x_width)) deallocate(state%x_width)
      if (allocated(state%y_width)) deallocate(state%y_width)
      if (allocated(state%z_width)) deallocate(state%z_width)
   end subroutine destroy_cycle_state

   !******************************************************************************!

   logical function same_word(left, right)
      implicit none
      character(len = *), intent(in) :: left, right
      character(len = max(len(left), len(right))) :: lhs, rhs
      integer :: i, c

      lhs = adjustl(left)
      rhs = adjustl(right)
      do i = 1, len_trim(lhs)
         c = iachar(lhs(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) lhs(i:i) = achar(c - 32)
      end do
      do i = 1, len_trim(rhs)
         c = iachar(rhs(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) rhs(i:i) = achar(c - 32)
      end do
      same_word = trim(lhs) == trim(rhs)
   end function same_word

   !******************************************************************************!

   elemental logical function finite_value(value)
      implicit none
      real(dp), intent(in) :: value

      finite_value = abs(value) <= huge(value)
   end function finite_value

   !******************************************************************************!

   subroutine refuel_error(output_unit, terminal_output, message)
      implicit none
      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      character(len = *), intent(in) :: message

      write(output_unit, '(A)') '   ERROR: %REFU ' // trim(message)
      if (terminal_output) write(*, '(A)') '   ERROR: %REFU ' // trim(message)
      stop
   end subroutine refuel_error

end module refuel
