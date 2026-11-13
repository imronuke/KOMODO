module cycle

   use sdata, only: dp
   implicit none
   private

   public :: cycle_read_card, cycle_validate, cycle_write_state, cycle_prepare

   integer, parameter :: assembly_id_length = 32
   character(len = assembly_id_length), allocatable, save :: assembly_id(:, :)
   character(len = 256), save :: cycle_output_path = ''
   logical, save :: cycle_active = .false.

contains

   subroutine cycle_read_card(unit, nx, ny, output_unit, terminal_output)
      implicit none

      integer, intent(in) :: unit, nx, ny, output_unit
      logical, intent(in) :: terminal_output

      integer :: i, j, ios, line_number, pair_i, pair_j
      integer :: position, token_start, token_length, entry_count, line_length
      character(len = 2) :: indicator
      character(len = 200) :: line
      character(len = 200), dimension(nx) :: row_entry
      character(len = 200) :: token

      allocate(assembly_id(nx, ny))
      do j = ny, 1, -1
         read(unit, '(A2,I5,A200)', iostat = ios) indicator, line_number, line
         if (ios /= 0) call cycle_error(output_unit, terminal_output, &
            'could not read an assembly-map row in %CYCL')
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
            if (entry_count > nx) call cycle_error(output_unit, terminal_output, &
               'each %CYCL assembly-map row must contain exactly NX entries')
            if (token_length > assembly_id_length) call cycle_error(output_unit, terminal_output, &
               'assembly-map entry exceeds the 32-character ID limit')
            row_entry(entry_count) = line(token_start:position - 1)
         end do
         if (entry_count /= nx) call cycle_error(output_unit, terminal_output, &
            'each %CYCL assembly-map row must contain exactly NX entries')
         do i = 1, nx
            token = trim(adjustl(row_entry(i)))
            if (len_trim(token) == 0) call cycle_error(output_unit, terminal_output, &
               'assembly-map entries cannot be blank')
            assembly_id(i, j) = token
            if (same_id(token, '0')) assembly_id(i, j) = '0'
         end do
      end do

      do j = 1, ny
         do i = 1, nx
            if (same_id(assembly_id(i, j), '0')) cycle
            do pair_j = 1, j
               do pair_i = 1, nx
                  if (pair_j == j .and. pair_i >= i) exit
                  if (same_id(assembly_id(i, j), assembly_id(pair_i, pair_j))) then
                     call cycle_error(output_unit, terminal_output, &
                        'each nonzero %CYCL assembly-map ID must be unique')
                  end if
               end do
            end do
         end do
      end do

      read(unit, '(A2,I5,A200)', iostat = ios) indicator, line_number, line
      if (ios /= 0) call cycle_error(output_unit, terminal_output, &
         '%CYCL is missing the cycle-state output filename after the assembly map')
      cycle_output_path = trim(adjustl(line))
      if (len_trim(cycle_output_path) == 0) call cycle_error(output_unit, terminal_output, &
         '%CYCL cycle-state output filename cannot be blank')

      cycle_active = .true.
      write(output_unit, '(A)') ' %CYCL assembly map and cycle-state output file are read.'
      if (terminal_output) write(*, '(A)') ' %CYCL assembly map and cycle-state output file are read.'
   end subroutine cycle_read_card

   !******************************************************************************!

   subroutine cycle_validate(output_unit, terminal_output, fuel_material)
      use sdata, only: nx, ny, nxx, nyy, nzz, nnod, nmat, mat, xdiv, ydiv, xdel, ydel, ix, iy
      implicit none

      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      logical, dimension(:), intent(in) :: fuel_material

      integer :: i, j, node, ref_i, ref_j, global, idx
      integer, dimension(:), allocatable :: x_owner, y_owner
      logical, dimension(:, :), allocatable :: has_fuel
      logical :: have_reference

      if (.not. cycle_active) return
      if (size(assembly_id, 1) /= nx .or. size(assembly_id, 2) /= ny) then
         call cycle_error(output_unit, terminal_output, '%CYCL assembly-map dimensions do not match %GEOM')
      end if
      if (size(fuel_material) /= nmat) call cycle_error(output_unit, terminal_output, &
         '%CYCL fuel-material map does not match NMAT')

      if (sum(xdiv) /= nxx) call cycle_error(output_unit, terminal_output, &
         '%CYCL x node layout does not match the geometry')
      if (sum(ydiv) /= nyy .or. nzz < 1) call cycle_error(output_unit, terminal_output, &
         '%CYCL geometry node layout is inconsistent')

      allocate(x_owner(nxx), y_owner(nyy))
      call build_axis_map(nx, xdiv, x_owner)
      call build_axis_map(ny, ydiv, y_owner)

      allocate(has_fuel(nx, ny))
      has_fuel = .false.
      do node = 1, nnod
         if (.not. fuel_material(mat(node))) cycle
         i = x_owner(ix(node))
         j = y_owner(iy(node))
         has_fuel(i, j) = .true.
         if (same_id(assembly_id(i, j), '0')) then
            call cycle_error(output_unit, terminal_output, &
               'a fuel-containing assembly position cannot be 0 in the %CYCL map')
         end if
      end do

      do j = 1, ny
         do i = 1, nx
            if (same_id(assembly_id(i, j), '0')) cycle
            if (.not. has_fuel(i, j)) then
               call cycle_error(output_unit, terminal_output, &
                  'a nonzero %CYCL assembly ID must identify a fuel-containing assembly')
            end if
         end do
      end do

      ! A cycle file stores each node in assembly-local coordinates. Require the
      ! same unrotated x/y subdivision and node widths for every fuel assembly.
      have_reference = .false.
      ref_i = 0
      ref_j = 0
      do j = 1, ny
         do i = 1, nx
            if (.not. has_fuel(i, j)) cycle
            if (.not. have_reference) then
               ref_i = i
               ref_j = j
               have_reference = .true.
               cycle
            end if
            if (xdiv(i) /= xdiv(ref_i) .or. ydiv(j) /= ydiv(ref_j)) then
               call cycle_error(output_unit, terminal_output, &
                  'all fuel assemblies in %CYCL must use the same node layout and orientation')
            end if
            idx = sum(xdiv(1:i - 1)) + 1
            global = sum(xdiv(1:ref_i - 1)) + 1
            if (.not. all(widths_match(xdel(idx:idx + xdiv(i) - 1), &
               xdel(global:global + xdiv(ref_i) - 1)))) then
               call cycle_error(output_unit, terminal_output, &
                  'all fuel assemblies in %CYCL must use the same node layout and orientation')
            end if
            idx = sum(ydiv(1:j - 1)) + 1
            global = sum(ydiv(1:ref_j - 1)) + 1
            if (.not. all(widths_match(ydel(idx:idx + ydiv(j) - 1), &
               ydel(global:global + ydiv(ref_j) - 1)))) then
               call cycle_error(output_unit, terminal_output, &
                  'all fuel assemblies in %CYCL must use the same node layout and orientation')
            end if
         end do
      end do

      deallocate(x_owner, y_owner, has_fuel)
   end subroutine cycle_validate

   !******************************************************************************!

   subroutine cycle_prepare(cycle_enabled, refuel_enabled, output_unit, terminal_output)
      use sdata, only: nmat, nnod
      use dlib, only: dlib_is_fuel, dlib_get_material_name, dlib_get_node_burnup, &
         dlib_get_node_history, dlib_set_node_state
      use refuel, only: refuel_apply
      implicit none

      logical, intent(in) :: cycle_enabled, refuel_enabled, terminal_output
      integer, intent(in) :: output_unit
      logical, dimension(:), allocatable :: fuel_material
      character(len = 64), dimension(:), allocatable :: material_name
      real(dp), dimension(:), allocatable :: burnup, history_boron
      real(dp), dimension(:), allocatable :: history_moderator_temperature, history_fuel_temperature
      integer :: material_id

      if (.not. cycle_enabled .and. .not. refuel_enabled) return

      allocate(fuel_material(nmat), material_name(nmat))
      do material_id = 1, nmat
         fuel_material(material_id) = dlib_is_fuel(material_id)
         call dlib_get_material_name(material_id, material_name(material_id))
      end do
      if (cycle_enabled) call cycle_validate(output_unit, terminal_output, fuel_material)

      if (refuel_enabled) then
         allocate(burnup(nnod), history_boron(nnod), &
            history_moderator_temperature(nnod), history_fuel_temperature(nnod))
         call dlib_get_node_burnup(burnup)
         call dlib_get_node_history(history_boron, history_moderator_temperature, history_fuel_temperature)
         call refuel_apply(output_unit, terminal_output, fuel_material, material_name, burnup, history_boron, &
            history_moderator_temperature, history_fuel_temperature)
         call dlib_set_node_state(burnup, history_boron, &
            history_moderator_temperature, history_fuel_temperature)
      end if
   end subroutine cycle_prepare

   !******************************************************************************!

   subroutine cycle_write_state(output_unit, terminal_output, fuel_material, hm_density_by_material, &
      material_name_by_material, burnup, hboron, htmod, htfuel)
      use sdata, only: nx, ny, nxx, nyy, nzz, nnod, nmat, ix, iy, iz, mat, &
         xdiv, ydiv, xdel, ydel, zdel, vdel
      implicit none

      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      logical, dimension(:), intent(in) :: fuel_material
      real(dp), dimension(:), intent(in) :: hm_density_by_material, burnup, hboron, htmod, htfuel
      character(len = *), dimension(:), intent(in) :: material_name_by_material

      integer :: i, j, node, unit, ios, local_x, local_y, asm_x, asm_y, ref_i, ref_j, global
      integer, dimension(:), allocatable :: x_owner, x_local, y_owner, y_local
      real(dp), dimension(:, :), allocatable :: hm_sum, bu_sum
      real(dp) :: hm_mass, node_burnup

      if (.not. cycle_active) return
      if (size(fuel_material) /= nmat .or. size(hm_density_by_material) /= nmat .or. &
         size(material_name_by_material) /= nmat .or. size(burnup) /= nnod .or. &
         size(hboron) /= nnod .or. size(htmod) /= nnod .or. size(htfuel) /= nnod) then
         call cycle_error(output_unit, terminal_output, '%CYCL node-state arrays do not match the geometry')
      end if
      allocate(x_owner(nxx), x_local(nxx), y_owner(nyy), y_local(nyy))
      call build_axis_map(nx, xdiv, x_owner, x_local)
      call build_axis_map(ny, ydiv, y_owner, y_local)
      allocate(hm_sum(nx, ny), bu_sum(nx, ny))
      hm_sum = 0._dp
      bu_sum = 0._dp

      do node = 1, nnod
         if (.not. fuel_material(mat(node))) cycle
         asm_x = x_owner(ix(node))
         asm_y = y_owner(iy(node))
         hm_mass = hm_density_by_material(mat(node)) * vdel(node) * 1.e-6_dp
         if (hm_mass <= 0._dp) call cycle_error(output_unit, terminal_output, &
            'fuel node has no positive HM mass while writing %CYCL state')
         hm_sum(asm_x, asm_y) = hm_sum(asm_x, asm_y) + hm_mass
         bu_sum(asm_x, asm_y) = bu_sum(asm_x, asm_y) + hm_mass * burnup(node)
      end do
      do j = 1, ny
         do i = 1, nx
            if (hm_sum(i, j) > 0._dp) then
               bu_sum(i, j) = bu_sum(i, j) / hm_sum(i, j)
            end if
         end do
      end do

      call print_assembly_map(output_unit, nx, ny, 'Assembly IDs')
      call print_burnup_map(output_unit, nx, ny, bu_sum)
      if (terminal_output) then
         call print_assembly_map(6, nx, ny, 'Assembly IDs')
         call print_burnup_map(6, nx, ny, bu_sum)
      end if

      open(newunit = unit, file = trim(cycle_output_path), status = 'replace', &
         action = 'write', iostat = ios)
      if (ios /= 0) call cycle_error(output_unit, terminal_output, &
         'could not open %CYCL output file: ' // trim(cycle_output_path))
      write(unit, '(A)') 'KOMODO_CYCLE 1'
      write(unit, '(A,4(1X,I0))') 'DIMENSIONS', nx, ny, nzz, nnod
      write(unit, '(A)') 'BURNUP_UNIT MWd/kgHM'
      write(unit, '(A)') 'HISTORY_UNITS boron_ppm moderator_temperature_K fuel_temperature_K'
      ref_i = 0
      ref_j = 0
      do j = 1, ny
         do i = 1, nx
            if (same_id(assembly_id(i, j), '0')) cycle
            ref_i = i
            ref_j = j
            exit
         end do
         if (ref_i /= 0) exit
      end do
      if (ref_i == 0) call cycle_error(output_unit, terminal_output, &
         'no fuel assembly is available for cycle node-layout metadata')
      write(unit, '(A,2(1X,I0))') 'ASSEMBLY_NODE_LAYOUT', xdiv(ref_i), ydiv(ref_j)
      write(unit, '(A)', advance = 'no') 'X_NODE_WIDTHS'
      do local_x = 1, xdiv(ref_i)
         global = sum(xdiv(1:ref_i - 1)) + local_x
         write(unit, '(1X,ES24.16E3)', advance = 'no') xdel(global)
      end do
      write(unit, '(A)') ''
      write(unit, '(A)', advance = 'no') 'Y_NODE_WIDTHS'
      do local_y = 1, ydiv(ref_j)
         global = sum(ydiv(1:ref_j - 1)) + local_y
         write(unit, '(1X,ES24.16E3)', advance = 'no') ydel(global)
      end do
      write(unit, '(A)') ''
      write(unit, '(A)', advance = 'no') 'Z_NODE_WIDTHS'
      do global = 1, nzz
         write(unit, '(1X,ES24.16E3)', advance = 'no') zdel(global)
      end do
      write(unit, '(A)') ''
      write(unit, '(A)') 'ASSEMBLIES i j id hm_weighted_axial_average_burnup'
      do j = ny, 1, -1
         do i = 1, nx
            write(unit, '(A,2(1X,I0),1X,A,1X,ES24.16E3)') &
               'ASSEMBLY', i, j, trim(assembly_id(i, j)), bu_sum(i, j)
         end do
      end do
      write(unit, '(A)') &
         'NODE_STATE node assembly_i assembly_j local_i local_j axial_k assembly_id material_name ' // &
         'burnup history_boron history_moderator_temperature history_fuel_temperature'
      do node = 1, nnod
         asm_x = x_owner(ix(node))
         asm_y = y_owner(iy(node))
         local_x = x_local(ix(node))
         local_y = y_local(iy(node))
         if (len_trim(material_name_by_material(mat(node))) == 0) then
            call cycle_error(output_unit, terminal_output, 'DLIB material name is blank for a node')
         end if
         node_burnup = burnup(node)
         if (.not. fuel_material(mat(node))) node_burnup = 0._dp
         write(unit, '(A,6(1X,I0),1X,A,1X,A,4(1X,ES24.16E3))') &
            'NODE', node, asm_x, asm_y, local_x, local_y, iz(node), &
            trim(assembly_id(asm_x, asm_y)), trim(material_name_by_material(mat(node))), node_burnup, &
            hboron(node), htmod(node), htfuel(node)
      end do
      write(unit, '(A)') 'END_CYCLE'
      close(unit)

      write(output_unit, '(A,A)') ' %CYCL state written to ', trim(cycle_output_path)
      if (terminal_output) write(*, '(A,A)') ' %CYCL state written to ', trim(cycle_output_path)
      deallocate(x_owner, x_local, y_owner, y_local, hm_sum, bu_sum)
   end subroutine cycle_write_state

   !******************************************************************************!

   subroutine build_axis_map(nassembly, divisions, owner, local_index)
      implicit none
      integer, intent(in) :: nassembly
      integer, dimension(:), intent(in) :: divisions
      integer, dimension(:), intent(out) :: owner
      integer, dimension(:), optional, intent(out) :: local_index
      integer :: a, local, global

      global = 0
      do a = 1, nassembly
         do local = 1, divisions(a)
            global = global + 1
            owner(global) = a
            if (present(local_index)) local_index(global) = local
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

   subroutine print_assembly_map(unit, nx, ny, title)
      implicit none
      integer, intent(in) :: unit, nx, ny
      character(len = *), intent(in) :: title
      integer :: i, j, field_width, max_id_length
      character(len = 32) :: number_text
      character(len = :), dimension(:, :), allocatable :: entry

      max_id_length = 1
      do j = 1, ny
         do i = 1, nx
            if (same_id(assembly_id(i, j), '0')) cycle
            max_id_length = max(max_id_length, len_trim(assembly_id(i, j)))
         end do
      end do
      write(number_text, '(I0)') nx
      field_width = max(max_id_length, len_trim(number_text)) + 1
      allocate(character(len = field_width) :: entry(nx, ny))
      entry = repeat(' ', field_width)
      do j = 1, ny
         do i = 1, nx
            if (.not. same_id(assembly_id(i, j), '0')) entry(i, j) = trim(assembly_id(i, j))
         end do
      end do
      write(unit, '(A)') ''
      write(unit, '(A)') ' ' // trim(title) // ' (fuel assemblies only):'
      call print_character_grid(unit, nx, ny, field_width, 5, 1, entry)
      deallocate(entry)
   end subroutine print_assembly_map

   !******************************************************************************!

   subroutine print_burnup_map(unit, nx, ny, burnup)
      implicit none
      integer, intent(in) :: unit, nx, ny
      real(dp), dimension(:, :), intent(in) :: burnup
      integer :: i, j, field_width, max_burnup_length, value_start
      character(len = 32) :: number_text, burnup_text
      character(len = :), dimension(:, :), allocatable :: entry

      max_burnup_length = 1
      do j = 1, ny
         do i = 1, nx
            if (same_id(assembly_id(i, j), '0')) cycle
            write(burnup_text, '(F24.3)') burnup(i, j)
            burnup_text = trim(adjustl(burnup_text))
            max_burnup_length = max(max_burnup_length, len_trim(burnup_text))
         end do
      end do
      write(number_text, '(I0)') nx
      field_width = max(max_burnup_length + 2, len_trim(number_text) + 2, 8)
      allocate(character(len = field_width) :: entry(nx, ny))
      entry = repeat(' ', field_width)
      do j = 1, ny
         do i = 1, nx
            if (same_id(assembly_id(i, j), '0')) cycle
            write(burnup_text, '(F24.3)') burnup(i, j)
            burnup_text = trim(adjustl(burnup_text))
            value_start = field_width - 1 - len_trim(burnup_text)
            entry(i, j)(value_start:field_width - 2) = trim(burnup_text)
         end do
      end do
      write(unit, '(A)') ' Axially averaged assembly burnup (MWd/kgHM; HM-mass weighted):'
      call print_character_grid(unit, nx, ny, field_width, 6, 3, entry)
      deallocate(entry)
   end subroutine print_burnup_map

   !******************************************************************************!

   subroutine print_character_grid(unit, nx, ny, field_width, prefix_width, row_label_start, entries)
      implicit none
      integer, intent(in) :: unit, nx, ny, field_width, prefix_width, row_label_start
      character(len = *), dimension(:, :), intent(in) :: entries

      integer :: i, j, start, line_width
      character(len = 32) :: number_text
      character(len = :), allocatable :: output_line

      line_width = prefix_width + nx * field_width
      allocate(character(len = line_width) :: output_line)
      output_line = repeat(' ', line_width)
      do i = 1, nx
         write(number_text, '(I0)') i
         start = prefix_width + 1 + (i - 1) * field_width
         output_line(start:start + len_trim(number_text) - 1) = trim(number_text)
      end do
      write(unit, '(A)') trim(output_line)

      do j = ny, 1, -1
         if (.not. any(entries(:, j) /= ' ')) cycle
         output_line = repeat(' ', line_width)
         write(number_text, '(I4)') j
         output_line(row_label_start:row_label_start + 3) = number_text(1:4)
         do i = 1, nx
            start = prefix_width + 1 + (i - 1) * field_width
            output_line(start:start + field_width - 1) = entries(i, j)
         end do
         write(unit, '(A)') trim(output_line)
      end do
      deallocate(output_line)
   end subroutine print_character_grid

   !******************************************************************************!

   logical function same_id(left, right)
      implicit none
      character(len = *), intent(in) :: left, right
      character(len = 32) :: lhs, rhs
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
      same_id = trim(lhs) == trim(rhs)
   end function same_id

   !******************************************************************************!

   subroutine cycle_error(output_unit, terminal_output, message)
      implicit none
      integer, intent(in) :: output_unit
      logical, intent(in) :: terminal_output
      character(len = *), intent(in) :: message

      write(output_unit, '(A)') '   ERROR: %CYCL ' // trim(message)
      if (terminal_output) write(*, '(A)') '   ERROR: %CYCL ' // trim(message)
      stop
   end subroutine cycle_error

end module cycle
