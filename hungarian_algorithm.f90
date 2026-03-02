!> Hungarian Algorithm (Kuhn-Munkres) module for optimal assignment problems.
!>
!> Algorithm Variant:
!>   This implements the "reduction + minimum line cover + adjustment" variant
!>   of the Hungarian algorithm, NOT the classical "starred/primed zeros"
!>   formulation (Munkres 1957). The steps are:
!>     1. Row reduction: subtract row minima
!>     2. Column reduction: subtract column minima
!>     3. Find minimum line cover via Konig's theorem (max matching = min cover)
!>     4. If cover < n, adjust matrix: subtract min uncovered from uncovered,
!>        add it to doubly covered elements
!>     5. Repeat from step 3 until cover == n (guaranteed in <= n iterations
!>        for well-conditioned input)
!>     6. Extract optimal assignment from the final zero structure
!>
!> Invariants:
!>   - After each adjustment step, the number of independent zeros
!>     (maximum matching) strictly increases.
!>   - The matrix remains non-negative.
!>   - Termination is guaranteed within O(n) adjustment iterations for
!>     exact arithmetic; 
!>
!> Memory Ownership (C API):
!>   - The caller owns all buffers passed via pointers.
!>   - c_matrix_ptr must point to n*n contiguous doubles in COLUMN-MAJOR order.
!>   - c_assignments_ptr must point to at least n contiguous int32_t elements.
!>   - c_cost_ptr must point to a valid single double.
!>   - All buffers must remain valid for the duration of the call.
!>   - The function is NOT thread-safe for concurrent calls with shared buffers.
module hungarian_mod
   use iso_c_binding
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   ! Public API
   public :: hungarian_algorithm, f90_hungarian_algorithm

   ! Define precision constants
   integer, parameter :: f64 = C_DOUBLE

contains

   ! ------------------------------------------------------------------
   ! Subroutine: Maximum Bipartite Matching
   ! Purpose:    Finds the maximum number of matches in a bipartite graph.
   ! Algorithm:  Iterative Depth-First Search (DFS) to find augmenting paths.
   ! ------------------------------------------------------------------
   subroutine maximum_bipartite_matching(adj, match_r_to_c)
      ! --- Inputs & Outputs using Assumed-Shape ---
      logical, intent(in) :: adj(:,:)
      integer, intent(out) :: match_r_to_c(:)

      ! --- Local Variables ---
      integer :: n
      ! Automatic arrays dimensioned via size()
      integer :: match_c_to_r(size(adj, 1))
      integer :: i, u, v, r, path_end_c, curr_c, prev_c_match, stack_top
      
      integer, allocatable :: stack(:), parent_row_of_col(:)
      logical, allocatable :: visited_cols(:)

      n = size(adj, 1)
      allocate(stack(n), visited_cols(n), parent_row_of_col(n))

      match_r_to_c = 0
      match_c_to_r = 0

      ! Try to find an augmenting path for every unmatched row 'i'
      do i = 1, n
         if (match_r_to_c(i) /= 0) cycle

         ! Initialize DFS
         stack_top = 1
         stack(stack_top) = i
         visited_cols = .false.
         parent_row_of_col = 0
         path_end_c = 0

         ! Iterative DFS Loop
         dfs_loop: do while (stack_top > 0)
            u = stack(stack_top)
            stack_top = stack_top - 1

            do v = 1, n
               ! If edge exists (zero in cost matrix) and column unvisited
               if (adj(u, v) .and. (.not. visited_cols(v))) then
                  visited_cols(v) = .true.
                  parent_row_of_col(v) = u

                  ! If column 'v' is unmatched, we found an augmenting path
                  if (match_c_to_r(v) == 0) then
                     path_end_c = v
                     exit dfs_loop
                  else
                     ! Stack depth is bounded mathematically by 'visited_cols'
                     stack_top = stack_top + 1
                     stack(stack_top) = match_c_to_r(v)
                  end if
               end if
            end do
         end do dfs_loop

         ! If path found, backtrack to flip edges and increase matching
         if (path_end_c /= 0) then
            curr_c = path_end_c
            do
               r = parent_row_of_col(curr_c)
               prev_c_match = match_r_to_c(r)
               
               match_c_to_r(curr_c) = r
               match_r_to_c(r) = curr_c
               
               if (r == i) exit
               curr_c = prev_c_match
            end do
         end if
      end do

      deallocate(stack, visited_cols, parent_row_of_col)
   end subroutine maximum_bipartite_matching

   ! ------------------------------------------------------------------
   ! Subroutine: Konig's Theorem (Minimum Line Cover)
   ! Purpose:    Finds the minimum number of lines needed to cover all zeros.
   ! Logic:      Min Vertex Cover in Bipartite Graph = Max Matching.
   ! ------------------------------------------------------------------
   subroutine minimum_line_cover(zeros_mat, lines_rows, lines_cols)
      ! --- Inputs & Outputs using Assumed-Shape ---
      logical, intent(in) :: zeros_mat(:,:)
      logical, intent(out) :: lines_rows(:), lines_cols(:)

      ! --- Local Variables ---
      integer :: n
      ! Automatic arrays dimensioned via size()
      integer :: match_r_to_c(size(zeros_mat, 1)), match_c_to_r(size(zeros_mat, 1))
      logical :: marked_rows(size(zeros_mat, 1)), marked_cols(size(zeros_mat, 1))
      
      integer, allocatable :: queue(:)
      integer :: q_head, q_tail, r, c, row, col, matched_row

      n = size(zeros_mat, 1)

      ! 1. Find Max Matching
      call maximum_bipartite_matching(zeros_mat, match_r_to_c)

      ! Build inverse lookup
      match_c_to_r = 0
      do r = 1, n
         c = match_r_to_c(r)
         if (c /= 0) match_c_to_r(c) = r
      end do

      ! 2. BFS for Konig's Construction
      allocate(queue(n))
      marked_rows = .false.
      marked_cols = .false.
      q_head = 1
      q_tail = 0

      ! Start BFS from all UNMATCHED rows
      do r = 1, n
         if (match_r_to_c(r) == 0) then
            marked_rows(r) = .true.
            q_tail = q_tail + 1
            queue(q_tail) = r
         end if
      end do

      do while (q_head <= q_tail)
         row = queue(q_head)
         q_head = q_head + 1

         do col = 1, n
            if (zeros_mat(row, col) .and. (.not. marked_cols(col))) then
               marked_cols(col) = .true.
               
               matched_row = match_c_to_r(col)
               if (matched_row /= 0 .and. (.not. marked_rows(matched_row))) then
                  marked_rows(matched_row) = .true.
                  q_tail = q_tail + 1
                  queue(q_tail) = matched_row
               end if
            end if
         end do
      end do

      ! 3. Construct Cover
      lines_rows = .not. marked_rows
      lines_cols = marked_cols

      deallocate(queue)
   end subroutine minimum_line_cover

   ! ==================================================================
   ! MAIN SOLVER: HUNGARIAN ALGORITHM
   ! ==================================================================
   subroutine hungarian_algorithm(cost_matrix, n, assignments, total_cost)
      real(f64), intent(in) :: cost_matrix(n, n)
      integer, intent(in) :: n
      integer, intent(out) :: assignments(n)
      real(f64), intent(out) :: total_cost

      real(f64), allocatable :: matrix(:,:)
      integer :: i, j
      real(f64) :: k, tol
      integer :: num_lines
      logical, allocatable :: zeros_bool(:,:)
      logical, allocatable :: lines_rows(:), lines_cols(:)
      integer, allocatable :: line_viz(:,:)
      logical, allocatable :: uncovered_mask(:,:), doubly_covered_mask(:,:)

      ! --- Input Validation ---
      if (n <= 0) then
         total_cost = ieee_value(1.0_f64, ieee_quiet_nan)
         return
      end if

      ! Check for NaN or Inf in the cost matrix
      if (any(.not. ieee_is_finite(cost_matrix))) then
         assignments = 0
         total_cost = ieee_value(1.0_f64, ieee_quiet_nan)
         return
      end if

      allocate(matrix(n, n))
      allocate(zeros_bool(n, n))
      allocate(lines_rows(n), lines_cols(n))
      allocate(line_viz(n, n))
      allocate(uncovered_mask(n, n), doubly_covered_mask(n, n))

      tol = epsilon(1.0_f64) * max(1.0_f64, maxval(abs(cost_matrix)))
      ! Create working copy
      matrix = cost_matrix

      ! --- Step 1: Subtract Row Minima ---
      do i = 1, n
         matrix(i, :) = matrix(i, :) - minval(matrix(i, :))
      end do

      ! --- Step 2: Subtract Column Minima ---
      do j = 1, n
         matrix(:, j) = matrix(:, j) - minval(matrix(:, j))
      end do

      ! --- Step 3: Iterative Reduction ---
      do
         zeros_bool = (abs(matrix) < tol)

         ! Cleaned up call: 'n' removed
         call minimum_line_cover(zeros_bool, lines_rows, lines_cols)
         num_lines = count(lines_rows) + count(lines_cols)

         if (num_lines == n) exit

         ! --- Step 4: Matrix Update ---
         line_viz = 0
         do i = 1, n
            if (lines_rows(i)) line_viz(i, :) = line_viz(i, :) + 1
            if (lines_cols(i)) line_viz(:, i) = line_viz(:, i) + 1
         end do

         uncovered_mask = (line_viz == 0)

         ! Find minimum value 'k' (Cache friendly loop order j, then i)
         k = huge(k)
         do j = 1, n
            do i = 1, n
               if (uncovered_mask(i, j)) then
                  if (matrix(i, j) < k) k = matrix(i, j)
               end if
            end do
         end do

         ! Subtract k from uncovered elements
         where (uncovered_mask)
            matrix = matrix - k
         end where

         ! Add k to doubly covered elements
         doubly_covered_mask = (line_viz == 2)
         where (doubly_covered_mask)
            matrix = matrix + k
         end where
      end do

      ! --- Step 5: Final Assignment ---
      zeros_bool = (abs(matrix) < tol)
      
      ! Cleaned up call: 'n' removed
      call maximum_bipartite_matching(zeros_bool, assignments)

      total_cost = 0.0_f64
      do i = 1, n
         if (assignments(i) /= 0) then
            total_cost = total_cost + cost_matrix(i, assignments(i))
         end if
      end do

      deallocate(matrix, zeros_bool, lines_rows, lines_cols)
      deallocate(line_viz, uncovered_mask, doubly_covered_mask)

   end subroutine hungarian_algorithm

   ! ==================================================================
   ! C WRAPPER
   ! ==================================================================
   subroutine f90_hungarian_algorithm(c_matrix_ptr, n, c_assignments_ptr, c_cost_ptr) bind(C, name="f90_hungarian_algorithm")
      type(c_ptr), value, intent(in) :: c_matrix_ptr
      integer(c_int32_t), value, intent(in) :: n
      type(c_ptr), value :: c_assignments_ptr
      type(c_ptr), value :: c_cost_ptr

      real(c_double), pointer :: ptr_matrix(:,:)
      integer(c_int32_t), pointer :: ptr_assign(:)
      real(c_double), pointer :: ptr_cost

      integer, allocatable :: assignments(:)
      real(f64) :: cost
      integer :: i, n_f

      n_f = int(n)

      ! Map Outputs First (so we can set them to NaN/-1 on error)
      call c_f_pointer(c_cost_ptr, ptr_cost)
      
      ! Handle edge case: N <= 0
      if (n_f <= 0) then
         ptr_cost = ieee_value(1.0_c_double, ieee_quiet_nan)
         return
      end if

      call c_f_pointer(c_assignments_ptr, ptr_assign, [n_f])

      ! 1. Zero-Copy Map
      call c_f_pointer(c_matrix_ptr, ptr_matrix, [n_f, n_f])

      allocate(assignments(n_f))

      ! 2. Solve 
      call hungarian_algorithm(ptr_matrix, n_f, assignments, cost)

      ! 3. Return Results
      ptr_cost = cost
      do i = 1, n_f
         if (assignments(i) /= 0 .and. ieee_is_finite(cost)) then
            ! Convert 1-based Fortran index to 0-based C index
            ptr_assign(i) = assignments(i) - 1
         else
            ptr_assign(i) = -1
         end if
      end do

      deallocate(assignments)

   end subroutine f90_hungarian_algorithm

end module hungarian_mod

! program test_hungarian
!    use hungarian_mod
!    use iso_c_binding
!    implicit none

!    ! Variables for 4x4 Test
!    integer, parameter :: n4 = 4
!    real(f64) :: cost_matrix4(n4, n4)
!    integer :: assignments4(n4)
!    real(f64) :: total_cost4

!    ! Variables for 5x5 Test
!    integer, parameter :: n5 = 5
!    real(f64) :: cost_matrix5(n5, n5)
!    integer :: assignments5(n5)
!    real(f64) :: total_cost5
!    integer :: i

!    print *, "========================================"
!    print *, "   TEST 1: 4x4 Matrix"
!    print *, "========================================"
   
!    cost_matrix4(1, :) = [82.0_f64, 83.0_f64, 69.0_f64, 92.0_f64]
!    cost_matrix4(2, :) = [77.0_f64, 37.0_f64, 49.0_f64, 92.0_f64]
!    cost_matrix4(3, :) = [11.0_f64, 69.0_f64,  5.0_f64, 86.0_f64]
!    cost_matrix4(4, :) = [ 8.0_f64,  9.0_f64, 98.0_f64, 23.0_f64]

!    call hungarian_algorithm(cost_matrix4, n4, assignments4, total_cost4)

!    print *, "Assignments (Index = Row, Value = Col):"
!    do i = 1, n4
!       print *, "  Row", i, " -> Col", assignments4(i), " (Cost: ", cost_matrix4(i, assignments4(i)), ")"
!    end do
!    print *, "Total Cost: ", total_cost4
!    print *, ""

!    print *, "========================================"
!    print *, "   TEST 2: 5x5 Matrix
!    print *, "========================================"
   
!    cost_matrix5(1, :) = [0.1_f64, 0.8_f64, 0.9_f64, 0.7_f64, 0.6_f64]
!    cost_matrix5(2, :) = [0.9_f64, 0.2_f64, 0.8_f64, 0.7_f64, 0.6_f64]
!    cost_matrix5(3, :) = [0.9_f64, 0.8_f64, 0.3_f64, 0.7_f64, 0.6_f64]
!    cost_matrix5(4, :) = [1.0_f64, 1.0_f64, 1.0_f64, 1.0_f64, 1.0_f64]
!    cost_matrix5(5, :) = [1.0_f64, 1.0_f64, 1.0_f64, 1.0_f64, 1.0_f64]

!    call hungarian_algorithm(cost_matrix5, n5, assignments5, total_cost5)

!    print *, "Assignments (Index = Row, Value = Col):"
!    do i = 1, n5
!       print *, "  Row", i, " -> Col", assignments5(i), " (Cost: ", cost_matrix5(i, assignments5(i)), ")"
!    end do
!    print *, "Total Cost: ", total_cost5
!    print *, "========================================"

! end program test_hungarian
