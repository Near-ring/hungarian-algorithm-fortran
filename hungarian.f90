module hungarian_mod
   use iso_c_binding
   implicit none

   ! Define precision constants
   integer, parameter :: f64 = C_DOUBLE

contains

   ! ------------------------------------------------------------------
   ! Subroutine: Maximum Bipartite Matching
   ! Purpose:    Finds the maximum number of matches in a bipartite graph.
   ! Algorithm:  Iterative Depth-First Search (DFS) to find augmenting paths.
   ! ------------------------------------------------------------------
   subroutine maximum_bipartite_matching(adj, n, match_r_to_c)
      ! --- Inputs ---
      logical, intent(in) :: adj(n, n)      ! Adjacency matrix (True where cost is 0)
      integer, intent(in) :: n              ! Dimension of the matrix

      ! --- Outputs ---
      integer, intent(out) :: match_r_to_c(n) ! Maps Row Index -> Matched Col Index (0 if unmatched)

      ! --- Local Variables ---
      integer :: match_c_to_r(n)
      integer :: i, u, v, r
      integer :: path_end_c, curr_c, prev_c_match

      ! Explicit stack for DFS to prevent stack overflow on large N
      integer, allocatable :: stack(:)
      logical, allocatable :: visited_cols(:)
      integer, allocatable :: parent(:)
      integer :: stack_top

      allocate(stack(n * n), visited_cols(n), parent(n))

      match_r_to_c = 0
      match_c_to_r = 0

      ! Try to find an augmenting path for every unmatched row 'i'
      do i = 1, n
         if (match_r_to_c(i) /= 0) cycle

         ! Initialize DFS
         stack_top = 1
         stack(stack_top) = i
         visited_cols = .false.
         parent = 0
         path_end_c = 0

         ! Iterative DFS Loop
         dfs_loop: do while (stack_top > 0)
            u = stack(stack_top)
            stack_top = stack_top - 1

            do v = 1, n
               ! If edge exists (zero in cost matrix) and column unvisited
               if (adj(u, v) .and. (.not. visited_cols(v))) then
                  visited_cols(v) = .true.
                  parent(v) = u

                  ! If column 'v' is unmatched, we found an augmenting path
                  if (match_c_to_r(v) == 0) then
                     path_end_c = v
                     exit dfs_loop
                  else
                     ! If matched, continue DFS from the row matching 'v'
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
               r = parent(curr_c)
               prev_c_match = match_r_to_c(r)

               match_c_to_r(curr_c) = r
               match_r_to_c(r) = curr_c

               if (r == i) exit
               curr_c = prev_c_match
            end do
         end if
      end do

      deallocate(stack, visited_cols, parent)
   end subroutine maximum_bipartite_matching

   ! ------------------------------------------------------------------
   ! Subroutine: Konig's Theorem (Minimum Line Cover)
   ! Purpose:    Finds the minimum number of lines (rows + cols) needed
   !             to cover all zeros in the matrix.
   ! Logic:      Min Vertex Cover in Bipartite Graph = Max Matching.
   !             1. Calculate Max Matching.
   !             2. Use BFS/DFS to mark reachable nodes from unmatched rows.
   !             3. Cover = (Unmarked Rows) + (Marked Cols).
   ! ------------------------------------------------------------------
   subroutine minimum_line_cover(zeros_mat, n, lines_rows, lines_cols)
      ! --- Inputs ---
      logical, intent(in) :: zeros_mat(n, n) ! Boolean matrix where True represents a Zero
      integer, intent(in) :: n

      ! --- Outputs ---
      logical, intent(out) :: lines_rows(n)  ! True if row is part of the cover
      logical, intent(out) :: lines_cols(n)  ! True if col is part of the cover

      ! --- Local Variables ---
      integer :: match_r_to_c(n), match_c_to_r(n)
      logical :: marked_rows(n), marked_cols(n)
      integer, allocatable :: queue(:)
      integer :: q_head, q_tail, r, c, row, col, matched_row

      ! 1. Find Max Matching on the zero structures
      call maximum_bipartite_matching(zeros_mat, n, match_r_to_c)

      ! Build inverse lookup (Col -> Row)
      match_c_to_r = 0
      do r = 1, n
         c = match_r_to_c(r)
         if (c /= 0) match_c_to_r(c) = r
      end do

      ! 2. BFS for Konig's Construction
      allocate(queue(n * n))
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
            ! Visit columns connected by a Zero edge
            if (zeros_mat(row, col) .and. (.not. marked_cols(col))) then
               marked_cols(col) = .true.

               matched_row = match_c_to_r(col)
               ! If this column is matched to a row, add that row to queue
               if (matched_row /= 0 .and. (.not. marked_rows(matched_row))) then
                  marked_rows(matched_row) = .true.
                  q_tail = q_tail + 1
                  queue(q_tail) = matched_row
               end if
            end if
         end do
      end do

      ! 3. Construct Cover
      ! Rows are covered if they are NOT marked (reachable from unmatched).
      ! Cols are covered if they ARE marked.
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
      integer, intent(out) :: assignments(n)  ! assignment(i) = j means Row i assigned to Col j
      real(f64), intent(out) :: total_cost

      real(f64), allocatable :: matrix(:,:)
      integer :: i, j
      real(f64) :: k
      integer :: num_lines, iteration
      logical, allocatable :: zeros_bool(:,:)
      logical, allocatable :: lines_rows(:), lines_cols(:)
      integer, allocatable :: line_viz(:,:)
      logical, allocatable :: uncovered_mask(:,:), doubly_covered_mask(:,:)

      allocate(matrix(n, n))
      allocate(zeros_bool(n, n))
      allocate(lines_rows(n), lines_cols(n))
      allocate(line_viz(n, n))
      allocate(uncovered_mask(n, n), doubly_covered_mask(n, n))

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
      iteration = 0
      do while (iteration < n*20) ! Safety break
         iteration = iteration + 1

         ! Identify zeros (using Epsilon for float stability)
         zeros_bool = (abs(matrix) < epsilon(1.0d0))

         ! Check if we have enough zeros to form a complete assignment
         call minimum_line_cover(zeros_bool, n, lines_rows, lines_cols)

         num_lines = count(lines_rows) + count(lines_cols)

         ! Optimization found?
         if (num_lines == n) exit

         ! --- Step 4: Matrix Update ---
         ! Build a map of covered lines
         ! 0 = Uncovered, 1 = Covered once, 2 = Covered twice (intersection)
         line_viz = 0
         do i = 1, n
            if (lines_rows(i)) line_viz(i, :) = line_viz(i, :) + 1
            if (lines_cols(i)) line_viz(:, i) = line_viz(:, i) + 1
         end do

         uncovered_mask = (line_viz == 0)

         ! Find minimum value 'k' in uncovered area
         k = huge(k)
         do i = 1, n
            do j = 1, n
               if (uncovered_mask(i, j)) then
                  if (matrix(i, j) < k) k = matrix(i, j)
               end if
            end do
         end do

         ! Subtract k from uncovered elements
         where (uncovered_mask)
            matrix = matrix - k
         end where

         ! Add k to doubly covered elements (intersections of lines)
         doubly_covered_mask = (line_viz == 2)
         where (doubly_covered_mask)
            matrix = matrix + k
         end where
      end do

      ! --- Step 5: Final Assignment ---
      zeros_bool = (abs(matrix) < epsilon(1.0d0))
      call maximum_bipartite_matching(zeros_bool, n, assignments)

      ! Calculate Total Cost based on original matrix
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
   ! NOTE on Memory Layout:
   ! C is Row-Major, Fortran is Column-Major.
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
      integer :: i

      ! 1. Zero-Copy Map (C array viewed as Fortran array)
      call c_f_pointer(c_matrix_ptr, ptr_matrix, [n, n])

      ! 2. Map Outputs
      call c_f_pointer(c_assignments_ptr, ptr_assign, [n])
      call c_f_pointer(c_cost_ptr, ptr_cost)

      allocate(assignments(n))

      ! 3. Solve
      call hungarian_algorithm(ptr_matrix, n, assignments, cost)

      ! 4. Return Results
      ptr_cost = cost
      do i = 1, n
         if (assignments(i) /= 0) then
            ! Convert 1-based Fortran index to 0-based C index
            ptr_assign(i) = assignments(i) - 1
         else
            ptr_assign(i) = -1
         end if
      end do

      deallocate(assignments)

   end subroutine f90_hungarian_algorithm

end module hungarian_mod
