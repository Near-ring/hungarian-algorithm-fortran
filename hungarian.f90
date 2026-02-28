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
!>   - The matrix remains non-negative throughout.
!>   - Termination is guaranteed within O(n) adjustment iterations for
!>     exact arithmetic; floating-point rounding may require additional
!>     iterations bounded by the scale-aware tolerance.
!>
!> Memory Ownership (C API):
!>   - The caller owns all buffers passed via pointers.
!>   - c_matrix_ptr must point to n*n contiguous doubles in COLUMN-MAJOR order.
!>   - c_assignments_ptr must point to at least n contiguous int32_t elements.
!>   - c_cost_ptr must point to a single double.
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

   ! Error codes returned by hungarian_algorithm via the info argument:
   !   0 = success
   !   1 = invalid input (n < 0, NaN/Inf in cost matrix)
   !   2 = allocation failure
   !   3 = algorithm did not converge within iteration limit
   !   4 = perfect matching not found (should not occur for valid input)
   integer, parameter, public :: HUNGARIAN_OK              = 0
   integer, parameter, public :: HUNGARIAN_ERR_INVALID     = 1
   integer, parameter, public :: HUNGARIAN_ERR_ALLOC       = 2
   integer, parameter, public :: HUNGARIAN_ERR_NO_CONVERGE = 3
   integer, parameter, public :: HUNGARIAN_ERR_NO_MATCH    = 4

contains

   ! ------------------------------------------------------------------
   ! Subroutine: Maximum Bipartite Matching
   ! Purpose:    Finds the maximum number of matches in a bipartite graph.
   ! Algorithm:  Iterative Depth-First Search (DFS) to find augmenting paths.
   ! ------------------------------------------------------------------
   subroutine maximum_bipartite_matching(adj, match_r_to_c, info)
      implicit none
      ! --- Inputs ---
      logical, intent(in) :: adj(:,:)            ! Adjacency matrix (True where cost is 0)

      ! --- Outputs ---
      integer, intent(out) :: match_r_to_c(:)    ! Maps Row Index -> Matched Col Index (0 if unmatched)
      integer, intent(out) :: info               ! 0 = success, nonzero = allocation failure

      ! --- Local Variables ---
      integer :: n
      integer, allocatable :: match_c_to_r(:)    ! Inverse mapping: Col -> matched Row
      integer :: i, u, v, r
      integer :: path_end_c, curr_c, prev_c_match

      ! Explicit stack for DFS (depth bounded by n, not n*n)
      integer, allocatable :: stack(:)
      logical, allocatable :: visited_cols(:)
      integer, allocatable :: parent_row_of_col(:)  ! parent_row_of_col(col) = row that discovered col
      integer :: stack_top, istat

      info = 0
      n = size(adj, 1)

      allocate(match_c_to_r(n), stack(n), visited_cols(n), parent_row_of_col(n), stat=istat)
      if (istat /= 0) then
         info = HUNGARIAN_ERR_ALLOC
         match_r_to_c = 0
         return
      end if

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
                     ! If matched, continue DFS from the row matching 'v'
                     ! Stack depth bounded by n (each row pushed at most once per DFS)
                     if (stack_top < n) then
                        stack_top = stack_top + 1
                        stack(stack_top) = match_c_to_r(v)
                     end if
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

      deallocate(match_c_to_r, stack, visited_cols, parent_row_of_col)
   end subroutine maximum_bipartite_matching

   ! ------------------------------------------------------------------
   ! Subroutine: Konig's Theorem (Minimum Line Cover)
   ! Purpose:    Finds the minimum number of lines (rows + cols) needed
   !             to cover all zeros in the matrix.
   ! Logic:      Min Vertex Cover in Bipartite Graph = Max Matching.
   !             1. Calculate Max Matching.
   !             2. Use BFS to mark reachable nodes from unmatched rows.
   !             3. Cover = (Unmarked Rows) + (Marked Cols).
   ! ------------------------------------------------------------------
   subroutine minimum_line_cover(zeros_mat, lines_rows, lines_cols, info)
      implicit none
      ! --- Inputs ---
      logical, intent(in) :: zeros_mat(:,:)      ! Boolean matrix where True represents a Zero

      ! --- Outputs ---
      logical, intent(out) :: lines_rows(:)      ! True if row is part of the cover
      logical, intent(out) :: lines_cols(:)       ! True if col is part of the cover
      integer, intent(out) :: info               ! 0 = success, nonzero = error

      ! --- Local Variables ---
      integer :: n
      integer, allocatable :: match_r_to_c(:), match_c_to_r(:)
      logical, allocatable :: marked_rows(:), marked_cols(:)
      integer, allocatable :: queue(:)
      integer :: q_head, q_tail, r, c, row, col, matched_row, istat

      info = 0
      n = size(zeros_mat, 1)

      allocate(match_r_to_c(n), match_c_to_r(n), marked_rows(n), marked_cols(n), &
               queue(n), stat=istat)
      if (istat /= 0) then
         info = HUNGARIAN_ERR_ALLOC
         lines_rows = .false.
         lines_cols = .false.
         return
      end if

      ! 1. Find Max Matching on the zero structure
      call maximum_bipartite_matching(zeros_mat, match_r_to_c, info)
      if (info /= 0) then
         deallocate(match_r_to_c, match_c_to_r, marked_rows, marked_cols, queue)
         lines_rows = .false.
         lines_cols = .false.
         return
      end if

      ! Build inverse lookup (Col -> Row)
      match_c_to_r = 0
      do r = 1, n
         c = match_r_to_c(r)
         if (c /= 0) match_c_to_r(c) = r
      end do

      ! 2. BFS for Konig's Construction
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

      deallocate(match_r_to_c, match_c_to_r, marked_rows, marked_cols, queue)
   end subroutine minimum_line_cover

   ! ==================================================================
   ! MAIN SOLVER: HUNGARIAN ALGORITHM
   !
   ! Solves the linear sum assignment problem for square n x n matrices.
   ! Finds the assignment that minimizes the total cost.
   !
   ! Arguments:
   !   cost_matrix  - Input n x n cost matrix (assumed-shape)
   !   assignments  - Output array of length n: assignments(i) = j means
   !                  row i is assigned to column j. 0 if unassigned.
   !   total_cost   - Output: sum of assigned costs from cost_matrix
   !   info         - Output error code:
   !                  0 = success (HUNGARIAN_OK)
   !                  1 = invalid input (HUNGARIAN_ERR_INVALID)
   !                  2 = allocation failure (HUNGARIAN_ERR_ALLOC)
   !                  3 = no convergence (HUNGARIAN_ERR_NO_CONVERGE)
   !                  4 = no perfect matching (HUNGARIAN_ERR_NO_MATCH)
   ! ==================================================================
   subroutine hungarian_algorithm(cost_matrix, assignments, total_cost, info)
      implicit none
      real(f64), intent(in) :: cost_matrix(:,:)
      integer, intent(out) :: assignments(:)
      real(f64), intent(out) :: total_cost
      integer, intent(out) :: info

      real(f64), allocatable :: matrix(:,:)
      logical, allocatable :: zeros_bool(:,:)
      logical, allocatable :: lines_rows(:), lines_cols(:)
      integer :: n, i, j, istat
      real(f64) :: k, tol
      integer :: num_lines
      integer(8) :: max_iter, iteration

      ! Initialize outputs
      info = 0
      total_cost = 0.0_f64

      n = size(cost_matrix, 1)

      ! --- Input Validation (M-3: define n=0 behavior, H-3: NaN/Inf) ---
      if (n == 0) then
         ! n=0 is valid: trivially no assignments to make
         return
      end if

      if (n < 0) then
         info = HUNGARIAN_ERR_INVALID
         return
      end if

      ! Validate cost_matrix is square
      if (size(cost_matrix, 2) /= n) then
         info = HUNGARIAN_ERR_INVALID
         assignments = 0
         return
      end if

      ! Check for NaN/Inf in cost matrix
      if (any(.not. ieee_is_finite(cost_matrix))) then
         info = HUNGARIAN_ERR_INVALID
         assignments = 0
         return
      end if

      ! --- Allocation with stat= checks (M-1) ---
      allocate(matrix(n, n), zeros_bool(n, n), lines_rows(n), lines_cols(n), stat=istat)
      if (istat /= 0) then
         info = HUNGARIAN_ERR_ALLOC
         assignments = 0
         return
      end if

      ! Create working copy
      matrix = cost_matrix

      ! Compute scale-aware tolerance (H-2)
      tol = 1.0d-12 * max(1.0_f64, maxval(abs(cost_matrix)))

      ! --- Step 1: Subtract Row Minima ---
      do i = 1, n
         matrix(i, :) = matrix(i, :) - minval(matrix(i, :))
      end do

      ! --- Step 2: Subtract Column Minima ---
      do j = 1, n
         matrix(:, j) = matrix(:, j) - minval(matrix(:, j))
      end do

      ! --- Step 3: Iterative Reduction ---
      ! Use int64 for iteration limit to avoid overflow (M-4)
      max_iter = int(n, 8) * 20_8
      iteration = 0

      do while (iteration < max_iter)
         iteration = iteration + 1

         ! Identify zeros using scale-aware tolerance (H-2)
         zeros_bool = (abs(matrix) <= tol)

         ! Check if we have enough zeros to form a complete assignment
         call minimum_line_cover(zeros_bool, lines_rows, lines_cols, info)
         if (info /= 0) then
            assignments = 0
            total_cost = 0.0_f64
            deallocate(matrix, zeros_bool, lines_rows, lines_cols)
            return
         end if

         num_lines = count(lines_rows) + count(lines_cols)

         ! Optimization found?
         if (num_lines == n) exit

         ! --- Step 4: Matrix Update ---
         ! Directly compute min uncovered value and apply adjustments
         ! without n×n temporary arrays (H-11, M-6, M-7)
         k = huge(1.0_f64)
         do j = 1, n
            if (lines_cols(j)) cycle
            do i = 1, n
               if (lines_rows(i)) cycle
               k = min(k, matrix(i, j))
            end do
         end do

         ! Subtract k from uncovered, add k to doubly covered
         do j = 1, n
            do i = 1, n
               if (.not. lines_rows(i) .and. .not. lines_cols(j)) then
                  matrix(i, j) = matrix(i, j) - k
               else if (lines_rows(i) .and. lines_cols(j)) then
                  matrix(i, j) = matrix(i, j) + k
               end if
            end do
         end do
      end do

      ! Check convergence (C-3)
      if (iteration >= max_iter .and. num_lines /= n) then
         info = HUNGARIAN_ERR_NO_CONVERGE
         assignments = 0
         total_cost = 0.0_f64
         deallocate(matrix, zeros_bool, lines_rows, lines_cols)
         return
      end if

      ! --- Step 5: Final Assignment ---
      zeros_bool = (abs(matrix) <= tol)
      call maximum_bipartite_matching(zeros_bool, assignments, info)
      if (info /= 0) then
         assignments = 0
         total_cost = 0.0_f64
         deallocate(matrix, zeros_bool, lines_rows, lines_cols)
         return
      end if

      ! Validate perfect matching (H-1)
      if (count(assignments /= 0) /= n) then
         info = HUNGARIAN_ERR_NO_MATCH
         total_cost = 0.0_f64
         deallocate(matrix, zeros_bool, lines_rows, lines_cols)
         return
      end if

      ! Calculate Total Cost based on original matrix
      total_cost = 0.0_f64
      do i = 1, n
         total_cost = total_cost + cost_matrix(i, assignments(i))
      end do

      deallocate(matrix, zeros_bool, lines_rows, lines_cols)

   end subroutine hungarian_algorithm

   ! ==================================================================
   ! C WRAPPER
   !
   ! NOTE on Memory Layout:
   !   C is Row-Major, Fortran is Column-Major.
   !   The input matrix MUST be in COLUMN-MAJOR order. If using row-major
   !   C arrays, transpose before calling. Eigen matrices are column-major
   !   by default and can be passed directly via .data().
   !
   ! Memory Ownership Contract (M-11):
   !   - Caller owns all buffers. They must remain valid during the call.
   !   - c_matrix_ptr: n*n contiguous doubles, column-major layout.
   !   - c_assignments_ptr: at least n contiguous int32_t elements.
   !   - c_cost_ptr: pointer to a single double.
   !   - Not thread-safe for shared buffers.
   !
   ! Return value via info pointer:
   !   0 = success, nonzero = error (see HUNGARIAN_ERR_* constants)
   ! ==================================================================
   subroutine f90_hungarian_algorithm(c_matrix_ptr, n, c_assignments_ptr, c_cost_ptr, c_info_ptr) &
      bind(C, name="f90_hungarian_algorithm")
      implicit none
      type(c_ptr), value, intent(in) :: c_matrix_ptr
      integer(c_int32_t), value, intent(in) :: n
      type(c_ptr), value, intent(in) :: c_assignments_ptr
      type(c_ptr), value, intent(in) :: c_cost_ptr
      type(c_ptr), value, intent(in) :: c_info_ptr

      real(c_double), pointer :: ptr_matrix(:,:)
      integer(c_int32_t), pointer :: ptr_assign(:)
      real(c_double), pointer :: ptr_cost
      integer(c_int32_t), pointer :: ptr_info

      integer, allocatable :: assignments(:)
      real(f64) :: cost
      integer :: i, n_f, info, istat

      ! Convert c_int32_t to default integer for internal use (C-1)
      n_f = int(n, kind(n_f))

      ! Validate pointers (C-2)
      if (.not. c_associated(c_matrix_ptr) .or. &
          .not. c_associated(c_assignments_ptr) .or. &
          .not. c_associated(c_cost_ptr) .or. &
          .not. c_associated(c_info_ptr)) then
         ! Cannot write error code if info pointer is null, just return
         if (c_associated(c_info_ptr)) then
            call c_f_pointer(c_info_ptr, ptr_info)
            ptr_info = int(HUNGARIAN_ERR_INVALID, c_int32_t)
         end if
         return
      end if

      ! Map info output first so we can report errors
      call c_f_pointer(c_info_ptr, ptr_info)

      ! Handle n <= 0
      if (n_f <= 0) then
         if (n_f == 0) then
            ptr_info = int(HUNGARIAN_OK, c_int32_t)
         else
            ptr_info = int(HUNGARIAN_ERR_INVALID, c_int32_t)
         end if
         return
      end if

      ! Map C pointers to Fortran arrays using default integer for shape (M-5)
      call c_f_pointer(c_matrix_ptr, ptr_matrix, [n_f, n_f])
      call c_f_pointer(c_assignments_ptr, ptr_assign, [n_f])
      call c_f_pointer(c_cost_ptr, ptr_cost)

      allocate(assignments(n_f), stat=istat)
      if (istat /= 0) then
         ptr_info = int(HUNGARIAN_ERR_ALLOC, c_int32_t)
         return
      end if

      ! Solve
      call hungarian_algorithm(ptr_matrix, assignments, cost, info)

      ! Report error code
      ptr_info = int(info, c_int32_t)

      if (info == HUNGARIAN_OK) then
         ! Return Results
         ptr_cost = cost
         do i = 1, n_f
            if (assignments(i) /= 0) then
               ! Convert 1-based Fortran index to 0-based C index
               ptr_assign(i) = int(assignments(i) - 1, c_int32_t)
            else
               ptr_assign(i) = -1_c_int32_t
            end if
         end do
      end if

      deallocate(assignments)

   end subroutine f90_hungarian_algorithm

end module hungarian_mod
