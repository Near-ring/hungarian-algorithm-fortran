# Hungarian Algorithm — Modern Fortran

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[中文文档](README.zh-CN.md)

A high-performance implementation of the **Hungarian Algorithm** (Kuhn-Munkres) for solving the linear sum assignment problem. Written in Modern Fortran with a C-compatible interface for C/C++ integration.

## Features

- Solves $N \times N$ cost matrix assignment in $O(n^3)$ time.
- C-interoperable API via `iso_c_binding` — call directly from C, C++, Python (ctypes), etc.
- Input validation: On invalid inputs (e.g., $N \le 0$, `NaN`/`Inf` in cost matrix), returns `NaN` for costs and `0`/`-1` for assignments .
- Zero external dependencies — single-file Fortran module.

## Prerequisites

- **Fortran compiler**: `gfortran` ≥ 7.0 (or any Fortran 2008-compliant compiler)


## Building

### Shared Library

The `Makefile` auto-detects your platform and builds the appropriate shared library:

| Platform        | Output               |
|-----------------|----------------------|
| Linux           | `libhungarian.so`    |
| macOS           | `libhungarian.dylib` |
| Windows (MinGW) | `libhungarian.dll`   |

**Manual build** (without Make):

```bash
# Linux (and wsl)
gfortran -shared -fPIC -O3 -march=native hungarian_algorithm.f90 -o libhungarian.so

# macOS
gfortran -shared -fPIC -O3 -march=native hungarian_algorithm.f90 -o libhungarian.dylib

# Windows (via MinGW)
gfortran -shared -O3 hungarian_algorithm.f90 -o libhungarian.dll
```

### Tests

```bash
make test
```

## API Reference

### Fortran API

```fortran
use hungarian_mod

call hungarian_algorithm(cost_matrix, n, assignments, total_cost)
```

| Argument      | Type                 | Intent | Description                                                           |
|---------------|----------------------|--------|-----------------------------------------------------------------------|
| `cost_matrix` | `real(f64) :: (n,n)` | `in`   | $N \times N$ cost matrix                                              |
| `n`           | `integer`            | `in`   | Dimension of the cost matrix                                          |
| `assignments` | `integer :: (n)`     | `out`  | Assignment result: `assignments(i) = j` means row $i \to$ column $j$  |
| `total_cost`  | `real(f64)`          | `out`  | Sum of assigned costs (`NaN` if input is invalid)                     |

**Example:**

```fortran
use hungarian_mod
use, intrinsic :: ieee_arithmetic

integer, parameter :: n = 4
real(f64) :: cost_matrix(n, n)
integer :: assignments(n)
real(f64) :: total_cost
integer :: i

cost_matrix(1, :) = [82.0_f64, 83.0_f64, 69.0_f64, 92.0_f64]
cost_matrix(2, :) = [77.0_f64, 37.0_f64, 49.0_f64, 92.0_f64]
cost_matrix(3, :) = [11.0_f64, 69.0_f64,  5.0_f64, 86.0_f64]
cost_matrix(4, :) = [ 8.0_f64,  9.0_f64, 98.0_f64, 23.0_f64]

call hungarian_algorithm(cost_matrix, n, assignments, total_cost)

if (.not. ieee_is_finite(total_cost)) then
   print *, 'Error: Invalid cost matrix provided.'
   stop 1
end if

print *, 'Total cost:', total_cost  ! 140.0
```

---

### C/C++ API

```c
void f90_hungarian_algorithm(
    const double* c_matrix_ptr,      // N×N cost matrix (column-major layout)
    int32_t       n,                 // Matrix dimension
    int32_t* c_assignments_ptr, // Output: N assignments (0-based indices, -1 on error)
    double* c_cost_ptr         // Output: total cost (NaN on error)
);
```

> [!IMPORTANT]
> The matrix must be in **column-major** layout (Fortran convention). Eigen matrices are column-major by default. For row-major C arrays, transpose the matrix before calling or store data as `matrix[col * n + row]`.

**Example (`main.cpp`):**

```cpp
#include <cstdint>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

extern "C" {
    void f90_hungarian_algorithm(const double* c_matrix_ptr, int32_t n,
                                 int32_t* c_assignments_ptr, double* c_cost_ptr);
}

int main() {
    int32_t n = 4;
    Eigen::MatrixXd cost(n, n);

    cost << 82, 83, 69, 92,
            77, 37, 49, 92,
            11, 69,  5, 86,
             8,  9, 98, 23;

    int32_t assignments[4];
    double total_cost = 0.0;

    f90_hungarian_algorithm(cost.data(), n, assignments, &total_cost);

    if (std::isnan(total_cost)) {
        std::cerr << "Error: Invalid cost matrix or dimension.\n";
        return 1;
    }

    std::cout << "Total Cost: " << total_cost << "\n";
    for (int32_t i = 0; i < n; ++i)
        std::cout << "Worker " << i << " -> Job " << assignments[i] << "\n";
    
    return 0;
}
```

## Algorithm

This implementation uses the **reduction + minimum line cover + adjustment** variant:

1. **Row reduction** — subtract each row's minimum.
2. **Column reduction** — subtract each column's minimum.
3. **Minimum line cover** — find minimum lines to cover all zeros. 
   
   *(This leverages König's theorem: the size of the maximum bipartite matching equals the size of the minimum vertex cover).*
4. **Adjustment** — if the number of covering lines is $< n$: subtract the minimum uncovered value from all uncovered cells, and add it to all doubly-covered cells.
5. Repeat steps 3–4 until cover $= n$.
6. **Extract assignment** from the final zero structure using Depth-First Search (DFS) for augmenting paths.

Convergence is guaranteed within $O(n)$ adjustment iterations for well-conditioned input. The implementation uses a scale-aware tolerance based on machine epsilon (`epsilon(1.0_f64) * max(1.0_f64, maxval(abs(cost_matrix)))`) for optimal floating-point stability.
