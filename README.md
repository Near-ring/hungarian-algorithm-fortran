# Hungarian Algorithm — Modern Fortran

[![CI](https://github.com/chuantian/hungarian-algorithm-fortran/actions/workflows/ci.yml/badge.svg)](https://github.com/chuantian/hungarian-algorithm-fortran/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[中文文档](README.zh-CN.md)

A high-performance implementation of the **Hungarian Algorithm** (Kuhn-Munkres) for solving the linear sum assignment problem. Written in Modern Fortran with a C-compatible interface for seamless C/C++ integration.

## Features

- Solves N×N cost matrix assignment in O(n³) time
- C-interoperable API via `iso_c_binding` — call directly from C, C++, Python (ctypes), etc.
- Robust error handling with typed error codes
- Input validation: rejects NaN, Inf, and invalid dimensions
- Scale-aware floating-point tolerance for numerical stability
- Zero external dependencies — single-file Fortran module

## Prerequisites

- **Fortran compiler**: `gfortran` ≥ 7.0 (or any Fortran 2008-compliant compiler)
- **Make** (optional, for build automation)
- **Eigen** (optional, only for the C++ example)

## Quick Start

```bash
# Build the shared library
make

# Build and run tests
make test

# Clean build artifacts
make clean
```

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
# Linux
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.so

# macOS
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.dylib

# Windows (MinGW)
gfortran -shared -O3 hungarian.f90 -o libhungarian.dll
```

### Tests

```bash
make test
```

This compiles `hungarian.f90` and `test_hungarian.f90` together, then runs the test suite. Expected output:

```text
=======================================
 Hungarian Algorithm Test Suite
=======================================
  PASS: 4x4 README example (cost=140)
  PASS: 1x1 matrix
  PASS: 2x2 simple (cost=5)
  PASS: 3x3 known (cost=17)
  PASS: 5x5 known (cost=72)
  PASS: n=0 (trivial, no-op)
  PASS: Negative costs (cost=-20)
  PASS: Identical costs (cost=21)
  PASS: Large costs (1e12 scale)
  PASS: Zero matrix (cost=0)
  PASS: NaN input rejected
  PASS: Inf input rejected
  PASS: Brute force 2x2
  PASS: Brute force 3x3
  PASS: Brute force 4x4
=======================================
 Results: 15/15 passed
 ALL TESTS PASSED
=======================================
```

## API Reference

### Fortran API

```fortran
use hungarian_mod

call hungarian_algorithm(cost_matrix, assignments, total_cost, info)
```

| Argument      | Type                 | Intent | Description                                                           |
|---------------|----------------------|--------|-----------------------------------------------------------------------|
| `cost_matrix` | `real(f64) :: (:,:)` | `in`   | N×N cost matrix                                                       |
| `assignments` | `integer :: (:)`     | `out`  | Assignment result: `assignments(i) = j` means row i → column j        |
| `total_cost`  | `real(f64)`          | `out`  | Sum of assigned costs                                                 |
| `info`        | `integer`            | `out`  | Error code (see below)                                                |

**Example:**

```fortran
use hungarian_mod

real(f64) :: cost_matrix(4, 4)
integer :: assignments(4), info
real(f64) :: total_cost

cost_matrix(1, :) = [82.0_f64, 83.0_f64, 69.0_f64, 92.0_f64]
cost_matrix(2, :) = [77.0_f64, 37.0_f64, 49.0_f64, 92.0_f64]
cost_matrix(3, :) = [11.0_f64, 69.0_f64,  5.0_f64, 86.0_f64]
cost_matrix(4, :) = [ 8.0_f64,  9.0_f64, 98.0_f64, 23.0_f64]

call hungarian_algorithm(cost_matrix, assignments, total_cost, info)

if (info /= HUNGARIAN_OK) then
   print *, 'Error:', info
   stop 1
end if

print *, 'Total cost:', total_cost  ! 140.0
```

### C/C++ API

```c
void f90_hungarian_algorithm(
    const double* c_matrix_ptr,     // N×N cost matrix (column-major!)
    int32_t       n,                // Matrix dimension
    int32_t*      c_assignments_ptr, // Output: N assignments (0-based indices)
    double*       c_cost_ptr,       // Output: total cost
    int32_t*      c_info_ptr        // Output: error code
);
```

> **⚠️ Column-Major Order Required:** The matrix must be in column-major layout (Fortran convention). Eigen matrices are column-major by default. For row-major C arrays, transpose before calling or store as `matrix[col * n + row]`.

**Example (`main.cpp`):**

```cpp
#include <cstdint>
#include <iostream>
#include <Eigen/Dense>

extern "C" {
    void f90_hungarian_algorithm(const double* c_matrix_ptr, int32_t n,
                                 int32_t* c_assignments_ptr, double* c_cost_ptr,
                                 int32_t* c_info_ptr);
}

int main() {
    int32_t n = 4;
    Eigen::MatrixXd cost(n, n);

    cost << 82, 83, 69, 92,
            77, 37, 49, 92,
            11, 69,  5, 86,
             8,  9, 98, 23;

    int32_t assignments[4], info = 0;
    double total_cost = 0.0;

    f90_hungarian_algorithm(cost.data(), n, assignments, &total_cost, &info);

    if (info != 0) {
        std::cerr << "Error: " << info << "\n";
        return 1;
    }

    std::cout << "Total Cost: " << total_cost << "\n";
    for (int32_t i = 0; i < n; ++i)
        std::cout << "Worker " << i << " -> Job " << assignments[i] << "\n";
}
```

```bash
g++ -O3 main.cpp -L. -lhungarian -Wl,-rpath,. -o example && ./example
```

**Output:**

```text
Total Cost: 140
Worker 0 -> Job 2
Worker 1 -> Job 1
Worker 2 -> Job 0
Worker 3 -> Job 3
```

### Error Codes

| Constant                     | Value | Description                                             |
|------------------------------|-------|---------------------------------------------------------|
| `HUNGARIAN_OK`               | 0     | Success                                                 |
| `HUNGARIAN_ERR_INVALID`      | 1     | Invalid input (n < 0, NaN/Inf in matrix, null pointer)  |
| `HUNGARIAN_ERR_ALLOC`        | 2     | Memory allocation failure                               |
| `HUNGARIAN_ERR_NO_CONVERGE`  | 3     | Algorithm did not converge within iteration limit       |
| `HUNGARIAN_ERR_NO_MATCH`     | 4     | Perfect matching not found                              |

## Algorithm

This implementation uses the **reduction + minimum line cover + adjustment** variant:

1. **Row reduction** — subtract each row's minimum
2. **Column reduction** — subtract each column's minimum
3. **Minimum line cover** — find minimum lines to cover all zeros (via König's theorem: max bipartite matching = min vertex cover)
4. **Adjustment** — if cover < n: subtract minimum uncovered value from uncovered cells, add it to doubly-covered cells
5. Repeat steps 3–4 until cover = n
6. **Extract assignment** from the final zero structure

Convergence is guaranteed within O(n) adjustment iterations for well-conditioned input. The implementation uses a scale-aware tolerance (`1e-12 × max|cost|`) for floating-point stability.

## Project Structure

```text
hungarian-algorithm-fortran/
├── hungarian.f90          # Core module: algorithm + C wrapper
├── test_hungarian.f90     # Test suite (15 tests)
├── Makefile               # Cross-platform build system
├── .github/
│   └── workflows/
│       └── ci.yml         # GitHub Actions CI
├── LICENSE                # MIT License
└── README.md
```

## License

[MIT](LICENSE)
