# Fortran-C++ Hungarian Algorithm Solver

A high-performance implementation of the Hungarian Algorithm. It is written in Modern Fortran and exposes a C-interface for C/C++ integration.

## Building the Shared Library

Compile the core Fortran module into a shared library. We highly recommend using `-O3` for maximum performance and `-march=native` to leverage modern CPU vectorization instructions (AVX2, etc.).

**Linux:**
```bash
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.so
```

**macOS:**
```bash
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.dylib
```

**Windows (MinGW):**
```bash
gfortran -shared -O3 hungarian.f90 -o libhungarian.dll
```

Alternatively, use the provided `Makefile`:
```bash
make        # builds shared library for current platform
make test   # builds and runs tests
make clean  # removes build artifacts
```

---

## Usage: Fortran users

Fortran users can use the module directly without the C-bindings.

Call subroutine `hungarian_algorithm` with an N × N cost matrix. The subroutine returns the optimal assignment array, minimum cost, and an error code (`info`):

```fortran
use hungarian_mod

real(f64) :: cost_matrix(4, 4)
integer :: assignments(4), info
real(f64) :: total_cost

! ... fill cost_matrix ...

call hungarian_algorithm(cost_matrix, assignments, total_cost, info)

if (info /= HUNGARIAN_OK) then
   ! Handle error: HUNGARIAN_ERR_INVALID, HUNGARIAN_ERR_ALLOC,
   !               HUNGARIAN_ERR_NO_CONVERGE, HUNGARIAN_ERR_NO_MATCH
end if
```

---

## Usage: C/C++ users

**⚠️ Important: Column-Major Order Required**

The Fortran backend uses column-major memory layout. You must pass your matrix in **column-major order**:
- **Eigen** (C++): Default layout is column-major — pass `.data()` directly.
- **Plain C/C++ arrays** (row-major): You must transpose the matrix before calling, or store it in column-major order (i.e., `matrix[col * n + row]`).

### Example (`main.cpp`) — using Eigen (column-major by default)
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
    Eigen::MatrixXd cost_matrix(n, n);
    
    cost_matrix << 82, 83, 69, 92,
                   77, 37, 49, 92,
                   11, 69,  5, 86,
                    8,  9, 98, 23;
                    
    int32_t assignments[4];
    double total_cost = 0.0;
    int32_t info = 0;

    f90_hungarian_algorithm(cost_matrix.data(), n, assignments, &total_cost, &info);

    if (info != 0) {
        std::cerr << "Error: " << info << "\n";
        return 1;
    }

    std::cout << "Total Cost: " << total_cost << "\n";
    for (int32_t i = 0; i < n; ++i) {
        std::cout << "Worker " << i << " -> Job " << assignments[i] << "\n";
    }

    return 0;
}
```

**Compile and Run:**
```bash
g++ -O3 -march=native main.cpp -I /usr/include/eigen3 -L. -lhungarian -Wl,-rpath=. -o test_cpp
./test_cpp
```

---

### Expected Output from the above example should be:
```text
Total Cost: 140
Worker 0 -> Job 2
Worker 1 -> Job 1
Worker 2 -> Job 0
Worker 3 -> Job 3
```
