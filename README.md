# Fortran-C++ Hungarian Algorithm Solver

A high-performance implementation of the Hungarian Algorithmm. It is written in Modern Fortran and exposes a C-interface for C/C++ integration.

## Building the Shared Library

Compile the core Fortran module into a shared library. We highly recommend using `-O3` for maximum performance and `-march=x86-64-v3` (or `-march=native`) to leverage modern CPU vectorization instructions (AVX2, etc.).

```bash
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.so
```

---

## Usage: Fortran users

Fortran users can use the module directly without the C-bindings.

call subroutine `hungarian_algorithm` with N x N cost matrix, subroutine will calculate the optimal assignment array and minimum cost.

---

## Usage: C++ users (using Eigen input as an example)

Eigen uses Column-Major memory layout by default, allowing you to pass the `.data()` pointer directly to the Fortran backend without transposing or copying the matrix.

### Example (`main.cpp`)
```cpp
#include <iostream>
#include <Eigen/Dense>

extern "C" {
    void f90_hungarian_algorithm(const double* c_matrix_ptr, int n, int* c_assignments_ptr, double* c_cost_ptr);
}

int main() {
    int n = 4;
    Eigen::MatrixXd cost_matrix(n, n);
    
    cost_matrix << 82, 83, 69, 92,
                   77, 37, 49, 92,
                   11, 69,  5, 86,
                    8,  9, 98, 23;
                    
    int assignments[4];
    double total_cost = 0.0;

    f90_hungarian_algorithm(cost_matrix.data(), n, assignments, &total_cost);

    std::cout << "Total Cost: " << total_cost << "\n";
    for (int i = 0; i < n; ++i) {
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
Worker 0 -> Job 2 (Cost: 69)
Worker 1 -> Job 1 (Cost: 37)
Worker 2 -> Job 0 (Cost: 11)
Worker 3 -> Job 3 (Cost: 23)
```
