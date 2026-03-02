# 匈牙利算法 — Modern Fortran 实现

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[English](README.md)

基于 **匈牙利算法**（Kuhn-Munkres）的高性能线性和分配问题求解器。使用 Modern Fortran 编写，通过 C 兼容接口实现与 C/C++ 的无缝集成。

## 特性

- O(n³) 时间复杂度求解 $N \times N$ 代价矩阵分配问题。
- 基于 `iso_c_binding` 的 C 互操作 API — 可从 C、C++、Python (ctypes) 等调用。
- 输入验证：对于无效输入（例如 $N \le 0$、代价矩阵中包含 `NaN`/`Inf`），总代价返回 `NaN`，分配结果返回 `0`/`-1`。
- 零外部依赖 — 单文件 Fortran 模块。

## 环境要求

- **Fortran 编译器**：`gfortran` ≥ 7.0（或任何符合 Fortran 2008 标准的编译器）

## 构建

### 共享库

`Makefile` 会自动检测平台并构建对应的共享库：

| 平台            | 输出文件               |
|-----------------|------------------------|
| Linux           | `libhungarian.so`      |
| macOS           | `libhungarian.dylib`   |
| Windows (MinGW) | `libhungarian.dll`     |

**手动构建**（不使用 Make）：

```bash
# Linux
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.so

# macOS
gfortran -shared -fPIC -O3 -march=native hungarian.f90 -o libhungarian.dylib

# Windows (MinGW)
gfortran -shared -O3 hungarian.f90 -o libhungarian.dll
```

### 测试

```bash
make test
```

## API 参考

### Fortran API

```fortran
use hungarian_mod

call hungarian_algorithm(cost_matrix, n, assignments, total_cost)
```

| 参数          | 类型                 | Intent | 说明                                                                  |
|---------------|----------------------|--------|-----------------------------------------------------------------------|
| `cost_matrix` | `real(f64) :: (n,n)` | `in`   | $N \times N$ 代价矩阵                                                 |
| `n`           | `integer`            | `in`   | 矩阵维度                                                              |
| `assignments` | `integer :: (n)`     | `out`  | 分配结果：`assignments(i) = j` 表示行 $i \to$ 列 $j$                  |
| `total_cost`  | `real(f64)`          | `out`  | 分配代价总和（输入无效时为 `NaN`）                                    |

**示例：**

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
    const double* c_matrix_ptr,      // N×N 代价矩阵（列主序！）
    int32_t       n,                 // 矩阵维度
    int32_t* c_assignments_ptr, // 输出：N 个分配结果（0 基索引，错误时为 -1）
    double* c_cost_ptr         // 输出：总代价（错误时为 NaN）
);
```

> [!IMPORTANT]
> 矩阵必须以 **列主序**（Fortran 约定）传入。Eigen 矩阵默认即为列主序。对于行主序的 C 数组，需在调用前转置，或按 `matrix[col * n + row]` 方式存储。

**示例 (`main.cpp`)：**

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

## 算法

本实现采用 **约减 + 最小线覆盖 + 调整** 变体：

1. **行约减** — 减去每行的最小值。
2. **列约减** — 减去每列的最小值。
3. **最小线覆盖** — 找到覆盖所有零元素的最少行列线。
   
   *（这利用了 König 定理：二部图最大匹配的大小等于最小顶点覆盖的大小）。*
4. **调整** — 若覆盖线数 $< n$：从未覆盖单元格中减去最小值，并将其加到所有双重覆盖的单元格上。
5. 重复步骤 3–4，直到覆盖线数 $= n$。
6. **提取分配** — 使用深度优先搜索（DFS）寻找增广路径，从最终的零元素结构中提取最优分配。

对于良态输入，保证在 $O(n)$ 次调整迭代内收敛。本实现使用基于机器精度（`epsilon(1.0_f64) * max(1.0_f64, maxval(abs(cost_matrix)))`）的尺度感知容差，以确保浮点稳定性。