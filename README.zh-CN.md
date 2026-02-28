# 匈牙利算法 — Modern Fortran 实现

[![CI](https://github.com/chuantian/hungarian-algorithm-fortran/actions/workflows/ci.yml/badge.svg)](https://github.com/chuantian/hungarian-algorithm-fortran/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[English](README.md)

基于 **匈牙利算法**（Kuhn-Munkres）的高性能线性和分配问题求解器。使用 Modern Fortran 编写，通过 C 兼容接口实现与 C/C++ 的无缝集成。

## 特性

- O(n³) 时间复杂度求解 N×N 代价矩阵分配问题
- 基于 `iso_c_binding` 的 C 互操作 API — 可直接从 C、C++、Python (ctypes) 等调用
- 完善的错误处理，提供类型化错误码
- 输入验证：拒绝 NaN、Inf 及无效维度
- 尺度感知的浮点容差，确保数值稳定性
- 零外部依赖 — 单文件 Fortran 模块

## 环境要求

- **Fortran 编译器**：`gfortran` ≥ 7.0（或任何符合 Fortran 2008 标准的编译器）
- **Make**（可选，用于自动化构建）
- **Eigen**（可选，仅 C++ 示例需要）

## 快速开始

```bash
# 构建共享库
make

# 构建并运行测试
make test

# 清理构建产物
make clean
```

## 构建

### 共享库

`Makefile` 会自动检测平台并构建对应的共享库：

| 平台 | 输出文件 |
| --- | --- |
| Linux | `libhungarian.so` |
| macOS | `libhungarian.dylib` |
| Windows (MinGW) | `libhungarian.dll` |

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

此命令会编译 `hungarian.f90` 和 `test_hungarian.f90`，然后运行测试套件。预期输出：

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

## API 参考

### Fortran API

```fortran
use hungarian_mod

call hungarian_algorithm(cost_matrix, assignments, total_cost, info)
```

| 参数 | 类型 | Intent | 说明 |
| --- | --- | --- | --- |
| `cost_matrix` | `real(f64) :: (:,:)` | `in` | N×N 代价矩阵 |
| `assignments` | `integer :: (:)` | `out` | 分配结果：`assignments(i) = j` 表示行 i 分配到列 j |
| `total_cost` | `real(f64)` | `out` | 分配代价总和 |
| `info` | `integer` | `out` | 错误码（见下表） |

**示例：**

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
    const double* c_matrix_ptr,     // N×N 代价矩阵（列主序！）
    int32_t       n,                // 矩阵维度
    int32_t*      c_assignments_ptr, // 输出：N 个分配结果（0 基索引）
    double*       c_cost_ptr,       // 输出：总代价
    int32_t*      c_info_ptr        // 输出：错误码
);
```

> **⚠️ 必须使用列主序：** 矩阵必须以列主序（Fortran 约定）传入。Eigen 矩阵默认即为列主序。对于行主序的 C 数组，需在调用前转置，或按 `matrix[col * n + row]` 方式存储。

**示例 (`main.cpp`)：**

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

**输出：**

```text
Total Cost: 140
Worker 0 -> Job 2
Worker 1 -> Job 1
Worker 2 -> Job 0
Worker 3 -> Job 3
```

### 错误码

| 常量 | 值 | 说明 |
| --- | --- | --- |
| `HUNGARIAN_OK` | 0 | 成功 |
| `HUNGARIAN_ERR_INVALID` | 1 | 无效输入（n < 0、矩阵含 NaN/Inf、空指针） |
| `HUNGARIAN_ERR_ALLOC` | 2 | 内存分配失败 |
| `HUNGARIAN_ERR_NO_CONVERGE` | 3 | 算法在迭代上限内未收敛 |
| `HUNGARIAN_ERR_NO_MATCH` | 4 | 未找到完美匹配 |

## 算法

本实现采用 **约减 + 最小线覆盖 + 调整** 变体：

1. **行约减** — 减去每行的最小值
2. **列约减** — 减去每列的最小值
3. **最小线覆盖** — 通过 König 定理（二部图最大匹配 = 最小顶点覆盖）找到覆盖所有零元素的最少行列线
4. **调整** — 若覆盖线数 < n：从未覆盖单元格中减去最小值，加到双重覆盖单元格上
5. 重复步骤 3–4，直到覆盖线数 = n
6. 从最终零元素结构中 **提取最优分配**

对良态输入，收敛保证在 O(n) 次调整迭代内完成。实现使用尺度感知容差（`1e-12 × max|cost|`）确保浮点稳定性。

## 项目结构

```text
hungarian-algorithm-fortran/
├── hungarian.f90          # 核心模块：算法 + C 包装器
├── test_hungarian.f90     # 测试套件（15 个测试）
├── Makefile               # 跨平台构建系统
├── .github/
│   └── workflows/
│       └── ci.yml         # GitHub Actions CI
├── LICENSE                # MIT 许可证
└── README.md
```

## 许可证

[MIT](LICENSE)
