# AGENTS.md

This file provides guidance to AI coding agents when working with code in this repository.

## Project Overview

Elemental is a C++ library for distributed-memory dense and sparse-direct linear algebra, conic optimization, and lattice reduction. It uses MPI for multi-node execution. All public API lives in namespace `El`. Version 0.88-dev, BSD 2-Clause licensed.

## Build Commands

Out-of-source builds are required. Dependencies (BLAS, LAPACK, METIS) are auto-downloaded if missing.

```bash
# Basic build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~/Install
make -j$(nproc)

# Build with tests and examples
cmake .. -DEL_TESTS=ON -DEL_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run all tests
ctest --output-on-failure

# Run a specific test (e.g. LU)
ctest -R Tests/lapack_like/LU

# Install
make install
```

Key CMake options: `EL_C_INTERFACE` (ON), `EL_USE_64BIT_INTS` (OFF), `EL_HYBRID` (OFF, enables OpenMP+MPI), `BUILD_SHARED_LIBS` (ON), `INSTALL_PYTHON_PACKAGE` (OFF), `EL_USE_QT5` (OFF).

## Architecture

### Matrix Type Hierarchy

The core abstraction is a two-level matrix system:

- `Matrix<Ring>` — local sequential column-major matrix
- `AbstractDistMatrix<Ring>` — base for all distributed matrices
  - `ElementalMatrix<Ring>` — elementwise (cyclic) distributions
    - `DistMatrix<Ring, ColDist, RowDist>` — 15 distribution combinations using `MC, MR, MD, STAR, VC, VR, CIRC`
  - `BlockMatrix<Ring>` — block-cyclic distributions (ScaLAPACK compatibility)
    - `DistMatrix<Ring, ColDist, RowDist, BLOCK>`

The `Grid` class wraps an MPI communicator into a 2D process grid with sub-communicators (`MCComm()`, `MRComm()`, etc.).

### Source Layout

- `include/El.hpp` — master include; `El-lite.hpp` for core+BLAS only; `El.h` for C interface
- `include/El/core/` — `Matrix`, `DistMatrix`, `Grid`, `SparseMatrix`, `DistSparseMatrix`, MPI/BLAS/LAPACK wrappers
- `include/El/blas_like/` — level 1/2/3 BLAS-style operations
- `include/El/lapack_like/` — factorizations, eigensolvers, SVD, linear solvers
- `include/El/optimization/` — LP/QP/SOCP, ADMM, proximal methods
- `src/` — implementation `.cpp` files mirroring the include structure
- `external/` — bundled pmrrr (parallel eigensolver) and suite_sparse (AMD/LDL)

### Template Patterns

Algorithms are templated on `Field` or `Ring` with trait-based SFINAE:
- `Base<Field>` extracts the real type (e.g., `Base<Complex<double>>` = `double`)
- `IsScalar<T>`, `IsField<T>`, `IsReal<T>` — type predicates
- `EnableIf`/`DisableIf` — constraint wrappers
- Same function name works for both `Matrix<T>` and `AbstractDistMatrix<T>` overloads

Supported scalar types: `float`, `double`, `Complex<float>`, `Complex<double>`, plus optional `Quad` (__float128), `DoubleDouble`, `QuadDouble` (QD library), `BigFloat`/`BigInt` (MPFR).

### C Interface

Parallel C API in `*-C.cpp` files using the `CReflect` pattern to convert between C structs and C++ objects. Used by the Python ctypes interface in `python/`.

### Common Code Patterns

- **RAII init**: `El::Environment env(argc, argv)` initializes MPI+Elemental, finalizes on destruction
- **Input parsing**: `El::Input("--flag", "description", default_value)` then `El::ProcessInput()`
- **Algorithm control**: dedicated control structs (e.g., `BPCtrl<Real>`) parameterize algorithms
- **FLAME partitioning**: `include/El/core/FlamePart.hpp` for blocked algorithm partitioning
- **Debug mode**: `EL_DEBUG_ONLY(...)` guards, `LogicError()`/`RuntimeError()` for exceptions

### Tests and Examples

Tests in `tests/{core,blas_like,lapack_like,optimization}/` — each `.cpp` is a standalone MPI program registered with CTest. Examples follow the same pattern in `examples/`. The `sandbox/test.cpp` is always built as a minimal sanity check.
