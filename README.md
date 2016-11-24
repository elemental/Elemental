<p align="left" style="padding: 20px">
<img src="http://libelemental.org/_static/elemental.png">
</p>

[![Build Status](https://api.travis-ci.org/elemental/Elemental.svg?branch=master)](https://travis-ci.org/elemental/Elemental)
[![Join the chat at https://gitter.im/elemental/chat](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/elemental/chat?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

**Elemental** is a modern C++ library for distributed-memory dense and
sparse-direct linear algebra, conic optimization, and lattice reduction.
The library was initially released in
[Elemental: A new framework for distributed memory dense linear algebra](https://dl.acm.org/citation.cfm?doid=2427023.2427030)
and absorbed, then greatly expanded upon, the functionality from the 
sparse-direct solver [Clique](http://www.github.com/poulson/Clique.git), which 
was originally released during a project on [Parallel Sweeping Preconditioners](http://epubs.siam.org/doi/abs/10.1137/120871985).

### Documentation

The (now outdated) [documentation for Elemental](http://libelemental.org/documentation) is built using [Sphinx](http://sphinx.pocoo.org) and the [Read the Docs Theme](http://docs.readthedocs.org/en/latest/theme.html)

### Unique features

Elemental supports a wide collection of sequential and distributed-memory
functionality, including sequential and distributed-memory support for the
datatypes:

- `float`, `El::Complex<float>`
- `double`, `El::Complex<double>`
- `El::DoubleDouble`, `El::Complex<El::DoubleDouble>` (on top of QD's *dd_real*)
- `El::QuadDouble`, `El::Complex<El::QuadDouble>` (on top of QD's *qd_real*)
- `El::Quad`, `El::Complex<El::Quad>` (on top of GCC's *__float128*)
- `El::BigFloat`, `El::Complex<El::BigFloat>` (on top of MPFR's *mpfr_t* and MPC's *mpc_t*)

**Linear algebra**:
* Dense and sparse-direct (generalized) Least Squares
  problems
    - Least Squares / Minimum Length
    - Tikhonov (and ridge) regression
    - Equality-constrained Least Squares
    - General (Gauss-Markov) Linear Models
* High-performance pseudospectral computation and visualization
* Aggressive Early Deflation Schur decompositions (currently sequential only)
* Blocked column-pivoted QR via Johnson-Lindenstrauss
* Quadratic-time low-rank Cholesky and LU modifications
* Bunch-Kaufman and Bunch-Parlett for accurate symmetric
  factorization
* LU and Cholesky with full pivoting
* Column-pivoted QR and interpolative/skeleton decompositions
* Quadratically Weighted Dynamic Halley iteration for the polar decomposition
* Many algorithms for Singular-Value soft-Thresholding (SVT)
* Tall-skinny QR decompositions
* Hermitian matrix functions
* Prototype Spectral Divide and Conquer Schur decomposition and Hermitian EVD
* Sign-based Lyapunov/Ricatti/Sylvester solvers
* Arbitrary-precision distributed SVD (QR and D&C support), (generalized) Hermitian EVPs (QR and D&C support), and Schur decompositions (e.g., via Aggressive Early Deflation)

**Convex optimization**:
* Dense and sparse Interior Point Methods for
  Linear, Quadratic, and Second-Order Cone Programs (**Note: Scalability for sparse IPMs will be lacking until more general sparse matrix distributions are introduced into Elemental**)
    - Basis Pursuit
    - Chebyshev Points
    - Dantzig selectors
    - LASSO / Basis Pursuit Denoising
    - Least Absolute Value regression
    - Non-negative Least Squares
    - Support Vector Machines
    - (1D) Total Variation
* Jordan algebras over products of Second-Order Cones
* Various prototype dense Alternating Direction Method of Multipliers routines
    - Sparse inverse covariance selection
    - Robust Principal Component Analysis
* Prototype alternating direction Non-negative Matrix Factorization

**Lattice reduction**:
* An extension of [Householder-based LLL](http://perso.ens-lyon.fr/damien.stehle/HLLL.html) to real and complex linearly-dependent bases (currently sequential only)
* Generalizations of [BKZ 2.0](http://link.springer.com/chapter/10.1007%2F978-3-642-25385-0_1) to complex bases (currently sequential only)
 incorporating ["y-sparse" enumeration](https://eprint.iacr.org/2014/980)
* Integer images/kernels and relation-finding (currently sequential only)

### The current development roadmap

**Core data structures**:
* (1a) Eliminate `DistMultiVec` in favor of the newly extended `DistMatrix`
* (1b) Extend `DistSparseMatrix` to support elementwise and blockwise 2D distributions

**Linear algebra**:
* (2a) Distributed iterative refinement tailored to two right-hand sides \[weakly depends on (1a)\]
* (2b) Extend black-box iterative refinement to `DistMatrix`
* (2c) Incorporate iterative refinement into linear solvers via optional control
  structure \[weakly depends upon (2b)\]
* (2d) Support for the Fix-Heiberger method for accurate generalized Hermitian-definite EVPs

**Convex optimization**:
* (3a) Add support for homogeneous self-dual embeddings \[weakly depends on (2a)\]
* (3b) Enhance sparse scalability via low edge-degree plus low-rank 
  decompositions \[depends on (1b); weakly depends on (1a)\]
* (3c) Distributed sparse semidefinite programs via chordal decompositions \[weakly depends on (3b)\]

### License

The vast majority of Elemental is distributed under the terms of the
[New BSD License](http://www.opensource.org/licenses/bsd-license.php).
Please see the [debian/copyright](https://github.com/elemental/Elemental/blob/master/debian/copyright) file for an overview of the copyrights and licenses for
the files in the library.

The optional external dependency
[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
is distributed under the (equally permissive)
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html),
though
[ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
can only be used for research purposes (and can be easily disabled).
[libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) is 
distributed under the terms of the [GNU Lesser General Public License, version 2.1 or later](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html),
while,
[QD](http://crd-legacy.lbl.gov/~dhbailey/mpdist/) is distributed under the
terms of the [LBNL-BSD-License](http://crd.lbl.gov/~dhbailey/mpdist/LBNL-BSD-License.doc).

### Dependencies

**Intranodal linear algebra**

* [BLAS](http://netlib.org/blas)
* [LAPACK](http://netlib.org/lapack)
* [libflame](http://www.cs.utexas.edu/~flame/web/libFLAME.html) (optional for faster bidiagonal SVDs)
* Elemental is packed with a greatly modified version of the Alternating
  Minimum Degree (AMD) reordering and unblocked sparse LDL factorization from
  [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) (**Note:** The used portions of SuiteSparse are licensed under the [GNU Lesser General Public License, version 2.1 or later](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html))

[OpenBLAS](http://www.openblas.net) is automatically downloaded and installed
if no vendor/tuned BLAS/LAPACK is detected.

**Intranodal graph partitioning**

* [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
* [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) (**Note:** commercial users must disable this option during configuration)

If ParMETIS is not disabled and cannot be found (including access to internal APIs), then it is automatically downloaded and installed;
otherwise, if METIS support is not detected, METIS is downloaded and installed.

**Internodal linear algebra**

* [Parallel MRRR](https://code.google.com/p/pmrrr/) (packaged with Elemental)
* [ScaLAPACK](http://netlib.org/scalapack) (optional for benchmarking)

If [ScaLAPACK](http://www.netlib.org/scalapack) support is not explicitly 
disabled, then Elemental looks for a previous installation and, failing that,
attempts to automatically download and install the library.

**Internodal communication**

* MPI2 (typically [MPICH](http://www.mpich.org/), [MVAPICH](http://mvapich.cse.ohio-state.edu/), or [OpenMPI](http://www.open-mpi.org/))

**Auxiliary libraries**

* [QD](http://crd-legacy.lbl.gov/~dhbailey/mpdist/) for efficient software analogues of 128-bit and 256-bit floating-point arithmetic (**Note:** QD is licensed under the [LBNL-BSD-License](http://crd.lbl.gov/~dhbailey/mpdist/LBNL-BSD-License.doc), which is a slight modification of the BSD License)

* [libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) for quad-precision support (especially for iterative refinement). (**Note:** libquadmath is licensed under the [GNU Lesser General Public License, version 2.1 or later](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html))

* [MPFR](http://www.mpfr.org/) for arbitrary-precision real arithmetic. (**Note:** MPFR is licensed under the [GNU Lesser General Public License, v3 or later](http://www.gnu.org/copyleft/lesser.html))

* [MPC](http://www.multiprecision.org/index.php?prog=mpc) for arbitrary-precision complex arithmetic. (**Note:** MPC is licensed under the [GNU Lesser General Public License, v3 or later](http://www.gnu.org/copyleft/lesser.html))

**Python interface**

* [matplotlib](http://matplotlib.org/) (optional for Python matrix visualization)
* [NetworkX](https://networkx.github.io/) (optional for Python graph visualization)
* [NumPy](http://www.numpy.org/)

**C++ visualization**

* [Qt5](http://qt-project.org/qt5) (optional for visualization from C++)

**Build system**

* [CMake >= 2.8.12](http://www.cmake.org/)

### Third-party interfaces

In addition to the C++11, C, and Python interfaces included within the project,
three external interfaces are currently being externally developed:

* [R-El](https://github.com/rocanale/R-Elemental) is an [R](http://www.r-project.org) interface to Elemental developed by [Rodrigo Canales](https://github.com/rocanale) and [Paolo Bientinesi](http://hpac.rwth-aachen.de/~pauldj/)

* [Elemental.jl](https://github.com/JuliaParallel/Elemental.jl) is an (in-progress) [Julia](http://julialang.org) interface to Elemental being developed by [Jake Bolewski](https://github.com/jakebolewski), [Jiahao Chen](https://jiahao.github.io), and [Andreas Noack](http://andreasnoack.github.io/academiccv.html).

* [CVXPY](https://github.com/cvxgrp/cvxpy) is a Python-embedded modeling language for convex optimization problems with an (in-progress) interface to Elemental's distributed Interior Point Methods. This effort is being led by [Steven Diamond](http://web.stanford.edu/~stevend2/).

### Related open-source projects

**Distributed dense linear algebra**:

* [ELPA](http://elpa.rzg.mpg.de)
* [CANDMC](https://github.com/solomonik/CANDMC)
* [PaRSEC/DPLASMA](http://icl.eecs.utk.edu/projectsdev/parsec/index.html)
* [PLAPACK](http://www.cs.utexas.edu/~plapack)
* [ScaLAPACK](http://www.netlib.org/scalapack)

**Distributed sparse-direct linear algebra**:

* [DSCPACK](http://www.cse.psu.edu/~raghavan/Dscpack/)
* [MUMPS](http://mumps.enseeiht.fr/)
* [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)

**Distributed linear algebra Frameworks**

* [PETSc](https://www.mcs.anl.gov/petsc/)
* [Trilinos](http://trilinos.sandia.gov)

**Convex optimization**

* [CVXOPT](http://cvxopt.org/)
* [ECOS](https://github.com/embotech/ecos)
* [L1-MAGIC](http://users.ece.gatech.edu/~justin/l1magic/)
* [SDPA](http://sdpa.sourceforge.net/index.html)

**Lattice reduction and number theory**

* [FPLLL](https://github.com/dstehle/fplll)
* [NTL](http://www.shoup.net/ntl/)
