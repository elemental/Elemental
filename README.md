<p align="left" style="padding: 20px">
<img src="http://libelemental.org/_static/elemental.png">
</p>

![Travis CI](https://travis-ci.org/elemental/Elemental.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/elemental/Elemental/badge.svg?branch=master)](https://coveralls.io/r/elemental/Elemental?branch=master)

**Elemental** is a modern C++ library for distributed-memory dense and
sparse-direct linear algebra and optimization.
The library was initially released in
[Elemental: A new framework for distributed memory dense linear algebra](https://dl.acm.org/citation.cfm?doid=2427023.2427030)
and has recently absorbed the functionality from the sparse-direct solver
[Clique](http://www.github.com/poulson/Clique.git).

Please visit [the download page](http://libelemental.org/download/) for
details about recent and upcoming releases.

### Documentation

The [documentation for Elemental](http://libelemental.org/documentation) is built using [Sphinx](http://sphinx.pocoo.org) and the [Read the Docs Theme](http://docs.readthedocs.org/en/latest/theme.html)

### Unique features

Elemental supports a wide collection of distributed-memory functionality,
including:

**Convex optimization**:
* Dense and sparse Interior Point Methods for Linear and Quadratic Programs
    - Basis Pursuit
    - Chebyshev Points
    - Dantzig selectors
    - LASSO / Basis Pursuit Denoising
    - Least Absolute Value regression
    - Non-negative Least Squares
    - Support Vector Machines
    - (1D) Total Variation
* Various prototype dense Alternating Direction Method of Multipliers routines
    - Sparse inverse covariance selection
    - Robust Principal Component Analysis
* Prototype alternating direction Non-negative Matrix Factorization

**Linear algebra**:
* Dense and sparse-direct (generalized) Least Squares problems
    - Least Squares / Minimum Length
    - Tikhonov (and ridge) regression
    - Equality-constrained Least Squares
    - General (Gauss-Markov) Linear Models
* High-performance pseudospectral computation and visualization
* Quadratic-time low-rank Cholesky and LU modifications
* Bunch-Kaufman and Bunch-Parlett for accurate symmetric factorization
* LU and Cholesky with full pivoting
* Column-pivoted QR and interpolative/skeleton decompositions
* Quadratically Weighted Dynamic Halley iteration for the polar decomposition
* Many algorithms for Singular-Value soft-Thresholding (SVT)
* Tall-skinny QR decompositions
* Hermitian matrix functions
* Prototype Spectral Divide and Conquer Schur decomposition and Hermitian EVD
* Sign-based Lyapunov/Ricatti/Sylvester solvers

### License

All of Elemental is distributed under the terms of the
[New BSD License](http://www.opensource.org/licenses/bsd-license.php),
with the exception of
[METIS 5.1.0](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview),
which is distributed under the (equally permissive)
[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html).

### Related open-source projects

**Distributed dense linear algebra**:

* [ELPA](http://elpa.rzg.mpg.de)
* [NuLAB](https://github.com/solomonik/NuLAB)
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

### Dependencies

**Intranodal linear algebra**

* [BLAS](http://netlib.org/blas)
* [LAPACK](http://netlib.org/lapack)
* [libflame](http://www.cs.utexas.edu/~flame/web/libFLAME.html) (optional for faster bidiagonal SVDs)

**Intranodal graph partitioning**

* [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
* [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) (**Note:** support is currently **not** enabled due to the restrictive license)

**Internodal linear algebra**

* [Parallel MRRR](https://code.google.com/p/pmrrr/) (packaged with Elemental)
* [ScaLAPACK](http://netlib.org/scalapack) (optional for Hessenberg QR algorithm)

**Internodal communication**

* MPI2 (typically [MPICH](http://www.mpich.org/), [MVAPICH](http://mvapich.cse.ohio-state.edu/), or [OpenMPI](http://www.open-mpi.org/))

**Python interface**

* [matplotlib](http://matplotlib.org/) (optional for python matrix visualization)
* [NetworkX](https://networkx.github.io/) (optional for python graph visualization)
* [NumPy](http://www.numpy.org/)

**C++ visualization**

* [Qt5](http://qt-project.org/qt5) (optional for visualization from C++)

**Build system**

* [CMake >= 2.8.12](http://www.cmake.org/)
