# Elemental

A modern C++ library for distributed-memory dense linear algebra.
The library was initially released in
[Elemental: A new framework for distributed memory dense linear algebra](https://dl.acm.org/citation.cfm?doid=2427023.2427030) 
and is the key building block for the distributed-memory sparse-direct solver 
[Clique](http://www.github.com/poulson/Clique.git).

Please visit [the download page](http://libelemental.org/download/) for 
download instructions.

### Documentation

The [documentation for the development version of Elemental](http://elemental.github.com/Elemental) is built using [Sphinx](http://sphinx.pocoo.org).

### Related open-source packages

Implementations:

1. [DPLASMA](http://icl.eecs.utk.edu/dplasma/)
2. [PLAPACK](http://www.cs.utexas.edu/~plapack)
3. [ScaLAPACK](http://www.netlib.org/scalapack) (and the add-on, [ELPA](http://elpa.rzg.mpg.de/))

Wrappers:

1. [PETSc](https://www.mcs.anl.gov/petsc/)
2. [Trilinos](http://trilinos.sandia.gov)

Note that [PETSc](https://www.mcs.anl.gov/petsc/) contains interfaces for both 
[Elemental](http://github.com/elemental/Elemental.git) and
[Clique](http://github.com/poulson/Clique.git).
