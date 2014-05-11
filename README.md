<p align="left" style="padding: 20px">
<img src="http://libelemental.org/_static/elemental.png">
</p>

**Elemental** is a modern C++ library for distributed-memory dense linear 
algebra.
The library was initially released in
[Elemental: A new framework for distributed memory dense linear algebra](https://dl.acm.org/citation.cfm?doid=2427023.2427030) 
and is the key building block for the distributed-memory sparse-direct solver 
[Clique](http://www.github.com/poulson/Clique.git).

Please visit [the download page](http://libelemental.org/download/) for 
download instructions.

### Documentation

The [documentation for the development version of Elemental](http://libelemental.org/documentation) is built using [Sphinx](http://sphinx.pocoo.org).

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

### Elemental's root directory

This is the root directory of the entire project, and it contains:

-  `AUTHORS`: the list of source code contributors
-  `cmake/`: auxiliary files for CMake configuration
-  `CMakeLists.txt`: the CMake configuration file
-  `doc/`: Sphinx documentation 
-  `examples/`: various concise examples of Elemental's functionality
-  `experimental/`: experimental code which is not yet library quality
-  `external/`: non-standard external code which Elemental builds on top of
-  `include/`: Elemental's header files; most of the library resides here
-  `LICENSE`: the New BSD License file
-  `octave/`: pedagogical versions of algorithms used in Elemental
-  `PUBLICATIONS`: publications directly related to this source code
-  `README.md`: this file
-  `REFERENCES`: some publications referenced in the creation of this software
-  `src/`: Elemental's source files; a small portion of the library is here
-  `tests/`: programs meant to test the accuracy of Elemental
-  `TODO`: a list of near/long-term goals of the project
-  `vagrant/`: lightweight preconfigured virtual-machines for Elemental
