<p align="left" style="padding: 20px">
<img src="http://libelemental.org/_static/elemental.png">
</p>

**Elemental** is a modern C++ library for distributed-memory dense and 
sparse-direct linear algebra and optimization.
The library was initially released in
[Elemental: A new framework for distributed memory dense linear algebra](https://dl.acm.org/citation.cfm?doid=2427023.2427030) 
and has recently absorbed the functionality from the sparse-direct solver
[Clique](http://www.github.com/poulson/Clique.git).

Please visit [the download page](http://libelemental.org/download/) for 
download instructions, but note that 0.85 wil be the first release containing sparse-direct support.

### Documentation

The [documentation for the development version of Elemental](http://libelemental.org/documentation) is built using [Sphinx](http://sphinx.pocoo.org).

### Related open-source packages

Implementations (in alphabetical order):

* [ELPA](http://elpa.rzg.mpg.de)
* [NuLAB](https://github.com/solomonik/NuLAB)
* [PaRSEC/DPLASMA](http://icl.eecs.utk.edu/projectsdev/parsec/index.html)
* [PLAPACK](http://www.cs.utexas.edu/~plapack)
* [ScaLAPACK](http://www.netlib.org/scalapack) 

Wrappers:

* [PETSc](https://www.mcs.anl.gov/petsc/)
* [Trilinos](http://trilinos.sandia.gov)

Note that [PETSc](https://www.mcs.anl.gov/petsc/) contains an interface for
some of Elemental's functionality.

### Elemental's root directory

This is the root directory of the entire project, and it contains:

-  `AUTHORS`: the list of source code contributors
-  `cmake/`: auxiliary files for CMake configuration
-  `CMakeLists.txt`: the CMake configuration file
-  `examples/`: various concise examples of Elemental's functionality
-  `experimental/`: experimental code which is not yet library quality
-  `external/`: external code which Elemental (optionally) builds on top of
-  `include/`: Elemental's C++ and C header files, as well as a python interface
-  `LICENSE`: the New BSD License file
-  `octave/`: pedagogical versions of algorithms used in Elemental (at the moment, just a rank-one modification of an LU factorization)
-  `PUBLICATIONS`: publications directly related to this source code
-  `README.md`: this file
-  `REFERENCES`: some publications referenced in the creation of this software
-  `sanbox/`: A makefile and test driver meant to be modified for testing an installation of the library
-  `src/`: Elemental's source files
-  `tests/`: programs meant to test the accuracy/correctness of Elemental
-  `TODO`: an incomplete list of near/long-term goals of the project
-  `vagrant/`: lightweight preconfigured virtual-machines for Elemental
