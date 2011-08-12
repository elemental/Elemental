Elemental's build system relies on `CMake <http://www.cmake.org>`_ 
in order to manage a large number of configuration options in a 
platform-independent manner; it can be easily configured to build on Linux and 
Unix environments (including Darwin) as well as various versions of 
Microsoft Windows.

Installing CMake
================
Elemental uses several new CMake modules, so it is important to ensure that 
version 2.8.5 or later is installed. Thankfully the 
`installation process <http://www.cmake.org/cmake/help/install.html>`_
is extremely straightforward: either download a platform-specific binary from
the `downloads page <http://www.cmake.org/cmake/resources/software.html>`_,
or instead grab the most recent stable tarball and have CMake bootstrap itself.
In the simplest case, the bootstrap process is as simple as running the 
following commands::

    ./bootstrap
    make
    make install

There are two important issues to consider:

1. By default, ``make install`` attempts a system-wide installation 
   (e.g., into ``/usr/bin``) and will likely require administrative privileges.
   A different installation folder may be specified with the ``--prefix`` 
   option to the ``bootstrap`` script, e.g.,::

    ./bootstrap --prefix=/home/your_username
    make
    make install

   Afterwards, it is a good idea to make sure that the environment variable 
   ``PATH`` includes the ``bin`` subdirectory of the installation folder, e.g.,
   ``/home/your_username/bin``.

2. Some highly optimizing compilers will not correctly build CMake, but the GNU
   compilers nearly always work. You can specify which compilers to use by
   setting the environment variables ``CC`` and ``CXX`` to the full paths to 
   your preferred C and C++ compilers before running the ``bootstrap`` script.

Working with CMake
==================
Though many configuration utilities, like 
`autoconf <http://www.gnu.org/software/autoconf/>`_, are designed such that
the user need only invoke ``./configure && make && make install`` from the
top-level source directory, CMake targets *out-of-source* builds, which is 
to say that the build process occurs away from the source code. The 
out-of-source build approach is ideal for projects that offer several 
different build modes, as each version of the project can be built in a 
separate folder.

A common approach is to create a folder named ``build`` in the top-level of 
the source directory and to invoke CMake from within it::

    mkdir build
    cd build
    cmake ..

The last line calls the command line version of CMake, ``cmake``,
and tells it that it should look in the parent directory for the configuration
instructions, which should be in a file named ``CMakeLists.txt``. Users that 
would prefer a graphical interface from the terminal (through ``curses``)
should instead use ``ccmake`` (on Unix platforms) or ``CMakeSetup`` 
(on Windows platforms). In addition, a GUI version is available through 
``cmake-gui``. 

Though running ``make clean`` will remove all files generated from running 
``make``, it will not remove configuration files. Thus, the best approach for
completely cleaning a build is to remove the entire build folder. On \*nix 
machines, this is most easily accomplished with::

    cd .. 
    rm -rf build

This is a better habit than simply running ``rm -rf *`` since, 
if accidentally run from the wrong directory, the former will most likely fail.

Getting Elemental's source 
==========================
There are two basic approaches:

1. Download a tarball of the most recent version from 
   `http://code.google.com/p/elemental/downloads/list 
   <http://code.google.com/p/elemental/downloads/list>`_. A new version is 
   released roughly once a month, on average.

2. Install `Mercurial <http://mercurial.selenic.com>`_ and check out a copy of 
   the repository by running ::

    hg clone https://elemental.googlecode.com/hg elemental

Building Elemental
==================
On \*nix machines with `BLAS <http://www.netlib.org/blas/>`_, 
`LAPACK <http://www.netlib.org/lapack/>`_, and 
`MPI <http://en.wikipedia.org/wiki/Message_Passing_Interface>`_ installed in 
standard locations, building Elemental can be as simple as::

    cd elemental
    mkdir build
    cd build
    cmake ..
    make
    make install

As with the installation of CMake, the default install location is 
system-wide, e.g., ``/usr/local``. The installation directory can be changed
at any time by running::

    cmake -D CMAKE_INSTALL_PREFIX=/your/desired/install/path ..
    make install


Though the above instructions will work on many systems, it is common to need
to manually specify several build options, especially when multiple versions of
libraries or several different compilers are available on your system. For 
instance, any C++, C, or Fortran compiler can respectively be set with the 
``CMAKE_CXX_COMPILER``, ``CMAKE_C_COMPILER``, and ``CMAKE_Fortran_COMPILER`` 
variables, e.g., ::

    cmake -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
          -D CMAKE_C_COMPILER=/usr/bin/gcc   \
          -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran ..
    
It is also common to need to specify which libraries need to be linked in order
to provide serial BLAS and LAPACK routines. The ``MATH_LIBS`` variable was
introduced for this purpose and an example usage would be ::

    cmake -D MATH_LIBS="-L/usr/lib -llapack -lblas -lm" ..

It is important to ensure that if library A depends upon library B, A should 
be specified to the left of B.

Discuss all of the available options, like optionally building PMRRR, as well
as the different build modes. Also give examples of building PureRelease versus
HybridDebug.

Using Elemental
===============
Give an example of linking against Elemental.
