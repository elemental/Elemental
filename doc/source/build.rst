==================
Building Elemental
==================
Elemental's build system heavily relies on `CMake <http://www.cmake.org>`_ 
in order to manage a large number of configuration options in a 
platform-independent manner.

----------------
Installing CMake
----------------
The `installation process <http://www.cmake.org/cmake/help/install.html>`_
is extremely straightforward: either download a platform-specific binary from
the `downloads page <http://www.cmake.org/cmake/resources/software.html>`_
or download the source code and run::

    ./bootstrap
    make

Though the binaries can be run from within the build directory, it might be 
more convenient to install them by subsequently running::

    make install

------------------
Working with CMake
------------------
Describe the process of creating a separate build tree that must be cleaned
before each configuration step.

---------------------
Configuring Elemental
---------------------
Discuss all of the available options, like optionally building PMRRR, as well
as the different build modes. Also give examples of building PureRelease versus
HybridDebug.

--------------------
Installing Elemental
--------------------
Discuss the differences between running "make" and "make install", and how
to choose the installation directory.

---------------
Using Elemental
---------------
Discuss how to link against libelemental and which headers need to be included, etc.
