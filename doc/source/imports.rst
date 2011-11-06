Imported library routines
*************************
Since one of the goals of Elemental is to provide high-performance 
datatype-independent parallel routines, yet Elemental's dependencies are 
datatype-dependent, it is convenient to first build a thin datatype-independent
abstraction on top of the necessary routines from BLAS, LAPACK, and MPI. The 
"first-class" datatypes are ``float``, ``double``, ``std::complex<float>``, and 
``std::complex<double>``, but ``int`` and ``byte`` (``unsigned char``) are 
supported for many cases, and support for higher precision arithmetic is in the
works.

.. toctree::
   :maxdepth: 2

   imports/blas
   imports/lapack
   imports/mpi
   imports/plcg
   imports/pmrrr
