Core functionality
******************

The Matrix class
================
This is the basic building block of the library: its purpose it to provide 
convenient mechanisms for performing basic matrix manipulation operations, 
such as setting and querying individual matrix entries, without giving up 
compatibility with interfaces such as BLAS and LAPACK, which assume column-major
storage.

.. note:: 

   There are essentially three constraints on the underlying datatype, say 
   ``T``:

    1. The datatype must be self-contained in the sense that all of 
       the relevant data of an instance of the datatype is contained within a 
       contiguous buffer of size less than or equal to ``sizeof(T)``.
    2. It must be sensical to cast the decimal literals ``0`` and ``1`` into 
       type ``T``, e.g., ``static_cast<T>(0)`` and ``static_cast<T>(1)``. 
       This requirement exists so that a routine can be provided to initialize
       a matrix to the identity operator.
    3. Setting all ``sizeof(T)`` bytes of an instance of type ``T`` to zero must
       be logically equivalent to setting the instance to ``static_cast<T>(0)``.
       This constraint is enforced so that ``std::memset`` can be used to 
       quickly set matrices to zero.

    It is important to note that these requirements are not strictly weaker or 
    strictly stronger than requiring ``T`` to be POD (Plain Old Data), as 
    requirements (2) and (3) are not necessary for a type to be POD, whereas 
    POD types are not allowed to have user-defined constructors. Any complex
    variable (defined using ``std::complex``) fits our requirements, but is not
    considered POD since it has several custom constructors.

An example of generating an :math:`m \times n` matrix of real double-precision 
numbers where the :math:`(i,j)` entry is equal to :math:`i-j` would be:

  .. code-block:: cpp

     #include "elemental.hpp"
     ...
     elemental::Matrix<double> A( m, n );
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             A.Set( i, j, (double)i-j );
     
The underlying data storage is simply a contiguous buffer that stores entries 
in a column-major fashion with an arbitrary leading dimension. For modifiable
instances of the ``Matrix`` class, the routine
``T* elemental::Matrix<T>::Buffer()`` returns a pointer to the underlying 
buffer, while ``int elemental::Matrix<T>::LDim() const`` returns the leading 
dimension; these two routines could be used to directly perform the equivalent
of the first code sample as follows:

  .. code-block:: cpp
     
     #include "elemental.hpp"
     ...
     elemental::Matrix<double> A( m, n );
     double* buffer = A.Buffer();
     const int ldim = A.LDim();
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             buffer[i+j*ldim] = (double)i-j;

For constant instances of the ``Matrix`` class, a ``const`` pointer
to the underlying data can similarly be returned with a call to 
``A.LockedBuffer()``. In addition, a (``const``) pointer to the place in the 
(``const``) buffer where entry :math:`(i,j)` resides can be easily retrieved
with a call to ``A.Buffer(i,j)`` or ``A.LockedBuffer(i,j)``.

It is also important to be able to create matrices which are simply *views* 
of existing (sub)matrices. For example, if ``A`` is a :math:`10 \times 10` 
matrix of complex doubles, then a matrix :math:`A_{BR}` can easily be created 
to view the bottom-right :math:`6 \times 7` submatrix using

  .. code-block:: cpp

     #include "elemental.hpp"
     ...
     elemental::Matrix<std::complex<double> > ABR;
     ABR.View( A, 4, 3, 6, 7 );

since the bottom-right :math:`6 \times 7` submatrix beings at index 
:math:`(4,3)`. In general, to view the :math:`M \times N` submatrix starting
at entry :math:`(i,j)`, one would call ``ABR.View( A, i, j, M, N );``.

The DistMatrix class
====================
The ``DistMatrix`` class is meant to provide a distributed-memory analogue of 
the ``Matrix`` class. Rather than simply supporting a single type of 
two-dimensional matrix distributions, it templates over the manner in which
individual columns and rows of the matrix should be distributed 
(in addition to templating over the underlying datatype).
Since it is crucial to know not only how many 
processes to distribute the data over, but *which* processes, and in what 
manner they should be decomposed into a logical two-dimensional grid, an 
instance of the ``Grid`` class must be passed into the constructor of 
the ``DistMatrix`` class.

For example, ...

The ``Axpy`` interface
======================

Environment routines
====================
