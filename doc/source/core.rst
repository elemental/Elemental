Core functionality
******************

The Matrix class
================
This is the basic building block of the library: its purpose it to provide 
convenient matrix abstractions over a contiguous chunk of memory, e.g., the
ability to set and/or get individual entries according to their matrix indices.
An important property of the ``Matrix`` class is that it is templated over 
the underlying datatype; it can just as easily make use of ``float``, 
``double``, and ``int``, as well as their complex counterparts. There are 
essentially three constraints on the underlying datatype, say ``T``:

1. The datatype must be self-contained in the sense that all of 
   the relevant data of an instance of the datatype is contained within a 
   contiguous buffer of size less than or equal to ``sizeof(T)``.
2. It must be sensical to cast the decimal literals ``0`` and ``1`` into 
   type ``T``, e.g., ``static_cast<T>(0)`` and ``static_cast<T>(1)``.
3. Setting all ``sizeof(T)`` bytes of an instance of type ``T`` to zero must
   be logically equivalent to setting the instance to ``static_cast<T>(0)``.

It is important to note that these requirements are not strictly weaker or 
strictly stronger than requiring ``T`` to be POD (Plain Old Data), as 
requirements (2) and (3) are not necessary for a type to be POD, and POD 
types are not allowed to have user-defined constructors. A good example of these
differences is any complex datatype, e.g., ``std::complex<double>``, which fits
all three of our requirements but is not POD since it has several custom 
constructors.

An example of generating an :math:`m \times n` matrix of real double-precision 
numbers where the :math:`(i,j)` entry is equal to :math:`i-j` would be:

  .. code-block:: cpp

     #include "elemental.hpp"
     ...
     elemental::Matrix<double> A( m, n );
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             A.Set( i, j, (double)i-j );
     
In order to retain compatibility with BLAS and LAPACK interfaces, the
underlying data storage is simply a contiguous buffer that stores entries 
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

It is important to note that the ``Matrix`` class supports much more general 
datatypes than just ``double``...**left off here**

The DistMatrix class
====================

The ``Axpy`` interface
======================

Environment routines
====================
