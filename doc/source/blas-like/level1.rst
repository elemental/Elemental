Level 1
=======

The prototypes for the following routines can be found at 
`include/elemental/blas-like.hpp <../../../../include/elemental/blas-like.hpp>`_, while the
implementations are in `include/elemental/blas-like/level1/ <../../../../include/elemental/blas-like/level1/>`_.

Adjoint
-------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`B := A^H`. 

.. cpp:function:: void Adjoint( const Matrix<T>& A, Matrix<T>& B )

   The serial version (templated over the datatype).

.. cpp:function:: void Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed version (templated over the datatype and the individual 
   distributions of :math:`A` and :math:`B`).

Axpy
----
Performs :math:`Y := \alpha X + Y` (hence the name *axpy*).

.. cpp:function:: void Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )

   The distributed implementation (templated over the datatype and the shared
   distribution of :math:`A` and :math:`B`).

Conjugate
---------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`A := \bar A`. For real datatypes, this is a no-op.

.. cpp:function:: void Conjugate( Matrix<T>& A )

   The serial version (templated over datatype).

.. cpp:function:: void Conjugate( DistMatrix<T,U,V>& A )

   The distributed version (templated over the datatype and the distribution of
   :math:`A`).

:math:`B := \bar A`.

.. cpp:function:: void Conjugate( const Matrix<T>& A, Matrix<T>& B )

   The serial version (templated over the datatype).

.. cpp:function:: void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed version (templated over the datatype and the individual 
   distributions of :math:`A` and :math:`B`).


Copy
----
Sets :math:`Y := X`.

.. cpp:function:: void Copy( const Matrix<T>& X, Matrix<T>& Y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed implementation (templated over the datatype and the
   individual distributions of :math:`A` and :math:`B`).

DiagonalScale
-------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Performs either :math:`X := \mbox{op}(D) X` or :math:`X := X \mbox{op}(D)`, 
where :math:`op(D)` equals :math:`D=D^T`, or :math:`D^H=\bar D`, where
:math:`D = \mbox{diag}(d)` and :math:`d` is a column vector.

.. cpp:function:: void DiagonalScale( LeftOrRight side, Orientation orientation, const Matrix<T>& d, Matrix<T>& X )

   The serial implementation (templated over the datatype).

.. cpp:function:: void DiagonalScale( LeftOrRight side, Orientation orientation, const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`d` and :math:`X`).

DiagonalSolve
-------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Performs either :math:`X := \mbox{op}(D)^{-1} X` or 
:math:`X := X \mbox{op}(D)^{-1}`, where :math:`D = \mbox{diag}(d)` and :math:`d`
is a column vector.

.. cpp:function:: void DiagonalSolve( LeftOrRight side, Orientation orientation, const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular=false )

   The serial implementation (templated over the datatype).

.. cpp:function:: void DiagonalSolve( LeftOrRight side, Orientation orientation, const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X, bool checkIfSingular=false )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`d` and :math:`X`).

Dot
---
Returns :math:`(x,y) = x^H y`. :math:`x` and :math:`y` are both allowed to be 
stored as column or row vectors, but will be interpreted as column vectors.

.. cpp:function:: T Dot( const Matrix<T>& x, const Matrix<T>& y )

   The serial implementation (templated over the datatype). 

.. cpp:function:: T Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`x` and :math:`y`).

Dotc
----
Same as ``Dot``. This routine name is provided since it is the usual 
BLAS naming convention.

.. cpp:function:: T Dotc( const Matrix<T>& x, const Matrix<T>& y )

   The serial implementation (templated over the datatype). 

.. cpp:function:: T Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`x` and :math:`y`).

Dotu
----
Returns :math:`x^T y`, which is **not** an inner product.

.. cpp:function:: T Dotu( const Matrix<T>& x, const Matrix<T>& y )

   The serial implementation (templated over the datatype). 

.. cpp:function:: T Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`x` and :math:`y`).

MakeTrapezoidal
---------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Sets all entries outside of the specified trapezoidal submatrix to zero.
The diagonal of the trapezoidal matrix is defined relative to either the 
upper-left or bottom-right corner of the matrix, depending on the 
value of ``side``; whether or not the trapezoid is upper or lower
(analogous to an upper or lower-triangular matrix) is determined by the 
``uplo`` parameter, and the last diagonal is defined with the ``offset`` 
integer.

.. cpp:function:: void MakeTrapezoidal( LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A )

   The serial implementation.

.. cpp:function:: void MakeTrapezoidal( LeftOrRight side, UpperOrLower uplo, int offset, DistMatrix<T,U,V>& A )

   The distributed implementation.

Nrm2
----
Returns :math:`||x||_2 = \sqrt{(x,x)} = \sqrt{x^H x}`. As with most other 
routines, even if :math:`x` is stored as a row vector, it will be interpreted
as a column vector.

.. cpp:function:: typename Base<F>::type Nrm2( const Matrix<F>& x )
.. cpp:function:: typename Base<F>::type Nrm2( const DistMatrix<F>& x )

Scal
----
:math:`X := \alpha X`.

.. cpp:function:: void Scal( T alpha, Matrix<T>& X )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Scal( T alpha, DistMatrix<T,U,V>& X )

   The distributed implementation (templated over the datatype and the 
   distribution of :math:`X`).

ScaleTrapezoid
--------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Scales the entries within the specified trapezoid of a general matrix.
The parameter conventions follow those of ``MakeTrapezoidal`` described above.

.. cpp:function:: void ScaleTrapezoid( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A )

   The serial implementation.

.. cpp:function:: void ScaleTrapezoid( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, DistMatrix<T,U,V>& A )

   The distributed implementation.

Transpose
---------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`B := A^T`. 

.. cpp:function:: void Transpose( const Matrix<T>& A, Matrix<T>& B )

   The serial version (templated over the datatype).

.. cpp:function:: void Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed version (templated over the datatype and the individual 
   distributions of :math:`A` and :math:`B`).

Zero
----
.. note::
   
   This is not a standard BLAS routine, but it is BLAS-like.

Sets all of the entries of the input matrix to zero.

.. cpp:function:: void Zero( Matrix<T>& A )

   The serial implementation.

.. cpp:function:: void Zero( DistMatrix<T,U,V>& A )

   The distributed implementation.


