Level 1
=======

The prototypes for the following routines can be found at 
`include/elemental/basic.hpp <../../../../include/elemental/basic.hpp>`_, while the
implementations are in `include/elemental/basic/level1/ <../../../../include/elemental/basic/level1/>`_.

Adjoint
-------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`B := A^H`. 

.. cpp:function:: void basic::Adjoint( const Matrix<T>& A, Matrix<T>& B )

   The serial version (templated over the datatype).

.. cpp:function:: void basic::Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed version (templated over the datatype and the individual 
   distributions of :math:`A` and :math:`B`).

Axpy
----
Performs :math:`Y := \alpha X + Y` (hence the name *axpy*).

.. cpp:function:: void basic::Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )

   The distributed implementation (templated over the datatype and the shared
   distribution of :math:`A` and :math:`B`).

Conjugate
---------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`A := \bar A`. For real datatypes, this is a no-op.

.. cpp:function:: void basic::Conjugate( Matrix<T>& A )

   The serial version (templated over datatype).

.. cpp:function:: void basic::Conjugate( DistMatrix<T,U,V>& A )

   The distributed version (templated over the datatype and the distribution of
   :math:`A`).

:math:`B := \bar A`.

.. cpp:function:: void basic::Conjugate( const Matrix<T>& A, Matrix<T>& B )

   The serial version (templated over the datatype).

.. cpp:function:: void basic::Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed version (templated over the datatype and the individual 
   distributions of :math:`A` and :math:`B`).


Copy
----
Sets :math:`Y := X`.

.. cpp:function:: void basic::Copy( const Matrix<T>& X, Matrix<T>& Y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed implementation (templated over the datatype and the
   individual distributions of :math:`A` and :math:`B`).

DiagonalScale
-------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Performs either :math:`X := \mbox{op}(D) X` or :math:`X := X \mbox{op}(D)`, 
where :math:`op(D)` equals :math:`D=D^T`, or :math:`D^H=\bar D`, where
:math:`D = \mbox{diag}(d)` and :math:`d` is a column vector.

.. cpp:function:: void basic::DiagonalScale( Side side, Orientation orientation, const Matrix<T>& d, Matrix<T>& X )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::DiagonalScale( Side side, Orientation orientation, const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`d` and :math:`X`).

DiagonalSolve
-------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Performs either :math:`X := \mbox{op}(D)^{-1} X` or 
:math:`X := X \mbox{op}(D)^{-1}`, where :math:`D = \mbox{diag}(d)` and :math:`d`
is a column vector.

.. cpp:function:: void basic::DiagonalSolve( Side side, Orientation orientation, const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular=false )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::DiagonalSolve( Side side, Orientation orientation, const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X, bool checkIfSingular=false )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`d` and :math:`X`).

Dot
---
Returns :math:`(x,y) = x^H y`. :math:`x` and :math:`y` are both allowed to be 
stored as column or row vectors, but will be interpreted as column vectors.

.. cpp:function:: T basic::Dot( const Matrix<T>& x, const Matrix<T>& y )

   The serial implementation (templated over the datatype). 

.. cpp:function:: T basic::Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`x` and :math:`y`).

Dotc
----
Same as ``basic::Dot``. This routine name is provided since it is the usual 
BLAS naming convention.

.. cpp:function:: T basic::Dotc( const Matrix<T>& x, const Matrix<T>& y )

   The serial implementation (templated over the datatype). 

.. cpp:function:: T basic::Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`x` and :math:`y`).

Dotu
----
Returns :math:`x^T y`, which is **not** an inner product.

.. cpp:function:: T basic::Dotu( const Matrix<T>& x, const Matrix<T>& y )

   The serial implementation (templated over the datatype). 

.. cpp:function:: T basic::Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )

   The distributed implementation (templated over the datatype and the 
   individual distributions of :math:`x` and :math:`y`).

Nrm2
----
Returns :math:`||x||_2 = \sqrt{(x,x)} = \sqrt{x^H x}`. As with most other 
routines, even if :math:`x` is stored as a row vector, it will be interpreted
as a column vector.

.. cpp:function:: R basic::Nrm2( const Matrix<R>& x )

   Serial version for real datatypes.

.. cpp:function:: R basic::Nrm2( const Matrix<std::complex<R> >& x )

   Serial version for complex datatypes.

.. cpp:function:: R basic::Nrm2( const DistMatrix<R,MC,MR>& x )

   Distributed version for real datatypes.

.. cpp:function:: R basic::Nrm2( const DistMatrix<std::complex<R>,MC,MR>& x )

   Distributed version for complex datatypes.

Scal
----
:math:`X := \alpha X`.

.. cpp:function:: void basic::Scal( T alpha, Matrix<T>& X )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Scal( T alpha, DistMatrix<T,U,V>& X )

   The distributed implementation (templated over the datatype and the 
   distribution of :math:`X`).

Transpose
---------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`B := A^T`. 

.. cpp:function:: void basic::Transpose( const Matrix<T>& A, Matrix<T>& B )

   The serial version (templated over the datatype).

.. cpp:function:: void basic::Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

   The distributed version (templated over the datatype and the individual 
   distributions of :math:`A` and :math:`B`).

