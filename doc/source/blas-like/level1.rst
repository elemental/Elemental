Level 1
=======

The prototypes for the following routines can be found at 
`include/elemental/blas-like_decl.hpp <https://github.com/elemental/Elemental/tree/master/include/elemental/blas-like_decl.hpp>`_, while the
implementations are in `include/elemental/blas-like/level1/ <https://github.com/elemental/Elemental/tree/master/include/elemental/blas-like/level1>`_.

Adjoint
-------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`B := A^H`. 

.. cpp:function:: void Adjoint( const Matrix<T>& A, Matrix<T>& B )
.. cpp:function:: void Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

Axpy
----
Performs :math:`Y := \alpha X + Y` (hence the name *axpy*).

.. cpp:function:: void Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )
.. cpp:function:: void Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )

Conjugate
---------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`A := \bar A`. For real datatypes, this is a no-op.

.. cpp:function:: void Conjugate( Matrix<T>& A )
.. cpp:function:: void Conjugate( DistMatrix<T,U,V>& A )

:math:`B := \bar A`.

.. cpp:function:: void Conjugate( const Matrix<T>& A, Matrix<T>& B )
.. cpp:function:: void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

Copy
----
Sets :math:`Y := X`.

.. cpp:function:: void Copy( const Matrix<T>& X, Matrix<T>& Y )
.. cpp:function:: void Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

DiagonalScale
-------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Performs either :math:`X := \mbox{op}(D) X` or :math:`X := X \mbox{op}(D)`, 
where :math:`op(D)` equals :math:`D=D^T`, or :math:`D^H=\bar D`, where
:math:`D = \mbox{diag}(d)` and :math:`d` is a column vector.

.. cpp:function:: void DiagonalScale( LeftOrRight side, Orientation orientation, const Matrix<T>& d, Matrix<T>& X )
.. cpp:function:: void DiagonalScale( LeftOrRight side, Orientation orientation, const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )

DiagonalSolve
-------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Performs either :math:`X := \mbox{op}(D)^{-1} X` or 
:math:`X := X \mbox{op}(D)^{-1}`, where :math:`D = \mbox{diag}(d)` and :math:`d`
is a column vector.

.. cpp:function:: void DiagonalSolve( LeftOrRight side, Orientation orientation, const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular=false )
.. cpp:function:: void DiagonalSolve( LeftOrRight side, Orientation orientation, const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X, bool checkIfSingular=false )

Dot
---
Returns :math:`(x,y) = x^H y`. :math:`x` and :math:`y` are both allowed to be 
stored as column or row vectors, but will be interpreted as column vectors.

.. cpp:function:: T Dot( const Matrix<T>& x, const Matrix<T>& y )
.. cpp:function:: T Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,U,V>& y )

Dotc
----
Same as ``Dot``. This routine name is provided since it is the usual 
BLAS naming convention.

.. cpp:function:: T Dotc( const Matrix<T>& x, const Matrix<T>& y )
.. cpp:function:: T Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,U,V>& y )

Dotu
----
Returns :math:`x^T y`, which is **not** an inner product.

.. cpp:function:: T Dotu( const Matrix<T>& x, const Matrix<T>& y )
.. cpp:function:: T Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,U,V>& y )

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

.. cpp:function:: void MakeTrapezoidal( UpperOrLower uplo, Matrix<T>& A, int offset=0, LeftOrRight side=LEFT )
.. cpp:function:: void MakeTrapezoidal( UpperOrLower uplo, DistMatrix<T,U,V>& A, int offset=0, LeftOrRight side=LEFT )

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
.. cpp:function:: void Scal( T alpha, DistMatrix<T,U,V>& X )

ScaleTrapezoid
--------------
.. note::

   This is not a standard BLAS routine, but it is BLAS-like.

Scales the entries within the specified trapezoid of a general matrix.
The parameter conventions follow those of ``MakeTrapezoidal`` described above.

.. cpp:function:: void ScaleTrapezoid( T alpha, UpperOrLower uplo, Matrix<T>& A, int offset=0, LeftOrRight side=LEFT )
.. cpp:function:: void ScaleTrapezoid( T alpha, UpperOrLower uplo, DistMatrix<T,U,V>& A, int offset=0, LeftOrRight side=LEFT )

Transpose
---------
.. note:: 

   This is not a standard BLAS routine, but it is BLAS-like.

:math:`B := A^T` or :math:`B := A^H`. 

.. cpp:function:: void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate=false )
.. cpp:function:: void Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )

Zero
----
.. note::
   
   This is not a standard BLAS routine, but it is BLAS-like.

Sets all of the entries of the input matrix to zero.

.. cpp:function:: void Zero( Matrix<T>& A )
.. cpp:function:: void Zero( DistMatrix<T,U,V>& A )

SetDiagonal
-----------
.. note::
   
   This is not a standard BLAS routine.

Sets all of the diagonal entries of a matrix to a given value.

.. cpp:function:: void SetDiagonal( Matrix<T>& A, T alpha )
.. cpp:function:: void SetDiagonal( DistMatrix<T,U,V>& A, T alpha )

.. cpp:function:: void SetDiagonal( Matrix<T>& A, T alpha, int offset=0, LeftOrRight side=LEFT )
.. cpp:function:: void SetDiagonal( DistMatrix<T,U,V>& A, T alpha, int offset=0, LeftOrRight side=LEFT )

UpdateDiagonal
--------------
.. note::
   
   This is not a standard BLAS routine.

Adds a given value to the diagonal of a matrix.

.. cpp:function:: void UpdateDiagonal( Matrix<T>& A, T alpha )
.. cpp:function:: void UpdateDiagonal( DistMatrix<T,U,V>& A, T alpha )

.. cpp:function:: void UpdateDiagonal( Matrix<T>& A, T alpha, int offset=0, LeftOrRight side=LEFT )
.. cpp:function:: void UpdateDiagonal( DistMatrix<T,U,V>& A, T alpha, int offset=0, LeftOrRight side=LEFT )
