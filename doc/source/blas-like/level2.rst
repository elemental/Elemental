Level 2
=======

The prototypes for the following routines can be found at          
`include/elemental/blas-like.hpp <../../../../include/elemental/blas-like.hpp>`_, while the
implementations are in `include/elemental/blas-like/level2/ <../../../../include/elemental/b
asic/level2/>`_.

Gemv
----
General matrix-vector multiply:
:math:`y := \alpha \mbox{op}(A) x + \beta y`,
where :math:`\mbox{op}(A)` can be :math:`A`, :math:`A^T`, or :math:`A^H`.
Whether or not :math:`x` and :math:`y` are stored as row vectors, they will
be interpreted as column vectors.

.. cpp:function:: void Gemv( Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )

   Serial implementation (templated over the datatype).

.. cpp:function:: void Gemv( Orientation orientation, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y )

   Distributed implementation (templated over the datatype).

Ger
---
General rank-one update: :math:`A := \alpha x y^H + A`. :math:`x` and :math:`y`
are free to be stored as either row or column vectors, but they will be 
interpreted as column vectors.

.. cpp:function:: void Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Ger( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype). 

Gerc
----
This is the same as ``Ger``, but the name is provided because it exists
in the BLAS.

.. cpp:function:: void Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Gerc( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype). 

Geru
----
General rank-one update (unconjugated): :math:`A := \alpha x y^T + A`. :math:`x` and :math:`y`
are free to be stored as either row or column vectors, but they will be 
interpreted as column vectors.

.. cpp:function:: void Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Geru( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype). 

Hemv
----
Hermitian matrix-vector multiply: :math:`y := \alpha A x + \beta y`, where 
:math:`A` is Hermitian.

.. cpp:function:: void Hemv( UpperOrLower uplo, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Hemv( UpperOrLower uplo, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y )

   The distributed implementation (templated over the datatype).

Please see ``SetLocalHemvBlocksize<T>( int blocksize )`` and 
``int LocalHemvBlocksize<T>()`` in the *Tuning parameters* section for 
information on tuning the distributed ``Hemv``.

Her
---
Hermitian rank-one update: implicitly performs :math:`A := \alpha x x^H + A`, 
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Her( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype).

Her2
----
Hermitian rank-two update: implicitly performs 
:math:`A := \alpha ( x y^H + y x^H ) + A`,
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Her2( UpperOrLower uplo, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Her2( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype).

Symv
----
Symmetric matrix-vector multiply: :math:`y := \alpha A x + \beta y`, where 
:math:`A` is symmetric.

.. cpp:function:: void Symv( UpperOrLower uplo, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Symv( UpperOrLower uplo, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y )

   The distributed implementation (templated over the datatype).

Please see ``SetLocalSymvBlocksize<T>( int blocksize )`` and 
``int LocalSymvBlocksize<T>()`` in the *Tuning parameters* section for 
information on tuning the distributed ``Symv``.

Syr
---
Symmetric rank-one update: implicitly performs :math:`A := \alpha x x^T + A`, 
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Syr( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Syr( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype).

Syr2
----
Symmetric rank-two update: implicitly performs 
:math:`A := \alpha ( x y^T + y x^T ) + A`,
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Syr2( UpperOrLower uplo, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Syr2( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

   The distributed implementation (templated over the datatype).

Trmv
----
Not yet written. Please call Trmm.

Trsv
----
Triangular solve with a vector: computes
:math:`x := \mbox{op}(A)^{-1} x`, where :math:`\mbox{op}(A)` is either 
:math:`A`, :math:`A^T`, or :math:`A^H`, and :math:`A` is treated an either a 
lower or upper triangular matrix, depending upon `uplo`. :math:`A` can also be 
treated as implicitly having a unit-diagonal if `diag` is set to ``UNIT``.

.. cpp:function:: void Trsv( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, const Matrix<F>& A, Matrix<F>& x )

   The serial implementation (templated over the datatype).

.. cpp:function:: void Trsv( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, const DistMatrix<F>& A, DistMatrix<F>& x )

   The distributed implementation (templated over the datatype).

