Level 2
=======

The prototypes for the following routines can be found at          
`include/elemental/blas-like_decl.hpp <https://github.com/elemental/Elemental/tree/master/include/elemental/blas-like_decl.hpp>`_, while the
implementations are in `include/elemental/blas-like/level2/ <https://github.com/elemental/Elemental/tree/master/include/elemental/blas-like/level2>`_.

Gemv
----
General matrix-vector multiply:
:math:`y := \alpha \mbox{op}(A) x + \beta y`,
where :math:`\mbox{op}(A)` can be :math:`A`, :math:`A^T`, or :math:`A^H`.
Whether or not :math:`x` and :math:`y` are stored as row vectors, they will
be interpreted as column vectors.

.. cpp:function:: void Gemv( Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
.. cpp:function:: void Gemv( Orientation orientation, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y )

Ger
---
General rank-one update: :math:`A := \alpha x y^H + A`. :math:`x` and :math:`y`
are free to be stored as either row or column vectors, but they will be 
interpreted as column vectors.

.. cpp:function:: void Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
.. cpp:function:: void Ger( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

Gerc
----
This is the same as :cpp:func:`Ger`, but the name is provided because it exists
in the BLAS.

.. cpp:function:: void Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
.. cpp:function:: void Gerc( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

Geru
----
General rank-one update (unconjugated): :math:`A := \alpha x y^T + A`. :math:`x`
and :math:`y` are free to be stored as either row or column vectors, but they 
will be interpreted as column vectors.

.. cpp:function:: void Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
.. cpp:function:: void Geru( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

Hemv
----
Hermitian matrix-vector multiply: :math:`y := \alpha A x + \beta y`, where 
:math:`A` is Hermitian.

.. cpp:function:: void Hemv( UpperOrLower uplo, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
.. cpp:function:: void Hemv( UpperOrLower uplo, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y )

Please see :cpp:func:`SetLocalSymvBlocksize\<T>` and 
:cpp:func:`LocalSymvBlocksize\<T>` in the :ref:`blas-tuning` section for 
information on tuning the distributed :cpp:func:`Hemv`.

Her
---
Hermitian rank-one update: implicitly performs :math:`A := \alpha x x^H + A`, 
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
.. cpp:function:: void Her( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A )

Her2
----
Hermitian rank-two update: implicitly performs 
:math:`A := \alpha ( x y^H + y x^H ) + A`,
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Her2( UpperOrLower uplo, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
.. cpp:function:: void Her2( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A )

Symv
----
Symmetric matrix-vector multiply: :math:`y := \alpha A x + \beta y`, where 
:math:`A` is symmetric.

.. cpp:function:: void Symv( UpperOrLower uplo, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y, bool conjugate=false )
.. cpp:function:: void Symv( UpperOrLower uplo, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y, bool conjugate=false )

Please see :cpp:func:`SetLocalSymvBlocksize\<T>` and 
:cpp:func:`LocalSymvBlocksize\<T>` in the :ref:`blas-tuning` section for 
information on tuning the distributed :cpp:func:`Symv`.

Syr
---
Symmetric rank-one update: implicitly performs :math:`A := \alpha x x^T + A`, 
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Syr( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A, bool conjugate=false )
.. cpp:function:: void Syr( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A, bool conjugate=false )

Syr2
----
Symmetric rank-two update: implicitly performs 
:math:`A := \alpha ( x y^T + y x^T ) + A`,
where only the triangle of :math:`A` specified by `uplo` is updated.

.. cpp:function:: void Syr2( UpperOrLower uplo, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A, bool conjugate=false )
.. cpp:function:: void Syr2( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A, bool conjugate=false )

Trmv
----
Not yet written. Please call :cpp:func:`Trmm` for now.

Trsv
----
Triangular solve with a vector: computes
:math:`x := \mbox{op}(A)^{-1} x`, where :math:`\mbox{op}(A)` is either 
:math:`A`, :math:`A^T`, or :math:`A^H`, and :math:`A` is treated an either a 
lower or upper triangular matrix, depending upon `uplo`. :math:`A` can also be 
treated as implicitly having a unit-diagonal if `diag` is set to ``UNIT``.

.. cpp:function:: void Trsv( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, const Matrix<F>& A, Matrix<F>& x )
.. cpp:function:: void Trsv( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, const DistMatrix<F>& A, DistMatrix<F>& x )
