Level 2
=======

Gemv
----
General matrix-vector multiply:
:math:`y := \alpha \mbox{op}(A) x + \beta y`,
where :math:`\mbox{op}(A)` can be :math:`A`, :math:`A^T`, or :math:`A^H`.
Whether or not :math:`x` and :math:`y` are stored as row vectors, they will
be interpreted as column vectors.

.. cpp:function:: void basic::Gemv( Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )

   Serial implementation (templated over the datatype).

.. cpp:function:: void basic::Gemv( Orientation orientation, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x, T beta, DistMatrix<T,MC,MR>& y )

   Distributed implementation (templated over the datatype).

Ger
---
General rank-one update: :math:`A := \alpha x y^H + A`. :math:`x` and :math:`y`
are free to be stored as either row or column vectors, but they will be 
interpreted as column vectors.

.. cpp:function:: void basic::Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Ger( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype). 

Gerc
----
This is the same as ``basic::Ger``, but the name is provided because it exists
in the BLAS.

.. cpp:function:: void basic::Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Gerc( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype). 

Geru
----
General rank-one update (unconjugated): :math:`A := \alpha x y^T + A`. :math:`x` and :math:`y`
are free to be stored as either row or column vectors, but they will be 
interpreted as column vectors.

.. cpp:function:: void basic::Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Geru( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype). 

Hemv
----
Hermitian matrix-vector multiply: :math:`y := \alpha A x + \beta y`, where 
:math:`A` is Hermitian.

.. cpp:function:: void basic::Hemv( Shape shape, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Hemv( Shape shape, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x, T beta, DistMatrix<T,MC,MR>& y )

   The distributed implementation (templated over the datatype).

Please see ``basic::SetLocalHemvBlocksize<T>( int blocksize )`` and 
``int basic::LocalHemvBlocksize<T>()`` in the *Tuning parameters* section for 
information on tuning the distributed ``basic::Hemv``.

Her
---
Hermitian rank-one update: implicitly performs :math:`A := \alpha x x^H + A`, 
where only the triangle of :math:`A` specified by `shape` is updated.

.. cpp:function:: void basic::Her( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Her( Shape shape, T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype).

Her2
----
Hermitian rank-two update: implicitly performs 
:math:`A := \alpha ( x y^H + y x^H ) + A`,
where only the triangle of :math:`A` specified by `shape` is updated.

.. cpp:function:: void basic::Her2( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Her2( Shape shape, T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype).

Symv
----
Symmetric matrix-vector multiply: :math:`y := \alpha A x + \beta y`, where 
:math:`A` is symmetric.

.. cpp:function:: void basic::Symv( Shape shape, T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Symv( Shape shape, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x, T beta, DistMatrix<T,MC,MR>& y )

   The distributed implementation (templated over the datatype).

Please see ``basic::SetLocalSymvBlocksize<T>( int blocksize )`` and 
``int basic::LocalSymvBlocksize<T>()`` in the *Tuning parameters* section for 
information on tuning the distributed ``basic::Symv``.

Syr
---
Symmetric rank-one update: implicitly performs :math:`A := \alpha x x^T + A`, 
where only the triangle of :math:`A` specified by `shape` is updated.

.. cpp:function:: void basic::Syr( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Syr( Shape shape, T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype).

Syr2
----
Symmetric rank-two update: implicitly performs 
:math:`A := \alpha ( x y^T + y x^T ) + A`,
where only the triangle of :math:`A` specified by `shape` is updated.

.. cpp:function:: void basic::Syr2( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Syr2( Shape shape, T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype).

Trmv
----
Not yet written. Please call Trmm.

Trsv
----
Triangular solve with a vector: computes
:math:`x := \mbox{op}(A)^{-1} x`, where :math:`\mbox{op}(A)` is either 
:math:`A`, :math:`A^T`, or :math:`A^H`, and :math:`A` is treated an either a 
lower or upper triangular matrix, depending upon `shape`. :math:`A` can also be 
treated as implicitly having a unit diagonal if `diagonal` is set to ``UNIT``.

.. cpp:function:: void basic::Trsv( Shape shape, Orientation orientation, Diagonal diagonal, const Matrix<F>& A, Matrix<F>& x )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Trsv( Shape shape, Orientation orientation, Diagonal diagonal, const DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& x )

   The distributed implementation (templated over the datatype).

