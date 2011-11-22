Level 3
=======

Gemm
----
General matrix-matrix multiplication: updates
:math:`C := \alpha \mbox{op}_A(A) \mbox{op}_B(B) + \beta C`,
where :math:`\mbox{op}_A(M)` and :math:`\mbox{op}_B(M)` can each be chosen from 
:math:`M`, :math:`M^T`, and :math:`M^H`.

.. cpp:function:: void basic::Gemm( Orientation orientationOfA, Orientation orientationOfB, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Gemm( Orientation orientationOfA, Orientation orientationOfB, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Hemm
----
Hermitian matrix-matrix multiplication: updates
:math:`C := \alpha A B + \beta C`, or 
:math:`C := \alpha B A + \beta C`, depending upon whether `side` is set to 
``LEFT`` or ``RIGHT``, respectively. In both of these types of updates, 
:math:`A` is implicitly Hermitian and only the triangle specified by `shape` is 
accessed.

.. cpp:function:: void basic::Hemm( Side side, Shape shape, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Hemm( Side side, Shape shape, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Her2k
-----
Hermitian rank-2K update: updates
:math:`C := \alpha (A B^H + B A^H) + \beta C`, or 
:math:`C := \alpha (A^H B + B^H A) + \beta C`, depending upon whether 
`orientation` is set to ``NORMAL`` or ``ADJOINT``, respectively. Only the 
triangle of :math:`C` specified by the `shape` parameter is modified.

.. cpp:function:: void basic::Her2k( Shape shape, Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Her2k( Shape shape, Orientation orientation, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Please see ``basic::SetLocalTrr2kBlocksize<T>( int blocksize )`` 
and ``int basic::LocalTrr2kBlocksize<T>()`` in the 
*Tuning parameters* section for information on tuning the distributed 
``basic::Her2k``.

Herk
----
Hermitian rank-K update: updates
:math:`C := \alpha A A^H + \beta C`, or 
:math:`C := \alpha A^H A + \beta C`, depending upon whether `orientation` is
set to ``NORMAL`` or ``ADJOINT``, respectively. Only the triangle of :math:`C` 
specified by the `shape` parameter is modified.

.. cpp:function:: void basic::Herk( Shape shape, Orientation orientation, T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Herk( Shape shape, Orientation orientation, T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Please see ``basic::SetLocalTrrkBlocksize<T>( int blocksize )`` 
and ``int basic::LocalTrrkBlocksize<T>()`` in the *Tuning parameters*
section for information on tuning the distributed ``basic::Herk``.

Hetrmm
------
.. note:: 

   This routine directly corresponds with the LAPACK routines ?lauum, but it 
   only involves matrix-matrix multiplication, so it is lumped in with the 
   BLAS-like routines in Elemental.

Hermitian triangular matrix-matrix multiply: performs 
:math:`L := L^H L` or :math:`U := U U^H`, depending upon the choice of the 
`shape` parameter. 

.. cpp:function:: void basic::Hetrmm( Shape shape, Matrix<T>& A )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Hetrmm( Shape shape, DistMatrix<T,MC,MR>& A )

   The distributed implementation (templated over the datatype).

Symm
----
Symmetric matrix-matrix multiplication: updates
:math:`C := \alpha A B + \beta C`, or 
:math:`C := \alpha B A + \beta C`, depending upon whether `side` is set to 
``LEFT`` or ``RIGHT``, respectively. In both of these types of updates, 
:math:`A` is implicitly symmetric and only the triangle specified by `shape` 
is accessed.

.. cpp:function:: void basic::Symm( Side side, Shape shape, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Symm( Side side, Shape shape, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Syr2k
-----
Symmetric rank-2K update: updates
:math:`C := \alpha (A B^T + B A^T) + \beta C`, or 
:math:`C := \alpha (A^T B + B^T A) + \beta C`, depending upon whether 
`orientation` is set to ``NORMAL`` or ``TRANSPOSE``, respectively. Only the 
triangle of :math:`C` specified by the `shape` parameter is modified.

.. cpp:function:: void basic::Syr2k( Shape shape, Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Syr2k( Shape shape, Orientation orientation, T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Please see ``basic::SetLocalTrr2kBlocksize<T>( int blocksize )`` 
and ``int basic::LocalTrr2kBlocksize<T>()`` in the 
*Tuning parameters* section for information on tuning the distributed 
``basic::Syr2k``.

Syrk
----
Symmetric rank-K update: updates
:math:`C := \alpha A A^T + \beta C`, or 
:math:`C := \alpha A^T A + \beta C`, depending upon whether `orientation` is
set to ``NORMAL`` or ``TRANSPOSE``, respectively. Only the triangle of :math:`C`
specified by the `shape` parameter is modified.

.. cpp:function:: void basic::Syrk( Shape shape, Orientation orientation, T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Syrk( Shape shape, Orientation orientation, T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C )

   The distributed implementation (templated over the datatype).

Please see ``basic::SetLocalTrrkBlocksize<T>( int blocksize )`` 
and ``int basic::LocalTrrkBlocksize<T>()`` in the *Tuning parameters*
section for information on tuning the distributed ``basic::Syrk``.

Trmm
----
Triangular matrix-matrix multiplication: performs
:math:`C := \alpha \mbox{op}(A) B`, or 
:math:`C := \alpha B \mbox{op}(A)`, depending upon whether `side` was chosen
to be ``LEFT`` or ``RIGHT``, respectively. Whether :math:`A` is treated as 
lower or upper triangular is determined by `shape`, and :math:`\mbox{op}(A)` 
can be any of :math:`A`, :math:`A^T`, and :math:`A^H` (and `diagonal` determines
whether :math:`A` is treated as unit diagonal or not).

.. cpp:function:: void basic::Trmm( Side side, Shape shape, Orientation orientation, Diagonal diagonal, T alpha, const Matrix<T>& A, Matrix<T>& B )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Trmm( Side side, Shape shape, Orientation orientation, Diagonal diagonal, T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )

   The distributed implementation (templated over the datatype).

Trsm
----
Triangular solve with multiple right-hand sides: performs
:math:`C := \alpha \mbox{op}(A)^{-1} B`, or 
:math:`C := \alpha B \mbox{op}(A)^{-1}`, depending upon whether `side` was 
chosen to be ``LEFT`` or ``RIGHT``, respectively. Whether :math:`A` is treated 
as lower or upper triangular is determined by `shape`, and :math:`\mbox{op}(A)` 
can be any of :math:`A`, :math:`A^T`, and :math:`A^H` (and `diagonal` determines
whether :math:`A` is treated as unit diagonal or not).

.. cpp:function:: void basic::Trsm( Side side, Shape shape, Orientation orientation, Diagonal diagonal, T alpha, const Matrix<T>& A, Matrix<T>& B )

   The serial implementation (templated over the datatype).

.. cpp:function:: void basic::Trsm( Side side, Shape shape, Orientation orientation, Diagonal diagonal, T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )

   The distributed implementation (templated over the datatype).
