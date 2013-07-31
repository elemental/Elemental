Level 3
=======

The prototypes for the following routines can be found at          
`include/elemental/blas-like_decl.hpp <https://github.com/elemental/Elemental/tree/master/include/elemental/blas-like_decl.hpp>`_, while the
implementations are in 
`include/elemental/blas-like/level3/ <https://github.com/elemental/Elemental/tree/master/include/elemental/blas-like/level3>`_.

Gemm
----
General matrix-matrix multiplication: updates
:math:`C := \alpha \mbox{op}_A(A) \mbox{op}_B(B) + \beta C`,
where :math:`\mbox{op}_A(M)` and :math:`\mbox{op}_B(M)` can each be chosen from 
:math:`M`, :math:`M^T`, and :math:`M^H`.

.. cpp:function:: void Gemm( Orientation orientationOfA, Orientation orientationOfB, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
.. cpp:function:: void Gemm( Orientation orientationOfA, Orientation orientationOfB, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C )

Hemm
----
Hermitian matrix-matrix multiplication: updates
:math:`C := \alpha A B + \beta C`, or 
:math:`C := \alpha B A + \beta C`, depending upon whether `side` is set to 
``LEFT`` or ``RIGHT``, respectively. In both of these types of updates, 
:math:`A` is implicitly Hermitian and only the triangle specified by `uplo` is 
accessed.

.. cpp:function:: void Hemm( LeftOrRight side, UpperOrLower uplo, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
.. cpp:function:: void Hemm( LeftOrRight side, UpperOrLower uplo, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C )

Her2k
-----
Hermitian rank-2K update: updates
:math:`C := \alpha (A B^H + B A^H) + \beta C`, or 
:math:`C := \alpha (A^H B + B^H A) + \beta C`, depending upon whether 
`orientation` is set to ``NORMAL`` or ``ADJOINT``, respectively. Only the 
triangle of :math:`C` specified by the `uplo` parameter is modified.

.. cpp:function:: void Her2k( UpperOrLower uplo, Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
.. cpp:function:: void Her2k( UpperOrLower uplo, Orientation orientation, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C )

Please see :cpp:func:`SetLocalTrr2kBlocksize\<T>` 
and :cpp:func:`LocalTrr2kBlocksize\<T>` in the 
:ref:`blas-tuning` section for information on tuning the distributed 
:cpp:func:`Her2k`.

Herk
----
Hermitian rank-K update: updates
:math:`C := \alpha A A^H + \beta C`, or 
:math:`C := \alpha A^H A + \beta C`, depending upon whether `orientation` is
set to ``NORMAL`` or ``ADJOINT``, respectively. Only the triangle of :math:`C` 
specified by the `uplo` parameter is modified.

.. cpp:function:: void Herk( UpperOrLower uplo, Orientation orientation, T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
.. cpp:function:: void Herk( UpperOrLower uplo, Orientation orientation, T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C )

Please see :cpp:func:`SetLocalTrrkBlocksize\<T>` 
and :cpp:func:`LocalTrrkBlocksize\<T>` in the :ref:`blas-tuning`
section for information on tuning the distributed :cpp:func:`Herk`.

Symm
----
Symmetric matrix-matrix multiplication: updates
:math:`C := \alpha A B + \beta C`, or 
:math:`C := \alpha B A + \beta C`, depending upon whether `side` is set to 
``LEFT`` or ``RIGHT``, respectively. In both of these types of updates, 
:math:`A` is implicitly symmetric and only the triangle specified by `uplo` 
is accessed.

.. cpp:function:: void Symm( LeftOrRight side, UpperOrLower uplo, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C, bool conjugate=false )
.. cpp:function:: void Symm( LeftOrRight side, UpperOrLower uplo, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C, bool conjugate=false )

Syr2k
-----
Symmetric rank-2K update: updates
:math:`C := \alpha (A B^T + B A^T) + \beta C`, or 
:math:`C := \alpha (A^T B + B^T A) + \beta C`, depending upon whether 
`orientation` is set to ``NORMAL`` or ``TRANSPOSE``, respectively. Only the 
triangle of :math:`C` specified by the `uplo` parameter is modified.

.. cpp:function:: void Syr2k( UpperOrLower uplo, Orientation orientation, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
.. cpp:function:: void Syr2k( UpperOrLower uplo, Orientation orientation, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C )

Please see :cpp:func:`SetLocalTrr2kBlocksize\<T>` 
and :cpp:func:`LocalTrr2kBlocksize\<T>` in the 
:ref:`blas-tuning` section for information on tuning the distributed 
:cpp:func:`Syr2k`.

Syrk
----
Symmetric rank-K update: updates
:math:`C := \alpha A A^T + \beta C`, or 
:math:`C := \alpha A^T A + \beta C`, depending upon whether `orientation` is
set to ``NORMAL`` or ``TRANSPOSE``, respectively. Only the triangle of :math:`C`
specified by the `uplo` parameter is modified.

.. cpp:function:: void Syrk( UpperOrLower uplo, Orientation orientation, T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
.. cpp:function:: void Syrk( UpperOrLower uplo, Orientation orientation, T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C )

Please see :cpp:func:`SetLocalTrrkBlocksize\<T>` 
and :cpp:func:`LocalTrrkBlocksize\<T>` in the :ref:`blas-tuning`
section for information on tuning the distributed :cpp:func:`Syrk`.

Trmm
----
Triangular matrix-matrix multiplication: performs
:math:`C := \alpha \mbox{op}(A) B`, or 
:math:`C := \alpha B \mbox{op}(A)`, depending upon whether `side` was chosen
to be ``LEFT`` or ``RIGHT``, respectively. Whether :math:`A` is treated as 
lower or upper triangular is determined by `uplo`, and :math:`\mbox{op}(A)` 
can be any of :math:`A`, :math:`A^T`, and :math:`A^H` (and `diag` determines
whether :math:`A` is treated as unit-diagonal or not).

.. cpp:function:: void Trmm( LeftOrRight side, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, T alpha, const Matrix<T>& A, Matrix<T>& B )
.. cpp:function:: void Trmm( LeftOrRight side, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, T alpha, const DistMatrix<T>& A, DistMatrix<T>& B )

Trr2k
-----
Triangular rank-2k update: performs 
:math:`E := \alpha ( \mbox{op}(A) \mbox{op}(B) + \mbox{op}(C) \mbox{op}(D) ) + \beta E`,
where only the triangle of `E` specified by `uplo` is modified, and 
:math:`\mbox{op}(X)` is determined by `orientationOfX`, for each 
:math:`X \in \left\{A,B,C,D\right\}`.

.. note::

   There is no corresponding BLAS routine, but it is a natural generalization
   of "symmetric" and "Hermitian" updates.

.. cpp:function:: void Trr2k( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB, Orientation orientationOfC, Orientation orientationOfD, T alpha, const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C, const Matrix<T>& D, T beta, Matrix<T>& E )
.. cpp:function:: void Trr2k( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB, Orientation orientationOfC, Orientation orientationOfD, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, const DistMatrix<T>& C, const DistMatrix<T>& D, T beta, DistMatrix<T>& E )

Trrk
----
Triangular rank-k update: performs 
:math:`C := \alpha \mbox{op}(A) \mbox{op}(B) + \beta C`, where only the 
triangle of `C` specified by `uplo` is modified, and :math:`\mbox{op}(A)` and 
:math:`\mbox{op}(B)` are determined by `orientationOfA` and `orientationOfB`, 
respectively.

.. note::

   There is no corresponding BLAS routine, but this type of update is frequently
   encountered, even in serial. For instance, the symmetric rank-k update 
   performed during an LDL factorization is symmetric but one of the 
   two update matrices is scaled by D.

.. cpp:function:: void Trrk( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB, T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
.. cpp:function:: void Trrk( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C )

Trtrmm
------
.. note:: 

   This routine loosely corresponds with the LAPACK routines ?lauum.

Symmetric/Hermitian triangular matrix-matrix multiply: performs
:math:`L := L^T L`, :math:`L := L^H L`, :math:`U := U U^T`, or 
:math:`U := U U^H`, depending upon the choice of the `orientation` and 
`uplo` parameters. 

.. cpp:function:: void Trtrmm( Orientation orientation, UpperOrLower uplo, Matrix<T>& A )
.. cpp:function:: void Trtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<T>& A )

Trdtrmm
-------
.. note:: 

   This is a modification of Trtrmm for LDL factorizations.

Symmetric/Hermitian triangular matrix-matrix multiply (with diagonal scaling): 
performs :math:`L := L^T D^{-1} L`, :math:`L := L^H D^{-1} L`, 
:math:`U := U D^{-1} U^T`, or :math:`U := U D^{-1} U^H`, 
depending upon the choice of the `orientation` and `uplo` parameters. 
Note that :math:`L` and :math:`U` are unit-diagonal and their diagonal is 
overwritten with :math:`D`.

.. cpp:function:: void Trdtrmm( Orientation orientation, UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void Trdtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<F>& A )

Trsm
----
Triangular solve with multiple right-hand sides: performs
:math:`C := \alpha \mbox{op}(A)^{-1} B`, or 
:math:`C := \alpha B \mbox{op}(A)^{-1}`, depending upon whether `side` was 
chosen to be ``LEFT`` or ``RIGHT``, respectively. Whether :math:`A` is treated 
as lower or upper triangular is determined by `uplo`, and :math:`\mbox{op}(A)` 
can be any of :math:`A`, :math:`A^T`, and :math:`A^H` (and `diag` determines
whether :math:`A` is treated as unit-diagonal or not).

.. cpp:function:: void Trsm( LeftOrRight side, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, F alpha, const Matrix<F>& A, Matrix<F>& B )

.. cpp:function:: void Trsm( LeftOrRight side, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag, F alpha, const DistMatrix<F>& A, DistMatrix<F>& B )

Two-sided Trmm
--------------
Performs a two-sided triangular multiplication with multiple right-hand sides 
which preserves the symmetry of the input matrix, 
either :math:`A := L^H A L` or :math:`A := U A U^H`.

.. cpp:function:: void TwoSidedTrmm( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<T>& A, const Matrix<T>& B )

.. cpp:function:: void TwoSidedTrmm( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<T>& A, const DistMatrix<T>& B )

Two-sided Trsm
--------------
Performs a two-sided triangular solves with multiple right-hand sides which 
preserves the symmetry of the input matrix, 
either :math:`A := L^{-1} A L^{-H}` or :math:`A := U^{-H} A U^{-1}`.

.. cpp:function:: void TwoSidedTrsm( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B )
.. cpp:function:: void TwoSidedTrsm( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& B )
