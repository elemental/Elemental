Linear solvers
==============

HPD solve
---------
Solves either :math:`AX=B` or :math:`A^T X=B` for :math:`X` given Hermitian 
positive-definite (HPD) :math:`A` and right-hand side matrix :math:`B`. 
The solution is computed by first finding the Cholesky factorization of 
:math:`A` and then performing two successive triangular solves against 
:math:`B`.

.. cpp:function:: void HPDSolve( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, Matrix<F>& B )
.. cpp:function:: void HPDSolve( UpperOrLower uplo, Orientation orientation, DistMatrix<F>& A, DistMatrix<F>& B )

   Overwrite `B` with the solution to :math:`AX=B` or :math:`A^T X=B`, 
   where `A` is Hermitian positive-definite and only the triangle of `A` 
   specified by `uplo` is accessed.

Gaussian elimination
--------------------
Solves :math:`AX=B` for :math:`X` given a general square nonsingular matrix 
:math:`A` and right-hand side matrix :math:`B`. The solution is computed through
(partially pivoted) Gaussian elimination.

.. cpp:function:: void GaussianElimination( Matrix<F>& A, Matrix<F>& B )
.. cpp:function:: void GaussianElimination( DistMatrix<F>& A, DistMatrix<F>& B )

   Upon completion, :math:`A` will have been overwritten with Gaussian 
   elimination and :math:`B` will be overwritten with :math:`X`.

Least-squares
-------------
Solves :math:`AX=B` or :math:`A^H X = B` for :math:`X` in a least-squares sense 
given a general full-rank matrix :math:`A \in \mathbb{F}^{m \times n}`. 
If :math:`m \ge n`, then the first step is to form the QR factorization of 
:math:`A`, otherwise the LQ factorization is computed. 

* If solving :math:`AX=B`, then either :math:`X=R^{-1} Q^H B` or 
  :math:`X=Q^H L^{-1} B`.

* If solving :math:`A^H X=B`, then either :math:`X=Q R^{-H} B` or 
  :math:`X=L^{-H} Q B`.

.. cpp:function:: void LeastSquares( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, Matrix<F>& X )
.. cpp:function:: void LeastSquares( Orientation orientation, DistMatrix<F>& A, const DistMatrix<F>& B, DistMatrix<F>& X )

   If `orientation` is set to ``NORMAL``, then solve :math:`AX=B`, otherwise 
   `orientation` must be equal to ``ADJOINT`` and :math:`A^H X=B` will 
   be solved. Upon completion, :math:`A` is overwritten with its QR or LQ 
   factorization, and :math:`X` is overwritten with the solution.

Solve after Cholesky
--------------------
Uses an existing in-place Cholesky factorization to solve against one or more 
right-hand sides.

.. cpp:function:: void cholesky::SolveAfter( UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
.. cpp:function:: void cholesky::SolveAfter( UpperOrLower uplo, Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B )

   Update :math:`B := A^{-1} B`, :math:`B := A^{-T} B`, or 
   :math:`B := A^{-H} B`, where one triangle of :math:`A` has been overwritten 
   with its Cholesky factor.

Solve after LU
--------------
Uses an existing in-place LU factorization (with or without partial pivoting) 
to solve against one or more right-hand sides.

.. cpp:function:: void lu::SolveAfter( Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
.. cpp:function:: void lu::SolveAfter( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B )

   Update :math:`B := A^{-1} B`, :math:`B := A^{-T} B`, or 
   :math:`B := A^{-H} B`, where :math:`A` has been overwritten with its LU 
   factors (without partial pivoting).

.. cpp:function:: void lu::SolveAfter( Orientation orientation, const Matrix<F>& A, const Matrix<int>& p, Matrix<F>& B )
.. cpp:function:: void lu::SolveAfter( Orientation orientation, const DistMatrix<F>& A, const DistMatrix<int,VC,STAR>& p, DistMatrix<F>& B )

   Update :math:`B := A^{-1} B`, :math:`B := A^{-T} B`, or 
   :math:`B := A^{-H} B`, where :math:`A` has been overwritten with 
   its LU factors with partial pivoting, which satisfy :math:`P A = L U`, where
   the permutation matrix :math:`P` is represented by the pivot vector ``p``.
