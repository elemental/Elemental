Linear solvers
==============

Cholesky solve
--------------
Solves :math:`AX=B` for :math:`X` given Hermitian positive-definite (HPD) 
:math:`A` and right-hand side matrix :math:`B`. The solution is computed by 
first finding the Cholesky factorization of :math:`A` and then performing two
successive triangular solves against :math:`B`:

.. math::

   B := A^{-1} B = (L L^H)^{-1} B = L^{-H} L^{-1} B


.. cpp:function:: void CholeskySolve( UpperOrLower uplo, DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B )

   Overwrite `B` with the solution to :math:`AX=B`, where `A` is Hermitian 
   positive-definite and only the triangle of `A` specified by `uplo` is 
   accessed.

Gaussian elimination
--------------------
Solves :math:`AX=B` for :math:`X` given a general square nonsingular matrix 
:math:`A` and right-hand side matrix :math:`B`. The solution is computed through
(partially pivoted) Gaussian elimination.

.. cpp:function:: void GaussianElimination( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B )

   Upon completion, :math:`A` will have been overwritten with Gaussian elimination
   and :math:`B` will be overwritten with :math:`X`.

Householder solve
-----------------
Solves :math:`AX=B` or :math:`A^H X = B` for :math:`X` in a least-squares sense 
given a general full-rank matrix :math:`A \in \mathbb{F}^{m \times n}`. 
If :math:`m \ge n`, then the first step is to form the QR factorization of 
:math:`A`, otherwise the LQ factorization is computed. 

* If solving :math:`AX=B`, then either :math:`X=R^{-1} Q^H B` or 
  :math:`X=Q^H L^{-1} B`.

* If solving :math:`A^H X=B`, then either :math:`X=Q R^{-H} B` or 
  :math:`X=L^{-H} Q B`.

.. cpp:function:: void HouseholderSolve( Orientation orientation, DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& B, DistMatrix<F,MC,MR>& X )

   If `orientation` is set to ``NORMAL``, then solve :math:`AX=B`, otherwise 
   `orientation` must be equal to ``ADJOINT`` and :math:`A^H X=B` will 
   be solved. Upon completion, :math:`A` is overwritten with its QR or LQ 
   factorization, and :math:`X` is overwritten with the solution.
