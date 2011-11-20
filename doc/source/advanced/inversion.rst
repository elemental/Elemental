In-place inversion
==================

General inversion
-----------------
Not yet written, but relatively trivial and planned.

HPD inversion
-------------
This routine uses a custom algorithm for computing the inverse of a
Hermitian positive-definite matrix :math:`A` as

.. math::
   :nowrap:

   \[
   A^{-1} = (L L^H)^{-1} = L^{-H} L^{-1}, 
   \]

where :math:`L` is the lower Cholesky factor of :math:`A` (the upper Cholesky
factor is computed in the case of upper-triangular storage). Rather than 
performing Cholesky factorization, triangular inversion, and then the Hermitian
triangular outer product in sequence, this routine makes use of the single-sweep 
algorithm described in Bientinesi et al.'s "Families of algorithms related to the 
inversion of a symmetric positive definite matrix", in particular, the variant 2
algorithm from Fig. 9. 

If the matrix is found to not be HPD, then a ``NonHPDMatrixException`` is thrown.

.. cpp:function:: void advanced::HPDInverse( Shape shape, Matrix<F>& A )

   Overwrite the `shape` triangle of the HPD matrix `A` with the same 
   triangle of the inverse of `A`.

.. cpp:function:: void advanced::HPDInverse( Shape shape, DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.


Triangular inversion
--------------------
Inverts a (possibly unit-diagonal) triangular matrix in-place.

.. cpp:function:: void advanced::TriangularInverse( Shape shape, Diagonal diagonal, Matrix<F>& A )

   Inverts the triangle of `A` specified by the parameter `shape`; 
   if `diagonal` is set to `UNIT`, then `A` is treated as unit-diagonal.

.. cpp:function:: void advanced::TriangularInverse( Shape shape, Diagonal diagonal, DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.
