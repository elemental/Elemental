Factorization-based inversion
=============================

General inversion
-----------------
This routine computes the in-place inverse of a general fully-populated 
(invertible) matrix :math:`A` as

.. math::

   A^{-1} = U^{-1} L^{-1} P,

where :math:`PA=LU` is the result of LU factorization with partial pivoting.
The algorithm essentially factors :math:`A`, inverts :math:`U` in place, 
solves against :math:`L` one block column at a time, and then applies the 
row pivots in reverse order to the columns of the result.

.. cpp:function:: void Inverse( Matrix<F>& A )
.. cpp:function:: void Inverse( DistMatrix<F>& A )

   Overwrites the general matrix `A` with its inverse.

HPD inversion
-------------
This routine uses a custom algorithm for computing the inverse of a
Hermitian positive-definite matrix :math:`A` as

.. math::

   A^{-1} = (L L^H)^{-1} = L^{-H} L^{-1}, 

where :math:`L` is the lower Cholesky factor of :math:`A` (the upper Cholesky
factor is computed in the case of upper-triangular storage). Rather than 
performing Cholesky factorization, triangular inversion, and then the Hermitian
triangular outer product in sequence, this routine makes use of the single-sweep
algorithm described in Bientinesi et al.'s "Families of algorithms related to 
the inversion of a symmetric positive definite matrix", in particular, the 
variant 2 algorithm from Fig. 9. 

If the matrix is found to not be HPD, then a :cpp:type:`NonHPDMatrixException`
is thrown.

.. cpp:function:: void HPDInverse( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void HPDInverse( UpperOrLower uplo, DistMatrix<F>& A )

   Overwrite the `uplo` triangle of the HPD matrix `A` with the same 
   triangle of the inverse of `A`.

Triangular inversion
--------------------
Inverts a (possibly unit-diagonal) triangular matrix in-place.

.. cpp:function:: void TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A )
.. cpp:function:: void TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A )

   Inverts the triangle of `A` specified by the parameter `uplo`; 
   if `diag` is set to `UNIT`, then `A` is treated as unit-diagonal.
