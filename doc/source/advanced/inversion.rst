In-place inversion
==================

General inversion
-----------------
This routine computes the in-place inverse of a general fully-populated 
(invertible) matrix :math:`A` as

.. math::
   :nowrap:

   \[
   A^{-1} = U^{-1} L^{-1} P,
   \]

where :math:`PA=LU` is the result of LU factorization with partial pivoting.
The algorithm essentially factors :math:`A`, inverts :math:`U` in place, 
solves against :math:`L` one block column at a time, and then applies the 
row pivots in reverse order to the columns of the result.

.. cpp:function:: void Inverse( Matrix<F>& A )

   Overwrites the general matrix `A` with its inverse.

.. cpp:function:: void Inverse( DistMatrix<F,MC,MR>& A )

   The same as above, but for distributed matrices.

Hermitian pseudoinverse
-----------------------
Computes the pseudoinverse of a Hermitian matrix through a customized version of
``RealHermitianFunction`` which used the eigenvalue mapping function

.. math::
   :nowrap:

   \[
   f(\omega_i) = \left\{\begin{array}{cc} 
     1/\omega_i, & |\omega_i| \ge \epsilon \, n \, ||A||_2 \\
         0,      & \mbox{otherwise}
   \end{array}\right.,
   \]

where :math:`\epsilon` is the relative machine precision,
:math:`n` is the height of :math:`A`, and :math:`||A||_2` can be computed
as the maximum absolute value of the eigenvalues of :math:`A`.

.. cpp:function:: HermitianPseudoinverse( UpperOrLower uplo, DistMatrix<R,MC,MR>& A )

   Computes the pseudoinverse of a distributed real symmetric matrix with data
   stored in the `uplo` triangle.

.. cpp:function:: HermitianPseudoinverse( UpperOrLower uplo, DistMatrix<std::complex<R>,MC,MR>& A )

   Computes the pseudoinverse of a distributed complex Hermitian matrix with
   data stored in the `uplo` triangle.

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

.. cpp:function:: void HPDInverse( UpperOrLower uplo, Matrix<F>& A )

   Overwrite the `uplo` triangle of the HPD matrix `A` with the same 
   triangle of the inverse of `A`.

.. cpp:function:: void HPDInverse( UpperOrLower uplo, DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.


Triangular inversion
--------------------
Inverts a (possibly unit-diagonal) triangular matrix in-place.

.. cpp:function:: void TriangularInverse( UpperOrLower uplo, Diagonal diagonal, Matrix<F>& A )

   Inverts the triangle of `A` specified by the parameter `uplo`; 
   if `diagonal` is set to `UNIT`, then `A` is treated as unit-diagonal.

.. cpp:function:: void TriangularInverse( UpperOrLower uplo, Diagonal diagonal, DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.
