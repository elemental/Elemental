Matrix functions
================

Hermitian functions
-------------------
Reform the matrix with the eigenvalues modified by a user-defined function. 
When the user-defined function is real-valued, the result will remain Hermitian,
but when the function is complex-valued, the result is best characterized as 
normal. 

When the user-defined function, say :math:`f`, is analytic, we can say much
more about the result: if the eigenvalue decomposition of the 
Hermitian matrix :math:`A` is :math:`A=Z \Omega Z^H`, then

.. math::

   f(A) = f(Z \Omega Z^H) = Z f(\Omega) Z^H.

Two important special cases are :math:`f(\lambda) = \exp(\lambda)` and 
:math:`f(\lambda)=\exp(i \lambda)`, where the former results in a Hermitian 
matrix and the latter in a normal (in fact, unitary) matrix.

.. note:: 

   Since Elemental currently depends on PMRRR for its tridiagonal 
   eigensolver, only double-precision results are supported as of now.

.. cpp:function:: void RealHermitianFunction( UpperOrLower uplo, DistMatrix<F>& A, const RealFunctor& f )

   Modifies the eigenvalues of the passed-in Hermitian matrix by replacing 
   each eigenvalue :math:`\omega_i` with :math:`f(\omega_i) \in \mathbb{R}`. 
   ``RealFunctor`` is any 
   class which has the member function ``R operator()( R omega ) const``.
   See `examples/lapack-like/RealSymmetricFunction.cpp <../../../../examples/lapack-like/RealHermitianFunction.cpp>`_ for an example usage.

.. cpp:function:: void ComplexHermitianFunction( UpperOrLower uplo, DistMatrix<Complex<R> >& A, const ComplexFunctor& f )

   Modifies the eigenvalues of the passed-in complex Hermitian matrix by
   replacing each eigenvalue :math:`\omega_i` with 
   :math:`f(\omega_i) \in \mathbb{C}`. ``ComplexFunctor`` can be any class
   which has the member function ``Complex<R> operator()( R omega ) const``.
   See `examples/lapack-like/ComplexHermitianFunction.cpp <../../../../examples/lapack-like/ComplexHermitianFunction.cpp>`_ for an example usage.

**TODO: A version of ComplexHermitianFunction which begins with a real matrix**

Pseudoinverse
-------------

.. cpp:function:: Pseudoinverse( DistMatrix<F>& A )

   Computes the pseudoinverse of a general matrix through computing its SVD,
   modifying the singular values with the function

   .. math::
      :nowrap:

      \[
      f(\sigma) = \left\{\begin{array}{cc} 
        1/\sigma, & \sigma \ge \epsilon \, n \, ||A||_2 \\
            0,      & \mbox{otherwise}
      \end{array}\right.,
      \]

   where :math:`\epsilon` is the relative machine precision,
   :math:`n` is the height of :math:`A`, and :math:`\|A\|_2` is the maximum
   singular value.

.. cpp:function:: HermitianPseudoinverse( UpperOrLower uplo, DistMatrix<F>& A )

   Computes the pseudoinverse of a Hermitian matrix through a customized version of ``RealHermitianFunction`` which used the eigenvalue mapping function

   .. math::
      :nowrap:

      \[
      f(\omega) = \left\{\begin{array}{cc} 
        1/\omega, & |\omega| \ge \epsilon \, n \, ||A||_2 \\
            0,      & \mbox{otherwise}
      \end{array}\right.,
      \]

   where :math:`\epsilon` is the relative machine precision,
   :math:`n` is the height of :math:`A`, and :math:`||A||_2` can be computed
   as the maximum absolute value of the eigenvalues of :math:`A`.

Square root
-----------
A matrix :math:`B` satisfying 

.. math::

   B^2 = A

is referred to as the *square-root* of the matrix :math:`A`. Such a matrix
is guaranteed to exist as long as :math:`A` is diagonalizable: if 
:math:`A = X \Lambda X^{-1}`, then we may put 

.. math::

   B = X \sqrt{\Lambda} X^{-1},

where each eigenvalue :math:`\lambda = r e^{i\theta}` maps to
:math:`\sqrt{\lambda} = \sqrt{r} e^{i\theta/2}`. 

.. cpp:function:: void HPSDSquareRoot( UpperOrLower uplo, DistMatrix<F>& A )

   Hermitian matrices with non-negative eigenvalues have a natural matrix square root which remains Hermitian. This routine attempts to overwrite a matrix with its square root and throws a ``NonHPSDMatrixException`` if any sufficiently negative eigenvalues are computed.

**TODO: HermitianSquareRoot**

Semi-definite Cholesky
----------------------

It is possible to compute the Cholesky factor of a Hermitian positive
semi-definite (HPSD) matrix through its eigenvalue decomposition, though it
is significantly more expensive than the HPD case: Let :math:`A = U \Lambda U^H`
be the eigenvalue decomposition of :math:`A`, where all entries of 
:math:`\Lambda` are non-negative. Then :math:`B = U \sqrt \Lambda U^H` is the 
matrix square root of :math:`A`, i.e., :math:`B B = A`, and it follows that the 
QR and LQ factorizations of :math:`B` yield Cholesky factors of :math:`A`:

.. math::
   A = B B = B^H B = (Q R)^H (Q R) = R^H Q^H Q R = R^H R,

and

.. math::
   A = B B = B B^H = (L Q) (L Q)^H = L Q Q^H L^H = L L^H.

If :math:`A` is found to have eigenvalues less than :math:`-n \epsilon ||A||_2`,
then a ``NonHPSDMatrixException`` will be thrown.

.. cpp:function:: void HPSDCholesky( UpperOrLower uplo, DistMatrix<F>& A )

   Overwrite the `uplo` triangle of the potentially singular matrix `A` with its Cholesky factor.

