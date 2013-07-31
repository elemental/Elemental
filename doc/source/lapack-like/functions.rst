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
Hermitian matrix :math:`A` is :math:`A=Z \Lambda Z^H`, then

.. math::

   f(A) = f(Z \Lambda Z^H) = Z f(\Lambda) Z^H.

Two important special cases are :math:`f(\lambda) = \exp(\lambda)` and 
:math:`f(\lambda)=\exp(i \lambda)`, where the former results in a Hermitian 
matrix and the latter in a normal (in fact, unitary) matrix.

.. note:: 

   Since Elemental currently depends on PMRRR for its tridiagonal 
   eigensolver, only double-precision results are supported as of now.

.. cpp:function:: void RealHermitianFunction( UpperOrLower uplo, Matrix<F>& A, const RealFunctor& f )
.. cpp:function:: void RealHermitianFunction( UpperOrLower uplo, DistMatrix<F>& A, const RealFunctor& f )

   Modifies the eigenvalues of the passed-in Hermitian matrix by replacing 
   each eigenvalue :math:`\lambda_i` with :math:`f(\lambda_i) \in \mathbb{R}`. 
   ``RealFunctor`` is any 
   class which has the member function ``R operator()( R omega ) const``.
   See `examples/lapack-like/RealSymmetricFunction.cpp <https://github.com/elemental/Elemental/tree/master/examples/lapack-like/RealHermitianFunction.cpp>`_ for an example usage.

.. cpp:function:: void ComplexHermitianFunction( UpperOrLower uplo, Matrix<Complex<R> >& A, const ComplexFunctor& f )
.. cpp:function:: void ComplexHermitianFunction( UpperOrLower uplo, DistMatrix<Complex<R> >& A, const ComplexFunctor& f )

   Modifies the eigenvalues of the passed-in complex Hermitian matrix by
   replacing each eigenvalue :math:`\lambda_i` with 
   :math:`f(\lambda_i) \in \mathbb{C}`. ``ComplexFunctor`` can be any class
   which has the member function ``Complex<R> operator()( R omega ) const``.
   See `examples/lapack-like/ComplexHermitianFunction.cpp <https://github.com/elemental/Elemental/tree/master/examples/lapack-like/ComplexHermitianFunction.cpp>`_ for an example usage.

**TODO: A version of ComplexHermitianFunction which begins with a real matrix**

Pseudoinverse
-------------

.. cpp:function:: Pseudoinverse( Matrix<F>& A, typename Base<F>::type tolerance=0 )
.. cpp:function:: Pseudoinverse( DistMatrix<F>& A, typename Base<F>::type tolerance=0 )

   Computes the pseudoinverse of a general matrix through computing its SVD,
   modifying the singular values with the function

   .. math::

      f(\sigma) = \left\{\begin{array}{cc} 
        1/\sigma, & \sigma \ge \epsilon \, n \, \| A \|_2 \\
            0,      & \mbox{otherwise}
      \end{array}\right.,

   where :math:`\epsilon` is the relative machine precision,
   :math:`n` is the height of :math:`A`, and :math:`\| A \|_2` is the maximum
   singular value.
   If a nonzero value for `tolerance` was specified, it is used instead of 
   :math:`\epsilon n \| A \|_2`.

.. cpp:function:: HermitianPseudoinverse( UpperOrLower uplo, Matrix<F>& A, typename Base<F>::type tolerance=0 )
.. cpp:function:: HermitianPseudoinverse( UpperOrLower uplo, DistMatrix<F>& A, typename Base<F>::type tolerance=0 )

   Computes the pseudoinverse of a Hermitian matrix through a customized version
   of :cpp:func:`RealHermitianFunction` which used the eigenvalue mapping 
   function

   .. math::

      f(\omega) = \left\{\begin{array}{cc} 
        1/\omega, & |\omega| \ge \epsilon \, n \, \| A \|_2 \\
            0,      & \mbox{otherwise}
      \end{array}\right.,

   where :math:`\epsilon` is the relative machine precision,
   :math:`n` is the height of :math:`A`, and :math:`\| A \|_2` can be computed
   as the maximum absolute value of the eigenvalues of :math:`A`.
   If a nonzero value for `tolerance` is specified, it is used instead of
   :math:`\epsilon n \| A \|_2`.

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

.. cpp:function:: void SquareRoot( Matrix<F>& A )
.. cpp:function:: void SquareRoot( DistMatrix<F>& A )

   Currently uses a Newton iteration to compute the general matrix square-root. 
   See ``square_root::Newton`` for the more detailed interface.

.. cpp:function:: void HPSDSquareRoot( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void HPSDSquareRoot( UpperOrLower uplo, DistMatrix<F>& A )

   Computes the Hermitian EVD, square-roots the eigenvalues, and then 
   reforms the matrix. If any of the eigenvalues were sufficiently negative,
   a :cpp:type:`NonHPSDMatrixException` is thrown.

**TODO: HermitianSquareRoot**

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:function:: int square_root::Newton( Matrix<F>& A, int maxIts=100, typename Base<F>::tol=0 )
.. cpp:function:: int square_root::Newton( DistMatrix<F>& A, int maxIts=100, typename Base<F>::tol=0 )

   Performs at most ``maxIts`` Newton steps in an attempt to compute the 
   matrix square-root within the specified tolerance, which defaults to 
   :math:`n \epsilon`, where :math:`n` is the matrix height and :math:`\epsilon`
   is the machine precision.

Sign
----
The matrix sign function can be written as

.. math::
   \text{sgn}(A) = A(A^2)^{-1/2},

as long as :math:`A` does not have any pure-imaginary eigenvalues.

.. cpp:function:: void Sign( Matrix<F>& A )
.. cpp:function:: void Sign( DistMatrix<F>& A )
.. cpp:function:: void Sign( Matrix<F>& A, Matrix<F>& N )
.. cpp:function:: void Sign( DistMatrix<F>& A, DistMatrix<F>& N )

   Compute the matrix sign through a globally-convergent Newton iteration
   scaled with the Frobenius norm of the iterate and its inverse.
   Optionally return the full decomposition, :math:`A=S N`, where :math:`A`
   is overwritten by :math:`S`.

.. cpp:function:: void HermitianSign( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void HermitianSign( UpperOrLower uplo, DistMatrix<F>& A )
.. cpp:function:: void HermitianSign( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N )
.. cpp:function:: void HermitianSign( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& N )

   Compute the Hermitian EVD, replace the eigenvalues with their sign, and then
   reform the matrix. Optionally return the full decomposition, :math:`A=SN`,
   where :math:`A` is overwritten by :math:`S`. Note that this will also be 
   a polar decomposition.

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:type:: sign::Scaling

   An enum for specifying the scaling strategy to be used for the Newton 
   iteration for the matrix sign function. It must be either ``NONE``, 
   ``DETERMINANT``, or ``FROB_NORM`` (the default).

.. cpp:function:: int sign::Newton( Matrix<F>& A, sign::Scaling scaling=FROB_NORM, int maxIts=100, typename Base<F>::type tol=0 )
.. cpp:function:: int sign::Newton( DistMatrix<F>& A, sign::Scaling scaling=FROB_NORM, int maxIts=100, typename Base<F>::type tol=0 )

   Runs a (scaled) Newton iteration for at most ``maxIts`` iterations with 
   the specified tolerance, which, if undefined, is set to :math:`n \epsilon`,
   where :math:`n` is the matrix dimension and :math:`\epsilon` is the 
   machine epsilon. The return value is the number of performed iterations.
