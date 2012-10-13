Norms
=====

Several matrix norm routines are provided for general, Hermitian, and symmetric 
(distributed) matrices; each of the following routines can return either
:math:`\|A\|_1`, :math:`\|A\|_\infty`, :math:`\|A\|_F` (the Frobenius norm),
the maximum entrywise norm, :math:`\|A\|_2`, or :math:`\|A\|_*` 
(the nuclear/trace norm). 

Norm
----

For computing norms of fully-populated matrices.

.. cpp:function:: typename Base<F>::type Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM )

HermitianNorm
-------------

Same as :cpp:func:`Norm`, but the (distributed) matrix is implicitly 
Hermitian with the data stored in the triangle specified by 
:cpp:type:`UpperOrLower`. 
Also, while :cpp:func:`Norm` supports every type of distribution, 
:cpp:func:`HermitianNorm` currently only supports the standard matrix 
distribution.

.. cpp:function:: typename Base<F>::type HermitianNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type HermitianNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )

SymmetricNorm
-------------

Same as :cpp:func:`Norm`, but the (distributed) matrix is implicitly 
symmetric with the data stored in the triangle specified by 
:cpp:type:`UpperOrLower`. 
Also, while :cpp:func:`Norm` supports every type of distribution, 
:cpp:func:`SymmetricNorm` currently only supports the standard matrix 
distribution.

.. cpp:function:: typename Base<F>::type SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type SymmetricNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )

Two-norm estimates
------------------

Since the two-norm is extremely useful, but expensive to compute, it is useful
to be able to compute rough lower and upper bounds for it. The following 
routines provide cheap, rough estimates. The ability to compute sharper 
estimates will likely be added later.

.. cpp:function:: typename Base<F>::type TwoNormLowerBound( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type TwoNormLowerBound( const DistMatrix<F>& A )

   Return the tightest lower bound on :math:`\|A\|_2` implied by the following inequalities:

   .. math::

      \|A\|_2 \ge \|A\|_{\mathrm{max}},

   .. math::

      \|A\|_2 \ge \frac{1}{\sqrt{n}} \|A\|_{\infty},

   .. math::

      \|A\|_2 \ge \frac{1}{\sqrt{m}} \|A\|_1,\;\;\mathrm{and}

   .. math::

      \|A\|_2 \ge \frac{1}{\mathrm{min}(m,n)} \|A\|_F.

.. cpp:function:: typename Base<F>::type TwoNormUpperBound( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type TwoNormUpperBound( const DistMatrix<F>& A )

   Return the tightest upper bound on :math:`\|A\|_2` implied by the following inequalities:

   .. math::

      \|A\|_2 \le \sqrt{m n} \|A\|_{\mathrm{max}},

   .. math::

      \|A\|_2 \le \sqrt{m} \|A\|_{\infty},

   .. math::

      \|A\|_2 \le \sqrt{n} \|A\|_1,\;\;\mathrm{and}

   .. math::

      \|A\|_2 \le \sqrt{ \|A\|_1 \|A\|_{\infty} }.
