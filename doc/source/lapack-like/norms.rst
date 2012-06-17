Norms
=====

Several matrix norm routines are provided for general, Hermitian, and symmetric 
(distributed) matrices; each of the following routines can return either
:math:`\|A\|_1`, :math:`\|A\|_\infty`, :math:`\|A\|_F` (the Frobenius norm), or 
the maximum entrywise norm. The matrix two-norm is quite expensive to directly 
compute, so a probabilistic algorithm (based upon Dixon's approach) will be 
added in the near future.

.. cpp:type:: NormType

   An enum that can be set to either

   * ``FROBENIUS_NORM``:

     .. math::

        \|A\|_F = \sqrt{\sum_{i=0}^{m-1} \sum_{j=0}^{n-1} |\alpha_{i,j}|^2}

   * ``INFINITY_NORM``: 

     .. math:: 
        :nowrap:

        \[
        \|A\|_{\infty} = \max_{\|x\|_{\infty}=1} \|Ax\|_{\infty} 
                       = \max_i \sum_{j=0}^{n-1} |\alpha_{i,j}|
        \]

   * ``ONE_NORM``: 
     
     .. math:: 
        :nowrap:

        \[
        \|A\|_1 = \max_{\|x\|_1=1} \|Ax\|_1 
                = \max_j \sum_{i=0}^{m-1} |\alpha_{i,j}|
        \]

   * ``MAX_NORM``: 
     
     .. math::
     
        \|A\|_{\mbox{max}} = \max_{i,j} |\alpha_{i,j}|

Norm
----

For computing norms of fully-populated matrices.

.. cpp:function:: typename Base<F>::type Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM )

HermitianNorm
-------------

Same as ``Norm``, but the (distributed) matrix is implicitly Hermitian 
with the data stored in the triangle specified by ``UpperOrLower``. Also, 
while ``Norm`` supports every type of distribution, ``HermitianNorm`` currently
only supports the standard matrix distribution.

SymmetricNorm
-------------

Same as ``Norm``, but the (distributed) matrix is implicitly symmetric
with the data stored in the triangle specified by ``UpperOrLower``.

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
