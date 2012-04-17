Norms
=====

Several matrix norm routines are provided for general, Hermitian, and symmetric 
(distributed) matrices; each of the following routines can return either
:math:`||A||_1`, :math:`||A||_\infty`, :math:`||A||_F` (the Frobenius norm), or 
the maximum entrywise norm. The matrix two-norm is quite expensive to directly 
compute, so a probabilistic algorithm (based upon Dixon's approach) will be 
added in the near future.

.. cpp:type:: NormType

   An enum that can be set to either

   * ``FROBENIUS_NORM``:

     .. math::

        ||A||_F = \sqrt{\sum_{i,j=0}^{n-1} |\alpha_{i,j}|^2}

   * ``INFINITY_NORM``: 

     .. math:: 
        :nowrap:

        \[
        ||A||_{\infty} = \max_{||x||_{\infty}=1} ||Ax||_{\infty} 
                       = \max_i \sum_{j=0}^{n-1} |\alpha_{i,j}|
        \]

   * ``ONE_NORM``: 
     
     .. math:: 
        :nowrap:

        \[
        ||A||_1 = \max_{||x||_1=1} ||Ax||_1 
                = \max_j \sum_{i=0}^{n-1} |\alpha_{i,j}|
        \]

   * ``MAX_NORM``: 
     
     .. math::
     
        ||A||_{\mbox{max}} = \max_{i,j} |\alpha_{i,j}|

Norm
----

For computing norms of fully-populated matrices.

.. cpp:function:: R Norm( const Matrix<R>& A, NormType type=FROBENIUS_NORM )

   Return the norm of the fully-populated real matrix `A`.

.. cpp:function:: R Norm( const DistMatrix<R,MC,MR>& A, NormType type=FROBENIUS_NORM )

   Return the norm of the fully-populated real distributed matrix `A`.

.. cpp:function:: R Norm( const Matrix<Complex<R> >& A, NormType type=FROBENIUS_NORM )

   Return the norm of the fully-populated complex matrix `A`.

.. cpp:function:: R Norm( const DistMatrix<Complex<R>,MC,MR>& A, NormType type=FROBENIUS_NORM )

   Return the norm of the fully-populated complex distributed matrix `A`.

HermitianNorm
-------------

Same as ``Norm``, but the (distributed) matrix is implicitly Hermitian 
with the data stored in the triangle specified by ``UpperOrLower``.

SymmetricNorm
-------------

Same as ``Norm``, but the (distributed) matrix is implicitly symmetric
with the data stored in the triangle specified by ``UpperOrLower``.

