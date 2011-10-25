Advanced linear algebra
***********************

This chapter describes all of the linear algebra operations which are not basic 
enough to fall within the BLAS (Basic Linear Algebra Subprograms) framework, and 
are therefore in the domain of LAPACK-like software.

Norms
=====

Several matrix norm routines are provided for general, Hermitian, and symmetric 
(distributed) matrices; each of the following routines can return either
:math:`||A||_1`, :math:`||A||_\infty`, :math:`||A||_F` (the Frobenius norm), or 
the maximum entrywise norm. The matrix two-norm is quite expensive to directly 
compute, so a probabilistic algorithm (based upon Dixon's approach) will be added 
in the near future.

Norm
----

.. cpp:function:: R advanced::Norm( const Matrix<R>& A, NormType type=FROBENIUS_NORM )

   Return the norm of the real matrix ``A``.

.. cpp:function:: R advanced::Norm( const DistMatrix<R,MC,MR>& A, NormType type=FROBENIUS_NORM )

   Return the norm of the real distributed matrix ``A``.

.. cpp:function:: R advanced::Norm( const Matrix<std::complex<R> >& A, NormType type=FROBENIUS_NORM )

   Return the norm of the complex matrix ``A``.

.. cpp:function:: R advanced::Norm( const DistMatrix<std::complex<R>,MC,MR>& A, NormType type=FROBENIUS_NORM )

   Return the norm of the complex distributed matrix ``A``.

HermitianNorm
-------------
Same as ``advanced::Norm``, but the (distributed) matrix is implicitly Hermitian 
with the data stored in the triangle specified by ``shape``.

SymmetricNorm
-------------
Same as ``advanced::Norm``, but the (distributed) matrix is implicitly symmetric
with the data stored in the triangle specified by ``shape``.

Factorizations
==============

Cholesky factorization
----------------------
It is well-known that Hermitian positive-definite (HPD) matrices can be decomposed
into the form :math:`A = L L^H` or :math:`A = U^H U`, where :math:`L=U^H` is lower
triangular, and Cholesky factorization provides such an :math:`L` (or :math:`U`) 
given an HPD :math:`A`.

.. cpp:function:: void advanced::Cholesky( Shape shape, Matrix<F>& A )

   Overwrite the ``shape`` triangle of the HPD matrix `A` with its Cholesky factor.

.. cpp:function:: void advanced::Cholesky( Shape shape, DistMatrix<F,MC,MR>& A )

   Overwrite the ``shape`` triangle of the distributed HPD matrix ``A`` with its 
   Cholesky factor.

:math:`LDL^H` factorization
---------------------------
Though the Cholesky factorization is ideal for most HPD matrices, there exist 
many Hermitian matrices whose eigenvalues are not all positive. The 
:math:`LDL^H` factorization exists as slight relaxation of the Cholesky 
factorization, i.e., it computes lower-triangular (with unit diagonal) :math:`L` 
and diagonal :math:`D` such that :math:`A = L D L^H`.

.. cpp:function:: void advanced::LDLH( Matrix<F>& A, Matrix<F>& d )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`, and then also return the diagonal of :math:`D` 
   in the vector ``d``. 

   .. warning::

      No pivoting is currently performed, so please use with caution.

.. cpp:function:: void advanced::LDLH( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but for distributed matrices.

:math:`LDL^T` factorization
---------------------------
While the :math:`LDL^H` factorization targets Hermitian matrices, the 
:math:`LDL^T` factorization targets symmetric matrices.

.. cpp:function:: void advanced::LDLT( Matrix<F>& A, Matrix<F>& d )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`, and then also return the diagonal of :math:`D` 
   in the vector ``d``. 

   .. warning::
      
      No pivoting is currently performed, so please use with caution.

.. cpp:function:: void advanced::LDLT( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but for distributed matrices.

:math:`LU` factorization
------------------------
**TODO:** Describe ``advanced::LU`` here.

:math:`LQ` factorization
------------------------
**TODO:** Describe ``advanced::LQ`` here.

:math:`QR` factorization
------------------------
**TODO:** Describe ``advanced::QR`` here.

Linear solvers
==============

Cholesky solve
--------------
**TODO:** Describe ``advanced::CholeskySolve`` here.

Gaussian Elimination
--------------------
**TODO:** Describe ``advanced::GaussianElimination`` here.

Householder solve
-----------------
**TODO:** Describe ``advanced::HouseholderSolve``. 
Solves a general overdetermined or underdetermined linear systems using 
a :math:`QR` or :math:`LQ` factorization, respectively. 

Direct inversion
================

HPD inversion
-------------
**TODO:** Describe ``advanced::HPDInverse`` here.

Triangular inversion
--------------------
**TODO:** Describe ``advanced::TriangularInverse`` here.

Reduction to Condensed Form
===========================

Hermitian to tridiagonal
------------------------
**TODO:** Describe ``advanced::HermitianTridiag`` here.

General to Hessenberg
---------------------
Not yet written.

General to bidiagonal
---------------------
Not yet written.

Eigensolvers and SVD
====================

Hermitian eigensolver
---------------------
**TODO:** Describe :math:`Ax=\lambda x` and ``advanced::HermitianEig`` here.

Skew-Hermitian eigensolver
--------------------------
**TODO:** Describe :math:`Gx=\lambda x` and ``advanced::SkewHermitianEig`` here.

Hermitian generalized-definite eigensolvers
-------------------------------------------
**TODO:** Describe :math:`ABx=\lambda x`, :math:`BAx=\lambda x`, and 
:math:`Ax=\lambda Bx` cases and ``advanced::HermitianGenDefiniteEig``.

Unitary eigensolver
-------------------
Not yet written.

Normal eigensolver
------------------
Not yet written and probably not in the near future.

General eigensolver
-------------------
Not yet written and probably not in the near future.

Schur decomposition
-------------------
Not yet written and probably not in the near future.

Hermitian SVD
-------------
Not yet written, but planned. Note that an SVD of a Hermitian matrix can easily be computed from the eigenvalue decomposition.

General SVD
-----------
Not yet written, but planned. 

Utilities
=========

Householder reflectors
----------------------
**TODO:** Describe major difference from LAPACK's conventions (i.e., we do not 
treat the identity matrix as a Householder transform since it requires the 
:math:`u` in :math:`H-I^2uu'` to have norm zero rather than one). 

Reduction of Hermitian generalized-definite EVPs
------------------------------------------------
**TODO:** Describe the reduction steps of :math:`ABx=\lambda x`, 
:math:`BAx=\lambda x`, and :math:`Ax=\lambda Bx` using the operations 
:math:`A := L^H A L` and :math:`A := L^{-1} A L^{-H}`.

Applying packed Householder transforms
--------------------------------------
**TODO:** Describe ``advanced::ApplyPackedReflectors`` here.

Tuning parameters
=================
**TODO:** Describe the tuning parameters.
