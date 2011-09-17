Advanced linear algebra
***********************

Norms
=====

Norm
--------------
Describe :math:`||A||_1`, :math:`||A||_\infty`, and
:math:`||A||_F`, and the maximum norm.

HermitianNorm
-------------
Same, but with added Shape parameter. Give examples.

SymmetricNorm
-------------
Meant for symmetric matrices. Currenty a wrapper around ``HermitianNorm``.

Factorizations
==============

Cholesky factorization
----------------------
Describe ``advanced::Cholesky`` here.

:math:`LDL^H` factorization
---------------------------
Describe ``advanced::LDLH`` here.

:math:`LDL^T` factorization
---------------------------
Describe ``advanced::LDLT`` here.

:math:`LU` factorization
------------------------
Describe ``advanced::LU`` here.

:math:`LQ` factorization
------------------------
Describe ``advanced::LQ`` here.

:math:`QR` factorization
------------------------
Describe ``advanced::QR`` here.

Linear solvers
==============

Cholesky solve
--------------
Describe ``advanced::CholeskySolve`` here.

Gaussian Elimination
--------------------
Describe ``advanced::GaussianElimination`` here.

Householder solve
-----------------
Solves a general overdetermined or underdetermined linear systems using 
a :math:`QR` or :math:`LQ` factorization, respectively. The routine is named
``advanced::HouseholderSolve``.

Direct inversion
================

HPD inversion
-------------
Describe ``advanced::HPDInverse`` here.

Triangular inversion
--------------------
Describe ``advanced::TriangularInverse`` here.

Reduction to Condensed Form
===========================

Hermitian to tridiagonal
------------------------
Describe ``advanced::HermitianTridiag`` here.

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
Describe :math:`Ax=\lambda x` and ``advanced::HermitianEig`` here.

Skew-Hermitian eigensolver
--------------------------
Describe :math:`Gx=\lambda x` and ``advanced::SkewHermitianEig`` here.

Hermitian generalized-definite eigensolvers
-------------------------------------------
Describe :math:`ABx=\lambda x`, :math:`BAx=\lambda x`, and 
:math:`Ax=\lambda Bx` cases and ``advanced::HermitianGenDefiniteEig``.

Non-Hermitian eigensolver
-------------------------
Not yet written and probably not in the near future.

SVD
---
Not yet written but planned. Note that an SVD of a Hermitian matrix can easily be computed from the eigenvalue decomposition.

Utilities
=========

Householder reflectors
----------------------
Describe major difference from LAPACK's conventions (i.e., we do not treat
the identity matrix as a Householder transform since it requires the 
:math:`u` in :math:`H-I^2uu'` to have norm zero rather than one). Describe 
``advanced::Reflector`` here.

Reduction of Hermitian generalized-definite EVPs
------------------------------------------------
Describe the reduction steps of :math:`ABx=\lambda x`, :math:`BAx=\lambda x`, 
and :math:`Ax=\lambda Bx` using the operations :math:`A := L^H A L` and 
:math:`A := L^{-1} A L^{-H}`.

Applying packed Householder transforms
--------------------------------------
Describe ``advanced::ApplyPackedReflectors`` here.

Environment routines
====================
