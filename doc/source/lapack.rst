=====
Norms
=====

--------
One norm
--------
Describe :math:`||A||_1`, ``elemental::lapack::OneNorm``, and 
``elemental::lapack::HermitianOneNorm`` here.

-------------
Infinity norm
-------------
Describe :math:`||A||_\infty`, ``elemental::lapack::InfinityNorm``, and 
``elemental::lapack::HermitianInfinityNorm`` here.

--------------
Frobenius norm
--------------
Describe :math:`||A||_F`, ``elemental::lapack::FrobeniusNorm``, and 
``elemental::lapack::HermitianFrobeniusNorm`` here.

==============
Factorizations
==============

----------------------
Cholesky factorization
----------------------
Describe ``elemental::lapack::Chol`` here.

----------------
LU factorization
----------------
Describe ``elemental::lapack::LU`` here.

----------------
QR factorization
----------------
Describe ``elemental::lapack::QR`` here.

==============
Linear Solvers
==============

--------------------
Gaussian Elimination
--------------------
Describe ``elemental::lapack::GaussElim`` here.

---------
QR Solver
---------
Not yet written but trivial.

===========================
Reduction to Condensed Form
===========================

------------------------------
Householder tridiagonalization
------------------------------
Describe ``elemental::lapack::Tridiag`` here.

-----------------------------
Householder bidiagonalization
-----------------------------
Not yet written.

----------------------------------------
Householder reduction to Hessenberg form
----------------------------------------
Not yet written.

====================
Eigensolvers and SVD
====================

---------------------
Hermitian eigensolver
---------------------
Describe :math:`Ax=\lambda x` and ``elemental::lapack::HermitianEig`` here.

--------------------------
Skew-Hermitian eigensolver
--------------------------
Describe :math:`Gx=\lambda x` and ``elemental::lapack::SkewHermitianEig`` here.

-------------------------------------------
Generalized Hermitian-definite eigensolvers
-------------------------------------------
Describe :math:`ABx=\lambda x`, :math:`BAx=\lambda x`, and 
:math:`Ax=\lambda Bx` cases and ``elemental::lapack::GeneralizedHermitianEig``.

-------------------------
Non-Hermitian eigensolver
-------------------------
Not yet written and probably not in the near future.

---
SVD
---
Not yet written but planned.

=========
Utilities
=========

----------------------
Householder reflectors
----------------------
Describe major difference from LAPACK's conventions (i.e., we do not treat
the identity matrix as a Householder transform since it requires the 
:math:`u` in :math:`H=I-2uu'` to have norm zero rather than one). Describe 
``elemental::lapack::Reflector`` here.

---------------------------------------------------
Reduction of the Generalized Hermitian-Definite EVP
---------------------------------------------------
Describe the reduction steps of :math:`ABx=\lambda x`, :math:`BAx=\lambda x`, 
and :math:`Ax=\lambda Bx` using the operations :math:`A := L^H A L` and 
:math:`A := L^{-1} A L^{-H}`.

-------------
UT transforms
-------------
Describe ``elemental::lapack::UT`` here.

