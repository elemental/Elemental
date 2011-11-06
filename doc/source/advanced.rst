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
compute, so a probabilistic algorithm (based upon Dixon's approach) will be 
added in the near future.

.. cpp:type:: enum NormType

   Can be set to either ``FROBENIUS_NORM``, ``INFINITY_NORM``, ``MAX_NORM``, 
   or ``ONE_NORM``.

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

Invariants
==========

Determinant
-----------
Though there are many different possible definitions of the determinant of a 
matrix :math:`A \in \mathbb{F}^{n \times n}`, the simplest one is in terms of 
the product of the eigenvalues (including multiplicity):

.. math::
   :nowrap:

   \[
   \mbox{det}(A) = \prod_{i=0}^{n-1} \lambda_i.
   \]

Since :math:`\mbox{det}(AB)=\mbox{det}(A)\mbox{det}(B)`, we can compute the 
determinant of an arbitrary matrix in :math:`\mathcal{O}(n^3)` work by 
computing its LU decomposition (with partial pivoting), :math:`PA=LU`, 
recognizing that :math:`\mbox{det}(P)=\pm 1` 
(the *signature* of the permutation), and computing

.. math::
   :nowrap:

   \[
   \mbox{det}(A) = \mbox{det}(P)\mbox{det}(L)\mbox{det}(U) 
                 = \mbox{det}(P) \prod_{i=0}^{n-1} \upsilon_{i,i}
                 = \pm \prod_{i=0}^{n-1} \upsilon_{i,i},
   \]

where :math:`\upsilon_{i,i}` is the i'th diagonal entry of :math:`U`.

.. note:: 

   The following functions overwrite the input matrix with its LU factorization
   in order to efficiently compute the determinant.

.. cpp:function:: F advanced::Determinant( Matrix<F>& A )

   Returns the determinant of the square matrix ``A``, which is overwritten 
   during the computation.

.. cpp:function:: F advanced::Determinant( DistMatrix<F,MC,MR>& A )

   Returns the determinant of the square distributed matrix ``A``, which is 
   overwritten during the computation.

.. cpp:type:: struct SafeProduct<F>

   Represents the product of ``n`` values as :math:`\rho \exp(\kappa n)`, 
   where :math:`|\rho| \le 1` and :math:`\kappa \in \mathbb{R}`.

   .. cpp:member:: F rho

      For nonzero values, ``rho`` is the modulus and lies *on* the unit 
      circle; in order to represent a value that is precisely zero, ``rho`` 
      is set to zero.

   .. cpp:member:: typename RealBase<F>::type kappa

      ``kappa`` can be an arbitrary real number.

   .. cpp:member:: int n

      The number of values in the product.

.. cpp:function:: SafeProduct<F> advanced::SafeDeterminant( Matrix<F>& A )

   Returns the determinant of the square matrix ``A`` in an expanded form 
   which is less likely to over/under-flow.

.. cpp:function:: SafeProduct<F> advanced::SafeDeterminant( DistMatrix<F,MC,MR>& A )

   Returns the determinant of the square distributed matrix ``A`` in an 
   expanded form which is less likely to over/under-flow.

Trace
-----
The two equally useful definitions of the trace of a square matrix 
:math:`A \in \mathbb{F}^{n \times n}` are

.. math::
   :nowrap:

   \[
   \mbox{tr}(A) = \sum_{i=0}^{n-1} \alpha_{i,i} = \sum_{i=0}^{n-1} \lambda_i,
   \]

where :math:`\alpha_{i,i}` is the i'th diagonal entry of :math:`A` and 
:math:`\lambda_i` is the i'th eigenvalue (counting multiplicity) of :math:`A`.

Clearly the former equation is easier to compute, but the latter is an 
important characterization.

.. cpp:function:: F advanced::Trace( const Matrix<F>& A )

   Return the trace of the square matrix ``A``.

.. cpp:function:: F advanced::Trace( const DistMatrix<F,MC,MR>& A )

   Return the trace of the square distributed matrix ``A``.

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

   .. warning::

      The following routines do not pivot, so please use with caution.

.. cpp:function:: void advanced::LDLH( Matrix<F>& A, Matrix<F>& d )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`, and then also return the diagonal of :math:`D` 
   in the vector ``d``. 

.. cpp:function:: void advanced::LDLH( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but for distributed matrices.

:math:`LDL^T` factorization
---------------------------
While the :math:`LDL^H` factorization targets Hermitian matrices, the 
:math:`LDL^T` factorization targets symmetric matrices.

   .. warning::

      The following routines do not pivot, so please use with caution.

.. cpp:function:: void advanced::LDLT( Matrix<F>& A, Matrix<F>& d )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`, and then also return the diagonal of :math:`D` 
   in the vector ``d``. 

.. cpp:function:: void advanced::LDLT( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but for distributed matrices.

:math:`LU` factorization
------------------------
Given :math:`A \in \mathbb{F}^{m \times n}`, an LU factorization 
(without pivoting) finds a unit lower-trapezoidal 
:math:`L \in \mathbb{F}^{m \times \mbox{min}(m,n)}` and upper-trapezoidal 
:math:`U \in \mathbb{F}^{\mbox{min}(m,n) \times n}` such that :math:`A=LU`. 
Since :math:`L` is required to have its diaganal entries set to one: the upper 
portion of :math:`A` can be overwritten with `U`, and the strictly lower 
portion of :math:`A` can be overwritten with the strictly lower portion of 
:math:`L`.

.. cpp:function:: void advanced::LU( Matrix<F>& A )

   Overwrites :math:`A` with its LU decomposition.

.. cpp:function:: void advanced::LU( DistMatrix<F,MC,MR>& A )

   Overwrites :math:`A` with its LU decomposition.

Since LU factorization without pivoting is known to be unstable for general 
matrices, it is standard practice to pivot the rows of :math:`A` during the 
factorization (this is called partial pivoting since the columns are not also 
pivoted). An LU factorization with partial pivoting therefore computes 
:math:`P`, :math:`L`, and :math:`U` such that :math:`PA=LU`, where :math:`L` 
and :math:`U` are as described above and :math:`P` is a permutation matrix.

.. cpp:function:: void advanced::LU( Matrix<F>& A, Matrix<int>& p )

   Ovewrites :math:`A` with the LU decomposition of :math:`PA`, where 
   :math:`P` is represented by the pivot vector ``p``.

.. cpp:function:: void advanced::LU( DistMatrix<F,MC,MR>& A, DistMatrix<F,VC,STAR>& p )

   Overwrites the distributed matrix :math:`A` with the LU decomposition of 
   :math:`PA`, where :math:`P` is represented by the pivot vector ``p``.

:math:`LQ` factorization
------------------------
Given :math:`A \in \mathbb{F}^{m \times n}`, an LQ factorization typically 
computes an implicit unitary matrix :math:`\hat Q \in \mathbb{F}^{n \times n}` 
such that :math:`\hat L \equiv A\hat Q^H` is lower trapezoidal. One can then 
form the thin factors :math:`L \in \mathbb{F}^{m \times \mbox{min}(m,n)}` and 
:math:`Q \in \mathbb{F}^{\mbox{min}(m,n) \times n}` by setting 
:math:`L` and :math:`Q` to first :math:`\mbox{min}(m,n)` columns and rows of 
:math:`\hat L` and :math:`\hat Q`, respectively. Upon completion :math:`L` is 
stored in the lower trapezoid of :math:`A` and the Householder reflectors 
representing :math:`\hat Q` are stored within the rows of the strictly upper 
trapezoid.

.. cpp:function:: void advanced::LQ( DistMatrix<R,MC,MR>& A )

   Overwrite the real distributed matrix :math:`A` with :math:`L` and the 
   Householder reflectors representing :math:`\hat Q`.

.. cpp:function:: void advanced::LQ( DistMatrix<std::complex<R>,MC,MR>& A, DistMatrix<std::complex<R>,MD,STAR>& t )

   Overwrite the complex distributed matrix :math:`A` with :math:`L` and the 
   Householder reflectors representing :math:`\hat Q`; unlike the real case, 
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector ``t``.

:math:`QR` factorization
------------------------
Given :math:`A \in \mathbb{F}^{m \times n}`, a QR factorization typically 
computes an implicit unitary matrix :math:`\hat Q \in \mathbb{F}^{m \times m}` 
such that :math:`\hat R \equiv \hat Q^H A` is upper trapezoidal. One can then 
form the thin factors :math:`Q \in \mathbb{F}^{m \times \mbox{min}(m,n)}` and
:math:`R \in \mathbb{F}^{\mbox{min}(m,n) \times n}` by setting 
:math:`Q` and :math:`R` to first :math:`\mbox{min}(m,n)` columns and rows of 
:math:`\hat Q` and :math:`\hat R`, respectively. Upon completion :math:`R` is 
stored in the upper trapezoid of :math:`A` and the Householder reflectors 
representing :math:`\hat Q` are stored within the columns of the strictly lower 
trapezoid.

.. cpp:function:: void advanced::QR( DistMatrix<R,MC,MR>& A )

   Overwrite the real distributed matrix :math:`A` with :math:`R` and the 
   Householder reflectors representing :math:`\hat Q`.

.. cpp:function:: void advanced::QR( DistMatrix<std::complex<R>,MC,MR>& A, DistMatrix<std::complex<R>,MD,STAR>& t )

   Overwrite the complex distributed matrix :math:`A` with :math:`R` and the 
   Householder reflectors representing :math:`\hat Q`; unlike the real case, 
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector ``t``.

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

**TODO:** Describe ``advanced::SortEig`` here.

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
:math:`u` in :math:`H=I-2uu'` to have norm zero rather than one). 

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
