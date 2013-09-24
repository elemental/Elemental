Eigensolvers and SVD
====================

Hermitian eigensolver
---------------------
Elemental provides a collection of routines for both full and partial 
solution of the Hermitian eigenvalue problem 

.. math::

   A Z = Z \Lambda,

where `A` is the given Hermitian matrix, and unitary `Z` and real diagonal 
:math:`\Lambda` are sought. In particular, with the eigenvalues and 
corresponding eigenpairs labeled in non-decreasing order, the three basic 
modes are:

1. Compute all eigenvalues or eigenpairs, :math:`\{\lambda_i\}_{i=0}^{n-1}` or 
   :math:`\{(z_i,\lambda_i)\}_{i=0}^{n-1}`.
2. Compute the eigenvalues or eigenpairs with a given range of indices, say  
   :math:`\{\lambda_i\}_{i=a}^b` or :math:`\{(z_i,\lambda_i)\}_{i=a}^b`, 
   with :math:`0 \le a \le b < n`.
3. Compute all eigenpairs (or just eigenvalues) with eigenvalues lying in a 
   particular half-open interval, either
   :math:`\{\lambda_i \;|\; \lambda_i \in (a,b] \}` or 
   :math:`\{ (z_i,\lambda_i) \;|\; \lambda_i \in (a,b] \}`.

As of now, all three approaches start with Householder tridiagonalization 
(ala :cpp:func:`HermitianTridiag`) and then call Matthias Petschow and 
Paolo Bientinesi's PMRRR for the tridiagonal eigenvalue problem.

.. note::

   Unfortunately, PMRRR currently only supports double-precision problems, and 
   so the parallel versions of these routines are limited to real and complex 
   double-precision matrices.

.. note:: 

   Please see the :ref:`lapack-tuning` section for information on optimizing
   the reduction to tridiagonal form, as it is the dominant cost in all of 
   Elemental's Hermitian eigensolvers.

Full spectrum computation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianEig( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& w, SortType sort=UNSORTED )
.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, SortType sort=UNSORTED )

   Compute the full set of eigenvalues of the Hermitian matrix `A`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& w, Matrix<F>& Z, SortType sort=UNSORTED )
.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<F>& Z, SortType sort=UNSORTED )

   Compute the full set of eigenpairs of the Hermitian matrix `A`.

Index-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianEig( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& w, int a, int b, SortType sort=UNSORTED )
.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, int a, int b, SortType sort=UNSORTED )

   Compute the eigenvalues of a Hermitian matrix `A` with indices in the range 
   :math:`a,a+1,...,b`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& w, Matrix<F>& Z, SortType sort=UNSORTED )
.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<F>& Z, SortType sort=UNSORTED )

   Compute the eigenpairs of a Hermitian matrix `A` with indices in the range 
   :math:`a,a+1,...,b`.

Range-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianEig( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& w, typename Base<F>::type a, typename Base<F>::type b, SortType sort=UNSORTED )
.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, typename Base<F>::type a, typename Base<F>::type b, SortType sort=UNSORTED )

   Compute the eigenvalues of a Hermitian matrix `A` lying in the half-open 
   interval :math:`(a,b]`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& w, Matrix<F>& Z, SortType sort=UNSORTED )
.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<F>& Z, SortType sort=UNSORTED )

   Compute the eigenpairs of a Hermitian matrix `A` with eigenvalues lying in 
   the half-open interval :math:`(a,b]`.

Spectral divide and conquer
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary references for this approach is Demmel et al.'s *Fast linear algebra
is stable* and Nakatsukasa et al.'s *Stable and efficient spectral divide and conquer algorithms for the symmetric eigenvalue problem*.

.. cpp:function:: void hermitian_eig::SDC( Matrix<F>& A, Matrix<typename Base<F>::type>& w, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )
.. cpp:function:: void hermitian_eig::SDC( DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )

   Compute the eigenvalues of the matrix :math:`A` via a QDWH-based spectral 
   divide and conquer process. 

   The cutoff controls when the problem is sufficiently small to switch to 
   a standard algorithm, the number of inner iterations is how many attempts 
   to make with the same randomized URV decomposition, and the number of outer 
   iterations is how many random Mobius transformations to try for each spectral
   split before giving up.

.. cpp:function:: void hermitian_eig::SDC( Matrix<F>& A, Matrix<typename Base<F>::type>& w, Matrix<F>& Q, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )
.. cpp:function:: void hermitian_eig::SDC( DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<F>& Q, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )

   Attempt to also compute the eigenvectors.

Skew-Hermitian eigensolver
--------------------------
Essentially identical to the Hermitian eigensolver, :cpp:func:`HermitianEig`;
for any skew-Hermitian matrix :math:`G`, :math:`iG` is Hermitian, as 

.. math::

   (iG)^H = -iG^H = iG.

This fact implies a fast method for solving skew-Hermitian eigenvalue problems:

1. Form :math:`iG` in :math:`O(n^2)` work 
   (switching to complex arithmetic in the real case)
2. Run a Hermitian eigensolve on :math:`iG`, yielding :math:`iG=Z \Lambda Z^H`.
3. Recognize that :math:`G=Z (-i \Lambda) Z^H` provides an EVD of :math:`G`.

Please see the :cpp:func:`HermitianEig` documentation for more details.

.. note::

   Unfortunately, PMRRR currently only supports double-precision problems, and 
   so the parallel versions of these routines are limited to real and complex 
   double-precision matrices.

.. note:: 

   Please see the :ref:`lapack-tuning` section for information on optimizing
   the reduction to tridiagonal form, as it is the dominant cost in all of 
   Elemental's Hermitian eigensolvers.

Full spectrum computation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, Matrix<F>& G, Matrix<typename Base<F>::type>& wImag, SortType sort=UNSORTED )
.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<F>& G, DistMatrix<typename Base<F>::type,VR,STAR>& wImag, SortType sort=UNSORTED )

   Compute the full set of eigenvalues of the skew-Hermitian matrix `G`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, Matrix<F>& G, Matrix<typename Base<F>::type>& wImag, Matrix<Complex<typename Base<F>::type> >& Z, SortType sort=UNSORTED )
.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<F>& G, DistMatrix<typename Base<F>::type,VR,STAR>& wImag, DistMatrix<Complex<typename Base<F>::type> >& Z, SortType sort=UNSORTED )

   Compute the full set of eigenpairs of the skew-Hermitian matrix `G`.

Index-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, Matrix<F>& G, Matrix<typename Base<F>::type>& wImag, int a, int b, SortType sort=UNSORTED )
.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<F>& G, DistMatrix<typename Base<F>::type,VR,STAR>& wImag, int a, int b, SortType sort=UNSORTED )

   Compute the eigenvalues of a skew-Hermitian matrix `G` with
   indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, Matrix<F>& G, Matrix<typename Base<F>::type>& wImag, Matrix<Complex<typename Base<F>::type> >& Z, SortType sort=UNSORTED )
.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<F>& G, DistMatrix<typename Base<F>::type,VR,STAR>& wImag, DistMatrix<Complex<typename Base<F>::type> >& Z, SortType sort=UNSORTED )

   Compute the eigenpairs of a skew-Hermitian matrix `G` with 
   indices in the range :math:`a,a+1,...,b`.

Range-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, Matrix<F>& G, Matrix<typename Base<F>::type>& wImag, typename Base<F>::type a, typename Base<F>::type b, SortType sort=UNSORTED )
.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<F>& G, DistMatrix<typename Base<F>::type,VR,STAR>& wImag, typename Base<F>::type a, typename Base<F>::type b, SortType sort=UNSORTED )

   Compute the eigenvalues of a skew-Hermitian matrix `G` 
   lying in the half-open interval :math:`(a,b]i`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, Matrix<F>& G, Matrix<typename Base<F>::type>& wImag, Matrix<F>& Z, SortType sort=UNSORTED )
.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<F>& G, DistMatrix<typename Base<F>::type,VR,STAR>& wImag, DistMatrix<F>& Z, SortType sort=UNSORTED )

   Compute the eigenpairs of a skew-Hermitian matrix `G` with 
   eigenvalues lying in the half-open interval :math:`(a,b]i`.

Hermitian generalized-definite eigensolvers
-------------------------------------------
The following Hermitian generalized-definite eigenvalue problems frequently 
appear, where both :math:`A` and :math:`B` are Hermitian, and :math:`B` is 
additionally positive-definite:

.. math::

   ABx = \lambda x,

which is denoted with the value ``ABX`` via the 
:cpp:type:`HermitianGenDefiniteEigType` enum,

.. math::

   BAx = \lambda x,

which uses the ``BAX`` value, and finally

.. math::

   Ax = \lambda B x,

which uses the ``AXBX`` enum value.

.. cpp:type:: HermitianGenDefiniteEigType

   An enum for specifying either the ``ABX``, ``BAX``, or ``AXBX`` 
   generalized eigenvalue problems (described above).

Full spectrum computation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B, Matrix<typename Base<F>::type>& w, SortType sort=UNSORTED )
.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<typename Base<F>::type,VR,STAR>& w, SortType sort=UNSORTED )

   Compute the full set of eigenvalues of a generalized EVP involving the 
   Hermitian matrices `A` and `B`, where `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B, Matrix<typename Base<F>::type>& w, Matrix<typename Base<F>::type>& Z, SortType sort=UNSORTED )
.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<typename Base<F>::type>& Z, SortType sort=UNSORTED )

   Compute the full set of eigenpairs of a generalized EVP involving the 
   Hermitian matrices `A` and `B`, where `B` is also positive-definite.

Index-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B, Matrix<typename Base<F>::type>& w, int a, int b, SortType sort=UNSORTED )
.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<typename Base<F>::type,VR,STAR>& w, int a, int b, SortType sort=UNSORTED )

   Compute the eigenvalues with indices in the range :math:`a,a+1,...,b` of a 
   generalized EVP involving the Hermitian matrices `A` and `B`, where `B` is 
   also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B, Matrix<typename Base<F>::type>& w, Matrix<F>& Z, SortType sort=UNSORTED )
.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<F>& Z, SortType sort=UNSORTED )

   Compute the eigenpairs with indices in the range :math:`a,a+1,...,b` of a 
   generalized EVP involving the Hermitian matrices `A` and `B`, where `B` is 
   also positive-definite.

Range-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B, Matrix<typename Base<F>::type>& w, typename Base<F>::type a, typename Base<F>::type b, SortType sort=UNSORTED )
.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<typename Base<F>::type,VR,STAR>& w, typename Base<F>::type a, typename Base<F>::type b, SortType sort=UNSORTED )

   Compute the eigenvalues lying in the half-open interval :math:`(a,b]` of a 
   generalized EVP involving the Hermitian matrices `A` and `B`, where `B` is 
   also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B, Matrix<typename Base<F>::type>& w, Matrix<F>& Z, SortType sort=UNSORTED )
.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<typename Base<F>::type,VR,STAR>& w, DistMatrix<F>& Z, SortType sort=UNSORTED )

   Compute the eigenpairs whose eigenvalues lie in the half-open interval 
   :math:`(a,b]` of a generalized EVP involving the Hermitian matrices `A` and 
   `B`, where `B` is also positive-definite.

Unitary eigensolver
-------------------
Not yet written, will likely be based on Ming Gu's unitary Divide and Conquer 
algorithm for unitary Hessenberg matrices. The spectral divide and conquer 
technique described below should suffice in the meantime.

Normal eigensolver
------------------
Not yet written, will likely be based on Angelika Bunse-Gerstner et al.'s 
Jacobi-like method for simultaneous diagonalization of the commuting Hermitian 
and skew-Hermitian portions of the matrix.
The spectral divide and conquer scheme described below should suffice in the 
meantime.

Schur decomposition
-------------------
Only a prototype spectral divide and conquer implementation is currently 
available, though Elemental will eventually also include an implementation of
Granat et al.'s parallel QR algorithm.

.. cpp:function:: void Schur( Matrix<F>& A )
.. cpp:function:: void Schur( Matrix<F>& A, Matrix<F>& Q )

   Currently defaults to the sequential Hessenberg QR algorithm.

.. cpp:function:: void Schur( DistMatrix<F>& A )
.. cpp:function:: void Schur( DistMatrix<F>& A, DistMatrix<F>& Q )

   Currently defaults to the prototype spectral divide and conquer approach.

Hessenberg QR algorithm
^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void schur::QR( Matrix<F>& A, Matrix<Complex<typename Base<F>::type>>& w )
.. cpp:function:: void schur::QR( Matrix<F>& A, Matrix<Complex<typename Base<F>::type>>& w, Matrix<F>& Q )

   Use a sequential QR algorithm to compute the Schur decomposition.

Spectral divide and conquer
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary reference for this approach is Demmel et al.'s *Fast linear algebra
is stable*. While the current implementation needs a large number of algorithmic
improvements, especially with respect to choosing the Mobius transformations,
it tends to succeed on random matrices.

.. cpp:function:: void schur::SDC( Matrix<F>& A, Matrix<Complex<typename Base<F>::type>>& w, bool formATR=false, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )
.. cpp:function:: void schur::SDC( DistMatrix<F>& A, DistMatrix<Complex<typename Base<F>::type>,VR,STAR>& w, bool formATR=false, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )

   Compute the eigenvalues of the matrix :math:`A` via a spectral divide and
   conquer process. On exit, the eigenvalues of :math:`A` will be stored on its
   diagonal, and, if ``formATR`` was set to true, the upper triangle of 
   :math:`A` will be its corresponding upper-triangular Schur factor.

   The cutoff controls when the problem is sufficiently small to switch to 
   a sequential Hessenberg QR algorithm, the number of inner iterations is 
   how many attempts to make with the same randomized URV decomposition, and 
   the number of outer iterations is how many random Mobius transformations to
   try for each spectral split before giving up.

.. cpp:function:: void schur::SDC( Matrix<F>& A, Matrix<Complex<typename Base<F>::type>>& w, Matrix<F>& Q, bool formATR=true, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )
.. cpp:function:: void schur::SDC( DistMatrix<F>& A, DistMatrix<Complex<typename Base<F>::type>,VR,STAR>& w, DistMatrix<F>& Q, bool formATR=true, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, typename Base<F>::type relTol=0 )

   Attempt to also compute the Schur vectors.

Hermitian SVD
-------------
Given an eigenvalue decomposition of a Hermitian matrix :math:`A`, say

.. math::

   A = V \Lambda V^H,

where :math:`V` is unitary and :math:`\Lambda` is diagonal and real. 
Then an SVD of :math:`A` can easily be computed as

.. math::

   A = U |\Lambda| V^H,

where the columns of :math:`U` equal the columns of :math:`V`, modulo sign 
flips introduced by negative eigenvalues.

.. cpp:function:: void HermitianSVD( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& U, Matrix<F>& V )
.. cpp:function:: void HermitianSVD( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, DistMatrix<F>& U, DistMatrix<F>& V )

   Return a vector of singular values, :math:`s`, and the left and right 
   singular vector matrices, :math:`U` and :math:`V`, such that 
   :math:`A=U \mathrm{diag}(s) V^H`.

.. cpp:function:: void HermitianSVD( UpperOrLower uplo, Matrix<F>& A, Matrix<typename Base<F>::type>& s )
.. cpp:function:: void HermitianSVD( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )

   Return the singular values of :math:`A` in `s`. Note that the appropriate 
   triangle of `A` is overwritten during computation.

Polar decomposition
-------------------
Every matrix :math:`A` can be written as :math:`A=QP`, where :math:`Q` is 
unitary and :math:`P` is Hermitian and positive semi-definite. This is known as
the *polar decomposition* of :math:`A` and can be constructed as 
:math:`Q := U V^H` and :math:`P := V \Sigma V^H`, where 
:math:`A = U \Sigma V^H` is the SVD of :math:`A`. Alternatively, it can be 
computed through (a dynamically-weighted) Halley iteration.

.. cpp:function:: void Polar( Matrix<F>& A )
.. cpp:function:: void Polar( DistMatrix<F>& A )
.. cpp:function:: void Polar( Matrix<F>& A, Matrix<F>& P )
.. cpp:function:: void Polar( DistMatrix<F>& A, DistMatrix<F>& P )

   Compute the polar decomposition of :math:`A`, :math:`A=QP`, returning 
   :math:`Q` within `A` and :math:`P` within `P`. The current implementation
   first computes the SVD.

.. cpp:function:: void HermitianPolar( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void HermitianPolar( UpperOrLower uplo, DistMatrix<F>& A )
.. cpp:function:: void HermitianPolar( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P )
.. cpp:function:: void HermitianPolar( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P )

   Compute the polar decomposition through a Hermitian EVD. Since this is 
   equivalent to a Hermitian sign decomposition (if :math:`\text{sgn}(0)` is 
   set to 1), these routines are equivalent to `HermitianSign`.

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:function:: int polar::QDWH( Matrix<F>& A, bool colPiv=false, int maxits=20 )
.. cpp:function:: int polar::QDWH( DistMatrix<F>& A, bool colPiv=false, int maxIts=20 )
.. cpp:function:: int hermitian_polar::QDWH( UpperOrLower uplo, Matrix<F>& A, bool colPiv=false, int maxits=20 )
.. cpp:function:: int hermitian_polar::QDWH( UpperOrLower uplo, DistMatrix<F>& A, bool colPiv=false, int maxIts=20 )

   Overwrites :math:`A` with the :math:`Q` from the polar decomposition using 
   a QR-based dynamically weighted Halley iteration. The number of iterations
   used is returned upon completion.
   **TODO: reference to Yuji's paper**

.. cpp:function:: int polar::QDWH( Matrix<F>& A, Matrix<F>& P, bool colPiv=false, int maxits=20 )
.. cpp:function:: int polar::QDWH( DistMatrix<F>& A, DistMatrix<F>& P, bool colPiv=false, int maxIts=20 )
.. cpp:function:: int hermitian_polar::QDWH( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, bool colPiv=false, int maxits=20 )
.. cpp:function:: int hermitian_polar::QDWH( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P, bool colPiv=false, int maxIts=20 )

   Return the full polar decomposition, where :math:`P` is HPD.

SVD
---
Given a general matrix :math:`A`, the *Singular Value Decomposition* is the 
triplet :math:`(U,\Sigma,V)` such that

.. math::

   A = U \Sigma V^H,

where :math:`U` and :math:`V` are unitary, and :math:`\Sigma` is diagonal with 
non-negative entries.

.. cpp:function:: void SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V )

.. cpp:function:: void SVD( DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, DistMatrix<F>& V )

   Overwrites `A` with :math:`U`, `s` with the diagonal entries of :math:`\Sigma`, and `V` with :math:`V`. 

.. cpp:function:: void SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s )

.. cpp:function:: void SVD( DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )

   Forms the singular values of :math:`A` in `s`. Note that `A` is overwritten in order to compute the singular values.

