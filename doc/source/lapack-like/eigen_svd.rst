Eigensolvers and SVD
====================

Hermitian eigensolver
---------------------
Elemental provides a collection of routines for both full and partial 
solution of the Hermitian eigenvalue problem 

.. math::

   A Z = Z \Omega,

where `A` is the given Hermitian matrix, and unitary `Z` and real diagonal 
:math:`\Omega` are sought. In particular, with the eigenvalues and 
corresponding eigenpairs labeled in non-decreasing order, the three basic 
modes are:

1. Compute all eigenvalues or eigenpairs, :math:`\{\omega_i\}_{i=0}^{n-1}` or 
   :math:`\{(x_i,\omega_i)\}_{i=0}^{n-1}`.
2. Compute the eigenvalues or eigenpairs with a given range of indices, say  
   :math:`\{\omega_i\}_{i=a}^b` or :math:`\{(x_i,\omega_i)\}_{i=a}^b`, 
   with :math:`0 \le a \le b < n`.
3. Compute all eigenpairs (or just eigenvalues) with eigenvalues lying in a 
   particular half-open interval, either
   :math:`\{\omega_i \;|\; \omega_i \in (a,b] \}` or 
   :math:`\{ (x_i,\omega_i) \;|\; \omega_i \in (a,b] \}`.

As of now, all three approaches start with Householder tridiagonalization 
(ala ``HermitianTridiag``) and then call Matthias Petschow and 
Paolo Bientinesi's PMRRR for the tridiagonal eigenvalue problem.

.. note:: 

   Please see the *Tuning parameters* section for information on optimizing 
   the reduction to tridiagonal form, as it is the dominant cost in all of 
   Elemental's Hermitian eigensolvers.

Full spectrum computation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,VR,STAR>& w )

   Compute the full set of eigenvalues of the double-precision real symmetric 
   distributed matrix `A`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<double,VR,STAR>& w )

   Compute the full set of eigenvalues of the double-precision complex
   Hermitian distributed matrix `A`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the full set of eigenpairs of the double-precision real symmetric 
   distributed matrix `A`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the full set of eigenpairs of the double-precision complex Hermitian
   distributed matrix `A`.

Index-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,VR,STAR>& w, int a, int b )

   Compute the eigenvalues of a double-precision real symmetric distributed 
   matrix `A` with indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<double,VR,STAR>& w, int a, int b )

   Compute the eigenvalues of a double-precision complex Hermitian distributed 
   matrix `A` with indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z, int a, int b )

   Compute the eigenpairs of a double-precision real symmetric distributed 
   matrix `A` with indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the eigenpairs of a double-precision complex Hermitian distributed 
   matrix `A` with indices in the range :math:`a,a+1,...,b`.

Range-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,VR,STAR>& w, double a, double b )

   Compute the eigenvalues of a double-precision real symmetric distributed 
   matrix `A` lying in the half-open interval :math:`(a,b]`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<double,VR,STAR>& w, double a, double b )

   Compute the eigenvalues of a double-precision complex Hermitian distributed 
   matrix `A` lying in the half-open interval :math:`(a,b]`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z, double a, double b )

   Compute the eigenpairs of a double-precision real symmetric distributed 
   matrix `A` with eigenvalues lying in the half-open interval :math:`(a,b]`.

.. cpp:function:: void HermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the eigenpairs of a double-precision complex Hermitian distributed 
   matrix `A` with eigenvalues lying in the half-open interval :math:`(a,b]`.

Sorting the eigenvalues/eigenpairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Since extra time is required in order to sort the eigenvalues/eigenpairs, 
they are not sorted by default. However, this can be remedied by the appropriate
routine from the following list:

.. cpp:function:: void SortEig( DistMatrix<R,VR,STAR>& w )

   Sort a column-vector of real eigenvalues into non-decreasing order.

.. cpp:function:: void SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<R,MC,MR>& Z )

   Sort a set of real eigenpairs into non-decreasing order (based on the 
   eigenvalues).

.. cpp:function:: void SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<Complex<R>,MC,MR>& Z )

   Sort a set of real eigenvalues and complex eigenvectors into non-decreasing
   order (based on the eigenvalues).

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

.. cpp:function:: void RealHermitianFunction( UpperOrLower uplo, DistMatrix<F,MC,MR>& A, const RealFunctor& f )

   Modifies the eigenvalues of the passed-in Hermitian matrix by replacing 
   each eigenvalue :math:`\omega_i` with :math:`f(\omega_i) \in \mathbb{R}`. 
   ``RealFunctor`` is any 
   class which has the member function ``R operator()( R omega ) const``.
   See `examples/lapack-like/RealSymmetricFunction.cpp <../../../../examples/lapack-like/RealHermitianFunction.cpp>`_ for an example usage.

.. cpp:function:: void ComplexHermitianFunction( UpperOrLower uplo, DistMatrix<Complex<R>,MC,MR>& A, const ComplexFunctor& f )

   Modifies the eigenvalues of the passed-in complex Hermitian matrix by
   replacing each eigenvalue :math:`\omega_i` with 
   :math:`f(\omega_i) \in \mathbb{C}`. ``ComplexFunctor`` can be any class
   which has the member function ``Complex<R> operator()( R omega ) const``.
   See `examples/lapack-like/ComplexHermitianFunction.cpp <../../../../examples/lapack-like/ComplexHermitianFunction.cpp>`_ for an example usage.

Skew-Hermitian eigensolver
--------------------------
Essentially identical to the Hermitian eigensolver, ``HermitianEig``;
for any skew-Hermitian matrix :math:`G`, :math:`iG` is Hermitian, as 

.. math::

   (iG)^H = -iG^H = iG.

This fact implies a fast method for solving skew-Hermitian eigenvalue problems:

1. Form :math:`iG` in :math:`O(n^2)` work 
   (switching to complex arithmetic in the real case)
2. Run a Hermitian eigensolve on :math:`iG`, yielding :math:`iG=Z \Omega Z^H`.
3. Recognize that :math:`G=Z (-i \Omega) Z^H` provides an EVD of :math:`G`.

Please see the ``HermitianEig`` documentation for more details.

.. note:: 

   Please see the *Tuning parameters* section for information on optimizing 
   the reduction to tridiagonal form, as it is the dominant cost in all of 
   Elemental's Hermitian eigensolvers.

Full spectrum computation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag )

   Compute the full set of eigenvalues of the double-precision real 
   skew-symmetric distributed matrix `G`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag )

   Compute the full set of eigenvalues of the double-precision complex
   skew-Hermitian distributed matrix `G`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, DistMatrix<Complex<double>,MC,MR>& Z )

   Compute the full set of eigenpairs of the double-precision real 
   skew-symmetric distributed matrix `G`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, DistMatrix<Complex<double>,MC,MR>& Z )

   Compute the full set of eigenpairs of the double-precision complex 
   skew-Hermitian distributed matrix `G`.

Index-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, int a, int b )

   Compute the eigenvalues of a double-precision real skew-symmetric 
   distributed matrix `G` with indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, int a, int b )

   Compute the eigenvalues of a double-precision complex skew-Hermitian 
   distributed matrix `G` with indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, DistMatrix<Complex<double>,MC,MR>& Z, int a, int b )

   Compute the eigenpairs of a double-precision real skew-symmetric distributed 
   matrix `G` with indices in the range :math:`a,a+1,...,b`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, DistMatrix<Complex<double>,MC,MR>& Z )

   Compute the eigenpairs of a double-precision complex skew-Hermitian 
   distributed matrix `G` with indices in the range :math:`a,a+1,...,b`.

Range-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, double a, double b )

   Compute the eigenvalues of a double-precision real skew-symmetric distributed
   matrix `G` lying in the half-open interval :math:`(a,b]i`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, double a, double b )

   Compute the eigenvalues of a double-precision complex skew-Hermitian 
   distributed matrix `G` lying in the half-open interval :math:`(a,b]i`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<double,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, DistMatrix<Complex<double>,MC,MR>& Z, double a, double b )

   Compute the eigenpairs of a double-precision real skew-symmetric distributed 
   matrix `G` with eigenvalues lying in the half-open interval :math:`(a,b]i`.

.. cpp:function:: void SkewHermitianEig( UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& G, DistMatrix<double,VR,STAR>& wImag, DistMatrix<Complex<double>,MC,MR>& Z )

   Compute the eigenpairs of a double-precision complex skew-Hermitian 
   distributed matrix `G` with eigenvalues lying in the half-open interval 
   :math:`(a,b]i`.

Hermitian generalized-definite eigensolvers
-------------------------------------------
The following Hermitian generalized-definite eigenvalue problems frequently 
appear, where both :math:`A` and :math:`B` are Hermitian, and :math:`B` is 
additionally positive-definite:

.. math::

   ABx = \omega x,

which is denoted with the value ``ABX`` via the ``HermitianGenDefiniteEigType``
enum,

.. math::

   BAx = \omega x,

which uses the ``BAX`` value, and finally

.. math::

   Ax = \omega B x,

which uses the ``AXBX`` enum value.

.. cpp:type:: HermitianGenDefiniteEigType

   An enum for specifying either the ``ABX``, ``BAX``, or ``AXBX`` 
   generalized eigenvalue problems (described above).

Full spectrum computation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B, DistMatrix<double,VR,STAR>& w )

   Compute the full set of eigenvalues of a generalized EVP involving the 
   double-precision real symmetric distributed matrices `A` and `B`, where 
   `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<Complex<double>,MC,MR>& B, DistMatrix<double,VR,STAR>& w )

   Compute the full set of eigenvalues of a generalized EVP involving the 
   double-precision complex Hermitian distributed matrices `A` and `B`, where
   `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the full set of eigenpairs of a generalized EVP involving the 
   double-precision real symmetric distributed matrices `A` and `B`, where 
   `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<Complex<double>,MC,MR>& B, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the full set of eigenpairs of a generalized EVP involving the 
   double-precision complex Hermitian distributed matrices `A` and `B`, where 
   `B` is also positive-definite.

Index-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B, DistMatrix<double,VR,STAR>& w, int a, int b )

   Compute the eigenvalues with indices in the range :math:`a,a+1,...,b` of a 
   generalized EVP involving the double-precision real symmetric distributed 
   matrices `A` and `B`, where `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<Complex<double>,MC,MR>& B, DistMatrix<double,VR,STAR>& w, int a, int b )

   Compute the eigenvalues with indices in the range :math:`a,a+1,...,b` of a 
   generalized EVP involving the double-precision complex Hermitian distributed 
   matrices `A` and `B`, where `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z, int a, int b )

   Compute the eigenpairs with indices in the range :math:`a,a+1,...,b` of a 
   generalized EVP involving the double-precision real symmetric distributed 
   matrices `A` and `B`, where `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<Complex<double>,MC,MR>& B, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the eigenpairs with indices in the range :math:`a,a+1,...,b` of a 
   generalized EVP involving the double-precision complex Hermitian distributed 
   matrices `A` and `B`, where `B` is also positive-definite.

Range-based subset computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B, DistMatrix<double,VR,STAR>& w, double a, double b )

   Compute the eigenvalues lying in the half-open interval :math:`(a,b]` of a 
   generalized EVP involving the double-precision real symmetric distributed 
   matrices `A` and `B`, where `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<Complex<double>,MC,MR>& B, DistMatrix<double,VR,STAR>& w, double a, double b )

   Compute the eigenvalues lying in the half-open interval :math:`(a,b]` of a 
   generalized EVP involving the double-precision complex Hermitian distributed 
   matrices `A` and `B`, where `B` is also positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z, double a, double b )

   Compute the eigenpairs whose eigenvalues lie in the half-open 
   interval :math:`(a,b]` of a generalized EVP involving the double-precision 
   real symmetric distributed matrices `A` and `B`, where `B` is also 
   positive-definite.

.. cpp:function:: void HermitianGenDefiniteEig( HermitianGenDefiniteEigType type, UpperOrLower uplo, DistMatrix<Complex<double>,MC,MR>& A, DistMatrix<Complex<double>,MC,MR>& B, DistMatrix<double,VR,STAR>& w, DistMatrix<double,MC,MR>& Z )

   Compute the eigenpairs whose eigenvalues lie in the half-open interval 
   :math:`(a,b]` of a generalized EVP involving the double-precision complex 
   Hermitian distributed matrices `A` and `B`, where `B` is also 
   positive-definite.

Unitary eigensolver
-------------------
Not yet written, will likely be based on Ming Gu's unitary Divide and Conquer 
algorithm for unitary Hessenberg matrices.

Normal eigensolver
------------------
Not yet written, will likely be based on Angelika Bunse-Gerstner et al.'s 
Jacobi-like method for simultaneous diagonalization of the commuting Hermitian 
and skew-Hermitian portions of the matrix.

Schur decomposition
-------------------
Not yet written, will likely eventually include Greg Henry et al.'s and 
Robert Granat et al.'s approaches.

Hermitian SVD
-------------
Given an eigenvalue decomposition of a Hermitian matrix :math:`A`, say

.. math::

   A = V \Lambda V^H,

where :math:`V` is unitary and :math:`\Lambda` is diagonal and real. 
Then an SVD of :math:`A` can easily be computed as

.. math::

   A = U |\Lambda| V^H,

where the columns of :math:`U` equal the columns of :math:`V`, modulo sign flips introduced by negative eigenvalues.

.. cpp:function:: void HermitianSVD( UpperOrLower uplo, DistMatrix<F,MC,MR>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, DistMatrix<F,MC,MR>& U, DistMatrix<F,MC,MR>& V )

   Return a vector of singular values, :math:`s`, and the left and right singular vector matrices, :math:`U` and :math:`V`, such that :math:`A=U \mathrm{diag}(s) V^H`.

.. cpp:function:: void HermitianSingularValues( UpperOrLower uplo, DistMatrix<F,MC,MR>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )

   Return the singular values of :math:`A` in `s`. Note that the appropriate triangle of `A` is overwritten during computation.

Square root
-----------
Hermitian matrices with non-negative eigenvalues have a
natural matrix square root through their eigenvalue decomposition. This routine
attempts to compute said matrix square root and throws a
``NonHPSDMatrixException`` if any sufficiently negative eigenvalues are
computed.

.. cpp:function:: void HPSDSquareRoot( UpperOrLower uplo, DistMatrix<F,MC,MR>& A )

   Overwrites the Hermitian positive semi-definite distributed matrix `A` with
   its matrix square root.

SVD
---
Given a general matrix :math:`A`, the *Singular Value Decomposition* is the 
triplet :math:`(U,\Sigma,V)` such that

.. math::

   A = U \Sigma V^H,

where :math:`U` and :math:`V` are unitary, and :math:`\Sigma` is diagonal with 
non-negative entries.

.. cpp:function:: void SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V )

.. cpp:function:: void SVD( DistMatrix<F,MC,MR>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, DistMatrix<F,MC,MR>& V )

   Overwrites `A` with :math:`U`, `s` with the diagonal entries of :math:`\Sigma`, and `V` with :math:`V`. 

   .. note:: 

      This routine has not yet been optimized for cases where the matrix is 
      highly rectangular. It is also far from the most efficient 
      possible algorithm, as it currently uses the QR algorithm for the
      bidiagonal SVD (as opposed to a Divide and Conquer algorithm) and 
      does not apply Householder reflectors in the most efficient manner 
      possible.

.. cpp:function:: void SingularValues( Matrix<F>& A, Matrix<typename Base<F>::type>& s )

.. cpp:function:: void SingularValues( DistMatrix<F,MC,MR>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )

   Forms the singular values of :math:`A` in `s`. Note that `A` is overwritten in order to compute the singular values.

   .. note::

      This routine has not yet been optimized for cases where the matrix is 
      highly rectangular. It is also far from the most efficient
      possible algorithm, as it currently uses the QR algorithm to redundantly
      compute the singular values of the bidiagonal matrix (as opposed to 
      bisection).

Polar decomposition
-------------------
Every matrix :math:`A` can be written as :math:`A=QP`, where :math:`Q` is unitary and :math:`P` is Hermitian and positive semi-definite. This is known as the *polar decomposition* of :math:`A` and can be constructed as :math:`Q := U V^H` and :math:`P := V \Sigma V^H`, where :math:`A = U \Sigma V^H` is the SVD of :math:`A`.

.. cpp:function:: void Polar( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& P )

   Compute the polar decomposition of :math:`A`, :math:`A=QP`, returning 
   :math:`Q` within `A` and :math:`P` within `P`.
