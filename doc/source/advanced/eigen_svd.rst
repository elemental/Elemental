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

.. cpp:function:: void RealHermitianFunction( UpperOrLower uplo, DistMatrix<R,MC,MR>& A, const RealFunctor& f )

   Modifies the eigenvalues of the passed-in real symmetric matrix by replacing 
   each eigenvalue :math:`\omega_i` with :math:`f(\omega_i) \in \mathbb{R}`. 
   ``RealFunctor`` is any 
   class which has the member function ``R operator()( R omega ) const``.
   See `examples/advanced/RealSymmetricFunction.cpp <../../../../examples/advanced/RealSymmetricFunction.cpp>`_ for an example usage.

.. cpp:function:: void RealHermitianFunction( UpperOrLower uplo, DistMatrix<Complex<R>,MC,MR>& A, const RealFunctor& f )

   Modifies the eigenvalues of the passed-in complex Hermitian matrix by 
   replacing each eigenvalue :math:`\omega_i` with 
   :math:`f(\omega_i) \in \mathbb{R}`. 
   ``RealFunctor`` can be any class which has the member function 
   ``R operator()( R omega ) const``.
   See `examples/advanced/RealHermitianFunction.cpp <../../../../examples/advanced/RealHermitianFunction.cpp>`_ for an example usage.

.. cpp:function:: void ComplexHermitianFunction( UpperOrLower uplo, DistMatrix<Complex<R>,MC,MR>& A, const ComplexFunctor& f )

   Modifies the eigenvalues of the passed-in complex Hermitian matrix by
   replacing each eigenvalue :math:`\omega_i` with 
   :math:`f(\omega_i) \in \mathbb{C}`. ``ComplexFunctor`` can be any class
   which has the member function ``Complex<R> operator()( R omega ) const``.
   See `examples/advanced/ComplexHermitianFunction.cpp <../../../../examples/advanced/ComplexHermitianFunction.cpp>`_ for an example usage.

Skew-Hermitian eigensolver
--------------------------
**TODO:** Describe :math:`Gx=\lambda x` and ``SkewHermitianEig`` here.

Hermitian generalized-definite eigensolvers
-------------------------------------------
**TODO:** Describe :math:`ABx=\lambda x`, :math:`BAx=\lambda x`, and 
:math:`Ax=\lambda Bx` cases and ``HermitianGenDefiniteEig``.

Unitary eigensolver
-------------------
Not yet written, will likely be based on Ming Gu's unitary Divide and Conquer 
algorithm for unitary Hessenberg matrices.

Normal eigensolver
------------------
Not yet written, will likely be based on Angelika Bunse-Gerstner et al.'s 
Jacobi-like method for simultaneous diagonalization of the commuting Hermitian 
and skew-Hermitian portions of the matrix.

General eigensolver
-------------------
Not yet written, will likely eventually include Greg Henry et al.'s and 
Robert Granat et al.'s approaches.

Hermitian SVD
-------------
Not yet written, but relatively trivial, as the SVD of a Hermitian matrix can 
easily be computed from its eigenvalue decomposition.

General SVD
-----------
Not yet written; will likely be based on the Ming Gu's approach to 
Divide and Conquer algorithm for the bidiagonal SVD.
