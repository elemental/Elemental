Invariants, inner products, norms, etc.
=======================================

Condition number
----------------
The condition number of a matrix with respect to a particular norm is

.. math::

   \kappa(A) = \|A\| \|A^{-1}\|,

with the most common choice being the matrix two-norm.

.. cpp:function:: Base<F> Condition( const Matrix<F>& A, NormType type=TWO_NORM )
.. cpp:function:: Base<F> Condition( const DistMatrix<F,U,V>& A, NormType type=TWO_NORM )

   Returns the condition number with respect to the specified norm 
   (one, two, or Frobenius).

.. cpp:function:: Base<F> FrobeniusCondition( const Matrix<F>& A )
.. cpp:function:: Base<F> FrobeniusCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the Frobenius norm.

.. cpp:function:: Base<F> InfinityCondition( const Matrix<F>& A )
.. cpp:function:: Base<F> InfinityCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the infinity norm.

.. cpp:function:: Base<F> MaxCondition( const Matrix<F>& A )
.. cpp:function:: Base<F> MaxCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the entrywise maximum norm.

.. cpp:function:: Base<F> OneCondition( const Matrix<F>& A )
.. cpp:function:: Base<F> OneCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the one norm.

.. cpp:function:: Base<F> TwoCondition( const Matrix<F>& A )
.. cpp:function:: Base<F> TwoCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the two norm.

Determinant
-----------
Though there are many different possible definitions of the determinant of a 
matrix :math:`A \in \mathbb{F}^{n \times n}`, the simplest one is in terms of 
the product of the eigenvalues (including multiplicity):

.. math::

   \mbox{det}(A) = \prod_{i=0}^{n-1} \lambda_i.

Since :math:`\mbox{det}(AB)=\mbox{det}(A)\mbox{det}(B)`, we can compute the 
determinant of an arbitrary matrix in :math:`\mathcal{O}(n^3)` work by 
computing its LU decomposition (with partial pivoting), :math:`PA=LU`, 
recognizing that :math:`\mbox{det}(P)=\pm 1` 
(the *signature* of the permutation), and computing

.. math::

   \mbox{det}(A) = \mbox{det}(P)\mbox{det}(L)\mbox{det}(U) 
                 = \mbox{det}(P) \prod_{i=0}^{n-1} \upsilon_{i,i}
                 = \pm \prod_{i=0}^{n-1} \upsilon_{i,i},

where :math:`\upsilon_{i,i}` is the i'th diagonal entry of :math:`U`.

.. cpp:function:: F Determinant( const Matrix<F>& A )
.. cpp:function:: F Determinant( const DistMatrix<F>& A )
.. cpp:function:: F Determinant( Matrix<F>& A, bool canOverwrite=false )
.. cpp:function:: F Determinant( DistMatrix<F>& A, bool canOverwrite=false )

   The determinant of the (fully populated) square matrix `A`.
   Some of the variants allow for overwriting the input matrix in order to 
   avoid forming another temporary matrix.

.. cpp:type:: struct SafeProduct<F>

   Represents the product of `n` values as :math:`\rho \exp(\kappa n)`, 
   where :math:`|\rho| \le 1` and :math:`\kappa \in \mathbb{R}`.

   .. cpp:member:: F rho

      For nonzero values, `rho` is the modulus and lies *on* the unit 
      circle; in order to represent a value that is precisely zero, `rho` 
      is set to zero.

   .. cpp:member:: Base<F> kappa

      `kappa` can be an arbitrary real number.

   .. cpp:member:: int n

      The number of values in the product.

.. cpp:function:: SafeProduct<F> SafeDeterminant( const Matrix<F>& A )
.. cpp:function:: SafeProduct<F> SafeDeterminant( const DistMatrix<F>& A )
.. cpp:function:: SafeProduct<F> SafeDeterminant( Matrix<F>& A, bool canOverwrite=false )
.. cpp:function:: SafeProduct<F> SafeDeterminant( DistMatrix<F>& A, bool canOverwrite=false )

   The determinant of the square matrix `A` in an expanded form 
   which is less likely to over/under-flow.

HPDDeterminant
--------------
A version of the above determinant specialized for Hermitian positive-definite
matrices (which will therefore have all positive eigenvalues and a positive 
determinant).

.. cpp:function:: Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: Base<F> HPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
.. cpp:function:: Base<F> HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )

   The determinant of the (fully populated) Hermitian positive-definite
   matrix `A`.
   Some of the variants allow for overwriting the input matrix in order to 
   avoid forming another temporary matrix.

.. cpp:function:: SafeProduct<F> SafeHPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: SafeProduct<F> SafeHPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: SafeProduct<F> SafeHPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
.. cpp:function:: SafeProduct<F> SafeHPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )

   The determinant of the Hermitian positive-definite matrix `A` in an 
   expanded form which is less likely to over/under-flow.

Norm
----
The following routines can return either
:math:`\|A\|_1`, :math:`\|A\|_\infty`, :math:`\|A\|_F` (the Frobenius norm),
the maximum entrywise norm, :math:`\|A\|_2`, or :math:`\|A\|_*`
(the nuclear/trace norm) of fully-populated matrices.

.. cpp:function:: Base<F> Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: Base<F> Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM )

   Assumes that the input is a fully-specified matrix.

.. cpp:function:: Base<F> HermitianNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: Base<F> HermitianNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: Base<F> SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: Base<F> SymmetricNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )

   Same as :cpp:func:`Norm`, but the matrix is implicitly
   Hermitian/symmetric with the data stored in the triangle specified by
   :cpp:type:`UpperOrLower`.
   Also, while :cpp:func:`Norm` supports every type of distribution,
   :cpp:func:`HermitianNorm`/:cpp:func:`SymmetricNorm` currently only supports 
   the standard matrix distribution.

Alternatively, one may directly call the following routines (note that the entrywise, KyFan, and Schatten norms have an extra parameter and must be called 
directly).

.. cpp:function:: Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> EntrywiseNorm( const DistMatrix<F,U,V>& A, Base<F> p )
.. cpp:function:: Base<F> HermitianEntrywiseNorm( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> HermitianEntrywiseNorm( UpperOrLower uplo, const DistMatrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> SymmetricEntrywiseNorm( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> SymmetricEntrywiseNorm( UpperOrLower uplo, const DistMatrix<F>& A, Base<F> p )

   The :math:`\ell_p` norm of the columns of `A` stacked into a single vector. 
   Note that the Frobenius norm corresponds to the :math:`p=2` case.

.. cpp:function:: Base<F> EntrywiseOneNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> EntrywiseOneNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: Base<F> SymmetricEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The :math:`\ell_1` norm of the columns of `A` stacked into a single vector. 

.. cpp:function:: Base<F> FrobeniusNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> FrobeniusNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The :math:`\ell_2` norm of the singular values (the Schatten norm with 
   :math:`p=2`), which can be cheaply computed as the :math:`\ell_2` norm of 
   :math:`\text{vec}(A)`.

.. cpp:function:: Base<F> KyFanNorm( const Matrix<F>& A, int k )
.. cpp:function:: Base<F> KyFanNorm( const DistMatrix<F,U,V>& A, int k )
.. cpp:function:: Base<F> HermitianKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, int k )
.. cpp:function:: Base<F> HermitianKyFanNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, int k )
.. cpp:function:: Base<F> SymmetricKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, int k )
.. cpp:function:: Base<F> SymmetricKyFanNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, int k )

   The sum of the largest `k` singular values.

.. cpp:function:: Base<F> InfinityNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> InfinityNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianInfinityNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: Base<F> SymmetricInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricInfinityNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The maximum :math:`\ell_1` norm of the rows of `A`. In the symmetric and 
   Hermitian cases, this is equivalent to the :math:`\|\cdot \|_1` norm.

.. cpp:function:: Base<F> MaxNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> MaxNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: Base<F> SymmetricMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The maximum absolute value of the matrix entries.

.. cpp:function:: Base<F> NuclearNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> NuclearNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )

   The sum of the singular values. This is equivalent to both the KyFan 
   norm with :math:`k=n` and the Schatten norm with :math:`p=1`.
   Note that the nuclear norm is dual to the two-norm, which is the 
   Schatten norm with :math:`p=\infty`.

.. cpp:function:: Base<F> OneNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> OneNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: Base<F> SymmetricOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The maximum :math:`\ell_1` norm of the columns of `A`. In the symmetric and 
   Hermitian cases, this is equivalent to the :math:`\| \cdot \|_\infty` norm.

.. cpp:function:: Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> SchattenNorm( const DistMatrix<F,U,V>& A, Base<F> p )
.. cpp:function:: Base<F> HermitianSchattenNorm( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> HermitianSchattenNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p )
.. cpp:function:: Base<F> SymmetricSchattenNorm( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
.. cpp:function:: Base<F> SymmetricSchattenNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p )

   The :math:`\ell_p` norm of the singular values.

.. cpp:function:: Base<F> TwoNorm( const Matrix<F>& A )
.. cpp:function:: Base<F> TwoNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> HermitianTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
.. cpp:function:: Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: Base<F> SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )

   The maximum singular value. This is equivalent to the KyFan norm with `k` 
   equal to one and the Schatten norm with :math:`p=\infty`.

.. cpp:function:: int ZeroNorm( const Matrix<F>& A )
.. cpp:function:: int ZeroNorm( const DistMatrix<F>& A )
.. cpp:function:: int HermitianZeroNorm( const Matrix<F>& A )
.. cpp:function:: int HermitianZeroNorm( const DistMatrix<F>& A )
.. cpp:function:: int SymmetricZeroNorm( const Matrix<F>& A )
.. cpp:function:: int SymmetricZeroNorm( const DistMatrix<F>& A )

   Return the number of nonzero entries in the matrix.

Two-norm estimates
------------------

.. cpp:function:: Base<F> TwoNormEstimate( Matrix<F>& A, Base<F> tol=1e-6 )
.. cpp:function:: Base<F> TwoNormEstimate( DistMatrix<F>& A, Base<F> tol=1e-6 )
.. cpp:function:: Base<F> HermitianTwoNormEstimate( Matrix<F>& A, Base<F> tol=1e-6 )
.. cpp:function:: Base<F> HermitianTwoNormEstimate( DistMatrix<F>& A, Base<F> tol=1e-6 )
.. cpp:function:: Base<F> SymmetricTwoNormEstimate( Matrix<F>& A, Base<F> tol=1e-6 )
.. cpp:function:: Base<F> SymmetricTwoNormEstimate( DistMatrix<F>& A, Base<F> tol=1e-6 )

   Return an estimate for the two-norm which should be accurate within a 
   factor of :math:`n` times the specified tolerance.

Trace
-----
The two equally useful definitions of the trace of a square matrix 
:math:`A \in \mathbb{F}^{n \times n}` are

.. math::

   \mbox{tr}(A) = \sum_{i=0}^{n-1} \alpha_{i,i} = \sum_{i=0}^{n-1} \lambda_i,

where :math:`\alpha_{i,i}` is the i'th diagonal entry of :math:`A` and 
:math:`\lambda_i` is the i'th eigenvalue (counting multiplicity) of :math:`A`.

Clearly the former equation is easier to compute, but the latter is an 
important characterization.

.. cpp:function:: F Trace( const Matrix<F>& A )
.. cpp:function:: F Trace( const DistMatrix<F>& A )

   Return the trace of the square matrix `A`.

Hadamard
--------
The Hadamard product of two :math:`m \times n` matrices :math:`A` and 
:math:`B` is given entrywise by :math:`\alpha_{i,j} \beta_{i,j}` and denoted
by :math:`C = A \circ B`.

.. cpp:function:: void Hadamard( const Matrix<F>& A, const Matrix<F>& B, Matrix<F>& C )
.. cpp:function:: void Hadamard( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B, DistMatrix<F,U,V>& C )

HermitianInertia
----------------

.. cpp:function:: void HermitianInertia( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_PARLETT )
.. cpp:function:: void HermitianInertia( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_PARLETT )

   Returns the triplet containing the number of positive, negative, and 
   zero eigenvalues of the Hermitian matrix by analyzing the block diagonal
   resulting from a pivoted LDL factorization.

HilbertSchmidt
--------------
The Hilbert-Schmidt inner-product of two :math:`m \times n` matrices :math:`A`
and :math:`B` is :math:`\mbox{tr}(A^H B)`.

.. cpp:function:: F HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B )
.. cpp:function:: F HilbertSchmidt( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
