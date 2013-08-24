Invariants, inner products, norms, etc.
=======================================

Condition number
----------------
The condition number of a matrix with respect to a particular norm is

.. math::

   \kappa(A) = \|A\| \|A^{-1}\|,

with the most common choice being the matrix two-norm.

.. cpp:function:: typename Base<F>::type Condition( const Matrix<F>& A, NormType type=TWO_NORM )
.. cpp:function:: typename Base<F>::type Condition( const DistMatrix<F,U,V>& A, NormType type=TWO_NORM )

   Returns the condition number with respect to the specified norm 
   (one, two, or Frobenius).

.. cpp:function:: typename Base<F>::type FrobeniusCondition( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type FrobeniusCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the Frobenius norm.

.. cpp:function:: typename Base<F>::type InfinityCondition( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type InfinityCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the infinity norm.

.. cpp:function:: typename Base<F>::type MaxCondition( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type MaxCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the entrywise maximum norm.

.. cpp:function:: typename Base<F>::type OneCondition( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type OneCondition( const DistMatrix<F,U,V>& A )

   Returns the condition number with respect to the one norm.

.. cpp:function:: typename Base<F>::type TwoCondition( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type TwoCondition( const DistMatrix<F,U,V>& A )

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

   .. cpp:member:: typename Base<F>::type kappa

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

.. cpp:function:: typename Base<F>::type HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type HPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
.. cpp:function:: typename Base<F>::type HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )

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

.. cpp:function:: typename Base<F>::type Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM )

   Assumes that the input is a fully-specified matrix.

.. cpp:function:: typename Base<F>::type HermitianNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type HermitianNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
.. cpp:function:: typename Base<F>::type SymmetricNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )

   Same as :cpp:func:`Norm`, but the matrix is implicitly
   Hermitian/symmetric with the data stored in the triangle specified by
   :cpp:type:`UpperOrLower`.
   Also, while :cpp:func:`Norm` supports every type of distribution,
   :cpp:func:`HermitianNorm`/:cpp:func:`SymmetricNorm` currently only supports 
   the standard matrix distribution.

Alternatively, one may directly call the following routines (note that the entrywise, KyFan, and Schatten norms have an extra parameter and must be called 
directly).

.. cpp:function:: typename Base<F>::type EntrywiseNorm( const Matrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type EntrywiseNorm( const DistMatrix<F,U,V>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type HermitianEntrywiseNorm( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type HermitianEntrywiseNorm( UpperOrLower uplo, const DistMatrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type SymmetricEntrywiseNorm( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type SymmetricEntrywiseNorm( UpperOrLower uplo, const DistMatrix<F>& A, typename Base<F>::type p )

   The :math:`\ell_p` norm of the columns of `A` stacked into a single vector. 
   Note that the Frobenius norm corresponds to the :math:`p=2` case.

.. cpp:function:: typename Base<F>::type EntrywiseOneNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type EntrywiseOneNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The :math:`\ell_1` norm of the columns of `A` stacked into a single vector. 

.. cpp:function:: typename Base<F>::type FrobeniusNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type FrobeniusNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianFrobeniusNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricFrobeniusNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The :math:`\ell_2` norm of the singular values (the Schatten norm with 
   :math:`p=2`), which can be cheaply computed as the :math:`\ell_2` norm of 
   :math:`\text{vec}(A)`.

.. cpp:function:: typename Base<F>::type KyFanNorm( const Matrix<F>& A, int k )
.. cpp:function:: typename Base<F>::type KyFanNorm( const DistMatrix<F,U,V>& A, int k )
.. cpp:function:: typename Base<F>::type HermitianKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, int k )
.. cpp:function:: typename Base<F>::type HermitianKyFanNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, int k )
.. cpp:function:: typename Base<F>::type SymmetricKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, int k )
.. cpp:function:: typename Base<F>::type SymmetricKyFanNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, int k )

   The sum of the largest `k` singular values.

.. cpp:function:: typename Base<F>::type InfinityNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type InfinityNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianInfinityNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricInfinityNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The maximum :math:`\ell_1` norm of the rows of `A`. In the symmetric and 
   Hermitian cases, this is equivalent to the :math:`\|\cdot \|_1` norm.

.. cpp:function:: typename Base<F>::type MaxNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type MaxNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The maximum absolute value of the matrix entries.

.. cpp:function:: typename Base<F>::type NuclearNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type NuclearNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )

   The sum of the singular values. This is equivalent to both the KyFan 
   norm with :math:`k=n` and the Schatten norm with :math:`p=1`.
   Note that the nuclear norm is dual to the two-norm, which is the 
   Schatten norm with :math:`p=\infty`.

.. cpp:function:: typename Base<F>::type OneNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type OneNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricOneNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )

   The maximum :math:`\ell_1` norm of the columns of `A`. In the symmetric and 
   Hermitian cases, this is equivalent to the :math:`\| \cdot \|_\infty` norm.

.. cpp:function:: typename Base<F>::type SchattenNorm( const Matrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type SchattenNorm( const DistMatrix<F,U,V>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type HermitianSchattenNorm( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type HermitianSchattenNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type SymmetricSchattenNorm( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
.. cpp:function:: typename Base<F>::type SymmetricSchattenNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, typename Base<F>::type p )

   The :math:`\ell_p` norm of the singular values.

.. cpp:function:: typename Base<F>::type TwoNorm( const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type TwoNorm( const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type HermitianTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
.. cpp:function:: typename Base<F>::type SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )

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

.. cpp:function:: typename Base<F>::type TwoNormEstimate( Matrix<F>& A, typename Base<F>::type tol=1e-6 )
.. cpp:function:: typename Base<F>::type TwoNormEstimate( DistMatrix<F>& A, typename Base<F>::type tol=1e-6 )
.. cpp:function:: typename Base<F>::type HermitianTwoNormEstimate( Matrix<F>& A, typename Base<F>::type tol=1e-6 )
.. cpp:function:: typename Base<F>::type HermitianTwoNormEstimate( DistMatrix<F>& A, typename Base<F>::type tol=1e-6 )
.. cpp:function:: typename Base<F>::type SymmetricTwoNormEstimate( Matrix<F>& A, typename Base<F>::type tol=1e-6 )
.. cpp:function:: typename Base<F>::type SymmetricTwoNormEstimate( DistMatrix<F>& A, typename Base<F>::type tol=1e-6 )

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

HilbertSchmidt
--------------
The Hilbert-Schmidt inner-product of two :math:`m \times n` matrices :math:`A`
and :math:`B` is :math:`\mbox{tr}(A^H B)`.

.. cpp:function:: F HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B )
.. cpp:function:: F HilbertSchmidt( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
