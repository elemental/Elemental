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

.. cpp:function:: F Determinant( Matrix<F>& A )

   Returns the determinant of the (fully populated) square matrix `A`, which is 
   overwritten during the computation.

.. cpp:function:: F Determinant( DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.

.. cpp:type:: struct SafeProduct<F>

   Represents the product of `n` values as :math:`\rho \exp(\kappa n)`, 
   where :math:`|\rho| \le 1` and :math:`\kappa \in \mathbb{R}`.

   .. cpp:member:: F rho

      For nonzero values, `rho` is the modulus and lies *on* the unit 
      circle; in order to represent a value that is precisely zero, `rho` 
      is set to zero.

   .. cpp:member:: typename RealBase<F>::type kappa

      `kappa` can be an arbitrary real number.

   .. cpp:member:: int n

      The number of values in the product.

.. cpp:function:: SafeProduct<F> SafeDeterminant( Matrix<F>& A )

   Returns the determinant of the square matrix `A` in an expanded form 
   which is less likely to over/under-flow.

.. cpp:function:: SafeProduct<F> SafeDeterminant( DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.

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

.. cpp:function:: F Trace( const Matrix<F>& A )

   Return the trace of the square matrix `A`.

.. cpp:function:: F Trace( const DistMatrix<F,MC,MR>& A )

   Same as above, but for a distributed matrix.

