LAPACK
======
A handful of LAPACK routines are currently used by Elemental: a few
routines for querying floating point characteristics, serial Cholesky and LU 
factorization, triangular inversion, and a few other utilities. In addition, 
there are several BLAS-like routines which are technically part of LAPACK 
(e.g., ``csyr``) which were included in the BLAS imports section.

The prototypes can be found in
`include/elemental/imports/lapack.hpp <../../../../include/elemental/imports/lapack.hpp>`_,
while the implementations are in
`src/imports/lapack.cpp <../../../../src/imports/lapack.cpp>`_.

Machine information
-------------------

In all of the following functions, `R` can be equal to either `float` or
`double`.

.. cpp:function:: R lapack::MachineEpsilon<R>()

   Return the relative machine precision.

.. cpp:function:: R lapack::MachineSafeMin<R>()

   Return the minimum number which can be inverted without underflow.

.. cpp:function:: R lapack::MachinePrecision<R>()

   Return the relative machine precision multiplied by the base.

.. cpp:function:: R lapack::MachineUnderflowExponent<R>()

   Return the minimum exponent before (gradual) underflow occurs.

.. cpp:function:: R lapack::MachineUnderflowThreshold<R>()

   Return the underflow threshold: ``(base)^((underflow exponent)-1)``.

.. cpp:function:: R lapack::MachineOverflowExponent<R>()

   Return the largest exponent before overflow.
    
.. cpp:function:: R lapack::MachineOverflowThreshold<R>()

   Return the overflow threshold: 
   ``(1-rel. prec.)) * (base)^(overflow exponent)``.

Safe norms
----------

.. cpp:function:: R lapack::SafeNorm( R alpha, R beta )

   Return :math:`\sqrt{\alpha^2+\beta^2}` in a manner which avoids 
   under/overflow. `R` can be equal to either `float` or `double`.

.. cpp:function:: R lapack::SafeNorm( R alpha, R beta, R gamma )

   Return :math:`\sqrt{\alpha^2+\beta^2+\gamma^2}` in a manner which avoids
   under/overflow. `R` can be equal to either `float` or `double`.

Givens rotations
----------------

Given :math:`\phi, \gamma \in \mathbb{C}^{n \times n}`, carefully compute 
:math:`c \in \mathbb{R}` and :math:`s, \rho \in \mathbb{C}` such that 

.. math::
   :nowrap:

   \[
   \left[\begin{array}{cc}
     c       & s \\
     -\bar s & c \end{array}\right] 
   \left[ \begin{array}{c} \phi \\ \gamma \end{array} \right] = 
   \left[ \begin{array}{c} \rho \\ 0 \end{array} \right],
   \]

where :math:`c^2 + |s|^2 = 1` and the mapping from :math:`(\phi,\gamma) \rightarrow (c,s,\rho)` is "as continuous as possible", in the manner described by 
Kahan and Demmel's "On computing Givens rotations reliably and efficiently".

.. cpp:function:: void lapack::ComputeGivens( R phi, R gamma, R* c, R* s, R* rho )

   Computes a Givens rotation for real :math:`\phi` and :math:`\gamma`.

.. cpp:function:: void lapack::ComputeGivens( C phi, C gamma, R* c, C* s, C* rho )

   Computes a Givens rotation for complex :math:`\phi` and :math:`\gamma`.

Cholesky factorization
----------------------

.. cpp:function:: void lapack::Cholesky( char uplo, int n, const F* A, int lda )

   Perform a Cholesky factorization on :math:`A \in F^{n \times n}`, where 
   :math:`A(i,j)` can be accessed at ``A[i+j*lda]`` and :math:`A` is implicitly
   Hermitian, with the data stored in the lower triangle if `uplo` equals 
   'L', or in the upper triangle if `uplo` equals 'U'.

LU factorization
----------------

.. cpp:function:: void lapack::LU( int m, int n, F* A, int lda, int* p )

   Perform an LU factorization with partial pivoting on 
   :math:`A \in F^{m \times n}`, where :math:`A(i,j)` can be accessed at 
   ``A[i+j*lda]``. On exit, the pivots are stored in the vector `p`, which 
   should be at least as large as ``min(m,n)``.

Triangular inversion
--------------------

.. cpp:function:: void lapack::TriangularInverse( char uplo, char diag, int n, const F* A, int lda )

   Overwrite either the lower or upper triangle of :math:`A \in F^{n \times n}`
   with its inverse. Which triangle is accessed is determined by `uplo` ('L' for lower or 'U' for upper), and setting `diag` equal to 'U' results in the 
   triangular matrix being treated as unit diagonal (set `diag` to 'N' 
   otherwise).

Hegst
-----

.. cpp:function:: void lapack::Hegst( int itype, char uplo, int n, F* A, int lda, const F* B, int ldb )

   Reduce a well-conditioned generalized Hermitian-definite eigenvalue problem 
   to Hermitian standard form. **TODO:** Explain in more detail.


