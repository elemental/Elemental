LAPACK
======
A handful of LAPACK routines are currently used by Elemental: a few
routines for querying floating point characteristics, serial Cholesky and LU 
factorization, triangular inversion, and a few other utilities. In addition, 
there are several BLAS-like routines which are technically part of LAPACK 
(e.g., ``csyr``) which were included in the BLAS imports section.

Machine information
-------------------

In all of the following functions, ``R`` can be equal to either ``float`` or
``double``.

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

Factorizations
--------------

.. cpp:function:: void lapack::Cholesky( char uplo, int n, const F* A, int lda )

   Perform a Cholesky factorization on :math:`A \in F^{n \times n}`, where 
   :math:`A(i,j)` can be accessed at ``A[i+j*lda]`` and :math:`A` is implicitly
   Hermitian, with the data stored in the lower triangle if ``uplo`` equals 
   'L', or in the upper triangle if ``uplo`` equals 'U'.

.. cpp:function:: void lapack::LU( int m, int n, F* A, int lda, int* p )

   Perform an LU factorization with partial pivoting on 
   :math:`A \in F^{m \times n}`, where :math:`A(i,j)` can be accessed at 
   ``A[i+j*lda]``. On exit, the pivots are stored in the vector ``p``, which 
   should be at least as large as ``min(m,n)``.

Inversion
---------

.. cpp:function:: void lapack::TriangularInverse( char uplo, char diag, int n, const F* A, int lda )

   Overwrite either the lower or upper triangle of :math:`A \in F^{n \times n}`
   with its inverse. Which triangle is accessed is determined by ``uplo`` ('L' for lower or 'U' for upper), and setting ``diag`` equal to 'U' results in the 
   triangular matrix being treated as unit diagonal (set ``diag`` to 'N' 
   otherwise).

Utilities
---------

.. cpp:function:: void lapack::Hegst( int itype, char uplo, int n, F* A, int lda, const F* B, int ldb )

   Reduce a generalized Hermitian-definite eigenvalue problem to Hermitian 
   standard form. **TODO:** Explain in more detail.

.. cpp:function:: R lapack::SafeNorm( R alpha, R beta )

   Return :math:`\sqrt{\alpha^2+\beta^2}` in a manner which avoids 
   under/overflow. ``R`` can be equal to either ``float`` or ``double``.

.. cpp:function:: R lapack::SafeNorm( R alpha, R beta, R gamma )

   Return :math:`\sqrt{\alpha^2+\beta^2+\gamma^2}` in a manner which avoids
   under/overflow. ``R`` can be equal to either ``float`` or ``double``.

