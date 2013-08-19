LAPACK
------
A handful of LAPACK routines are currently used by Elemental: a few
routines for querying floating point characteristics and a few other utilities.
In addition, there are several BLAS-like routines which are technically part 
of LAPACK (e.g., ``csyr``) which were included in the BLAS imports section.

The prototypes can be found in
`include/elemental/core/imports/lapack.hpp <https://github.com/elemental/Elemental/tree/master/include/elemental/core/imports/lapack.hpp>`_,
while the implementations are in
`src/imports/lapack.cpp <https://github.com/elemental/Elemental/tree/master/src/imports/lapack.cpp>`_.

Machine information
^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^

.. cpp:function:: R lapack::SafeNorm( R alpha, R beta )

   Return :math:`\sqrt{\alpha^2+\beta^2}` in a manner which avoids 
   under/overflow. `R` can be equal to either `float` or `double`.

.. cpp:function:: R lapack::SafeNorm( R alpha, R beta, R gamma )

   Return :math:`\sqrt{\alpha^2+\beta^2+\gamma^2}` in a manner which avoids
   under/overflow. `R` can be equal to either `float` or `double`.

Givens rotations
^^^^^^^^^^^^^^^^

Given :math:`\phi, \gamma \in \mathbb{C}^{n \times n}`, carefully compute 
:math:`c \in \mathbb{R}` and :math:`s, \rho \in \mathbb{C}` such that 

.. math::

   \left[\begin{array}{cc}
     c       & s \\
     -\bar s & c \end{array}\right] 
   \left[ \begin{array}{c} \phi \\ \gamma \end{array} \right] = 
   \left[ \begin{array}{c} \rho \\ 0 \end{array} \right],

where :math:`c^2 + |s|^2 = 1` and the mapping from :math:`(\phi,\gamma) \rightarrow (c,s,\rho)` is "as continuous as possible", in the manner described by 
Kahan and Demmel's "On computing Givens rotations reliably and efficiently".

.. cpp:function:: void lapack::ComputeGivens( F phi, F gamma, typename Base<F>::type* c, F* s, F* rho )

   Computes a Givens rotation.

MRRR-based Hermitian EVP 
^^^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void lapack::HermitianEig( char job, char range, char uplo, int n, F* A, int lda, typename Base<F>::type vl, typename Base<F>::type vu, int il, int iu, typename Base<F>::type abstol, typename Base<F>::type* w, F* Z, int ldz )

   Computes the eigenvalue decomposition of a Hermitian matrix using MRRR.

QR- and DQDS-based SVD
^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void lapack::QRSVD( int m, int n, F* A, int lda, typename Base<F>::type* s, F* U, int ldu, F* VAdj, int ldva )

   Computes the singular value decomposition of a general matrix by running the 
   QR algorithm on the condensed bidiagonal matrix.

.. cpp:function:: void lapack::SVD( int m, int n, F* A, int lda, typename Base<F>::type* s )

   Computes the singular values of a general matrix by running DQDS on the 
   condensed bidiagonal matrix.

Divide-and-conquer SVD
^^^^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void lapack::DivideAndConquerSVD( int m, int n, F* A, int lda, typename Base<F>::type* s, F* U, int ldu, F* VAdj, int ldva )

   Computes the SVD of a general matrix using a divide-and-conquer algorithm on
   the condensed bidiagonal matrix.

Bidiagonal QR
^^^^^^^^^^^^^

.. cpp:function:: void lapack::BidiagQRAlg( char uplo, int n, int numColsVTrans, int numRowsU, typename Base<F>::type* d, typename Base<F>::type* e, F* VAdj, int ldva, F* U, int ldu )

   Computes the SVD of a bidiagonal matrix using the QR algorithm.

Hessenberg QR
^^^^^^^^^^^^^

.. cpp:function:: void lapack::HessenbergEig( int n, F* H, int ldh, Complex<typename Base<F>::type>* w )

   Computes the eigenvalues of an upper Hessenberg matrix using the QR 
   algorithm.

Schur decomposition
^^^^^^^^^^^^^^^^^^^

.. cpp:function:: void lapack::Eig( int n, F* A, int lda, Complex<typename Base<F>::type>* w, bool fullTriangle=false )

   Returns the eigenvalues of a square matrix using the QR algorithm.

.. cpp:function:: void lapack::Schur( int n, F* A, int lda, F* Q, int ldq, Complex<typename Base<F>::type>* w, bool fullTriangle=true )

   Returns the Schur decomposition of a square matrix using the QR algorithm.
