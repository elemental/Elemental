Imported library routines
*************************
Since one of the goals of Elemental is to provide high-performance 
datatype-independent parallel routines, yet Elemental's dependencies are 
datatype-dependent, it is convenient to first build a thin datatype-independent
abstraction on top of the necessary routines from BLAS, LAPACK, and MPI. The 
"first-class" datatypes are ``float``, ``double``, ``std::complex<float>``, and 
``std::complex<double>``, but ``int`` and ``byte`` (``unsigned char``) are 
supported for many cases, and support for higher precision arithmetic is in the
works.

BLAS
====
The Basic Linear Algebra Subprograms (BLAS) are heavily exploited within 
Elemental in order to achieve high performance whenever possible. Since the 
official BLAS interface uses different routine names for different datatypes, 
the following interfaces are built directly on top of the datatype-specific 
versions.

Level 1
-------

.. cpp:function:: void blas::Axpy( int n, T alpha, const T* x, int incx, T* y, int incy )

   Performs :math:`y := \alpha x + y` for vectors :math:`x,y \in T^n` and 
   scalar :math:`\alpha \in T`. ``x`` and ``y`` must be stored such that 
   :math:`x_i` occurs at ``x[i*incx]`` (and likewise for ``y``).

.. cpp:function:: T blas::Dot( int n, const T* x, int incx, T* y, int incy )

   Returns :math:`\alpha := x^H y`, where ``x`` and ``y`` are stored in the 
   same manner as in ``blas::Axpy``.

.. cpp:function:: T blas::Dotc( int n, const T* x, int incx, T* y, int incy )

   Equivalent to ``blas::Dot``, but this name is kept for historical purposes
   (the BLAS provide ``?dotc`` and ``?dotu`` for complex datatypes).

.. cpp:function:: T blas::Dotu( int n, const T* x, int incx, T* y, int incy )

   Similar to ``blas::Dot``, but this routine instead returns 
   :math:`\alpha := x^T y` (``x`` is not conjugated).

.. cpp:function:: RealBase<T>::type blas::Nrm2( int n, const T* x, int incx )

   Return the Euclidean two-norm of the vector ``x``, where
   :math:`||x||_2 = \sqrt{\sum_{i=0}^{n-1} |x_i|^2}`. Note that if ``T`` 
   represents a complex field, then the return type is the underlying real field
   (e.g., ``T=std::complex<double>`` results in a return type of ``double``), 
   otherwise ``T`` equals the return type.

.. cpp:function:: void blas::Scal( int n, T alpha, T* x, int incx )

   Performs :math:`x := \alpha x`, where :math:`x \in T^n` is stored in the 
   manner described in ``blas::Axpy``, and :math:`\alpha \in T`.

Level 2
-------

.. cpp:function:: void blas::Gemv( char trans, int m, int n, T alpha, const T* A, int lda, const T* x, int incx, T beta, T* y, int incy )

   Updates :math:`y := \alpha \mbox{op}(A) x + \beta y`, where 
   :math:`A \in T^{m \times n}` and 
   :math:`\mbox{op}(A) \in \left\{A,A^T,A^H\right\}` is chosen by choosing 
   ``trans`` from :math:`\{N,T,C\}`, respectively. Note that ``x`` is stored
   in the manner repeatedly described in the Level 1 routines, e.g., 
   ``blas::Axpy``, but ``A`` is stored such that :math:`A(i,j)` is located
   at ``A[i+j*lda]``.

.. cpp:function:: void blas::Ger( int m, int n, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Updates :math:`A := \alpha x y^H + A`, where :math:`A \in T^{m \times n}` and
   ``x``, ``y``, and ``A`` are stored in the manner described in ``blas::Gemv``.

.. cpp:function:: void blas::Gerc( int m, int n, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Equivalent to ``blas::Ger``, but the name is provided for historical 
   reasons (the BLAS provides ``?gerc`` and ``?geru`` for complex datatypes).

.. cpp:function:: void blas::Geru( int m, int n, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Same as ``blas::Ger``, but instead perform :math:`A := \alpha x y^T + A` 
   (``y`` is not conjugated).

.. cpp:function:: void blas::Hemv( char uplo, int m, T alpha, const T* A, int lda, const T* x, int incx, T beta, T* y, int incy )

   Performs :math:`y := \alpha A x + \beta y`, where 
   :math:`A \in T^{m \times n}` is assumed to be Hermitian with the data stored
   in either the lower or upper triangle of ``A`` (depending upon whether 
   ``uplo`` is equal to 'L' or 'U', respectively).

.. cpp:function:: void blas::Her( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda )

   Performs :math:`A := \alpha x x^H + A`, where :math:`A \in T^{m \times m}` 
   is assumed to be Hermitian, with the data stored in the triangle specified
   by ``uplo`` (depending upon whether ``uplo`` is equal to 'L' or 'U', 
   respectively).

.. cpp:function:: void blas::Her2( char uplo, int m, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Performs :math:`A := \alpha ( x y^H + y x^H ) + A`, where
   :math:`A \in T^{m \times m}` is assumed to be Hermitian, with the data 
   stored in the triangle specified by ``uplo`` (depending upon whether ``uplo``
   is equal to 'L' or 'U', respectively).

.. cpp:function:: void blas::Symv( char uplo, int m, T alpha, const T* A, int lda, const T* x, int incx, T beta, T* y, int incy )

   The same as ``blas::Hemv``, but :math:`A \in T^{m \times m}` is instead 
   assumed to be *symmetric*, and the update is 
   :math:`y := \alpha A x + \beta y`.

.. cpp:function:: void blas::Syr( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda )

   The same as ``blas::Her``, but :math:`A \in T^{m \times m}` is instead 
   assumed to be *symmetric*, and the update is :math:`A := \alpha x x^T + A`.

.. cpp:function:: void blas::Syr2( char uplo, int m, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   The same as ``blas::Her2``, but :math:`A \in T^{m \times m}` is instead
   assumed to be *symmetric*, and the update is 
   :math:`A := \alpha ( x y^T + y x^T ) + A`.

.. cpp:function:: void blas::Trmv( char uplo, char trans, char diag, int m, const T* A, int lda, T* x, int incx )

   Perform the update :math:`x := \alpha \mbox{op}(A) x`, 
   where :math:`A \in T^{m \times m}` is assumed to be either lower or upper
   triangular (depending on whether ``uplo`` is 'L' or 'U'), unit diagonal if 
   ``diag`` equals 'U', and :math:`\mbox{op}(A) \in \left\{A,A^T,A^H\right\}` 
   is determined by ``trans`` being chosen as 'N', 'T', or 'C', respectively.

.. cpp:function:: void blas::Trsv( char uplo, char trans, char diag, int m, const T* A, int lda, T* x, int incx )

   Perform the update :math:`x := \alpha \mbox{op}(A)^{-1} x`, 
   where :math:`A \in T^{m \times m}` is assumed to be either lower or upper
   triangular (depending on whether ``uplo`` is 'L' or 'U'), unit diagonal if 
   ``diag`` equals 'U', and :math:`\mbox{op}(A) \in \left\{A,A^T,A^H\right\}` 
   is determined by ``trans`` being chosen as 'N', 'T', or 'C', respectively.

Level 3
-------

..  cpp:function:: void blas::Gemm( char transA, char transB, int m, int n, int k, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

    Perform the update 
    :math:`C := \alpha \mbox{op}_A(A) \mbox{op}_B(B) + \beta C`, 
    where :math:`\mbox{op}_A` and :math:`\mbox{op}_B` are each determined 
    (according to ``transA`` and ``transB``) in the manner described for 
    ``blas::Trmv``; it is required that :math:`C \in T^{m \times n}` and that
    the inner dimension of :math:`\mbox{op}_A(A) \mbox{op}_B(B)` is ``k``.

.. cpp:function:: void blas::Hemm( char side, char uplo, int m, int n, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

    Perform either :math:`C := \alpha A B + \beta C` or 
    :math:`C := \alpha B A + \beta C` 
    (depending upon whether ``side`` is respectively 'L' or 'R') where 
    :math:`A` is assumed to be Hermitian with its data stored in either the
    lower or upper triangle (depending upon whether ``uplo`` is set to 'L' or 
    'U', respectively) and :math:`C \in T^{m \times n}`.

.. cpp:function:: void blas::Her2k( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha ( A B^H + B A^H ) \beta C` or 
   :math:`C := \alpha ( A^H B + B^H A ) \beta C` (depending upon whether 
   ``trans`` is respectively 'N' or 'C'), where :math:`C \in T^{n \times n}` 
   is assumed to be Hermitian, with the data stored in the triangle specified 
   by ``uplo`` (see ``blas::Hemv``) and the inner dimension of :math:`A B^H` or 
   :math:`A^H B` is equal to ``k``.

.. cpp:function:: void blas::Herk( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha A A^H + \beta C` or 
   :math:`C := \alpha A^H A + \beta C` (depending upon whether ``trans`` is 
   respectively 'N' or 'C'), where :math:`C \in T^{n \times n}` is assumed to
   be Hermitian with the data stored in the triangle specified by ``uplo``
   (see ``blas::Hemv``) and the inner dimension of :math:`A A^H` or 
   :math:`A^H A` equal to ``k``.

.. cpp:function:: void blas::Hetrmm( char uplo, int n, T* A, int lda )

   Form either :math:`A := L^H L` or :math:`A := U U^H`, depending upon the 
   choice of ``uplo``: if ``uplo`` equals 'L', then :math:`L \in T^{n \times n}`
   is equal to the lower triangle of ``A``, otherwise :math:`U` is read from 
   the upper triangle of ``A``. In both cases, the relevant triangle of ``A`` 
   is overwritten in order to store the Hermitian product.

.. cpp:function:: void blas::Symm( char side, char uplo, int m, int n, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

    Perform either :math:`C := \alpha A B + \beta C` or
    :math:`C := \alpha B A + \beta C`
    (depending upon whether ``side`` is respectively 'L' or 'R') where
    :math:`A` is assumed to be symmetric with its data stored in either the
    lower or upper triangle (depending upon whether ``uplo`` is set to 'L' or
    'U', respectively) and :math:`C \in T^{m \times n}`.

.. cpp:function:: void blas::Syr2k( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha ( A B^T + B A^T ) \beta C` or
   :math:`C := \alpha ( A^T B + B^T A ) \beta C` (depending upon whether
   ``trans`` is respectively 'N' or 'T'), where :math:`C \in T^{n \times n}`
   is assumed to be symmetric, with the data stored in the triangle specified
   by ``uplo`` (see ``blas::Symv``) and the inner dimension of :math:`A B^T` or
   :math:`A^T B` is equal to ``k``.

.. cpp:function:: void blas::Syrk( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha A A^T + \beta C` or
   :math:`C := \alpha A^T A + \beta C` (depending upon whether ``trans`` is
   respectively 'N' or 'T'), where :math:`C \in T^{n \times n}` is assumed to
   be symmetric with the data stored in the triangle specified by ``uplo``
   (see ``blas::Symv``) and the inner dimension of :math:`A A^T` or
   :math:`A^T A` equal to ``k``.

.. cpp:function:: void blas::Trmm( char side, char uplo, char trans, char unit, int m, int n, T alpha, const T* A, int lda, T* B, int ldb )

   Performs :math:`C := \alpha \mbox{op}(A) B` or 
   :math:`C := \alpha B \mbox{op}(A)`, depending upon whether ``side`` was 
   chosen as 'L' or 'R', respectively. Whether :math:`A` is treated as lower 
   or upper triangular is determined by whether ``uplo`` is 'L' or 'U' (setting
   ``unit`` equal to 'U' treats :math:`A` as unit diagonal, otherwise it should
   be set to 'N'). :math:`\mbox{op}` is determined in the same manner as in 
   ``blas::Trmv``.

.. cpp:function:: void blas::Trsm( char side, char uplo, char trans, char unit, int m, int n, T alpha, const T* A, int lda, T* B, int ldb )

   Performs :math:`C := \alpha \mbox{op}(A)^{-1} B` or 
   :math:`C := \alpha B \mbox{op}(A)^{-1}`, depending upon whether ``side`` was 
   chosen as 'L' or 'R', respectively. Whether :math:`A` is treated as lower 
   or upper triangular is determined by whether ``uplo`` is 'L' or 'U' (setting
   ``unit`` equal to 'U' treats :math:`A` as unit diagonal, otherwise it should
   be set to 'N'). :math:`\mbox{op}` is determined in the same manner as in 
   ``blas::Trmv``.


LAPACK
======
Only a handful of LAPACK routines are currently used by Elemental: a few
routines for querying floating point characteristics, serial Cholesky and LU 
factorization kernels, and a few random utilities.

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

.. cpp:function:: void lapack::TriangularInverse( char uplo, char diag, int n, const F* A, int lda )

   Overwrite either the lower or upper triangle of :math:`A \in F^{n \times n}`
   with its inverse. Which triangle is accessed is determined by ``uplo`` ('L' for lower or 'U' for upper), and setting ``diag`` equal to 'U' results in the 
   triangular matrix being treated as unit diagonal (set ``diag`` to 'N' 
   otherwise).

MPI
===
All communication within Elemental is built on top of the Message Passing 
Interface (MPI). Just like with BLAS and LAPACK, a minimal set of datatype 
independent abstractions has been built directly on top of the standard 
MPI interface. This has the added benefit of localizing the changes required
for porting Elemental to architectures that do not have full MPI 
implementations available.

Datatypes
---------

.. cpp:type:: mpi::Comm

   Equivalent to ``MPI_Comm``.

.. cpp:type:: mpi::Datatype

   Equivalent to ``MPI_Datatype``.

.. cpp:type:: mpi::ErrorHandler

   Equivalent to ``MPI_Errhandler``.

.. cpp:type:: mpi::Group

   Equivalent to ``MPI_Group``.

.. cpp:type:: mpi::Op

   Equivalent to ``MPI_Op``.

.. cpp:type:: mpi::Request

   Equivalent to ``MPI_Request``.

.. cpp:type:: mpi::Status

   Equivalent to ``MPI_Status``.

.. cpp:type:: mpi::UserFunction

   Equivalent to ``MPI_User_function``.

Constants
---------

.. cpp:member:: const int mpi::ANY_SOURCE

   Equivalent to ``MPI_ANY_SOURCE``.

.. cpp:member:: const int mpi::ANY_TAG

   Equivalent to ``MPI_ANY_TAG``.

.. cpp:member:: const int mpi::THREAD_SINGLE

   Equivalent to ``MPI_THREAD_SINGLE``.

.. cpp:member:: const int mpi::THREAD_FUNNELED

   Equivalent to ``MPI_THREAD_FUNNELED``.

.. cpp:member:: const int mpi::THREAD_SERIALIZED

   Equivalent to ``MPI_THREAD_SERIALIZED``.

.. cpp:member:: const int mpi::THREAD_MULTIPLE

   Equivalent to ``MPI_THREAD_MULTIPLE``.

.. cpp:member:: const int mpi::UNDEFINED

   Equivalent to ``MPI_UNDEFINED``.

.. cpp:member:: const mpi::Comm mpi::COMM_WORLD

   Equivalent to ``MPI_COMM_WORLD``.

.. cpp:member:: const mpi::ErrorHandler mpi::ERRORS_RETURN
   
   Equivalent to ``MPI_ERRORS_RETURN``.

.. cpp:member:: const mpi::ErrorHandler mpi::ERRORS_ARE_FATAL

   Equivalent to ``MPI_ERRORS_ARE_FATAL``.

.. cpp:member:: const mpi::Group mpi::GROUP_EMPTY

   Equivalent to ``MPI_GROUP_EMPTY``.

.. cpp:member:: const mpi::Request mpi::REQUEST_NULL

   Equivalent to ``MPI_REQUEST_NULL``.

.. cpp:member:: const mpi::Op mpi::MAX

   Equivalent to ``MPI_MAX``.

.. cpp:member:: const mpi::Op mpi::PROD

   Equivalent to ``MPI_PROD``.

.. cpp:member:: const mpi::Op mpi::SUM

   Equivalent to ``MPI_SUM``.

.. cpp:member:: const int mpi::MIN_COLL_MSG

   The minimum message size for collective communication, e.g., the minimum
   number of elements contributed by each process in an ``MPI_Allgather``. 
   By default, it is hardcoded to ``1`` in order to avoid problems with 
   MPI implementations that do not support the ``0`` corner case.

Routines
--------

.. rubric:: Environmental

.. cpp:function:: void mpi::Initialize( int& argc, char**& argv )

   Equivalent of ``MPI_Init`` 
   (but notice the difference in the calling convention).

   .. code-block:: cpp

      #include "elemental.hpp"
      using namespace elemental;

      int main( int argc, char* argv[] )
      {
          mpi::Initialize( argc, argv );
          ...
          mpi::Finalize();
          return 0;
      }

.. cpp:function:: int mpi::InitializeThread( int& argc, char**& argv, int required )

   The threaded equivalent of ``mpi::Initialize``; the return integer indicates
   the level of achieved threading support, e.g., ``mpi::THREAD_MULTIPLE``.

.. cpp:function:: void mpi::Finalize()

   Shut down the MPI environment, freeing all of the allocated resources.

.. cpp:function:: bool mpi::Initialized()

   Return whether or not MPI has been initialized.

.. cpp:function:: bool mpi::Finalized()

   Return whether or not MPI has been finalized.

.. cpp:function:: double mpi::Time()

   Return the current wall-time in seconds.

.. cpp:function:: void mpi::OpCreate( mpi::UserFunction* func, bool commutes, Op& op )

   Create a custom operation for use in reduction routines, e.g., 
   ``mpi::Reduce``, ``mpi::AllReduce``, and ``mpi::ReduceScatter``, where
   ``mpi::UserFunction`` could be defined as

   .. code-block:: cpp

      namespace mpi {
      typedef void (UserFunction) ( void* a, void* b, int* length, mpi::Datatype* datatype );
      }

   The ``commutes`` parameter is also important, as it specifies whether or not
   the operation ``b[i] = a[i] op b[i], for i=0,...,length-1``, can be 
   performed in an arbitrary order (for example, using a minimum spanning tree).

.. cpp:function:: void mpi::OpFree( mpi::Op& op )

   Free the specified MPI reduction operator.

.. rubric:: Communicator manipulation

.. cpp:function:: int mpi::CommRank( mpi::Comm comm )

   Return our rank in the specified communicator.

.. cpp:function:: int mpi::CommSize( mpi::Comm comm )

   Return the number of processes in the specified communicator.

.. cpp:function:: void mpi::CommCreate( mpi::Comm parentComm, mpi::Group subsetGroup, mpi::Comm& subsetComm )

   Create a communicator (``subsetComm``) which is a subset of ``parentComm`` 
   consisting of the processes specified by ``subsetGroup``.

.. cpp:function:: void mpi::CommDup( mpi::Comm original, mpi::Comm& duplicate )

   Create a copy of a communicator.

.. cpp:function:: void mpi::CommSplit( mpi::Comm comm, int color, int key, mpi::Comm& newComm )

   Split the communicator ``comm`` into different subcommunicators, where each 
   process specifies the ``color`` (unique integer) of the subcommunicator it 
   will reside in, as well as its ``key`` (rank) for the new subcommunicator.

.. cpp:function:: void mpi::CommFree( mpi::Comm& comm )

   Free the specified communicator.

.. cpp:function:: bool mpi::CongruentComms( mpi::Comm comm1, mpi::Comm comm2 )

   Return whether or not the two communicators consist of the same set of 
   processes (in the same order).

.. cpp:function:: void mpi::ErrorHandlerSet( mpi::Comm comm, mpi::ErrorHandler errorHandler )

   Modify the specified communicator to use the specified error-handling 
   approach.

.. rubric:: Cartesian communicator manipulation

.. cpp:function:: void mpi::CartCreate( mpi::Comm comm, int numDims, const int* dimensions, const int* periods, bool reorder, mpi::Comm& cartComm )

   Create a Cartesian communicator (``cartComm``) from the specified 
   communicator (``comm``), given the number of dimensions (``numDims``), 
   the sizes of each dimension (``dimensions``), whether or not each 
   dimension is periodic (``periods``), and whether or not the ordering of the 
   processes may be changed (``reorder``).

.. cpp:function:: void mpi::CartSub( mpi::Comm comm, const int* remainingDims, mpi::Comm& subComm )

   Create this process's subcommunicator of ``comm`` that results from only 
   keeping the specified dimensions (``0`` for ignoring and ``1`` for keeping).

.. rubric:: Group manipulation

.. cpp:function:: int mpi::GroupRank( mpi::Group group )

   Return our rank in the specified group.

.. cpp:function:: int mpi::GroupSize( mpi::Group group )

   Return the number of processes in the specified group.

.. cpp:function:: void mpi::CommGroup( mpi::Comm comm, mpi::Group& group )

   Extract the underlying group from the specified communicator.

.. cpp:function:: void mpi::GroupIncl( mpi::Group group, int n, const int* ranks, mpi::Group& subGroup )

   Create a subgroup of ``group`` that consists of the ``n`` processes whose 
   ranks are specified in the ``ranks`` array.

.. cpp:function:: void mpi::GroupDifference( mpi::Group parent, mpi::Group subset, mpi::Group& complement )

   Form a group (``complement``) out of the set of processes which are in 
   the ``parent`` communicator, but not in the ``subset`` communicator.

.. cpp:function:: void mpi::GroupFree( mpI::Group& group )

   Free the specified group.

.. cpp:function:: void mpi::GroupTranslateRanks( mpi::Group origGroup, int size, const int* origRanks, mpi::Group newGroup, int* newRanks )

   Return the ranks within ``newGroup`` of the ``size`` processes specified 
   by their ranks in the ``origGroup`` communicator using the ``origRanks`` 
   array. The result will be in the ``newRanks`` array, which must have been 
   preallocated to a length at least as large as ``size``.

.. rubric:: Utilities

.. cpp:function:: void mpi::Barrier( mpi::Comm comm )

   Pause until all processes within the ``comm`` communicator have called this
   routine.

.. cpp:function:: void mpi::Wait( mpi::Request& request )

   Pause until the specified request has completed.

.. cpp:function:: bool mpi::Test( mpi::Request& request )

   Return whether or not the specified request has completed.

.. cpp:function:: bool mpi::IProbe( int source, int tag, mpi::Comm comm, mpi::Status& status )

   Return whether or not there is a message ready which

   * is from the process with rank ``source`` in the communicator ``comm``
     (note that ``mpi::ANY_SOURCE`` is allowed)
   * had the integer tag ``tag``

   If ``true`` was returned, then ``status`` will have been filled with the 
   relevant information, e.g., the source's rank.

.. cpp:function:: int mpi::GetCount<T>( mpi::Status& status )

   Return the number of entries of the specified datatype which are ready to 
   be received.

.. rubric:: Point-to-point communication

.. cpp:function:: void mpi::Send( const T* buf, int count, int to, int tag, mpi::Comm comm )

   Send ``count`` entries of type ``T`` to the process with rank ``to`` in the 
   communicator ``comm``, and tag the message with the integer ``tag``.

.. cpp:function:: void mpi::ISend( const T* buf, int count, int to, int tag, mpi::Comm comm, mpi::Request& request )

   Same as ``mpi::Send``, but the call is non-blocking.

.. cpp:function:: void mpi::ISSend( const T* buf, int count, int to, int tag, mpi::Comm comm, mpi::Request& request )

   Same as ``mpi::ISend``, but the call is in synchronous mode.

.. cpp:function:: void mpi::Recv( T* buf, int count, int from, int tag, mpi::Comm comm )

   Receive ``count`` entries of type ``T`` from the process with rank ``from`` 
   in the communicator ``comm``, where the message must have been tagged with 
   the integer ``tag``.

.. cpp:function:: void mpi::IRecv( T* buf, int count, int from, int tag, mpi::Comm comm, mpi::Request& request )

   Same as ``mpi::Recv``, but the call is non-blocking.

.. cpp:function:: void mpi::SendRecv( const T* sendBuf, int sendCount, int to, int sendTag, T* recvBuf, int recvCount, int from, int recvTag, mpi::Comm comm )

   Send ``sendCount`` entries of type ``T`` to process ``to``, and 
   simultaneously receive ``recvCount`` entries of type ``T`` from process 
   ``from``.

.. rubric:: Collective communication

.. cpp:function:: void mpi::Broadcast( T* buf, int count, int root, mpi::Comm comm )

   The contents of ``buf`` (``count`` entries of type ``T``) on process ``root``
   are duplicated in the local buffers of every process in the communicator.

.. cpp:function:: void mpi::Gather( const T* sendBuf, int sendCount, T* recvBuf, int recvCount, int root, mpi::Comm comm )

   Each process sends an independent amount of data (i.e., ``sendCount`` 
   entries of type ``T``) to the process with rank ``root``; the ``root`` 
   process must specify the maximum number of entries sent from each process, 
   ``recvCount``, so that the data received from process ``i`` lies within the 
   ``[i*recvCount,(i+1)*recvCount)`` range of the receive buffer.

.. cpp:function:: void mpi::AllGather( const T* sendBuf, int sendCount, T* recvBuf, int recvCount, mpi::Comm comm )

   Same as ``mpi::Gather``, but every process receives the result.

.. cpp:function:: void mpi::Scatter( const T* sendBuf, int sendCount, T* recvBuf, int recvCount, int root, mpi::Comm comm )

   The same as ``mpi::Gather``, but in reverse: the root process starts with 
   an array of data and sends the ``[i*sendCount,(i+1)*sendCount)`` entries 
   to process ``i``. 

.. cpp:function:: void mpi::AllToAll( const T* sendBuf, int sendCount, T* recvBuf, int recvCount, mpi::Comm comm )

   This can be thought of as every process simultaneously scattering data: after
   completion, the ``[i*recvCount,(i+1)*recvCount)`` portion of the receive 
   buffer on process ``j`` will contain the ``[j*sendCount,(j+1)*sendCount)`` 
   portion of the send buffer on process ``i``, where ``sendCount`` refers to 
   the value specified on process ``i``, and ``recvCount`` refers to the value
   specified on process ``j``.

.. cpp:function:: void mpi::AllToAll( const T* sendBuf, const int* sendCounts, const int* sendDispls, T* recvBuf, const int* recvCounts, const int* recvDispls, mpi::Comm comm )

   Same as previous ``mpi::AllToAll``, but the amount of data sent to and 
   received from each process is allowed to vary; after completion, the 
   ``[recvDispls[i],recvDispls[i]+recvCounts[i])`` portion of the receive buffer
   on process ``j`` will contain the 
   ``[sendDispls[j],sendDispls[j]+sendCounts[j])`` portion of the send buffer
   on process ``i``.

.. cpp:function:: void mpi::Reduce( const T* sendBuf, T* recvBuf, int count, mpi::Op op, int root, mpi::Comm comm )

   The ``root`` process receives the result of performing 

   :math:`S_{p-1} + (S_{n-2} + \cdots (S_2 + (S_1 + S_0)) \cdots )`,
   where :math:`S_i` represents the send buffer of process ``i``, and :math:`+`
   represents the operation specified by ``op``.

.. cpp:function:: void mpi::AllReduce( const T* sendBuf, T* recvBuf, int count, mpi::Op op, mpi::Comm comm )

   Same as ``mpi::Reduce``, but every process receives the result.

.. cpp:function:: void mpi::ReduceScatter( const T* sendBuf, T* recvBuf, const int* recvCounts, mpi::Op op, mpi::Comm comm )

   Same as ``mpi::AllReduce``, but process ``0`` only receives the 
   ``[0,recvCounts[0])`` portion of the result, process ``1`` only receives the 
   ``[recvCounts[0],recvCounts[0]+recvCounts[1])`` portion of the result, 
   etc.

Parallel LCG
============
Since it is often necessary to generate a large matrix with pseudo-random 
entries in parallel, a method for ensuring that a large set of processes can 
each generate independent uniformly random samples is required. The purpose of
Parallel LCG (PLCG) is to provide a provably independent generalization of a
simple (but well-studied) Linear Congruential Generator. Knuth's constants from
The Art of Computer Programming Vol. 2 are used.

Datatypes
---------

.. cpp:type:: plcg::UInt32

   Since the vast majority of modern systems make use of ``unsigned`` for
   storing 32-bit unsigned integers, we simply hardcode the type. If your 
   system does not follow this convention, then this typedef will need to be
   changed!

.. cpp:type:: struct plcg::UInt64

   A custom 64-bit unsigned integer which is simply the concatenation of two 
   32-bit unsigned integers (``UInt32``).

.. cpp:type:: struct plcg::ExpandedUInt64

   A custom 64-bit unsigned integer which is stores each of the four 16-bit
   pieces within the first 16 bits of a 32-bit unsigned integer. This is done
   so that two such expanded 16-bit numbers can be multiplied without any 
   chance of overflow.

LCG primitives
--------------

.. cpp:function:: plcg::UInt32 plcg::Lower16Bits( plcg::UInt32 a )

   Return the lower 16 bits of ``a`` in the lower 16 bits of the returned 
   32-bit unsigned integer.

.. cpp:function:: plcg::UInt32 plcg::Upper16Bits( plcg::UInt32 a )

   Return the upper 16 bits of ``a`` in the lower 16 bits of the returned
   32-bit unsigned integer.

.. cpp:function:: plcg::ExpandedUInt64 plcg::Expand( plcg::UInt32 a )

   Expand a 32-bit unsigned integer into a 64-bit expanded representation.

.. cpp:function:: plcg::ExpandedUInt64 plcg::Expand( plcg::UInt64 a )

   Expand a 64-bit unsigned integer into a 64-bit expanded representation.

.. cpp:function:: plcg::UInt64 plcg::Deflate( plcg::ExpandedUInt64 a )

   Deflate an expanded 64-bit unsigned integer into the standard 64-bit form.

.. cpp:function:: void plcg::CarryUpper16Bits( plcg::ExpandedUInt64& a )

   Carry the results stored in the upper 16-bits of each of the four pieces 
   into the next lower 16 bits.

.. cpp:function:: plcg::ExpandedUInt64 plcg::AddWith64BitMod( plcg::ExpandedUInt64 a, plcg::ExpandedUInt64 b )

   Return :math:`a+b \mod 2^{64}`.

.. cpp:function:: plcg::ExpandedUInt64 plcg::MultiplyWith64BitMod( plcg::ExpandedUInt64 a, plcg::ExpandedUInt64 b )

   Return :math:`ab \mod 2^{64}`.

.. cpp:function:: plcg::ExpandedUInt64 plcg::IntegerPowerWith64BitMod( plcg::ExpandedUInt64 x, plcg::ExpandedUInt64 n )

   Return :math:`x^n \mod 2^{64}`.

.. cpp:function:: void plcg::Halve( plcg::ExpandedUInt64& a )

   :math:`a := a/2`.

.. cpp:function:: void plcg::SeedSerialLcg( plcg::UInt64 globalSeed )

   Set the initial state of the serial Linear Congruential Generator.

.. cpp:function:: void plcg::SeedParallelLcg( plcg::UInt32 rank, plcg::UInt32 commSize, plcg::UInt64 globalSeed )

   Have our process seed a separate LCG meant for parallel computation, where 
   the calling process has the given rank within a communicator of the 
   specified size.

.. cpp:function:: plcg::UInt64 plcg::SerialLcg()

   Return the current state of the serial LCG, and then advance to the next one.

.. cpp:function:: plcg::UInt64 plcg::ParallelLcg()

   Return the current state of our process's portion of the parallel LCG, 
   and then advance to our next local state.

.. cpp:function:: void plcg::ManualLcg( plcg::ExpandedUInt64 a, plcg::ExpandedUInt64 c, plcg::ExpandedUInt64& X )

   :math:`X := a X + c \mod 2^{64}`.

Sampling
--------

.. cpp:function:: R plcg::SerialUniform()

   Return a uniform sample from :math:`(0,1]` using the serial LCG.

.. cpp:function:: R plcg::ParallelUniform()

   Return a uniform sample from :math:`(0,1]` using the parallel LCG.

.. cpp:function:: void plcg::SerialBoxMuller( R& X, R& Y )

   Return two samples from a normal distribution with mean 0 and standard 
   deviation of 1 using the serial LCG.

.. cpp:function:: void plcg::ParallelBoxMuller( R& X, R& Y )

   Return two samples from a normal distribution with mean 0 and standard
   deviation 1, but using the parallel LCG.

.. cpp:function:: void plcg::SerialGaussianRandomVariable( R& X )

   Return a single sample from a normal distribution with mean 0 and 
   standard deviation 1 using the serial LCG.

.. cpp:function:: void plcg::ParallelGaussianRandomVariable( R& X )

   Return a single sample from a normal distribution with mean 0 and 
   standard deviation 1, but using the parallel LCG.
   
.. cpp:function:: void plcg::SerialGaussianRandomVariable( std::complex<R>& X )

   Return a single complex sample from a normal distribution with mean 0 and 
   standard deviation 1 using the serial LCG.

.. cpp:function:: void plcg::ParallelGaussianRandomVariable( std::complex<R>& X )

   Return a single complex sample from a normal distribution with mean 0 and 
   standard deviation 1, but using the parallel LCG.

PMRRR
=====
Rather than directly using Petschow and Bientinesi's parallel implementation of 
the Multiple Relatively Robust Representations (MRRR) algorithm, several 
simplified interfaces have been exposed.

Data structures
---------------

.. cpp:type:: struct pmrrr::Estimate

   For returning upper bounds on the number of local and global eigenvalues
   with eigenvalues lying in the specified interval, :math:`(a,b]`.

   .. cpp:member:: int numLocalEigenvalues

      The upper bound on the number of eigenvalues in the specified interval 
      that our process stores locally.

   .. cpp:member:: int numGlobalEigenvalues

      The upper bound on the number of eigenvalues in the specified interval.

.. cpp:type:: struct pmrrr::Info

   For returning information about the computed eigenvalues.

   .. cpp:member:: int numLocalEigenvalues

      The number of computed eigenvalues that our process locally stores.

   .. cpp:member:: int numGlobalEigenvalues

      The number of computed eigenvalues.

   .. cpp:member:: int firstLocalEigenvalue

      The index of the first eigenvalue stored locally on our process.

Compute all eigenvalues
-----------------------

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm )

   Compute all of the eigenvalues of the real symmetric tridiagonal matrix with 
   diagonal ``d`` and subdiagonal ``e``: the eigenvalues will be stored in 
   ``w`` and the work will be divided among the processors in ``comm``.

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, double* Z, int ldz, mpi::Comm comm )

   Same as above, but also compute the corresponding eigenvectors.

Compute eigenvalues within interval
-----------------------------------

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm, double a, double b )

   Only compute the eigenvalues lying within the interval :math:`(a,b]`.

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, double* Z, int ldz, mpi::Comm comm, double a, double b )

   Same as above, but also compute the corresponding eigenvectors.

.. cpp:function:: pmrrr::Estimate pmrrr::EigEstimate( int n, const double* d, const double* w, mpi::Comm comm, double a, double b )

   Return upper bounds on the number of local and global eigenvalues lying 
   within the specified interval.

Compute eigenvalues in index range
----------------------------------

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm, int a, int b )

   Only compute the ``a-b`` eigenvalues of the tridiagonal matrix, where 
   :math:`0 \le a \le b < n`.

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, double* Z, int ldz, mpi::Comm comm, int a, int b )

   Same as above, but also compute the corresponding eigenvectors.
