Environment
===========

Set up and clean up
-------------------

.. cpp:function:: void Initialize( int& argc, char**& argv )

   Initializes Elemental and (if necessary) MPI. The usage is very similar to 
   ``MPI_Init``, but the ``argc`` and ``argv`` can be directly passed in.

   .. code-block:: cpp

      #include "elemental.hpp"

      int
      main( int argc, char* argv[] )
      {
          elemental::Initialize( argc, argv );

          // ...

          elemental::Finalize();
          return 0;
      }

.. cpp:function:: void Finalize()

   Frees all resources allocated by Elemental and (if necessary) MPI.

.. cpp:function:: bool Initialized()

   Return whether or not Elemental is currently initialized.

Blocksize manipulation
----------------------

.. cpp:function:: int Blocksize()

   Return the currently chosen algorithmic blocksize. The optimal value 
   depends on the problem size, algorithm, and architecture; the default value
   is 128.

.. cpp:function:: void SetBlocksize( int blocksize )

   Change the algorithmic blocksize to the specified value.

.. cpp:function:: void PushBlocksizeStack( int blocksize )

   It is frequently useful to temporarily change the algorithmic blocksize, so 
   rather than having to manually store and reset the current state, one can 
   simply push a new value onto a stack 
   (and later pop the stack to reset the value).

.. cpp:function:: void PopBlocksizeStack() 

   Pops the stack of blocksizes. See above.

Scalar manipulation
-------------------

.. cpp:function:: Z Abs( Z alpha )

   Return the absolute value of non-complex :math:`\alpha`.

.. cpp:function:: Z Abs( std::complex<Z> alpha )

   Return the absolute value of complex :math:`\alpha`: 

   .. math::

      |\alpha| = \sqrt{|\mathcal{R}(\alpha)|^2+|\mathcal{I}(\alpha)|^2}

.. cpp:function:: Z FastAbs( Z alpha )

   Same as :cpp:function::Abs: for non-complex :math:`\alpha`.

.. cpp:function:: Z FastAbs( std::complex<Z> alpha )

   Return a cheaper norm of the complex :math:`\alpha`:

   .. math::
   
      |\alpha|_{\mbox{fast}} = |\mathcal{R}(\alpha)| + |\mathcal{I}(\alpha)|

.. cpp:function:: Z Conj( Z alpha )

   Return the complex conjugate of the non-complex variable :math:`\alpha`,
   which is simply :math:`\alpha`.

.. cpp:function:: std::complex<Z> Conj( std::complex<Z> alpha )

   Return the complex conjugate of the complex variable :math:`\alpha`,

   .. math::

      \bar \alpha = \mathcal{R}(\alpha) - \mathcal{I}(\alpha) i

.. cpp:class:: RealBase<F>

   .. cpp:type:: type

      The underlying real datatype of the (potentially complex) datatype ``F``.
      For example, ``typename RealBase< std::complex<double> >::type`` and 
      ``typename RealBase<double>::type`` are both equivalent to ``double``.
      This is often extremely useful in implementing routines which are 
      templated over real and complex datatypes but still make use of real 
      datatypes.

Custom datatypes
----------------

.. cpp:type:: byte

   ``typedef unsigned char byte;``

.. cpp:type:: scomplex

   ``typedef std::complex<float> scomplex;``

.. cpp:type:: dcomplex

   ``typedef std::complex<double> dcomplex;``

.. cpp:type:: Diagonal

   An enum specifying whether a triangular matrix's diagonal is explicitly stored 
   (``NON_UNIT``) or is implicitly unit diaognal (``UNIT``).

.. cpp:type:: Orientation

   An enum for specifying whether a matrix, say :math:`A`, should be implicitly 
   treated as :math:`A` (``NORMAL``), :math:`A^H` (``ADJOINT``), or :math:`A^T`
   (``TRANSPOSE``).

.. cpp:type:: Side

   An enum for specifying ``LEFT`` or ``RIGHT``.

.. cpp:type:: VectorDirection

   An enum for specifying ``VERTICAL`` or ``HORIZONTAL``.

.. cpp:type:: ForwardOrBackward

   An enum for specifying ``FORWARD`` or ``BACKWARD``.

.. cpp:type:: Distribution

   An enum for specifying the distribution of a row or column of a distributed
   matrix:

   * ``MC``: Column of a matrix distribution
   * ``MD``: Diagonal of a matrix distribution
   * ``MR``: Row of a matrix distribution
   * ``VC``: Column-major vector distribution
   * ``VR``: Row-major vector distribution
   * ``STAR``: Redundantly stored

.. cpp:type:: GridOrder

   An enum for specifying either ``ROW_MAJOR`` or ``COLUMN_MAJOR``.

Custom exceptions
-----------------

.. cpp:class:: SingularMatrixException

   An extension of ``std::runtime_error`` which is meant to be thrown when 
   a singular matrix is unexpectedly encountered.

   .. cpp:function:: SingularMatrixException( const char* msg="Matrix was singular" )

      Builds an instance of the exception which allows one to optionally 
      specify the error message.

   .. code-block:: cpp

      throw elemental::SingularMatrixException();

.. cpp:class:: NonHPDMatrixException 

   An extension of ``std::runtime_error`` which is meant to be thrown when
   a non positive-definite Hermitian matrix is unexpectedly encountered
   (e.g., during Cholesky factorization).

   .. cpp:function:: NonHPDMatrixException( const char* msg="Matrix was not HPD" )

      Builds an instance of the exception which allows one to optionally 
      specify the error message.

   .. code-block:: cpp

      throw elemental::NonHPDMatrixException();

Call stack manipulation
-----------------------

.. note::

   The following call stack manipulation routines are only available in 
   non-release builds (i.e., PureDebug and HybridDebug) and are meant to allow 
   for the call stack to be printed (ala DumpCallStack) when an 
   exception is caught.

.. cpp:function:: void PushCallStack( std::string s )

   Push the given routine name onto the call stack.

.. cpp:function:: void PopCallStack()

   Remove the routine name at the top of the call stack.

.. cpp:function:: void DumpCallStack()

   Print (and empty) the contents of the call stack.
   
Default process grid
--------------------

.. cpp:function:: Grid& DefaultGrid()

   Return a process grid built over :cpp:type:`mpi::COMM_WORLD`. This is 
   typically used as a means of allowing instances of the 
   :cpp:class:`DistMatrix\<T,MC,MR>` class to be constructed without having to 
   manually specify a process grid, e.g., 

   .. code-block:: cpp

      // Build a 10 x 10 distributed matrix over mpi::COMM_WORLD
      elemental::DistMatrix<T,MC,MR> A( 10, 10 );


