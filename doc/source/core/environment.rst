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
          elem::Initialize( argc, argv );

          // ...

          elem::Finalize();
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

.. cpp:function:: Z Abs( const Z& alpha )

   Return the absolute value of the real variable :math:`\alpha`.

.. cpp:function:: Z Abs( const Complex<Z>& alpha )

   Return the absolute value of complex :math:`\alpha`: 

   .. math::

      |\alpha| = \sqrt{|\mathcal{R}(\alpha)|^2+|\mathcal{I}(\alpha)|^2}

.. cpp:function:: Z FastAbs( const Z& alpha )

   Same as :cpp:function::Abs: for real :math:`\alpha`.

.. cpp:function:: Z FastAbs( const Complex<Z>& alpha )

   Return a cheaper norm of the complex :math:`\alpha`:

   .. math::
   
      |\alpha|_{\mbox{fast}} = |\mathcal{R}(\alpha)| + |\mathcal{I}(\alpha)|

.. cpp:function:: Z RealPart( const Z& alpha )

   Return the real part of the real variable :math:`\alpha`, which is 
   :math:`\alpha` itself. 

.. cpp:function:: Z RealPart( const Complex<Z>& alpha )

   Return the real part of the complex variable :math:`\alpha`.

.. cpp:function:: Z ImagPart( const Z& alpha )

   Return the imaginary part of the real variable :math:`\alpha`, which is 
   trivially zero.

.. cpp:function:: Z ImagPart( const Complex<Z>& alpha )

   Return the imaginary part of the complex variable :math:`\alpha`.

.. cpp:function:: Z Conj( const Z& alpha )

   Return the complex conjugate of the real variable :math:`\alpha`,
   which is simply :math:`\alpha`.

.. cpp:function:: Complex<Z> Conj( const Complex<Z>& alpha )

   Return the complex conjugate of the complex variable :math:`\alpha`,

   .. math::

      \bar \alpha = \mathcal{R}(\alpha) - \mathcal{I}(\alpha) i

.. cpp:class:: Base<F>

   .. cpp:type:: type

      The underlying real datatype of the (potentially complex) datatype ``F``.
      For example, ``typename Base<Complex<double> >::type`` and 
      ``typename Base<double>::type`` are both equivalent to ``double``.
      This is often extremely useful in implementing routines which are 
      templated over real and complex datatypes but still make use of real 
      datatypes.

Custom datatypes
----------------

.. cpp:type:: byte

   ``typedef unsigned char byte;``

.. cpp:class:: Complex<R>

   .. cpp:type:: R BaseType

   .. cpp:member:: R real  
   
      The real part of the complex number

   .. cpp:member:: R imag

      The imaginary part of the complex number

   .. cpp:function:: Complex()
 
      This default constructor is a no-op.

   .. cpp:function:: Complex( R a )

      Construction from a real value.

   .. cpp:function:: Complex( R a, R b )
   
      Construction from a complex value.

   .. cpp:function:: Complex( const std::complex<R>& alpha )

      Construction from an ``std::complex<R>`` instance.

   .. cpp:function:: Complex<R>& operator=( const R& alpha )

      Assignment from a real value.

   .. cpp:function:: Complex<R>& operator+=( const R& alpha )

      Increment with a real value.

   .. cpp:function:: Complex<R>& operator-=( const R& alpha )

      Decrement with a real value.

   .. cpp:function:: Complex<R>& operator*=( const R& alpha )

      Scale with a real value.

   .. cpp:function:: Complex<R>& operator/=( const R& alpha )

      Divide with a real value.

   .. cpp:function:: Complex<R>& operator=( const Complex<R>& alpha )

      Assignment from a complex value.

   .. cpp:function:: Complex<R>& operator+=( const Complex<R>& alpha )

      Increment with a complex value.

   .. cpp:function:: Complex<R>& operator-=( const Complex<R>& alpha )

      Decrement with a complex value.

   .. cpp:function:: Complex<R>& operator*=( const Complex<R>& alpha )

      Scale with a complex value.

   .. cpp:function:: Complex<R>& operator/=( const Complex<R>& alpha )

      Divide with a complex value.

.. cpp:function:: Complex<R> operator+( const Complex<R>& alpha, const Complex<R>& beta )

   (complex,complex) addition.

.. cpp:function:: Complex<R> operator+( const Complex<R>& alpha, const R& beta )

   (complex,real) addition.

.. cpp:function:: Complex<R> operator+( const R& alpha, const Complex<R>& beta )

   (real,complex) addition.

.. cpp:function:: Complex<R> operator-( const Complex<R>& alpha, const Complex<R>& beta )

   (complex,complex) subtraction.

.. cpp:function:: Complex<R> operator-( const Complex<R>& alpha, R& beta )

   (complex,real) subtraction.

.. cpp:function:: Complex<R> operator-( const R& alpha, const Complex<R>& beta )

   (real,complex) subtraction.

.. cpp:function:: Complex<R> operator*( const Complex<R>& alpha, const Complex<R>& beta )

   (complex,complex) multiplication.

.. cpp:function:: Complex<R> operator*( const Complex<R>& alpha, R& beta )

   (complex,real) multiplication.

.. cpp:function:: Complex<R> operator*( const R& alpha, const Complex<R>& beta )

   (real,complex) multiplication.

.. cpp:function:: Complex<R> operator/( const Complex<R>& alpha, const Complex<R>& beta )

   (complex,complex) division.

.. cpp:function:: Complex<R> operator/( const Complex<R>& alpha, const R& beta )

   (complex,real) division.

.. cpp:function:: Complex<R> operator/( const R& alpha, const Complex<R>& beta )

   (real,complex) division.

.. cpp:function:: Complex<R> operator+( const Complex<R>& alpha )

   Returns `alpha`.

.. cpp:function:: Complex<R> operator-( const Complex<R>& alpha )

   Returns negative `alpha`.

.. cpp:function:: bool operator==( const Complex<R>& alpha, const Complex<R>& beta )

   (complex,complex) equality check.

.. cpp:function:: bool operator==( const Complex<R>& alpha, const R& beta )

   (complex,real) equality check.

.. cpp:function:: bool operator==( const R& alpha, const Complex<R>& beta )

   (real,complex) equality check.

.. cpp:function:: bool operator!=( const Complex<R>& alpha, const Complex<R>& beta )

   (complex,complex) inequality check.

.. cpp:function:: bool operator!=( const Complex<R>& alpha, const R& beta )

   (complex,real) inequality check.

.. cpp:function:: bool operator!=( const R& alpha, const Complex<R>& beta )

   (real,complex) inequality check.

.. cpp:function:: std::ostream& operator<<( std::ostream& os, Complex<R> alpha )

   Pretty prints `alpha` in the form ``a+bi``.

.. cpp:type:: scomplex

   ``typedef Complex<float> scomplex;``

.. cpp:type:: dcomplex

   ``typedef Complex<double> dcomplex;``

.. cpp:type:: UnitOrNonUnit

   An enum for specifying either ``UNIT`` or ``NON_UNIT``; typically used for 
   stating whether or not a triangular matrix's diagonal is explicitly stored
   (``NON_UNIT``) or is implicitly unit-diagonal (``UNIT``).

.. cpp:type:: Orientation

   An enum for specifying whether a matrix, say :math:`A`, should be implicitly 
   treated as :math:`A` (``NORMAL``), :math:`A^H` (``ADJOINT``), or :math:`A^T`
   (``TRANSPOSE``).

.. cpp:type:: LeftOrRight

   An enum for specifying ``LEFT`` or ``RIGHT``.

.. cpp:type:: ForwardOrBackward

   An enum for specifying ``FORWARD`` or ``BACKWARD``.

.. cpp:type:: UpperOrLower

   An enum for specifying ``LOWER`` or ``UPPER`` (triangular).

.. cpp:type:: VerticalOrHorizontal

   An enum for specifying ``VERTICAL`` or ``HORIZONTAL``.

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

   An enum for specifying either a ``ROW_MAJOR`` or ``COLUMN_MAJOR`` ordering;
   it is used to tune one of the algorithms in ``HermitianTridiag``
   which requires building a smaller square process grid from a rectangular 
   process grid, as the ordering of the processes can greatly impact 
   performance. See ``SetHermitianTridiagGridOrder``.

Custom exceptions
-----------------

.. cpp:class:: SingularMatrixException

   An extension of ``std::runtime_error`` which is meant to be thrown when 
   a singular matrix is unexpectedly encountered.

   .. cpp:function:: SingularMatrixException( const char* msg="Matrix was singular" )

      Builds an instance of the exception which allows one to optionally 
      specify the error message.

   .. code-block:: cpp

      throw elem::SingularMatrixException();

.. cpp:class:: NonHPDMatrixException 

   An extension of ``std::runtime_error`` which is meant to be thrown when
   a non positive-definite Hermitian matrix is unexpectedly encountered
   (e.g., during Cholesky factorization).

   .. cpp:function:: NonHPDMatrixException( const char* msg="Matrix was not HPD" )

      Builds an instance of the exception which allows one to optionally 
      specify the error message.

   .. code-block:: cpp

      throw elem::NonHPDMatrixException();

.. cpp:class:: NonHPSDMatrixException 

   An extension of ``std::runtime_error`` which is meant to be thrown when
   a non positive semi-definite Hermitian matrix is unexpectedly encountered
   (e.g., during computation of the square root of a Hermitian matrix).

   .. cpp:function:: NonHPSDMatrixException( const char* msg="Matrix was not HPSD" )

      Builds an instance of the exception which allows one to optionally 
      specify the error message.

   .. code-block:: cpp

      throw elem::NonHPSDMatrixException();

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
      elem::DistMatrix<T,MC,MR> A( 10, 10 );


