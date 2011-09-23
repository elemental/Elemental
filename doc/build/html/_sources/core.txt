Core functionality
******************

The Matrix class
================
This is the basic building block of the library: its purpose it to provide 
convenient mechanisms for performing basic matrix manipulation operations, 
such as setting and querying individual matrix entries, without giving up 
compatibility with interfaces such as BLAS and LAPACK, which assume column-major
storage.

.. note:: 

   There are essentially four constraints on the underlying datatype, say 
   ``T``:

    1. The datatype must be self-contained in the sense that all of 
       the relevant data of an instance of the datatype is contained within a 
       contiguous buffer of size less than or equal to ``sizeof(T)``.
    2. It must be sensical to cast the decimal literals ``0`` and ``1`` into 
       type ``T``, e.g., ``static_cast<T>(0)`` and ``static_cast<T>(1)``. 
       This requirement exists so that a routine can be provided to initialize
       a matrix to the identity operator.
    3. Setting all ``sizeof(T)`` bytes of an instance of type ``T`` to zero must
       be logically equivalent to setting the instance to ``static_cast<T>(0)``.
       This constraint is enforced so that ``std::memset`` can be used to 
       quickly set matrices to zero.
    4. If the ``SetToRandom()`` member function is to be used, then the routine 
       ``T elemental::SampleUnitBall<T>()`` must exist. This is already 
       provided for ``int``, ``float``, ``double``, ``std::complex<float>``, 
       and ``std::complex<double>``, and could be easily extended to new 
       datatypes.

    It is important to note that these requirements are not strictly weaker or 
    strictly stronger than requiring ``T`` to be POD (Plain Old Data), as 
    requirements (2) and (3) are not necessary for a type to be POD, whereas 
    POD types are not allowed to have user-defined constructors. Any complex
    variable (defined using ``std::complex``) fits our requirements, but is not
    considered POD since it has several custom constructors.

An example of generating an :math:`m \times n` matrix of real double-precision 
numbers where the :math:`(i,j)` entry is equal to :math:`i-j` would be:

  .. code-block:: cpp

     #include "elemental.hpp"
     using namespace elemental;
     ...
     Matrix<double> A( m, n );
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             A.Set( i, j, (double)i-j );
     
The underlying data storage is simply a contiguous buffer that stores entries 
in a column-major fashion with an arbitrary leading dimension. For modifiable
instances of the ``Matrix`` class, the routine
``T* Matrix<T>::Buffer()`` returns a pointer to the underlying 
buffer, while ``int Matrix<T>::LDim() const`` returns the leading 
dimension; these two routines could be used to directly perform the equivalent
of the first code sample as follows:

  .. code-block:: cpp
     
     #include "elemental.hpp"
     using namespace elemental;
     ...
     Matrix<double> A( m, n );
     double* buffer = A.Buffer();
     const int ldim = A.LDim();
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             buffer[i+j*ldim] = (double)i-j;

For constant instances of the ``Matrix`` class, a ``const`` pointer
to the underlying data can similarly be returned with a call to 
``A.LockedBuffer()``. In addition, a (``const``) pointer to the place in the 
(``const``) buffer where entry :math:`(i,j)` resides can be easily retrieved
with a call to ``A.Buffer(i,j)`` or ``A.LockedBuffer(i,j)``.

It is also important to be able to create matrices which are simply *views* 
of existing (sub)matrices. For example, if ``A`` is a :math:`10 \times 10` 
matrix of complex doubles, then a matrix :math:`A_{BR}` can easily be created 
to view the bottom-right :math:`6 \times 7` submatrix using

  .. code-block:: cpp

     #include "elemental.hpp"
     ...
     Matrix<std::complex<double> > ABR;
     ABR.View( A, 4, 3, 6, 7 );

since the bottom-right :math:`6 \times 7` submatrix beings at index 
:math:`(4,3)`. In general, to view the :math:`M \times N` submatrix starting
at entry :math:`(i,j)`, one would call ``ABR.View( A, i, j, M, N );``.

.. cpp:class:: Matrix<T>

   .. rubric:: Constructors

   .. cpp:function:: Matrix()

      This simply creates a default :math:`0 \times 0` matrix with a leading 
      dimension of one (BLAS and LAPACK require positive leading dimensions).

   .. cpp:function:: Matrix( int height, int width )

      A *height* :math:`\times` *width* matrix is created with an unspecified
      leading dimension (though it is currently implemented as 
      ``std::max(height,1)``).

   .. cpp:function:: Matrix( int height, int width, int ldim )

      A *height* :math:`\times` *width* matrix is created with a leading 
      dimension equal to *ldim* (which must be greater than or equal 
      ``std::min(height,1)``).

   .. cpp:function:: Matrix( int height, int width, const T* buffer, int ldim )

      A matrix is built around column-major constant buffer ``const T* buffer`` 
      with the specified dimensions. The memory pointed to by ``buffer`` should
      not be freed until after the ``Matrix`` object is destructed.

   .. cpp:function:: Matrix( int height, int width, T* buffer, int ldim )

      A matrix is built around the column-major modifiable buffer ``T* buffer``
      with the specified dimensions. The memory pointed to by ``buffer`` should
      not be freed until after the ``Matrix`` object is destructed.

   .. cpp:function:: Matrix( const Matrix<T>& A )

      A copy (not a view) of the matrix :math:`A` is built.

   .. rubric:: Basic information

   .. cpp:function:: int Height() const

      Return the height of the matrix.

   .. cpp:function:: int Width() const

      Return the width of the matrix.

   .. cpp:function:: int DiagonalLength( int offset=0 ) const

      Return the length of the specified diagonal of the matrix: an offset of 
      :math:`0` refers to the main diagonal, an offset of :math:`1` refers to 
      the superdiagonal, an offset of :math:`-1` refers to the subdiagonal, 
      etc.

   .. cpp:function:: int LDim() const

      Return the leading dimension of the underlying buffer.

   .. cpp:function:: int MemorySize() const

      Return the number of entries of type ``T`` that this ``Matrix`` instance 
      has allocated space for.

   .. cpp:function:: T* Buffer()

      Return a pointer to the underlying buffer.

   .. cpp:function:: const T* LockedBuffer() const

      Return a pointer to the underlying buffer that does not allow for 
      modifying the data.

   .. cpp:function:: T* Buffer( int i, int j )

      Return a pointer to the portion of the buffer that holds entry 
      :math:`(i,j)`.

   .. cpp:function:: const T* LockedBuffer( int i, int j ) const

      Return a pointer to the portion of the buffer that holds entry
      :math:`(i,j)` that does not allow for modifying the data.

   .. cpp:function:: T* Buffer( int i, int j, int height, int width )

      Same as the version without *height* and *width*, but in **Debug** modes 
      it will ensure that the *height* :math:`\times` *width* submatrix starting
      at entry :math:`(i,j)` does not go out of bounds.

   .. cpp:function:: const T* LockedBuffer( int i, int j, int height, int width ) const

      Same as above, but the data cannot be modified using the returned pointer.

   .. rubric:: I/O

   .. cpp:function:: void Print( const std::string msg="" ) const

   The matrix is printed to standard output (``std::cout``) with the preceding
   message ``msg`` (which is empty if unspecified).

   .. cpp:function:: void Print( std::ostream& os, const std::string msg="" ) const

      The matrix is printed to the output stream ``os`` with the preceding 
      message ``msg`` (which is empty if unspecified).

   .. rubric:: Entry manipulation

   .. cpp:function:: T Get( int i, int j ) const

      Return entry :math:`(i,j)`.

   .. cpp:function:: void Set( int i, int j, T alpha )

      Set entry :math:`(i,j)` to :math:`\alpha`.

   .. cpp:function:: void Update( int i, int j, T alpha )

      Add :math:`\alpha` to entry :math:`(i,j)`.

   .. cpp:function:: void GetDiagonal( Matrix<T>& d, int offset=0 ) const

      Modify :math:`d` into a column-vector containing the entries lying on the 
      ``offset`` diagonal of our matrix (for instance, the main diagonal has 
      offset :math:`0`, the subdiagonal has offset :math:`-1`, and the 
      superdiagonal has offset :math:`+1`).

   .. cpp:function:: void SetDiagonal( const Matrix<T>& d, int offset=0 )

      Set the entries in the ``offset`` diagonal entries from the contents of 
      the column-vector :math:`d`.

   .. cpp:function:: void UpdateDiagonal( const Matrix<T>& d, int offset=0 )

      Add the contents of :math:`d` onto the entries in the ``offset`` diagonal.

   .. note::

      The remainder of this group is only valid for complex datatypes.

   .. cpp:function:: typename RealBase<T>::type GetReal( int i, int j ) const

      Return the real part of entry :math:`(i,j)`.

   .. cpp:function:: typename RealBase<T>::type GetImag( int i, int j ) const

      Return the imaginary part of entry :math:`(i,j)`.

   .. cpp:function:: void SetReal( int i, int j, typename RealBase<T>::type alpha )

      Set the real part of entry :math:`(i,j)` to :math:`\alpha`.

   .. cpp:function:: void SetImag( int i, int j, typename RealBase<T>::type alpha )

      Set the imaginary part of entry :math:`(i,j)` to :math:`\alpha`.

   .. cpp:function:: void UpdateReal( int i, int j, typename RealBase<T>::type alpha )

      Add :math:`\alpha` to the real part of entry :math:`(i,j)`.

   .. cpp:function:: void UpdateImag( int i, int j, typename RealBase<T>::type alpha ) 

      Add :math:`\alpha` to the imaginary part of entry :math:`(i,j)`.

   .. cpp:function:: void GetRealDiagonal( Matrix<typename RealBase<T>::type>& d, int offset=0 ) const

      Modify :math:`d` into a column-vector containing the real parts of the
      entries in the ``offset`` diagonal.

   .. cpp:function:: void GetImagDiagonal( Matrix<typename RealBase<T>::type>& d, int offset=0 ) const

      Modify :math:`d` into a column-vector containing the imaginary parts of 
      the entries in the ``offset`` diagonal.

   .. cpp:function:: void SetRealDiagonal( const Matrix<typename RealBase<T>::type>& d, int offset=0 )

      Set the real parts of the entries in the ``offset`` diagonal from the 
      contents of the column-vector :math:`d`.

   .. cpp:function:: void SetImagDiagonal( const Matrix<typename RealBase<T>::type>& d, int offset=0 )

      Set the imaginary parts of the entries in the ``offset`` diagonal from 
      the column-vector :math:`d`.

   .. cpp:function:: void UpdateRealDiagonal( const Matrix<typename RealBase<T>::type>& d, int offset=0 )

      Add the contents of the column-vector :math:`d` onto the real parts of the
      entries in the ``offset`` diagonal.

   .. cpp:function:: void UpdateImagDiagonal( const Matrix<typename RealBase<T>::type>& d, int offset=0 )

      Add the contents of the column-vector :math:`d` onto the imaginary parts 
      of the entries in the ``offset`` diagonal.

   .. rubric:: Views

   .. cpp:function:: bool Viewing() const

      Return whether or not this matrix is currently viewing another matrix.

   .. cpp:function:: bool LockedView() const

      Return whether or not we can modify the data we are viewing.

   .. cpp:function:: void View( int height, int width, T* buffer, int ldim )

      Reconfigure the matrix around the specified buffer.

   .. cpp:function:: void View( Matrix<T>& A )

      Reconfigure the matrix around the modifiable buffer underlying ``A``.

   .. cpp:function:: void LockedView( int height, int width, const T* buffer, int ldim )

      Reconfigure the matrix around the specified unmodifiable buffer.

   .. cpp:function:: void LockedView( const Matrix<T>& A )

      Reconfigure the matrix around the unmodifiable buffer underlying ``A``.

   .. cpp:function:: void View( Matrix<T>& A, int i, int j, int height, int width )

      Reconfigure the matrix around the modifiable buffer underlying ``A``, but
      only the portion that holds the *height* :math:`\times` *width* submatrix 
      starting at entry ``(i,j)``

   .. cpp:function:: void LockedView( const Matrix<T>& A, int i, int j, int height, int width )

      Same as above, but the resulting matrix data is unmodifiable.

   .. cpp:function:: void View1x2( Matrix<T>& AL, Matrix<T>& AR )

      Reconfigure the matrix to use the modifiable buffer that spans the 
      matrices :math:`A_L` and :math:`A_R` such that it behaves like 
      :math:`[A_L A_R]` (this routine requires that :math:`A_R`'s buffer begins 
      at the same memory location that an extra column of :math:`A_L` would 
      have).

   .. cpp:function:: void LockedView1x2( const Matrix<T>& AL, const Matrix<T>& AR )

      Same as above, but the resulting matrix data is unmodifiable.

   .. cpp:function:: void View2x1( Matrix<T>& AT, Matrix<T>& AB )

      Reconfigure the matrix to use the modifiable buffer that spans the 
      matrices :math:`A_T` and :math:`A_B` such that it behaves like 
      :math:`[A_T;A_B]` (this routine requires that :math:`A_B`'s buffer begins 
      at the same memory location that an extra row of :math:`A_T` would have).

   .. cpp:function:: void LockedView2x1( const Matrix<T>& AT, const Matrix<T>& AB )

      Same as above, but the resulting matrix data is unmodifiable.

   .. cpp:function:: void View2x2( Matrix<T>& ATL, Matrix<T>& ATR, Matrix<T>& ABL, Matrix<T>& ABR )

      Reconfigure the matrix to behave like 
      :math:`[A_{TL} A_{TR}; A_{BL} A_{BR}]`
      (the buffer requirements are similar to ``View1x2`` and ``View2x1``).

   .. cpp:function:: void LockedView2x2( const Matrix<T>& ATL, const Matrix<T>& ATR, const Matrix<T>& ABL, const Matrix<T>& ABR )

      Same as above, but the resulting matrix data is unmodifiable.

   .. rubric:: Utilities

   .. cpp:function:: const Matrix<T>& operator=( const Matrix<T>& A )

      Create a copy of matrix :math:`A`.

   .. cpp:function:: void Empty()

      Sets the matrix to :math:`0 \times 0` and frees the underlying buffer.

   .. cpp:function:: void ResizeTo( int height, int width )

      Reconfigures the matrix to be *height* :math:`\times` *width*.

   .. cpp:function:: void ResizeTo( int height, int width, int ldim )

      Reconfigures the matrix to be *height* :math:`\times` *width*, but with 
      leading dimension equal to *ldim* (which must be greater than or equal to 
      ``std::min(height,1)``).

   .. cpp:function:: void SetToIdentity()

      Sets the entire matrix to zero, with the exception of the main diagonal 
      being set to one. For square matrices, this corresponds to the identity 
      operator.

   .. cpp:function:: void SetToRandom()

      Sets each entry in the matrix to a uniform sample from the most natural 
      interpretation of the unit ball specified by the datatype.

   .. cpp:function:: void SetToZero()

      Sets every entry of the matrix to zero.

The DistMatrix class
====================
The ``DistMatrix`` class is meant to provide a distributed-memory analogue of 
the ``Matrix`` class. Similar to PLAPACK, roughly ten different matrix 
distributions are provided and it is trivial (in the programmability sense) to 
redistribute from one to another: in PLAPACK, one would simply call 
``PLA_Copy``, whereas, in Elemental, it is handled through overloading the 
:math:`=` operator.

Since it is crucial to know not only how many 
processes to distribute the data over, but *which* processes, and in what 
manner they should be decomposed into a logical two-dimensional grid, an 
instance of the ``Grid`` class must be passed into the constructor of 
the ``DistMatrix`` class.

.. note:: 
   
   Since the ``DistMatrix`` class makes use of MPI for message passing, 
   custom interfaces must be written for nonstandard datatypes. As of now, 
   the following datatypes are fully supported for ``DistMatrix``:
   ``int``, ``float``, ``double``, ``std::complex<float>``, and
   ``std::complex<double>``.

AbstractDistMatrix
------------------

This abstract class defines the list of member functions that are guaranteed 
to be available for all matrix distributions.

.. cpp:class:: AbstractDistMatrix<T>

   .. rubric:: Basic information

   .. cpp:function:: int Height() const

      Return the height of the matrix.

   .. cpp:function:: int Width() const

      Return the width of the matrix.

   .. cpp:function:: int LocalHeight() const

      Return the local height of the matrix.

   .. cpp:function:: int LocalWidth() const

      Return the local width of the matrix.

   .. cpp:function:: int LocalLDim() const

      Return the local leading dimension of the matrix.

   .. cpp:function:: size_t AllocatedMemory() const

      Return the number of entries of type ``T`` that we have locally allocated
      space for.

   .. cpp:function:: const elemental::Grid& Grid() const

      Return the grid that this distributed matrix is distributed over.

   .. cpp:function:: T* LocalBuffer( int iLocal=0, int jLocal=0 )

      Return a pointer to the portion of the local buffer that stores entry 
      ``(iLocal,jLocal)``.

   .. cpp:function:: const T* LockedLocalBuffer( int iLocal=0, int jLocal=0 ) const

      Return a pointer to the portion of the local buffer that stores entry
      ``(iLocal,jLocal)``, but do not allow for the data to be modified through
      the returned pointer.

   .. cpp:function:: Matrix<T>& LocalMatrix()

      Return a reference to the local matrix.

   .. cpp:function:: const Matrix<T>& LockedLocalMatrix() const

      Return an unmodifiable reference to the local matrix.

   .. rubric:: I/O

   .. cpp:function:: void Print( const std::string msg="" ) const

      Print the distributed matrix to standard output (``std::cout``).

   .. cpp:function:: void Print( std::ostream& os, const std::string msg="" ) const

      Print the distributed matrix to the output stream ``os``.

   .. cpp:function:: void Write( const std::string filename, const std::string msg="" ) const

      Print the distributed matrix to the file named ``filename``.

   .. rubric:: Alignments

   .. cpp:function:: void FreeAlignments()

      Free all alignment constaints.

   .. cpp:function:: bool ConstrainedColAlignment() const

      Return whether or not the column alignment is constrained.

   .. cpp:function:: bool ConstrainedRowAlignment() const

      Return whether or not the row alignment is constrained.

   .. cpp:function:: int ColAlignment() const

      Return the alignment of the columns of the matrix.

   .. cpp:function:: int RowAlignment() const

      Return the alignment of the rows of the matrix.

   .. cpp:function:: int ColShift() const

      Return the first global row that our process owns.

   .. cpp:function:: int RowShift() const

      Return the first global column that our process owns.

   .. rubric:: Entry manipulation

   .. cpp:function:: T Get( int i, int j ) const

      Return the ``(i,j)`` entry of the global matrix. This operation is 
      collective.

   .. cpp:function:: void Set( int i, int j, T alpha )

      Set the ``(i,j)`` entry of the global matrix to :math:`\alpha`. This 
      operation is collective.

   .. cpp:function:: void Update( int i, int j, T alpha )

      Add :math:`\alpha` to the ``(i,j)`` entry of the global matrix. This 
      operation is collective.

   .. cpp:function:: T GetLocalEntry( int iLocal, int jLocal ) const

      Return the ``(iLocal,jLocal)`` entry of our local matrix.

   .. cpp:function:: void SetLocalEntry( int iLocal, int jLocal, T alpha )

      Set the ``(iLocal,jLocal)`` entry of our local matrix to :math:`\alpha`.

   .. cpp:function:: void UpdateLocalEntry( int iLoca, int jLocal, T alpha )

      Add :math:`\alpha` to the ``(iLocal,jLocal)`` entry of our local matrix.

   .. note::

      The remainder of this group is only valid for complex datatypes.

   .. cpp:function:: typename RealBase<T>::type GetReal( int i, int j ) const

      Return the real part of the ``(i,j)`` entry of the global matrix. This
      operation is collective.

   .. cpp:function:: typename RealBase<T>::type GetImag( int i, int j ) const

      Return the imaginary part of the ``(i,j)`` entry of the global matrix. 
      This operation is collective.

   .. cpp:function:: void SetReal( int i, int j, typename RealBase<T>::type alpha )

      Set the real part of the ``(i,j)`` entry of the global matrix to 
      :math:`\alpha`.

   .. cpp:function:: void SetImag( int i, int j, typename RealBase<T>::type alpha )

      Set the imaginary part of the ``(i,j)`` entry of the global matrix to 
      :math:`\alpha`.

   .. cpp:function:: void UpdateReal( int i, int j, typename RealBase<T>::type alpha )

      Add :math:`\alpha` to the real part of the ``(i,j)`` entry of the global 
      matrix.

   .. cpp:function:: void UpdateImag( int i, int j, typename RealBase<T>::type alpha )

      Add :math:`\alpha` to the imaginary part of the ``(i,j)`` entry of the 
      global matrix.

   .. cpp:function:: typename RealBase<T>::type GetRealLocalEntry( int iLocal, int jLocal ) const

      Return the real part of the ``(iLocal,jLocal)`` entry of our local matrix.

   .. cpp:function:: typename RealBase<T>::type GetImagLocalEntry( int iLocal, int jLocal ) const

      Return the imaginary part of the ``(iLocal,jLocal)`` entry of our local 
      matrix.

   .. cpp:function:: void SetRealLocalEntry( int iLocal, int jLocal, typename RealBase<T>::type alpha )

      Set the real part of the ``(iLocal,jLocal)`` entry of our local matrix.

   .. cpp:function:: void SetImagLocalEntry( int iLocal, int jLocal, typename RealBase<T>::type alpha )

      Set the imaginary part of the ``(iLocal,jLocal)`` entry of our local 
      matrix.

   .. cpp:function:: void UpdateRealLocalEntry( int iLocal, int jLocal, typename RealBase<T>::type alpha )

      Add :math:`\alpha` to the real part of the ``(iLocal,jLocal)`` entry of 
      our local matrix.

   .. cpp:function:: void UpdateImagLocalEntry( int iLocal, int jLocal, typename RealBase<T>::type alpha )

      Add :math:`\alpha` to the imaginary part of the ``(iLocal,jLocal)`` entry 
      of our local matrix.

   .. rubric:: Viewing

   .. cpp:function:: bool Viewing() const

      Return whether or not this ``DistMatrix`` is viewing another.

   .. cpp:function:: bool LockedView() const

      Return whether or not this ``DistMatrix`` is viewing another in a manner
      that does not allow for modifying the viewed data.

   .. rubric:: Utilities

   .. cpp:function:: void Empty()

      Resize the distributed matrix so that it is :math:`0 \times 0` and free 
      all allocated storage.

   .. cpp:function:: void MakeTrapezoidal( Side side, Shape shape, int offset=0 )

      Explicitly introduce zeroes into the distributed matrix such that it is 
      trapezoidal with respect to the left or right diagonal (as chosen by the 
      ``side`` parameter). Whether or not the matrix is lower or upper 
      trapezoidal is determined by the ``shape`` parameter, and the diagonal 
      offset is chosen by the ``offset`` parameter (:math:`0` denotes the main 
      diagonal, :math:`-1` denotes the subdiagonal, and :math:`+1` denotes the 
      superdiagonal).

   .. cpp:function:: void ScaleTrapezoidal( T alpha, Side side, Shape shape, int offset=0 )

      Scale the portion of the matrix determined by the above discussion by the 
      scalar :math:`\alpha`.

   .. cpp:function:: void ResizeTo( int height, int width )

      Reconfigure the matrix so that it is *height* :math:`\times` *width*.

   .. cpp:function:: void SetToIdentity()

      Set the entire matrix to zero and then introduce ones onto the main 
      diagonal.

   .. cpp:function:: void SetToRandom()

      Independently draw each entry of the matrix from the uniform distribution
      over the unit ball.

   .. cpp:function:: void SetToRandomHermitian()

      Same as above, but the diagonal is forced to be real-valued
      (the rest of the symmetry is implicit).

   .. cpp:function:: void SetToRandomHPD()

      Same as above, but a sufficiently large constant is added to every 
      diagonal entry in order to ensure that the matrix is positive-definite.

   .. cpp:function:: void SetToZero()

      Set all entries of the distributed matrix to zero.

``[MC,MR]``
-----------

This is the standard matrix distribution... **left off here**

``[MC,* ]``
-----------

``[* ,MR]``
-----------

``[MR,MC]``
-----------

``[MR,* ]``
-----------

``[* ,MC]``
-----------

``[MD,* ]``
-----------

``[* ,MD]``
-----------

``[VC,* ]``
-----------

``[* ,VC]``
-----------

``[VR,* ]``
-----------

``[* ,VR]``
-----------

``[* ,* ]``
-----------

The Grid class
==============

Constructors
------------


The Axpy interface
======================

Environment routines
====================
