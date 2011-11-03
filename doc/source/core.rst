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

The Grid class
==============

This class is responsible for converting MPI communicators into a 
two-dimensional process grid meant for distributing matrices (ala the 
soon-to-be-discussed ``DistMatrix`` class).

.. cpp:class:: Grid

   .. rubric:: Basic constructors

   .. cpp:function:: Grid( mpi::Comm comm=mpi::COMM_WORLD )

      Construct a process grid over the specified communicator and let Elemental
      decide the process grid dimensions. If no communicator is specified, 
      mpi::COMM_WORLD is used.

   .. cpp:function:: Grid( mpi::Comm comm, int height, int width )

      Construct a process grid over the specified communicator with the 
      given dimensions. Note that the size of the communicator should be 
      *height* :math:`\times` *width*.

   .. rubric:: Basic information

   .. cpp:function:: bool InGrid() const

      Return whether or not our process is actively participating in the process
      grid.

   .. cpp:function:: int Size() const

      Return the number of active processes in the process grid. This number 
      is equal to *Height()* :math:`\times` *Width()*.

   .. cpp:function:: int Height() const

      Return the height of the process grid.

   .. cpp:function:: int Width() const

      Return the width of the process grid.

   .. cpp:function:: int GCD() const

      Return the greatest common denominator of the height and width of the 
      process grid.

   .. cpp:function:: int LCM() const

      Return the lowest common multiple of the height and width of the process
      grid.

   .. cpp:function:: int MCRank() const

      Return our process's rank in the *MC* (Matrix Column) communicator. This 
      corresponds to our row in the process grid.

   .. cpp:function:: int MRRank() const
     
      Return our process's rank in the *MR* (Matrix Row) communicator. This
      corresponds to our column in the process grid.

   .. cpp:function:: int VCRank() const

      Return our process's rank in the *VC* (Vector Column) communicator. This
      corresponds to our rank in a column-major ordering of the process grid.

   .. cpp:function:: int VRRank() const

      Return our process's rank in the *VR* (Vector Row) communicator. This 
      corresponds to our rank in a row-major ordering of the process grid.

   .. cpp:function:: mpi::Comm MCComm() const

      Return the *MC* (Matrix Column) communicator. This consists of the set
      of processes within our column of the grid (ordered top-to-bottom).

   .. cpp:function:: mpi::Comm MRComm() const

      Return the *MR* (Matrix Row) communicator. This consists of the set of
      processes within our row of the grid (ordered left-to-right).

   .. cpp:function:: mpi::Comm VCComm() const

      Return the *VC* (Vector Column) communicator. This consists of the entire
      set of processes in the grid, but ordered in a column-major fashion.

   .. cpp:function:: mpi::Comm VRComm() const

      Return the *VR* (Vector Row) communicator. This consists of the entire 
      set of processes in the grid, but ordered in a row-major fashion.

   .. rubric:: Advanced routines

   .. cpp:function:: Grid( mpi::Comm viewingComm, mpi::Group owningGroup )

      Construct a process grid where only a subset of the participating 
      processes should actively participate in the process grid. In particular,
      *viewingComm* should consist of the set of all processes constructing 
      this ``Grid`` instance, and *owningGroup* should define a subset of the
      processes in *viewingComm*. Elemental then chooses the grid dimensions. 
      Most users should not call this routine, as this type of grid is only 
      supported for a few ``DistMatrix`` types.

   .. cpp:function:: Grid( mpi::Comm viewingComm, mpi::Group owningGroup, int height, int width )

      This is the same as the previous routine, but the process grid dimensions
      are explicitly specified, and it is required that *height* :math:`\times`
      *width* equals the size of *owningGroup*. Most users should not call this
      routine, as it is only supported for a few ``DistMatrix`` types.

   .. cpp:function:: int OwningRank() const

      Return our process's rank within the set of processes that are actively
      participating in the grid.

   .. cpp:function:: int ViewingRank() const

      Return our process's rank within the entire set of processes that 
      constructed this grid.

   .. cpp:function:: int VCToViewingMap() const

      Map the given column-major grid rank to the rank in the (potentially)
      larger set of processes which constructed the grid.

   .. cpp:function:: mpi::Group OwningGroup() const

      Return the group of processes which is actively participating in the 
      grid.

   .. cpp:function:: mpi::Comm OwningComm() const

      Return the communicator for the set of processes actively participating
      in the grid. Note that this can only be valid if the calling process
      is an active member of the grid!

   .. cpp:function:: mpi::Comm ViewingComm() const

      Return the communicator for the entire set of processes which constructed
      the grid.

   .. cpp:function:: int DiagPath() const

      Return our unique diagonal index in an tesselation of the process grid.

   .. cpp:function:: int DiagPath( int vectorColRank ) const

      Return the unique diagonal index of the process with the given 
      column-major vector rank in an tesselation of the process grid.

   .. cpp:function:: int DiagPathRank() const

      Return our process's rank out of the set of processes lying in our 
      diagonal of the tesselation of the process grid.

   .. cpp:function:: int DiagPathRank( int vectorColRank ) const

      Return the rank of the given process out of the set of processes in its
      diagonal of the tesselation of the process grid.

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

   .. cpp:function:: void SetGrid( const elemental::Grid& grid )

      Clear the distributed matrix's contents and reconfigure for the new 
      process grid.

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

This is by far the most important matrix distribution in Elemental, as the vast
majority of parallel routines expect the input to be in this form. For a
:math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process grid,
individual entries would be owned by the following processes (assuming the 
column and row alignments are both 0):

.. math::
   :nowrap:

   \[
   \left(\begin{array}{cccccccccc}
     0 & 2 & 4 & 0 & 2 & 4 & 0 \\
     1 & 3 & 5 & 1 & 3 & 5 & 1 \\ 
     0 & 2 & 4 & 0 & 2 & 4 & 0 \\
     1 & 3 & 5 & 1 & 3 & 5 & 1 \\ 
     0 & 2 & 4 & 0 & 2 & 4 & 0 \\
     1 & 3 & 5 & 1 & 3 & 5 & 1 \\ 
     0 & 2 & 4 & 0 & 2 & 4 & 0  
   \end{array}\right)
   \]

Similarly, if the column alignment is kept at 0 and the row alignment is changed
to 2 (meaning that the third process column owns the first column of the 
matrix), the individual entries would be owned as follows:

.. math::
   :nowrap:

   \[
   \left(\begin{array}{cccccccccc}
     4 & 0 & 2 & 4 & 0 & 2 & 4 \\
     5 & 1 & 3 & 5 & 1 & 3 & 5 \\ 
     4 & 0 & 2 & 4 & 0 & 2 & 4 \\
     5 & 1 & 3 & 5 & 1 & 3 & 5 \\ 
     4 & 0 & 2 & 4 & 0 & 2 & 4 \\
     5 & 1 & 3 & 5 & 1 & 3 & 5 \\ 
     4 & 0 & 2 & 4 & 0 & 2 & 4 
   \end{array}\right)
   \]

.. cpp:class:: DistMatrix<T,MC,MR>

   .. rubric:: Constructors

   .. cpp:function:: DistMatrix( const elemental::Grid& grid=DefaultGrid() )
      
      Create a :math:`0 \times 0` distributed matrix over the specified grid.

   .. cpp:function:: DistMatrix( int height, int width, const elemental::Grid& grid=DefaultGrid() )

      Create a ``height`` :math:`\times` ``width`` distributed matrix over the
      specified grid.

   .. cpp:function:: DistMatrix( int height, int width, bool constrainedColAlignment, bool constrainedRowAlignment, int colAlignment, int rowAlignment, const elemental::Grid& grid )

      Create a ``height`` :math:`\times` ``width`` distributed matrix 
      distributed over the specified process grid, but with the top-left entry
      owned by the ``colAlignment`` process row and the ``rowAlignment`` 
      process column. Each of these alignments may be *constrained* to remain
      constant when redistributing data into this ``DistMatrix``.

   .. cpp:function:: DistMatrix( int height, int width, bool constrainedColAlignment, bool constrainedRowAlignment, int colAlignment, int rowAlignment, int ldim, const elemental::Grid& grid )

      Same as above, but the local leading dimension is also specified.

   .. cpp:function:: DistMatrix( int height, int width, int colAlignment, int rowAlignment, const T* buffer, int ldim, const elemental::Grid& grid )

      View a constant distributed matrix's buffer; the buffer must correspond 
      to the local portion of an elemental distributed matrix with the 
      specified row and column alignments and leading dimension, ``ldim``.

   .. cpp:function:: DistMatrix( int height, int width, int colAlignment, int rowAlignment, T* buffer, int ldim, const elemental::Grid& grid )

      Same as above, but the contents of the matrix are modifiable.

   .. cpp:function:: DistMatrix( const DistMatrix<T,U,V>& A )

      Build a copy of the distributed matrix ``A``, but force it to be in the
      ``[MC,MR]`` distribution.

   .. rubric:: Redistribution

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MC,MR>& A )

      If this matrix can be properly aligned with ``A``, then perform a local
      copy, otherwise perform an ``mpi::SendRecv`` permutation first.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MC,STAR>& A )

      Perform a local (filtered) copy to form an ``[MC,MR ]`` distribution and 
      then, if necessary, fix the alignment of the ``MC`` distribution via an 
      ``mpi::SendRecv`` within process columns.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,MR>& A )
       
      Perform a local (filtered) copy to form an ``[MC,MR ]`` distribution and 
      then, if necessary, fix the alignment of the ``MR`` distribution via an 
      ``mpi::SendRecv`` within process rows.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MD,STAR>& A )

      Since the ``[MD,STAR]`` distribution is defined such that its columns are
      distributed like a diagonal of an ``[MC,MR]`` distributed matrix, this 
      operation is not very common. 

      .. note::
         This redistribution routine is not yet implemented.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,MD>& A )

      .. note::
         This redistribution routine is not yet implemented.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MR,MC>& A )

      This routine serves to transpose the distribution of ``A[MR,MC]`` into 
      the standard matrix distribution, ``A[MC,MR]``. This redistribution is 
      implemented with four different approaches: one for matrices that are 
      taller than they are wide, one for matrices that are wider than they are 
      tall, one for column vectors, and one for row vectors.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MR,STAR>& A )

      This is similar to the above routine, but with each row of ``A`` being 
      undistributed, and only one approach is needed: 
      :math:`A[M_C,M_R] \leftarrow A[V_C,\star] \leftarrow A[V_R,\star] \leftarrow A[M_R,\star]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,MC>& A )

      This routine is dual to the :math:`A[M_C,M_R] \leftarrow A[M_R,\star]` 
      redistribution and is accomplished through the sequence: 
      :math:`A[M_C,M_R] \leftarrow A[\star,V_R] \leftarrow A[\star,V_C] \leftarrow A[\star,M_C]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,VC,STAR>& A )

      Perform an ``mpi::AllToAll`` within process rows in order to redistribute
      to the ``[MC,MR]`` distribution (an ``mpi::SendRecv`` within process 
      columns may be required for alignment).

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,VC>& A )

      Accomplished through the sequence 
      :math:`A[M_C,M_R] \leftarrow A[\star,V_R] \leftarrow A[\star,V_C]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,VR,STAR>& A )

      Accomplished through the sequence
      :math:`A[M_C,M_R] \leftarrow A[V_C,\star] \leftarrow A[V_R,\star]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,VR>& A )

      Perform an ``mpi::AllToAll`` within process columns in order to 
      redistribute to the ``[MC,MR]`` distribution (an ``mpi::SendRecv`` within
      process rows may be required for alignment).

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,STAR>& A )

      Perform an ``mpi::AllGather`` over the entire grid in order to give every
      process a full copy of ``A``.

   .. rubric:: Diagonal manipulation

   .. cpp:function:: void GetDiagonal( DistMatrix<T,MD,STAR>& d, int offset=0 ) const

      The :math:`[M_D,\star]` distribution is defined such that its columns 
      are distributed like diagonals of the standard matrix distribution, 
      `[M_C,M_R]`. Thus, ``d`` can be formed locally if the distribution can
      be aligned with that of the ``offset`` diagonal of :math:`A[M_C,M_R]`. 

   .. cpp:function:: void GetDiagonal( DistMatrix<T,STAR,MD>& d, int offset=0 ) const

      This is the same as above, but ``d`` is a row-vector instead of a 
      column-vector.

   .. cpp:function:: void SetDiagonal( const DistMatrix<T,MD,STAR>& d, int offset=0 )

      Same as ``GetDiagonal``, but in reverse.

   .. cpp:function:: void SetDiagonal( const DistMatrix<T,STAR,MD>& d, int offset=0 )

      Same as ``GetDiagonal``, but in reverse.

   .. note:: 

      The following are only valid for complex datatypes and are analogous to
      their general counterparts from above in the obvious manner.

   .. cpp:function:: void GetRealDiagonal( DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 ) const

   .. cpp:function:: void GetImagDiagonal( DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 ) const

   .. cpp:function:: void GetRealDiagonal( DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 ) const

   .. cpp:function:: void GetImagDiagonal( DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 ) const

   .. cpp:function:: void SetRealDiagonal( const DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 )

   .. cpp:function:: void SetImagDiagonal( const DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 )

   .. cpp:function:: void SetRealDiagonal( const DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 )

   .. cpp:function:: void SetImagDiagonal( const DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 )

   .. rubric:: Alignment

   All of the following clear the distributed matrix's contents and then 
   reconfigure the alignments as described.

   .. cpp:function:: void Align( int colAlignment, int rowAlignment )

      Specify the process row, ``colAlignment``, and process column,
      ``rowAlignment``, which own the top-left entry.

   .. cpp:function:: void AlignCols( int colAlignment )

      Specify the process row which owns the top-left entry.

   .. cpp:function:: void AlignRows( int rowAlignment )

      Specify the process column which owns the top-left entry.

   .. cpp:function:: void AlignWith( const DistMatrix<S,MC,MR>& A )

      Force the alignments to match those of ``A``.

   .. cpp:function:: void AlignWith( const DistMatrix<S,MC,STAR>& A )

      Force the column alignment to match that of ``A``.

   .. cpp:function:: void AlignWith( const DistMatrix<S,STAR,MR>& A )

      Force the row alignment to match that of ``A``.

   .. cpp:function:: void AlignWith( const DistMatrix<S,MR,MC>& A )

      Force the column alignment to match the row alignment of ``A`` (and 
      vice-versa).

   .. cpp:function:: void AlignWith( const DistMatrix<S,MR,STAR>& A )

      Force the row alignment to match the column alignment of ``A``.

   .. cpp:function:: void AlignWith( const DistMatrix<S,STAR,MC>& A )

      Force the column alignment to match the row alignment of ``A``.

   .. cpp:function:: void AlignWith( const DistMatrix<S,VC,STAR>& A )

      Force the column alignment to be equal to that of ``A`` (modulo 
      the number of process rows).

   .. cpp:function:: void AlignWith( const DistMatrix<S,STAR,VC>& A )

      Force the column alignment to equal the row alignment of ``A`` (modulo
      the number of process rows).

   .. cpp:function:: void AlignWith( const DistMatrix<S,VR,STAR>& A )

      Force the row alignment to equal the column alignment of ``A`` (modulo
      the number of process columns).

   .. cpp:function:: void AlignWith( const DistMatrix<S,STAR,VR>& A )

      Force the row alignment to equal the row alignment of ``A`` (modulo
      the number of process columns).

   .. cpp:function:: void AlignColsWith( const DistMatrix<S,MC,MR>& A )

      Force the column alignment to match that of ``A``.

   .. cpp:function:: void AlignColsWith( const DistMatrix<S,MC,STAR>& A )

      Force the column alignment to match that of ``A``.

   .. cpp:function:: void AlignColsWith( const DistMatrix<S,MR,MC>& A )

      Force the column alignment to match the row alignment of ``A``.

   .. cpp:function:: void AlignColsWith( const DistMatrix<S,STAR,MC>& A )

      Force the column alignment to match the row alignment of ``A``.

   .. cpp:function:: void AlignColsWith( const DistMatrix<S,VC,STAR>& A )

      Force the column alignment to match the column alignment of ``A`` 
      (modulo the number of process rows).

   .. cpp:function:: void AlignColsWith( const DistMatrix<S,STAR,VC>& A )

      Force the column alignment to match the row alignment of ``A`` 
      (modulo the number of process rows).

   .. cpp:function:: void AlignRowsWith( const DistMatrix<S,MC,MR>& A )

      Force the row alignment to match that of ``A``.

   .. cpp:function:: void AlignRowsWith( const DistMatrix<S,STAR,MR>& A )

      Force the row alignment to match that of ``A``.

   .. cpp:function:: void AlignRowsWith( const DistMatrix<S,MR,MC>& A )

      Force the row alignment to match the column alignment of ``A``.

   .. cpp:function:: void AlignRowsWith( const DistMatrix<S,MR,STAR>& A )

      Force the row alignment to match the column alignment of ``A``.

   .. cpp:function:: void AlignRowsWith( const DistMatrix<S,VR,STAR>& A )

      Force the row alignment to match the column alignment of ``A`` (modulo
      the number of process columns).

   .. cpp:function:: void AlignRowsWith( const DistMatrix<S,STAR,VR>& A )

      Force the row alignment to match the row alignment of ``A`` (modulo
      the number of process columns).

   .. rubric:: Views

   .. cpp:function:: void View( DistMatrix<T,MC,MR>& A )

      Reconfigure this matrix such that it is essentially a copy of the 
      distributed matrix ``A``, but the local data buffer simply points to 
      the one from ``A``.

   .. cpp:function:: void LockedView( const DistMatrix<T,MC,MR>& A )

      Same as above, but this matrix is "locked", meaning that it cannot 
      change the data from ``A`` that it points to.

   .. cpp:function:: void View( DistMatrix<T,MC,MR>& A, int i, int j, int height, int width )

      View a subset of ``A`` rather than the entire matrix. In particular, 
      reconfigure this matrix to behave like the submatrix defined from the 
      ``[i,i+height)`` rows and ``[j,j+width)`` columns of ``A``.

   .. cpp:function:: void LockedView( const DistMatrix<T,MC,MR>& A, int i, int j, int height, int width )

      Same as above, but this matrix is "locked", meaning that it cannot
      change the data from ``A`` that it points to.

   .. cpp:function:: void View( int height, int width, int colAlignment, int rowAlignment, T* buffer, int ldim, const elemental::Grid& grid )

      Reconfigure this distributed matrix around an implicit ``[M_C,M_R]`` 
      distributed matrix of the specified dimensions, alignments, local buffer, 
      local leading dimension, and process grid.

   .. cpp:function:: void LockedView( int height, int width, int colAlignment, int rowAlignment, const T* buffer, int ldim, const elemental::Grid& grid )

      Same as above, but the resulting matrix is "locked", meaning that it 
      cannot modify the underlying local data.

   .. note::

      The following functions have strict requirements on the input matrices 
      and must be used with care in ``PureRelease`` and ``HybridRelease`` modes.

   .. cpp:function:: void View1x2( DistMatrix<T,MC,MR>& AL, DistMatrix<T,MC,MR>& AR )

      Recombine two adjacent submatrices to form :math:`[A_L A_R]`. 

   .. cpp:function:: void LockedView1x2( const DistMatrix<T,MC,MR>& AL, const DistMatrix<T,MC,MR>& AR )

      Same as above, but the result is "locked" (the data is not modifiable).

   .. cpp:function:: void View2x1( DistMatrix<T,MC,MR>& AT, DistMatrix<T,MC,MR>& AB )

      Recombine two adjacent submatrices to form :math:`[A_T; A_B]`.

   .. cpp:function:: void LockedView2x1( const DistMatrix<T,MC,MR>& AT, const DistMatrix<T,MC,MR>& AB )

      Same as above, but the result is "locked" (the data is not modifiable).

   .. cpp:function:: void View2x2( DistMatrix<T,MC,MR>& ATL, DistMatrix<T,MC,MR>& ATR, DistMatrix<T,MC,MR>& ABL, DistMatrix<T,MC,MR>& ABR )

      Recombine four adjacent submatrices to form 
      :math:`[A_{TL} A_{TR}; A_{BL} A_{BR}]`.

   .. cpp:function:: void LockedView2x2( const DistMatrix<T,MC,MR>& ATL, const DistMatrix<T,MC,MR>& ATR, const DistMatrix<T,MC,MR>& ABL, const DistMatrix<T,MC,MR>& ABR )

      Same as above, but the result is "locked" (the data is not modifiable).

   .. rubric:: Custom communication routines

   The following routines primarily exist as a means of avoiding the poor 
   memory bandwidth which results from packing or unpacking large amounts of 
   data without a unit stride. PLAPACK noticed this issue and avoided the 
   problem by carefully (conjugate-)transposing matrices in strategic places,
   usually before a gather or after a scatter, and we follow suit.

   .. cpp:function:: void SumScatterFrom( const DistMatrix<T,MC,STAR>& A )

      Simultaneously sum :math:`A[M_C,\star]` within each process row and scatter 
      the entries in each row to form the result in an :math:`[M_C,M_R]` 
      distribution.

   .. cpp:function:: void SumScatterUpdate( T alpha, const DistMatrix<T,MC,STAR>& A )

      Same as above, but add :math:`\alpha` times the result onto the parent
      distributed matrix rather than simply assigning the result to it.

   .. cpp:function:: void SumScatterFrom( const DistMatrix<T,STAR,MR>& A )

      Simultaenously sum :math:`A[\star,M_R]` within each process column and 
      scatter the entries in each column to form the result in an 
      :math:`[M_C,M_R]` distribution.

   .. cpp:function:: void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,MR>& A )

      Same as above, but add :math:`\alpha` times the result onto the parent
      distributed matrix rather than simply assigning the result to it.

   .. cpp:function:: void SumScatterFrom( const DistMatrix<T,STAR,STAR>& A )

      Simultaneously sum :math:`A[\star,\star]` over the entire process grid and 
      scatter the entries in each row and column to form the result in an 
      :math:`[M_C,M_R]` distribution.

   .. cpp:function:: void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR>& A )

      Same as above, but add :math:`\alpha` times the result onto the parent
      distributed matrix rather than simply assigning the result to it.

   .. cpp:function:: void AdjointFrom( const DistMatrix<T,STAR,MC>& A )

      Set the parent matrix equal to the redistributed adjoint of 
      :math:`A[\star,M_C]`; in particular, 
      :math:`(A[\star,M_C])^H = A^H[M_C,\star]`, so perform an 
      :math:`[M_C,M_R] \leftarrow [M_C,\star]` redistribution on the adjoint of
      ``A``, which typically just consists of locally copying (and conjugating) 
      subsets of the data from :math:`A[\star,M_C]`.

   .. cpp:function:: void AdjointFrom( const DistMatrix<T,MR,STAR>& A )

      This routine is the dual of the above routine, and performs an
      :math:`[M_C,M_R] \leftarrow [\star,M_R]` redistribution on the adjoint of 
      ``A``.

   .. cpp:function:: void TransposeFrom( const DistMatrix<T,STAR,MC>& A )

      Same as the corresponding ``AdjointFrom``, but with no conjugation.

   .. cpp:function:: void TransposeFrom( const DistMatrix<T,MR,STAR>& A )

      Same as the corresponding ``AdjointFrom``, but with no conjugation.

``[MC,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[M_C,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[* ,MR]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,M_R]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[MR,MC]``
-----------
**TODO**, but not as high of a priority since the :math:`[M_R,M_C]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[MR,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[M_R,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[* ,MC]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,M_C]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[MD,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[M_D,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[* ,MD]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,M_D]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[VC,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[V_C,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[* ,VC]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,V_C]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[VR,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[V_R,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[* ,VR]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,V_R]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

``[* ,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

Partitioning
============
The following routines are slight tweaks of the FLAME project's 
(as well as PLAPACK's) approach to submatrix tracking; the difference is that 
they have "locked" versions, which are meant for creating partitionings where 
the submatrices cannot be modified.

.. cpp:function:: void PartitionUp( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, int heightAB=Blocksize() )

   **Left off here**

The Axpy interface
==================
The Axpy interface is Elemental's version of the PLAPACK Axpy interface, where 
*Axpy*  is derived from the BLAS shorthand for :math:`Y := \alpha X + Y` 
(Alpha X Plus Y). Rather than always requiring users to manually fill their 
distributed matrix, this interface provides a mechanism so that individual processes
can independently submit local submatrices which will be automatically redistributed 
and added onto the global distributed matrix 
(this would be ``LOCAL_TO_GLOBAL`` mode). The interface also allows for the reverse: 
each process may asynchronously request arbitrary subset of the global distributed 
matrix (``GLOBAL_TO_LOCAL`` mode).

.. note:: 
   
   The catch is that, in order for this behavior to be possible, all of the 
   processes that share a particular distributed matrix must synchronize at the 
   beginning and end of the Axpy interface usage (these synchronizations correspond 
   to the ``Attach`` and ``Detach`` member functions). The distributed matrix 
   should **not** be manually modified between the ``Attach`` and ``Detach`` calls.

An example usage might be:

.. code-block:: cpp

   #include "elemental.hpp"
   using namespace elemental;
   ...
   // Create an 8 x 8 distributed matrix over the given grid
   DistMatrix<double,MC,MR> A( 8, 8, grid );

   // Set every entry of A to zero
   A.SetToZero();

   // Print the original A
   A.Print("Original distributed A");

   // Open up a LOCAL_TO_GLOBAL interface to A 
   AxpyInterface<double> interface;
   interface.Attach( LOCAL_TO_GLOBAL, A );

   // If we are process 0, then create a 3 x 3 identity matrix, B,
   // and Axpy it into the bottom-right of A (using alpha=2)
   // NOTE: The bottom-right 3 x 3 submatrix starts at the (5,5) 
   //       entry of A.
   // NOTE: Every process is free to Axpy as many submatrices as they 
   //       desire at this point.
   if( grid.VCRank() == 0 )
   {
       Matrix<double> B( 3, 3 );
       B.SetToIdentity();
       interface.Axpy( 2.0, B, 5, 5 );
   }

   // Have all processes collectively detach from A
   interface.Detach();

   // Print the updated A
   A.Print("Updated distributed A");

   // Reattach to A, but in the GLOBAL_TO_LOCAL direction
   interface.Attach( GLOBAL_TO_LOCAL, A );

   // Have process 0 request a copy of the entire distributed matrix
   //
   // NOTE: Every process is free to Axpy as many submatrices as they 
   //       desire at this point.
   Matrix<double> C;
   if( grid.VCRank() == 0 )
   {
       C.ResizeTo( 8, 8 );
       C.SetToZero();
       interface.Axpy( 1.0, C, 0, 0 );
   }

   // Collectively detach in order to finish filling process 0's request
   interface.Detach();
   
   // Process 0 can now locally print its copy of A
   if( g.VCRank() == 0 )
       C.Print("Process 0's local copy of A");

The output would be ::

    Original distributed A
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0

    Updated distributed A
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 2 0 0
    0 0 0 0 0 0 2 0
    0 0 0 0 0 0 0 2

    Process 0's local copy of A
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0
    0 0 0 0 0 2 0 0
    0 0 0 0 0 0 2 0
    0 0 0 0 0 0 0 2

.. cpp:type:: enum AxpyType

   Can take on the value of either ``LOCAL_TO_GLOBAL`` or ``GLOBAL_TO_LOCAL``, 
   with the meanings described above.

.. cpp:class:: AxpyInterface<T>

   .. cpp:function:: AxpyInterface()

      Initialize a blank instance of the interface class. It will need to later be 
      attached to a distributed matrix before any Axpy's can occur.

   .. cpp:function:: AxpyInterface( AxpyType type, DistMatrix<T,MC,MR>& Z )

      Initialize an interface to the distributed matrix ``Z``, where ``type`` 
      can be either ``LOCAL_TO_GLOBAL`` or ``GLOBAL_TO_LOCAL``.

   .. cpp:function:: AxpyInterface( AxpyType type, const DistMatrix<T,MC,MR>& Z )

      Initialize an interface to the (unmodifiable) distributed matrix ``Z``; 
      since ``Z`` cannot be modified, the only sensical ``AxpyType`` is 
      ``GLOBAL_TO_LOCAL``. The ``AxpyType`` argument was kept in order to be 
      consistent with the previous routine.

   .. cpp:function:: void Attach( AxpyType type, DistMatrix<T,MC,MR>& Z )

      Attach to the distributed matrix ``Z``, where ``type`` can be either 
      ``LOCAL_TO_GLOBAL`` or ``GLOBAL_TO_LOCAL``.

   .. cpp:function:: void Attach( AxpyType type, const DistMatrix<T,MC,MR>& Z )

      Attach to the (unmodifiable) distributed matrix ``Z``; as mentioned above, 
      the only sensical value of ``type`` is ``GLOBAL_TO_LOCAL``, but the
      ``AxpyType`` argument was kept for consistency.

   .. cpp:function:: void Axpy( T alpha, Matrix<T>& Z, int i, int j )

      If the interface was previously attached in the ``LOCAL_TO_GLOBAL`` 
      direction, then the matrix ``\alpha Z`` will be added onto the associated
      distributed matrix starting at the :math:`(i,j)` global index; otherwise 
      :math:`\alpha` times the submatrix of the associated distributed matrix,
      which starts at index :math:`(i,j)` and is of the same size as ``Z``, will 
      be added onto ``Z``.

   .. cpp:function:: void Axpy( T alpha, const Matrix<T>& Z, int i, int j )

      Same as above, but since ``Z`` is unmodifiable, the attachment must have 
      been in the ``LOCAL_TO_GLOBAL`` direction.

   .. cpp:function:: void Detach()

      All processes collectively finish handling each others requests and then 
      detach from the associated distributed matrix.

Environment routines
====================
**TODO**
