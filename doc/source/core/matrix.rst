The Matrix class
================
This is the basic building block of the library: its purpose it to provide 
convenient mechanisms for performing basic matrix manipulation operations, 
such as setting and querying individual matrix entries, without giving up 
compatibility with interfaces such as BLAS and LAPACK, which assume column-major
storage.

An example of generating an :math:`m \times n` matrix of real double-precision 
numbers where the :math:`(i,j)` entry is equal to :math:`i-j` would be:

  .. code-block:: cpp

     #include "elemental.hpp"
     using namespace elem;
     ...
     Matrix<double> A( m, n );
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             A.Set( i, j, double(i-j) );
     
The underlying data storage is simply a contiguous buffer that stores entries 
in a column-major fashion with an arbitrary leading dimension. For modifiable
instances of the :cpp:type:`Matrix\<T>` class, the routine
:cpp:type:`Matrix\<T>::Buffer` returns a pointer to the underlying 
buffer, while :cpp:type:`Matrix\<T>::LDim` returns the leading 
dimension; these two routines could be used to directly perform the equivalent
of the first code sample as follows:

  .. code-block:: cpp
     
     #include "elemental.hpp"
     using namespace elem;
     ...
     Matrix<double> A( m, n );
     double* buffer = A.Buffer();
     const int ldim = A.LDim();
     for( int j=0; j<n; ++j )
         for( int i=0; i<m; ++i )
             buffer[i+j*ldim] = double(i-j);

For constant instances of the :cpp:type:`Matrix\<T>` class, a ``const`` pointer
to the underlying data can similarly be returned with a call to 
:cpp:func:`Matrix\<T>::LockedBuffer`.
In addition, a (``const``) pointer to the place in the 
(``const``) buffer where entry :math:`(i,j)` resides can be easily retrieved
with a call to :cpp:func:`Matrix\<T>::Buffer` or 
:cpp:func:`Matrix\<T>::LockedBuffer`.

It is also important to be able to create matrices which are simply *views* 
of existing (sub)matrices. For example, if `A` is a :math:`10 \times 10` 
matrix of complex doubles, then a matrix :math:`A_{BR}` can easily be created 
to view the bottom-right :math:`6 \times 7` submatrix using

  .. code-block:: cpp

     #include "elemental.hpp"
     ...
     auto ABR = View( A, 1, 2, 3, 4 );

since the bottom-right :math:`3 \times 4` submatrix beings at index 
:math:`(1,2)`. In general, to view the :math:`M \times N` submatrix starting
at entry :math:`(i,j)`, one would call ``View( ABR, A, i, j, M, N );``.

.. cpp:type:: class Matrix<T>

   The most general case, where the underlying datatype `T` is only assumed to 
   be a ring; that is, it supports multiplication and addition and has the 
   appropriate identities.

   .. rubric:: Constructors

   .. cpp:function:: Matrix()

      This simply creates a default :math:`0 \times 0` matrix with a leading 
      dimension of one (BLAS and LAPACK require positive leading dimensions).

   .. cpp:function:: Matrix( int height, int width )

      A `height` :math:`\times` `width` matrix is created with an unspecified
      leading dimension (though it is currently implemented as 
      ``std::max(height,1)``).

   .. cpp:function:: Matrix( int height, int width, int ldim )

      A `height` :math:`\times` `width` matrix is created with a leading 
      dimension equal to `ldim` (which must be greater than or equal 
      ``std::min(height,1)``).

   .. cpp:function:: Matrix( int height, int width, const T* buffer, int ldim )

      A matrix is built around column-major constant buffer ``const T* buffer`` 
      with the specified dimensions. The memory pointed to by `buffer` should
      not be freed until after the :cpp:type:`Matrix\<T>` object is destructed.

   .. cpp:function:: Matrix( int height, int width, T* buffer, int ldim )

      A matrix is built around the column-major modifiable buffer ``T* buffer``
      with the specified dimensions. The memory pointed to by `buffer` should
      not be freed until after the :cpp:type:`Matrix\<T>` object is destructed.

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

      Return the number of entries of type `T` that this :cpp:type:`Matrix\<T>`
      instance has allocated space for.

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

   .. rubric:: Entry manipulation

   .. cpp:function:: T Get( int i, int j ) const

      Return entry :math:`(i,j)`.

   .. cpp:function:: void Set( int i, int j, T alpha )

      Set entry :math:`(i,j)` to :math:`\alpha`.

   .. cpp:function:: void Update( int i, int j, T alpha )

      Add :math:`\alpha` to entry :math:`(i,j)`.

   .. cpp:function:: void GetDiagonal( Matrix<T>& d, int offset=0 ) const

      Modify :math:`d` into a column-vector containing the entries lying on the 
      `offset` diagonal of our matrix (for instance, the main diagonal has 
      offset :math:`0`, the subdiagonal has offset :math:`-1`, and the 
      superdiagonal has offset :math:`+1`).

   .. cpp:function:: Matrix<T> GetDiagonal( Int offset=0 ) const

      Efficiently construct and return the diagonal via C++11 move semantics.

   .. cpp:function:: void SetDiagonal( const Matrix<T>& d, int offset=0 )

      Set the entries in the `offset` diagonal entries from the contents of 
      the column-vector :math:`d`.

   .. cpp:function:: void UpdateDiagonal( const Matrix<T>& d, int offset=0 )

      Add the contents of :math:`d` onto the entries in the `offset` diagonal.

   .. note::

      Many of the following routines are only valid for complex datatypes.

   .. cpp:function:: typename Base<T>::type GetRealPart( int i, int j ) const
   .. cpp:function:: typename Base<T>::type GetImagPart( int i, int j ) const

      Return the real (imaginary) part of entry :math:`(i,j)`.

   .. cpp:function:: void SetRealPart( int i, int j, typename Base<T>::type alpha )
   .. cpp:function:: void SetImagPart( int i, int j, typename Base<T>::type alpha )

      Set the real (imaginary) part of entry :math:`(i,j)` to :math:`\alpha`.

   .. cpp:function:: void UpdateRealPart( int i, int j, typename Base<T>::type alpha )
   .. cpp:function:: void UpdateImagPart( int i, int j, typename Base<T>::type alpha ) 

      Add :math:`\alpha` to the real (imaginary) part of entry :math:`(i,j)`.

   .. cpp:function:: void GetRealPartOfDiagonal( Matrix<typename Base<T>::type>& d, int offset=0 ) const
   .. cpp:function:: void GetImagPartOfDiagonal( Matrix<typename Base<T>::type>& d, int offset=0 ) const

      Modify :math:`d` into a column-vector containing the real (imaginary) 
      parts of the entries in the `offset` diagonal.

   .. cpp:function:: Matrix<typename Base<T>::type> GetRealPartOfDiagonal( Int offset=0 ) const
   .. cpp:function:: Matrix<typename Base<T>::type> GetRealPartOfDiagonal( Int offset=0 ) const

      Use C++11 move semantics to return the real (imaginary) part of the 
      diagonal.

   .. cpp:function:: void SetRealPartOfDiagonal( const Matrix<typename Base<T>::type>& d, int offset=0 )
   .. cpp:function:: void SetImagPartOfDiagonal( const Matrix<typename Base<T>::type>& d, int offset=0 )

      Set the real (imaginary) parts of the entries in the `offset` diagonal 
      from the contents of the column-vector :math:`d`.

   .. cpp:function:: void UpdateRealPartOfDiagonal( const Matrix<typename Base<T>::type>& d, int offset=0 )
   .. cpp:function:: void UpdateImagPartOfDiagonal( const Matrix<typename Base<T>::type>& d, int offset=0 )

      Add the contents of the column-vector :math:`d` onto the real (imaginary)
      parts of the entries in the `offset` diagonal.

   .. rubric:: Views

   .. cpp:function:: bool Viewing() const

      Return whether or not this matrix is currently viewing another matrix.

   .. cpp:function:: bool Locked() const

      Return whether or not we can modify the data we are viewing.

   .. cpp:function:: void Attach( int height, int width, T* buffer, int ldim )

      Reconfigure the matrix around the specified buffer.

   .. cpp:function:: void LockedAttach( int height, int width, const T* buffer, int ldim )

      Reconfigure the matrix around the specified unmodifiable buffer.

   .. rubric:: Utilities

   .. cpp:function:: const Matrix<T>& operator=( const Matrix<T>& A )

      Create a copy of matrix :math:`A`.

   .. cpp:function:: void Empty()

      Sets the matrix to :math:`0 \times 0` and frees the underlying buffer.

   .. cpp:function:: void ResizeTo( int height, int width )

      Reconfigures the matrix to be `height` :math:`\times` `width`.

   .. cpp:function:: void ResizeTo( int height, int width, int ldim )

      Reconfigures the matrix to be `height` :math:`\times` `width`, but with 
      leading dimension equal to `ldim` (which must be greater than or equal to 
      ``std::min(height,1)``).

Special cases used in Elemental
-------------------------------
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`Matrix\<T>`.

.. cpp:type:: class Matrix<R>

   Used to denote that the underlying datatype `R` is real.

.. cpp:type:: class Matrix<Complex<R> >

   Used to denote that the underlying datatype :cpp:type:`Complex\<R>` is
   complex with base type `R`.

.. cpp:type:: class Matrix<F>

   Used to denote that the underlying datatype `F` is a field.

