The DistMatrix class
====================
The :cpp:type:`DistMatrix\<T,U,V>` class is meant to provide a 
distributed-memory analogue of the :cpp:type:`Matrix\<T>` class. 
Similarly to PLAPACK, roughly ten different matrix 
distributions are provided and it is trivial (in the programmability sense) to 
redistribute from one to another: in PLAPACK, one would simply call 
``PLA_Copy``, whereas, in Elemental, it is handled through overloading the 
:math:`=` operator.

Since it is crucial to know not only how many 
processes to distribute the data over, but *which* processes, and in what 
manner they should be decomposed into a logical two-dimensional grid, an 
instance of the :cpp:type:`Grid` class must be passed into the constructor of 
the :cpp:type:`DistMatrix\<T,U,V>` class.

.. note:: 
   
   Since the :cpp:type:`DistMatrix\<T,U,V>` class makes use of MPI for 
   message passing, custom interfaces must be written for nonstandard datatypes.
   As of now, the following datatypes are fully supported for 
   :cpp:type:`DistMatrix\<T,U,V>`:
   ``int``, ``float``, ``double``, ``Complex<float>``, and ``Complex<double>``.

.. cpp:type:: struct DistData

   .. cpp:member:: Distribution colDist
   
   .. cpp:member:: Distribution rowDist

   .. cpp:member:: int colAlignment

   .. cpp:member:: int rowAlignment

   .. cpp:member:: int diagPath

   .. cpp:member:: const Grid* grid

AbstractDistMatrix
------------------

This abstract class defines the list of member functions that are guaranteed 
to be available for all matrix distributions.

.. cpp:type:: class AbstractDistMatrix<T>

   The most general case, where the underlying datatype `T` is only assumed 
   to be a ring; that is, it supports multiplication and addition and has the 
   appropriate identities.

   .. rubric:: Basic information

   .. cpp:function:: int Height() const

      Return the height of the matrix.

   .. cpp:function:: int Width() const

      Return the width of the matrix.

   .. cpp:function:: int LocalHeight() const

      Return the local height of the matrix.

   .. cpp:function:: int LocalWidth() const

      Return the local width of the matrix.

   .. cpp:function:: int LDim() const

      Return the local leading dimension of the matrix.

   .. cpp:function:: size_t AllocatedMemory() const

      Return the number of entries of type `T` that we have locally allocated
      space for.

   .. cpp:function:: const Grid& Grid() const

      Return the grid that this distributed matrix is distributed over.

   .. cpp:function:: T* Buffer( int iLocal=0, int jLocal=0 )

      Return a pointer to the portion of the local buffer that stores entry 
      `(iLocal,jLocal)`.

   .. cpp:function:: const T* LockedBuffer( int iLocal=0, int jLocal=0 ) const

      Return a pointer to the portion of the local buffer that stores entry
      `(iLocal,jLocal)`, but do not allow for the data to be modified through
      the returned pointer.

   .. cpp:function:: Matrix<T>& Matrix()

      Return a reference to the local matrix.

   .. cpp:function:: const Matrix<T>& LockedMatrix() const

      Return an unmodifiable reference to the local matrix.

   .. rubric:: Distribution details

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

   .. cpp:function:: int ColStride() const

      Return the number of rows between locally owned entries.

   .. cpp:function:: int RowStride() const

      Return the number of columns between locally owned entries.

   .. cpp:function:: elem::DistData DistData() const

      Returns a description of the distribution and alignment information

   .. rubric:: Entry manipulation

   .. cpp:function:: T Get( int i, int j ) const

      Return the `(i,j)` entry of the global matrix. This operation is 
      collective.

   .. cpp:function:: void Set( int i, int j, T alpha )

      Set the `(i,j)` entry of the global matrix to :math:`\alpha`. This 
      operation is collective.

   .. cpp:function:: void Update( int i, int j, T alpha )

      Add :math:`\alpha` to the `(i,j)` entry of the global matrix. This 
      operation is collective.

   .. cpp:function:: T GetLocal( int iLocal, int jLocal ) const

      Return the `(iLocal,jLocal)` entry of our local matrix.

   .. cpp:function:: void SetLocal( int iLocal, int jLocal, T alpha )

      Set the `(iLocal,jLocal)` entry of our local matrix to :math:`\alpha`.

   .. cpp:function:: void UpdateLocal( int iLoca, int jLocal, T alpha )

      Add :math:`\alpha` to the `(iLocal,jLocal)` entry of our local matrix.

   .. note::

      Many of the following routines are only valid for complex datatypes.

   .. cpp:function:: typename Base<T>::type GetRealPart( int i, int j ) const
   .. cpp:function:: typename Base<T>::type GetImagPart( int i, int j ) const

      Return the real (imaginary) part of the `(i,j)` entry of the global 
      matrix. This operation is collective.

   .. cpp:function:: void SetRealPart( int i, int j, typename Base<T>::type alpha )
   .. cpp:function:: void SetImagPart( int i, int j, typename Base<T>::type alpha )

      Set the real (imaginary) part of the `(i,j)` entry of the global matrix to
      :math:`\alpha`.

   .. cpp:function:: void UpdateRealPart( int i, int j, typename Base<T>::type alpha )
   .. cpp:function:: void UpdateImagPart( int i, int j, typename Base<T>::type alpha )

      Add :math:`\alpha` to the real (imaginary) part of the `(i,j)` entry of 
      the global matrix.

   .. cpp:function:: typename Base<T>::type GetRealPartLocal( int iLocal, int jLocal ) const
   .. cpp:function:: typename Base<T>::type GetLocalImagPart( int iLocal, int jLocal ) const

      Return the real (imaginary) part of the `(iLocal,jLocal)` entry of our 
      local matrix.

   .. cpp:function:: void SetLocalRealPart( int iLocal, int jLocal, typename Base<T>::type alpha )
   .. cpp:function:: void SetLocalImagPart( int iLocal, int jLocal, typename Base<T>::type alpha )

      Set the real (imaginary) part of the `(iLocal,jLocal)` entry of our local 
      matrix.

   .. cpp:function:: void UpdateRealPartLocal( int iLocal, int jLocal, typename Base<T>::type alpha )
   .. cpp:function:: void UpdateLocalImagPart( int iLocal, int jLocal, typename Base<T>::type alpha )

      Add :math:`\alpha` to the real (imaginary) part of the `(iLocal,jLocal)` 
      entry of our local matrix.

   .. rubric:: Viewing

   .. cpp:function:: bool Viewing() const

      Return whether or not this matrix is viewing another.

   .. cpp:function:: bool Locked() const

      Return whether or not this matrix is viewing another in a manner
      that does not allow for modifying the viewed data.

   .. rubric:: Utilities

   .. cpp:function:: void Empty()

      Resize the distributed matrix so that it is :math:`0 \times 0` and free 
      all allocated storage.

   .. cpp:function:: void ResizeTo( int height, int width )

      Reconfigure the matrix so that it is `height` :math:`\times` `width`.

   .. cpp:function:: void ResizeTo( int height, int width, int ldim )

      Same as above, but the local leading dimension is also specified.

   .. cpp:function:: void SetGrid( const Grid& grid )

      Clear the distributed matrix's contents and reconfigure for the new 
      process grid.

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`AbstractDistMatrix\<T>`.

.. cpp:type:: class AbstractDistMatrix<R>

   Used to denote that the underlying datatype `R` is real.

.. cpp:type:: class AbstractDistMatrix<Complex<R> >

   Used to denote that the underlying datatype :cpp:type:`Complex\<R>` is 
   complex with base type `R`.

.. cpp:type:: class AbstractDistMatrix<F>

   Used to denote that the underlying datatype `F` is a field. 

DistMatrix
----------

.. cpp:type:: class DistMatrix<T,U,V>

   This templated class for manipulating distributed matrices is only defined
   for the following choices of the column and row :cpp:type:`Distribution`'s, 
   `U` and `V` (`T` is a ring in this case).

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special 
cases of :cpp:type:`DistMatrix\<T,U,V>`.

.. cpp:type:: class DistMatrix<double,U,V>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,U,V>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,U,V>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,U,V>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`. 

.. cpp:type:: class DistMatrix<F,U,V>

   The underlying datatype `F` is a field.

``[MC,MR]``
-----------

This is by far the most important matrix distribution in Elemental, as the vast
majority of parallel routines expect the input to be in this form. For a
:math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process grid,
individual entries would be owned by the following processes (assuming the 
column and row alignments are both 0):

.. math::

   \left(\begin{array}{ccccccc}
     0 & 2 & 4 & 0 & 2 & 4 & 0 \\
     1 & 3 & 5 & 1 & 3 & 5 & 1 \\ 
     0 & 2 & 4 & 0 & 2 & 4 & 0 \\
     1 & 3 & 5 & 1 & 3 & 5 & 1 \\ 
     0 & 2 & 4 & 0 & 2 & 4 & 0 \\
     1 & 3 & 5 & 1 & 3 & 5 & 1 \\ 
     0 & 2 & 4 & 0 & 2 & 4 & 0  
   \end{array}\right)

Similarly, if the column alignment is kept at 0 and the row alignment is changed
to 2 (meaning that the third process column owns the first column of the 
matrix), the individual entries would be owned as follows:

.. math::

   \left(\begin{array}{ccccccc}
     4 & 0 & 2 & 4 & 0 & 2 & 4 \\
     5 & 1 & 3 & 5 & 1 & 3 & 5 \\ 
     4 & 0 & 2 & 4 & 0 & 2 & 4 \\
     5 & 1 & 3 & 5 & 1 & 3 & 5 \\ 
     4 & 0 & 2 & 4 & 0 & 2 & 4 \\
     5 & 1 & 3 & 5 & 1 & 3 & 5 \\ 
     4 & 0 & 2 & 4 & 0 & 2 & 4 
   \end{array}\right)

It should also be noted that this is the default distribution format for the 
:cpp:type:`DistMatrix\<T,U,V>` class, as :cpp:type:`DistMatrix\<T>` defaults to
:cpp:type:`DistMatrix\<T,MC,MR>`.

.. cpp:type:: class DistMatrix<T>

.. cpp:type:: class DistMatrix<T,MC,MR>

   The most general case, where the underlying datatype `T` is only assumed 
   to be a ring.

   .. rubric:: Constructors

   .. cpp:function:: DistMatrix( const Grid& grid=DefaultGrid() )
      
      Create a :math:`0 \times 0` distributed matrix over the specified grid.

   .. cpp:function:: DistMatrix( int height, int width, const Grid& grid=DefaultGrid() )

      Create a `height` :math:`\times` `width` distributed matrix over the
      specified grid.

   .. cpp:function:: DistMatrix( int height, int width, int colAlignment, int rowAlignment, const Grid& grid )

      Create a `height` :math:`\times` `width` distributed matrix 
      distributed over the specified process grid, but with the top-left entry
      owned by the `colAlignment` process row and the `rowAlignment` 
      process column.

   .. cpp:function:: DistMatrix( int height, int width, int colAlignment, int rowAlignment, int ldim, const Grid& grid )

      Same as above, but the local leading dimension is also specified.

   .. cpp:function:: DistMatrix( int height, int width, int colAlignment, int rowAlignment, const T* buffer, int ldim, const Grid& grid )

      View a constant distributed matrix's buffer; the buffer must correspond 
      to the local portion of an elemental distributed matrix with the 
      specified row and column alignments and leading dimension, `ldim`.

   .. cpp:function:: DistMatrix( int height, int width, int colAlignment, int rowAlignment, T* buffer, int ldim, const Grid& grid )

      Same as above, but the contents of the matrix are modifiable.

   .. cpp:function:: DistMatrix( const DistMatrix<T,U,V>& A )

      Build a copy of the distributed matrix `A`, but force it to be in the
      ``[MC,MR]`` distribution.

   .. rubric:: Redistribution

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MC,MR>& A )

      If this matrix can be properly aligned with `A`, then perform a local
      copy, otherwise perform an :cpp:func:`mpi::SendRecv` permutation first.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,MC,STAR>& A )

      Perform a local (filtered) copy to form an ``[MC,MR ]`` distribution and 
      then, if necessary, fix the alignment of the ``MC`` distribution via an 
      :cpp:func:`mpi::SendRecv` within process columns.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,MR>& A )
       
      Perform a local (filtered) copy to form an ``[MC,MR ]`` distribution and 
      then, if necessary, fix the alignment of the ``MR`` distribution via an 
      :cpp:func:`mpi::SendRecv` within process rows.

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

      This is similar to the above routine, but with each row of `A` being 
      undistributed, and only one approach is needed: 
      :math:`A[M_C,M_R] \leftarrow A[V_C,\star] \leftarrow A[V_R,\star] \leftarrow A[M_R,\star]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,MC>& A )

      This routine is dual to the :math:`A[M_C,M_R] \leftarrow A[M_R,\star]` 
      redistribution and is accomplished through the sequence: 
      :math:`A[M_C,M_R] \leftarrow A[\star,V_R] \leftarrow A[\star,V_C] \leftarrow A[\star,M_C]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,VC,STAR>& A )

      Perform an :cpp:func:`mpi::AllToAll` within process rows in order to 
      redistribute to the ``[MC,MR]`` distribution 
      (an :cpp:func:`mpi::SendRecv` within process columns may be required for 
      alignment).

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,VC>& A )

      Accomplished through the sequence 
      :math:`A[M_C,M_R] \leftarrow A[\star,V_R] \leftarrow A[\star,V_C]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,VR,STAR>& A )

      Accomplished through the sequence
      :math:`A[M_C,M_R] \leftarrow A[V_C,\star] \leftarrow A[V_R,\star]`.

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,VR>& A )

      Perform an :cpp:func:`mpi::AllToAll` within process columns in order to 
      redistribute to the ``[MC,MR]`` distribution 
      (an :cpp:func:`mpi::SendRecv` within process rows may be required for 
      alignment).

   .. cpp:function:: const DistMatrix<T,MC,MR>& operator=( const DistMatrix<T,STAR,STAR>& A )

      Perform an :cpp:func:`mpi::AllGather` over the entire grid in order to 
      give every process a full copy of `A`.

   .. rubric:: Diagonal manipulation

   .. cpp:function:: void GetDiagonal( DistMatrix<T,MD,STAR>& d, int offset=0 ) const
   .. cpp:function:: void GetDiagonal( DistMatrix<T,STAR,MD>& d, int offset=0 ) const

      The :math:`[M_D,\star]` (:math:`[\star,M_D]`) distribution is defined 
      such that its columns (rows) are distributed like diagonals of the 
      standard matrix distribution, ``[MC,MR]``. 
      Thus, `d` can be formed locally if the distribution can
      be aligned with that of the `offset` diagonal of :math:`A[M_C,M_R]`. 

   .. cpp:function:: void SetDiagonal( const DistMatrix<T,MD,STAR>& d, int offset=0 )
   .. cpp:function:: void SetDiagonal( const DistMatrix<T,STAR,MD>& d, int offset=0 )

      Same as :cpp:func:`DistMatrix\<T>::GetDiagonal`, but in reverse.

   .. note:: 

      Many of the following routines are only valid for complex datatypes and
      are analogous to their general counterparts from above in the obvious 
      manner.

   .. cpp:function:: void GetRealPartOfDiagonal( DistMatrix<typename Base<T>::type,MD,STAR>& d, int offset=0 ) const

   .. cpp:function:: void GetImagPartOfDiagonal( DistMatrix<typename Base<T>::type,MD,STAR>& d, int offset=0 ) const

   .. cpp:function:: void GetRealPartOfDiagonal( DistMatrix<typename Base<T>::type,STAR,MD>& d, int offset=0 ) const

   .. cpp:function:: void GetImagPartOfDiagonal( DistMatrix<typename Base<T>::type,STAR,MD>& d, int offset=0 ) const

   .. cpp:function:: DistMatrix<typename Base<T>::type,MD,STAR> GetRealPartOfDiagonal( int offset=0 ) const
   .. cpp:function:: DistMatrix<typename Base<T>::type,MD,STAR> GetImagPartOfDiagonal( int offset=0 ) const

   .. cpp:function:: void SetRealPartOfDiagonal( const DistMatrix<typename Base<T>::type,MD,STAR>& d, int offset=0 )

   .. cpp:function:: void SetImagPartOfDiagonal( const DistMatrix<typename Base<T>::type,MD,STAR>& d, int offset=0 )

   .. cpp:function:: void SetRealPartOfDiagonal( const DistMatrix<typename Base<T>::type,STAR,MD>& d, int offset=0 )

   .. cpp:function:: void SetImagPartOfDiagonal( const DistMatrix<typename Base<T>::type,STAR,MD>& d, int offset=0 )

   .. rubric:: Alignment

   All of the following clear the distributed matrix's contents and then 
   reconfigure the alignments as described.

   .. cpp:function:: void AlignWith( const AbstractDistMatrix<T>& A )

      Force the alignments to match those of `A`.

   .. cpp:function:: void AlignWith( const elem::DistData& data )

      A mechanism for aligning with a distributed matrix of a different 
      datatype, via ``AlignWith( A.DistData() );``

   .. cpp:function:: void AlignColsWith( const AbstractDistMatrix<T>& A )

      Force the column alignment to match that of `A`.

   .. cpp:function:: void AlignColsWith( const elem::DistData& data )

      A mechanism for aligning with a distributed matrix of a different 
      datatype, via ``AlignColsWith( A.DistData() );``

   .. cpp:function:: void AlignRowsWith( const AbstractDistMatrix<T>& A )

      Force the row alignment to match that of `A`.

   .. cpp:function:: void AlignRowsWith( const elem::DistData& data )

      A mechanism for aligning with a distributed matrix of a different 
      datatype, via ``AlignRowsWith( A.DistData() );``

   .. rubric:: Views

   .. cpp:function:: void Attach( int height, int width, int colAlignment, int rowAlignment, T* buffer, int ldim, const Grid& grid )

      Reconfigure this distributed matrix around an implicit ``[MC,MR]`` 
      distributed matrix of the specified dimensions, alignments, local buffer, 
      local leading dimension, and process grid.

   .. cpp:function:: void LockedAttach( int height, int width, int colAlignment, int rowAlignment, const T* buffer, int ldim, const Grid& grid )

      Same as above, but the resulting matrix is "locked", meaning that it 
      cannot modify the underlying local data.

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
      `A`, which typically just consists of locally copying (and conjugating) 
      subsets of the data from :math:`A[\star,M_C]`.

   .. cpp:function:: void AdjointFrom( const DistMatrix<T,MR,STAR>& A )

      This routine is the dual of the above routine, and performs an
      :math:`[M_C,M_R] \leftarrow [\star,M_R]` redistribution on the adjoint of 
      `A`.

   .. cpp:function:: void TransposeFrom( const DistMatrix<T,STAR,MC>& A, bool conjugate=false )

      Same as the corresponding :cpp:func:`DistMatrix\<T>::AdjointFrom`, but 
      with no conjugation by default.

   .. cpp:function:: void TransposeFrom( const DistMatrix<T,MR,STAR>& A, bool conjugate=false )

      Same as the corresponding :cpp:func:`DistMatrix\<T>::AdjointFrom`, but 
      with no conjugation by default.

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special 
cases of :cpp:type:`DistMatrix\<T,MC,MR>` = :cpp:type:`DistMatrix\<T>`.

.. cpp:type:: class DistMatrix<double>

.. cpp:type:: class DistMatrix<double,MC,MR>

   The underlying datatype is the set of double-precision real numbers. 

.. cpp:type:: class DistMatrix<Complex<double>>

.. cpp:type:: class DistMatrix<Complex<double>,MC,MR>

   The underlying datatype is the set of double-precision complex numbers. 

.. cpp:type:: class DistMatrix<R>

.. cpp:type:: class DistMatrix<R,MC,MR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>>

.. cpp:type:: class DistMatrix<Complex<R>,MC,MR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`. 

.. cpp:type:: class DistMatrix<F>

.. cpp:type:: class DistMatrix<F,MC,MR>

   The underlying datatype `F` is a field.

``[MC,* ]``
-----------

This distribution is often used as part of matrix-matrix multiplication. For a
:math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process grid,
individual entries would be owned by the following processes (assuming the 
column alignment is 0):

.. math::

   \left(\begin{array}{ccccccc}
     \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & 
     \{0,2,4\} & \{0,2,4\} \\
     \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & 
     \{1,3,5\} & \{1,3,5\} \\ 
     \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & 
     \{0,2,4\} & \{0,2,4\} \\
     \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & 
     \{1,3,5\} & \{1,3,5\} \\ 
     \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & 
     \{0,2,4\} & \{0,2,4\} \\
     \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & \{1,3,5\} & 
     \{1,3,5\} & \{1,3,5\} \\ 
     \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & \{0,2,4\} & 
     \{0,2,4\} & \{0,2,4\} 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,MC,STAR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,MC,STAR>`.

.. cpp:type:: class DistMatrix<double,MC,STAR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,MC,STAR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,MC,STAR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,MC,STAR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,MC,STAR>

   The underlying datatype `F` is a field.

``[* ,MR]``
-----------
This distribution is also frequently used for matrix-matrix multiplication. 
For a :math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process 
grid, individual entries would be owned by the following processes (assuming 
the row alignment is 0):

.. math::

   \left(\begin{array}{ccccccc}
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} \\
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} \\
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} \\
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} \\
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} \\
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} \\
     \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} & \{2,3\} & \{4,5\} & \{0,1\} 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,STAR,MR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,STAR,MR>`.

.. cpp:type:: class DistMatrix<double,STAR,MR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,STAR,MR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,STAR,MR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,STAR,MR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,STAR,MR>

   The underlying datatype `F` is a field.

``[MR,MC]``
-----------
This is essentially the transpose of the standard matrix distribution, 
``[MC,MR]``. For a
:math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process grid,
individual entries would be owned by the following processes (assuming the 
column and row alignments are both 0):

.. math::

   \left(\begin{array}{ccccccc}
     0 & 1 & 0 & 1 & 0 & 1 & 0 \\
     2 & 3 & 2 & 3 & 2 & 3 & 2 \\
     4 & 5 & 4 & 5 & 4 & 5 & 4 \\
     0 & 1 & 0 & 1 & 0 & 1 & 0 \\
     2 & 3 & 2 & 3 & 2 & 3 & 2 \\
     4 & 5 & 4 & 5 & 4 & 5 & 4 \\
     0 & 1 & 0 & 1 & 0 & 1 & 0 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,MR,MC>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,MR,MC>`.

.. cpp:type:: class DistMatrix<double,MR,MC>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,MR,MC>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,MR,MC>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,MR,MC>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,MR,MC>

   The underlying datatype `F` is a field.
 
``[MR,* ]``
-----------
This is the transpose of the ``[* ,MR]`` distribution and is, like many of 
the previous distributions, useful for matrix-matrix multiplication.
For a :math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process 
grid, individual entries would be owned by the following processes (assuming 
the column alignment is 0):

.. math::

   \left(\begin{array}{ccccccc}
     \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} \\
     \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} \\
     \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} \\
     \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} \\
     \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} & \{2,3\} \\
     \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} & \{4,5\} \\
     \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} & \{0,1\} 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,MR,STAR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,MR,STAR>`.

.. cpp:type:: class DistMatrix<double,MR,STAR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,MR,STAR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,MR,STAR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,MR,STAR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,MR,STAR>

   The underlying datatype `F` is a field.

``[* ,MC]``
-----------
This is the transpose of the ``[MC,*]`` distribution and is, like many of 
the previous distributions, useful for matrix-matrix multiplication.
For a :math:`7 \times 7` matrix distributed over a :math:`2 \times 3` process 
grid, individual entries would be owned by the following processes (assuming 
the column alignment is 0):

.. math::

   \left(\begin{array}{ccccccc}
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} \\
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} \\
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} \\
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} \\
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} \\
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} \\
     \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & \{0,2,4\} & \{1,3,5\} & 
     \{0,2,4\} 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,STAR,MC>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,STAR,MC>`.

.. cpp:type:: class DistMatrix<double,STAR,MC>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,STAR,MC>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,STAR,MC>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,STAR,MC>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,STAR,MC>

   The underlying datatype `F` is a field.

``[MD,* ]``
-----------
**TODO**, but not as high of a priority since the :math:`[M_D,\star]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

.. cpp:type:: class DistMatrix<T,MD,STAR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,MD,STAR>`.

.. cpp:type:: class DistMatrix<double,MD,STAR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,MD,STAR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,MD,STAR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,MD,STAR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,MD,STAR>

   The underlying datatype `F` is a field.

``[* ,MD]``
-----------
**TODO**, but not as high of a priority since the :math:`[\star,M_D]` 
distribution is not as crucial for end users as many other details that have 
not yet been documented.

.. cpp:type:: class DistMatrix<T,STAR,MD>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,STAR,MD>`.

.. cpp:type:: class DistMatrix<double,STAR,MD>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,STAR,MD>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,STAR,MD>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,STAR,MD>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,STAR,MD>

   The underlying datatype `F` is a field.

``[VC,* ]``
-----------
This distribution makes use of a 1d distribution which uses a column-major 
ordering of the entire process grid. Since 1d distributions are useful for 
distributing *vectors*, and a *column-major* ordering is used, the distribution 
symbol is ``VC``. Again using the simple :math:`2 \times 3` process grid, 
with a zero column alignment, each entry of a :math:`7 \times 7` matrix 
would be owned by the following sets of processes:

.. math::

   \left(\begin{array}{ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     1 & 1 & 1 & 1 & 1 & 1 & 1 \\
     2 & 2 & 2 & 2 & 2 & 2 & 2 \\
     3 & 3 & 3 & 3 & 3 & 3 & 3 \\
     4 & 4 & 4 & 4 & 4 & 4 & 4 \\
     5 & 5 & 5 & 5 & 5 & 5 & 5 \\
     0 & 0 & 0 & 0 & 0 & 0 & 0
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,VC,STAR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,VC,STAR>`.

.. cpp:type:: class DistMatrix<double,VC,STAR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,VC,STAR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,VC,STAR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,VC,STAR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,VC,STAR>

   The underlying datatype `F` is a field.

``[* ,VC]``
-----------
This is the transpose of the above ``[VC,* ]`` distribution. On the standard
:math:`2 \times 3` process grid with a row alignment of zero, a 
:math:`7 \times 7` matrix would be distributed as:

.. math::

   \left(\begin{array}{ccccccc}
   0 & 1 & 2 & 3 & 4 & 5 & 0 \\
   0 & 1 & 2 & 3 & 4 & 5 & 0 \\
   0 & 1 & 2 & 3 & 4 & 5 & 0 \\
   0 & 1 & 2 & 3 & 4 & 5 & 0 \\
   0 & 1 & 2 & 3 & 4 & 5 & 0 \\
   0 & 1 & 2 & 3 & 4 & 5 & 0 \\
   0 & 1 & 2 & 3 & 4 & 5 & 0 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,STAR,VC>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,STAR,VC>`.

.. cpp:type:: class DistMatrix<double,STAR,VC>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,STAR,VC>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,STAR,VC>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,STAR,VC>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,STAR,VC>

   The underlying datatype `F` is a field.

``[VR,* ]``
-----------
This distribution makes use of a 1d distribution which uses a row-major 
ordering of the entire process grid. Since 1d distributions are useful for 
distributing *vectors*, and a *row-major* ordering is used, the distribution 
symbol is ``VR``. Again using the simple :math:`2 \times 3` process grid, 
with a zero column alignment, each entry of a :math:`7 \times 7` matrix 
would be owned by the following sets of processes:

.. math::

   \left(\begin{array}{ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     2 & 2 & 2 & 2 & 2 & 2 & 2 \\
     4 & 4 & 4 & 4 & 4 & 4 & 4 \\
     1 & 1 & 1 & 1 & 1 & 1 & 1 \\
     3 & 3 & 3 & 3 & 3 & 3 & 3 \\
     5 & 5 & 5 & 5 & 5 & 5 & 5 \\
     0 & 0 & 0 & 0 & 0 & 0 & 0
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,VR,STAR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,VR,STAR>`.

.. cpp:type:: class DistMatrix<double,VR,STAR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,VR,STAR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,VR,STAR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,VR,STAR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,VR,STAR>

   The underlying datatype `F` is a field.

``[* ,VR]``
-----------
This is the transpose of the above ``[VR,* ]`` distribution. On the standard
:math:`2 \times 3` process grid with a row alignment of zero, a 
:math:`7 \times 7` matrix would be distributed as:

.. math::

   \left(\begin{array}{ccccccc}
   0 & 2 & 4 & 1 & 3 & 5 & 0 \\
   0 & 2 & 4 & 1 & 3 & 5 & 0 \\
   0 & 2 & 4 & 1 & 3 & 5 & 0 \\
   0 & 2 & 4 & 1 & 3 & 5 & 0 \\
   0 & 2 & 4 & 1 & 3 & 5 & 0 \\
   0 & 2 & 4 & 1 & 3 & 5 & 0 \\
   0 & 2 & 4 & 1 & 3 & 5 & 0 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,STAR,VR>

   **TODO:** Add the member functions. 

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,STAR,VR>`.

.. cpp:type:: class DistMatrix<double,STAR,VR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,STAR,VR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,STAR,VR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,STAR,VR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,STAR,VR>

   The underlying datatype `F` is a field.

``[* ,* ]``
-----------
This "distribution" actually redundantly stores every entry of the associated
matrix on every process. Again using a :math:`2 \times 3` process grid, 
the entries of a :math:`7 \times 7` matrix would be owned by the following
sets of processes:

.. math::

   \left(\begin{array}{ccccccc}
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} \\
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} \\
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} \\
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} \\
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} \\
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} \\
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} & 
   \{0,1,...,5\} & \{0,1,...,5\} & \{0,1,...,5\} 
   \end{array}\right)

.. cpp:type:: class DistMatrix<T,STAR,STAR>

   **TODO:** Add the member functions. 

``[o ,o ]``
-----------
This ``distribution`` stores the entire matrix on a single process.

.. cpp:type:: class DistMatrix<T,CIRC,CIRC>

   .. cpp:function:: int Root()

      Returns the rank of the process owning the matrix.

   .. cpp:function:: void SetRoot( int root )

      Sets the rank of the process owning the matrix (and clears the current
      contents).

Special cases used in Elemental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This list of special cases is here to help clarify the notation used throughout
Elemental's source (as well as this documentation). These are all special
cases of :cpp:type:`DistMatrix\<T,STAR,STAR>`.

.. cpp:type:: class DistMatrix<double,STAR,STAR>

   The underlying datatype is the set of double-precision real numbers.

.. cpp:type:: class DistMatrix<Complex<double>,STAR,STAR>

   The underlying datatype is the set of double-precision complex numbers.

.. cpp:type:: class DistMatrix<R,STAR,STAR>

   The underlying datatype `R` is real.

.. cpp:type:: class DistMatrix<Complex<R>,STAR,STAR>

   The underlying datatype :cpp:type:`Complex\<R>` is complex with base type 
   `R`.

.. cpp:type:: class DistMatrix<F,STAR,STAR>

   The underlying datatype `F` is a field.

