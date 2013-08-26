Viewing
=======

View a full matrix
------------------

.. cpp:function:: void View( Matrix<T>& A, Matrix<T>& B )
.. cpp:function:: void View( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B )

   Make `A` a view of the matrix `B`.

.. cpp:function:: Matrix<T> View( Matrix<T>& B )
.. cpp:function:: DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B )

   Return a view of the matrix `B`.

.. cpp:function:: void LockedView( Matrix<T>& A, const Matrix<T>& B )
.. cpp:function:: void LockedView( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B )

   Make `A` a non-mutable view of the matrix `B`.

.. cpp:function:: Matrix<T> LockedView( const Matrix<T>& B )
.. cpp:function:: DistMatrix<T,U,V> LockedView( const DistMatrix<T,U,V>& B )

   Return a view of the matrix `B`.

View a submatrix
----------------

.. cpp:function:: void View( Matrix<T>& A, Matrix<T>& B, int i, int j, int height, int width )
.. cpp:function:: void View( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B, int i, int j, int height, int width )

   Make `A` a view of the `height x width` submatrix of `B` starting at 
   coordinate ``(i,j)``.

.. cpp:function:: Matrix<T> View( Matrix<T>& B, int i, int j, int height, int width )
.. cpp:function:: DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B, int i, int j, int height, int width )

   Return a view of the specified submatrix of `B`.

.. cpp:function:: void LockedView( Matrix<T>& A, const Matrix<T>& B, int i, int j, int height, int width )
.. cpp:function:: void LockedView( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B, int i, int j, int height, int width )

   Make `A` a non-mutable view of the `height x width` submatrix of `B` 
   starting at coordinate ``(i,j)``.

.. cpp:function:: Matrix<T> LockedView( const Matrix<T>& B, int i, int j, int height, int width )
.. cpp:function:: DistMatrix<T,U,V> LockedView( const DistMatrix<T,U,V>& B, int i, int j, int height, int width )

   Return an immutable view of the specified submatrix of `B`.

View 1x2 matrices
-----------------

.. cpp:function:: void View1x2( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR )
.. cpp:function:: void View1x2( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR )

   Make `A` a view of the matrix 
   :math:`\left(\begin{array}{cc} B_L & B_R \end{array}\right)`.

.. cpp:function:: Matrix<T> View1x2( Matrix<T>& BL, Matrix<T>& BR )
.. cpp:function:: DistMatrix<T,U,V> View1x2( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR )

   Return a view of the merged matrix.

.. cpp:function:: void LockedView1x2( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR )
.. cpp:function:: void LockedView1x2( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR )

   Make `A` a non-mutable view of the matrix 
   :math:`\left(\begin{array}{cc} B_L & B_R \end{array}\right)`.

.. cpp:function:: Matrix<T> LockedView1x2( const Matrix<T>& BL, const Matrix<T>& BR )
.. cpp:function:: DistMatrix<T,U,V> LockedView1x2( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR )

   Return an immutable view of the merged matrix.

View 2x1 matrices
-----------------

.. cpp:function:: void View2x1( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB )
.. cpp:function:: void View2x1( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB )

   Make `A` a view of the matrix 
   :math:`\left(\begin{array}{c} B_T \\ B_B \end{array}\right)`.

.. cpp:function:: Matrix<T> View2x1( Matrix<T>& BT, Matrix<T>& BB )
.. cpp:function:: DistMatrix<T,U,V> View2x1( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB )

   Return a view of the merged matrix.

.. cpp:function:: void LockedView2x1( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB )
.. cpp:function:: void LockedView2x1( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB )

   Make `A` a non-mutable view of the matrix 
   :math:`\left(\begin{array}{c} B_T \\ B_B \end{array}\right)`.

.. cpp:function:: Matrix<T> LockedView2x1( const Matrix<T>& BT, const Matrix<T>& BB )
.. cpp:function:: DistMatrix<T,U,V> LockedView2x1( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB )

   Return a view of the merged matrix.

View 2x2 matrices
-----------------

.. cpp:function:: void View2x2( Matrix<T>& A, Matrix<T>& BTL, Matrix<T>& BTR, Matrix<T>& BBL, Matrix<T>& BBR )
.. cpp:function:: void View2x2( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR, DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR )

   Make `A` a view of the matrix 
   :math:`\left(\begin{array}{cc} B_{TL} & B_{TR} \\ B_{BB} & B_{BR} \end{array}\right)`.

.. cpp:function:: Matrix<T> View2x2( Matrix<T>& BTL, Matrix<T>& BTR, Matrix<T>& BBL, Matrix<T>& BBR )
.. cpp:function:: DistMatrix<T,U,V> View2x2( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR, DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR )

   Return a view of the merged matrix.

.. cpp:function:: void LockedView2x2( Matrix<T>& A, const Matrix<T>& BTL, const Matrix<T>& BTR, const Matrix<T>& BBL, const Matrix<T>& BBR )
.. cpp:function:: void LockedView2x2( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR, const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR )

    Make `A` a non-mutable view of the matrix 
    :math:`\left(\begin{array}{cc} B_{TL} & B_{TR} \\ B_{BB} & B_{BR} \end{array}\right)`.

.. cpp:function:: Matrix<T> LockedView2x2( const Matrix<T>& BTL, const Matrix<T>& BTR, const Matrix<T>& BBL, const Matrix<T>& BBR )
.. cpp:function:: DistMatrix<T,U,V> LockedView2x2( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR, const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR )

   Return an immutable view of the merged matrix.
