Utilities
=========

Householder reflectors
----------------------
**TODO:** Describe major difference from LAPACK's conventions (i.e., we do not 
treat the identity matrix as a Householder transform since it requires the 
:math:`u` in :math:`H=I-2uu'` to have norm zero rather than one). 

.. cpp:function:: F Reflector( Matrix<F>& chi, Matrix<F>& x )
.. cpp:function:: F Reflector( DistMatrix<F>& chi, DistMatrix<F>& x )

.. cpp:function:: F reflector::Col( DistMatrix<F>& chi, DistMatrix<F>& x )
.. cpp:function:: F reflector::Row( DistMatrix<F>& chi, DistMatrix<F>& x )

.. cpp:function:: void ApplyPackedReflectors( LeftOrRight side, UpperOrLower uplo, VerticalOrHorizontal dir, ForwardOrBackward order, Conjugation conjugation, Int offset, const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
.. cpp:function:: void ApplyPackedReflectors( LeftOrRight side, UpperOrLower uplo, VerticalOrHorizontal dir, ForwardOrBackward order, Conjugation conjugation, Int offset, const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
.. cpp:function:: void ApplyPackedReflectors( LeftOrRight side, UpperOrLower uplo, VerticalOrHorizontal dir, ForwardOrBackward order, Conjugation conjugation, Int offset, const DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& A )

.. cpp:function:: void ExpandPackedReflectors( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation, Int offset, Matrix<F>& H, const Matrix<F>& t )
.. cpp:function:: void ExpandPackedReflectors( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation, Int offset, DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t )
.. cpp:function:: void ExpandPackedReflectors( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation, Int offset, DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t )

Applying pivots
---------------

.. cpp:function:: void ApplyColumnPivots( Matrix<F>& A, const Matrix<int>& p )
.. cpp:function:: void ApplyColumnPivots( DistMatrix<F>& A, const DistMatrix<int,U,V>& p )
.. cpp:function:: void ApplyInverseColumnPivots( Matrix<F>& A, const Matrix<int>& p )
.. cpp:function:: void ApplyInverseColumnPivots( DistMatrix<F>& A, const DistMatrix<int,U,V>& p )

.. cpp:function:: void ApplyRowPivots( Matrix<F>& A, const Matrix<int>& p )
.. cpp:function:: void ApplyRowPivots( DistMatrix<F>& A, const DistMatrix<int,U,V>& p )
.. cpp:function:: void ApplyInverseRowPivots( Matrix<F>& A, const Matrix<int>& p )
.. cpp:function:: void ApplyInverseRowPivots( DistMatrix<F>& A, const DistMatrix<int,U,V>& p )

Sorting
-------

.. cpp:function:: void Sort( Matrix<Real>& X, SortType sort=ASCENDING )
.. cpp:function:: void Sort( DistMatrix<Real,U,V>& X, SortType sort=ASCENDING )

.. cpp:function:: std::vector<ValueInt<Real> > TaggedSort( const Matrix<Real>& X, SortType sort=ASCENDING )
.. cpp:function:: std::vector<ValueInt<Real> > TaggedSort( const DistMatrix<Real,U,V>& X, SortType sort=ASCENDING )

.. cpp:function:: ValueInt<Real> Median( const Matrix<Real>& x )
.. cpp:function:: ValueInt<Real> Median( const DistMatrix<Real,U,V>& x )

