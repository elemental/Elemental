Reduction to condensed form
===========================

Hermitian to tridiagonal
------------------------
The currently best-known algorithms for computing eigenpairs of dense Hermitian 
matrices begin by performing a unitary similarity transformation which reduces 
the matrix to real symmetric tridiagonal form (usually through Householder 
transformations). This routine performs said reduction on a Hermitian matrix 
and stores the scaled Householder vectors in place of the introduced zeroes. 

.. cpp:function:: void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void HermitianTridiag( UpperOrLower uplo, DistMatrix<F>& A )
.. cpp:function:: void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t )
.. cpp:function:: void HermitianTridiag( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )

   Overwrites the main and sub (super) diagonal of the real matrix 
   `A` with its similar symmetric tridiagonal matrix and stores the scaled 
   Householder vectors below (above) its tridiagonal entries.
   Complex Hermitian reductions have the added complication of needing to 
   also store the phase information for the Householder vectors (the scaling can
   be inferred since the Householder vectors must be unit length) if they are 
   to be applied (in the column vector `t`). 

Please see the :ref:`lapack-tuning` section for extensive information on 
maximizing the performance of Householder tridiagonalization.

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:function:: void hermitian_tridiag::ApplyQ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
.. cpp:function:: void hermitian_tridiag::ApplyQ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
.. cpp:function:: void hermitian_tridiag::ApplyQ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )

   Apply (from the left or right) the implicitly defined unitary matrix 
   (or its adjoint) represented by the Householder transformations stored within
   the specified triangle of `A` and the phase information stored in `t`.

General to Hessenberg
---------------------
Not yet written, but it is planned and relatively straightforward after 
writing the reductions to tridiagonal and bidiagonal form.

General to bidiagonal
---------------------
Reduces a general fully-populated :math:`m \times n` matrix to bidiagonal form 
through two-sided Householder transformations; when the :math:`m \ge n`, the 
result is upper bidiagonal, otherwise it is lower bidiagonal. This routine is 
most commonly used as a preprocessing step in computing the SVD of a general
matrix.

.. cpp:function:: void Bidiag( Matrix<F>& A )
.. cpp:function:: void Bidiag( DistMatrix<F>& A )
.. cpp:function:: void Bidiag( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
.. cpp:function:: void Bidiag( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )

   Overwrites the main and sub (or super) diagonal of the real matrix `A` with 
   the resulting bidiagonal matrix and stores the scaled Householder vectors in 
   the remainder of the matrix.
   The complex case must also store the phases of the Householder 
   transformations (in `tP` and `tQ`) if they are to be applied.

   .. note:: The :math:`m < n` case is not yet supported for the distributed 
             version.

Detailed interface
^^^^^^^^^^^^^^^^^^
**TODO**

.. cpp:function:: void ApplyU( LeftOrRight side, Orientation orientation, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
.. cpp:function:: void ApplyU( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
.. cpp:function:: void ApplyU( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )

.. cpp:function:: void ApplyV( LeftOrRight side, Orientation orientation, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
.. cpp:function:: void ApplyV( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
.. cpp:function:: void ApplyV( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )

