Factorizations
==============

Cholesky factorization
----------------------
It is well-known that Hermitian positive-definite (HPD) matrices can be 
decomposed into the form :math:`A = L L^H` or :math:`A = U^H U`, where 
:math:`L=U^H` is lower triangular, and Cholesky factorization provides such an 
:math:`L` (or :math:`U`) given an HPD :math:`A`. If :math:`A` is found to be 
numerically indefinite, then a :cpp:type:`NonHPDMatrixException` will be 
thrown.

.. cpp:function:: void Cholesky( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void Cholesky( UpperOrLower uplo, DistMatrix<F>& A )

   Overwrite the `uplo` triangle of the HPD matrix `A` with its Cholesky factor.

It is possible to compute the Cholesky factor of a Hermitian positive
semi-definite (HPSD) matrix through its eigenvalue decomposition, though it
is significantly more expensive than the HPD case: Let :math:`A = U \Lambda U^H`
be the eigenvalue decomposition of :math:`A`, where all entries of
:math:`\Lambda` are non-negative. Then :math:`B = U \sqrt \Lambda U^H` is the
matrix square root of :math:`A`, i.e., :math:`B B = A`, and it follows that the
QR and LQ factorizations of :math:`B` yield Cholesky factors of :math:`A`:

.. math::
   A = B B = B^H B = (Q R)^H (Q R) = R^H Q^H Q R = R^H R,

and

.. math::
   A = B B = B B^H = (L Q) (L Q)^H = L Q Q^H L^H = L L^H.

If :math:`A` is found to have eigenvalues less than
:math:`-n \epsilon \| A \|_2`, then a :cpp:type:`NonHPSDMatrixException` will
be thrown.

.. cpp:function:: void HPSDCholesky( UpperOrLower uplo, Matrix<F>& A )
.. cpp:function:: void HPSDCholesky( UpperOrLower uplo, DistMatrix<F>& A )

   Overwrite the `uplo` triangle of the potentially singular matrix `A` with
   its Cholesky factor.

:math:`LDL^H` factorization
---------------------------
Though the Cholesky factorization is ideal for most HPD matrices, there exist 
many Hermitian matrices whose eigenvalues are not all positive. The 
:math:`LDL^H` factorization exists as slight relaxation of the Cholesky 
factorization, i.e., it computes lower-triangular (with unit diagonal) :math:`L`
and diagonal :math:`D` such that :math:`A = L D L^H`. If :math:`A` is found to 
be numerically singular, then a :cpp:type:`SingularMatrixException` will be 
thrown.

   .. warning::

      The following routines do not pivot, so please use with caution.

.. cpp:function:: void LDLH( Matrix<F>& A )
.. cpp:function:: void LDLH( DistMatrix<F>& A )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`.

.. cpp:function:: void LDLH( Matrix<F>& A, Matrix<F>& d )
.. cpp:function:: void LDLH( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but also return the diagonal in the column vector `d`.

:math:`LDL^T` factorization
---------------------------
While the :math:`LDL^H` factorization targets Hermitian matrices, the 
:math:`LDL^T` factorization targets symmetric matrices. If :math:`A` is found to 
be numerically singular, then a :cpp:type:`SingularMatrixException` will be 
thrown.

   .. warning::

      The following routines do not pivot, so please use with caution.

.. cpp:function:: void LDLT( Matrix<F>& A )
.. cpp:function:: void LDLT( DistMatrix<F>& A )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`.

.. cpp:function:: void LDLT( Matrix<F>& A, Matrix<F>& d )
.. cpp:function:: void LDLT( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but also return the diagonal in the vector `d`.

:math:`LU` factorization
------------------------
Given :math:`A \in \mathbb{F}^{m \times n}`, an LU factorization 
(without pivoting) finds a unit lower-trapezoidal 
:math:`L \in \mathbb{F}^{m \times \mbox{min}(m,n)}` and upper-trapezoidal 
:math:`U \in \mathbb{F}^{\mbox{min}(m,n) \times n}` such that :math:`A=LU`. 
Since :math:`L` is required to have its diaganal entries set to one: the upper 
portion of :math:`A` can be overwritten with `U`, and the strictly lower 
portion of :math:`A` can be overwritten with the strictly lower portion of 
:math:`L`. If :math:`A` is found to be numerically singular, then a 
:cpp:type:`SingularMatrixException` will be thrown.

.. cpp:function:: void LU( Matrix<F>& A )
.. cpp:function:: void LU( DistMatrix<F>& A )

   Overwrites :math:`A` with its LU decomposition.

Since LU factorization without pivoting is known to be unstable for general 
matrices, it is standard practice to pivot the rows of :math:`A` during the 
factorization (this is called partial pivoting since the columns are not also 
pivoted). An LU factorization with partial pivoting therefore computes 
:math:`P`, :math:`L`, and :math:`U` such that :math:`PA=LU`, where :math:`L` 
and :math:`U` are as described above and :math:`P` is a permutation matrix.

.. cpp:function:: void LU( Matrix<F>& A, Matrix<int>& p )
.. cpp:function:: void LU( DistMatrix<F>& A, DistMatrix<F,VC,STAR>& p )

   Overwrites the matrix :math:`A` with the LU decomposition of 
   :math:`PA`, where :math:`P` is represented by the pivot vector `p`.

.. cpp:function:: void LU( Matrix<F>& A, Matrix<int>& p, Matrix<int>& q )
.. cpp:function:: void LU( DistMatrix<F>& A, DistMatrix<F,VC,STAR>& p, DistMatrix<F,VC,STAR>& q )

   Overwrites the matrix :math:`A` with the LU decomposition of 
   :math:`PAQ`, where :math:`P` is represented by the pivot vector `p`, 
   and likewise for :math:`Q`.

:math:`LQ` factorization
------------------------
Given :math:`A \in \mathbb{F}^{m \times n}`, an LQ factorization typically 
computes an implicit unitary matrix :math:`\hat Q \in \mathbb{F}^{n \times n}` 
such that :math:`\hat L \equiv A\hat Q^H` is lower trapezoidal. One can then 
form the thin factors :math:`L \in \mathbb{F}^{m \times \mbox{min}(m,n)}` and 
:math:`Q \in \mathbb{F}^{\mbox{min}(m,n) \times n}` by setting 
:math:`L` and :math:`Q` to first :math:`\mbox{min}(m,n)` columns and rows of 
:math:`\hat L` and :math:`\hat Q`, respectively. Upon completion :math:`L` is 
stored in the lower trapezoid of :math:`A` and the Householder reflectors 
representing :math:`\hat Q` are stored within the rows of the strictly upper 
trapezoid.

.. cpp:function:: void LQ( Matrix<F>& A )
.. cpp:function:: void LQ( DistMatrix<F>& A )
.. cpp:function:: void LQ( Matrix<F>& A, Matrix<F>& t )
.. cpp:function:: void LQ( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )

   Overwrite the complex matrix :math:`A` with :math:`L` and the 
   Householder reflectors representing :math:`\hat Q`. In the complex case, 
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector `t`.

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:function:: void lq::ApplyQ( LeftOrRight side, Orientation orientation, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
.. cpp:function:: void lq::ApplyQ( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
.. cpp:function:: void lq::ApplyQ( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )

   Applies the implicitly-defined :math:`Q` (or its adjoint) stored within
   `A` and `t` from either the left or the right to :math:`B`.

:math:`QR` factorization
------------------------
Given :math:`A \in \mathbb{F}^{m \times n}`, a QR factorization typically 
computes an implicit unitary matrix :math:`\hat Q \in \mathbb{F}^{m \times m}` 
such that :math:`\hat R \equiv \hat Q^H A` is upper trapezoidal. One can then 
form the thin factors :math:`Q \in \mathbb{F}^{m \times \mbox{min}(m,n)}` and
:math:`R \in \mathbb{F}^{\mbox{min}(m,n) \times n}` by setting 
:math:`Q` and :math:`R` to first :math:`\mbox{min}(m,n)` columns and rows of 
:math:`\hat Q` and :math:`\hat R`, respectively. Upon completion :math:`R` is 
stored in the upper trapezoid of :math:`A` and the Householder reflectors 
representing :math:`\hat Q` are stored within the columns of the strictly lower 
trapezoid.

.. cpp:function:: void QR( Matrix<F>& A )
.. cpp:function:: void QR( DistMatrix<F>& A )
.. cpp:function:: void QR( Matrix<F>& A, Matrix<F>& t )
.. cpp:function:: void QR( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )

   Overwrite the complex matrix :math:`A` with :math:`R` and the 
   Householder reflectors representing :math:`\hat Q`. In the complex case,
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector `t`.

.. cpp:function:: void QR( Matrix<F>& A, Matrix<int>& p )
.. cpp:function:: void QR( DistMatrix<F>& A, DistMatrix<int,VR,STAR>& p )
.. cpp:function:: void QR( Matrix<F>& A, Matrix<F>& t, Matrix<int>& p )
.. cpp:function:: void QR( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<int,VR,STAR>& p )

   Column-pivoted QR factorization. The current implementation uses 
   Businger-Golub pivoting.

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:function:: void qr::Explicit( Matrix<F>& A, bool colPiv=false )
.. cpp:function:: void qr::Explicit( DistMatrix<F>& A, bool colPiv=false )

   Overwrite :math:`A` with the orthogonal matrix from its QR factorization
   (with or without column pivoting).

.. cpp:function:: void qr::Explicit( Matrix<F>& A, Matrix<F>& R, bool colPiv=false )
.. cpp:function:: void qr::Explicit( DistMatrix<F>& A, DistMatrix<F>& R, bool colPiv=false )

   Additionally explicitly return the :math:`R` from the QR factorization.

.. cpp:function:: void qr::ApplyQ( LeftOrRight side, Orientation orientation, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
.. cpp:function:: void qr::ApplyQ( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
.. cpp:function:: void qr::ApplyQ( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )

   Applies the implicitly-defined :math:`Q` (or its adjoint) stored within
   `A` and `t` from either the left or the right to :math:`B`.

.. cpp:function:: void qr::BusingerGolub( Matrix<F>& A, Matrix<int>& p )
.. cpp:function:: void qr::BusingerGolub( DistMatrix<F>& A, DistMatrix<int,VR,STAR>& p )
.. cpp:function:: void qr::BusingerGolub( Matrix<F>& A, Matrix<F>& t, Matrix<int>& p )
.. cpp:function:: void qr::BusingerGolub( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<int,VR,STAR>& p )

   Column-pivoted versions of the above routines which use the Businger/Golub 
   strategy, i.e., the pivot is chosen as the remaining column with maximum
   two norm.

.. cpp:function:: void qr::BusingerGolub( Matrix<F>& A, Matrix<int>& p, int numSteps )
.. cpp:function:: void qr::BusingerGolub( DistMatrix<F>& A, DistMatrix<int,VR,STAR>& p, int numSteps )
.. cpp:function:: void qr::BusingerGolub( Matrix<F>& A, Matrix<F>& t, Matrix<int>& p, int numSteps )
.. cpp:function:: void qr::BusingerGolub( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<int,VR,STAR>& p, int numSteps )

   Same as above, but only execute a fixed number of steps of the rank-revealing
   factorization.

.. cpp:function:: void qr::BusingerGolub( Matrix<F>& A, Matrix<int>& p, int maxSteps, R tol )
.. cpp:function:: void qr::BusingerGolub( DistMatrix<F>& A, DistMatrix<int,VR,STAR>& p, int maxSteps, R tol )
.. cpp:function:: void qr::BusingerGolub( Matrix<F>& A, Matrix<F>& t, Matrix<int>& p, int maxSteps, R tol )
.. cpp:function:: void qr::BusingerGolub( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<int,VR,STAR>& p, int maxSteps, R tol )

   Either execute `maxSteps` iterations or stop after the maximum remaining 
   column norm is less than or equal to `tol` times the maximum original column
   norm.

:math:`RQ` factorization
------------------------
Just like an LQ factorization, but the orthogonalization process starts from the bottom row and produces a 
much sparser triangular factor when the matrix is wider than it is tall.

.. cpp:function:: void RQ( Matrix<F>& A )
.. cpp:function:: void RQ( DistMatrix<F>& A )
.. cpp:function:: void RQ( Matrix<F>& A, Matrix<F>& t )
.. cpp:function:: void RQ( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )

   Overwrite the complex matrix :math:`A` with :math:`R` and the 
   Householder reflectors representing :math:`\hat Q`. In the complex case, 
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector `t`.

Detailed interface
^^^^^^^^^^^^^^^^^^

.. cpp:function:: void rq::ApplyQ( LeftOrRight side, Orientation orientation, const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
.. cpp:function:: void rq::ApplyQ( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
.. cpp:function:: void rq::ApplyQ( LeftOrRight side, Orientation orientation, const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )

   Applies the implicitly-defined :math:`Q` (or its adjoint) stored within
   `A` and `t` from either the left or the right to :math:`B`.

Interpolative Decomposition (ID)
--------------------------------
Interpolative Decompositions (ID's) are closely related to pivoted QR 
factorizations and are useful for representing (approximately) low-rank 
matrices in terms of linear combinations of a few of their columns, i.e., 

.. math::

   A P = \hat{A} \begin{pmatrix} I & Z \end{pmatrix},

where :math:`P` is a permutation matrix, :math:`\hat{A}` is a small set of 
columns of :math:`A`, and :math:`Z` is an interpolation matrix responsible for
representing the remaining columns in terms of the selected columns of 
:math:`A`.

.. cpp:function:: void ID( const Matrix<F>& A, Matrix<int>& p, Matrix<F>& Z, int numSteps )
.. cpp:function:: void ID( const DistMatrix<F>& A, DistMatrix<int,VR,STAR>& p, DistMatrix<F,STAR,VR>& Z, int numSteps )

   `numSteps` steps of a pivoted QR factorization are used to return an 
   Interpolative Decomposition of :math:`A`.

.. cpp:function:: void ID( const Matrix<F>& A, Matrix<int>& p, Matrix<F>& Z, int maxSteps, typename Base<F>::type tol )
.. cpp:function:: void ID( const DistMatrix<F>& A, DistMatrix<int,VR,STAR>& p, DistMatrix<F,STAR,VR>& Z, int maxSteps, typename Base<F>::type tol )

   Either `maxSteps` steps of a pivoted QR factorization are used, or 
   executation stopped after the maximum remaining column norm was less than or
   equal to `tol` times the maximum original column norm.

Skeleton decomposition
----------------------
Skeleton decompositions are essentially two-sided interpolative decompositions,
but the terminology is unfortunately extremely contested. We follow the 
convention that a skeleton decomposition is an approximation

.. math::

   A \approx A_C Z A_R,

where :math:`A_C` is a (small) selection of columns of :math:`A`, 
:math:`A_R` is a (small) selection of rows of :math:`A`, and :math:`Z` is a 
(small) square matrix. When :math:`Z` is allowed to be rectangular, it is more
common to call this a CUR decomposition.

.. cpp:function:: void Skeleton( const Matrix<F>& A, Matrix<int>& pR, Matrix<int>& pC, Matrix<F>& Z, int maxSteps, typename Base<F>::type tol )
.. cpp:function:: void Skeleton( const DistMatrix<F>& A, DistMatrix<int,VR,STAR>& pR, DistMatrix<int,VR,STAR>& pC, int maxSteps, typename Base<F>::type tol )

   Rather than returning :math:`A_R` and :math:`A_C`, the permutation matrices
   which implicitly define them are returned instead. At most `maxSteps` steps 
   of a pivoted QR decomposition will be used in order to generate the 
   row/column subsets, and less steps will be taken if a pivot norm is less 
   than or equal to `tolerance` times the first pivot norm.
