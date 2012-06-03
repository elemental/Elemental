Factorizations
==============

Cholesky factorization
----------------------
It is well-known that Hermitian positive-definite (HPD) matrices can be decomposed
into the form :math:`A = L L^H` or :math:`A = U^H U`, where :math:`L=U^H` is lower
triangular, and Cholesky factorization provides such an :math:`L` (or :math:`U`) 
given an HPD :math:`A`. If :math:`A` is found to be numerically indefinite, then 
a ``NonHPDMatrixException`` will be thrown.

.. cpp:function:: void Cholesky( UpperOrLower uplo, Matrix<F>& A )

   Overwrite the `uplo` triangle of the HPD matrix `A` with its Cholesky factor.

.. cpp:function:: void Cholesky( UpperOrLower uplo, DistMatrix<F>& A )

   Overwrite the `uplo` triangle of the distributed HPD matrix `A` with its 
   Cholesky factor.

It is also possible to form the Cholesky factors of Hermitian positive 
semi-definite (HPSD) matrices through their eigenvalue decomposition, though it 
is significantly more expensive than the HPD case: Let :math:`A = U \Lambda U^H`
be the eigenvalue decomposition of :math:`A`, where all entries of :math:`\Lambda`
are non-negative. Then :math:`B = U \sqrt \Lambda U^H` is the matrix square root
of :math:`A`, i.e., :math:`B B = A`, and it follows that the QR and LQ 
factorizations of :math:`B` yield Cholesky factors of :math:`A`:

.. math::
   A = B B = B^H B = (Q R)^H (Q R) = R^H Q^H Q R = R^H R,

and

.. math::
   A = B B = B B^H = (L Q) (L Q)^H = L Q Q^H L^H = L L^H.

If :math:`A` is found to have eigenvalues less than :math:`-n \epsilon ||A||_2`, 
then a ``NonHPSDMatrixException`` will be thrown.

.. cpp:function:: void HPSDCholesky( UpperOrLower uplo, DistMatrix<F>& A )

   Overwrite the `uplo` triangle of the distributed HPSD matrix `A` with its
   Cholesky factor.

:math:`LDL^H` factorization
---------------------------
Though the Cholesky factorization is ideal for most HPD matrices, there exist 
many Hermitian matrices whose eigenvalues are not all positive. The 
:math:`LDL^H` factorization exists as slight relaxation of the Cholesky 
factorization, i.e., it computes lower-triangular (with unit diagonal) :math:`L`
and diagonal :math:`D` such that :math:`A = L D L^H`. If :math:`A` is found to 
be numerically singular, then a ``SingularMatrixException`` will be thrown.

   .. warning::

      The following routines do not pivot, so please use with caution.

.. cpp:function:: void LDLH( Matrix<F>& A )
.. cpp:function:: void LDLH( DistMatrix<F>& A )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`.

.. cpp:function:: void LDLH( Matrix<F>& A, Matrix<F>& d )
.. cpp:function:: void LDLH( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but also return the diagonal in the column vector ``d``.

:math:`LDL^T` factorization
---------------------------
While the :math:`LDL^H` factorization targets Hermitian matrices, the 
:math:`LDL^T` factorization targets symmetric matrices. If :math:`A` is found to 
be numerically singular, then a ``SingularMatrixException`` will be thrown.

   .. warning::

      The following routines do not pivot, so please use with caution.

.. cpp:function:: void LDLT( Matrix<F>& A )
.. cpp:function:: void LDLT( DistMatrix<F>& A )

   Overwrite the strictly lower triangle of :math:`A` with the strictly lower 
   portion of :math:`L` (:math:`L` implicitly has ones on its diagonal) and 
   the diagonal with :math:`D`.

.. cpp:function:: void LDLT( Matrix<F>& A, Matrix<F>& d )
.. cpp:function:: void LDLT( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )

   Same as above, but also return the diagonal in the vector ``d``.

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
``SingularMatrixException`` will be thrown.

.. cpp:function:: void LU( Matrix<F>& A )

   Overwrites :math:`A` with its LU decomposition.

.. cpp:function:: void LU( DistMatrix<F>& A )

   Overwrites :math:`A` with its LU decomposition.

Since LU factorization without pivoting is known to be unstable for general 
matrices, it is standard practice to pivot the rows of :math:`A` during the 
factorization (this is called partial pivoting since the columns are not also 
pivoted). An LU factorization with partial pivoting therefore computes 
:math:`P`, :math:`L`, and :math:`U` such that :math:`PA=LU`, where :math:`L` 
and :math:`U` are as described above and :math:`P` is a permutation matrix.

.. cpp:function:: void LU( Matrix<F>& A, Matrix<int>& p )

   Ovewrites :math:`A` with the LU decomposition of :math:`PA`, where 
   :math:`P` is represented by the pivot vector `p`.

.. cpp:function:: void LU( DistMatrix<F>& A, DistMatrix<F,VC,STAR>& p )

   Overwrites the distributed matrix :math:`A` with the LU decomposition of 
   :math:`PA`, where :math:`P` is represented by the pivot vector `p`.

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

.. cpp:function:: void LQ( Matrix<R>& A )
.. cpp:function:: void LQ( DistMatrix<R>& A )

   Overwrite the real matrix :math:`A` with :math:`L` and the 
   Householder reflectors representing :math:`\hat Q`.

.. cpp:function:: void LQ( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
.. cpp:function:: void LQ( DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,MD,STAR>& t )

   Overwrite the complex matrix :math:`A` with :math:`L` and the 
   Householder reflectors representing :math:`\hat Q`; unlike the real case, 
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector `t`.

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

.. cpp:function:: void QR( Matrix<R>& A )
.. cpp:function:: void QR( DistMatrix<R>& A )

   Overwrite the real matrix :math:`A` with :math:`R` and the 
   Householder reflectors representing :math:`\hat Q`.

.. cpp:function:: void QR( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
.. cpp:function:: void QR( DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,MD,STAR>& t )

   Overwrite the complex matrix :math:`A` with :math:`R` and the 
   Householder reflectors representing :math:`\hat Q`; unlike the real case, 
   phase information is needed in order to define the (generalized) 
   Householder transformations and is stored in the column vector `t`.

