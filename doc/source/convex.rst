Convex optimization
*******************

LogBarrier
----------
Uses a careful calculation of the log of the determinant in order to return
the *log barrier* of a Hermitian positive-definite matrix `A`,
:math:`-\log(\mbox{det}(A))`.

.. cpp:function:: typename Base<F>::type LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
.. cpp:function:: typename Base<F>::type LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A )
.. cpp:function:: typename Base<F>::type LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
.. cpp:function:: typename Base<F>::type LogBarrier( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )

LogDetDivergence
----------------
The *log-det divergence* of a pair of :math:`n \times n` Hermitian
positive-definite matrices :math:`A` and :math:`B` is

.. math::

   D_{ld}(A,B) = \mbox{tr}(A B^{-1}) -\log(\mbox{det}(A B^{-1})) - n,

which can be greatly simplified using the Cholesky factors of :math:`A` and :math:`B`.
In particular, if we set :math:`Z = L_B^{-1} L_A`, where :math:`A=L_A L_A^H` and 
:math:`B=L_B L_B^H` are Cholesky factorizations, then

.. math::

   D_{ld}(A,B) = \| Z \|_F^2 - 2 \log(\mbox{det}(Z)) - n.

.. cpp:function:: typename Base<F>::type LogDetDivergence( UpperOrLower uplo, const Matrix<F>& A, const Matrix<F>& B )
.. cpp:function:: typename Base<F>::type LogDetDivergence( UpperOrLower uplo, const DistMatrix<F>& A, const DistMatrix<F>& B )

Singular-value soft-thresholding
--------------------------------
Overwrites :math:`A` with :math:`U S_{\tau}(\Sigma) V^H`, where :math:`U \Sigma V^H` is the singular-value decomposition of :math:`A` upon input and :math:`S_{\tau}` performs soft-thresholding with parameter :math:`\tau`.
The return value is the rank of the soft-thresholded matrix.

.. cpp:function:: int SVT( Matrix<F>& A, typename Base<F>::type tau )
.. cpp:function:: int SVT( DistMatrix<F>& A, typename Base<F>::type tau )

   Uses a thresholded cross-product SVD.

.. cpp:function:: int SVT( Matrix<F>& A, typename Base<F>::type tau, int numSteps )
.. cpp:function:: int SVT( DistMatrix<F>& A, typename Base<F>::type tau, int numSteps )

   Same as above, but run the thresholded cross-product SVD on the :math:`R` 
   from the partial :math:`QR` decomposition produced from `numSteps` iterations
   of (Businger/Golub) column-pivoted QR.

Soft-thresholding
-----------------
Overwrites each entry of :math:`A` with its soft-thresholded value.

.. cpp:function:: void SoftThreshold( Matrix<F>& A, typename Base<F>::type tau )
.. cpp:function:: void SoftThreshold( DistMatrix<F>& A, typename Base<F>::type tau )
