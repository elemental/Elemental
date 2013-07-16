Control theory
**************

The following algorithms draw heavily from the second chapter of 
Nicholas J. Higham's "Functions of Matrices: Theory and Computation" and depend
heavily on the matrix sign function. They have only undergone cursory testing.

Sylvester
---------
As long as both :math:`A` and :math:`B` only have eigenvalues in the open 
right-half plane, the following routines solve for :math:`X` in the *Sylvester
equation*

.. math::

   A X + X B = C

via computing :math:`\text{sgn}(W)`, where

.. math::

   W = \begin{pmatrix} A & -C \\ 0 & -B \end{pmatrix}.

.. cpp:function:: int Sylvester( const Matrix<F>& A, const Matrix<F>& B, const Matrix<F>& C, Matrix<F>& X )
.. cpp:function:: int Sylvester( const DistMatrix<F>& A, const DistMatrix<F>& B, const DistMatrix<F>& C, DistMatrix<F>& X )

One may also directly pass in :math:`W` in order to save space.

.. cpp:function:: int Sylvester( int m, Matrix<F>& W, Matrix<F>& X )
.. cpp:function:: int Sylvester( int m, DistMatrix<F>& W, DistMatrix<F>& X )

Lyapunov
--------
A special case of the Sylvester solver, where :math:`B = A^H`.

.. cpp:function:: int Lyapunov( const Matrix<F>& A, const Matrix<F>& C, Matrix<F>& X )
.. cpp:function:: int Lyapunov( const DistMatrix<F>& A, const DistMatrix<F>& C, DistMatrix<F>& X )


Algebraic Ricatti
-----------------
Under suitable conditions, the following routines solve for :math:`X` in the 
*algebraic Ricatti equation*

.. math::

   X K X - A^H X - X A = L,

where both :math:`K` and :math:`L` are Hermitian. In each case, the number
of Newton iterations required for convergence of :math:`\text{sgn}(W)`, where

.. math::

   W = \begin{pmatrix} A^H & L \\ K & -A \end{pmatrix},

is returned.

.. cpp:function:: int Ricatti( UpperOrLower uplo, const Matrix<F>& A, const Matrix<F>& K, const Matrix<F>& L, Matrix<F>& X )
.. cpp:function:: int Ricatti( UpperOrLower uplo, const DistMatrix<F>& A, const DistMatrix<F>& K, const DistMatrix<F>& L, DistMatrix<F>& X )

Alternatively, one may directly fill the matrix :math:`W`.

.. cpp:function:: int Ricatti( Matrix<F>& W, Matrix<F>& X )
.. cpp:function:: int Ricatti( DistMatrix<F>& W, DistMatrix<F>& X )
