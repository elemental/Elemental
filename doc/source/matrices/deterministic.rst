Deterministic
=============

Cauchy
------
An :math:`m \times n` matrix :math:`A` is called *Cauchy* if there exist 
vectors :math:`x` and :math:`y` such that 

.. math::

   \alpha_{i,j} = \frac{1}{\chi_i - \eta_j},

where :math:`\chi_i` is the :math:`i`'th entry of :math:`x` and :math:`\eta_j`
is the :math:`j`'th entry of :math:`y`.

.. cpp:function:: void Cauchy( const std::vector<F>& x, const std::vector<F>& y, Matrix<F>& A )
.. cpp:function:: void Cauchy( const std::vector<F>& x, const std::vector<F>& y, DistMatrix<F,U,V>& A )

   Generate a Cauchy matrix using the defining vectors, :math:`x` and :math:`y`. 

Cauchy-like
-----------
An :math:`m \times n` matrix :math:`A` is called *Cauchy-like* if there exist 
vectors :math:`r`, :math:`s`, :math:`x`, and :math:`y` such that 

.. math::

   \alpha_{i,j} = \frac{\rho_i \psi_j}{\chi_i - \eta_j},

where :math:`\rho_i` is the :math:`i`'th entry of :math:`r`, :math:`\psi_j` is the :math:`j`'th 
entry of :math:`s`, :math:`\chi_i` is the :math:`i`'th entry of :math:`x`, and :math:`\eta_j`
is the :math:`j`'th entry of :math:`y`.

.. cpp:function:: void CauchyLike( const std::vector<F>& r, const std::vector<F>& s, const std::vector<F>& x, const std::vector<F>& y, Matrix<F>& A )
.. cpp:function:: void CauchyLike( const std::vector<F>& r, const std::vector<F>& s, const std::vector<F>& x, const std::vector<F>& y, DistMatrix<F,U,V>& A )

   Generate a Cauchy-like matrix using the defining vectors: :math:`r`, :math:`s`, :math:`x`, and :math:`y`.

Circulant
---------
An :math:`n \times n` matrix :math:`A` is called *circulant* if there exists a vector :math:`b` 
such that 

.. math::

   \alpha_{i,j} = \beta_{(i-j) \bmod n},

where :math:`\beta_k` is the :math:`k`'th entry of vector :math:`b`.

.. cpp:function:: void Circulant( const std::vector<T>& a, Matrix<T>& A )
.. cpp:function:: void Circulant( const std::vector<T>& a, DistMatrix<T,U,V>& A )

   Generate a circulant matrix using the vector ``a``.

Diagonal
--------
An :math:`n \times n` matrix :math:`A` is called *diagonal* if each entry :math:`(i,j)`, where 
:math:`i \neq j`, is :math:`0`. They are therefore defined by the *diagonal* values, where 
:math:`i = j`.

.. cpp:function:: void Diagonal( const std::vector<T>& d, Matrix<T>& D )
.. cpp:function:: void Diagonal( const std::vector<T>& d, DistMatrix<T,U,V>& D )

   Construct a diagonal matrix from the vector of diagonal values, :math:`d`.

DiscreteFourier
---------------
The :math:`n \times n` *Discrete Fourier Transform* (DFT) matrix, say :math:`A`, is given by

.. math::

   \alpha_{i,j} = \frac{e^{-2\pi i j / n}}{\sqrt{n}}.

.. cpp:function:: void DiscreteFourier( int n, Matrix<Complex<R> >& A )
.. cpp:function:: void DiscreteFourier( int n, DistMatrix<Complex<R>,U,V>& A )

   Set the matrix ``A`` equal to the :math:`n \times n` DFT matrix.

.. cpp:function:: void MakeDiscreteFourier( Matrix<Complex<R> >& A )
.. cpp:function:: void MakeDiscreteFourier( DistMatrix<Complex<R>,U,V>& A )

   Turn the existing :math:`n \times n` matrix ``A`` into a discrete Fourier 
   matrix.

Hankel
------
An :math:`m \times n` matrix :math:`A` is called a *Hankel matrix* if there 
exists a vector :math:`b` such that

.. math::

   \alpha_{i,j} = \beta_{i+j},

where :math:`\alpha_{i,j}` is the :math:`(i,j)` entry of :math:`A` and 
:math:`\beta_k` is the :math:`k`'th entry of the vector :math:`b`.

.. cpp:function:: void Hankel( int m, int n, const std::vector<T>& b, Matrix<T>& A )
.. cpp:function:: void Hankel( int m, int n, const std::vector<T>& b, DistMatrix<T,U,V>& A )

   Create an :math:`m \times n` Hankel matrix from the generate vector, 
   :math:`b`.

Hilbert
-------
The Hilbert matrix of order :math:`n` is the :math:`n \times n` matrix where
entry :math:`(i,j)` is equal to :math:`1/(i+j+1)`.

.. cpp:function:: void Hilbert( int n, Matrix<F>& A )
.. cpp:function:: void Hilbert( int n, DistMatrix<F,U,V>& A )

   Generate the :math:`n \times n` Hilbert matrix ``A``.

.. cpp:function:: void MakeHilbert( Matrix<F>& A )
.. cpp:function:: void MakeHilbert( DistMatrix<F,U,V>& A )

   Turn the square matrix ``A`` into a Hilbert matrix.

Identity
--------
The :math:`n \times n` *identity matrix* is simply defined by setting entry 
:math:`(i,j)` to one if :math:`i = j`, and zero otherwise. For various 
reasons, we generalize this definition to nonsquare, :math:`m \times n`, 
matrices.

.. cpp:function:: void Identity( int m, int n, Matrix<T>& A )
.. cpp:function:: void Identity( int m, int n, DistMatrix<T,U,V>& A )

   Set the matrix ``A`` equal to the :math:`m \times n` identity(-like) matrix.

.. cpp:function:: void MakeIdentity( Matrix<T>& A )
.. cpp:function:: void MakeIdentity( DistMatrix<T,U,V>& A ) 

   Set the matrix ``A`` to be identity-like.

Kahan
-----
For any pair :math:`(\phi,\zeta)` such that :math:`|\phi|^2+|\zeta|^2=1`, 
the corresponding :math:`n \times n` Kahan matrix is given by:

.. math::

   K = \text{diag}(1,\phi,\ldots,\phi^{n-1}) \begin{pmatrix} 
   1      & -\zeta & -\zeta & \cdots & -\zeta \\
   0      & 1      & -\zeta & \cdots & -\zeta \\
          & \ddots &        & \vdots & \vdots \\
   \vdots &        &        & 1      & -\zeta \\
   0      &        & \cdots &        & 1 \end{pmatrix}

.. cpp:function:: void Kahan( F phi, int n, Matrix<F>& A )
.. cpp:function:: void Kahan( F phi, int n, DistMatrix<F>& A )

   Sets the matrix ``A`` equal to the :math:`n \times n` Kahan matrix with 
   the specified value for :math:`\phi`.

Legendre
--------
The :math:`n \times n` tridiagonal Jacobi matrix associated with the Legendre
polynomials. Its main diagonal is zero, and the off-diagonal terms are given 
by 

.. math::

   \beta_j = \frac{1}{2}\left(1-(2(j+1))^{-2}\right)^{-1/2},

where :math:`\beta_j` connects the :math:`j`'th degree of freedom to the 
:math:`j+1`'th degree of freedom, counting from zero.
The eigenvalues of this matrix lie in :math:`[-1,1]` and are the locations for 
Gaussian quadrature of order :math:`n`. The corresponding weights may be found 
by doubling the square of the first entry of the corresponding normalized 
eigenvector.

.. cpp:function:: void Legendre( int n, Matrix<F>& A )
.. cpp:function:: void Legendre( int n, DistMatrix<F,U,V>& A )

   Sets the matrix ``A`` equal to the :math:`n \times n` Jacobi matrix.

Ones
----
Create an :math:`m \times n` matrix of all ones.

.. cpp:function:: void Ones( int m, int n, Matrix<T>& A )
.. cpp:function:: void Ones( int m, int n, DistMatrix<T,U,V>& A )

   Set the matrix ``A`` to be an :math:`m \times n` matrix of all ones.

Change all entries of the matrix :math:`A` to one.

.. cpp:function:: void MakeOnes( Matrix<T>& A )
.. cpp:function:: void MakeOnes( DistMatrix<T,U,V>& A )

   Change the entries of the matrix to ones.

OneTwoOne
---------
A "1-2-1" matrix is tridiagonal with a diagonal of all twos and sub- and 
super-diagonals of all ones.

.. cpp:function:: void OneTwoOne( int n, Matrix<T>& A )
.. cpp:function:: void OneTwoOne( int n, DistMatrix<T,U,V>& A )

   Set ``A`` to a :math:`n \times n` "1-2-1" matrix.

.. cpp:function:: void MakeOneTwoOne( Matrix<T>& A )
.. cpp:function:: void MakeOneTwoOne( DistMatrix<T,U,V>& A )

   Modify the entries of the square matrix ``A`` to be "1-2-1".

Toeplitz
--------
An :math:`m \times n` matrix is *Toeplitz* if there exists a vector :math:`b` such that, for each entry :math:`\alpha_{i,j}` of :math:`A`,

.. math::

   \alpha_{i,j} = \beta_{i-j+(n-1)},

where :math:`\beta_k` is the :math:`k`'th entry of :math:`b`.

.. cpp:function:: void Toeplitz( int m, int n, const std::vector<T>& b, Matrix<T>& A )
.. cpp:function:: void Toeplitz( int m, int n, const std::vector<T>& b, DistMatrix<T,U,V>& A )

   Build the matrix ``A`` using the generating vector :math:`b`.

.. cpp:function:: void MakeToeplitz( const std::vector<T>& b, Matrix<T>& A )
.. cpp:function:: void MakeToeplitz( const std::vector<T>& b, DistMatrix<T,U,V>& A )

   Turn the matrix ``A`` into a Toeplitz matrix defined from the generating 
   vector :math:`b`.

Walsh
-----
The Walsh matrix of order :math:`k` is a :math:`2^k \times 2^k` matrix, where

.. math::

   W_1 = \left(\begin{array}{cc} 1 & 1 \\ 1 & -1 \end{array}\right),

and 

.. math::

   W_k = \left(\begin{array}{cc} W_{k-1} & W_{k-1} \\ W_{k-1} & -W_{k-1} 
               \end{array}\right).

A *binary* Walsh matrix changes the bottom-right entry of :math:`W_1` from 
:math:`-1` to :math:`0`.

.. cpp:function:: void Walsh( int k, Matrix<T>& W, bool binary=false )
.. cpp:function:: void Walsh( int k, DistMatrix<T,U,V>& W, bool binary=false )

   Set the matrix :math:`W` equal to the :math:`k`'th (possibly binary) Walsh 
   matrix.

Wilkinson
---------
A *Wilkinson matrix* of order :math:`k` is a tridiagonal matrix with diagonal

.. math::

   [k,k-1,k-2,...,1,0,1,...,k-2,k-1,k],

and sub- and super-diagonals of all ones.

.. cpp:function:: void Wilkinson( int k, Matrix<T>& W )
.. cpp:function:: void Wilkinson( int k, DistMatrix<T,U,V>& W )

   Set the matrix :math:`W` equal to the :math:`k`'th Wilkinson matrix.

Zeros
-----
Create an :math:`m \times n` matrix of all zeros.

.. cpp:function:: void Zeros( int m, int n, Matrix<T>& A )
.. cpp:function:: void Zeros( int m, int n, DistMatrix<T,U,V>& A )

   Set the matrix ``A`` to be an :math:`m \times n` matrix of all zeros. 

Change all entries of the matrix :math:`A` to zero.

.. cpp:function:: void MakeZeros( Matrix<T>& A )
.. cpp:function:: void MakeZeros( DistMatrix<T,U,V>& A )

   Change the entries of the matrix to zero.
