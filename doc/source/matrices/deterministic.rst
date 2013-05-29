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

.. cpp:function:: void Cauchy( Matrix<F>& A, const std::vector<F>& x, const std::vector<F>& y )
.. cpp:function:: void Cauchy( DistMatrix<F,U,V>& A, const std::vector<F>& x, const std::vector<F>& y )

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

.. cpp:function:: void CauchyLike( Matrix<F>& A, const std::vector<F>& r, const std::vector<F>& s, const std::vector<F>& x, const std::vector<F>& y )
.. cpp:function:: void CauchyLike( DistMatrix<F,U,V>& A, const std::vector<F>& r, const std::vector<F>& s, const std::vector<F>& x, const std::vector<F>& y )

   Generate a Cauchy-like matrix using the defining vectors: :math:`r`, :math:`s`, :math:`x`, and :math:`y`.

Circulant
---------
An :math:`n \times n` matrix :math:`A` is called *circulant* if there exists a vector :math:`b` 
such that 

.. math::

   \alpha_{i,j} = \beta_{(i-j) \bmod n},

where :math:`\beta_k` is the :math:`k`'th entry of vector :math:`b`.

.. cpp:function:: void Circulant( Matrix<T>& A, const std::vector<T>& a )
.. cpp:function:: void Circulant( DistMatrix<T,U,V>& A, const std::vector<T>& a )

   Generate a circulant matrix using the vector ``a``.

Diagonal
--------
An :math:`n \times n` matrix :math:`A` is called *diagonal* if each entry :math:`(i,j)`, where 
:math:`i \neq j`, is :math:`0`. They are therefore defined by the *diagonal* values, where 
:math:`i = j`.

.. cpp:function:: void Diagonal( Matrix<T>& D, const std::vector<T>& d )
.. cpp:function:: void Diagonal( DistMatrix<T,U,V>& D, const std::vector<T>& d )

   Construct a diagonal matrix from the vector of diagonal values, :math:`d`.

Egorov
------
**TODO**

Extended Kahan
--------------
**TODO**

Fiedler
-------

Forsythe
--------
**TODO**

Fourier
-------
The :math:`n \times n` *Discrete Fourier Transform* (DFT) matrix, say :math:`A`, is given by

.. math::

   \alpha_{i,j} = \frac{e^{-2\pi i j / n}}{\sqrt{n}}.

.. cpp:function:: void Fourier( Matrix<Complex<R> >& A, int n )
.. cpp:function:: void Fourier( DistMatrix<Complex<R>,U,V>& A, int n )

   Set the matrix ``A`` equal to the :math:`n \times n` DFT matrix.

.. cpp:function:: void MakeFourier( Matrix<Complex<R> >& A )
.. cpp:function:: void MakeFourier( DistMatrix<Complex<R>,U,V>& A )

   Turn the existing :math:`n \times n` matrix ``A`` into a DFT matrix.

GCDMatrix
---------
**TODO**

Gear
----
**TODO**

Golub/Klema/Stewart 
-------------------
**TODO**

Grcar
-----
**TODO**

Hankel
------
An :math:`m \times n` matrix :math:`A` is called a *Hankel matrix* if there 
exists a vector :math:`b` such that

.. math::

   \alpha_{i,j} = \beta_{i+j},

where :math:`\alpha_{i,j}` is the :math:`(i,j)` entry of :math:`A` and 
:math:`\beta_k` is the :math:`k`'th entry of the vector :math:`b`.

.. cpp:function:: void Hankel( Matrix<T>& A, int m, int n, const std::vector<T>& b )
.. cpp:function:: void Hankel( DistMatrix<T,U,V>& A, int m, int n, const std::vector<T>& b )

   Create an :math:`m \times n` Hankel matrix from the generate vector, 
   :math:`b`.

Hanowa
------
**TODO**

Helmholtz
---------
**TODO**

Hilbert
-------
The Hilbert matrix of order :math:`n` is the :math:`n \times n` matrix where
entry :math:`(i,j)` is equal to :math:`1/(i+j+1)`.

.. cpp:function:: void Hilbert( Matrix<F>& A, int n )
.. cpp:function:: void Hilbert( DistMatrix<F,U,V>& A, int n )

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

.. cpp:function:: void Identity( Matrix<T>& A, int m, int n )
.. cpp:function:: void Identity( DistMatrix<T,U,V>& A, int m, int n )

   Set the matrix ``A`` equal to the :math:`m \times n` identity(-like) matrix.

.. cpp:function:: void MakeIdentity( Matrix<T>& A )
.. cpp:function:: void MakeIdentity( DistMatrix<T,U,V>& A ) 

   Set the matrix ``A`` to be identity-like.

Jordan
------
**TODO**

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

.. cpp:function:: void Kahan( Matrix<F>& A, int n, F phi )
.. cpp:function:: void Kahan( DistMatrix<F>& A, int n, F phi )

   Sets the matrix ``A`` equal to the :math:`n \times n` Kahan matrix with 
   the specified value for :math:`\phi`.

KMS
---
**TODO**

Laplacian
---------
**TODO**

Lauchli
-------
**TODO**

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

.. cpp:function:: void Legendre( Matrix<F>& A, int n )
.. cpp:function:: void Legendre( DistMatrix<F,U,V>& A, int n )

   Sets the matrix ``A`` equal to the :math:`n \times n` Jacobi matrix.

Lehmer
------
**TODO**

Lotkin
------
**TODO**

MinIJ
-----
**TODO**

Ones
----
Create an :math:`m \times n` matrix of all ones.

.. cpp:function:: void Ones( Matrix<T>& A, int m, int n )
.. cpp:function:: void Ones( DistMatrix<T,U,V>& A, int m, int n )

   Set the matrix ``A`` to be an :math:`m \times n` matrix of all ones.

Change all entries of the matrix :math:`A` to one.

.. cpp:function:: void MakeOnes( Matrix<T>& A )
.. cpp:function:: void MakeOnes( DistMatrix<T,U,V>& A )

   Change the entries of the matrix to ones.

OneTwoOne
---------
A "1-2-1" matrix is tridiagonal with a diagonal of all twos and sub- and 
super-diagonals of all ones.

.. cpp:function:: void OneTwoOne( Matrix<T>& A, int n )
.. cpp:function:: void OneTwoOne( DistMatrix<T,U,V>& A, int n )

   Set ``A`` to a :math:`n \times n` "1-2-1" matrix.

.. cpp:function:: void MakeOneTwoOne( Matrix<T>& A )
.. cpp:function:: void MakeOneTwoOne( DistMatrix<T,U,V>& A )

   Modify the entries of the square matrix ``A`` to be "1-2-1".

Parter
------
**TODO**

Pei
---
**TODO**

Redheffer
---------
**TODO**

Riemann
-------
**TODO**

Ris
---
**TODO**

Toeplitz
--------
An :math:`m \times n` matrix is *Toeplitz* if there exists a vector :math:`b` such that, for each entry :math:`\alpha_{i,j}` of :math:`A`,

.. math::

   \alpha_{i,j} = \beta_{i-j+(n-1)},

where :math:`\beta_k` is the :math:`k`'th entry of :math:`b`.

.. cpp:function:: void Toeplitz( Matrix<T>& A, int m, int n, const std::vector<T>& b )
.. cpp:function:: void Toeplitz( DistMatrix<T,U,V>& A, int m, int n, const std::vector<T>& b )

   Build the matrix ``A`` using the generating vector :math:`b`.

TriW
----
**TODO**

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

.. cpp:function:: void Walsh( Matrix<T>& W, int k, bool binary=false )
.. cpp:function:: void Walsh( DistMatrix<T,U,V>& W, int k, bool binary=false )

   Set the matrix :math:`W` equal to the :math:`k`'th (possibly binary) Walsh 
   matrix.

Wilkinson
---------
A *Wilkinson matrix* of order :math:`k` is a tridiagonal matrix with diagonal

.. math::

   [k,k-1,k-2,...,1,0,1,...,k-2,k-1,k],

and sub- and super-diagonals of all ones.

.. cpp:function:: void Wilkinson( Matrix<T>& W, int k )
.. cpp:function:: void Wilkinson( DistMatrix<T,U,V>& W, int k )

   Set the matrix :math:`W` equal to the :math:`k`'th Wilkinson matrix.

Zeros
-----
Create an :math:`m \times n` matrix of all zeros.

.. cpp:function:: void Zeros( Matrix<T>& A, int m, int n )
.. cpp:function:: void Zeros( DistMatrix<T,U,V>& A, int m, int n )

   Set the matrix ``A`` to be an :math:`m \times n` matrix of all zeros. 

Change all entries of the matrix :math:`A` to zero.

.. cpp:function:: void MakeZeros( Matrix<T>& A )
.. cpp:function:: void MakeZeros( DistMatrix<T,U,V>& A )

   Change the entries of the matrix to zero.
