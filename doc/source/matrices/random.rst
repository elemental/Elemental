Random
======

Gaussian
--------
An :math:`m \times n` matrix is Gaussian if each entry is independently drawn
from a normal distribution.

.. cpp:function:: void Gaussian( Matrix<T>& A, int m, int n, T mean=0, typename Base<T>::type stddev=1 )
.. cpp:function:: void Gaussian( DistMatrix<T,U,V>& A, int m, int n, T mean=0, typename Base<T>::type stddev=1 )

   Sets the matrix ``A`` to an :math:`m \times n` Gaussian matrix with the
   specified mean and standard deviation.

.. cpp:function:: void MakeGaussian( Matrix<T>& A, T mean=0, typename Base<T>::type stddev=1 )
.. cpp:function:: void MakeGaussian( DistMatrix<T,U,V>& A, T mean=0, typename Base<T>::type stddev=1 )

   Changes each entry to an independent sample from the specified normal
   distribution.

Wigner
------
A Hermitian matrix whose entries in one triangle are all independent samples
from a normal distribution. The spectra of these matrices are well-studied.

.. cpp:function:: void Wigner( Matrix<T>& A, int n, T mean=0, typename Base<T>::type stddev=1 )
.. cpp:function:: void Wigner( DistMatrix<T,U,V>& A, int n, T mean=0, typename Base<T>::type stddev=1 )

   Sets the matrix ``A`` to an :math:`n \times n` Wigner matrix with the
   specified mean and standard deviation.

Haar
----
The Haar distribution is the uniform distribution over the space of real or 
complex unitary matrices. 

.. cpp:function:: void Haar( Matrix<F>& A, int n )
.. cpp:function:: void Haar( DistMatrix<F>& A, int n )

   Draws ``A`` from the Haar distribution. The current scheme performs a QR
   factorization of a Gaussian matrix, but Stewart introduced a well-known 
   scheme which only requires quadratic work for the implicit representation 
   as a product of random Householder reflectors.

.. cpp:function:: void ImplicitHaar( Matrix<F>& A, Matrix<F>& t, int n )
.. cpp:function:: void ImplicitHaar( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, int n )

   Sets ``A`` to a set of Householder reflectors with the same structure as
   the result of a QR decomposition. The product of these reflectors is a 
   sample from the Haar distribution.

Uniform
-------
We call an :math:`m \times n` matrix is uniformly random if each entry is drawn 
from a uniform distribution over some ball :math:`B_r(x)`, which is centered 
around some point :math:`x` and of radius :math:`r`.

.. cpp:function:: void Uniform( Matrix<T>& A, int m, int n, T center=0, typename Base<T>::type radius=1 )
.. cpp:function:: void Uniform( DistMatrix<T,U,V>& A, int m, int n, T center=0, typename Base<T>::type radius=1 )

   Set the matrix ``A`` to an :math:`m \times n` matrix with each entry sampled from the uniform distribution centered at `center` with radius `radius`.

.. cpp:function:: void MakeUniform( Matrix<T>& A, T center=0, typename Base<T>::type radius=1 )
.. cpp:function:: void MakeUniform( DistMatrix<T,U,V>& A, T center=0, typename Base<T>::type radius=1 )

   Sample each entry of ``A`` from :math:`U(B_r(x))`, where :math:`r` is given by ``radius`` and :math:`x` is given by ``center``.

HermitianUniformSpectrum
------------------------
These routines sample a diagonal matrix from the specified interval of the 
real line and then perform a similarity transformation using a random 
Householder transform.

.. cpp:function:: void HermitianUniformSpectrum( Matrix<F>& A, int n, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )
.. cpp:function:: void HermitianUniformSpectrum( DistMatrix<F,U,V>& A, int n, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )

   Build the :math:`n \times n` matrix ``A`` with a spectrum sampled uniformly 
   from the interval :math:`(lower,upper]`.

.. cpp:function:: void MakeHermitianUniformSpectrum( Matrix<F>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )
.. cpp:function:: void MakeHermitianUniformSpectrum( DistMatrix<F,U,V>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )

   Sample the entries of the square matrix ``A`` from the interval 
   :math:`(lower,upper]`.

NormalUniformSpectrum
---------------------
These routines sample a diagonal matrix from the specified ball in the 
complex plane and then perform a similarity transformation using a random 
Householder transform.

.. cpp:function:: void NormalUniformSpectrum( Matrix<Complex<R> >& A, int n, Complex<R> center=0, R radius=1 )
.. cpp:function:: void NormalUniformSpectrum( DistMatrix<Complex<R>,U,V>& A, int n, Complex<R> center=0, R radius=1 )

   Build the :math:`n \times n` matrix ``A`` with a spectrum sampled uniformly 
   from the ball :math:`B_{\mathrm{radius}}(\mathrm{center})`.

.. cpp:function:: void MakeNormalUniformSpectrum( Matrix<Complex<R> >& A, Complex<R> center=0, R radius=1 )
.. cpp:function:: void MakeNormalUniformSpectrum( DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 )

   Sample the entries of the square matrix ``A`` from the ball in the complex 
   plane centered at ``center`` with radius ``radius``.
