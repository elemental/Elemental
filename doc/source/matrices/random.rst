Random
======

Uniform
-------
We call an :math:`m \times n` matrix uniformly random if each entry is drawn 
from a uniform distribution over some ball :math:`B_r(x)`, which is centered 
around some point :math:`x` and of radius :math:`r`.

.. cpp:function:: void Uniform( int m, int n, Matrix<T>& A, T center=0, typename Base<T>::type radius=1 )

   Set the serial matrix ``A`` to an :math:`m \times n` matrix with each entry sampled from the uniform distribution centered at `center` with radius `radius` (templated over datatype, `T`).

.. cpp:function:: void Uniform( int m, int n, DistMatrix<T,U,V>& A, T center=0, typename Base<T>::type radius=1 )

   Set the distributed matrix ``A`` to an :math:`m \times n` matrix with each entry sampled from the uniform distribution centered at `center` with radius `radius` (templated over datatype, `T`, and distribution scheme, `(U,V)`).

.. cpp:function:: void MakeUniform( Matrix<T>& A, T center=0, typename Base<T>::type radius=1 )

   Sample each entry of ``A`` from :math:`U(B_r(x))`, where :math:`r` is given by ``radius`` and :math:`x` is given by ``center`` (templated over the datatype, `T`).

.. cpp:function:: void MakeUniform( DistMatrix<T,U,V>& A, T center=0, typename Base<T>::type radius=1 )

   Sample each entry of ``A`` from :math:`U(B_r(x))`, where :math:`r` is given by ``radius`` and :math:`x` is given by ``center`` (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

HermitianUniformSpectrum
------------------------
These routines sample a diagonal matrix from the specified interval of the 
real line and then perform a similarity transformation using a random 
Householder transform.

.. cpp:function:: void HermitianUniformSpectrum( int n, Matrix<F>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )

   Build the :math:`n \times n` serial matrix ``A`` with a spectrum sampled uniformly from the interval :math:`(lower,upper]` (templated over the datatype, `F`).

.. cpp:function:: void HermitianUniformSpectrum( int n, DistMatrix<F,U,V>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )

   Build the :math:`n \times n` distributed matrix ``A`` with a spectrum sampled uniformly from the interval :math:`(lower,upper]` (templated over the datatype, `F`, which must be a field, and the distribution scheme, `(U,V)`).

.. cpp:function:: void MakeHermitianUniformSpectrum( Matrix<F>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )

   Sample the entries of the square serial matrix ``A`` from the interval :math:`(lower,upper]` (templated over the datatype, `F`).

.. cpp:function:: void MakeHermitianUniformSpectrum( DistMatrix<F,U,V>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )

   Sample the entries of the square distributed matrix ``A`` from the interval :math:`(lower,upper]` (templated over the datatype, `F`, and the distribution scheme, `(U,V)`).

NormalUniformSpectrum
---------------------
These routines sample a diagonal matrix from the specified ball in the 
complex plane and then perform a similarity transformation using a random 
Householder transform.

.. cpp:function:: void NormalUniformSpectrum( int n, Matrix<Complex<R> >& A, Complex<R> center=0, R radius=1 )

   Build the :math:`n \times n` serial matrix ``A`` with a spectrum sampled uniformly from the ball :math:`B_{\mathrm{radius}}(\mathrm{center})` (templated over the real datatype, `R`).

.. cpp:function:: void NormalUniformSpectrum( int n, DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 )

   Build the :math:`n \times n` distributed matrix ``A`` with a spectrum sampled uniformly from the ball :math:`B_{\mathrm{radius}}(\mathrm{center})` (templated over the real datatype, `R`, and the distribution scheme, `(U,V)`).

.. cpp:function:: void MakeNormalUniformSpectrum( Matrix<Complex<R> >& A, Complex<R> center=0, R radius=1 )

   Sample the entries of the square serial matrix ``A`` from the ball in the complex plane centered at ``center`` with radius ``radius`` (templated over the real datatype, `R`).

.. cpp:function:: void MakeNormalUniformSpectrum( DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 )

   Sample the entries of the square distributed matrix ``A`` from the ball in the complex plane centered at ``center`` with radius ``radius`` (templated over the real datatype, `R`, and the distribution scheme, `(U,V)`).
