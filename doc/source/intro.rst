Elemental is a high-performance framework for distributed-memory dense linear 
algebra that is inspired by 
`PLAPACK <http://www.cs.utexas.edu/users/plapack/new/using.html>`_
and, more recently, 
`FLAME <http://zold.cs.utexas.edu/wiki/flame.wiki/FrontPage>`_. The main goal 
of the project is to provide a convenient interface to efficient distributed 
versions of as much of BLAS and LAPACK-level functionality as possible. 
Elemental templates routines over the most general datatypes possible 
(e.g., Cholesky factorization is implemented in terms of a general field),
so the datatype-specific parameter of BLAS routines can be removed, e.g., 
``{S,D,C,Z}GEMM`` simply become ``Gemm``. Due to several fundamental 
differences in the implementation details of LAPACK-level routines, no attempt 
has been made to mirror LAPACK naming conventions. For instance, the LAPACK 
Cholesky factorization routines are named ``{S,D,C,Z}POTRF``, for 
"POsitive-definite TRiangular Factorization", but in Elemental they are simply 
named ``Chol``.

**Need to discuss different matrix distributions here**

Elemental is made available under the 
`New BSD license <http://www.opensource.org/licenses/bsd-license.php>`_.
If this license is not sufficient for your purposes, please contact Jack Poulson
at ``jack.poulson@gmail.com``.
