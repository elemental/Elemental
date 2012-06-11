High-level linear algebra
*************************

This chapter describes all of the linear algebra operations which are not basic
enough to fall within the BLAS (Basic Linear Algebra Subprograms) framework.
In particular, algorithms which would traditionally have fallen into the 
domain of LAPACK (Linear Algebra PACKage), such as factorizations and 
eigensolvers, are placed here.

.. toctree::
   :maxdepth: 2

   lapack-like/norms
   lapack-like/invariants
   lapack-like/factorizations
   lapack-like/solvers
   lapack-like/inversion
   lapack-like/condensed
   lapack-like/eigen_svd
   lapack-like/functions
   lapack-like/util
   lapack-like/tuning
