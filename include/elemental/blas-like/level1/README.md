### `include/elemental/blas-like/level1/`

The level-1 BLAS-like routines implemented in this folder are:

-  `Adjoint.hpp`: form the adjoint (conjugate-transpose) of a matrix
-  `Axpy.hpp`: Y becomes alpha X + Y (AXPY)
-  `AxpyTriangle.hpp`: Axpy a triangular matrix onto a triangular matrix
-  `Conjugate.hpp`: form the conjugate of a matrix
-  `Copy.hpp`: form a copy of a matrix
-  `DiagonalScale.hpp`: apply a diagonal matrix
-  `DiagonalSolve.hpp`: apply the inverse of a diagonal matrix
-  `Dot.hpp`: form the dot product of two vectors
-  `Dotu.hpp`: form the unconjugated dot product of two vectors
-  `MakeHermitian.hpp`: force a matrix to be Hermitian
-  `MakeReal.hpp`: force a matrix to be real
-  `MakeSymmetric.hpp`: force a matrix to be symmetric
-  `MakeTrapezoidal.hpp`: force a matrix to be trapezoidal
-  `MakeTriangular.hpp`: force a matrix to be triangular
-  `Max.hpp`: return the maximum entry in a vector or matrix
-  `Nrm2.hpp`: compute the two-norm of a vector
-  `QuasiDiagonalScale.hpp`: apply a quasi-diagonal matrix
-  `QuasiDiagonalSolve.hpp`: apply the inverse of a quasi-diagonal matrix
-  `Scale.hpp`: scale a matrix or vector 
-  `ScaleTrapezoid.hpp`: scale a trapezoid of a matrix or vector
-  `SetDiagonal.hpp`: set the entire diagonal of a matrix to a particular value
-  `Swap.hpp`: swap the contents of two matrices, or swap two rows or columns of
   a particular matrix
-  `Symmetric2x2Inv.hpp`: form the inverse of a symmetric 2x2 matrix
-  `Symmetric2x2Scale.hpp`: apply a symmetric 2x2 matrix
-  `Symmetric2x2Solve.hpp`: apply the inverse of a symmetric 2x2 matrix
-  `Transpose.hpp`: form the transpose of a matrix
-  `UpdateDiagonal.hpp`: add a constant value to every diagonal entry of a 
   matrix
-  `Zero.hpp`: set every entry in a matrix to zero
