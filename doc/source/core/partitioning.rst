Partitioning
============
The following routines are slight tweaks of the FLAME project's 
(as well as PLAPACK's) approach to submatrix tracking; the difference is that 
they have "locked" versions, which are meant for creating partitionings where 
the submatrices cannot be modified.

.. cpp:function:: void PartitionUp( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, int heightAB=Blocksize() )

   **Left off here**
