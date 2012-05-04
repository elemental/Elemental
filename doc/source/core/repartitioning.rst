Repartitioning
==============

RepartitionUp
-------------
Given the partition

.. math::

   A = \left(\begin{array}{c} A_T \\ \hline A_B \end{array}\right),

and a blocksize, :math:`n_b`, turn the two-way partition into the three-way
partition 

.. math::

   \left(\begin{array}{c} A_T \\ \hline A_B \end{array}\right) = 
   \left(\begin{array}{c} A_0 \\ A_1 \\ \hline A_2 \end{array}\right),

where :math:`A_1` is of height :math:`n_b` and :math:`A_2 = A_B`.

.. cpp:function:: void RepartitionUp( Matrix<T>& AT, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& AB, Matrix<T>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionUp( const Matrix<T>& AT, Matrix<T>& A0, Matrix<T>& A1, const Matrix<T>& AB, Matrix<T>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`.

.. cpp:function:: void RepartitionUp( DistMatrix<T,U,V>& AT, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& AB, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionUp( const DistMatrix<T,U,V>& AT, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, const DistMatrix<T,U,V>& AB, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   RepartitionUp( AT,  A0,
                       A1,
                 /**/ /**/
                  AB,  A2, blocksize );

RepartitionDown
---------------
Given the partition

.. math::

   A = \left(\begin{array}{c} A_T \\ \hline A_B \end{array}\right),

and a blocksize, :math:`n_b`, turn the two-way partition into the three-way
partition 

.. math::

   \left(\begin{array}{c} A_T \\ \hline A_B \end{array}\right) = 
   \left(\begin{array}{c} A_0 \\ \hline A_1 \\ A_2 \end{array}\right),

where :math:`A_1` is of height :math:`n_b` and :math:`A_0 = A_T`.

.. cpp:function:: void RepartitionDown( Matrix<T>& AT, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& AB, Matrix<T>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionDown( const Matrix<T>& AT, Matrix<T>& A0, Matrix<T>& A1, const Matrix<T>& AB, Matrix<T>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`.

.. cpp:function:: void RepartitionDown( DistMatrix<T,U,V>& AT, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& AB, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionDown( const DistMatrix<T,U,V>& AT, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, const DistMatrix<T,U,V>& AB, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   RepartitionDown( AT,  A0,
                   /**/ /**/
                         A1,
                    AB,  A2, blocksize );

RepartitionLeft
---------------
Given the partition

.. math::

   A = \left(\begin{array}{c|c} A_L & A_R \end{array}\right),

and a blocksize, :math:`n_b`, turn the two-way partition into the three-way 
partition

.. math::

   \left(\begin{array}{c|c} A_L & A_R \end{array}\right) = 
   \left(\begin{array}{cc|c} A_0 & A_1 & A_2 \end{array}\right),

where :math:`A_1` is of width :math:`n_b` and :math:`A_2=A_R`.

.. cpp:function:: void RepartitionLeft( Matrix<T>& AL, Matrix<T>& AR, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionLeft( const Matrix<T>& AL, const Matrix<T>& AR, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`.

.. cpp:function:: void RepartitionLeft( DistMatrix<T,U,V>& AL, DistMatrix<T,U,V>& AR, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionLeft( const DistMatrix<T,U,V>& AL, const DistMatrix<T,U,V>& AR, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   RepartitionLeft( AL,     /**/ AR,
                    A0, A1, /**/ A2, blocksize );

RepartitionRight
----------------
Given the partition

.. math::

   A = \left(\begin{array}{c|c} A_L & A_R \end{array}\right),

and a blocksize, :math:`n_b`, turn the two-way partition into the three-way 
partition

.. math::

   \left(\begin{array}{c|c} A_L & A_R \end{array}\right) = 
   \left(\begin{array}{c|cc} A_0 & A_1 & A_2 \end{array}\right),

where :math:`A_1` is of width :math:`n_b` and :math:`A_0=A_L`.

.. cpp:function:: void RepartitionRight( Matrix<T>& AL, Matrix<T>& AR, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionRight( const Matrix<T>& AL, const Matrix<T>& AR, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`.

.. cpp:function:: void RepartitionRight( DistMatrix<T,U,V>& AL, DistMatrix<T,U,V>& AR, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionRight( const DistMatrix<T,U,V>& AL, const DistMatrix<T,U,V>& AR, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& A2, int bsize=Blocksize() )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   RepartitionRight( AL, /**/ AR,
                     A0, /**/ A1, A2, blocksize );

RepartitionUpDiagonal
---------------------
Given the partition

.. math::

   A = \left(\begin{array}{c|c} A_{TL} & A_{TR} \\ \hline A_{BL} & A_{BR}
             \end{array}\right),

turn the two-by-two partition into the three-by-three partition

.. math::

   \left(\begin{array}{c|c} A_{TL} & A_{TR} \\ 
                            \hline
                            A_{BL} & A_{BR} \end{array}\right) = 
   \left(\begin{array}{cc|c} A_{00} & A_{01} & A_{02} \\ 
                             A_{10} & A_{11} & A_{12} \\
                             \hline
                             A_{20} & A_{21} & A_{22} \end{array}\right),

where :math:`A_{11}` is :math:`n_b \times n_b` and the corresponding quadrants are equivalent.

.. cpp:function:: void RepartitionUpDiagonal( Matrix<T>& ATL, Matrix<T>& ATR, Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, Matrix<T>& ABL, Matrix<T>& ABR, Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionUpDiagonal( const Matrix<T>& ATL, const Matrix<T>& ATR, Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, const Matrix<T>& ABL, const Matrix<T>& ABR, Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, int bsize=Blocksize() )

   Templated over the datatype, `T`.

.. cpp:function:: void RepartitionUpDiagonal( DistMatrix<T,U,V>& ATL, DistMatrix<T,U,V>& ATR, DistMatrix<T,U,V>& A00, DistMatrix<T,U,V>& A01, DistMatrix<T,U,V>& A02, DistMatrix<T,U,V>& A10, DistMatrix<T,U,V>& A11, DistMatrix<T,U,V>& A12, DistMatrix<T,U,V>& ABL, DistMatrix<T,U,V>& ABR, DistMatrix<T,U,V>& A20, DistMatrix<T,U,V>& A21, DistMatrix<T,U,V>& A22, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionUpDiagonal( const DistMatrix<T,U,V>& ATL, const DistMatrix<T,U,V>& ATR, DistMatrix<T,U,V>& A00, DistMatrix<T,U,V>& A01, DistMatrix<T,U,V>& A02, DistMatrix<T,U,V>& A10, DistMatrix<T,U,V>& A11, DistMatrix<T,U,V>& A12, const DistMatrix<T,U,V>& ABL, const DistMatrix<T,U,V>& ABR, DistMatrix<T,U,V>& A20, DistMatrix<T,U,V>& A21, DistMatrix<T,U,V>& A22, int bsize=Blocksize() )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   RepartitionUpDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                               /**/       A10, A11, /**/ A12,
                         /*************/ /******************/
                          ABL, /**/ ABR,  A20, A21, /**/ A22, blocksize );

RepartitionDownDiagonal
-----------------------
Given the partition

.. math::

   A = \left(\begin{array}{c|c} A_{TL} & A_{TR} \\ \hline A_{BL} & A_{BR}
             \end{array}\right),

turn the two-by-two partition into the three-by-three partition

.. math::

   \left(\begin{array}{c|c} A_{TL} & A_{TR} \\ 
                            \hline
                            A_{BL} & A_{BR} \end{array}\right) = 
   \left(\begin{array}{c|cc} A_{00} & A_{01} & A_{02} \\ 
                             \hline
                             A_{10} & A_{11} & A_{12} \\
                             A_{20} & A_{21} & A_{22} \end{array}\right),

where :math:`A_{11}` is :math:`n_b \times n_b` and the corresponding quadrants are equivalent.

.. cpp:function:: void RepartitionDownDiagonal( Matrix<T>& ATL, Matrix<T>& ATR, Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, Matrix<T>& ABL, Matrix<T>& ABR, Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionDownDiagonal( const Matrix<T>& ATL, const Matrix<T>& ATR, Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, const Matrix<T>& ABL, const Matrix<T>& ABR, Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, int bsize=Blocksize() )

   Templated over the datatype, `T`.

.. cpp:function:: void RepartitionDownDiagonal( DistMatrix<T,U,V>& ATL, DistMatrix<T,U,V>& ATR, DistMatrix<T,U,V>& A00, DistMatrix<T,U,V>& A01, DistMatrix<T,U,V>& A02, DistMatrix<T,U,V>& A10, DistMatrix<T,U,V>& A11, DistMatrix<T,U,V>& A12, DistMatrix<T,U,V>& ABL, DistMatrix<T,U,V>& ABR, DistMatrix<T,U,V>& A20, DistMatrix<T,U,V>& A21, DistMatrix<T,U,V>& A22, int bsize=Blocksize() )

.. cpp:function:: void LockedRepartitionDownDiagonal( const DistMatrix<T,U,V>& ATL, const DistMatrix<T,U,V>& ATR, DistMatrix<T,U,V>& A00, DistMatrix<T,U,V>& A01, DistMatrix<T,U,V>& A02, DistMatrix<T,U,V>& A10, DistMatrix<T,U,V>& A11, DistMatrix<T,U,V>& A12, const DistMatrix<T,U,V>& ABL, const DistMatrix<T,U,V>& ABR, DistMatrix<T,U,V>& A20, DistMatrix<T,U,V>& A21, DistMatrix<T,U,V>& A22, int bsize=Blocksize() )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   RepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                           /*************/ /******************/
                                 /**/       A10, /**/ A11, A12,
                            ABL, /**/ ABR,  A20, /**/ A21, A22, blocksize );

