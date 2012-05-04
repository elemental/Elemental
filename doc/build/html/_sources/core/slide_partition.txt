Sliding partitions
==================

SlidePartitionUp
----------------
Simultaneously slide and merge the partition

.. math::

   A = \left(\begin{array}{c} A_0 \\ A_1 \\ \hline A_2 \end{array}\right),

into

.. math::

   \left(\begin{array}{c} A_T \\ \hline A_B \end{array}\right) = 
   \left(\begin{array}{c} A_0 \\ \hline A_1 \\ A_2 \end{array}\right).

.. cpp:function:: void SlidePartitionUp( Matrix<T>& AT, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& AB, Matrix<T>& A2 )

.. cpp:function:: void SlideLockedPartitionUp( Matrix<T>& AT, const Matrix<T>& A0, const Matrix<T>& A1, Matrix<T>& AB, const Matrix<T>& A2 )

   Templated over the datatype, `T`.

.. cpp:function:: void SlidePartitionUp( DistMatrix<T,U,V>& AT, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& AB, DistMatrix<T,U,V>& A2 )

.. cpp:function:: void SlideLockedPartitionUp( DistMatrix<T,U,V>& AT, const DistMatrix<T,U,V>& A0, const DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& AB, const DistMatrix<T,U,V>& A2 )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   SlidePartitionUp( AT,  A0,
                    /**/ /**/
                          A1,
                     AB,  A2 );

SlidePartitionDown
------------------
Simultaneously slide and merge the partition

.. math::

   A = \left(\begin{array}{c} A_0 \\ \hline A_1 \\ A_2 \end{array}\right),

into

.. math::

   \left(\begin{array}{c} A_T \\ \hline A_B \end{array}\right) = 
   \left(\begin{array}{c} A_0 \\ A_1 \\ \hline A_2 \end{array}\right).

.. cpp:function:: void SlidePartitionDown( Matrix<T>& AT, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& AB, Matrix<T>& A2 )

.. cpp:function:: void SlideLockedPartitionDown( Matrix<T>& AT, const Matrix<T>& A0, const Matrix<T>& A1, Matrix<T>& AB, const Matrix<T>& A2 )

   Templated over the datatype, `T`.

.. cpp:function:: void SlidePartitionDown( DistMatrix<T,U,V>& AT, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& AB, DistMatrix<T,U,V>& A2 )

.. cpp:function:: void SlideLockedPartitionDown( DistMatrix<T,U,V>& AT, const DistMatrix<T,U,V>& A0, const DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& AB, const DistMatrix<T,U,V>& A2 )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   SlidePartitionDown( AT,  A0,
                            A1,
                      /**/ /**/
                       AB,  A2 );

SlidePartitionLeft
------------------
Simultaneously slide and merge the partition

.. math::

   A = \left(\begin{array}{cc|c} A_0 & A_1 & A_2 \end{array}\right)

into

.. math::

   \left(\begin{array}{c|c} A_L & A_R \end{array}\right) = 
   \left(\begin{array}{c|cc} A_0 & A_1 & A_2 \end{array}\right).

.. cpp:function:: void SlidePartitionLeft( Matrix<T>& AL, Matrix<T>& AR, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 )

.. cpp:function:: void SlidePartitionLeft( DistMatrix<T,U,V>& AL, DistMatrix<T,U,V>& AR, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& A2 )

   Templated over the datatype, `T`.

.. cpp:function:: void SlideLockedPartitionLeft( Matrix<T>& AL, Matrix<T>& AR, const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 )

.. cpp:function:: void SlideLockedPartitionLeft( DistMatrix<T,U,V>& AL, DistMatrix<T,U,V>& AR, const DistMatrix<T,U,V>& A0, const DistMatrix<T,U,V>& A1, const DistMatrix<T,U,V>& A2 )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   SlidePartitionLeft( AL, /**/ AR,
                       A0, /**/ A1, A2 );

SlidePartitionRight
-------------------
Simultaneously slide and merge the partition

.. math::

   A = \left(\begin{array}{c|cc} A_0 & A_1 & A_2 \end{array}\right)

into

.. math::

   \left(\begin{array}{c|c} A_L & A_R \end{array}\right) = 
   \left(\begin{array}{cc|c} A_0 & A_1 & A_2 \end{array}\right).

.. cpp:function:: void SlidePartitionRight( Matrix<T>& AL, Matrix<T>& AR, Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 )

.. cpp:function:: void SlidePartitionRight( DistMatrix<T,U,V>& AL, DistMatrix<T,U,V>& AR, DistMatrix<T,U,V>& A0, DistMatrix<T,U,V>& A1, DistMatrix<T,U,V>& A2 )

   Templated over the datatype, `T`.

.. cpp:function:: void SlideLockedPartitionRight( Matrix<T>& AL, Matrix<T>& AR, const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 )

.. cpp:function:: void SlideLockedPartitionRight( DistMatrix<T,U,V>& AL, DistMatrix<T,U,V>& AR, const DistMatrix<T,U,V>& A0, const DistMatrix<T,U,V>& A1, const DistMatrix<T,U,V>& A2 )

   Templated over the datatype, `T`, and distribution scheme, `(U,V)`.

Note that each of the above routines is meant to be used in a manner similar 
to the following:

.. code-block:: cpp

   SlidePartitionRight( AL,     /**/ AR,
                        A0, A1, /**/ A2 );

SlidePartitionUpDiagonal
------------------------
Simultaneously slide and merge the partition

.. math::

   A = \left(\begin{array}{cc|c} A_{00} & A_{01} & A_{02} \\
                                 A_{10} & A_{11} & A_{12} \\
                                 \hline
                                 A_{20} & A_{21} & A_{22} \end{array}\right)

into

.. math::

   \left(\begin{array}{c|c} A_{TL} & A_{TR} \\
                            \hline
                            A_{BL} & A_{BR} \end{array}\right) = 
   \left(\begin{array}{c|cc} A_{00} & A_{01} & A_{02} \\
                             \hline
                             A_{10} & A_{11} & A_{12} \\
                             A_{20} & A_{21} & A_{22} \end{array}\right).

Note that the above routines are meant to be used as:

.. code-block:: cpp

   SlidePartitionUpDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                            /*************/ /******************/
                                  /**/       A10, /**/ A11, A12,
                             ABL, /**/ ABR,  A20, /**/ A21, A22 );

SlidePartitionDownDiagonal
--------------------------
Simultaneously slide and merge the partition

.. math::

   A = \left(\begin{array}{c|cc} A_{00} & A_{01} & A_{02} \\
                                 \hline
                                 A_{10} & A_{11} & A_{12} \\
                                 A_{20} & A_{21} & A_{22} \end{array}\right)

into

.. math::

   \left(\begin{array}{c|c} A_{TL} & A_{TR} \\
                            \hline
                            A_{BL} & A_{BR} \end{array}\right) = 
   \left(\begin{array}{cc|c} A_{00} & A_{01} & A_{02} \\
                             A_{10} & A_{11} & A_{12} \\
                             \hline
                             A_{20} & A_{21} & A_{22} \end{array}\right).

Note that the above routines are meant to be used as:

.. code-block:: cpp

   SlidePartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                    /**/       A10, A11, /**/ A12,
                              /*************/ /******************/
                               ABL, /**/ ABR,  A20, A21, /**/ A22 );


