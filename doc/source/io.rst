Input/output
************

Display
-------

.. cpp:function:: void Display( const Matrix<T>& A, std::string title="Default" )
.. cpp:function:: void Display( const DistMatrix<T,U,V>& A, std::string title="Default" )

   If Qt5 was detected during configuration, display the matrix on screen.
   Otherwise, print it to the console.

Print
-----

.. cpp:function:: void Print( const Matrix<T>& A, std::string title="", std::ostream& os=std::cout )
.. cpp:function:: void Print( const DistMatrix<T,U,V>& A, std::string title="", std::ostream& os=std::cout )

   Prints the matrix to the console.

Spy
---

.. cpp:function:: void Spy( const Matrix<T>& A, std::string title="Default", typename Base<T>::type tol=0 )
.. cpp:function:: void Spy( const DistMatrix<T,U,V>& A, std::string title="Default", typename Base<T>::type tol=0 )

   Only available if Qt5 was detected during configuration. A spy plot of 
   the elements with absolute values greater than or equal to ``tol`` is 
   displayed on screen.

Write
-----

.. cpp:function:: void Write( const Matrix<T>& A, std::string title="", std::string filename="Matrix" )
.. cpp:function:: void Write( const DistMatrix<T>& A, std::string title="", std::string filename="DistMatrix" )

   The matrix is printed to the file with the specified name.
