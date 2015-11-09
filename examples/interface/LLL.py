import El;

# Form 
# 1 -1 3
# 1  0 5
# 1  2 6
A = El.Matrix()
El.Zeros(A,3,3)

A.Set(0,0,1)
A.Set(1,0,1)
A.Set(2,0,1)

A.Set(0,1,-1)
A.Set(1,1, 0)
A.Set(2,1, 2)

A.Set(0,2,3)
A.Set(1,2,5)
A.Set(2,2,6)

El.Print(A)

B=El.Matrix()

El.Copy( A, B )
El.LLL(B,0.75,0.5001,0.,0.001)
El.Print(B)

El.Copy( A, B )
El.LLL(B,0.5,0.5001,0.,0.001)
El.Print(B)
