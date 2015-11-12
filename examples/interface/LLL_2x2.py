import El

# Form 
# 1 4
# 2 5
# 3 6
A = El.Matrix()
El.Zeros(A,3,2)
A.Set(0,0,1)
A.Set(1,0,2)
A.Set(2,0,3)
A.Set(0,1,4)
A.Set(1,1,5)
A.Set(2,1,6)
El.Print(A)

B=El.Matrix()
El.Copy( A, B )
QR, numBacktrack = El.LLL(B,0.75)
El.Print(B)
