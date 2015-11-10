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

El.Print(A,"A")

B=El.Matrix()

El.Copy( A, B )
QR = El.LLL(B,0.75)
El.Print(B,"B(0.75)")
delta = El.LLLDelta(QR)
print "delta=", delta
print ""

El.Copy( A, B )
QR = El.LLL(B,0.5)
El.Print(B,"B(0.5)")
delta = El.LLLDelta(QR)
print "delta=", delta

