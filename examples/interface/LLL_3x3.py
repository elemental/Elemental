import El

loopTol = 0.
progress = False
output = True

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

if output:
  El.Print(A,"A")

B=El.Matrix()

for presorted, smallestFirst in (True,False), (True,True), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    El.Copy( A, B )
    QR, numBacktrack = \
      El.LLL(B,deltaLower,loopTol,presorted,smallestFirst,progress)
    if output:
      El.Print(B,"B(0.75)")
    delta = El.LLLDelta(QR)
    print "delta=", delta
    print "num backtracks=", numBacktrack
    print ""
