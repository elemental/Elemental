import El

# Form 
# -6   9 -15 -18 
#  4  -6  10  12 
# 10 -15  18  35 
#-24  36 -46 -82 
A = El.Matrix()
El.Zeros(A,4,4)
A.Set(0,0,-6)
A.Set(1,0,+4)
A.Set(2,0,10)
A.Set(3,0,-24)
A.Set(0,1,9)
A.Set(1,1,-6)
A.Set(2,1,15)
A.Set(3,1,36)
A.Set(2,1,-15)
A.Set(0,2,-15)
A.Set(1,2,10)
A.Set(2,2,18)
A.Set(3,2,-46)
A.Set(0,3,-18)
A.Set(1,3,12)
A.Set(2,3,35)
A.Set(3,3,-82)
El.Print(A)

B=El.Matrix()

for presorted, smallestFirst in (True,True), (True,False), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    print "Testing with presorted=%r, smallestFirst=%r, deltaLower=%f" % \
      (presorted,smallestFirst,deltaLower)
    El.Copy( A, B )
    El.LLL(B,deltaLower,0.,presorted,smallestFirst)
