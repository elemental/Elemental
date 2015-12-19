import El

# Form 
# -6   9 -15 -18 
#  4  -6  10  12 
# 10 -15  18  35 
#-24  36 -46 -82 
A = El.Matrix()
n=4
El.Zeros(A,n,n)
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

loopTol=0.
full=True
weak=False
progress=False
output=True

if output:
  El.Print( A, "B" )

B=El.Matrix()
for presort, smallestFirst in (True,True), (True,False), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    print "Testing with presort=%r, smallestFirst=%r, deltaLower=%f" % \
      (presort,smallestFirst,deltaLower)
    El.Copy( A, B )
    startTime = El.mpi.Time()
    U, UInv, QR, nullity, numBacktracks = \
      El.LLL(B,deltaLower,loopTol,full=full,weak=weak,presort=presort,
        smallestFirst=smallestFirst,progress=progress)
    runTime = El.mpi.Time() - startTime
    print "  runtime: %f seconds" % runTime
    b1Norm = El.FrobeniusNorm(B[:,0])
    print "  || b_1 ||_2 = ",b1Norm

    if output:
      El.Print( U, "U" )
      El.Print( UInv, "UInv" )
      El.Print( B, "BNew" )
      El.Print( QR, "QR" )

    numLower = 0
    minInd = 0
    minNorm = b1Norm
    for j in xrange(1,n):
      bjNorm = El.FrobeniusNorm(B[:,j])
      if bjNorm < b1Norm:
        numLower = numLower + 1
      if bjNorm < minNorm:
        minNorm = bjNorm
        minInd = j
    if numLower == 0:
      print "  b1 was the smallest column"
    else:
      print "  %d columns were smaller than b1, with || b_%d ||_2 = %f the smallest" % (numLower,minInd,minNorm)

    print "  nullity: ", nullity
    print "  num backtracks: ", numBacktracks

    delta = El.LLLDelta(QR,weak) 
    print "  delta=", delta
    print ""
