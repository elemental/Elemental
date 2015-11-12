import El

m = 2000
n = 1000
loopTol = 0.
progress = False
output = False

A = El.Matrix(El.dTag)
El.Uniform( A, m, n, 0., 10.)
El.Round( A )
if output:
  El.Print( A, "A" )

B=El.Matrix(El.dTag)

for presorted, smallestFirst in (True,False), (True,True), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    El.Copy( A, B )
    print "Testing with presorted=%r, smallestFirst=%r, deltaLower=%f" % \
      (presorted,smallestFirst,deltaLower)
    startTime = El.mpi.Time()
    QR, numBacktrack = \
      El.LLL(B,deltaLower,loopTol,presorted,smallestFirst,progress)
    runTime = El.mpi.Time() - startTime
    print "  runtime: %f seconds" % runTime 
    b1Norm = El.FrobeniusNorm(B[:,0])
    print "  || b_1 ||_2 =", b1Norm

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

    if output:
      El.Print(B,"B")
    print "  num backtracks: ", numBacktrack
    delta = El.LLLDelta(QR)
    print "  delta=", delta
    print ""
