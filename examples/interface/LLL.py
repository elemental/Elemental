import El

m = 50
n = 50

progress = False
timeLLL = False
output = False

BOrig = El.Matrix(El.dTag)
El.Uniform( BOrig, m, n, 0., 10.)
El.Round( BOrig )
if output:
  El.Print( BOrig, "B" )

ctrl = El.LLLCtrl_d()
ctrl.progress = progress
ctrl.time = timeLLL

B=El.Matrix(El.dTag)
for presort, smallestFirst in (True,False), (True,True), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    for variant in El.LLL_WEAK, El.LLL_NORMAL, El.LLL_DEEP, El.LLL_DEEP_REDUCE:

      print "variant=%d, presort=%r, smallest1st=%r, deltaLower=%f:" % \
        (variant,presort,smallestFirst,deltaLower)

      El.Copy( BOrig, B )

      ctrl.delta = deltaLower
      ctrl.variant = variant
      ctrl.presort = presort
      ctrl.smallestFirst = smallestFirst

      # Run the LLL reduction
      startTime = El.mpi.Time()
      mode=El.LLL_FULL
      U, R, info = El.LLL(B,mode,ctrl)
      runTime = El.mpi.Time() - startTime
      print "  runtime: %f seconds" % runTime 
      print "  delta=", info.delta
      print "  eta  =", info.eta
      print "  nullity: ", info.nullity
      print "  num swaps: ", info.numSwaps
      if output:
        El.Print( U, "U" );
        El.Print( B, "BNew" )
        El.Print( R, "R" )

      # Test how small the first column is compared to the others
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
      print ""
