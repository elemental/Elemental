import El

n=4
progress=False
timeLLL=False
output=True

# Form 
# -6   9 -15 -18 
#  4  -6  10  12 
# 10 -15  18  35 
#-24  36 -46 -82 
BOrig = El.Matrix()
El.Zeros(BOrig,n,n)
BOrig.Set(0,0,-6)
BOrig.Set(1,0,+4)
BOrig.Set(2,0,10)
BOrig.Set(3,0,-24)
BOrig.Set(0,1,9)
BOrig.Set(1,1,-6)
BOrig.Set(2,1,15)
BOrig.Set(3,1,36)
BOrig.Set(2,1,-15)
BOrig.Set(0,2,-15)
BOrig.Set(1,2,10)
BOrig.Set(2,2,18)
BOrig.Set(3,2,-46)
BOrig.Set(0,3,-18)
BOrig.Set(1,3,12)
BOrig.Set(2,3,35)
BOrig.Set(3,3,-82)
if output:
  El.Print( BOrig, "B" )

ctrl=El.LLLCtrl_d()
ctrl.progress = progress
ctrl.time = timeLLL

B=El.Matrix()
for presort, smallestFirst in (True,True), (True,False), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    for variant in El.LLL_WEAK, El.LLL_NORMAL, El.LLL_DEEP, El.LLL_DEEP_REDUCE:

      print "variant=%d, presort=%r, smallestFirst=%r, deltaLower=%f" %\
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
        El.Print( U, "U" )
        El.Print( B, "BNew" )
        El.Print( R, "R" )

      # Test how small the first column is compared to the others
      b1Norm = El.FrobeniusNorm(B[:,0])
      print "  || b_1 ||_2 = ",b1Norm
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
