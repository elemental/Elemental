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
    for variant in El.LLL_NORMAL, El.LLL_DEEP, El.LLL_DEEP_REDUCE:

      print('variant={}, presort={}, smallest1st={}, deltaLower={}'.format( \
        variant,presort,smallestFirst,deltaLower))

      El.Copy( BOrig, B )

      ctrl.delta = deltaLower
      ctrl.variant = variant
      ctrl.presort = presort
      ctrl.smallestFirst = smallestFirst

      # Compute the image and kernel
      startTime = El.mpi.Time()
      M, K = El.LatticeImageAndKernel(B,ctrl)
      runTime = El.mpi.Time() - startTime
      print('  runtime: {} seconds'.format(runTime))
      print('  nullity: {}'.format(K.Width()))
      if output:
        El.Print( M, "Image" )
        El.Print( K, "Kernel" )

El.Finalize()
