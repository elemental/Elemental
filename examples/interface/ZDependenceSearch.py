import El

# NOTE:
# This is *not* yet functioning for modest-sized n!
# There appears to be coefficient
# explosion and I am unsure of whether or not it is inherent to LLL
# or due to me not using high-enough precision. This will soon be investigated.

n=75
tag=El.zTag
progress=True
timeLLL=False
output=True

# Initially draw z out of a uniform ball of radius five about 10+0i
z = El.Matrix(tag)
El.Uniform(z,n,1,10.,5.)

# Compute a (hidden) Gaussian integer vector aHidden
aHidden = El.Matrix(tag)
El.Uniform(aHidden,n-1,1,0.,5.)
El.Round(aHidden)

# Force the last entry of z to be the above integer linear combination of the
# first n-1 entries
zetaLast=El.Dotu(aHidden,z[0:n-1,:])
z.Set(n-1,0,zetaLast)

if output:
  El.Print( aHidden, "aHidden" )
  El.Print( z, "z" )

ctrl=El.LLLCtrl_d()
ctrl.progress = progress
ctrl.time = timeLLL

NSqrt=10000.

B=El.Matrix()
for presort, smallestFirst in (True,True), (True,False), (False,False):
  for deltaLower in 0.5, 0.75, 0.95, 0.98:
    for weak in True, False:

      print "weak=%r, presort=%r, smallestFirst=%r, deltaLower=%f" % \
        (weak,presort,smallestFirst,deltaLower)

      ctrl.delta = deltaLower
      ctrl.weak = weak 
      ctrl.presort = presort
      ctrl.smallestFirst = smallestFirst

      # Search for the linear dependence
      startTime = El.mpi.Time()
      numFound, B, U = El.ZDependenceSearch(z,NSqrt,ctrl)
      runTime = El.mpi.Time() - startTime
      print "  runtime: %f seconds" % runTime
      print "  num found: ", numFound
      if output:
        El.Print( B, "B" )
        El.Print( U, "U" )
