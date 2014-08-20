	## check constaints
	#constraintsOK == 0
	#for ci in b2c[index] :
	#    constraintsOK = constraints[ci].satisfied()
	#    if constraintsOK == 0 : break
	#if constraintsOK == 0 : continue
#import os, sys
#print os.environ["PYTHONPATH"]
#print os.environ["LD_LIBRARY_PATH"]
#print sys.argv
#sys.exit(0)

from geom import Grid
from builders import Builder, VecInt, VecFloat, VecVecFloat, SimpleBuilder

pts = VecVecFloat()
radii = VecFloat()
clashExclInds, builders, constraints = [], [], []
b2c = {} #maps bldr index to indices of constraints to be tested after that builder

## a build strategy, another could be divide-and-rule
def build(grid, index) :
    cb, cbip, cbop = builders[index], builders[index].getIP(), builders[index].getOP()
    for nt in range(50) :
	print "BUILDER index %d trial %d" % (index, nt)
	## just build
	cb.build()
	## clashcheck using grid
	if g.add(cbop, clashExclInds[index]) == 0 : continue
	## call next builder
	if len(builders) == index+1 : return 1
	nextBldrOp = build(grid, index+1)
	if nextBldrOp : return 1
	## remove added points from grid
	for ai in list(cbop): grid.remove(ai)
    return 0

if __name__ == "__main__" :
    ## prepare all points and radii, init starting points ##############
    pts.push_back( VecFloat([0,2,0]) )
    pts.push_back( VecFloat([0,0,0]) )
    pts.push_back( VecFloat([2,0,0]) )
    pts.push_back( VecFloat([-999,-999,-999]) )
    pts.push_back( VecFloat([-999,-999,-999]) )
    for i in range(pts.size()) : radii.push_back(1)
    ####################################################################

    ## prepare builders #################################################
    builders.append( SimpleBuilder(pts, VecInt([0,1,2]), VecInt([3]), 2, 135, 0) )
    clashExclInds.append( [2] )
    builders.append( SimpleBuilder(pts, VecInt([1,2,3]), VecInt([4]), 2, 135, 0) )
    clashExclInds.append( [3] )
    for ci in range(len(clashExclInds)) : clashExclInds[ci] = VecInt( clashExclInds[ci] )
    #####################################################################

    ## start model building
    nmodels = 10000
    for nm in range(nmodels) :
	print "starting MODEL ", nm
	## add starting points to grid unconditionally
	g = Grid(1, pts, radii)
	g.justAdd(builders[0].getIP())
	## start building
	if build(g, 0) :
	    print "MODEL %d succeeded" % nm
	    for i in range(pts.size()) : print list(pts[i])
