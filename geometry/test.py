from geom import Grid, FloatVector
from random import random

def main() :
    grid = Grid(1)
    points = []
    for i in range(1000) :
	ri = -1*random() * 5
	rj = random() * 5
	rk = random() * 5
	print "Start Addition"
	id = grid.add(FloatVector([ri,rj,rk]), 1)
	print "id", id
	if id >= 0 : points.append( [ri,rj,rk] )
	else :
	    clash = 0
	    for pi in range(len(points)) :
		P, Q = points[pi], [ri,rj,rk]
		dist = (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2])
		if dist < 4 : clash = 1
	    assert clash == 1
    for pi in range(len(points)) :
	for pj in range(pi+1,len(points)) :
	    P, Q = points[pi], points[pj]
	    dist = (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2])
	    assert dist > 4

if __name__ == "__main__" :
    main()
