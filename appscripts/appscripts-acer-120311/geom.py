#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################


import math, re, string, sys

def cross_product(v1, v2):
	x1 = v1[0]
	y1 = v1[1]
	z1 = v1[2]
	x2 = v2[0]
	y2 = v2[1]
	z2 = v2[2]
	cp = [0,0,0]
	cp[0] = y1*z2 - y2*z1
	cp[1] = (x1*z2 - x2*z1) * (-1)
	cp[2] = x1*y2 - x2*y1
	return cp

def normalize(v):
	mag = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
	if mag == 0 : mag = 1
	return [v[0]/mag, v[1]/mag, v[2]/mag]

def dot_product(v1, v2) :
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def vec_angle_deg(x,y) :
	#print x, v1, v2
	#print y, v3, v2
	x,y = normalize(x), normalize(y)
	dp = dot_product(x,y)
	if dp > 0.9999 : dp = 0.9999
	if dp < -0.9999 : dp = -0.9999
	return 180./math.pi * math.acos(dp)

def vec_angle_deg3(v1, v2, v3):
	x = vec_diff(v1,v2)
	y = vec_diff(v3,v2)
	return vec_angle_deg(x,y)

def vec_dihed_deg(v1, v2, v3, v4):
	## check angle v2-v3-v4 first. if colinear, return 0
	ang = vec_angle_deg3(v2,v3,v4)
	if math.fabs(180-ang) < 1 or math.fabs(ang) < 1 : return 0
	ang = vec_angle_deg3(v1,v2,v3)
	if math.fabs(180-ang) < 1 or math.fabs(ang) < 1 : return 0
	## now find actual dihedral
	n21 = normalize(vec_diff(v1,v2))
	n23 = normalize(vec_diff(v3,v2))
	n2 = normalize(cross_product(n21, n23))
	n32 = normalize(vec_diff(v2,v3))
	n34 = normalize(vec_diff(v4,v3))
	n3 = normalize(cross_product(n32, n34))
	dp = dot_product(n2,n3)
	if dp > 0.9999 : dp = 0.9999
	if dp < -0.9999 : dp = -0.9999
	ang = 180./math.pi * math.acos(dp)
	if dot_product(cross_product(n2,n3), vec_diff(v3, v2)) < 0. : ang = -1. * ang
	return ang

def vec_diff(v1, v2):
	return [ v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2] ]

def vec_dist(v1, v2):
	vd = vec_diff(v1,v2)
	return math.sqrt( vd[0]*vd[0] + vd[1]*vd[1] + vd[2]*vd[2] )

def vec_add(v1, v2):
	return [ v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2] ]

def vec_scale(v, scale):
	return [ v[0]*scale, v[1]*scale, v[2]*scale ]


#def poly_vol_sa(vert):
#	cmd = "qconvex FA"
#	lines_in = []
#	lines_in.append("3")
#	lines_in.append("%d" % len(vert))
#	for v in vert : lines_in.append("%f %f %f" % (v[0], v[1], v[2]))
#	#print lines_in
#    import procrun
#	lines_out = procrun.proc_run_exitOnError(cmd, lines_in)[1]
#	v, a = 0.,0.
#	for l in lines_out :
#		#print l
#		if re.compile(" volume:").search(l) :
#			v = string.atof(re.sub(".* volume:", "", l))
#		if re.compile(" area:").search(l) :
#			a = string.atof(re.sub(".* area:", "", l))
#	return v, a
#	assert(0)

def triangle_area(v) :
	assert( len(v) == 3 )
	a = [0.,0.,0.]
	for p in range(3) :
		q = p + 1
		if q == 3 : q = 0
		for i in range(3) :
			d = v[p][i] - v[q][i]
			a[p] = a[p] + d * d
		a[p] = math.sqrt(a[p])
	s = (a[0] + a[1] + a[2]) / 2.
	return math.sqrt(s * (s-a[0]) * (s-a[1]) * (s-a[2]))

# assuming convex ordered polygon and adding triangle areas
def poly_area(v):
	assert(len(v) > 2)
	area = 0.
	for i in range(1,len(v)-1) :
		area = area + triangle_area((v[0],v[i],v[i+1]))
	return area

def poly_area_old(vert):
	cp = cross_product( vec_diff(vert[1], vert[2]), vec_diff(vert[0], vert[2]) )
	cp = normalize(cp)

	cmd = "qconvex FA"
	lines_in = []
	lines_in.append("3")
	lines_in.append("%d" % (len(vert)*2))
	for v in vert :
		lines_in.append("%f %f %f" % (v[0], v[1], v[2]))
		lines_in.append("%f %f %f" % (v[0]+cp[0], v[1]+cp[1], v[2]+cp[2]))
	lines_out = proc_run(cmd, lines_in)
	for l in lines_out :
		#print l
		if re.compile(" volume:").search(l) :
			return string.atof(re.sub(".* volume:", "", l))
	assert(0)

def find4thPoint(a, b, c, dist, ang, dihed) :
	ang1 = 180 - ang
	cosA, sinA = math.cos( math.pi * ang1/180. ), math.sin( math.pi * ang1/180. )
	bcNorm = normalize( vec_diff(c,b) )
	proj = -1 * dot_product( vec_diff(a,b), bcNorm )
	yAxis = vec_add ( vec_scale(bcNorm,proj) , vec_diff(a,b) )
	#print bcNorm, yAxis, dot_product(bcNorm,yAxis)
	yAxis = normalize ( yAxis )
	#print bcNorm, yAxis, dot_product(bcNorm,yAxis)
	assert math.fabs(dot_product(bcNorm,yAxis)) < 1e-3
	pt = vec_add( vec_scale(bcNorm, dist*cosA), vec_scale(yAxis, dist*sinA) )
	pt1 = rotate(pt, bcNorm, dihed)
	pt = vec_add(c, pt1)
	return pt
	#print "CHECK", dist, ang, dihed, vec_dist(c,pt), vec_angle_deg3(b,c,pt), vec_dihed_deg(a,b,c,pt)
	#retpt = rotate(pt, bcNorm, dihed)
	print "CHECK", dist, ang, dihed, vec_dist(c,retpt), vec_angle_deg3(b,c,retpt), vec_dihed_deg(a,b,c,retpt)
	return retpt


## reflect pt abt plane given by points pl1, pl2, pl3
## first find the projection of pt on plane
def reflect(pt, pl1, pl2, pl3) :
	d = vec_dist(pl3, pt)
	an = vec_angle_deg3(pl2, pl3, pt)
	dh = vec_dihed_deg(pl1, pl2, pl3, pt)
	return find4thPoint(pl1, pl2, pl3, d, an, -1*dh)

## rotate abt axis A through D degrees
def rotate(pt, axis, angle_deg) :
	angle = math.pi/180. * angle_deg
	c = math.cos(angle)
	s = math.sin(angle)
	R = [ [0.,0.,0.], [0.,0.,0.], [0.,0.,0.], ]
	axis = normalize(axis)
	x,y,z = axis[0], axis[1], axis[2]
	t = 1. - c
	R[0][0],R[0][1],R[0][2] = t*x*x + c,    t*x*y - s*z,   t*x*z + s*y
	R[1][0],R[1][1],R[1][2] = t*x*y + s*z,  t*y*y + c,     t*y*z - s*x
	R[2][0],R[2][1],R[2][2] = t*x*z - s*y,  t*y*z + s*x,   t*z*z + c
	ret = [0.,0.,0.]
	ret[0] = R[0][0]*pt[0] + R[0][1]*pt[1] + R[0][2]*pt[2]
	ret[1] = R[1][0]*pt[0] + R[1][1]*pt[1] + R[1][2]*pt[2]
	ret[2] = R[2][0]*pt[0] + R[2][1]*pt[1] + R[2][2]*pt[2]
	return ret

def dihedDiff(ang1,ang2) : ## absolute difference between dihedral angles in degrees
	assert -180 <= ang1 and ang1 <= 180
	assert -180 <= ang2 and ang2 <= 180
	diff = math.fabs(ang1 - ang2)
	if diff > 180 : diff = 360-diff
	return diff

if __name__=="__main__" :
	axis = (2,1,3)
	pt = (0,2,0)
	print pt, rotate(pt, axis, 90)
	print pt, rotate(pt, axis, -90)
	print pt, rotate(pt, axis, 0)
	print pt, rotate(pt, axis, 180)
	sys.exit(0)
	print dihedDiff(100,90)
	print dihedDiff(100,-90)
	print dihedDiff(-100,90)
	print dihedDiff(0,0)
	print reflect([1,1,1], [1,0,0], [0,1,0], [0,0,0])
