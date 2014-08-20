from geom import VecFloat, calcDihed
import string

n, ca, c, cb = None, None, None, None
for l in open("/home/swanand/downloads/1A1F.pdb", 'r').readlines():
#ATOM    577  N   ALA A 173      10.105  15.455  36.454  1.00 28.85      DSNA N
#012345678901234567890123456789012345678901234567890123456789
#000000000011111111112222222222333333333344444444445555555555
    if l[0:4] != 'ATOM' : continue
    crd = [
	string.atof(l[30:38]),
	string.atof(l[38:46]),
	string.atof(l[46:54]),
    ]
    if l[12:16] == " N  " :
	n, ca, c, cb = crd, None, None, None
    elif l[12:16] == " CA " : ca = crd
    elif l[12:16] == " C  " : c = crd
    elif l[12:16] == " CB " :
	cb = crd
	#print n, ca, c, cb
	if n and ca and c and cb :
	    print l[17:20], calcDihed(
		VecFloat(n), VecFloat(c),
		VecFloat(ca), VecFloat(cb),
	    )
	n, ca, c, cb = None, None, None, None
