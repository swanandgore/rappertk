import os, sys, popen2

if "RTKROOT" not in os.environ.keys() or len(os.environ["RTKROOT"]) == 0 :
    print "RTKROOT is not defined, define and try again."
    sys.exit(1)

modulesInUse = ['scep', 'builders', 'restraints', 'samplers', 'geometry', 'misc',]

envvars = ["LD_LIBRARY_PATH", "PYTHONPATH"]
for v in envvars :
    if v not in os.environ : os.environ[v] = ''
    for m in modulesInUse :
	os.environ[v] = os.environ[v] + ":" + os.environ["RTKROOT"] + os.sep + m
    os.environ[v] = os.environ[v] + ":" + os.environ["RTKROOT"]
 
#print os.environ["LD_LIBRARY_PATH"]
#print os.environ["PYTHONPATH"]

args = ['rappertk'] + [os.environ["RTKROOT"] + os.sep + "appscripts-new" + os.sep + sys.argv[1] + ".py"] + sys.argv[2:] #padding 0th arg

#print
#print

#print "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
#print "  Welcome to RapperTK - a versatile engine for conformational sampling of macromolecules   "
#print "        Developed by Swanand Gore and Anjum Karmali at Dept of Biochemistry, University of Cambridge         "
#print "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"


#print args
#sys.exit(0)
os.execvpe("python", args, os.environ)

#c = popen2.Popen3(cmd)
#c.wait()
#for l in c.fromchild.readlines() : print l
