#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################

import popen2, re, sys

def proc_run(cmd, input_lines=[], acc_exit_status = []) :
	c = popen2.Popen3(cmd, capturestderr=0)
	print "EXECUTING", cmd
	for line in input_lines :
		c.tochild.write(line + "\n")
	c.tochild.close()
	out_lines, err_lines = [], []
	for l in c.fromchild.readlines() :
		out_lines.append(re.sub("\n", "", l))
	c.fromchild.close()
	exitStatus = c.wait()
	#c.childerr.close()
#	print "EXIT STAT", exitStatus, acc_exit_status
	if exitStatus != 0 and exitStatus not in acc_exit_status :
		print "****exitStatus of -", cmd, "- was nonzero****", exitStatus
		for l in out_lines : print "STDOUT", re.sub('\n', '', l)
		for l in err_lines : print "STDERR", re.sub('\n', '', l)
		##print "Is the given structure factor file an MTZ file ? If not, please give a MTZ format file!"
		#import sys ; sys.exit()
	return exitStatus, out_lines, err_lines

def proc_run_exitOnError(cmd, input_lines=[], acc_exit_status=[]) :
	exitStatus, ol, el = proc_run(cmd, input_lines, acc_exit_status)
	if exitStatus != 0 and exitStatus not in acc_exit_status : sys.exit(exitStatus)
	return exitStatus, ol, el

if __name__=="__main__" :
	exitStatus, olines, elines = proc_run("ls /tmp/devices/devils", [])
	for l in olines : print "STDOUT", re.sub('\n', '', l)
	exitStatus, olines, elines = proc_run("ls", [])
	for l in olines : print "STDOUT", re.sub('\n', '', l)

def fileOpen(filename, mode) :
	f = None
	try : f = open(filename, mode)
	except IOError : print "Cannot open %s in %s mode" % (filename, mode)
	if f : return f
	return None
