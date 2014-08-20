import re, string, sys, os, sys

from protinfo import *
import pdbinfo
from geom import vec_dist

SCWRL = "scwrl"
#SCWRL = "/software/scwrl/scwrl3"

nl = '\n'

def findMissingPatters(patterns, strings) :
	missing_pats = []
	for pat in patterns :
		absent = 1
		for s in strings :
			if re.compile(pat).search(s) :
				absent = 0
				break
		if absent == 1 : missing_pats.append(pat)
	if len(missing_pats) > 0 : return missing_pats
	return None

def isPdbAtomLine(line) :
	return re.compile("^ATOM").search(line) or re.compile("^HETATM").search(line)
def isPdbNUline(line) :
	cats = getResidueCategory(line2resn(line))
	return 'rna' in cats or 'dna' in cats
def isPdbAAline(line) :
	return getResidueCategory(line2resn(line)) == ['aa']
def isPdbHetLine(line) :
	if not isPdbAtomLine(line) : assert(0==1)
	if line[0:6] == "HETATM" : return 1
	if getResidueCategory(line2resn(line)) == ['het'] : return 1
	return None
def isPdbHOHline(line) :
	return line[pdbinfo.resn[0]:pdbinfo.resn[1]] == 'HOH'
def isPdbHydrogenLine(line) :
	an = line2atomname(line)
	return an[0] == 'H' or an[1] == 'H'
def line2crd(line) :
	return ( string.atof(line[pdbinfo.xcrd[0]:pdbinfo.xcrd[1]]), string.atof(line[pdbinfo.ycrd[0]:pdbinfo.ycrd[1]]), string.atof(line[pdbinfo.zcrd[0]:pdbinfo.zcrd[1]]) )
def line2atomname(line) :
	return line[pdbinfo.atmn[0]:pdbinfo.atmn[1]]
def line2atomid(line) :
	return line2resid(line) + line2atomname(line)
def line2resn(line) :
	return line[pdbinfo.resn[0]:pdbinfo.resn[1]]
def line2resn1(line) :
	return AA31[line[pdbinfo.resn[0]:pdbinfo.resn[1]]]
def line2resnum(line) :
	return line[pdbinfo.resi[0]:pdbinfo.resi[1]]
def line2resic(line) :
	return line[pdbinfo.ic[0]:pdbinfo.ic[1]]
def makeResid(resn, chid, resnum, inscode) :
	return "%3s%1s%4s%1s" % (resn, chid, resnum, inscode)
def line2resid(line) :
	return line[pdbinfo.resn[0]:pdbinfo.resn[1]] + line[pdbinfo.chid[0]:pdbinfo.chid[1]] + line[pdbinfo.resi[0]:pdbinfo.resi[1]] + line[pdbinfo.ic[0]:pdbinfo.ic[1]]
def changeResid(line, newl) :
	assert(line2resn(line) == line2resn(newl))
	return line[0:pdbinfo.chid[0]] + newl[pdbinfo.chid[0]:pdbinfo.chid[1]] + line[pdbinfo.chid[1]:pdbinfo.resi[0]] + newl[pdbinfo.resi[0]:pdbinfo.resi[1]] + line[pdbinfo.resi[1]:]
def line2chid(line) :
	return line[pdbinfo.chid[0]:pdbinfo.chid[1]]
def missingMCatoms(reslines) :
	atomnames = []
	for rl in reslines : atomnames.append(line2atomname(rl))
	return findMissingPatters(mc_atomnames, atomnames)
def missingSCatoms(reslines) :
	atomnames = []
	for rl in reslines : atomnames.append(line2atomname(rl))
	return findMissingPatters(sc_atomnames[line2resn(reslines[0])], atomnames)
def sameChain(line1, line2) :
	if line2chid(line1) == line2chid(line2) : return 1
	return None
def sameResidue(pdbline1, pdbline2) :
	if line2resid(pdbline1) == line2resid(pdbline2) : return 1
	return None
def covBonded(pdbline1,pdbline2) :
	if sameResidue(pdbline1,pdbline2) and sameResidueCovbonded(line2atomname(pdbline1), line2atomname(pdbline2), line2resn(pdbline1)) : return 1
	return None
def phosSugarO3P(pdbline1, pdbline2) :
	if not isPdbNUline(pdbline1) or not isPdbNUline(pdbline2) : return None
	resn1,resn2 = line2resn(pdbline1), line2resn(pdbline2)
	if line2atomname(pdbline1) != ' O3*' and not line2atomname(pdbline1) in whichResinfo(resn1)[0].atsyn[' O3*'] : return None
	if line2atomname(pdbline2) != ' P  ' and not line2atomname(pdbline2) in whichResinfo(resn2)[0].atsyn[' P  '] : return None
	return 1
def peptideCN(pdbline1, pdbline2) :
	if not isPdbAAline(pdbline1) and not isPdbAAline(pdbline2) : return None
	resn1,resn2 = line2resn(pdbline1), line2resn(pdbline2)
	if line2atomname(pdbline1) != ' C  ' and not line2atomname(pdbline1) in whichResinfo(resn1)[0].atsyn[' C  '] : return None
	if line2atomname(pdbline2) != ' N  ' and not line2atomname(pdbline2) in whichResinfo(resn2)[0].atsyn[' N  '] : return None
	return 1
## special handling of gromacs O1,O2 instead of simple O, take O1 as O and ignore O2
def gromacsO1O2(l) :
	if isPdbAAline(l) and line2atomname(l) == ' O2 ' : return 'ignore this line please'
	if isPdbAAline(l) and line2atomname(l) == ' O1 ' : return l[0:pdbinfo.atmn[0]] + ' O  ' + l[pdbinfo.atmn[1]:]
	return l
def ignoreOXT(l) :
	if isPdbAAline(l) and line2atomname(l) == ' OXT' : return l[0:pdbinfo.atmn[0]] + ' O  ' + l[pdbinfo.atmn[1]:]
	return l
def change_O_OW_HOH(l) :
	if isPdbHOHline(l) and line2atomname(l)[1] == 'O' : return l[0:pdbinfo.atmn[0]] + ' OW ' + l[pdbinfo.atmn[1]:]
	return l

class protein :
	def __init__(self, filenameORlines, read_hydrogens, read_waters, read_hets, remove_duplicate_water=0, remove_duplicates=1, read_aa_only=0) :
		if 'list' in "%s" % type(filenameORlines) :
			self.name = "read_from_lines"
			print "Reading protein from %d lines" % len(filenameORlines)
			self.readFromLines(filenameORlines, read_hydrogens, read_waters, read_hets, remove_duplicate_water, remove_duplicates, read_aa_only)
		else :
			self.name = filenameORlines
			print "Reading protein from", filenameORlines
			f = open(filenameORlines, 'r')
			if not f : sys.exit(1)
			filelines = f.readlines()
			f.close()
			self.readFromLines(filelines, read_hydrogens, read_waters, read_hets, remove_duplicate_water, remove_duplicates, read_aa_only)

	def readFromLines(self, filelines, read_hydrogens, read_waters, read_hets, remove_duplicate_water, remove_duplicates, read_aa_only) :
		## remove all duplicate atoms in a verbose manner
		atomlines, atomids = [], {}
		for l in filelines :
			l = re.sub(nl, '', l)
			l = ignoreOXT(l)
			l = gromacsO1O2(l)
			l = change_O_OW_HOH(l)
			if re.compile('UNK').search(l) : continue
			if not isPdbAtomLine(l) : continue
			if read_hydrogens == 0 and isPdbHydrogenLine(l) : continue
			if read_waters == 0 and isPdbHOHline(l) : continue
			if read_hets == 0 and isPdbHetLine(l) : continue
			if read_aa_only==1 and not isPdbAAline(l) : continue
			if remove_duplicates == 1 :
				in_atomids = 0
				try :
					if atomids[line2atomid(l)] == 1 : in_atomids = 1
				except KeyError : pass
				if in_atomids == 0 : ##line2atomid(l) not in atomids :
					atomlines.append(l)
					atomids[line2atomid(l)] = 1
				elif isPdbHOHline(l) and remove_duplicate_water == 0: atomlines.append(l)
				else : print "WARNING : detected,ignored duplicate atom :", l
			else : atomlines.append(l)
		self.atomlines = atomlines
		reslines = []
		prevResid = ''
		for i in range(len(atomlines)) :
			resid = line2resid(atomlines[i])
			if i==0 :
				prevResid = resid
				continue
			if i > 0 and resid != prevResid :
				prevResid = resid
				lastResEnd = 0
				if len(reslines) > 0 : lastResEnd = reslines[ len(reslines)-1 ][1]
				reslines.append((lastResEnd, i))
		if len(reslines) == 0 : reslines.append( (0, len(atomlines)) )
		else : reslines.append( (reslines[len(reslines)-1][1], len(atomlines)) )
		self.reslines = reslines

	def atomindex2resindex(self, ai) :
		for ri in range(len(self.reslines)) :
			start, end = self.reslines[ri]
			if start <= ai and ai < end : return ri
		return None
	def allcrds(self) :
		crds = []
		for l in self.atomlines :
			crds.append(line2crd(l))
		return crds
	def findMissingAtoms(self, mode) :
		missing = 0
		seq = ''
		for (start,end) in self.reslines :
			if not isPdbAAline(self.atomlines[start]) : continue
			reslines = []
			for i in range(start,end) : reslines.append(self.atomlines[i])
			if mode == 'mc' : missing_atoms = missingMCatoms(reslines)
			elif mode == 'sc' : missing_atoms = missingSCatoms(reslines)
			else : assert(0)
			if missing_atoms == None :
				seq = seq + string.lower(line2resn1(reslines[0]))
			else :
				print "WARNING : Missing Atoms in", line2resid(reslines[0]), missing_atoms
				missing = 1
				seq = seq + line2resn1(reslines[0])
		if missing == 0 : return None
		return seq
	def info(self) :
		print "--------------------------------------"
		print "PROTEIN derived from :", self.name
		chNres, hets = {}, []
		nAA, nHOH, nHET, nNU = 0,0,0,0
		for (start, end) in self.reslines :
			line = self.atomlines[start]
			if isPdbHOHline(line) : nHOH = nHOH + 1
			elif isPdbHetLine(line) :
				nHET = nHET + 1
				hets.append(line2resn(line))
			elif isPdbAAline(line) : nAA = nAA + 1
			elif isPdbNUline(line) : nNU = nNU + 1
			else : assert(0)
			chain = line2chid(line)
			if chain not in chNres.keys() : chNres[chain] = 0
			chNres[chain] = chNres[chain] + 1
			atomnames = []
			for i in range(start,end) : atomnames.append(line2atomname(self.atomlines[i]))
			print line2resid(line), line2resn(line), getResidueCategory(line2resn(line)), getResidueCategory(line2resn(line), atomnames), 'has', end-start, 'atoms'
		for ch,n in chNres.items() : print "Chain -", ch, "- contains", n, "residues"
		print "nAA, nNU, nHOH, nHET = ", nAA, nNU, nHOH, nHET
		hetstr = "HETS : "
		for h in hets : hetstr = hetstr + " " + h
		print hetstr
		self.findMissingAtoms('mc')
		self.findMissingAtoms('sc')
		print "--------------------------------------"
	def nonHOHsampleLines(self) :
		retlines = []
		for (start, end) in self.reslines :
			line = self.atomlines[start]
			if isPdbHOHline(line) : continue
			retlines.append(line)
		return tuple(retlines)
	def changeNonHOHnumbering(self, reslines) :
		myreslines = self.nonHOHsampleLines()
		allWell = 0
		if len(myreslines) == len(reslines) :
			allWell = 1
			for i in range(len(reslines)) :
				if line2resn(reslines[i]) != line2resn(myreslines[i]) : allWell = 0
		if allWell == 0 :
			print "WARNING : renumbering cannot be done, lengths %d %d" % (len(myreslines), len(reslines))
			return None
		count = -1
		newatomlines = []
		for (start, end) in self.reslines :
			line = self.atomlines[start]
			if isPdbHOHline(line) :
				for i in range(start,end) :
					newatomlines.append(self.atomlines[i])
			else :
				count = count + 1
				assert( line2resn(reslines[count]) == line2resn(myreslines[count]) )
				for i in range(start,end) :
					newatomlines.append( changeResid(self.atomlines[i],reslines[count]) )
		assert( len(self.atomlines) == len(newatomlines) )
		self.atomlines = newatomlines
		return 1
	def writeNonAA(self, filename) :
		f = open(filename, 'w')
		if not f : return
		for l in self.atomlines :
			if isPdbAAline(l) : continue
			f.write(l + nl)
		f.close()
	def write_noSC(self, filename, chain='') : ## write everything except aa sidechains
		f = open(filename, 'w')
		if not f : return
		for l in self.atomlines :
			if len(chain) == 1 and line2chid(l) != chain : continue
			if isPdbAAline(l) and not line2atomname(l) in mc_atomnames : continue
			f.write(l + nl)
		f.close()
	def writeAAchainAndHetHoh(self, filename, chain='') :
		f = open(filename, 'w')
		if not f : return
		for l in self.atomlines :
			if isPdbHydrogenLine(l) : continue
			if len(chain) == 1 and line2chid(l) != chain : continue
			f.write(l + nl)
		f.close()
	def writeAAonly(self, filename, chain='') :
		f = open(filename, 'w')
		if not f : return
		for l in self.atomlines :
			if isPdbHydrogenLine(l) or isPdbHetLine(l) or isPdbHOHline(l) : continue
			if len(chain) == 1 and line2chid(l) != chain : continue
			if isPdbAAline(l) : f.write(l + nl)
		f.close()
	def write(self, filename, nonHOHnonHET=1, het=0, water=1, hydrogen=0) :
		f = open(filename, 'w')
		if not f : return
		for l in self.atomlines :
			if hydrogen==0 and isPdbHydrogenLine(l) : continue
			if nonHOHnonHET==1 and not isPdbHetLine(l) and not isPdbHOHline(l) : f.write(l + nl)
			if het==1 and isPdbHetLine(l) : f.write(l + nl)
			if water==1 and isPdbHOHline(l) : f.write(l + nl)
		f.close()
	def makeSites(self, metasiteDefs) :
		msnames, msframes, msSlaves, mschems, sitecrds, collapseMSLater = [],[],[],[],[],[]
		resns, resids, atomnames, atomids, atomline2metasites = [], [], [], [], []
		for ai in range(len(self.atomlines)) :
			sitecrds.append(line2crd(self.atomlines[ai]))
			resns.append(line2resn(self.atomlines[ai]))
			resids.append(line2resid(self.atomlines[ai]))
			atomnames.append(line2atomname(self.atomlines[ai]))
			atomids.append(line2atomid(self.atomlines[ai]))
			atomline2metasites.append([])
		for ms in metasiteDefs :
			ms.info()
			hasFrame = None
			if ms.frA and ms.frB and ms.frC and len(ms.frA) >= 4 and len(ms.frB) >= 4 and len(ms.frC) >= 4 : hasFrame = 1
			nextresAtoms = []
			for sla in ms.slaveatoms :
				if re.compile("^next").search(sla) : nextresAtoms.append(re.sub("^next", "", sla))
			for ri in range(len(self.reslines)) : ## look for metasites per residue
				start,end, atominds = self.reslines[ri][0],self.reslines[ri][1], []
				if applicableMetasiteResidue(self.atomlines[start],ms) == None : continue
				## other than exact atomname match, matches for following keywords are searched :
				## all_atoms_of_the_residue, next....
				## in case of next, next residue is searched, this is for peptide grps
				## if all reqd atoms are not found, metasite is not defined
				if ms.restype == 'HOH' :
					atominds.append(start) ## make shortcut for time-consuming waters
				elif "all_atoms_of_the_residue" in ms.slaveatoms :
					for k in range(start,end) : atominds.append(k)
				elif "all_sc_atoms" in ms.slaveatoms :
					assert( isPdbAAline(self.atomlines[start]) )
					for k in range(start,end) :
						if not isMCatom(atomnames[k], resns[k]) : atominds.append(k)
				elif "all_mc_atoms" in ms.slaveatoms :
					assert( isPdbAAline(self.atomlines[start]) )
					for k in range(start,end) :
						if isMCatom(atomnames[k], resns[k]) : atominds.append(k)
				elif "all_sugar_atoms" in ms.slaveatoms :
					assert( isPdbNUline(self.atomlines[start]) )
					for k in range(start,end) :
						if isSugarAtom(atomnames[k], resns[k]) : atominds.append(k)
				elif "all_base_atoms" in ms.slaveatoms :
					assert( isPdbNUline(self.atomlines[start]) )
					for k in range(start,end) :
						if isBaseAtom(atomnames[k], resns[k]) : atominds.append(k)
				elif "all_phosphate_atoms" in ms.slaveatoms :
					assert( isPdbNUline(self.atomlines[start]) )
					for k in range(start,end) :
						if isPhosphateAtom(atomnames[k], resns[k]) : atominds.append(k)
				else :
					if len(nextresAtoms) > 0 :
						if end == len(self.atomlines) : continue
						if not sameChain(self.atomlines[start],self.atomlines[end]) : continue
						for k in range(self.reslines[ri+1][0],self.reslines[ri+1][1]) :
							if atomnames[k] in nextresAtoms : atominds.append(k)
					for k in range(start,end) :
						if atomnames[k] in ms.slaveatoms : atominds.append(k)
				if len(atominds) > 0 :
					msnames.append("[" + resids[start] + "]_" + ms.name)
					if ms.style == 'atomic' : collapseMSLater.append(0)
					else : collapseMSLater.append(1)
					msSlaves.append(atominds)
					mschems.append(ms.physchem)
					if not hasFrame : msframes.append(None)
					else :
						a,b,c = None, None, None
						for k in range(start, end) :
							if atomnames[k] == ms.frA : a = sitecrds[k]
							if atomnames[k] == ms.frB : b = sitecrds[k]
							if atomnames[k] == ms.frC : c = sitecrds[k]
						assert a
						assert b
						assert c
						assert a != b
						assert c != b
#(a and b and c and (not a==b) and (not b==c))
						msframes.append((tuple(a),tuple(b),tuple(c)))
					#print "making metasite", ms.name, "for", resids[start], "with :"
					for i in atominds :
						atomline2metasites[i].append(len(msnames)-1)
						#print line2atomid(self.atomlines[i])
		for li in range(len(atomline2metasites)) : ## each atom is in at most one metasite
			if len(atomline2metasites[li]) == 0 :
				print "WARNING :", line2atomid(self.atomlines[li]), "has no metasite"
			elif len(atomline2metasites[li]) > 1 :
				print "Multiple metasites, fatal !! :", li, msnames[atomline2metasites[li][0]], msnames[atomline2metasites[li][1]]
				assert(0)
		newscrds, newMSslaves, newSlavenames = [],[],[]
		for msi in range(len(msnames)) :	
			newMSslaves.append([])
			newSlavenames.append([])
			for si in msSlaves[msi] :
				newSlavenames[msi].append(self.atomlines[si])
			if collapseMSLater[msi] == 0 :
				for si in msSlaves[msi] :
					newMSslaves[msi].append(len(newscrds))
					newscrds.append(sitecrds[si])
			else :
				meancrd = [.0,.0,.0]
				for si in msSlaves[msi] :
					meancrd[0] = meancrd[0] + sitecrds[si][0]
					meancrd[1] = meancrd[1] + sitecrds[si][1]
					meancrd[2] = meancrd[2] + sitecrds[si][2]
				if len(msSlaves[msi]) > 1 :
					meancrd[0] = meancrd[0] / len(msSlaves[msi])
					meancrd[1] = meancrd[1] / len(msSlaves[msi])
					meancrd[2] = meancrd[2] / len(msSlaves[msi])
				newMSslaves[msi].append(len(newscrds))
				newscrds.append(tuple(meancrd))
		assert len(msframes) == len(msnames)
		return msnames, newMSslaves, mschems, newscrds, newSlavenames, msframes

## checks if metasite residue is same as residue name, if not looks for keyword matches
## keywords are AA, HOH, all_nonHOH_nonAA, any
def applicableMetasiteResidue(atomline, ms) :
	if ms.restype == "any" : return 1
	if ms.restype == line2resn(atomline) : return 1
	if ms.restype == "AA" and isPdbAAline(atomline) : return 1
	if ms.restype == "NU" and isPdbNUline(atomline) : return 1
	if ms.restype == "HOH" and isPdbHOHline(atomline) : return 1
	if ms.restype == "all_nonHOH_nonAA" and (not isPdbHOHline(atomline)) and (not isPdbAAline(atomline)) : return 1
	if ms.restype == "all_nonHOH_nonAA_nonNU" and (not isPdbHOHline(atomline)) and (not isPdbAAline(atomline)) and (not isPdbNUline(atomline)) : return 1
	return None

## check whether consecutive CA atoms in same chain are less than 5A apart from e.o., return None if farther
def checkMCcont(prot) :
	for ri in range(len(prot.reslines)-1) :
		s1,e1 = prot.reslines[ri]
		s2,e2 = prot.reslines[ri+1]
		if not isPdbAAline(prot.atomlines[s1]) or not isPdbAAline(prot.atomlines[s2]) : continue
		if line2chid(prot.atomlines[s1]) != line2chid(prot.atomlines[s2]) : continue
		ca1, ca2 = -1, -1
		for i in range(s1,e1) :
			if line2atomname(prot.atomlines[i]) == ' CA ' :
				ca1 = i
				break
		for i in range(s2,e2) :
			if line2atomname(prot.atomlines[i]) == ' CA ' :
				ca2 = i
				break
		if ca1 == -1 or ca2 == -1 : continue ## missing CA atoms is none of this function's business
		caD = vec_dist( line2crd(prot.atomlines[ca1]), line2crd(prot.atomlines[ca2]) )
		if caD > 5. :
			print "-%s-%s-" % (line2atomid(prot.atomlines[ca1]), line2atomid(prot.atomlines[ca2])), "are farther than expected for consecutive CA atoms", caD
			return None
	return 1

# abort is there are missing mainchain atoms
# run scwrl if there are missing sidechain atoms
# write output file as pdbfile_out
# hydrogens are ignored
def checkRepairMissingSCatoms(prot, pdbfile_out, tempdir='/tmp/') :
	seq = prot.findMissingAtoms('sc')
	if not seq :
		prot.write(pdbfile_out, hydrogen=0, nonHOHnonHET=1, water=1, het=1)
		#proc_run_exitOnError("cp %s %s" % (pdbfile, pdbfile_out))
		return 1
	## we have missing sidechains
	scwrl_basefile = tempdir + '/' + "scwrl"
	#scwrl_basefile = tempdir + '/' + re.sub("\..pdb", "", os.path.basename(pdbfile))
	scwrl_in = "%s.in" % scwrl_basefile
	scwrl_out = "%s.out" % scwrl_basefile
	scwrl_frame = "%s.frame" % scwrl_basefile
	scwrl_seq = "%s.seq" % scwrl_basefile
	prot.writeAAonly(scwrl_in)
	prot.writeNonAA(scwrl_frame)
	f = open(scwrl_seq, 'w')
	f.write(seq + '\n')
	f.close()
	proc_run_exitOnError("%s -i %s -o %s -s %s -f %s" % (SCWRL, scwrl_in, scwrl_out, scwrl_seq, scwrl_frame),[],[256,35840])
	proc_run_exitOnError("cp %s %s" % (scwrl_out, pdbfile_out))
	proc_run_exitOnError("cat %s >> %s" % (scwrl_frame, pdbfile_out))
	#proc_run_exitOnError("rm -Rf %s %s %s %s" % (scwrl_in, scwrl_out, scwrl_seq, scwrl_frame))

if __name__ == "__main__" :
	checkRepairMissingAtoms(sys.argv[1], sys.argv[2])
