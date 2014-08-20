#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################

import string, re

class Resinfo :
	def __init__(self, name3, name1, atomnames, covconn, type, atomname_synonyms, name3_synonyms) :
		assert( type == 'aa' or type == 'het' or type == 'rna' or type == 'dna')
		assert( len(atomnames) > 0 )
		assert( len(covconn) > 0 )
		assert( len(name1) == 1 )
		assert( len(name3) == 3 )
		self.name3 = name3
		self.covconn = covconn
		self.name1 = name1
		#print self.name3, self.covconn
		self.atomnames = atomnames
		self.atsyn = {}
		for atn in atomnames :
			if atn in atomname_synonyms.keys() : self.atsyn[atn] = atomname_synonyms[atn]
			else : self.atsyn[atn] = []
		for i in range(len(self.atomnames)) : ## check all atomnames are distinct
			for k in range(i+1,len(self.atomnames)) :
				assert (self.atomnames[i] != self.atomnames[k])
		for (a,b) in covconn : ## check covalent connctivity is good
			assert(a in self.atomnames)
			assert(b in self.atomnames)
		for an in self.atsyn.keys() : assert(an in self.atomnames)
		self.name3syn = name3_synonyms
		self.type = type
		for atn1 in self.atomnames : ## check whether atomnames and their synonyms dont clash
			atoms1 = []
			if atn1 in self.atsyn.keys() : atoms1 = [atn1] + self.atsyn[atn1]
			else : atoms1 = [atn1]
			for atn2 in self.atomnames :
				if atn1 == atn2 : continue
				if atn2 in self.atsyn.keys() : atoms2 = [atn2] + self.atsyn[atn2]
				else : atoms2 = [atn2]
				for at in atoms1 : assert not at in atoms2
	def containsAtom(self, atomname) :
		for at in self.atomnames :
			if atomname == at or atomname in self.atsyn[at] : return 1
		return None
	def consistsOfAtoms(self, atomnames) :
		if len(self.atomnames) != len(atomnames) : return None
		for atn in atomnames :
			if not self.containsAtom(atn) :
				print self.name3, 'does not contain', atn
				return None
		return 1
	def bondedAtoms(self, atn1, atn2) :
		a1, a2 = '', ''
		for at in self.atomnames :
			if atn1 == at or atn1 in self.atsyn[at] : a1 = at
			if atn2 == at or atn2 in self.atsyn[at] : a2 = at
		if a1 == '' or a2 == '' : return None
		if (a1,a2) in self.covconn or (a2,a1) in self.covconn : return 1
		return None

knownResinfo = []

AA20 = [
"PHE", "TYR", "TRP", "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "ASP",
"GLU", "HIS", "LYS", "ARG", "CYS", "MET", "MSE", "SER", "THR", "ASN", "GLN", ]

	
AA31 = {
"PHE":'F',
"TYR":'Y',
"TRP":'W',
"GLY":'G',
"ALA":'A',
"VAL":'V',
"LEU":'L',
"ILE":'I',
"PRO":'P',
"ASP":'D',
"GLU":'E',
"HIS":'H',
"LYS":'K',
"ARG":'R',
"CYS":'C',
"MET":'M',
"MSE":'O',
"SER":'S',
"THR":'T',
"ASN":'N',
"GLN":'Q',
}

AA13 = {}
for k,v in AA31.items() : AA13[v] = k

mc_atomnames = [' N  ', ' CA ', ' C  ', ' O  ']

sc_atomnames = {
"PHE":[' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ '],
"TYR":[' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' OH '],
"TRP":[' CB ', ' CG ', ' CD1', ' CD2', ' NE1', ' CE2', ' CE3', ' CZ2', ' CZ3', ' CH2'],
"GLY":[],
"ALA":[' CB '],
"VAL":[' CB ', ' CG1', ' CG2'],
"LEU":[' CB ', ' CG ', ' CD1', ' CD2'],
"ILE":[' CB ', ' CG1', ' CG2', ' CD1'],
"PRO":[' CB ', ' CG ', ' CD '],
"ASP":[' CB ', ' CG ', ' OD1', ' OD2'],
"GLU":[' CB ', ' CG ', ' CD ', ' OE1', ' OE2'],
"HIS":[' CB ', ' CG ', ' ND1', ' CD2', ' CE1', ' NE2'],
"LYS":[' CB ', ' CG ', ' CD ', ' CE ', ' NZ '],
"ARG":[' CB ', ' CG ', ' CD ', ' NE ', ' CZ ', ' NH1', ' NH2'],
"CYS":[' CB ', ' SG '],
"MET":[' CB ', ' CG ', ' SD ', ' CE '],
"MSE":[' CB ', ' CG ', 'SE  ', ' CE '],
"SER":[' CB ', ' OG '],
"THR":[' CB ', ' OG1', ' CG2'],
"ASN":[' CB ', ' CG ', ' OD1', ' ND2'],
"GLN":[' CB ', ' CG ', ' CD ', ' OE1', ' NE2'],
}

### residue-wise list of covalent connectivity
resAAconn = {}

resAAconn['GLY'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
]
resAAconn['ALA'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
]
resAAconn['VAL'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG1'),
( ' CB ', ' CG2'),
]
resAAconn['LEU'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD1'),
( ' CG ', ' CD2'),
]
resAAconn['ILE'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG1'),
( ' CB ', ' CG2'),
( ' CG1', ' CD1'),
]
resAAconn['PRO'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD '),
( ' CD ', ' N  '),
]
resAAconn['PHE'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD1'),
( ' CG ', ' CD2'),
( ' CD1', ' CE1'),
( ' CD2', ' CE2'),
( ' CE1', ' CZ '),
( ' CE2', ' CZ '),
]
resAAconn['TYR'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD1'),
( ' CG ', ' CD2'),
( ' CD1', ' CE1'),
( ' CD2', ' CE2'),
( ' CE1', ' CZ '),
( ' CE2', ' CZ '),
( ' CZ ', ' OH '),
]
resAAconn['TRP'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD1'),
( ' CG ', ' CD2'),
( ' CD1', ' NE1'),
( ' CD2', ' CE2'),
( ' CD2', ' CE3'),
( ' NE1', ' CE2'),
( ' CE2', ' CZ2'),
( ' CE3', ' CZ3'),
( ' CZ2', ' CH2'),
( ' CZ3', ' CH2'),
]
resAAconn['ASP'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' OD1'),
( ' CG ', ' OD2'),
]
resAAconn['GLU'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD '),
( ' CD ', ' OE1'),
( ' CD ', ' OE2'),
]
resAAconn['HIS'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' ND1'),
( ' CG ', ' CD2'),
( ' ND1', ' CE1'),
( ' CD2', ' NE2'),
( ' CE1', ' NE2'),
]
resAAconn['LYS'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD '),
( ' CD ', ' CE '),
( ' CE ', ' NZ '),
]
resAAconn['ARG'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD '),
( ' CD ', ' NE '),
( ' NE ', ' CZ '),
( ' CZ ', ' NH1'),
( ' CZ ', ' NH2'),
]
resAAconn['CYS'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' SG '),
]
resAAconn['MET'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' SD '),
( ' SD ', ' CE '),
]
resAAconn['MSE'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', 'SE  '),
( 'SE  ', ' CE '),
]
resAAconn['SER'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' OG '),
]
resAAconn['THR'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' OG1'),
( ' CB ', ' CG2'),
]
resAAconn['ASN'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' OD1'),
( ' CG ', ' ND2'),
]
resAAconn['GLN'] = [
( ' N  ', ' CA '),
( ' CA ', ' C  '),
( ' C  ', ' O  '),
( ' CA ', ' CB '),
( ' CB ', ' CG '),
( ' CG ', ' CD '),
( ' CD ', ' OE1'),
( ' CD ', ' NE2'),
]

O_synonyms = {' O  ' : [' O1 ', ' O2 ', ' OXT',]}
for aa in AA20 :
	#print 'OSYN :', O_synonyms
	knownResinfo.append(Resinfo(aa, AA31[aa], mc_atomnames + sc_atomnames[aa], resAAconn[aa], 'aa', O_synonyms, []))

NU5 = ['  A', '  T', '  C', '  G', '  U']
nu_atomnames, nu_covconn = {}, {}
nu_synonyms = {
' OP1':[' O1P', ' OA '],
' OP2':[' O2P', ' OB '],
}
phosphate_atomnames = [' OP1', ' OP2', ' P  ']
phosphate_covconn = [(' P  ',' OP1'), (' P  ',' OP2')]
deoxyribose_sugar_atomnames = [ ' C1*', ' C2*', ' C3*', ' O3*', ' C4*', ' O4*', ' C5*', ' O5*',]
deoxyribose_sugar_covconn = [
(' O5*', ' C5*'),
(' C4*', ' C5*'),
(' C4*', ' O4*'),
(' C1*', ' O4*'),
(' C4*', ' C3*'),
(' O3*', ' C3*'),
(' C2*', ' C3*'),
(' C2*', ' C1*'),
]
ribose_sugar_atomnames = deoxyribose_sugar_atomnames + [' O2*']
ribose_sugar_covconn = deoxyribose_sugar_covconn + [(' C2*', ' O2*')]
AG_atomnames = [ ' N1 ', ' C2 ', ' N3 ', ' C4 ', ' C5 ', ' C6 ', ' N7 ', ' C8 ', ' N9 ', ]
nu_atomnames['  A'] = AG_atomnames + [' N6 ']
nu_atomnames['  G'] = AG_atomnames + [' O6 ', ' N2 ']
AG_covconn = [
(' N1 ', ' C2 '),
(' N3 ', ' C2 '),
(' N3 ', ' C4 '),
(' C5 ', ' C4 '),
(' C5 ', ' C6 '),
(' N1 ', ' C6 '),
(' C5 ', ' N7 '),
(' C8 ', ' N7 '),
(' C8 ', ' N9 '),
(' C4 ', ' N9 '),
]
nu_covconn['  A'] = AG_covconn + [(' C6 ', ' N6 ')]
nu_covconn['  G'] = AG_covconn + [(' C2 ', ' N2 '), (' C6 ', ' O6 ')]
TCU_atomnames = [ ' N1 ', ' C2 ', ' O2 ', ' N3 ', ' C4 ', ' C5 ', ' C6 ', ]
TCU_covconn = [
(' N1 ', ' C2 '),
(' N3 ', ' C2 '),
(' N3 ', ' C4 '),
(' C4 ', ' C5 '),
(' C6 ', ' C5 '),
(' C6 ', ' N1 '),
]
nu_atomnames['  T'] = TCU_atomnames + [' O4 ', ' C5M']
nu_covconn['  T'] = TCU_covconn + [(' C5 ', ' C5M'), (' C4 ', ' O4 ')]
nu_atomnames['  C'] = TCU_atomnames + [' N4 ']
nu_covconn['  C'] = TCU_covconn + [(' C4 ', ' N4 ')]
nu_atomnames['  U'] = TCU_atomnames + [' O4 ']
nu_covconn['  U'] = TCU_covconn + [(' C4 ', ' O4 ')]
P_sugar_covconn = [(' P  ', ' O5*')]
sugar_base_covconn = {}
sugar_base_covconn['  A'] = [(' N9 ', ' C1*')]
sugar_base_covconn['  G'] = sugar_base_covconn['  A']
sugar_base_covconn['  T'] = [(' N1 ', ' C1*')]
sugar_base_covconn['  C'] = sugar_base_covconn['  T']
sugar_base_covconn['  U'] = sugar_base_covconn['  T']

for nu in NU5 :
	knownResinfo.append(Resinfo(nu, nu[2],
		nu_atomnames[nu] + phosphate_atomnames + ribose_sugar_atomnames,
		nu_covconn[nu] + phosphate_covconn + ribose_sugar_covconn + P_sugar_covconn + sugar_base_covconn[nu],
		'rna', nu_synonyms, [' ' + 'R' + nu[2]])) # rna
	knownResinfo.append(Resinfo(nu, nu[2],
		nu_atomnames[nu] + phosphate_atomnames + deoxyribose_sugar_atomnames,
		nu_covconn[nu] + phosphate_covconn + deoxyribose_sugar_covconn + P_sugar_covconn + sugar_base_covconn[nu],
		'dna', nu_synonyms, [' ' + 'D' + nu[2]])) # dna

def getResidueCategory(resn, atomnames=None) :
	if len(resn) > 3 : assert 0==1
	if len(resn) == 1 : resn = ' ' + ' ' + resn
	elif len(resn) == 2 : resn = ' ' + resn
	ret = []
	for ri in knownResinfo :
		if not resn == ri.name3 and not resn in ri.name3syn : continue
		if atomnames != None :
			if ri.consistsOfAtoms(atomnames) : ret.append(ri.type)
		else : ret.append(ri.type)
	if ret == [] : ret = ['het']
	return ret

def isKnownResidue(resn) :
	retri = []
	for ri in knownResinfo :
		if resn == ri.name3 or resn in ri.name3syn : retri.append(ri.name3)
	if len(retri) > 0 : return retri
	return None
	
def sameResidueCovbonded(atomname1, atomname2, resn) :
	assert(len(atomname1) == 4)
	assert(len(atomname2) == 4)
	knownRes = whichResinfo(resn)
	if not knownRes :
		print 'unknown residue -%s-' % resn
		return None
	for kr in knownRes :
		if not kr.containsAtom(atomname1) : continue
		if not kr.containsAtom(atomname2) : continue
		if not kr.bondedAtoms(atomname1,atomname2) : continue
		return 1
	return None

def whichResinfo(resn):
	retri = []
	for ri in knownResinfo :
		if resn == ri.name3 or resn in ri.name3syn : retri.append(ri)
	if len(retri) > 0 : return retri
	return None

def isPhosphateAtom(atomname, resn) :
	cats = getResidueCategory(resn)
	if not 'rna' in cats and not 'dna' in cats : return None
	RIs = whichResinfo(resn)
	rind = -1
	for ri in range(len(RIs)) :
		if RIs[ri].type == 'dna' or RIs[ri].type == 'rna' :
			rind = ri
			break
	assert(rind >= 0)
	if atomname == ' P  ' or atomname in RIs[rind].atsyn[' P  '] : return 1
	if atomname == ' OP1' or atomname in RIs[rind].atsyn[' OP1'] : return 1
	if atomname == ' OP2' or atomname in RIs[rind].atsyn[' OP2'] : return 1
	return None

def isSugarAtom(atomname, resn) :
	cats = getResidueCategory(resn)
	if not 'rna' in cats and not 'dna' in cats : return None
	RIs = whichResinfo(resn)
	rind = -1
	for ri in range(len(RIs)) :
		if RIs[ri].type == 'dna' or RIs[ri].type == 'rna' :
			rind = ri
			break
	assert(rind >= 0)
	for atn in RIs[rind].atomnames :
		if atn[3] != '*' : continue
		if atomname == atn or atomname in RIs[rind].atsyn[atn] : return 1
	return None

def isBaseAtom(atomname, resn) :
	cats = getResidueCategory(resn)
	if not 'rna' in cats and not 'dna' in cats : return None
	if not isSugarAtom(atomname, resn) and not isPhosphateAtom(atomname,resn) : return 1
	return None

def isMCatom(atomname, resn) :
	if 'aa' not in getResidueCategory(resn) : return None
	RIs = whichResinfo(resn)
	assert(len(RIs) == 1)
	if atomname == ' N  ' or atomname in RIs[0].atsyn[' N  '] : return 1
	if atomname == ' CA ' or atomname in RIs[0].atsyn[' CA '] : return 1
	if atomname == ' C  ' or atomname in RIs[0].atsyn[' C  '] : return 1
	if atomname == ' O  ' or atomname in RIs[0].atsyn[' O  '] : return 1
	return None

def main() :
	for ri in knownResinfo :
		print 'KnownResidue', ri.name3, ri.name3syn

if __name__=='__main__' :
	main()
