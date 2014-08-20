from samplers import OmegaSampler, PhipsiSampler, InformedPhipsiSampler, PRLsampler, BBdepChiSampler, NAsuiteSampler, SCLsampler
from pdbr import makeResid, isAAres
from builders import VecInt, VecFloat, VecBuilder, SCLbuilder, VecVecFloat
from builders import PeptideBridgeBuilder, NanchorBuilder, PeptideBuilder, ProtResBuilder, InformedPeptideBuilder, ChiBuilder, Rotator, TransRotator, BuilderGroup, NAsuiteBuilder, CBbuilder
from data import consts, resAtoms
from restraints import SphPosRestr, CentroidPosRestraint, DistanceRestraint, AngleRestraint, EDrestraint, DihedRestraint
from geometry import calcDist

import os, math, string

cbDatapath = os.environ["RTKROOT"] + "/data/"
fwdRatDataPath = cbDatapath + "/rat/fwd.0.000001"
bwdRatDataPath = cbDatapath + "/rat/bwd.0.000001"
nas = NAsuiteSampler("%s/rnaSuite" % cbDatapath, 1)

def addBridgeRestraints(irlist, brgind, res, resids) :
    bridgeapaprt = 4.
    irlist.append( DistanceRestraint(VecInt([res[brgind-1][' CA '],res[brgind+1][' CA ']]), "BridgeApart [%s][%s] %f" % (resids[brgind-1],resids[brgind+1],bridgeapaprt), bridgeapaprt, 1e5) )
    irlist.append( AngleRestraint(VecInt([res[brgind-2][' CA '],res[brgind-1][' CA '],res[brgind][' CA ']]), "BridgeAngle", 70, 150) )
    irlist.append( AngleRestraint(VecInt([res[brgind][' CA '],res[brgind+1][' CA '],res[brgind+2][' CA ']]), "BridgeAngle", 70, 150) )

def buildAcb(res, resids, resnums, resns, chids, incsodes, pts, ri) :
    if not isAAres(resns[ri]) or resns[ri] == "GLY" : return
    if not res[ri].has_key(" CB ") : pts.append([0.,0.,0.]) ; res[ri][" CB "] = len(pts)-1
    newpts, ipinds, opinds = [], VecInt([0,1,2]), VecInt([3,])
    newpts.append( pts[ res[ri][' N  '] ] )
    newpts.append( pts[ res[ri][' C  '] ] )
    newpts.append( pts[ res[ri][' CA '] ] )
    newpts.append( pts[ res[ri][' CB '] ] )
    newpts = VecVecFloat(newpts)
    CBbuilder(ipinds, opinds, consts, "CB builder", resns[ri]).build (newpts)
    pts[res[ri][' CB ']] = [ newpts[3][0], newpts[3][1], newpts[3][2] ]

def buildCBs(pdbfilename) :
    from pdbr import protein
    import prot2res
    from peptidebuild import ModelRenderer
    prot = protein(pdbfilename, read_hydrogens=1, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    for ri in resids.keys() :
        if not isAAres(resns[ri]) or resns[ri] == "GLY" : continue
        buildAcb(res, resids, resnums, resns, chids, inscodes, pts, ri)
    ModelRenderer(res, resns, chids, resnums, inscodes, [], pdbfilename).render(pts)


def printResAtoms(res, resids) :
    print "\n\nResidues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print "%5d:"%k, resids[k], res[k]
    print "------------------------------------------------\n\n\n"

## except for ignoreResidues (the loop residues) make mc and sc restraints on all residues
class CArestraintGenerator :
    def __init__(s, ignoreResidues, caRad, scRad) :
        s.ignoreResidues, s.caRad, s.scRad = ignoreResidues, caRad, scRad
    def generate(s, blist, aiSC, res, resids, pts) : ## blist has no role really
        rlist = []
        for ri in res.keys() :
            if resids[ri] in s.ignoreResidues : continue
            if not " CA " in res[ri].keys() : continue ## a mild check for being amino acid
            if aiSC[ res[ri][" CA "] ] != -1 :
                rlist.append( makeSphPosRestr(ri, ri, 0., s.caRad, res, resids, pts, ' CA ') )
            isdummy = None
            for an,ai in res[ri].items() : ## dont put pos restr on dummies please
                if aiSC[ai] == -1 : isdummy == 1 ; break
            if isdummy == None :
                ar = makeSCcentroidRestr(res, pts, resids, s.scRad, ri)
                if ar != None : rlist.append( ar )
        return rlist, []
        atomri = {}
        for ri in res.keys() :
            for an,ai in res[ri].items() : atomri[ai] = (an,ri)
        rlist, optional = [], []
        for b in blist :
            bop = b.getOP()
            for i in range(bop.size()) :
                pi = bop[i]
                if atomri[pi][0] != " CA " or aiSC[pi] == -1 : continue
                rlist.append( makeSphPosRestr(atomri[pi][1], atomri[pi][1], 0, s.caRad, res, resids, pts) )
        return rlist, []

class XrayRestraintsGenerator :
    def __init__(s, mtzfn, folabel="Fo", fclabel="Fc", philabel="PHIC", maptype="2F1-F2", min=.25, max=2., mean=1., exclBldr=[], makeOptionalRestraints=None) :
        s.makeOptionalRestraints = makeOptionalRestraints
        if folabel == "map" :
            s.type = "map"
            s.mapfn = mtzfn
        else :
            s.type = "mtz"
            if maptype == '2F1-F2' : maptype = 0
            elif maptype == 'F1' : maptype = 1
            else : print "unknown maptype ", maptype ; sys.exit(1)
            s.mtzfn, s.folabel, s.fclabel, s.philabel, s.maptype = mtzfn, folabel, fclabel, philabel, maptype 
        s.min, s.max, s.mean = min, max, mean
        s.excludedBuilders = exclBldr
    def generate(s, blist, aiSC, res, resids, pts) :
        if s.type == "mtz" : return s.generateMTZ(blist, aiSC)
        elif s.type == "map" : return s.generateMAP(blist, aiSC)
        else : assert None
    def isOptional(s,b) :
        if s.makeOptionalRestraints : return 1
        else : return None
        s.compulsory = ["LigandBuilder",] #"PeptideBridgeBuilder","PeptideBuilder"]
        for bn in s.compulsory :
            if bn in b.name() :
                print b.name(), "is compulsory, not optional" ; return None
        return 1
    def ignoreBuilder(s, b) :
        ignB = None
        for exn in s.excludedBuilders :
            if exn in b.name() : ignB = 1
        return ignB
    def generateMTZ(s, blist, aiSC) :
        xrlist, optional = [], []
        for b in blist :
            if s.ignoreBuilder(b) : continue
            bop = b.getOP()
            edrIP = []
            for i in range(bop.size()) :
                if aiSC[bop[i]] != -1 : edrIP.append(bop[i])
            if len(edrIP) > 0 :
                xrlist.append( EDrestraint.makeEDrestraintFromMTZ(VecInt(edrIP), "EDrestraint %4.2f %4.2f %4.2f on o/p of %s" % (s.min,s.max,s.mean,b.name()), \
                    s.mtzfn, s.folabel, s.fclabel, s.philabel, s.maptype, s.min, s.max, s.mean) )
                if s.isOptional(b) : optional.append(len(xrlist)-1)
        return xrlist, optional
    def generateMAP(s, blist, aiSC) : ##TODO
        print s.mapfn, s.min, s.max, s.mean
        erlist, optional = [], []
        for b in blist :
            if s.ignoreBuilder(b) : continue
            bop = b.getOP()
            edrIP = []
            for i in range(bop.size()) :
                if aiSC[bop[i]] != -1 : edrIP.append(bop[i])
            rstr = "XEnvelopeRestraint %4.2f %4.2f %4.2f on op of %s" % (s.min, s.max, s.mean, b.name())
            print s.mapfn, s.min, s.max, s.mean
            erlist.append( EDrestraint.makeEDrestraintFromMap(VecInt(edrIP), rstr, s.mapfn, s.min, s.max, s.mean) )
            if s.isOptional(b) : optional.append(len(erlist)-1)
        return erlist, optional

class PRL_BB_SC_SamplerProvider :
    def __init__(self) :
        self.samplers, self.newsamplers = {}, []
        import protinfo, copy
        self.aaMap = copy.deepcopy(protinfo.AA31)
        self.aaMap["MSE"] = "M"
    def get(self, resn, addsample=None) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = PRLsampler( "%s/richardson.lib" % cbDatapath, self.aaMap[resn] )
        if addsample != None :
            newsam = self.samplers[resn].makeCopy()
            assert len(addsample) == 2 and len(addsample[1]) >= 4
            newsam.addSample(addsample[0], VecFloat(addsample[1]))
            self.newsamplers.append(newsam)
            return newsam
        return self.samplers[resn]

class SCL_BB_SC_SamplerProvider :
    def __init__(self, resolution) :
        self.samplers, self.newsamplers, self.minSamples = {}, [], {}
        self.crdfile  = "%s/SCL/scl-B30-occ1.0-rmsd%3.1f-prop20.0.pdb" % (cbDatapath, resolution)
        self.propfile = "%s/SCL/scl-B30-occ1.0-rmsd%3.1f-prop20.0.dat" % (cbDatapath, resolution)
        import protinfo, copy
        self.aaMap = copy.deepcopy(protinfo.AA31)
        self.aaMap["MSE"] = "M"
    def get(self, resn, addsample=None) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = SCLsampler( self.crdfile, self.propfile, self.aaMap[resn] )
        if addsample != None :
            newsam = self.samplers[resn].makeCopy()
            assert len(addsample) == 2
            newsam.addSample(addsample[0], VecVecFloat(addsample[1]))
            self.newsamplers.append(newsam)
            return newsam
        return self.samplers[resn]
    def findMinNumSamples(self, resn) :
        if not resn in self.minSamples.keys() :
            sclsampler = self.get(resn)
            retpts, maxRI = None, -1000000 ##XXX
            for phi in range(-180,180,40) :
                for psi in range(-180,180,40) :
                    nrotPS = -1
                    for ri in range(100000) :
                        tmp = sclsampler.sample(phi, psi, retpts, ri)
                        if tmp < 0 or tmp != ri : break
                        nrotPS = tmp
                    #print phi, psi, resn, maxRI, nrotPS
                    nrotPS += 1
                    if nrotPS <= 0 : continue
                    if nrotPS > maxRI : maxRI = nrotPS ##XXX
            print "SCL-min-numsamples", resn , maxRI
            assert maxRI > 0 and maxRI < 1000
            self.minSamples[resn] = maxRI
        return self.minSamples[resn]

class Dunbrack_BB_SC_SamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = BBdepChiSampler( "%s/BBdepSCrot/sc%s" % (cbDatapath,resn) )
        return self.samplers[resn]

class RapSamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = RapSampler( cbDatapath + "/rapPSO", resn )
            print "read RapSampler", resn
        return self.samplers[resn]

class PhipsiSamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn, ss="coil") :
        if ss == "coil" : sstype = 0
        elif ss == "alpha" : sstype = 1
        elif ss == "beta" : sstype = 2
        else : assert None
        sstype = 0
        if not resn in self.samplers.keys() : self.samplers[resn] = {}
        if not sstype in self.samplers[resn].keys() :
            self.samplers[resn][sstype] = PhipsiSampler( cbDatapath + "/PhipsiWeightedProp/ps%s" % resn , sstype )
        return self.samplers[resn][sstype]

class OmegaSamplerProvider :
    def __init__(self) :
        self.omegaSampler = OmegaSampler(0.99)
        self.preproOmegaSampler = OmegaSampler(0.90)
    def get(self, preOmegaResn, postOmegaResn) :
        if postOmegaResn == 'PRO' : return self.preproOmegaSampler
        return self.omegaSampler

def makeSphPosRestr(onIndex, duetoIndex, Rlow, Rup, res, resids, pts, atomname=' CA ') :
    descr = "SphPosRestraint on [%s]%s %f %f due to [%s]" % (resids[onIndex], atomname, Rlow, Rup, resids[duetoIndex])
    return SphPosRestr( VecInt([res[onIndex][atomname]]), descr, VecFloat(pts[res[duetoIndex][atomname]]), Rlow, Rup )


def makeCBbuilder(index, res, resns, resids) :
    IPs, OPs = [ res[index][' N  '],res[index][' C  '],res[index][' CA '] ], [ res[index][' CB '] ]
    return CBbuilder(VecInt(IPs), VecInt(OPs), consts, "CBbuilder %s" % (resids[index]), resns[index])

## make Nanchor-builder at beginning and CO builder in the end
psp = PhipsiSamplerProvider() ## make peptide sampler
omsp = OmegaSamplerProvider() ## make peptide sampler

class PrepareChain :
    def __init__(s, sctype) :
        s.scl, s.prl = None, None
        if sctype[0:3] == "PRL" :
            s.prl = PRL_BB_SC_SamplerProvider() ## make sidechain sampler with Lovell library
        elif sctype[0:3] == "SCL" :
            s.scl = SCL_BB_SC_SamplerProvider(string.atof(sctype[3:])) ## make sidechain sampler with Shetty library
        s.isamps = []
    def makeChiBuilder(s, i, j, k, res, resns, resids, addsample=None) :
        IPs, OPs = [ res[i][' C  '],res[j][' N  '],res[j][' CA '],res[j][' C  '],res[k][' N  '],res[j][' CB '] ], []
        for name in resAtoms[resns[j]]:
            if name in [ ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ] : continue
            OPs.append(res[j][name])

        if s.scl :
            if addsample != None :
                scsampler = s.scl.get(resns[j], addsample)
            else : scsampler = s.scl.get(resns[j])
        elif s.prl :
            if addsample != None :
                if len(addsample[1]) == 1 : addsample[1].append(-999.)
                if len(addsample[1]) == 2 : addsample[1].append(-999.)
                if len(addsample[1]) == 3 : addsample[1].append(-999.)
                scsampler = s.prl.get(resns[j], addsample)
            else : scsampler = s.prl.get(resns[j])
        else : assert None
        
        if s.scl :
            return SCLbuilder(VecInt(IPs), VecInt(OPs), consts, "SCLbuilder %s" % (resids[j]), scsampler)
        elif s.prl :
            return ChiBuilder(VecInt(IPs), VecInt(OPs), consts, "ChiBuilder %s" % (resids[j]), resns[j], scsampler)
        else : assert None

    ## i < j < k if fwd; i < j < k is bwd.
    ## k'th CA built in fwd mode, i'th in bwd mode
    def makeBgrp(s, i, j, k, pts, res, resns, resids, mconly, fwd, informedRadius) :
        '''make peptide builder that builds k'th CA, j'th CB and sidechain, assuming i,j,k are one after other indices into res'''
        assert mconly == 1 or mconly == None
        assert fwd == 1 or fwd == 0
        bgrp = []

        if fwd == 1 : IPs, OPs = [ res[i][' C  '],res[j][' N  '],res[j][' CA '] ], [ res[j][' C  '],res[j][' O  '],res[k][' N  '],res[k][' CA '] ]
        else : IPs, OPs = [ res[k][' N  '],res[j][' C  '],res[j][' CA '] ], [ res[j][' N  '],res[i][' C  '],res[i][' O  '],res[i][' CA '] ]
    
        #print resids[j], res[j]
        if mconly != 1 and resns[j] != "GLY" : OPs.append(res[j][' CB '])
    
        if fwd == 1 : bstr = "PeptideBuilder [%s]" % resids[k]
        else : bstr = "R-PeptideBuilder [%s]" % resids[i]
        if informedRadius : bstr = "Inf-" + bstr
    
        if not mconly : bstr = "%s sc%s" % (bstr, resns[j])
    
        if informedRadius :
            if fwd == 1 :
                s.isamps.append( InformedPhipsiSampler(VecFloat(pts[res[k][' CA ']]), informedRadius, resns[j], fwdRatDataPath, 1) )
                bgrp.append( InformedPeptideBuilder( VecInt(IPs), VecInt(OPs), consts, resns[j], bstr, psp.get(resns[j]), omsp.get(resns[j],resns[k]), s.isamps[len(s.isamps)-1], 1 ) )
            else :
                s.isamps.append( InformedPhipsiSampler(VecFloat(pts[res[i][' CA ']]), informedRadius, resns[j], bwdRatDataPath, 0) )
                bgrp.append( InformedPeptideBuilder( VecInt(IPs), VecInt(OPs), consts, resns[j], bstr, psp.get(resns[j]), omsp.get(resns[i],resns[j]), s.isamps[len(s.isamps)-1], 0 ) )
        else :
            if fwd == 1 : bgrp.append( PeptideBuilder( VecInt(IPs), VecInt(OPs), consts, resns[j], bstr, psp.get(resns[j]), omsp.get(resns[j],resns[k]), 1 ) )
            else : bgrp.append( PeptideBuilder( VecInt(IPs), VecInt(OPs), consts, resns[j], bstr, psp.get(resns[j]), omsp.get(resns[i],resns[j]), 0 ) )
    
        if not mconly and resns[j] not in ['GLY','ALA'] : # can build sidechain
            bgrp.append( s.makeChiBuilder(i,j,k, res, resns, resids) )
        return bgrp

    def preparePeptideChain(s, chid, res, resids, resnums, resns, chids, inscodes, pts, mconly, buildfwd, guidedSamplingRadius, caRad, scRad, scMissInds=None) :
        assert mconly == 1 or mconly == None
        ## for a chain in pdbfile, rebuild using some CA-sphere restraint radius
        ## assume all else in pdbfile as just steric obstructions
        assert len(chid) == 1
    
        if scMissInds == None and mconly == None : scMissInds = incompleteSCcorrection(res, resns, pts)
    
        loop_started, index, startindex, endindex = None, 0, None, None
        keys = res.keys() ; keys.sort()
        for index in keys :
            if chids[index] == chid and not loop_started : loop_started, startindex = 1, index
            if chids[index] != chid and loop_started : loop_started, endindex = None, index-1
        if startindex != None and endindex == None : endindex = keys[len(keys)-1]
        assert startindex != None and endindex != None and startindex < endindex
        for index in range(startindex,endindex+1) : assert chid == chids[index]
    
        dummies,firstindex,lastindex = addNCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
    
        blist, bfoll, numtrials, rlist = [], {}, [], []
        useSCfollowers = 1
    
        assert buildfwd == 1 or buildfwd == 0
        if buildfwd == 1 :
            IPs, OPs = [], [ res[startindex][' CA '], res[startindex][' C  '], res[startindex][' O  '], res[startindex+1][' N  '], res[startindex+1][' CA '], ]
            blist.append( NanchorBuilder( VecInt(IPs), VecInt(OPs), consts, "NanchorBuilder [%s][%s]"%(resids[startindex],resids[startindex+1]), VecFloat(pts[res[startindex][' CA ']]), caRad, VecFloat(pts[res[startindex+1][' CA ']]), caRad ) )
            rest1, rest2 = makeSCcentroidRestr(res,pts,resids,scRad,startindex), makeSCcentroidRestr(res,pts,resids,scRad,startindex+1)
            if rest1 and not startindex in scMissInds : rlist.append(rest1)
            if rest2 and not startindex+1 in scMissInds : rlist.append(rest2)
            bgrp = s.makeBgrp(firstindex, startindex, startindex+1, pts, res, resns, resids, mconly, 0, None)
            for bl in bgrp : blist.append(bl)
            if len(bgrp) > 1 and useSCfollowers : #XXX
                bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
            blist1, bfoll1, numtrials1, rlist1 = s.prepareChainTerminal("Cterm", startindex+2, endindex, lastindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, guidedSamplingRadius, scMissInds)
            Rlow, Rup = 3.81-caRad, 3.81+caRad
            if Rlow < 0 : Rlow = 0
            rlist1.append( makeSphPosRestr(startindex+1, startindex+2, Rlow, Rup, res, resids, pts) )
        else :
            IPs, OPs = [], [ res[endindex-1][' CA '], res[endindex-1][' C  '], res[endindex-1][' O  '], res[endindex][' N  '], res[endindex][' CA '], ]
            blist.append( NanchorBuilder( VecInt(IPs), VecInt(OPs), consts, "NanchorBuilder [%s][%s]"%(resids[endindex-1],resids[endindex]), VecFloat(pts[res[endindex-1][' CA ']]), caRad, VecFloat(pts[res[endindex][' CA ']]), caRad ) )
            rest1, rest2 = makeSCcentroidRestr(res,pts,resids,scRad,endindex-1), makeSCcentroidRestr(res,pts,resids,scRad,endindex)
            if rest1 and not endindex-1 in scMissInds : rlist.append(rest1)
            if rest2 and not endindex in scMissInds : rlist.append(rest2)
            bgrp = s.makeBgrp(endindex-1, endindex, lastindex, pts, res, resns, resids, mconly, 1, None)
            for bl in bgrp : blist.append(bl)
            if len(bgrp) > 1 and useSCfollowers : #XXX
                bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
            blist1, bfoll1, numtrials1, rlist1 = s.prepareChainTerminal("Nterm", startindex, endindex-2, firstindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, guidedSamplingRadius, scMissInds)
            Rlow, Rup = 3.81-caRad, 3.81+caRad
            if Rlow < 0 : Rlow = 0
            rlist1.append( makeSphPosRestr(endindex-1, endindex-2, Rlow, Rup, res, resids, pts) )
    
        numtrials = s.makeNumtrials(blist)
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)

        return blist, numtrials, bfoll, rlist, dummies
    def prepareChainTerminal(s, NCterm, startindex, endindex, specialIndex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, guidedSamplingRadius, scMissInds=[]) :
        assert NCterm in ["Nterm", "Cterm"]
        assert startindex <= endindex
        #print "prepareChainTerminal", NCterm, startindex, endindex, specialIndex, "[%s] [%s] [%s]" % (resids[startindex], resids[endindex], resids[specialIndex]), res[specialIndex]
        blist, bfoll = [], {}
        useSCfollowers = 1
        if NCterm == "Cterm" :
            lastindex = specialIndex
            # for index, make PeptideBuilder that builds index'th CA
            # for index, make CBbuilder and ChiBuilder that builds index-1'th CB and sidechain
            # startindex : startindex + 1 : .... endindex : lastindex
            indices = range(startindex, endindex+1) + [lastindex]
            for index in indices :
                i, j, k, informed = index-2, index-1, index, guidedSamplingRadius
                #print i, j, k
                if index == lastindex : i,j,k, informed = endindex-1, endindex, lastindex, None
                bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 1, informed)
                for b in bgrp : blist.append(b)
                if len(bgrp) > 1 and useSCfollowers :
                    if index == endindex and resns[index] in ["ALA","GLY"] : pass
                    else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
        else :
            #tempIndex = startindex ; startindex = endindex ; endindex = tempIndex
            firstindex = specialIndex
            # for index, make PeptideBuilder that builds index'th CA
            # for index, make CBbuilder and ChiBuilder that builds index+1'th CB and sidechain
            # endindex : endindex - 1 : ... : startindex : firstindex
            indices = range(endindex, startindex-1, -1) + [firstindex]
            for index in indices :
                i, j, k, informed = index, index+1, index+2, guidedSamplingRadius
                if index == firstindex : i,j,k,informed = firstindex, startindex, startindex+1, None
                bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 0, informed)
                for b in bgrp : blist.append(b)
                if len(bgrp) > 1 and useSCfollowers :
                    if index == startindex and resns[index] in ["ALA","GLY"] : pass
                    else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
    
        numtrials = s.makeNumtrials(blist)
    
        ## make CA-tracing restraints. common for both fwd and bkwd builds
        irlist = []
        for index in range(startindex, endindex+1) :
            ## CA restraint
            irlist.append( makeSphPosRestr(index, index, 0, caRad, res, resids, pts) )
            if NCterm == "Cterm" :
                indmax = endindex
                if index < endindex - 5 : indmax = index + 5
                proprange = range(index+1,indmax+1, 1)
            else :
                indmin = startindex
                if index > startindex + 5 : indmin = index-5
                proprange = range(index-1,indmin-1, -1)
            for ind in proprange : # backpropagate CA restraints in the at-most-so-far and at-least-this-far sense
                #d_low, d_up = 0., abs(index-ind) * 3.81 + caRad
                d_low, d_up = 0., ((abs(index-ind)/2)*2*math.sin(math.pi*70./180.) + abs(index-ind)%2) * 3.81 + caRad
                if abs(index-ind) == 1 : #and guidedSamplingRadius :
                    dCA = calcDist( VecFloat(pts[res[index][' CA ']]) , VecFloat(pts[res[ind][' CA ']]) )
                    if abs(dCA-3.81) < abs(dCA-2.8) :
                        if 3.81-caRad > 0 : d_low = 3.81-caRad
                        if guidedSamplingRadius : d_up = 3.81+caRad
                    else :
                        if 2.8-caRad > 0 : d_low = 2.8-caRad
                        if guidedSamplingRadius : d_up = 2.8+caRad
                irlist.append( makeSphPosRestr(index, ind, d_low, d_up, res, resids, pts) )
            ## sidechain centroid restraints
            IPs, mean, count = [], [0.,0.,0.], 0
            for k,v in res[index].items() :
                if k in (' N  ',' CA ',' C  ',' O  ',) : continue
                IPs.append(v) ; count = count + 1
                mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
            if mconly == 1 or count == 0 or index in scMissInds : continue
            mean[0] /= count ; mean[1] /= count ; mean[2] /= count
            irlist.append( CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[index],scRad), VecFloat(mean), scRad) )
    
        return blist, bfoll, numtrials, irlist
    def makeNumtrials(s, blist) :
        numtrials = []
        from data import numRotPRL, numRotSCL
        for b in blist :
            btype = "%s" % type(b)
            if "PeptideBuilder" in btype : numtrials.append(1000)
            elif "PeptideBridgeBuilder" in btype :
                numtrials.append(1000)
                from builders import PeptideBridgeBuilder ; PeptideBridgeBuilder.setCTtrials(250); PeptideBridgeBuilder.setThetastep(5);
                from data import consts ; consts.set("TAU_QUALITY", 30.)
            elif "NanchorBuilder" in btype : numtrials.append(1000)
            elif "ChiBuilder" in btype : numtrials.append(numRotPRL[b.name()[11:14]])
            elif "SCLbuilder" in btype : numtrials.append( s.scl.findMinNumSamples(b.name()[11:14]) )
            elif "CBbuilder" in btype : numtrials.append(1)
            else : print btype ; assert None
        return numtrials
    def preparePeptideMidloop1(s, startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, scMissInds=[], guidedRadius=None) :
        if not mconly in [1,None] : print mconly ; assert None
        blist, bfoll = [], {}
        useSCfollowers = 1
        brgind = int( math.ceil((startindex + endindex) / 2.) )
        ## make a fwd and a bwd peptide builder, leading to the bridge residue
        for index in range(startindex, brgind) :
            i, j, k = index-2, index-1, index
            bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 1, guidedRadius)
            for b in bgrp : blist.append(b)
            if len(bgrp) > 1 and useSCfollowers :
                if index == endindex and resns[index] in ["ALA","GLY"] : pass
                else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
            ind = endindex - (index-startindex)
            if ind == brgind : continue
            i, j, k = ind, ind+1, ind+2
            bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 0, guidedRadius)
            for b in bgrp : blist.append(b)
            if len(bgrp) > 1 and useSCfollowers :
                if index == endindex and resns[index] in ["ALA","GLY"] : pass
                else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
        ## bridge
        IPs = [ res[brgind-2][' C  '], res[brgind-1][' N  '], res[brgind-1][' CA '], res[brgind+1][' CA '], res[brgind+1][' C  '], ]
        OPs = [ res[brgind-1][' C  '], res[brgind-1][' O  '], res[brgind+0][' N  '], res[brgind+0][' CA '],
                res[brgind+0][' C  '], res[brgind+0][' O  '], res[brgind+1][' N  '], ]
        for ei in range(brgind-1,brgind+2) :
            if not mconly and resns[ei] != "GLY" : OPs += [ res[ei][' CB '], ]
        bridgeB = PeptideBridgeBuilder(VecInt(IPs), VecInt(OPs), consts, "PeptideBridgeBuilder [%s][%s]"%(resids[brgind-1],resids[brgind]),
            resns[brgind-1], resns[brgind], resns[brgind+1], fwdRatDataPath)
        blist.append(bridgeB) ; bridgeIndex = len(blist)-1 ; bfoll[bridgeIndex] = []
        if not mconly :
            for ei in [brgind, brgind-1, brgind+1] :
                if resns[ei] != "GLY" :
                    if resns[ei] != "ALA" :
                        blist.append( s.makeChiBuilder(ei-1, ei, ei+1, res, resns, resids) )
                        bfoll[ bridgeIndex ].append( len(blist)-1 )
        if len( bfoll[bridgeIndex] ) == 0 : del bfoll[bridgeIndex]
        numtrials = s.makeNumtrials(blist)
    
        irlist = []
        for index in range(startindex, brgind+1) :
            for pind in range(index,brgind+1) : irlist.append( makeSphPosRestr(index, pind, 0, (pind-index)*3.81+caRad, res, resids, pts) )
            ind1 = endindex + 1 - (index-startindex)
            if ind1 > index :
                propdist = 3.81*(ind1-index)+0.01 #     ( 2*math.sin(math.pi*70./180.)*((ind1-index)/2) + (ind1-index)%2 )
                irlist.append( DistanceRestraint(VecInt([res[ind1][' CA '],res[index][' CA ']]), "LoopClosure [%s][%s] %f" % (resids[ind1],resids[index],propdist), 2.8, propdist) )
            endi = endindex - (index-startindex)
            for pind in range(brgind, endi+1) : irlist.append( makeSphPosRestr(endi, pind, 0, (endi-pind)*3.81+caRad, res, resids, pts) )
            if endi > index :
                propdist = 3.81*(endi-index)+0.01 #   ( 2*math.sin(math.pi*70./180.)*((endi-index)/2) + (endi-index)%2 )
                irlist.append( DistanceRestraint(VecInt([res[endi][' CA '],res[index][' CA ']]), "LoopClosure [%s][%s] %f" % (resids[endi],resids[index],propdist), 2.8, propdist) )
        addBridgeRestraints(irlist, brgind, res, resids)
        if not mconly :
            for index in range(startindex-1, endindex+2) :
                if index in scMissInds : continue
                if resns[index] in ["ALA","GLY"] : continue
                IPs, mean, count = [], [0.,0.,0.], 0
                for k,v in res[index].items() :
                    if k in (' N  ',' CA ',' C  ',' O  ',) : continue
                    IPs.append(v) ; count = count + 1
                    mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
                mean[0] /= count ; mean[1] /= count ; mean[2] /= count
                irlist.append( CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[index],scRad), VecFloat(mean), scRad) )
        return blist, bfoll, numtrials, irlist
    def preparePeptideMidloop(s, startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, scMissInds=[], guidedRadius=None) :
        if not mconly in [1,None] : print mconly ; assert None
        blist, bfoll = [], {}
        useSCfollowers = 1
        ## make a fwd and a bwd peptide builder, leading to the bridge residue
        builderDone, lastfwdind = [], None
        for index in range(startindex, endindex) :
            if index in builderDone : break
            i, j, k = index-2, index-1, index
            bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 1, guidedRadius)
            for b in bgrp : blist.append(b) ; builderDone.append(index) ; lastfwdind = index
            if len(bgrp) > 1 and useSCfollowers :
                if index == endindex and resns[index] in ["ALA","GLY"] : pass
                else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
            ind = endindex - (index-startindex)
            if ind in builderDone : break
            i, j, k = ind, ind+1, ind+2
            bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 0, guidedRadius)
            for b in bgrp : blist.append(b) ; builderDone.append(ind)
            if len(bgrp) > 1 and useSCfollowers :
                if index == endindex and resns[index] in ["ALA","GLY"] : pass
                else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
        brgind = lastfwdind + 1
        ## bridge
        IPs = [ res[brgind-2][' C  '], res[brgind-1][' N  '], res[brgind-1][' CA '], ]
        OPs = [ res[brgind-1][' C  '], res[brgind-1][' O  '], res[brgind+0][' N  '], ]
        blist.append( PeptideBuilder( VecInt(IPs), VecInt(OPs), consts, resns[brgind-1], "PeptideCONBuilder [%s] sc%s sc%s"%(resids[brgind],resns[brgind-1],resns[brgind]),
            psp.get(resns[brgind-1]), omsp.get(resns[brgind-1],resns[brgind]), 1 ) )
        brind = len(blist)-1 ; bfoll[brind] = []
        if not mconly and resns[brgind-1] != "GLY" :
            blist.append( makeCBbuilder(brgind-1, res, resns, resids) )
            bfoll[brind].append( len(blist) - 1 )
        if not mconly and resns[brgind] != "GLY" :
            blist.append( makeCBbuilder(brgind, res, resns, resids) )
            bfoll[brind].append( len(blist)-1 )
        if not mconly and not resns[brgind-1] in ["GLY","ALA"] :
            blist.append( s.makeChiBuilder(brgind-2, brgind-1, brgind, res, resns, resids) )
            bfoll[brind].append( len(blist)-1 )
        if not mconly and not resns[brgind] in ["GLY","ALA"] :
            blist.append( s.makeChiBuilder(brgind-1, brgind, brgind+1, res, resns, resids) )
            bfoll[brind].append( len(blist)-1 )

        numtrials = s.makeNumtrials(blist)
    
        irlist = []
        for index in range(startindex, brgind+1) :
            for pind in range(index,brgind+1) : irlist.append( makeSphPosRestr(index, pind, 0, (pind-index)*3.81+caRad, res, resids, pts) )
            ind1 = endindex + 1 - (index-startindex)
            if ind1 > index :
                propdist = 3.81*(ind1-index)+0.01 #     ( 2*math.sin(math.pi*70./180.)*((ind1-index)/2) + (ind1-index)%2 )
                irlist.append( DistanceRestraint(VecInt([res[ind1][' CA '],res[index][' CA ']]), "LoopClosure [%s][%s] %f" % (resids[ind1],resids[index],propdist), 2.8, propdist) )
            endi = endindex - (index-startindex)
            for pind in range(brgind, endi+1) : irlist.append( makeSphPosRestr(endi, pind, 0, (endi-pind)*3.81+caRad, res, resids, pts) )
            if endi > index :
                propdist = 3.81*(endi-index)+0.01 #   ( 2*math.sin(math.pi*70./180.)*((endi-index)/2) + (endi-index)%2 )
                irlist.append( DistanceRestraint(VecInt([res[endi][' CA '],res[index][' CA ']]), "LoopClosure [%s][%s] %f" % (resids[endi],resids[index],propdist), 2.8, propdist) )

        irlist.append( DistanceRestraint(VecInt([res[brgind][' N  '],res[brgind][' CA ']]), "LoopClosure N-CA [%s] 1-2" % (resids[brgind]), .5, 2.) )
        irlist.append( AngleRestraint(VecInt([res[brgind][' N  '],res[brgind][' CA '],res[brgind][' C  ']]), "LoopClosure N-CA-C [%s] 80-140" % (resids[brgind]), 80., 140.) )
        irlist.append( AngleRestraint(VecInt([res[brgind-1][' C  '],res[brgind][' N  '],res[brgind][' CA ']]), "LoopClosure C-N-CA [%s] 90-150" % (resids[brgind]), 90., 150.) )
        irlist.append( DihedRestraint(VecInt([res[brgind-1][' CA '],res[brgind-1][' C  '],res[brgind][' N  '],res[brgind][' CA ']]),
            "LoopClosure CA-C-N-CA [%s] cis/trans" % (resids[brgind]), VecFloat([0,-180]), VecFloat(30,30)) )

        if not mconly :
            for index in range(startindex-1, endindex+2) :
                if index in scMissInds : continue
                if resns[index] in ["ALA","GLY"] : continue
                IPs, mean, count = [], [0.,0.,0.], 0
                for k,v in res[index].items() :
                    if k in (' N  ',' CA ',' C  ',' O  ',) : continue
                    IPs.append(v) ; count = count + 1
                    mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
                mean[0] /= count ; mean[1] /= count ; mean[2] /= count
                irlist.append( CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[index],scRad), VecFloat(mean), scRad) )
        return blist, bfoll, numtrials, irlist
    def preparePeptideLoop1(s, startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, scMissInds=[], guidedRadius=None, taperRestraints=None) :
        if not mconly in [1,None] : print mconly ; assert None
        blist, bfoll = [], {}
        useSCfollowers = 1
        ## make loop before bridge
        for index in range(startindex, endindex) :
            i, j, k = index-2, index-1, index
            bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 1, guidedRadius)
            for b in bgrp : blist.append(b)
            if len(bgrp) > 1 and useSCfollowers :
                if index == endindex and resns[index] in ["ALA","GLY"] : pass
                else : bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
        ## bridge
        IPs = [ res[endindex-2][' C  '], res[endindex-1][' N  '], res[endindex-1][' CA '], res[endindex+1][' CA '], res[endindex+1][' C  '], ]
        OPs = [ res[endindex-1][' C  '], res[endindex-1][' O  '], res[endindex+0][' N  '], res[endindex+0][' CA '],
                res[endindex+0][' C  '], res[endindex+0][' O  '], res[endindex+1][' N  '], ]
        for ei in range(endindex-1,endindex+2) :
            if not mconly and resns[ei] != "GLY" : OPs += [ res[ei][' CB '], ]
        bridgeB = PeptideBridgeBuilder(VecInt(IPs), VecInt(OPs), consts, "PeptideBridgeBuilder [%s][%s]"%(resids[endindex-1],resids[endindex]),
            resns[endindex-1], resns[endindex], resns[endindex+1], fwdRatDataPath)
        blist.append(bridgeB) ; bridgeIndex = len(blist)-1 ; bfoll[bridgeIndex] = []
        if not mconly :
            for ei in [endindex, endindex-1, endindex+1] :
                if resns[ei] != "GLY" :
                    if resns[ei] != "ALA" :
                        blist.append( s.makeChiBuilder(ei-1, ei, ei+1, res, resns, resids) )
                        bfoll[ bridgeIndex ].append( len(blist)-1 )
        if len( bfoll[bridgeIndex] ) == 0 : del bfoll[bridgeIndex]
        numtrials = s.makeNumtrials(blist)
        ## starting from at least 0.5 & 1, caRad & scRad have to be scaled up to be maximum at the center and tapering towards known ends 
        gradCA = (caRad+0.) / ((endindex+1 - startindex+1) / 2.)
        irlist = []
        for index in range(startindex, endindex+1) :
            indDiff = index - startindex + 1
            if endindex + 1 - index < indDiff : indDiff = endindex - index + 1
            rad = caRad
            if taperRestraints : rad = 0.5 + gradCA * indDiff
            irlist.append( makeSphPosRestr(index, index, 0, rad, res, resids, pts) )
            ## min and max dist w.r.t. next CA restraint, esp imp for guided mode
            minR = 3.81 - caRad ;
            if minR < 0 : minR = 0
            irlist.append( makeSphPosRestr(index, index+1, minR, 3.81+caRad, res, resids, pts) )
            #print "adding CA restraint", resids[index], rad, gradCA
            propdist = 3.81*( 2*math.sin(math.pi*70./180.)*((endindex+1-index)/2) + (endindex+1-index)%2 )
            irlist.append( DistanceRestraint(VecInt([res[endindex+1][' CA '],res[index][' CA ']]), "LoopClosure [%s][%s] %f" % (resids[endindex+1],resids[index],propdist), 2.8, propdist) )
        addBridgeRestraints(irlist, endindex, res, resids)
        if not mconly :
            gradSC = (scRad-0.) / ((endindex+2 - startindex+2) / 2)
            for index in range(startindex-1, endindex+2) :
                if index in scMissInds : continue
                if resns[index] in ["ALA","GLY"] : continue
                IPs, mean, count = [], [0.,0.,0.], 0
                for k,v in res[index].items() :
                    if k in (' N  ',' CA ',' C  ',' O  ',) : continue
                    IPs.append(v) ; count = count + 1
                    mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
                mean[0] /= count ; mean[1] /= count ; mean[2] /= count
                indDiff = index - startindex + 2
                if endindex + 2 - index < indDiff : indDiff = endindex - index + 2
                rad = scRad
                if taperRestraints : rad = 1 + gradSC*indDiff
                irlist.append( CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[index],rad), VecFloat(mean), scRad) )
                #print "adding Centroid restraint", resids[index], rad, gradSC
        return blist, bfoll, numtrials, irlist
    def preparePeptideLoop(s, startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, scMissInds=[], guidedRadius=None, taperRestraints=None) :
        if not mconly in [1,None] : print mconly ; assert None
        blist, bfoll = [], {}
        useSCfollowers = 1
        lastpeptbuilder = None
        ## make loop before bridge
        for index in range(startindex, endindex+2) :
            i, j, k = index-2, index-1, index
            bgrp = s.makeBgrp(i,j,k, pts, res, resns, resids, mconly, 1, guidedRadius)
            for b in bgrp : blist.append(b) ; lastpeptbuilder = index
            if len(bgrp) > 1 and useSCfollowers :
                bfoll[len(blist) - len(bgrp)] = range(len(blist)-len(bgrp)+1, len(blist))
        numtrials = s.makeNumtrials(blist)
        irlist = []
        caradii = {}
        for index in range(startindex,endindex+1) : caradii[index] = caRad
        caradii[endindex+1] = 0.5 ; caradii[endindex+2] = 0.25
        for index in range(startindex, endindex+2) :
            irlist.append( makeSphPosRestr(index, index, 0, caradii[index], res, resids, pts) )
            ## min and max dist w.r.t. next CA restraint, esp imp for guided mode
            for nextind in range(index+1,endindex+2) :
                irlist.append( makeSphPosRestr(index, nextind, 0, 3.81*(nextind-index)+caradii[nextind], res, resids, pts) )
            propdist = 3.81*( 2*math.sin(math.pi*70./180.)*((endindex+2-index)/2) + (endindex+2-index)%2 ) + 0.01
            irlist.append( DistanceRestraint(VecInt([res[endindex+2][' CA '],res[index][' CA ']]), "LoopClosure [%s][%s] %f" % (resids[endindex+2],resids[index],propdist), 2.8, propdist) )
        irlist.append( makeSphPosRestr(endindex+1, endindex+1, 0, caradii[endindex+1], res, resids, pts, ' N  ') )
        if not mconly :
            for index in range(startindex-1, endindex+2) :
                if index in scMissInds : continue
                if resns[index] in ["ALA","GLY"] : continue
                IPs, mean, count = [], [0.,0.,0.], 0
                for k,v in res[index].items() :
                    if k in (' N  ',' CA ',' C  ',' O  ',) : continue
                    IPs.append(v) ; count = count + 1
                    mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
                mean[0] /= count ; mean[1] /= count ; mean[2] /= count
                rad = scRad
                irlist.append( CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[index],rad), VecFloat(mean), scRad) )
        return blist, bfoll, numtrials, irlist



def mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1, dum1=None, dum2=None) :
    assert len(blist) == len(numtrials)
    assert len(blist1) == len(numtrials1)
    blistsize = len(blist)
    for bi in range(len(blist1)) :
        blist.append( blist1[bi] )
        numtrials.append( numtrials1[bi] )
        if bfoll1.has_key(bi) :
            bfoll[blistsize+bi] = []
            for bi1 in bfoll1[bi] : bfoll[blistsize + bi].append( blistsize + bi1 )
    for r in rlist1 : rlist.append(r)
    #print "DUMMIES", dum1,
    if dum1!=None and dum2!=None :
        for k,v in dum2.items() : dum1[k] = v
    #print dum1

dummyrn = 1000
def getNextDummyResnum() :
    global dummyrn
    dummyrn = dummyrn-1
    return "%3d" % dummyrn

## make N or C terminal of the chain

def makeSCcentroidRestr(res, pts, resids, scRad, index) :
    IPs = []
    mean, count = [0.,0.,0.], 0
    for k,v in res[index].items() :
        if k in (' N  ',' CA ',' C  ',' O  ',) : continue
        IPs.append(v) ; count = count + 1
        mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
    if count == 0 : return None
    mean[0] /= count ; mean[1] /= count ; mean[2] /= count
    return CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[index],scRad), VecFloat(mean), scRad)

def prepareSConly(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, scRad, scMissInds=[]) :
    blist, rlist = [], []
    for ind in range(startindex,endindex+1) :
        if resns[ind] in ["ALA","GLY"] : continue
        blist.append( makeChiBuilder(ind-1, ind, ind+1, res, resns, resids) )
        IPs, mean, count = [], [0.,0.,0.], 0
        for k,v in res[ind].items() :
            if k in (' N  ',' CA ',' C  ',' O  ',) : continue
            IPs.append(v) ; count = count + 1
            mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
        mean[0] /= count ; mean[1] /= count ; mean[2] /= count
        rlist.append( CentroidPosRestraint(VecInt(IPs), "CentroidPosRestr %s sc %f" % (resids[ind],scRad), VecFloat(mean), scRad) )
    return blist, {}, makeNumtrials(blist), rlist



def addNdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts) :
    ## add dummy residues, this is independent of direction of building
    dummies = {}
    keys = list(res.keys())
    keys.sort()
    ## insert a dummy GLY residue with atoms C,O,CA at the beginning of traced chain for last CO and sidechain
    newkey = keys[0]-1
    pts.append([0,0,0]) ; pts.append([0,0,0]) ; pts.append([0,0,0])
    res[newkey] = { ' C  ': len(pts)-1, ' O  ': len(pts)-2, ' CA ': len(pts)-3 }
    resnums[newkey] = getNextDummyResnum() ; inscodes[newkey] = ' ' ; resns[newkey] = 'GLY' ; chids[newkey] = chid
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    dummies[newkey] = [' C  ',' O  ',' CA ']
    firstindex = newkey
    #print dummies
    return dummies, firstindex
def addCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts) :
    ## add dummy residues, this is independent of direction of building
    dummies = {}
    keys = list(res.keys())
    keys.sort()
    ## insert a dummy GLY residue with atoms N,CA at the end of traced chain for last CO and sidechain
    newkey = keys[len(keys)-1] + 1
    pts.append([0,0,0]) ; pts.append([0,0,0])
    res[newkey] = { ' N  ': len(pts)-2, ' CA ': len(pts)-1 }
    resnums[newkey] = getNextDummyResnum() ; inscodes[newkey] = ' ' ; resns[newkey] = 'GLY' ; chids[newkey] = chid
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    dummies[newkey] = [' N  ',' CA ']
    lastindex = newkey
    return dummies, lastindex
def addNCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts) :
    dumN, firstindex = addNdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
    dumC, lastindex = addCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
    for k,v in dumC.items() : dumN[k] = v
    return dumN, firstindex, lastindex

## if sidechain atoms are missing, add them and change res,pts accordingly.
## return res indices where sidechain atoms were missing
def incompleteSCcorrection(res, resns, pts) :
    scMissInds = []
    for ind in res.keys() :
        if not resAtoms.has_key(resns[ind]) : continue
        scmiss = None
        for an in resAtoms[resns[ind]] :
            if an in [' N  ',' CA ',' C  ',' O  '] : continue
            if an in res[ind].keys() : continue
            scmiss = 1
            res[ind][an] = len(pts)
            pts.append( (0.,0.,0.) )
        if scmiss : scMissInds.append(ind)
    return scMissInds

def makeRNAbuilder(res, pts, resns, resids, ri, ri1) :
    IPs, OPs = [], []
    IPs.append( res[ri][' C5*'] ) ; IPs.append( res[ri][' C4*'] ) ; IPs.append( res[ri][' C3*'] )

    #print resids[ri], resids[ri1]
    OPs.append( res[ri][' O3*'] )
    OPs.append( res[ri1][' P  '] )
    OPs.append( res[ri1][' O5*'] )
    OPs.append( res[ri1][' C5*'] )
    OPs.append( res[ri1][' C4*'] )
    OPs.append( res[ri1][' C3*'] )
    OPs.append( res[ri1][' O1P'] )
    OPs.append( res[ri1][' O2P'] )

    OPs.append( res[ri][' O4*'] ) ; OPs.append( res[ri][' C1*'] ) ; OPs.append( res[ri][' C2*'] ) ## sugar
    dna = 0
    if ' O2*' in res[ri].keys() : dna = 0
    if dna==0 : OPs.append( res[ri][' O2*'] ) ## extra oxygen for RNA

    if resns[ri] in [ '  T', '  C', '  U', ' DT', ' DC', ' DU' ] : ## TCU bases pyrimidines
        OPs.append(res[ri][' N1 ']) ; OPs.append(res[ri][' C2 ']) ; OPs.append(res[ri][' O2 ']); OPs.append(res[ri][' N3 '])
        OPs.append(res[ri][' C4 ']) ; OPs.append(res[ri][' C5 ']) ; OPs.append(res[ri][' C6 '])
        if resns[ri][2] == 'T' :
            OPs.append(res[ri][' O4 ']) ; OPs.append(res[ri][' C5M'])
        elif resns[ri][2] == 'C' :
            OPs.append(res[ri][' N4 '])
        elif resns[ri][2] == 'U' :
            OPs.append(res[ri][' O4 '])
        else : assert 0
    elif resns[ri] in [ '  A', '  G', ' DA', ' DG' ] : ## AG bases purines
        OPs.append( res[ri][' N1 ']) ; OPs.append( res[ri][' C2 ']) ; OPs.append( res[ri][' N3 '])
        OPs.append( res[ri][' C4 ']) ; OPs.append( res[ri][' C5 ']) ; OPs.append( res[ri][' C6 '])
        OPs.append( res[ri][' N7 ']) ; OPs.append( res[ri][' C8 ']) ; OPs.append( res[ri][' N9 '])
        if resns[ri][2] == 'A' :
            OPs.append( res[ri][' N6 '] )
        elif resns[ri][2] == 'G' :
            OPs.append( res[ri][' O6 '] ) ; OPs.append( res[ri][' N2 '] )
        else : assert 0
    else : assert 0
    return NAsuiteBuilder(VecInt(IPs), VecInt(OPs), consts, "NAsuiteBuilder %s" % resids[ri], nas, resns[ri][2], dna, "C3'-endo", "C3'-endo")

def appendDummyEndResRNA(res, pts, resids, resns, resnums, chids, inscodes, chid, newresnum) :
    keys = list(res.keys()) ; keys.sort()
    newkey = keys[len(keys)-1] + 1
    res[newkey] = {}
    dummies = {newkey:[]}
    for an in [' P  ',' O5*',' C5*',' C4*',' C3*',' O1P',' O2P'] :
        pts.append([0,0,0])
        res[newkey][an] = len(pts)-1
        dummies[newkey].append(an)
    resns[newkey] = '  A'
    resnums[newkey] = newresnum
    chids[newkey] = chid
    inscodes[newkey] = ' '
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    return dummies

componentBuilders = [] #just a hack to stop crashing when component builders in BuilderGroup are called
def prepareRNAchain(chid, res, resids, resnums, resns, chids, inscodes, pts, rnaRad=2.) :
    blist, numtrials, irlist = [], [], []

    maxEndIndex = res.keys()[ len(res.keys())-1 ]
    startindex, endindex = None, res.keys()[ len(res.keys())-1 ]
    for ind in res.keys() :
        if startindex == None and chids[ind] == chid : startindex = ind
        if chids[ind] != chid and startindex != None and endindex == maxEndIndex :
            endindex = ind-1 ; break
    assert startindex < endindex

    ## add additional partial end residue (with atoms P,O5*,C5*,C4*,C3*,O1P,O2P) to each rna chain so that RNASuiteBuilder works
    dummies = appendDummyEndResRNA(res, pts, resids, resns, resnums, chids, inscodes, chid, getNextDummyResnum())

    ## bootstrapping of RNA chain is done by finding C5*,C4* positions within 1A of given ones and translating P,O5*,C3* along.
    ## then C3* is rotated around C5*,C4* until it is with 1A of its given position, same rotation is given to P,O5* too
    ## after that normal rna suite builders are constructed
    IPs, OPs = [], []
    for an in [' C5*',' C4*',' C3*'] : print an, resids[startindex]; print res[startindex].keys(); OPs.append(res[startindex][an]) # compulsorily have to present
    for an in [' P  ',' O1P',' O2P',' O5*'] : # optional
        if an in res[startindex].keys() : OPs.append(res[startindex][an])
    b1 = TransRotator(VecInt(IPs), VecInt(OPs), None, "RNAbootstrap-1", VecFloat(pts[res[startindex][' C5*']]), rnaRad, VecFloat(pts[res[startindex][' C4*']]), rnaRad)
    IPs, OPs = [res[startindex][' C5*'], res[startindex][' C4*']], [res[startindex][' C3*']]
    for an in [' P  ',' O5*',' O1P',' O2P'] : # optional
        if an in res[startindex].keys() : OPs.append(res[startindex][an])
    b2 = Rotator(VecInt(IPs), VecInt(OPs), None, "RNAbootstrap-2", -180,180,10)
    blist.append( BuilderGroup(VecBuilder([b1,b2]),"RNAbootstrap") )
    componentBuilders.append(b1) ; componentBuilders.append(b2)
    numtrials.append(100)
    irlist.append( SphPosRestr(VecInt([res[startindex][' C3*']]), "SphPosRestraint %s C3*" % resids[startindex], VecFloat(pts[res[startindex][' C3*']]), 0, rnaRad) )
    for index in range(startindex,endindex+1) :
        #print resids[index]
        ri0,ri1 = index,index+1
        if ri0 == endindex : ri1 = dummies.keys()[0]
        #print ri0, ri1, endindex, res[ri0], res[ri1]
        blist.append( makeRNAbuilder(res, pts, resns, resids, ri0, ri1) )
        irlist.append( SphPosRestr(VecInt([res[ri0][' C3*']]), "SphPosRestraint %s C3*" % resids[ri0], VecFloat(pts[res[ri0][' C3*']]), 0, rnaRad) )
        numtrials.append(100)
    return blist, numtrials, {}, irlist, dummies #no builder-followers here

def removeSC(res, pts, forIndices) :
    newpts = []
    for index in forIndices :
        for an in res[index].keys() :
            if not an in [ ' N  ', ' CA ', ' C  ', ' O  ' ] : del res[index][an]
    for index in forIndices :
        for an,pi in res[index].items() :
            newpts.append( pts[pi] )
            res[index][an] = len(newpts)-1
    return res, newpts

def makeCApropRestraints(CAind, proprange, limDOWN, limUP, caRad, res, resids, pts) :
    rlist = []
    for index in range(CAind-proprange, CAind+proprange+1) :
        if index < limDOWN or index > limUP : continue
        caRR = 3.81 * abs(CAind-index) + caRad
        rlist.append( makeSphPosRestr(CAind, index, 0, caRR, res, resids, pts) )
    return rlist

## given an N and CO, make N-O distance and C-O-N angle restraint for hbonding
def makeHbondRestraints(ni, oi, res, resids) :
    rlist = []
    ONhbd = (1.5, 3.5) # distance
    CONhbd = (100, 200) # angle limits
    rlist.append( DistanceRestraint(VecInt([res[ni][' N  '],res[oi][' O  ']]), "betaHbond-NO [%s] [%s]" % (resids[ni],resids[oi]), ONhbd[0], ONhbd[1]) )
    rlist.append( AngleRestraint(VecInt([res[ni][' N  '],res[oi][' O  '],res[oi][' C  ']]), "betaHbond-NOC [%s] [%s]" % (resids[ni],resids[oi]), CONhbd[0], CONhbd[1]) )
    return rlist

## can be extended to pi helix also later
## puts hbonding restraints on O(i),N(i+4)
## CA(i),CA(i+2) (5.1 - 5.8, in place of alpha-specific phi-psi sampling)
## alpha and beta specific sampling certainly a TODO
def prepareAlphaHelix(helixtype, res, resids, resnums, resns, chids, inscodes, pts, caRad, start, end, bootstrap, mconly, informed) :
    useSCfollowers, blist, bfoll, numtrials, rlist = 1, [], {}, [], []
    if bootstrap :
        b1, b2 = bootstrap
        IPs, OPs = [], [ res[b1][' CA '], res[b1][' C  '], res[b1][' O  '], res[b2][' N  '], res[b2][' CA '], ]
        blist.append( NanchorBuilder( VecInt(IPs), VecInt(OPs), consts, "NanchorBuilder [%s][%s]"%(resids[b1],resids[b2]),
            VecFloat(pts[res[b1][' CA ']]), caRad, VecFloat(pts[res[b2][' CA ']]), caRad ) )
    if start < end :
        for ind in range(start,end+1) :
            bgrp = makeBgrp(ind-2, ind-1, ind, pts, res, resns, resids, mconly, 1, informed)
            blistsize = len(blist)
            blist = blist + bgrp
            if len(bgrp) > 1 and useSCfollowers : bfoll[blistsize] = range(blistsize+1,len(blist))
            if helixtype == "4_13" and ind <= end-4 : rlist += makeHbondRestraints(ind+4, ind, res, resids)
            if helixtype == "3_10" and ind <= end-3 : rlist += makeHbondRestraints(ind+3, ind, res, resids)
            if ind <= end-2 : rlist.append( DistanceRestraint(VecInt([res[ind][' CA '],res[ind+2][' CA ']]), "helix-CA-i-i+2 dist-restraint", 5.0, 6.0) ) # 5.1-5.8
            rlist.append( makeSphPosRestr(ind, ind, 0, caRad, res, resids, pts) )
    else :
        for ind in range(start,end-1,-1) : # bcz end is start and start is end
            #print ind, ind+1, ind+2
            bgrp = makeBgrp(ind, ind+1, ind+2, pts, res, resns, resids, mconly, 0, informed)
            blistsize = len(blist)
            blist = blist + bgrp
            if len(bgrp) > 1 and useSCfollowers : bfoll[blistsize] = range(blistsize+1,len(blist))
            if helixtype == "4_13" and ind >= start+4 : rlist += makeHbondRestraints(ind, ind-4, res, resids)
            if helixtype == "3_10" and ind >= start+3 : rlist += makeHbondRestraints(ind, ind-3, res, resids)
            if ind >= start+2 : rlist.append( DistanceRestraint(VecInt([res[ind][' CA '],res[ind-2][' CA ']]), "helix-CA-i-i+2 dist-restraint", 5.0, 6.0) ) # 5.1-5.8
            rlist.append( makeSphPosRestr(ind, ind, 0, caRad, res, resids, pts) )
    ## positional restraints propagated
    if start < end : proprange = range(start,end+1)
    else : proprange = range(end,start+1)
    for index in proprange :
        rlist = rlist + makeCApropRestraints(index, 3, proprange[0], proprange[len(proprange)-1], caRad, res, resids, pts)
    numtrials = makeNumtrials(blist)
    return blist, bfoll, numtrials, rlist

## for parallel beta sheets, Ni - Oj-1, Oi - Nj+1, every alternate residue
## for antiparallel beta sheets, Ni - Oj, Oi - Nj, every alternate residue
## corrInds[si] is a pair of corresponding indices in strand si and si+1
def prepareBetaSheet(res, resids, resnums, resns, chids, inscodes, pts, caRad, starts, ends, sheetBS, directions, corrInds, mconly, informed) :
    for d in directions : assert d == "fwd" or d == "bkwd"
    for si in range(len(starts)) : assert starts[si] < ends[si]
    # make steps and ladder
    ladder = {} # ladder[residue_index] gives step index
    steps = [] # steps[step_index][strand_index] gives peer in that strand
    for si in range(len(starts)-1) : # need not access last strand
        for index in range(starts[si],ends[si]+1) :
            if not index in ladder.keys() :
                steps.append( [None] * len(starts) )
                ladder[index] = len(steps)-1
                steps[ladder[index]][si] = index
            if directions[si] == directions[si+1] : # parallel
                cind = corrInds[si][1] + index - corrInds[si][0]
            else : # antiparallel
                cind = corrInds[si][1] - index + corrInds[si][0]
            if not cind in range(starts[si+1],ends[si+1]+1) : continue
            if not cind in ladder.keys() : ladder[cind] = ladder[index]
            steps[ladder[cind]][si+1] = cind
    for st in steps :
        steplen = 0
        for sti in st :
            if sti != None : steplen += 1
        if steplen < 2 : assert None
    # order the steps
    stepsOrder = range(len(steps))
    for ii in range(len(stepsOrder)) :
        for jj in range(ii+1,len(stepsOrder)) :
            i,j = stepsOrder[ii], stepsOrder[jj]
            jLarger = 0 ; ncomparisons = 0
            for si in range(len(steps[i])) :
                if steps[i][si] == None or steps[j][si] == None : continue
                if directions[si] == "fwd" :
                    if steps[i][si] < steps[j][si] :
                        assert jLarger in (1,0) ; jLarger = 1
                    else :
                        assert jLarger in (None,0) ; jLarger = None
                else :
                    if steps[i][si] > steps[j][si] :
                        assert jLarger in (1,0) ; jLarger = 1
                    else :
                        assert jLarger in (None,0) ; jLarger = None
                ncomparisons += 1
            if ncomparisons == 0 : jLarger, ncomparisons = None, 1
            assert ncomparisons > 0 and jLarger in (1,None)
            if jLarger :
                temp = stepsOrder[ii]
                stepsOrder[ii] = stepsOrder[jj]
                stepsOrder[jj] = temp
    newsteps = []
    for i in stepsOrder : newsteps.insert(0, steps[i])
    steps = newsteps
    for st in steps :
        print "ORDSTEP",
        for sti in st :
            if sti != None : print "[%s]" % resids[sti],
            else : print "----None---",
        print ''
    # make builders according to steps
    # for a strand-residue in a step, if strand is fwd, make a fwd builder for that CA and bkwd otherwise
    blist, bfoll, useSCfollowers = [], {}, 1
    strandBootstrapped = [None] * len(steps[0])
    for step in steps :
        for si in range(len(step)) :
            if step[si] == None : continue
            if sheetBS[si] and not strandBootstrapped[si] :
                strandBootstrapped[si], b1,b2 = 1, sheetBS[si][0], sheetBS[si][1]
                IPs, OPs = [], [ res[b1][' CA '], res[b1][' C  '], res[b1][' O  '], res[b2][' N  '], res[b2][' CA '], ]
                blist.append( NanchorBuilder( VecInt(IPs), VecInt(OPs), consts, "NanchorBuilder [%s][%s]"%(resids[b1],resids[b2]),
                    VecFloat(pts[res[b1][' CA ']]), caRad, VecFloat(pts[res[b2][' CA ']]), caRad ) )
            if directions[si] == "fwd" : bgrp = makeBgrp(step[si]-2, step[si]-1, step[si], pts, res, resns, resids, mconly, 1, informed)
            else : bgrp = makeBgrp(step[si], step[si]+1, step[si]+2, pts, res, resns, resids, mconly, 0, informed)
            blistsize = len(blist)
            blist = blist + bgrp
            if len(bgrp) > 1 and useSCfollowers : bfoll[blistsize] = range(blistsize+1,len(blist))
    builtAtomInds = []
    for b in blist :
        bop = b.getOP()
        for i in range(bop.size()) : builtAtomInds.append(bop[i])
    ## make hbond restraints between strands si , si+1. make the restraint only if this function is returning builders for both atoms
    builtAtomInds = builtAtomInds
    rlist = []
    for si in range(len(starts)-1) :
        for index in range(starts[si],ends[si]+1) :
            if (index - corrInds[si][0]) % 2 == 1 : continue ## alternate residues only
            if directions[si] == directions[si+1] : # parallel
                cind = corrInds[si][1] + index - corrInds[si][0]
                if res[index][' N  '] in builtAtomInds and res[cind-1][' O  '] in builtAtomInds :
                    rlist += makeHbondRestraints(index, cind-1, res, resids)
                if res[index][' O  '] in builtAtomInds and res[cind+1][' N  '] in builtAtomInds :
                    rlist += makeHbondRestraints(cind+1, index, res, resids)
            else : # antiparallel
                cind = corrInds[si][1] + corrInds[si][0] - index
                if res[index][' N  '] in builtAtomInds and res[cind][' O  '] in builtAtomInds :
                    rlist += makeHbondRestraints(index, cind, res, resids)
                if res[index][' O  '] in builtAtomInds and res[cind][' N  '] in builtAtomInds :
                    rlist += makeHbondRestraints(cind, index, res, resids)

    ## positional restraints and their propagation
    for si in range(len(starts)) :
        assert starts[si] < ends[si]
        for index in range(starts[si], ends[si]+1) :
            rlist += makeCApropRestraints(index, 3, starts[si], ends[si], caRad, res, resids, pts)

    numtrials = makeNumtrials(blist)
    return blist, bfoll, numtrials, rlist
