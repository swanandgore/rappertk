import re
from pdbr import protein, isAAres
from data import consts
import prot2res
import data
from data import resAtoms, vdwr
from prepareChain import PrepareChain, printResAtoms
from builders import VecInt, VecVecInt, VecVecFloat, VecFloat, CBbuilder
from scep import SCenergyProvider, findBC, VecBuilder
from geometry import CAtraceGH, Grid
from evalCAtrace import findChis

def getAllAtomCrds(res,pts,ri,resn) :
    retpts = []
    for an in resAtoms[resn] :
        retpts.append( pts[ res[ri][an] ] )
    return retpts

class SCplacement :
    def __init__(s, pdbfile, scReduction, outpdb, dotfile, withDEE, mtzfn, fplabel, fclabel, philabel, maptype, esmin, esmax, addSC,
                    considerGivenRotamers, residsToChange, prepC) :
        s.mt = {"2F1-F2":0, "F1":1}
        s.makeConfidentDecisions = None
        s.pdbfile, s.scReduction, s.outpdb, s.dotfile, s.withDEE, s.mtzfn, s.fplabel, s.fclabel, s.philabel, s.maptype, s.addSC, s.residsToChange \
             = pdbfile, scReduction, outpdb, dotfile, withDEE, mtzfn, fplabel, fclabel, philabel, maptype, addSC, residsToChange
        s.esmin, s.esmax = esmin, esmax
        s.considerGivenRotamers = considerGivenRotamers
        if prepC == None : s.prepC = PrepareChain("PRL")
        else : s.prepC = prepC

    def run(s) :
        prot = protein(s.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    
        # addSC == 1 means that sidechain in AA residues have to be added.
        # add those atoms in our data and set coordinates to 0,0,0. build CB where reqd
        if s.addSC == 1 :
            prot2res.addSC(res, resns, pts)
            for ri in res.keys() :
                if not isAAres(resns[ri]) or resns[ri] == "GLY" : continue
                newpts, ipinds, opinds = [], VecInt([0,1,2]), VecInt([3,])
                newpts.append( pts[ res[ri][' N  '] ] )
                newpts.append( pts[ res[ri][' C  '] ] )
                newpts.append( pts[ res[ri][' CA '] ] )
                newpts.append( pts[ res[ri][' CB '] ] )
                newpts = VecVecFloat(newpts)
                CBbuilder(ipinds, opinds, consts, "CB builder", resns[ri]).build (newpts)
                pts[res[ri][' CB ']] = [ newpts[3][0], newpts[3][1], newpts[3][2] ]
    
        #printResAtoms(res, resids)

        ## are all sidechain heavy-atoms present? and which all are not-sidechain atoms
        knownPositions, scbldrs, scepBtags = [], [], []
        for ri in res.keys() :
            scAssignable = None
            if resAtoms.has_key(resns[ri]) :
                scAssignable = 1
                for an in resAtoms[resns[ri]] :
                    if not res[ri].has_key(an) : scAssignable = None
            if resns[ri] in ["ALA","GLY"] : scAssignable = None
            if s.residsToChange != None and not resids[ri] in s.residsToChange : scAssignable = None
            for an,pi in res[ri].items() :
                if scAssignable and not an in [' N  ',' CA ',' C  ',' O  ',' CB '] : continue
                knownPositions.append(pi)
            if not scAssignable : continue
            print "[%s]" % resids[ri], "will be sc-assigned"
            ri0, ri1 = ri-1, ri+1
            if not ri0 in resns.keys() or not resns[ri0] in resAtoms.keys() : ri0 = ri1
            if not ri1 in resns.keys() or not resns[ri1] in resAtoms.keys() : ri1 = ri0
            addsample = None
            if s.considerGivenRotamers != None and s.prepC.prl : addsample = [0.1, findChis(res,pts,ri,resns[ri])]
            elif s.considerGivenRotamers != None and s.prepC.scl : addsample = [0.1, getAllAtomCrds(res,pts,ri,resns[ri])]
            #print addsample
            scbldrs.append( s.prepC.makeChiBuilder(ri0,ri,ri1, res, resns, resids, addsample) )
            scepBtags.append( "[%s]" % resids[ri] )
    
        radii = [0] * len(pts) ## make radii for all atoms including dummies
        for index, val in res.items() :
            for aname, pi in val.items() :
                if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
                else :
                    try : radii[pi] = vdwr['XXX'][aname[1]]
                    except KeyError : radii[pi] = vdwr['XXX']['C']
    
        ai2pos, aiSC = [-999] * len(pts), [-999] * len(pts) ## flag atoms
        for index,val in res.items() :
            for name, ai in val.items() :
                ai2pos[ai] = index
                if name == ' N  ' : aiSC[ai] = 0
                elif name == ' CA ' : aiSC[ai] = 1
                elif name == ' C  ' : aiSC[ai] = 2
                elif name == ' O  ' : aiSC[ai] = 3
                elif name == ' SG ' and resns[index] == 'CYS' : aiSC[ai] = 5 # for disulphide vdw reduction
                elif name == ' CD ' and resns[index] == 'PRO' : aiSC[ai] = 6 # PRO CD shdnt be clash-checked with prev res's C
                else : aiSC[ai] = 4
        ssReduction = data.PROBE_DISULFIDE_OVERLAP_MARGIN
        gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), s.scReduction, ssReduction)
    
        radii = VecFloat(radii)
        pts = VecVecFloat(pts)

        grid, SCep, scbldrs1 = None, None, None

        while 1 :
            grid = Grid(3, pts, radii, gridHelper)
            grid.justAdd(VecInt(knownPositions))
            scbldrs1 = VecBuilder(scbldrs)
            SCep = SCenergyProvider(scbldrs1, pts, grid)

            if re.compile("\.mtz$").search(s.mtzfn) : SCep.addMTZinfo(s.mtzfn, s.fplabel, s.fclabel, s.philabel, s.mt[s.maptype], s.esmin, s.esmax)
            elif re.compile("\.map$").search(s.mtzfn) : SCep.addMAPinfo(s.mtzfn, s.esmin, s.esmax)
            else : print "Cant decide whether this is map or mtz", mtzfn ; assert None
            givenE = SCep.scoreGivenRotamers()
            print "Egiven", givenE ;# sys.exit(0)
            badbldris = []
            for i in range(len(scbldrs)) :
                if SCep.numRot(i) == 0 : badbldris.append(i)
            if len(badbldris) > 0 :
                newbldrs = []
                for i in range(len(scbldrs)) :
                    if not i in badbldris :
                        newbldrs.append( scbldrs[i] )
                    else :
                        bop = scbldrs[i].getOP()
                        for ai in range(bop.size()) : knownPositions.append( bop[ai] )
                scbldrs = newbldrs ; continue
            for i in range(len(scbldrs)) :
                print i, "[%s]" % scepBtags[i], SCep.numRot(i),
                for k in range(SCep.numRot(i)) : print "%6.3f"%SCep.selfEn(i,k),
                print ''
            break

        #import sys; sys.exit(0)
        ## find correlated sc-pairs except when members are far-off
        sciPairs = []
        for i in range(len(scbldrs)) :
            for k in range(i+1,len(scbldrs)) :
                #if SCep.interCBdist(i,k) > SCep.getMaxCBDist(i) + SCep.getMaxCBDist(k) : continue
                if SCep.interact(i,k) == 0 : continue
                #print len(sciPairs), i, k #,  SCep.interCBdist(i,k), SCep.getMaxCBDist(i), SCep.getMaxCBDist(k) 
                sciPairs.append( (i,k) )
        print sciPairs
    
        if s.withDEE :
            for akey in range(len(scbldrs)) :
                print "dee", akey, SCep.numRot(akey)
                SCep.deeGoldstein(akey)
        #import sys; sys.exit(0)
    #    sciPairs = []
    #    sciPairs.append( (0, 1) )
    #    sciPairs.append( (1, 2) )
    #    sciPairs.append( (1, 5) )
    #    sciPairs.append( (1, 6) )
    #    sciPairs.append( (2, 3) )
    #    sciPairs.append( (2, 4) )
    #    sciPairs.append( (3, 4) )
    #    sciPairs.append( (3, 7) )
    #    sciPairs.append( (7, 8) )
    #    sciPairs.append( (7, 10) )
    #    sciPairs.append( (8, 9) )
    
    
        sciPairs = VecVecInt(sciPairs)
        components, artipts = VecVecInt(), VecInt()
        findBC(sciPairs, components, artipts)
        
        for k in range(components.size()) :
            print "COMPONENT", k, "size", len(components[k]), "fallsinto", scepBtags[artipts[k]]
            for i in range(len(components[k])) : print scepBtags[components[k][i]]
        #import sys; sys.exit(0)
    
        if s.dotfile :
            grf = open(s.dotfile, 'w')
            print >> grf, "digraph G {"
            for ai in range(artipts.size()) :
                print >> grf, "%d [color=red]" % artipts[ai]
                if ai > 0 :
                    print >> grf, "%d -> %d [color=red]" % (artipts[ai-1], artipts[ai])
            for abi in range(sciPairs.size()) :
                a,b = sciPairs[abi][0], sciPairs[abi][1]
                gstr = "%d -> %d [dir=none]" % (a,b)
                for ci in range(len(components)) :
                    if a in components[ci] and b in components[ci] : gstr += " [label = %d]" % ci
                print >> grf, gstr
            print >> grf, "}"
        #return
        for i in range(artipts.size()) :
            for k in range(i+1) : print k, "size", len(components[k]), components[k], "fallsinto", artipts[k]
            SCep.collapse( components[i], artipts[i] )
    
        finalAssig = VecInt()
        SCep.join(finalAssig)
        SCep.buildAssig(finalAssig)
        from peptidebuild import ModelRenderer
        modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], s.outpdb)
        modelRenderer.render(pts)
    

        evalResult = None
        if evalResult :
            SCep.resetEnergy()
            if s.withDEE :
                for akey in range(len(scbldrs)) :
                    print "dee", akey, SCep.numRot(akey)
                    SCep.deeGoldstein(akey)
    
            allsc = VecInt( range(finalAssig.size()) )
            ejoin = SCep.findEsofar(allsc, allsc.size()-1, finalAssig, 1)
            print "Ejoin", ejoin,
            for i in range(finalAssig.size()) : print finalAssig[i],
            print ''
        return SCep.scoreGivenRotamers()


def callmain() :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--dotfile", action='store', type='string', dest='dotfile', help='graphviz output file', default=None)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw radius in case of sidechain clashchecks. DEFAULT 1 ie no reduction', default=1)
    parser.add_option("--outpdb", action='store', type='str', dest='outpdb', help='file to write the models to, must be set to valid filepath')
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--dee", action='store', type='int', dest='dee', help='use Goldstein DEE, no by default', default=None)
    parser.add_option("--considerGivenRotamers", action='store', type='int', dest='considerGivenRotamers', help='use present rotamers also in addition to rotamer lib', default=None)

    parser.add_option("--mtz", action='store', type='string', dest='mtzfn', help='map, or mtz, recognized by file extension', default=None)
    parser.add_option("--folabel", action='store', type='string', dest='folabel', help='Fo label', default="FP")
    parser.add_option("--fclabel", action='store', type='string', dest='fclabel', help='Fc label', default="FC")
    parser.add_option("--philabel", action='store', type='string', dest='philabel', help='PHIC label', default="PHIC")
    parser.add_option("--esmin", action='store', type='float', dest='esmin', help='min sigma level below which xray density is penalized (0.25)', default=0.25)
    parser.add_option("--esmax", action='store', type='float', dest='esmax', help='max sigma level above which xray density remains flat (2.)', default=2.)

    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)

    SCplacement(options.pdbfile, options.scReduction, options.outpdb, options.dotfile, options.dee,
        options.mtzfn, options.folabel, options.fclabel, options.philabel, options.esmin, options.esmax,
        None, options.considerGivenRotamers, None).run()


if __name__ == "__main__" : callmain()
