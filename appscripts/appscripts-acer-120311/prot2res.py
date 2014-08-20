
from pdbr import protein, line2crd, line2resnum, line2resn, line2crd, line2resid, line2chid, line2resic, line2atomname, isAAres , line2atomnumber , line2bfac , isHetAtomLine

from data import resAtoms

## find the AA residues with incomplete or no sidechain, and add missing coordintaes set to 0,0,0
def addSC(res, resns, pts) :
    for ri in res.keys() :
        if not isAAres(resns[ri]) : continue
        missingSCatoms = list(resAtoms[resns[ri]])
        for an in [' N  ',' CA ',' C  ',' O  '] : missingSCatoms.remove(an)
        for an in res[ri].keys() :
            if an in missingSCatoms : missingSCatoms.remove(an)
        for an in missingSCatoms :
            res[ri][an] = len(pts) ; pts.append([0.,0.,0.])

## wherever resn has changed in the AA, remove old sc and add new sc.
## CB coordinate is also erased during this!
def changeSeq(res, resns, pts, newresns) :
    print res
    newpts = []
    for ri in res.keys() :
        if resns[ri] == newresns[ri] or not isAAres(resns[ri]) :
            for an,ai in res[ri].items() : newtpts.append( pts[ai] )
            continue
        scatomsOld = set(resAtoms[resns[ri]])
        scatomsNew = set(resAtoms[newresns[ri]])
        for an in scatomsNew.intersection(scatomsOld) :
            newpts.append(pts[res[ri][an]]) ; res[ri][an] = len(newpts)-1
        for an in scatomsNew.difference(scatomsOld) :
            newpts.append([0.,0.,0.]) ; res[ri][an] = len(newpts)
    for i in range(len(pts)) : pts.pop()
    for p in newpts : pts.append(p)


def readProtRes(prot) :
    res, resids, resnums, resnames, chids, resics, crds = {}, {}, {}, {}, {}, {}, []
    
    
    for ri in range(len(prot.reslines)) :
        start, stop = prot.reslines[ri]
        if start == 0 and  stop == 0:
            print "No ATOM records found in file"
            continue
        resids[ri] = line2resid( prot.atomlines[start] )
        resnames[ri] = line2resn( prot.atomlines[start] )
        resnums[ri] = line2resnum( prot.atomlines[start] )
        chids[ri] = line2chid( prot.atomlines[start] )
        resics[ri] = line2resic( prot.atomlines[start] )
        res[ri] = {}
        for ai in range(start, stop) :
            crd = line2crd( prot.atomlines[ai] )
            aname = line2atomname( prot.atomlines[ai] )
            crds.append(crd)
            res[ri][aname] = len(crds)-1
    if len(res.keys()) == 0  and start ==0 and stop ==0: 
        print "Coordinate file corrupted, pleas check file format. Unable to read coordinae file"
        import sys ; sys.exit()
    return res, resids, resnums, resnames, chids, resics, crds


def readProtRes3(prot) :
    res, resids, resnums, resnames, chids, resics, crds, hetids = {}, {}, {}, {}, {}, {}, [] , {}
    for ri in range(len(prot.reslines)) :
#        print prot.reslines[ri]
        start, stop = prot.reslines[ri]
#        print prot.atomlines[start]
        resids[ri] = line2resid( prot.atomlines[start] )
        resnames[ri] = line2resn( prot.atomlines[start] )
        resnums[ri] = line2resnum( prot.atomlines[start] )
        chids[ri] = line2chid( prot.atomlines[start] )
        resics[ri] = line2resic( prot.atomlines[start] )
        if isHetAtomLine( prot.atomlines[start] ) :
            hetids[ri] = 1
        else :

            hetids[ri] = 0
        res[ri] = {}

        for ai in range(start, stop) :
            crd = line2crd( prot.atomlines[ai] )
            aname = line2atomname( prot.atomlines[ai] )
            crds.append(crd)
            res[ri][aname] = len(crds)-1
    return res, resids, resnums, resnames, chids, resics, crds , hetids


def readProtRes2(prot) :
    res, resids, resnums, resnames, chids, resics, crds,atomids, bfac = {}, {}, {}, {}, {}, {}, [], {} , {} 
    for ri in range(len(prot.reslines)) :
        start, stop = prot.reslines[ri]
        resids[ri] = line2resid( prot.atomlines[start] )
        resnames[ri] = line2resn( prot.atomlines[start] )
        resnums[ri] = line2resnum( prot.atomlines[start] )
        chids[ri] = line2chid( prot.atomlines[start] )
        resics[ri] = line2resic( prot.atomlines[start] )
        res[ri] = {}
        atomids[ri] = {}
        bfac[ri] = {}
        for ai in range(start, stop) :
            crd = line2crd( prot.atomlines[ai] )
            aname = line2atomname( prot.atomlines[ai] )
            crds.append(crd)
            res[ri][aname] = len(crds)-1
            atomids[ri][aname] = line2atomnumber( prot.atomlines[ai] )
            bfac[ri][aname] = line2bfac( prot.atomlines[ai] )
    return res, resids, resnums, resnames, chids, resics, crds , atomids , bfac


if __name__ == "__main__" :
    res = {0:{" N  ":0, " CA ":1, " C  ":2, " O  ":3, " CB ":4}}
    resns = {0:"ALA"}; newresns = {0:"GLU"}
    pts = []
    for i in range(5) : pts.append( [0.,0.,0.] )
    changeSeq(res, resns, pts, newresns)
    print res
    print pts
    sys.exit(0)
    prot = protein("/tmp/2db3.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resnames, chids, resics, crds = readProtRes(prot)
    assert res.keys() == resids.keys()
    assert res.keys() == resnums.keys()
    assert res.keys() == resnames.keys()
    assert res.keys() == chids.keys()
    assert res.keys() == resics.keys()
    print res