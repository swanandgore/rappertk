
from pdbr import protein, line2crd, line2resnum, line2resn, line2crd, line2resid, line2chid, line2resic, line2atomname

def readProtRes(prot) :
    res, resids, resnums, resnames, chids, resics, crds = {}, {}, {}, {}, {}, {}, []
    for ri in range(len(prot.reslines)) :
        start, stop = prot.reslines[ri]
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
    return res, resids, resnums, resnames, chids, resics, crds


if __name__ == "__main__" :
    prot = protein("/tmp/2db3.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resnames, chids, resics, crds = readProtRes(prot)
    assert res.keys() == resids.keys()
    assert res.keys() == resnums.keys()
    assert res.keys() == resnames.keys()
    assert res.keys() == chids.keys()
    assert res.keys() == resics.keys()
    print res
