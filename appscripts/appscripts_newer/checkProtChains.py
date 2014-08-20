from pdbr import protein
import prot2res
from data import resAtoms
from geometry import calcDist, VecFloat

def check(res, resids, resnums, resns, chids, inscodes, pts, cadistCutoff=4.) :
    chains = set(chids.values())

    mcmiss, scmiss, chbr = set(), set(), set()
    for chid in chains :
        keys = res.keys() ; keys.sort()
        start, end = None, keys[len(keys)-1]
        print "CHAIN", chid,
        for k in keys :
            if chids[k] == chid and not start : start = k
            elif chids[k] != chid and start and end == keys[len(keys)-1] : end = k-1
        print "------------------ [%s] to [%s] ----------------" % (resids[start], resids[end]) ; assert start and end
        for k in range(start,end+1) :
            assert chids[k] == chid
            mcok, scok = [], []
            for an in resAtoms[resns[k]] :
                if not res[k].has_key(an) :
                    if an in [' N  ',' CA ',' C  ',' O  '] : mcok.append(an)
                    else : scok.append(an)
            if len(mcok) > 0 : mcmiss.add(k) ; print "Missing mainchain atoms in", resids[k], mcok
            if len(scok) > 0 : scmiss.add(k) ; print "Missing sidechain atoms in", resids[k], scok
            if k == end : continue
            if not res[k].has_key(' CA ') or not res[k+1].has_key(' CA ') :
                chbr.add(k) ; print "Cant determine chainbreak due to missing CA atoms", resids[k], resids[k+1]
            if calcDist( VecFloat(pts[res[k][' CA ']]), VecFloat(pts[res[k+1][' CA ']]) ) > cadistCutoff :
                chbr.add(k) ; print "Chainbreak between", resids[k], resids[k+1]
    return mcmiss, scmiss, chbr

if __name__ == "__main__" :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--cacaCutoff", action='store', type='float', dest='pdbfile', help='min dist between adjacent CA to detect a chain-break', default=4.)
    (options, args) = parser.parse_args()
    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    check(res, resids, resnums, resns, chids, inscodes, pts, options.cacaCutoff)
