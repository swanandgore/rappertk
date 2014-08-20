from pdbr import protein
import prot2res, string
from peptidebuild import ModelRenderer

if __name__ == "__main__" :
    import sys ; print sys.argv
    halfWidth = string.atof(sys.argv[3]) / 2.
    prot = protein(sys.argv[1], read_hydrogens=1, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    ptis = []
    for ind in res.keys() :
        ocrd = pts[res[ind]['OH2 ']]
        print resids[ind], ocrd
        if -1*halfWidth > ocrd[0] or halfWidth < ocrd[0] : continue
        if -1*halfWidth > ocrd[1] or halfWidth < ocrd[1] : continue
        if -1*halfWidth > ocrd[2] or halfWidth < ocrd[2] : continue
        ptis = ptis + res[ind].values()
        print "IN", resids[ind], len(ptis)
    print ptis
    mr = ModelRenderer(res, resns, chids, resnums, inscodes, sys.argv[2])
    mr.render(pts, ptis)
