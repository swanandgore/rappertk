from pdbr import protein
import prot2res
from peptidebuild import ModelRenderer
prot = protein("complete.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

reskeys = res.keys()
reskeys.sort()

closeCA = {}
for ri in reskeys :
    for rk in reskeys :
        if rk <= ri : continue
        pi = pts[ res[ri][" CA "] ]
        pk = pts[ res[rk][" CA "] ]
        dist = 0.
        for i in range(3) : dist += (pi[i]-pk[i])*(pi[i]-pk[i])
        if dist < 3.*3. :
            if not ri in closeCA.keys() : closeCA[ri] = []
            if not rk in closeCA.keys() : closeCA[rk] = []
            closeCA[ri].append(rk)
            closeCA[rk].append(ri)


for ri in reskeys :
    print "[%s]" % resids[ri],
    if not ri in closeCA.keys() or len( closeCA[ri] ) == 0 : print 0
    else : print 1
