from pdbr import protein
import prot2res
from peptidebuild import ModelRenderer
prot = protein("orig1SG1.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

reskeys = res.keys()
reskeys.sort()

for ri in reskeys :
    resnums[ri] = "%4d" % (ri+1)

ModelRenderer(res, resns, chids, resnums, inscodes, [], "orig1SG1.pdb").render(pts)
