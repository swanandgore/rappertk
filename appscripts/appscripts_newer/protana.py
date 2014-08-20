import sys, os
from pdbr import protein
import prot2res
from geometry import calcDihed, VecFloat
from samplers import PhipsiSampler

cbDatapath = os.environ["RTKROOT"] + "/data/"
class PhipsiSamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = PhipsiSampler( cbDatapath + "/PhipsiWeightedProp/ps%s" % resn )
        return self.samplers[resn]

def calcDihed1(a, b, c, d) :
    return calcDihed(VecFloat(a), VecFloat(b), VecFloat(c), VecFloat(d))

if __name__ == "__main__" :
    prot = protein(sys.argv[1], read_hydrogens=0, read_waters=0, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    psp = PhipsiSamplerProvider()

    keys = res.keys() ; keys.sort()
    startindex, endindex = -9999, -9999
    for k in keys :
        if chids[k] == sys.argv[2]  and startindex == -9999 : startindex = k
        if startindex != -9999 and endindex == -9999 and chids[k] != sys.argv[2] : endindex = k
    if startindex != -9999 and endindex == -9999 : endindex = keys[len(keys)-1]
    print startindex, endindex

    for k in range(startindex+1, endindex) :
        ps = psp.get(resns[k])
        omega = calcDihed1( pts[res[k-1][' CA ']], pts[res[k-1][' C  ']], pts[res[k][' N  ']], pts[res[k][' CA ']] )
        phi = calcDihed1( pts[res[k-1][' C  ']], pts[res[k][' N  ']], pts[res[k][' CA ']], pts[res[k][' C  ']] )
        psi = calcDihed1( pts[res[k][' N  ']], pts[res[k][' CA ']], pts[res[k][' C  ']], pts[res[k+1][' N  ']] )
        print resids[k], "\t%6.1f\t%6.1f\t%6.1f\t" % (omega, phi, psi), ps.findProb(phi, psi) * 5184
