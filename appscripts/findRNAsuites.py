from samplers import NAsuiteSampler, SuiteVal
from pdbr import protein, line2crd, line2resnum, line2resn, line2crd, line2resid, line2chid, line2resic, line2atomname
import os, string, sys
from geometry import calcDihed, VecFloat
import prot2res

cbDatapath = os.environ["RTKROOT"] + "/data/"

def calcDihed1(a, b, c, d) :
    return int( round (calcDihed(VecFloat(a), VecFloat(b), VecFloat(c), VecFloat(d))) )

def main() :
    print sys.argv
    prot = protein(sys.argv[1], read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resnames, chids, resics, crds = prot2res.readProtRes(prot)

    
    nas = NAsuiteSampler( "%s/rnaSuite" % cbDatapath , 0 )
    for i in range(len(res)) :
        try :
            delta0 = calcDihed1( crds[res[i][' C5*']], crds[res[i][' C4*']], crds[res[i][' C3*']], crds[res[i][' O3*']] )
            epsilon = calcDihed1( crds[res[i][' C4*']], crds[res[i][' C3*']], crds[res[i][' O3*']], crds[res[i+1][' P  ']] )
            zeta = calcDihed1( crds[res[i][' C3*']], crds[res[i][' O3*']], crds[res[i+1][' P  ']], crds[res[i+1][' O5*']] )
            alpha = calcDihed1( crds[res[i][' O3*']], crds[res[i+1][' P  ']], crds[res[i+1][' O5*']], crds[res[i+1][' C5*']] )
            beta = calcDihed1( crds[res[i+1][' P  ']], crds[res[i+1][' O5*']], crds[res[i+1][' C5*']], crds[res[i+1][' C4*']] )
            gamma = calcDihed1( crds[res[i+1][' O5*']], crds[res[i+1][' C5*']], crds[res[i+1][' C4*']], crds[res[i+1][' C3*']] )
            delta1 = calcDihed1( crds[res[i+1][' C5*']], crds[res[i+1][' C4*']], crds[res[i+1][' C3*']], crds[res[i+1][' O3*']] )
        except KeyError : continue

        print "[%s]\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (resids[i], delta0, epsilon, zeta, alpha, beta, gamma, delta1)
        sv = SuiteVal(delta0, epsilon, zeta, alpha, beta, gamma, delta1)
        #sv = SuiteVal(84, -150, 60, 65, 180, 55, 84)
        nas.closestSuiteKey(sv)
 

if __name__ == "__main__" : main()
