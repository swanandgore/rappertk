import math
from geometry import findSuperpositionTransform, TRTtransform, calcDist
from builders import VecVecFloat, VecFloat
from pdbr import protein
import prot2res

def nicestr(p) :
    return "(%8.3f %8.3f %8.3f)" % (p[0],p[1],p[2])

def simpleRMSD() :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='reference PDB file', default=None)
    parser.add_option("--multipdb", action='store', type='string', dest='multipdb', help='many models in this PDB file', default=None)
    parser.add_option("--outpdb", action='store', type='string', dest='outpdb', help='multipdb + RMSD remarks', default=None)
    (options, args) = parser.parse_args()

    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res1, resids1, resnums1, resns1, chids1, inscodes1, pts1 = prot2res.readProtRes(prot)
    pointsFrom = []
    for ri in res1.keys() :
        pointsFrom.append( pts1[ res1[ri][' CA '] ] )

    ofp = open(options.multipdb, 'r') ; multilines = ofp.readlines() ; ofp.close()
    ofp = open(options.outpdb, 'w')
    modelLines = None
    for l in multilines :
        if l[0:5] == "MODEL" : modelLines = [l]
        elif l[0:6] == "ENDMDL" :
            modelLines.append(l)
            prot = protein(modelLines, read_hydrogens=0, read_waters=0, read_hets=0)
            res2, resids2, resnums2, resns2, chids2, inscodes2, pts2 = prot2res.readProtRes(prot)
            pointsOnto = []
            for ri in res1.keys() :
                assert resids1[ri] == resids2[ri]
                pointsOnto.append( pts2[ res2[ri][' CA '] ] )
            assert len(pointsOnto) == len(pointsFrom)
            rmsd, numDiffPts = 0., 0
            for i in range(len(pointsFrom)) :
                ptdist = calcDist( VecFloat(pointsFrom[i]), VecFloat(pointsOnto[i]) )
                if ptdist < 1e-6 : continue
                rmsd += ptdist*ptdist ; numDiffPts += 1
            if numDiffPts == 0 :
                print "RMSD error because all points are too close to corresponding points" ; sys.exit(1)
            rmsd = math.sqrt(rmsd / numDiffPts)
            print "RMSD", rmsd, "numpts", numDiffPts
            for l in modelLines :
                if "RMSD" in l : continue
                if l[0:5] == "MODEL" : ofp.write(l) ; ofp.write("REMARK         RMSD %f over %d atoms\n" % (rmsd,numDiffPts))
                else : ofp.write(l)
        elif modelLines != None : modelLines.append(l)
    ofp.close()

def superMulti() :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='reference PDB file', default=None)
    parser.add_option("--multipdb", action='store', type='string', dest='multipdb', help='many models in this PDB file', default=None)
    parser.add_option("--outpdb", action='store', type='string', dest='outpdb', help='multipdb + RMSD remarks', default=None)
    (options, args) = parser.parse_args()

    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res1, resids1, resnums1, resns1, chids1, inscodes1, pts1 = prot2res.readProtRes(prot)
    pointsFrom = []
    for ri in res1.keys() :
        pointsFrom.append( pts1[ res1[ri][' CA '] ] )
    pointsFrom = VecVecFloat(pointsFrom)

    ofp = open(options.multipdb, 'r') ; multilines = ofp.readlines() ; ofp.close()
    ofp = open(options.outpdb, 'w')
    rot = VecVecFloat([ [0,0,0], [0,0,0], [0,0,0] ])
    t1,t2 = VecFloat( [0,0,0] ), VecFloat( [0,0,0] )
    modelLines = None
    for l in multilines :
        if l[0:5] == "MODEL" : modelLines = [l]
        elif l[0:6] == "ENDMDL" :
            modelLines.append(l)
            prot = protein(modelLines, read_hydrogens=0, read_waters=0, read_hets=0)
            res2, resids2, resnums2, resns2, chids2, inscodes2, pts2 = prot2res.readProtRes(prot)
            pointsOnto = []
            for ri in res1.keys() :
                assert resids1[ri] == resids2[ri]
                pointsOnto.append( pts2[ res2[ri][' CA '] ] )
            pointsOnto = VecVecFloat(pointsOnto)
            rmsd = findSuperpositionTransform(pointsFrom, pointsOnto, t1, rot, t2)
            print "RMSD", rmsd
            for l in modelLines :
                if "RMSD" in l : continue
                if l[0:5] == "MODEL" : ofp.write(l) ; ofp.write("REMARK         RMSD %f\n" % rmsd)
                else : ofp.write(l)
        elif modelLines != None : modelLines.append(l)
    ofp.close()

def supermain() :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb1", action='store', type='string', dest='pdbfile1', help='transform coordinates in this PDB file', default=None)
    parser.add_option("--pdb2", action='store', type='string', dest='pdbfile2', help='to superpose onto coordinates in this PDB file', default=None)
    parser.add_option("--outpdb", action='store', type='string', dest='outpdb', help='and write resulting coordinates into this file in PDB format', default=None)
    (options, args) = parser.parse_args()

    prot = protein(options.pdbfile1, read_hydrogens=1, read_waters=1, read_hets=1)
    res1, resids1, resnums1, resns1, chids1, inscodes1, pts1 = prot2res.readProtRes(prot)
    prot = protein(options.pdbfile2, read_hydrogens=1, read_waters=1, read_hets=1)
    res2, resids2, resnums2, resns2, chids2, inscodes2, pts2 = prot2res.readProtRes(prot)

    pointsFrom, pointsOnto = [], []
    for ri in res1.keys() :
        assert resids1[ri] == resids2[ri]
        pointsFrom.append( pts1[ res1[ri][' CA '] ] )
        pointsOnto.append( pts2[ res2[ri][' CA '] ] )
    pointsFrom, pointsOnto = VecVecFloat(pointsFrom), VecVecFloat(pointsOnto)

    rot = VecVecFloat([ [0,0,0], [0,0,0], [0,0,0] ])
    t1,t2 = VecFloat( [0,0,0] ), VecFloat( [0,0,0] )
    rmsd = findSuperpositionTransform(pointsFrom, pointsOnto, t1, rot, t2)
    print "RMSD", rmsd

    pts1 = VecVecFloat(pts1)
    TRTtransform(pts1, t1,rot,t2)

    from peptidebuild import ModelRenderer
    ModelRenderer(res1, resns1, chids1, resnums1, inscodes1, [], options.outpdb).render(pts1)

def supermain_test() :
    pointsFrom = [[16.674, 23.342, 23.285], [15.705, 23.079, 24.32], [15.069, 21.706, 24.131]]
    pointsOnto = [[-9.81934, 21.5857, -6.42715], [-8.59078, 22.358, -6.28567], [-8.3812, 22.7826, -4.83613]]

    #pointsFrom = [ [10,50,0], [0,0,0], [-10,20,0] ]
    #pointsOnto = [ [-10,50,0], [0,0,0], [10,20,0] ]
    fromP = VecVecFloat(pointsFrom)
    ontoP = VecVecFloat(pointsOnto)

    for pi in range(fromP.size()) : print "BEFORE", nicestr(fromP[pi]), nicestr(ontoP[pi])
    rot = VecVecFloat([ [0,0,0], [0,0,0], [0,0,0] ])
    t1,t2 = VecFloat( [0,0,0] ), VecFloat( [0,0,0] )
    rmsd = findSuperpositionTransform(fromP, ontoP, t1, rot, t2)
    print "RMSD", rmsd
    print "T1", t1[0], t1[1], t1[2]
    print "rot", rot[0][0], rot[0][1], rot[0][2]
    print "rot", rot[1][0], rot[1][1], rot[1][2]
    print "rot", rot[2][0], rot[2][1], rot[2][2]
    print "T2", t2[0], t2[1], t2[2]
    for pi in range(fromP.size()) : print "AFTER ", nicestr(fromP[pi]), nicestr(ontoP[pi])
    TRTtransform(fromP, t1, rot, t2)
    for pi in range(fromP.size()) : print "AFTER ", nicestr(fromP[pi]), nicestr(ontoP[pi])

if __name__ == "__main__" : superMulti(); #supermain()
