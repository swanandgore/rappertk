import os, sys
cbDatapath = os.environ["RTKROOT"] + "/data/"
from data import resAtoms, consts
from samplers import PhipsiSampler, OmegaSampler
from builders import NanchorBuilder, VecVecFloat, VecInt, VecFloat
from geometry import calcDist, calcAngle, calcDihed
from peptidebuild import ModelRenderer

if __name__ == "__main__" :
    IPs, OPs = [], []

    res, resns, resnums, chids, inscodes = {}, {}, {}, {}, {}
    pts = []

    resns[0] = "THR" ; resnums[0] = "  92" ; chids[0] = "A" ; inscodes[0] = ' ' ; res[0] = {}
    pts.append( [0,0,0] ) ; res[0][' N  '] = len(pts)-1
    pts.append( [0,0,0] ) ; res[0][' CA '] = len(pts)-1
    pts.append( [0,0,0] ) ; res[0][' C  '] = len(pts)-1
    pts.append( [0,0,0] ) ; res[0][' O  '] = len(pts)-1

    resns[1] = "THR" ; resnums[1] = "  93" ; chids[1] = "A" ; inscodes[1] = ' ' ; res[1] = {}
    pts.append( [0,0,0] ) ; res[1][' N  '] = len(pts)-1
    pts.append( [0,0,0] ) ; res[1][' CA '] = len(pts)-1

    pts = VecVecFloat(pts)

    OPs = [
            res[0][' N  '], res[0][' CA '], res[0][' C  '],
            res[0][' O  '], res[0][' N  '], res[0][' CA '],
        ]
    nab = NanchorBuilder(VecInt(IPs), VecInt(OPs), consts, "NanchorBuilder",
            VecFloat([0,0,0]), .5, VecFloat([3,0,0]), .2)

    nab.build(pts)
    mr = ModelRenderer(res, resns, chids, resnums, inscodes)
    mr.render(pts)
